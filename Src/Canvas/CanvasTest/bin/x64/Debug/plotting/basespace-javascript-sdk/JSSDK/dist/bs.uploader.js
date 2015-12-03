/// <reference path="bs.sdk.ts" />

var bs;
(function (bs) {
    

    var SliceJob = (function () {
        function SliceJob(slice) {
            this.slice = slice;
            this.xhr = new XMLHttpRequest();
        }
        return SliceJob;
    })();
    bs.SliceJob = SliceJob;

    var FileMetaData = (function () {
        function FileMetaData() {
            this.parts = {};
            this.bytesUploaded = 0;
        }
        return FileMetaData;
    })();
    bs.FileMetaData = FileMetaData;

    var SliceUploaderManager = (function () {
        function SliceUploaderManager(customOptions) {
            var _this = this;
            this.retryCount = 0;
            this.queue = [];
            this.filesInProgress = [];
            this.totalSlices = 0;
            this.totalSlicesUploaded = 0;
            this.onUploadSliceRetryCallBack = function (event) {
                setTimeout(function () {
                    if (event.retryCount < SliceUploaderManager.MAX_RETRY) {
                        event.retryCount++;
                        _this.uploadSlice(event.sliceJob, event.fileMetaData, event.retryCount);
                        console.log('retrying slice upload');
                    } else {
                        $(_this).trigger("failed");
                        console.log('retry ' + SliceUploaderManager.MAX_RETRY + ' and failed');
                    }
                }, SliceUploaderManager.RETRY_TIMEOUT * event.retryCount);
            };
            this.onSliceUploadComplete = function (event, responseText) {
                var partNumber = event.partResponse.Response.Number;
                var part = event.fileMetaData.parts[partNumber];

                event.fileMetaData.bytesUploaded += part.size;

                _this.removeCompletedSlice(partNumber);
                _this.totalSlicesUploaded++;
                $(_this).trigger("progressChanged", responseText);

                if (event.fileMetaData.bytesUploaded == event.fileMetaData.totalSize) {
                    _this.onMultiPartUploadComplete(part, responseText);
                }
            };
            this.bsApiCustomOptions = customOptions;
        }
        SliceUploaderManager.prototype.addSliceToQueue = function (slice, fileMetaData) {
            this.queue.push(slice);
            this.fileMetaData = fileMetaData;
            this.totalSlices++;

            this.ensureTimerRunning();
        };

        SliceUploaderManager.prototype.ensureTimerRunning = function () {
            var _this = this;
            if (this.queueTimer == null) {
                this.queueTimer = setInterval(function () {
                    _this.onQueueTimer();
                }, 1000);
            }
        };

        SliceUploaderManager.prototype.stopTimerIfNoWork = function () {
            if (this.queue.length == 0) {
                clearInterval(this.queueTimer);
                this.queueTimer = null;
            }
        };

        SliceUploaderManager.prototype.onQueueTimer = function () {
            this.updateQueue();
        };

        SliceUploaderManager.prototype.updateQueue = function () {
            if (this.filesInProgress.length < SliceUploaderManager.MAX_CONCURRENT_UPLOADS && this.queue.length > 0) {
                var slice = this.queue.pop();
                var sliceJob = new SliceJob(slice);
                this.filesInProgress.push(sliceJob);

                sliceJob.slice = this.fileMetaData.parts[slice.partNumber];
                this.uploadSlice(sliceJob, this.fileMetaData);
            }

            this.stopTimerIfNoWork();
        };

        SliceUploaderManager.prototype.removeCompletedSlice = function (partNumber) {
            var _this = this;
            this.filesInProgress.forEach(function (sliceJob, index) {
                if (sliceJob.slice.partNumber == partNumber) {
                    _this.filesInProgress = _this.filesInProgress.slice(index, 1);
                }
            });
        };

        SliceUploaderManager.prototype.cancelUpload = function () {
            this.queue = [];
            var deferred = $.Deferred();
            var xhrs = [];
            this.filesInProgress.forEach(function (sliceJob) {
                xhrs.push(sliceJob.xhr.abort());
            });

            $.when(xhrs).always(function () {
                deferred.resolve();
            });

            this.filesInProgress = [];

            return deferred;
        };

        SliceUploaderManager.prototype.uploadSlice = function (sliceJob, fileMetaData, retryCount) {
            var _this = this;
            if (typeof retryCount === "undefined") { retryCount = 0; }
            var xhr = sliceJob.xhr;

            xhr.open("PUT", this.bsApiCustomOptions.baseUri + "files/" + sliceJob.slice.id + "/parts/" + sliceJob.slice.partNumber, true);
            xhr.withCredentials = true;
            xhr.setRequestHeader("Content-Type", sliceJob.slice.type);

            sliceJob.xhr.onreadystatechange = function () {
                _this.handleUploadSliceStateChange(xhr.responseText, xhr.readyState, xhr.status, sliceJob, fileMetaData, retryCount, _this.onSliceUploadComplete, _this.onUploadSliceRetryCallBack);
            };

            sliceJob.xhr.send(sliceJob.slice);
        };

        SliceUploaderManager.prototype.handleUploadSliceStateChange = function (responseText, readyState, status, sliceJob, fileMetaData, retryCount, onCompleteCallBack, onRetryCallBack) {
            try  {
                if (readyState === 4) {
                    var partResponse = JSON.parse(responseText);
                    var partNumber = partResponse.Response.Number;
                    var part = fileMetaData.parts[partNumber];

                    if (status === 200) {
                        onCompleteCallBack({ partResponse: partResponse, fileMetaData: fileMetaData });
                    } else {
                        onRetryCallBack({ sliceJob: sliceJob, fileMetaData: fileMetaData, retryCount: retryCount });
                    }
                }
            } catch (ex) {
                onRetryCallBack({ sliceJob: sliceJob, fileMetaData: fileMetaData, retryCount: retryCount });
                console.log(ex);
            }
        };

        SliceUploaderManager.prototype.onMultiPartUploadComplete = function (part, responseText) {
            var _this = this;
            $.ajax({
                type: "post",
                url: this.bsApiCustomOptions.baseUri + "files/" + part.id + "?uploadstatus=complete",
                xhrFields: {
                    withCredentials: $.support.cors
                },
                crossDomain: $.support.cors,
                dataType: "application/json"
            }).always(function (updateFileResponse) {
                $(_this).trigger("multiPartUploadComplete", [updateFileResponse, updateFileResponse.responseText]);
                $(_this).trigger("progressChanged", [updateFileResponse, updateFileResponse.responseText]);
            });
        };
        SliceUploaderManager.MAX_CONCURRENT_UPLOADS = 5;
        SliceUploaderManager.MAX_RETRY = 5;
        SliceUploaderManager.RETRY_TIMEOUT = 2000;
        return SliceUploaderManager;
    })();
    bs.SliceUploaderManager = SliceUploaderManager;

    var Uploader = (function () {
        function Uploader(customOptions) {
            var _this = this;
            this.isMultiPart = false;
            this.fileMetaData = {};
            this.onSingleFileUploadComplete = function (responseText) {
                _this.onProgressChangedCallBack(responseText);
                _this.onCompleteCallBack(responseText);
            };
            this.onSingleFileUploadRetry = function (responseText) {
            };
            this.fileAPISupportCheck();
            this.bsApi = new bs.Api(customOptions);
            this.bsApiCustomOptions = customOptions;
            this.sliceUploaderManager = new SliceUploaderManager(customOptions);

            $(this.sliceUploaderManager).on("progressChanged", function (event, responseText) {
                _this.onProgressChangedCallBack(responseText);
            });
            $(this.sliceUploaderManager).on("multiPartUploadComplete", function (event) {
                var xhr = [];
                for (var _i = 0; _i < (arguments.length - 1); _i++) {
                    xhr[_i] = arguments[_i + 1];
                }
                _this.onCompleteCallBack(xhr[0].responseText);
            });
            $(this.sliceUploaderManager).on("failed", function (event) {
                var xhr = [];
                for (var _i = 0; _i < (arguments.length - 1); _i++) {
                    xhr[_i] = arguments[_i + 1];
                }
                _this.onErrorCallBack(xhr[0].responseText);
            });
        }
        Uploader.prototype.fileAPISupportCheck = function () {
            if ((typeof (File) !== 'undefined') && (typeof (Blob) !== 'undefined') && (typeof (FileList) !== 'undefined') && (!!Blob.prototype.slice || false))
                return true;
            return false;
        };

        Uploader.prototype.uploadFileToAppResult = function (appresultId, file, callBack) {
            var fileMetaData = new FileMetaData();

            fileMetaData.totalSize = file.size;
            fileMetaData.appResultId = appresultId;
            this.fileMetaData[file.name] = fileMetaData;
            this.callBack = callBack;

            // DT: where is filetype coming from?
            //get the file type (if available)
            var filetype = filetype || 'application/octet-stream';
            var megabyte = 1024 * 1024;
            var minFileSize = file.size < 100 * megabyte ? 5 * megabyte : 25 * megabyte;

            if (minFileSize <= file.size) {
                this.multiPartUpload(file, filetype, minFileSize, appresultId);
            } else {
                this.singleFileUpload(file, filetype, appresultId);
            }
        };

        Uploader.prototype.singleFileUpload = function (file, fileType, appresultId) {
            var _this = this;
            this.isMultiPart = false;
            this.singleFileXhr = new XMLHttpRequest();
            this.singleFileXhr.withCredentials = true;
            this.singleFileXhr.open("POST", this.bsApiCustomOptions.baseUri + "/appresults/" + appresultId + "/files?name=" + encodeURI(file.name) + "&multipart=false", true);
            this.singleFileXhr.setRequestHeader("Content-Type", fileType);
            this.singleFileXhr.send(file);
            this.singleFileXhr.onreadystatechange = function () {
                _this.onSingleFileUploadStateChange(_this.singleFileXhr.responseText, _this.singleFileXhr.readyState, _this.singleFileXhr.status, _this.onSingleFileUploadComplete, _this.onSingleFileUploadRetry);
            };
        };

        Uploader.prototype.onSingleFileUploadStateChange = function (responseText, readyState, status, onSingleFileCompleteCallBack, onSingleFileRetryCallBack) {
            if (readyState === 4) {
                if (status === 201) {
                    onSingleFileCompleteCallBack(responseText);
                } else {
                    onSingleFileRetryCallBack(responseText);
                }
            }
        };

        Uploader.prototype.cancelUpload = function () {
            var _this = this;
            var deferred = $.Deferred();
            if (this.isMultiPart) {
                $.ajax({
                    type: "delete",
                    url: this.bsApiCustomOptions.baseUri + "files/" + this.fileId,
                    xhrFields: {
                        withCredentials: $.support.cors
                    },
                    crossDomain: $.support.cors,
                    dataType: "application/json"
                }).always(function () {
                    _this.sliceUploaderManager.cancelUpload().always(function () {
                        deferred.resolve();
                    });
                });
            } else {
                if (this.singleFileXhr != null)
                    this.singleFileXhr.abort();
                deferred.resolve();
            }
            return deferred;
        };

        Uploader.prototype.multiPartUpload = function (file, fileType, minFileSize, appresultId) {
            var _this = this;
            this.isMultiPart = true;

            $.ajax({
                type: "post",
                url: this.bsApiCustomOptions.baseUri + "/appresults/" + appresultId + "/files?name=" + encodeURI(file.name) + "&multipart=true",
                xhrFields: {
                    withCredentials: $.support.cors
                },
                crossDomain: $.support.cors,
                headers: {
                    "Content-Type": fileType
                }
            }).then(function (fileResponse) {
                var response = fileResponse;
                var fileId = response.Response.Id;
                _this.fileId = fileId;

                var fileSlices = _this.sliceFiles(file, fileId, fileType, minFileSize);
                var fileMeta = _this.fileMetaData[file.name];

                fileSlices.forEach(function (fileSlice) {
                    fileMeta.parts[fileSlice.partNumber] = fileSlice;
                    _this.sliceUploaderManager.addSliceToQueue(fileSlice, fileMeta);
                });
            });
        };

        Uploader.prototype.sliceFiles = function (file, fileId, fileType, minFileSize) {
            var slices = new Array();

            //get the slices for this file
            var chunkSize = minFileSize;
            var beginPosition = 0;
            var endPosition = chunkSize;
            var partNumber = 1;

            do {
                var fileslice = file.slice(beginPosition, endPosition, fileType);
                fileslice.partNumber = partNumber;
                fileslice.id = fileId;
                fileslice.name = file.name;
                fileslice.type = fileType;

                slices.push(fileslice);

                beginPosition = endPosition;
                endPosition += chunkSize;
                partNumber++;
            } while(beginPosition < file.size);

            return slices;
        };

        Uploader.prototype.onProgressChangedCallBack = function (responseText) {
            if (this.isMultiPart) {
                var percentComplete = this.sliceUploaderManager.totalSlicesUploaded / this.sliceUploaderManager.totalSlices;
                this.callBack.onProgress(responseText, percentComplete, this.sliceUploaderManager.totalSlices, this.sliceUploaderManager.totalSlicesUploaded);
            } else {
                // single file will return 1 for percent complete since there are no slices
                this.callBack.onProgress(responseText, 1);
            }
        };

        Uploader.prototype.onCompleteCallBack = function (responseText) {
            if (this.isMultiPart)
                this.callBack.onComplete(responseText, this.sliceUploaderManager.totalSlices, this.sliceUploaderManager.totalSlicesUploaded);
            else
                this.callBack.onComplete(responseText);
        };

        Uploader.prototype.onErrorCallBack = function (responseText) {
            this.callBack.onError(responseText);
        };
        return Uploader;
    })();
    bs.Uploader = Uploader;
})(bs || (bs = {}));
