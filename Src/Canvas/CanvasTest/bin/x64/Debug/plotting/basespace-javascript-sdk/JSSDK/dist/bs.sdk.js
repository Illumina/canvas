/// <reference path="../lib/jquery.d.ts" />
/// <reference path="../lib/moment.d.ts" />
var bs;
(function (bs) {
    bs.xhrCache = bs.xhrCache || {};

    var Api = (function () {
        function Api(customOptions) {
            /*
            *    This isis code should be removed after isis is made a native app
            */
            this.fetchIsisAnalysisLowPercentageChart = function (analysisId, successCallback, failCallback, queryParameters) {
                this.apiRequest("GET", "isis-analysis/" + analysisId + "/low-percentage/data", queryParameters, successCallback, failCallback);
            };
            this.fetchIsisAnalysisMismatchChart = function (analysisId, successCallback, failCallback, queryParameters) {
                this.apiRequest("GET", "isis-analysis/" + analysisId + "/mismatch/data", queryParameters, successCallback, failCallback);
            };
            this.fetchIsisAnalysisHighPercentageChart = function (analysisId, successCallback, failCallback, queryParameters) {
                this.apiRequest("GET", "isis-analysis/" + analysisId + "/high-percentage/data", queryParameters, successCallback, failCallback);
            };
            this.fetchIsisAnalysisClustersChart = function (analysisId, successCallback, failCallback, queryParameters) {
                this.apiRequest("GET", "isis-analysis/" + analysisId + "/clusters/data", queryParameters, successCallback, failCallback);
            };
            this.fetchIsisAnalysisTrimmedLengthChart = function (analysisId, successCallback, failCallback, queryParameters) {
                this.apiRequest("GET", "isis-analysis/" + analysisId + "/trimmed-length/data", queryParameters, successCallback, failCallback);
            };
            this.checkCors();

            this.options = {
                baseUri: "https://api.basespace.illumina.com/v1pre3/",
                accessToken: null,
                baseUriProxyFallback: "https://api.basespace.illumina.com/api/proxy/",
                useXhrCache: false,
                cacheTimeout: Api.DEFAULT_CACHE_TIMEOUT
            };

            if (customOptions) {
                $.extend(true, this.options, customOptions);
            }

            // default to cookie-based authentication.
            // if token is provided, use token for authentication.
            if (!this.options.accessToken) {
                this.authType = Api.COOKIE_CREDENTIALS;
            } else {
                this.authType = Api.ACCESS_TOKEN;
            }
        }
        Api.prototype.checkCors = function () {
            var xhr = new XMLHttpRequest();
            $.support.cors = ("withCredentials" in xhr);
            xhr = null;
            if (!$.support.cors)
                console.log("CORS isn't supported.");
        };

        Api.prototype.getBaseUri = function () {
            return this.options.baseUri;
        };
        Api.prototype.setBaseUri = function (s) {
            this.options.baseUri = s;
        };

        Api.prototype.fetchApplication = function (appId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "applications/" + appId, queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchAppSessions = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "users/current/appsessions", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchAppSessionsByProject = function (projectId, successCallback, failCallback, queryParameters) {
            queryParameters = $.extend(true, { "output.projects": projectId }, queryParameters);
            this.apiRequest("GET", "users/current/appsessions", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchAppSession = function (appSessionId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "appsessions/" + appSessionId, queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchCurrentUser = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "users/current", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchProjects = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "users/current/projects", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchProjectName = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "users/current/projects", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchProject = function (projectId, successCallback, failCallback) {
            this.apiRequest("GET", "projects/" + projectId, null, successCallback, failCallback);
        };

        Api.prototype.fetchProjectSamples = function (projectId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "projects/" + projectId + "/samples", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchAppSessionsDeep = function (successCallback, failCallback, queryParameters) {
            var _this = this;
            var fetchChildrenListMethod = function (success, fail, queryParams) {
                return _this.fetchAppSessions(success, fail, queryParams);
            };
            var fetchChildDetailsMethod = this.fetchAppSession;
            return this.apiRequestDeep(fetchChildrenListMethod, fetchChildDetailsMethod, successCallback, failCallback, queryParameters);
        };

        Api.prototype.fetchProjectSamplesDeep = function (projectId, successCallback, failCallback, queryParameters) {
            var _this = this;
            var fetchChildrenListMethod = function (success, fail, queryParams) {
                return _this.fetchProjectSamples(projectId, success, fail, queryParams);
            };
            var fetchChildDetailsMethod = this.fetchSample;
            return this.apiRequestDeep(fetchChildrenListMethod, fetchChildDetailsMethod, successCallback, failCallback, queryParameters);
        };

        Api.prototype.fetchProjectAppResultsDeep = function (projectId, successCallback, failCallback, queryParameters) {
            var _this = this;
            var fetchChildrenListMethod = function (success, fail, queryParams) {
                return _this.fetchProjectAppResults(projectId, success, fail, queryParams);
            };
            var fetchChildDetailsMethod = this.fetchAppResult;
            return this.apiRequestDeep(fetchChildrenListMethod, fetchChildDetailsMethod, successCallback, failCallback, queryParameters);
        };

        Api.prototype.fetchGenomesDeep = function (successCallback, failCallback, queryParameters) {
            return this.apiRequestDeep(this.fetchGenomes, this.fetchGenome, successCallback, failCallback, queryParameters);
        };

        Api.prototype.fetchProjectAppResults = function (projectId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "projects/" + projectId + "/appresults", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchProjectFiles = function (projectId, successCallback, failCallback, queryParameters) {
            var _this = this;
            this.apiRequest("GET", "projects/" + projectId + "/appresults", {}, function (json) {
                var appResults = json.Response.Items;
                var queriedSoFar = 0;
                var appResultsToQuery = appResults.length;

                var combinedResponse = {
                    Items: []
                };

                if (appResultsToQuery < 1) {
                    successCallback({ 'Response': combinedResponse });
                }

                for (var i = 0; i < appResults.length; i++) {
                    var appResult = appResults[i];

                    _this.fetchAppResultFiles(appResult.Id, function (json) {
                        var appResultFilesResponse = json.Response;
                        var currentItems = appResultFilesResponse.Items;

                        for (var j = 0; j < currentItems.length; j++) {
                            var extraDetailsAboutThisFile = {
                                AppResultName: appResult.Name,
                                AppResultId: appResult.Id,
                                AppResultStatus: appResult.Status,
                                UserOwnedBy: appResult.UserOwnedBy
                            };

                            $.extend(currentItems[j], extraDetailsAboutThisFile);
                        }
                        combinedResponse.Items = combinedResponse.Items.concat(currentItems);
                        queriedSoFar++;

                        if (queriedSoFar == appResultsToQuery) {
                            successCallback({ 'Response': combinedResponse });
                        }
                    }, failCallback, queryParameters);
                }
            }, failCallback);
        };

        Api.prototype.fetchSample = function (sampleId, successCallback, failCallback) {
            this.apiRequest("GET", "samples/" + sampleId, null, successCallback, failCallback);
        };

        Api.prototype.fetchSampleFiles = function (sampleId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "samples/" + sampleId + "/files", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchSampleAppResultFiles = function (sampleId, successCallback, failCallback, queryParameters) {
            var _this = this;
            this.apiRequest("GET", "samples/" + sampleId, {}, function (json) {
                var response = json.Response;
                var references = response.References;
                var queriedSoFar = 0;
                var appResultsToQuery = references.length;

                var combinedResponse = {
                    Items: []
                };

                if (appResultsToQuery < 1) {
                    successCallback({ 'Response': combinedResponse });
                }

                for (var i = 0; i < references.length; i++) {
                    var reference = references[i];

                    if (reference.Rel == Api.REFERENCE_REL_USED_BY && reference.Type == Api.REFERENCE_TYPE_APP_RESULT) {
                        var appResultId = reference.Content.Id;
                        var appResultName = reference.Content.Name;
                        var appResultStatus = reference.Content.Status;
                        var appResultUserOwnedBy = reference.Content.UserOwnedBy;

                        _this.fetchAppResultFiles(appResultId, function (json) {
                            var appResultFilesResponse = json.Response;
                            var currentItems = appResultFilesResponse.Items;

                            for (var j = 0; j < currentItems.length; j++) {
                                var extraDetailsAboutThisFile = {
                                    AppResultName: appResultName,
                                    AppResultId: appResultId,
                                    AppResultStatus: appResultStatus,
                                    UserOwnedBy: appResultUserOwnedBy
                                };

                                $.extend(currentItems[j], extraDetailsAboutThisFile);
                            }
                            combinedResponse.Items = combinedResponse.Items.concat(currentItems);
                            queriedSoFar++;

                            if (queriedSoFar == appResultsToQuery) {
                                successCallback({ 'Response': combinedResponse });
                            }
                        }, failCallback, queryParameters);
                    }
                }
            }, failCallback);
        };

        Api.prototype.fetchAppResult = function (appResultId, successCallback, failCallback) {
            this.apiRequest("GET", "appresults/" + appResultId, null, successCallback, failCallback);
        };

        Api.prototype.fetchAppResultsByAppSession = function (appSessionId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "appsessions/" + appSessionId + "/appresults", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchSampleSource = function (sampleSourceId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "sources/" + sampleSourceId, queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchSampleSources = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "users/current/sources", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchLibraries = function (sampleSourceId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "sources/" + sampleSourceId + "/libraries", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchLibrary = function (sampleLibraryId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "libraries/" + sampleLibraryId, queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchAppResultFiles = function (appResultId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "appresults/" + appResultId + "/files", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchGenomes = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "genomes", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchGenome = function (genomeId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "genomes/" + genomeId, null, successCallback, failCallback);
        };

        Api.prototype.fetchRuns = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "users/current/runs", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchRun = function (runId, successCallback, failCallback) {
            this.apiRequest("GET", "runs/" + runId, null, successCallback, failCallback);
        };

        Api.prototype.fetchFile = function (fileId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "files/" + fileId, queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchFileContent = function (fileId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "files/" + fileId + "/content", queryParameters, successCallback, failCallback);
        };

        /*
        QC plots
        */
        Api.prototype.fetchRunCyclePlot = function (runId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "runs/" + runId + "/cycleplot/data", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchRunLanePlot = function (runId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "runs/" + runId + "/laneplot/data", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchRunQScoreHistogram = function (runId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "runs/" + runId + "/qscorehistogram/data", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchRunQScoreHeatmap = function (runId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "runs/" + runId + "/qscoreheatmap/data", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchRunFlowcellChart = function (runId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "runs/" + runId + "/flowcellchart/data", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchRunSampleQCChart = function (runId, successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "runs/" + runId + "/sampleqcchart/data", queryParameters, successCallback, failCallback);
        };

        Api.prototype.fetchSpecies = function (successCallback, failCallback, queryParameters) {
            this.apiRequest("GET", "species", queryParameters, successCallback, failCallback);
        };

        // ***************************
        Api.prototype.createProject = function (queryParameters, successCallback, failCallback) {
            var maxNameLength = 100;
            if (queryParameters.Name.length > maxNameLength) {
                queryParameters.Name = queryParameters.Name.substring(0, maxNameLength);
            }

            this.apiRequest("POST", "projects", queryParameters, successCallback, failCallback);
        };

        // ****************************
        Api.prototype.search = function (request, queryOverride, successCallBack, failCallback) {
            if (queryOverride != undefined)
                request.Query = queryOverride;

            if (!request.Query)
                request.Query = "*";

            this.apiRequest("GET", "search", request, successCallBack, failCallback);
        };

        Api.prototype.makePlural = function (entityType) {
            var plural;
            switch (entityType) {
                case Api.SAMPLELIBRARY:
                    plural = "libraries";
                    break;
                default:
                    plural = entityType + "s";
            }
            return plural;
        };

        Api.prototype.fetchProperties = function (entityType, entityId, queryParameters, successCallback, failCallback) {
            // make plural with "+s" might not work with library
            this.apiRequest("GET", this.makePlural(entityType) + "/" + entityId + "/properties", queryParameters, successCallback, failCallback);
        };

        Api.prototype.isCachedXhrExpired = function (url) {
            if (bs.xhrCache.hasOwnProperty(url)) {
                var now = new Date();
                var lastFetched = bs.xhrCache[url].lastFetched;

                if (now - lastFetched < this.options.cacheTimeout)
                    return false;
            }

            return true;
        };

        Api.prototype.apiRequest = function (method, path, data, successCallback, failCallback) {
            var url, dataType, xhr, urlWithQueryParams;

            if ($.support.cors) {
                url = this.options.baseUri + path;
            } else if (!$.support.cors && this.options.baseUriProxyFallback != null) {
                url = this.options.baseUriProxyFallback + path;
            } else {
                $.error("CORS not supported, and no fallback provided.");
                return;
            }

            urlWithQueryParams = url + (data ? "?" + $.param(data) : "");

            // if useXhrCache is enabled on this instance of the Api client, try to get
            // the cached xhr object, and attach callbacks to it, instead of firing
            // a new ajax call.
            if (this.options.useXhrCache && bs.xhrCache.hasOwnProperty(urlWithQueryParams) && !this.isCachedXhrExpired(urlWithQueryParams)) {
                xhr = bs.xhrCache[urlWithQueryParams].xhr;

                xhr.then(function () {
                    var reasons = [];
                    for (var _i = 0; _i < (arguments.length - 0); _i++) {
                        reasons[_i] = arguments[_i + 0];
                    }
                    var data = reasons[0], textStatus = reasons[1], jqXHR = reasons[2];

                    successCallback(data, textStatus, jqXHR);
                }, function () {
                    var reasons = [];
                    for (var _i = 0; _i < (arguments.length - 0); _i++) {
                        reasons[_i] = arguments[_i + 0];
                    }
                    var jqXHR = reasons[0], textStatus = reasons[1], errorThrown = reasons[2];

                    failCallback(jqXHR, textStatus, errorThrown);
                });

                return xhr;
            }

            var ajaxSettings = {
                type: method,
                url: url,
                xhrFields: {
                    withCredentials: $.support.cors
                },
                crossDomain: $.support.cors,
                data: data,
                dataType: "json"
            };

            if (this.authType == Api.COOKIE_CREDENTIALS) {
                xhr = $.ajax(ajaxSettings);

                xhr.then(function () {
                    var reasons = [];
                    for (var _i = 0; _i < (arguments.length - 0); _i++) {
                        reasons[_i] = arguments[_i + 0];
                    }
                    var data = reasons[0], textStatus = reasons[1], jqXHR = reasons[2];

                    successCallback(data, textStatus, jqXHR);
                }, function () {
                    var reasons = [];
                    for (var _i = 0; _i < (arguments.length - 0); _i++) {
                        reasons[_i] = arguments[_i + 0];
                    }
                    var jqXHR = reasons[0], textStatus = reasons[1], errorThrown = reasons[2];

                    failCallback(jqXHR, textStatus, errorThrown);
                });

                // this is the first time making this ajax request, so cache it.
                if (this.options.useXhrCache) {
                    bs.xhrCache[urlWithQueryParams] = {
                        xhr: xhr,
                        lastFetched: new Date()
                    };
                }
            }

            //else if (this.authType == Api.ACCESS_TOKEN) {
            //    //TODO: Token-based Api calls, to support cross-domain requests.
            //}
            return xhr;
        };

        Api.prototype.apiRequestDeep = function (fetchChildrenListMethod, fetchChildDetailsMethod, successCallback, failCallback, queryParameters) {
            var _this = this;
            fetchChildrenListMethod.call(this, function (json) {
                var response = json.Response;
                var children = response.Items;
                var queriedSoFar = 0;
                var numToQuery = children.length;
                var combinedResponse = {
                    Items: children,
                    DisplayedCount: response.DisplayedCount,
                    TotalCount: response.TotalCount,
                    Offset: response.Offset,
                    Limit: response.Limit,
                    SortDir: response.SortDir,
                    SortBy: response.SortBy
                };

                if (numToQuery < 1) {
                    successCallback(json);
                    return;
                }

                for (var i = 0; i < numToQuery; i++) {
                    fetchChildDetailsMethod.call(_this, children[i].Id, function (json) {
                        var childDeepResponse = json.Response;
                        var childId = childDeepResponse.Id;
                        for (var j = 0; j < combinedResponse.Items.length; j++) {
                            if (combinedResponse.Items[j].Id == childId) {
                                combinedResponse.Items[j] = childDeepResponse;
                                console.log(childDeepResponse);
                                break;
                            }
                        }
                        queriedSoFar++;
                        if (queriedSoFar == numToQuery) {
                            successCallback({ 'Response': combinedResponse });
                        }
                    }, failCallback);
                }
            }, failCallback, queryParameters);
        };

        Api.prototype.filterResponseItems = function (items, json) {
            var matches = [];

            for (var i = 0; i < items.length; i++) {
                var item = items[i];

                var mergedObj = {};

                $.extend(mergedObj, item, json);

                if (mergedObj == item) {
                    matches.push[item];
                }
            }

            return matches;
        };
        Api.PROJECT = 'project';
        Api.PROJECTFILE = 'projectfile';
        Api.SAMPLE = 'sample';
        Api.APPSESSION = 'appsession';
        Api.APPRESULT = 'appresult';
        Api.SAMPLEFILE = 'samplefile';
        Api.FILE = 'file';
        Api.SAMPLESOURCE = 'samplesource';
        Api.SAMPLELIBRARY = 'library';
        Api.GENOME = 'genome';
        Api.RUN = 'run';
        Api.SPECIES = 'species';

        Api.PERMISSIONS_OWNER = 'owner';
        Api.PERMISSIONS_COLLABORATOR = 'collaborator';

        Api.COOKIE_CREDENTIALS = 'credentials';
        Api.ACCESS_TOKEN = 'access-token';

        Api.REFERENCE_REL_USED_BY = 'UsedBy';
        Api.REFERENCE_TYPE_APP_RESULT = 'AppResult';

        Api.DEFAULT_CACHE_TIMEOUT = 15000;
        return Api;
    })();
    bs.Api = Api;

    var SearchScopes = (function () {
        function SearchScopes() {
        }
        SearchScopes.fromResourceType = function (resourceType) {
            switch (resourceType.toLowerCase()) {
                case Api.PROJECT:
                    return this.PROJECTS;
                case Api.APPRESULT:
                    return this.APPRESULTS;
                case Api.SAMPLE:
                    return this.SAMPLES;
                case Api.APPSESSION:
                    return this.APPSESSIONS;
                case Api.RUN:
                    return this.RUNS;
                case Api.FILE:
                    return this.APPRESULT_FILES;
            }
        };
        SearchScopes.ALL = "ALL";
        SearchScopes.RUNS = "RUNS";
        SearchScopes.SAMPLES = "SAMPLES";
        SearchScopes.APPRESULTS = "APPRESULTS";
        SearchScopes.APPS = "APPS";
        SearchScopes.PROJECTS = "PROJECTS";
        SearchScopes.APPSESSIONS = "APPSESSIONS";
        SearchScopes.SAMPLE_FILES = "SAMPLE_FILES";
        SearchScopes.APPRESULT_FILES = "APPRESULT_FILES";
        SearchScopes.USERS = "USERS";
        return SearchScopes;
    })();
    bs.SearchScopes = SearchScopes;

    var SearchResourceTypes = (function () {
        function SearchResourceTypes() {
        }
        SearchResourceTypes.fromResourceType = function (resourceType) {
            switch (resourceType.toLowerCase()) {
                case Api.PROJECT:
                    return this.PROJECT;
                case Api.APPRESULT:
                    return this.APPRESULT;
                case Api.SAMPLE:
                    return this.SAMPLE;
                case Api.APPSESSION:
                    return this.APPSESSION;
                case Api.RUN:
                    return this.RUN;
                case Api.FILE:
                    return this.APPRESULT_FILE;
            }
        };
        SearchResourceTypes.RUN = "Run";
        SearchResourceTypes.PROJECT = "Project";
        SearchResourceTypes.APPRESULT = "AppResult";
        SearchResourceTypes.SAMPLE_FILE = "Sample_File";
        SearchResourceTypes.APPRESULT_FILE = "AppResult_File";
        SearchResourceTypes.APPSESSION = "AppSession";
        SearchResourceTypes.SAMPLE = "Sample";
        return SearchResourceTypes;
    })();
    bs.SearchResourceTypes = SearchResourceTypes;
})(bs || (bs = {}));

var bs;
(function (bs) {
    var Util = (function () {
        function Util() {
        }
        // *****************
        Util.escapeHtml = function (str) {
            var entityMap = {
                "&": "&amp;",
                "<": "&lt;",
                ">": "&gt;",
                '"': '&quot;',
                "'": '&#39;',
                "/": '&#x2F;'
            };

            return String(str).replace(/[&<>"'\/]/g, function (s) {
                return entityMap[s];
            });
        };

        // *****************
        Util.getSizeInGigabases = function (sample) {
            return ((parseInt(sample.Read1) + parseInt(sample.Read2)) * parseInt(sample.NumReadsPF)) / 1000000000;
        };

        // ******************
        Util.bytesToSize = function (bytes, precision) {
            var kilobyte = 1024;
            var megabyte = kilobyte * 1024;
            var gigabyte = megabyte * 1024;
            var terabyte = gigabyte * 1024;

            if ((bytes >= 0) && (bytes < kilobyte)) {
                return bytes + ' B';
            } else if ((bytes >= kilobyte) && (bytes < megabyte)) {
                return (bytes / kilobyte).toFixed(precision) + ' KB';
            } else if ((bytes >= megabyte) && (bytes < gigabyte)) {
                return (bytes / megabyte).toFixed(precision) + ' MB';
            } else if ((bytes >= gigabyte) && (bytes < terabyte)) {
                return (bytes / gigabyte).toFixed(precision) + ' GB';
            } else if (bytes >= terabyte) {
                return (bytes / terabyte).toFixed(precision) + ' TB';
            } else {
                return bytes + ' B';
            }
        };

        Util.stringToBytes = function (str) {
            str = str.toLowerCase();
            var bytes = 1;
            var kilobyte = 1024;
            var megabyte = kilobyte * 1024;
            var gigabyte = megabyte * 1024;
            var terabyte = gigabyte * 1024;

            var multiplier = bytes;

            // if the string is less than three characters long,
            // return it.
            if (str.length < 3)
                return parseInt(str);
            else {
                var lastTwoChars = str.substring(str.length - 2, str.length);
            }

            if (lastTwoChars == "kb")
                multiplier = kilobyte;
            else if (lastTwoChars == "mb")
                multiplier = megabyte;
            else if (lastTwoChars == "gb")
                multiplier = gigabyte;

            var num = parseFloat(str.replace(lastTwoChars, "").trim());

            return num * multiplier;
        };

        // ***************************
        Util.getContent = function (ajaxOptions) {
            if ($.support.cors) {
                $.ajax(ajaxOptions);
            } else {
                var xdr = new XDomainRequest();

                // xdr.onprogress = progres;
                // xdr.timeout = tbTO.value;
                xdr.onerror = ajaxOptions.error;
                xdr.onload = ajaxOptions.success;
                xdr.open("get", ajaxOptions.url);
                xdr.send();
            }
        };

        // ***************************
        Util.getFileNameExtension = function (filename) {
            var ext = filename.split('.').pop();
            if (ext)
                return ext.toLowerCase();
            else
                return "";
        };

        Util.getFileCategoryByExtension = function (ext) {
            if (Util.TEXT_FILE_EXTENSIONS.indexOf(ext) > -1) {
                return Util.TEXT_FILE;
            } else if (Util.IMAGE_FILE_EXTENSIONS.indexOf(ext) > -1) {
                return Util.IMAGE_FILE;
            } else if (Util.PDF_FILE_EXTENSIONS.indexOf(ext) > -1) {
                return Util.PDF_FILE;
            } else {
                return Util.OTHER_FILE;
            }
        };

        // ***************************
        // Note: deprecated. Use the other date functions (below) which use moment.js.
        Util.getFriendlyDate = function (unfriendlyDate) {
            var monthNames = [];
            monthNames["01"] = "January";
            monthNames["02"] = "February";
            monthNames["03"] = "March";
            monthNames["04"] = "April";
            monthNames["05"] = "May";
            monthNames["06"] = "June";
            monthNames["07"] = "July";
            monthNames["08"] = "August";
            monthNames["09"] = "September";
            monthNames["10"] = "October";
            monthNames["11"] = "November";
            monthNames["12"] = "December";

            var yyyymmdd = unfriendlyDate.substr(0, 10);
            var yyyy = unfriendlyDate.substr(0, 4);
            var mm = unfriendlyDate.substr(5, 2);
            var dd = unfriendlyDate.substr(8, 2);
            var date = new Date(yyyymmdd);
            var friendlyMonth = monthNames[mm];
            var friendlyDate = friendlyMonth + " " + dd + ", " + yyyy;

            return friendlyDate;
        };

        // ***************************
        Util.getDateConcise = function (unfriendlyUTCDate, toLocal, twitterTimeToday) {
            if (typeof toLocal === "undefined") { toLocal = false; }
            if (typeof twitterTimeToday === "undefined") { twitterTimeToday = true; }
            if (!unfriendlyUTCDate)
                return "";

            if (twitterTimeToday) {
                // check if the date is in the past, relative to the current moment.
                //      - if so, check if the date occurred within the last 24 hours.
                //          - if so, get the duration and append "ago"
                var now = moment.utc();
                var diff = now.diff(moment.utc(unfriendlyUTCDate));

                if (diff > 0 && diff < 86400000) {
                    var duration = moment.duration(diff);
                    var time = {
                        hour: duration.hours(),
                        minute: duration.minutes(),
                        second: duration.seconds()
                    };
                    var timeStr = " ago";
                    if (time.hour > 0) {
                        timeStr = time.hour + " hour" + (time.hour > 1 ? "s" : "") + timeStr;
                    } else if (time.minute > 0) {
                        timeStr = time.minute + " minute" + (time.minute > 1 ? "s" : "") + timeStr;
                    } else if (time.second > 0) {
                        timeStr = time.second + " second" + (time.second > 1 ? "s" : "") + timeStr;
                    }
                    return timeStr;
                }
            }

            if (toLocal)
                return moment.utc(unfriendlyUTCDate).local().format(Util.FORMAT_FRIENDLY_DATE_CONCISE);
            else
                return moment.utc(unfriendlyUTCDate).format(Util.FORMAT_FRIENDLY_DATE_CONCISE);
        };

        // ***************************
        Util.getDatetimeVerbose = function (unfriendlyUTCDate, toLocal) {
            if (typeof toLocal === "undefined") { toLocal = false; }
            if (!unfriendlyUTCDate)
                return "";

            if (toLocal)
                return moment.utc(unfriendlyUTCDate).local().format(Util.FORMAT_FRIENDLY_DATETIME_VERBOSE);
            else
                return moment.utc(unfriendlyUTCDate).format(Util.FORMAT_FRIENDLY_DATETIME_VERBOSE);
        };

        // ***************************
        Util.getHumanizedDuration = function (startDateUnformatted, endDateUnformatted) {
            if (!startDateUnformatted || !endDateUnformatted)
                return "";

            var diff = moment(endDateUnformatted).diff(moment(startDateUnformatted));
            return moment.duration(diff).humanize();
        };

        // ***************************
        Util.getDuration = function (startDateUnformatted, endDateUnformatted, smallestUnit) {
            if (!startDateUnformatted || !endDateUnformatted)
                return "";

            var diff = moment.utc(endDateUnformatted).diff(moment.utc(startDateUnformatted));
            var duration = moment.duration(diff);

            var time = {
                year: duration.years(),
                month: duration.months(),
                day: duration.days(),
                hour: duration.hours(),
                minute: duration.minutes(),
                second: duration.seconds()
            };

            var durationString = "";

            for (var key in time) {
                var value = time[key];

                if (value > 0 || key === smallestUnit) {
                    durationString = durationString + value + " " + key + (value != 1 ? "s" : "") + " ";
                }
                if (key === smallestUnit)
                    break;
            }
            durationString = durationString.trim();
            return durationString;
        };

        // ***************************
        Util.getDateConciseHtml = function (unfriendlyUTCDate, toLocal, twitterTimeToday) {
            if (typeof toLocal === "undefined") { toLocal = false; }
            if (typeof twitterTimeToday === "undefined") { twitterTimeToday = true; }
            return "<time datetime=\"" + unfriendlyUTCDate + "\" rel=\"tooltip\" title=\"" + bs.Util.getDatetimeVerbose(unfriendlyUTCDate, toLocal) + (toLocal ? " (Local timezone)" : " (UTC)") + "\">" + bs.Util.getDateConcise(unfriendlyUTCDate, toLocal, twitterTimeToday) + "</time>";
        };

        // ***************************
        Util.getDatetimeVerboseHtml = function (unfriendlyUTCDate, toLocal) {
            if (typeof toLocal === "undefined") { toLocal = false; }
            return "<time datetime=\"" + unfriendlyUTCDate + "\" rel=\"tooltip\" title=\"" + (toLocal ? "Local timezone" : "UTC") + "\">" + bs.Util.getDatetimeVerbose(unfriendlyUTCDate, toLocal) + "</time>";
        };

        // ****************************
        Util.bindEventsOnce = function (eventName, selector, fn) {
            console.log("Binding events once.");
            var eventBoundClass = eventName + '-event-bound';
            $(document).on(eventName, selector + ':not(.' + eventBoundClass + ')', fn);
            $(selector).addClass(eventBoundClass);
        };

        // ****************************
        Util.uniquifySelector = function (selector) {
            if ($(selector).length < 1)
                return selector;
            else {
                var i = 1;
                do {
                    var newSelector = selector + '-' + (i++);
                } while($(newSelector).length > 0);

                return newSelector;
            }
        };

        Util.formatNumberWithCommas = function (x) {
            if (x.toString().indexOf('.') < 0)
                return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
            return x.toString();
        };

        Util.trimAndAddUnits = function (str) {
            var realStr = str.toString();
            if (this.endsWith(realStr, "000000000"))
                return realStr.replace(new RegExp('000000000$'), ' G');
            if (this.endsWith(realStr, "000000"))
                return realStr.replace(new RegExp('000000$'), ' M');
            if (this.endsWith(realStr, "000"))
                return realStr.replace(new RegExp('000$'), ' K');
            return realStr;
        };

        Util.trimAndAddUnitsV2 = function (str) {
            var realStr = str.toString().replace(/(\d)(?=(\d\d\d)+(?!\d))/g, "$1,").replace(/,\d\d\d/g, ",000");
            if (this.endsWith(realStr, ",000,000,000"))
                return realStr.replace(new RegExp(',000,000,000$'), ' G');
            if (this.endsWith(realStr, ",000,000"))
                return realStr.replace(new RegExp(',000,000$'), ' M');
            if (this.endsWith(realStr, ",000"))
                return realStr.replace(new RegExp(',000$'), ' K');
            return realStr;
        };

        Util.formatThreeDigitPercision = function (str) {
            var realStr = str.toString();
            if (realStr.indexOf('.') < 0) {
                //return trimAndAddUnitsV2(str);
            } else {
                return parseFloat(realStr).toPrecision(3).replace(/^0./, ".");
            }
            return str;
        };

        Util.endsWith = function (str, key) {
            var realStr = str.toString();
            return realStr.indexOf(key, realStr.length - key.length) !== -1;
        };
        Util.IMAGE_FILE = "image";
        Util.TEXT_FILE = "text";
        Util.OTHER_FILE = "other";
        Util.PDF_FILE = "pdf";

        Util.FORMAT_FRIENDLY_DATETIME_VERBOSE = "dddd, MMMM Do YYYY, h:mm:ss a";
        Util.FORMAT_FRIENDLY_DATE_CONCISE = "MMM DD YYYY";

        Util.TEXT_FILE_EXTENSIONS = [
            "txt",
            "asc",
            "bed",
            "gff",
            "gff3",
            "fastq",
            "vcf",
            "csv",
            "xml",
            "tab",
            "tsv",
            "html",
            "xml",
            "xhtml",
            "log",
            "json",
            "sql",
            "yaml"
        ];

        Util.IMAGE_FILE_EXTENSIONS = [
            "jpg",
            "jpeg",
            "gif",
            "png"
        ];

        Util.PDF_FILE_EXTENSIONS = ["pdf"];
        return Util;
    })();
    bs.Util = Util;
})(bs || (bs = {}));

var bs;
(function (bs) {
    (function (MathUtils) {
        // see http://rosettacode.org/wiki/Averages/Simple_moving_average for description of behavior
        var MovingAverager = (function () {
            function MovingAverager(period) {
                this.nums = [];
                this.period = period;
            }
            MovingAverager.prototype.add = function (num) {
                this.nums.push(num);
                if (this.nums.length > this.period) {
                    this.nums.splice(0, 1); // remove the first element of the array
                }
                var sum = 0;
                for (var i = 0; i < this.nums.length; i++) {
                    sum += this.nums[i];
                }
                var divisor = (this.nums.length < this.period) ? this.nums.length : this.period;
                return (sum / divisor);
            };
            return MovingAverager;
        })();
        MathUtils.MovingAverager = MovingAverager;
    })(bs.MathUtils || (bs.MathUtils = {}));
    var MathUtils = bs.MathUtils;
})(bs || (bs = {}));

// **********************************************
// IE8 Shims.
// Here's the Trim Shim.
if (typeof String.prototype.trim !== 'function') {
    String.prototype.trim = function () {
        return this.replace(/^\s+|\s+$/g, '');
    };
}
