var __extends = this.__extends || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    __.prototype = b.prototype;
    d.prototype = new __();
};
var bs;
(function (bs) {
    var HistogramChart = (function (_super) {
        __extends(HistogramChart, _super);
        function HistogramChart(responseData, customConfig) {
            var _this = this;
            _super.call(this, responseData, $.extend({}, {
                scaleXMaxPadding: 1
            }, customConfig));

            if (this.config.contentType == '') {
                this.config.contentType = 'qscore';
            }

            var labelIndex = this.responseData.Series[0].Fields.label;

            if (!this.config.xAxis.tickFormat && this.responseData.Series[0].Data.length > 0 && this.responseData.Series[0].Data[0][labelIndex]) {
                this.config.xAxis.tickFormat = function (d) {
                    return _this.formatXAxis(d);
                };
            }

            this.gradientPfId = bs.Util.uniquifySelector('#gradient-passing-filter');
            this.gradientBelowId = bs.Util.uniquifySelector('#gradient-below-filter');

            this.addGradients(this.config.color);

            this.update(responseData);
        }
        HistogramChart.prototype.update = function (updatedData, updatedSettings) {
            var _this = this;
            _super.prototype.update.call(this, updatedData, updatedSettings);
            if (this.legend == null)
                this.legend = d3.select(this.config.chartContainerSelector).append('div').attr('class', 'legend');
            this.dataLookupIndex = this.responseData.Series[0].Fields;

            var histogram = this.svg.selectAll('g.histogram').data([this.responseData]);
            histogram.enter().append('g').attr('class', 'histogram');

            var series = histogram.selectAll('g.series').data(function (d, i) {
                return d.Series;
            });
            series.enter().append('g').attr('class', 'series');

            var setTitle = function (d) {
                var title = "";

                title += d[_this.dataLookupIndex.label] || (_this.responseData.Axes.X.Title + ": " + d[_this.dataLookupIndex.x]);

                title += "\n" + _this.responseData.Axes.Y.Title + ": " + bs.Util.formatNumberWithCommas(d[_this.dataLookupIndex.height]);

                return title;
            };

            var rect = series.selectAll('rect.bar').data(function (d) {
                var areBarsCentered = false;
                if (d.Options)
                    areBarsCentered = (d.Options.indexOf("Centered") > -1) ? true : false;
                for (var j = 0; j < d.Data.length; j++)
                    d.Data[j].push(areBarsCentered);
                return d.Data;
            });

            var setRectFill = function (d) {
                if (_this.config.contentType == 'qscore')
                    return "url(" + (d[_this.dataLookupIndex.x] >= HistogramChart.PASSING_FILTER_THRESHOLD ? _this.gradientPfId : _this.gradientBelowId) + ")";
                else {
                    return "url(" + _this.gradientPfId + ")";
                }
            };

            rect.enter().append('rect').attr('class', 'bar').attr('title', setTitle).attr('fill', setRectFill).attr('x', function (d) {
                var areBarsCentered = d[d.length - 1];
                return _this.scaleX(d[_this.dataLookupIndex.x]) - (areBarsCentered ? _this.scaleBarWidth(d[_this.dataLookupIndex.width]) / 2 : 0);
            }).attr('width', function (d) {
                return _this.scaleBarWidth(d[_this.dataLookupIndex.width]);
            }).attr('y', function (d) {
                return d3.max(_this.scaleY.range()) - _this.scaleBarHeight(d[_this.dataLookupIndex.height]);
            }).attr('height', function (d) {
                return _this.scaleBarHeight(d[_this.dataLookupIndex.height]);
            });

            rect.attr('title', setTitle).attr('x', function (d) {
                var areBarsCentered = d[d.length - 1];
                return _this.scaleX(d[_this.dataLookupIndex.x]) - (areBarsCentered ? _this.scaleBarWidth(d[_this.dataLookupIndex.width]) / 2 : 0);
            }).attr('width', function (d) {
                return _this.scaleBarWidth(d[_this.dataLookupIndex.width]);
            }).transition().duration(this.config.transitionSpeed).attr('y', function (d) {
                return d3.max(_this.scaleY.range()) - _this.scaleBarHeight(d[_this.dataLookupIndex.height]);
            }).attr('height', function (d) {
                return _this.scaleBarHeight(d[_this.dataLookupIndex.height]);
            }).attr('fill', setRectFill);

            rect.exit().remove();

            if (this.config.labelBarValues) {
                var label = series.selectAll('text.bar-label').data(function (d) {
                    var areBarsCentered = false;
                    if (d.Options)
                        areBarsCentered = (d.Options.indexOf("Centered") > -1) ? true : false;
                    for (var j = 0; j < d.Data.length; j++)
                        d.Data[j].push(areBarsCentered);
                    return d.Data;
                });

                var postFix = this.config.labelBarPostfix || '';

                label.enter().append('text').attr('class', 'bar-label').attr('x', function (d) {
                    var areBarsCentered = d[d.length - 1];
                    return _this.scaleX(d[_this.dataLookupIndex.x]) + (!areBarsCentered ? _this.scaleBarWidth(d[_this.dataLookupIndex.width]) / 2 : 0);
                }).attr('y', function (d) {
                    return d3.max(_this.scaleY.range()) - _this.scaleBarHeight(d[_this.dataLookupIndex.height]) - 3;
                }).text(function (d) {
                    return bs.Util.formatThreeDigitPercision(d[_this.dataLookupIndex.height]) + postFix;
                }).style("text-anchor", "middle");
                ;

                label.attr('x', function (d) {
                    var areBarsCentered = d[d.length - 1];
                    return _this.scaleX(d[_this.dataLookupIndex.x]) + (!areBarsCentered ? _this.scaleBarWidth(d[_this.dataLookupIndex.width]) / 2 : 0);
                }).attr('y', function (d) {
                    return d3.max(_this.scaleY.range()) - _this.scaleBarHeight(d[_this.dataLookupIndex.height]) - 3;
                }).text(function (d) {
                    return bs.Util.formatThreeDigitPercision(d[_this.dataLookupIndex.height]) + postFix;
                }).style("text-anchor", "middle").transition().duration(this.config.transitionSpeed);

                label.exit().remove();
            }

            if (this.config.contentType == 'qscore')
                this.updateLegend();
            else if (this.legend) {
                this.legend.remove();
            }
        };

        HistogramChart.prototype.addGradients = function (color) {
            if (typeof color === "undefined") { color = '#98cc78'; }
            var gradientPassingFilter = this.svg.append('defs').append('linearGradient').attr('id', this.gradientPfId.replace('#', '')).attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

            gradientPassingFilter.append('stop').attr('offset', '0%').attr('style', 'stop-color: ' + color + '; stop-opacity: .5;');

            gradientPassingFilter.append('stop').attr('offset', '50%').attr('style', 'stop-color: ' + color + '; stop-opacity: 1');

            var gradientBelowFilter = this.svg.append('defs').append('linearGradient').attr('id', this.gradientBelowId.replace('#', '')).attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

            gradientBelowFilter.append('stop').attr('offset', '0%').attr('style', 'stop-color: #2c53a1; stop-opacity: .5;');

            gradientBelowFilter.append('stop').attr('offset', '50%').attr('style', 'stop-color: #2c53a1; stop-opacity: 1');
        };

        HistogramChart.prototype.updateLegend = function () {
            var _this = this;
            var legendList = this.legend.selectAll('ul').data([this.responseData]);
            legendList.enter().append('ul').append('li').text('% >= Q30').style('color', 'green').style('font-weight', 'bold');
            legendList.exit().remove();
            var sumAboveQ30 = 0;
            var sumAll = 0;

            var liSize = legendList.selectAll('li.size').data(function (d) {
                for (var i = 0; i < _this.responseData.Series.length; i++) {
                    var series = _this.responseData.Series[i];
                    for (var k = 0; k < series.Data.length; k++) {
                        var height = series.Data[k][_this.dataLookupIndex.height];
                        var x = series.Data[k][_this.dataLookupIndex.x];
                        sumAll += height;

                        if (x >= 30) {
                            sumAboveQ30 += height;
                        }
                    }
                }
                var sumMultiplied = parseFloat((sumAboveQ30 / 1000).toFixed(1));
                return [sumMultiplied + "G"];
            });

            liSize.enter().append('li').attr('class', 'size').text(function (size) {
                return size;
            });
            liSize.exit().remove();
            liSize.attr('class', 'size').text(function (size) {
                return size;
            });

            var percentAboveQ30 = sumAboveQ30 / sumAll;

            var liSumAboveQ30 = legendList.selectAll('li.sum').data([(percentAboveQ30 * 100).toFixed(1) + '%']);

            liSumAboveQ30.enter().append('li').attr('class', 'sum').text(function (sum) {
                return sum;
            });
            liSumAboveQ30.exit().remove();
            liSumAboveQ30.attr('class', 'sum').text(function (sum) {
                return sum;
            });
        };

        HistogramChart.prototype.formatXAxis = function (d) {
            var labelIndex = this.responseData.Series[0].Fields.label;
            var data = this.responseData.Series[0].Data[d];
            return data ? this.trim(data[labelIndex], 10) : null;
        };

        HistogramChart.prototype.trim = function (str, cnt) {
            str = str || '';
            return str.substring(0, Math.min(cnt, str.length));
        };

        HistogramChart.prototype.updateScales = function () {
            var that = this;

            this.scaleX = d3.scale.linear();
            this.scaleX.domain([this.responseData.Axes.X.Min + this.config.scaleXMinPadding, this.responseData.Axes.X.Max + this.config.scaleXMaxPadding]);
            this.scaleX.range([that.config.margins.left, that.config.width - that.config.margins.right]);

            this.scaleY = d3.scale.linear();
            this.scaleY.domain([this.responseData.Axes.Y.Min, this.responseData.Axes.Y.Max]);
            this.scaleY.range([that.config.height - that.config.margins.bottom, that.config.margins.top]);

            this.scaleBarWidth = d3.scale.linear();
            this.scaleBarWidth.domain([0, Math.abs((this.responseData.Axes.X.Max + this.config.scaleXMaxPadding) - (this.responseData.Axes.X.Min + this.config.scaleXMinPadding))]);
            this.scaleBarWidth.range([0, that.config.width - (that.config.margins.left + that.config.margins.right)]);

            this.scaleBarHeight = d3.scale.linear();
            this.scaleBarHeight.domain([this.responseData.Axes.Y.Min, this.responseData.Axes.Y.Max]);
            this.scaleBarHeight.range([0, that.config.height - (that.config.margins.top + that.config.margins.bottom)]);
        };
        HistogramChart.CHART_NAME = "Histogram";

        HistogramChart.PASSING_FILTER_THRESHOLD = 30;
        return HistogramChart;
    })(bs.Chart);
    bs.HistogramChart = HistogramChart;
})(bs || (bs = {}));
