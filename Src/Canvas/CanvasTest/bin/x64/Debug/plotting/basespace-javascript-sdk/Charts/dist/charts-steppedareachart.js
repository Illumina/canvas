var __extends = this.__extends || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    __.prototype = b.prototype;
    d.prototype = new __();
};
var bs;
(function (bs) {
    var SteppedAreaChart = (function (_super) {
        __extends(SteppedAreaChart, _super);
        function SteppedAreaChart(responseData, customConfig) {
            _super.call(this, responseData, customConfig);
            this.groupLabel2Number = {};
            this.groupDisplayLabels = [];
            this.groupLabelPixelsPerChar = 10;
            this.xTickValues = [];

            this.update(responseData, customConfig);
        }
        SteppedAreaChart.prototype.update = function (updatedData, updatedSettings) {
            this.updateWidth(updatedData, updatedSettings);

            _super.prototype.update.call(this, updatedData, updatedSettings);

            var steppedAreaChart = this.svg.selectAll('g.steppedareachart').data([this.responseData]);
            steppedAreaChart.enter().append('g').attr('class', 'steppedareachart');
            steppedAreaChart.exit().remove();

            this.updateSteppedAreaSeries(steppedAreaChart, this.responseData.Series);

            this.updateTitle();
            this.updateLegend();

            this.callTooltips("[title]");

            this.bindMouseEvents();
        };

        SteppedAreaChart.prototype.updateSteppedAreaSeries = function (steppedAreaChart, steppedAreaSeries) {
            var _this = this;
            var lookupIndex;
            var color;

            var series = steppedAreaChart.selectAll('g.series.steppedarea').data(steppedAreaSeries);
            series.enter().append('g');
            series.exit().remove();
            series.datum(function (d) {
                lookupIndex = d.Fields;
                color = d.Color;
                return d;
            }).attr('class', 'series steppedarea');

            var rect = series.selectAll('rect').data(function (d) {
                return d.Data;
            });
            rect.enter().append('rect');
            rect.exit().remove();
            rect.attr("data-powertip", function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "powertip" in parentData.Accessors) {
                    return parentData.Accessors.powertip(d);
                } else {
                    return d.Title;
                }
            }).attr('title', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                return "(" + d[lookupIndex.tooltipTitle] + ", " + parentData.Title + "): " + d[lookupIndex.y];
            }).transition().duration(this.config.transitionSpeed).attr('width', this.rectWidth).attr('height', function (d) {
                return Math.abs(_this.scaleY(d[lookupIndex.y]) - _this.scaleY(0));
            }).attr('x', function (d) {
                return _this.scaleX(_this.groupLabel2Number[d[lookupIndex.x]]) - _this.rectWidth / 2;
            }).attr('y', function (d) {
                return _this.scaleY(Math.max(0, d[lookupIndex.y]));
            });
            rect.style('fill', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "fill" in parentData.Accessors) {
                    return parentData.Accessors.fill(d);
                } else {
                    return parentData.Color.toLowerCase();
                }
            }).style('fill-opacity', SteppedAreaChart.DEFAULT_BAR_FILL_OPACITY).style('stroke', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "fill" in parentData.Accessors) {
                    return parentData.Accessors.fill(d);
                } else {
                    return parentData.Color.toLowerCase();
                }
            }).style('stroke-width', SteppedAreaChart.DEFAULT_BAR_STROKE_WIDTH).on("mouseenter", function (d) {
                $(this).css("stroke-width", SteppedAreaChart.DEFAULT_BAR_HOVER_STROKE_WIDTH).css("fill-opacity", SteppedAreaChart.DEFAULT_BAR_HOVER_FILL_OPACITY);
            }).on("mouseleave", function (d) {
                $(this).css("stroke-width", SteppedAreaChart.DEFAULT_BAR_STROKE_WIDTH).css("fill-opacity", SteppedAreaChart.DEFAULT_BAR_FILL_OPACITY);
            });
        };

        SteppedAreaChart.prototype.updateWidth = function (updatedData, updatedSettings) {
            if (updatedData != null)
                this.responseData = updatedData;
            if (updatedSettings != null)
                $.extend(true, this.config, updatedSettings);
            this.updateTitle();
            this.updateLegend();
            var plotWidth = 0;
            if (this.config.stepSize > 0) {
                plotWidth = (2 + this.responseData.Axes.X.GroupLabels.length) * this.config.stepSize;
            }
            var minWidth = this.config.margins.left + this.config.margins.right + (this.config.fontSize == null ? SteppedAreaChart.DEFAULT_FONT_SIZE : this.config.fontSize) + Math.max(this.titleWidth, this.legendWidth, plotWidth);

            var width = Math.max(minWidth, this.config.width);
            if (updatedSettings != null) {
                updatedSettings.width = width;
            }
            $(this.config.chartContainerSelector).width(width);
        };

        SteppedAreaChart.prototype.updateAxes = function () {
            var _this = this;
            this.svg.selectAll('.axis-x, .axis-y').remove();

            var fontSize = this.config.fontSize != null ? this.config.fontSize : SteppedAreaChart.DEFAULT_FONT_SIZE;

            this.yAxis = d3.svg.axis();
            this.yAxis.scale(this.scaleY);
            this.yAxis.tickSize(-(this.config.width - ((this.config.margins.left) + this.config.margins.right)), 3, 0);

            this.configureAxis(this.yAxis, this.config.yAxis);

            this.svg.append('g').attr('class', 'axis-y').attr('transform', 'translate(' + this.config.margins.left + ', 0)').call(this.yAxis).selectAll('text').style('font-size', fontSize + 'px');

            this.groupLabel2Number = {};
            this.groupDisplayLabels = [];
            this.maxGroupDisplayLabelLength = Math.floor((this.config.margins.bottom - 20) / this.groupLabelPixelsPerChar);
            this.xTickValues = [];
            for (var i = 0; i < this.responseData.Axes.X.GroupLabels.length; i++) {
                var label = this.responseData.Axes.X.GroupLabels[i];
                var displayLabel = label;
                if (displayLabel.length > this.maxGroupDisplayLabelLength) {
                    displayLabel = displayLabel.substr(0, this.maxGroupDisplayLabelLength - 3) + "...";
                }
                this.groupDisplayLabels.push(displayLabel);
                this.groupLabel2Number[label] = i + 1;
                this.xTickValues.push(i + 1);
            }
            this.xAxis = d3.svg.axis();
            this.xAxis.scale(this.scaleX);
            this.xAxis.orient('bottom');
            this.xAxis.tickSize(0, 3, 0);
            this.xAxis.tickValues(this.xTickValues);
            this.xAxis.tickFormat(function (v) {
                return _this.groupDisplayLabels[v - 1];
            });

            this.svg.append('g').attr('class', 'axis-x').attr('transform', 'translate(0, ' + (this.config.height - this.config.margins.bottom) + ')').call(this.xAxis).selectAll('text').style('font-size', fontSize + 'px').style("text-anchor", "end").attr("dx", "-5px").attr("dy", "-1px").attr("transform", "rotate(-90)").append("title").text(function (d, i) {
                return _this.responseData.Axes.X.GroupLabels[i];
            });

            var nPixels = 0, nChars = 0;
            var texts = this.svg.selectAll('g.axis-x text')[0];
            for (var i = 0; i < texts.length; i++) {
                nChars += this.groupDisplayLabels[i].length;
                nPixels += texts[i].getBBox().width;
            }
            this.groupLabelPixelsPerChar = Math.ceil(nPixels / nChars);

            this.svg.select('.axis-y').append('text').style('font-size', fontSize + 'px').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'rotate(-90 0 0), translate(' + (-(this.config.height - this.config.margins.top - this.config.margins.bottom) / 2 - this.config.margins.top) + ', ' + (-this.config.margins.left + 20) + ')').text(this.responseData.Axes.Y.Title);

            this.svg.select('.axis-x').append('text').style('font-size', fontSize + 'px').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'translate(' + ((this.config.width - this.config.margins.left - this.config.margins.right) / 2 + this.config.margins.left) + ', ' + (this.config.margins.bottom - 20) + ')').text(this.responseData.Axes.X.Title);
        };

        SteppedAreaChart.prototype.updateTitle = function () {
            if (this.svgTitle == null) {
                this.svgTitle = this.svg.append('g').attr('class', 'title');
            }
            this.svgTitle.selectAll('text').remove();
            var fontSize = this.config.fontSize != null ? this.config.fontSize : SteppedAreaChart.DEFAULT_FONT_SIZE;
            this.svgTitle.attr('transform', 'translate(' + (this.config.margins.left + fontSize) + ', ' + (2 * fontSize) + ')');
            var text = this.svgTitle.append("text").attr("x", 0).attr("y", 0).attr('alignment-baseline', 'middle').style('font-size', fontSize + 'px').text(this.responseData.Title);
            this.titleWidth = text[0][0].getBBox().width;
        };

        SteppedAreaChart.prototype.updateLegend = function () {
            if (this.svgLegend == null) {
                this.svgLegend = this.svg.append('g').attr('class', 'legend');
            }

            if (!this.config.legend) {
                this.svgLegend.style('display', 'none');
            } else {
                this.svgLegend.style('display', 'block');
            }

            this.svgLegend.selectAll('text').remove();
            this.svgLegend.selectAll('rect').remove();
            var fontSize = this.config.fontSize != null ? this.config.fontSize : SteppedAreaChart.DEFAULT_FONT_SIZE;
            this.svgLegend.attr('transform', 'translate(' + (this.config.margins.left + fontSize) + ', ' + (4 * fontSize) + ')');

            for (var legendXOffset = 0, legendYOffset = 0, i = 0; i < this.responseData.Series.length; i++) {
                var series = this.responseData.Series[i];
                this.svgLegend.append("rect").attr("class", "legend-rect").attr("x", legendXOffset).attr("y", legendYOffset - fontSize).attr("width", fontSize).attr("height", fontSize).style("fill", series.Color.toLowerCase()).style("fill-opacity", SteppedAreaChart.DEFAULT_BAR_FILL_OPACITY).style("stroke", series.Color.toLowerCase()).style("stroke-width", SteppedAreaChart.DEFAULT_BAR_STROKE_WIDTH);
                legendXOffset += (fontSize + 5);
                var text = this.svgLegend.append("text").attr("x", legendXOffset).attr("y", legendYOffset).attr('alignment-baseline', 'middle').style('font-size', fontSize + 'px').style('fill', series.Color.toLowerCase()).text(series.Title);
                legendXOffset += (text[0][0].getBBox().width + fontSize);
            }
            this.legendWidth = legendXOffset;
        };

        SteppedAreaChart.prototype.addLegend = function () {
        };

        SteppedAreaChart.prototype.updateScales = function () {
            this.scaleX = d3.scale.linear();
            this.scaleX.domain([0, this.responseData.Axes.X.GroupLabels.length + 1]);
            this.scaleX.range([this.config.margins.left, this.config.width - this.config.margins.right]);
            this.scaleX.nice();
            this.rectWidth = Math.floor(this.scaleX(0.45) - this.scaleX(0)) * 2;

            var yMin = 0, yMax = 0;
            for (var i = 0; i < this.responseData.Series.length; i++) {
                for (var j = 0; j < this.responseData.Series[i].Data.length; j++) {
                    var y = this.responseData.Series[i].Data[j][this.responseData.Series[i].Fields.y];
                    yMin = Math.min(yMin, y);
                    yMax = Math.max(yMax, y);
                }
            }
            this.responseData.Axes.Y.Min = Math.min(this.responseData.Axes.Y.Min, yMin);
            this.responseData.Axes.Y.Max = Math.max(this.responseData.Axes.Y.Max, yMax);

            this.scaleY = this.config.scaleY();
            if (this.config.scaleY === d3.scale.log) {
                this.scaleY.domain([Math.max(bs.Chart.MIN_LOG_INPUT, this.responseData.Axes.Y.Min), this.responseData.Axes.Y.Max]);
            } else {
                this.scaleY.domain([this.responseData.Axes.Y.Min, this.responseData.Axes.Y.Max]);
            }
            this.scaleY.range([this.config.height - this.config.margins.bottom, this.config.margins.top]);
            if (this.config.scaleY && this.config.scaleY.nice === true) {
                this.scaleY.nice();
            }
        };

        SteppedAreaChart.prototype.bindMouseEvents = function () {
            var that = this;

            $(this.config.chartContainerSelector + ' .series.steppedarea rect').hover(function () {
                var parentData = d3.select(this.parentNode).datum();
                $(that.config.chartContainerSelector + ' g.legend text').each(function () {
                    if (parentData.Title == $(this).text()) {
                        $(this).css('text-decoration', 'underline');
                    } else {
                        $(this).css('text-decoration', '');
                    }
                });
            }, function () {
                $(that.config.chartContainerSelector + ' g.legend text').css('text-decoration', '');
            });
            $(this.config.chartContainerSelector + ' g.legend text').hover(function () {
                var title = $(this).text();
                $(that.config.chartContainerSelector + ' .series.steppedarea rect').each(function () {
                    var parentData = d3.select(this.parentNode).datum();
                    if (parentData.Title == title) {
                        $(this).css("fill-opacity", SteppedAreaChart.DEFAULT_BAR_HOVER_FILL_OPACITY);
                        $(this).css("stroke-width", SteppedAreaChart.DEFAULT_BAR_HOVER_STROKE_WIDTH);
                    }
                });
            }, function () {
                $(that.config.chartContainerSelector + ' .series.steppedarea rect').each(function () {
                    var parentData = d3.select(this.parentNode).datum();
                    $(this).css("fill-opacity", SteppedAreaChart.DEFAULT_BAR_FILL_OPACITY);
                    $(this).css("stroke-width", SteppedAreaChart.DEFAULT_BAR_STROKE_WIDTH);
                });
            });
        };
        SteppedAreaChart.CHART_NAME = "SteppedAreaChart";

        SteppedAreaChart.DEFAULT_FONT_SIZE = 10;
        SteppedAreaChart.DEFAULT_BAR_FILL_OPACITY = 0.5;
        SteppedAreaChart.DEFAULT_BAR_STROKE_WIDTH = 2;
        SteppedAreaChart.DEFAULT_BAR_HOVER_FILL_OPACITY = 0.8;
        SteppedAreaChart.DEFAULT_BAR_HOVER_STROKE_WIDTH = 3;

        SteppedAreaChart.DEFAULT_FIELD_INDICES = {
            tooltipTitle: 0,
            x: 0,
            y: 1
        };
        return SteppedAreaChart;
    })(bs.Chart);
    bs.SteppedAreaChart = SteppedAreaChart;
})(bs || (bs = {}));
