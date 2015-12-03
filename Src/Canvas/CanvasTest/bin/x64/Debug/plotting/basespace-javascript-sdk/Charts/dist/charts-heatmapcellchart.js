var __extends = this.__extends || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    __.prototype = b.prototype;
    d.prototype = new __();
};
var bs;
(function (bs) {
    var HeatmapCellChart = (function (_super) {
        __extends(HeatmapCellChart, _super);
        function HeatmapCellChart(responseData, customConfig) {
            _super.call(this, responseData, customConfig);
            this.xGroupLabel2Number = {};
            this.xGroupDisplayLabels = [];
            this.xTickValues = [];
            this.yGroupLabel2Number = {};
            this.yGroupDisplayLabels = [];
            this.yTickValues = [];
            this.groupLabelPixelsPerChar = 10;

            var config = {
                valueAsHeight: false,
                colorScale: true
            };
            $.extend(true, config, this.config);
            this.config = config;

            this.update(responseData, customConfig);
        }
        HeatmapCellChart.prototype.update = function (updatedData, updatedSettings) {
            this.updateDefs();

            _super.prototype.update.call(this, updatedData, updatedSettings);

            var heatmapCellChart = this.svg.selectAll('g.heatmapcellchart').data([this.responseData]);
            heatmapCellChart.enter().append('g').attr('class', 'heatmapcellchart');
            heatmapCellChart.exit().remove();

            this.updateHeatmapCellData(heatmapCellChart);

            this.callTooltips("[title]");

            this.bindMouseEvents();
        };

        HeatmapCellChart.prototype.updateHeatmapCellData = function (heatmapCellChart) {
            var _this = this;
            var heatmap = heatmapCellChart.selectAll('g.heatmapcell').data(function (d) {
                return [d];
            });
            heatmap.enter().append('g').attr('class', 'heatmapcell');
            heatmap.exit().remove();

            var rects = heatmap.selectAll('rect').data(function (d) {
                var rectData = [];
                for (var i = 0; i < d.Data.length; i++) {
                    for (var j = 0; j < d.Data[i].length; j++) {
                        rectData.push({
                            value: d.Data[i][j],
                            x: j + 1,
                            y: i + 1
                        });
                    }
                }
                return rectData;
            });
            rects.enter().append('rect');
            rects.exit().remove();

            rects.transition().duration(this.config.transitionSpeed).attr('width', this.rectWidth).attr('height', function (d) {
                if (_this.config.valueAsHeight) {
                    return _this.heightScale(d.value);
                } else {
                    return _this.rectHeight;
                }
            }).attr('x', function (d) {
                return _this.scaleX(d.x) - _this.rectWidth / 2;
            }).attr('y', function (d) {
                if (_this.config.valueAsHeight) {
                    return _this.scaleY(d.y) + _this.rectHeight / 2 - _this.heightScale(d.value);
                } else {
                    return _this.scaleY(d.y) - _this.rectHeight / 2;
                }
            });

            rects.style('fill', function (d) {
                return _this.colorScale(d.value);
            }).style('fill-opacity', HeatmapCellChart.DEFAULT_RECT_FILL_OPACITY).style('stroke', HeatmapCellChart.DEFAULT_RECT_STROKE).style('stroke-width', HeatmapCellChart.DEFAULT_RECT_STROKE_WIDTH).on('mouseenter', function (d) {
                var strokeWidth = $(this).css("stroke-width");
                $(this).css('stroke-width', HeatmapCellChart.DEFAULT_RECT_HOVER_STROKE_WIDTH);
                var that = this;
                setTimeout(function () {
                    $(that).css("stroke-width", strokeWidth);
                }, 200);
            }).on('click', function (d) {
                heatmap.selectAll("rect").style("stroke", HeatmapCellChart.DEFAULT_RECT_STROKE).style("stroke-width", HeatmapCellChart.DEFAULT_RECT_STROKE_WIDTH);

                $(this).css("stroke", HeatmapCellChart.DEFAULT_RECT_CLICK_STROKE).css("stroke-width", HeatmapCellChart.DEFAULT_RECT_CLICK_STROKE_WIDTH);
            });

            rects.append('title').text(function (d) {
                return "(" + _this.responseData.Axes.X.GroupLabels[d.x - 1] + ", " + _this.responseData.Axes.Y.GroupLabels[d.y - 1] + "): " + d.value;
            });
        };

        HeatmapCellChart.prototype.renderXYAxisTickLabels = function () {
            var _this = this;
            this.svg.selectAll('.axis-x, .axis-y').remove();

            var fontSize = this.config.fontSize != null ? this.config.fontSize : HeatmapCellChart.DEFAULT_FONT_SIZE;

            this.yGroupLabel2Number = {};
            this.yGroupDisplayLabels = [];
            this.maxYGroupDisplayLabelLength = Math.floor((this.config.margins.left - 2 * fontSize) / this.groupLabelPixelsPerChar);
            this.yTickValues = [];
            for (var i = 0; i < this.responseData.Axes.Y.GroupLabels.length; i++) {
                var label = this.responseData.Axes.Y.GroupLabels[i];
                var displayLabel = label;
                if (displayLabel.length > this.maxYGroupDisplayLabelLength) {
                    displayLabel = displayLabel.substr(0, this.maxYGroupDisplayLabelLength - 3) + "...";
                }
                this.yGroupDisplayLabels.push(displayLabel);
                this.yGroupLabel2Number[label] = i + 1;
                this.yTickValues.push(i + 1);
            }
            this.yAxis = d3.svg.axis();
            this.yAxis.scale(this.scaleY);
            this.yAxis.orient('left');
            this.yAxis.tickSize(0, 3, 0);
            this.yAxis.tickValues(this.yTickValues);
            this.yAxis.tickFormat(function (v) {
                return _this.yGroupDisplayLabels[v - 1];
            });

            this.svg.append('g').attr('class', 'axis-y').attr('transform', 'translate(' + this.config.margins.left + ', 0)').call(this.yAxis).selectAll('text').style('font-size', fontSize + 'px').style("text-anchor", "end").append("title").text(function (d, i) {
                return _this.responseData.Axes.Y.GroupLabels[i];
            });

            this.xGroupLabel2Number = {};
            this.xGroupDisplayLabels = [];
            this.maxXGroupDisplayLabelLength = Math.floor((this.config.margins.top - 2 * fontSize) / this.groupLabelPixelsPerChar);
            this.xTickValues = [];
            for (var i = 0; i < this.responseData.Axes.X.GroupLabels.length; i++) {
                var label = this.responseData.Axes.X.GroupLabels[i];
                var displayLabel = label;
                if (displayLabel.length > this.maxXGroupDisplayLabelLength) {
                    displayLabel = displayLabel.substr(0, this.maxXGroupDisplayLabelLength - 3) + "...";
                }
                this.xGroupDisplayLabels.push(displayLabel);
                this.xGroupLabel2Number[label] = i + 1;
                this.xTickValues.push(i + 1);
            }
            this.xAxis = d3.svg.axis();
            this.xAxis.scale(this.scaleX);
            this.xAxis.orient('top');
            this.xAxis.tickSize(0, 3, 0);
            this.xAxis.tickValues(this.xTickValues);
            this.xAxis.tickFormat(function (v) {
                return _this.xGroupDisplayLabels[v - 1];
            });

            this.svg.append('g').attr('class', 'axis-x').attr('transform', 'translate(0, ' + (this.config.margins.top) + ')').call(this.xAxis).selectAll('text').style('font-size', fontSize + 'px').style("text-anchor", "start").attr("dy", "7px").attr("transform", "rotate(-90)").append("title").text(function (d, i) {
                return _this.responseData.Axes.X.GroupLabels[i];
            });
        };

        HeatmapCellChart.prototype.calculateGroupLabelPixelPerChar = function () {
            var nPixels = 0, nChars = 0;
            var texts = this.svg.selectAll('g.axis-x text')[0];
            for (var i = 0; i < texts.length; i++) {
                nChars += this.xGroupDisplayLabels[i].length;
                nPixels += texts[i].getBBox().width;
            }
            texts = this.svg.selectAll('g.axis-y text')[0];
            for (var i = 0; i < texts.length; i++) {
                nChars += this.yGroupDisplayLabels[i].length;
                nPixels += texts[i].getBBox().width;
            }
            this.groupLabelPixelsPerChar = Math.ceil(nPixels / nChars);
        };

        HeatmapCellChart.prototype.updateAxes = function () {
            this.renderXYAxisTickLabels();

            this.calculateGroupLabelPixelPerChar();

            this.renderXYAxisTickLabels();

            var fontSize = this.config.fontSize != null ? this.config.fontSize : HeatmapCellChart.DEFAULT_FONT_SIZE;

            this.svg.select('.axis-y').append('text').style('font-size', fontSize + 'px').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'rotate(-90 0 0), translate(' + (-(this.config.height - this.config.margins.top - this.config.margins.bottom) / 2 - this.config.margins.top) + ', ' + (-this.config.margins.left + fontSize * 2) + ')').text(this.responseData.Axes.Y.Title);

            this.svg.select('.axis-x').append('text').style('font-size', fontSize + 'px').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'translate(' + ((this.config.width - this.config.margins.left - this.config.margins.right) / 2 + this.config.margins.left) + ', ' + (-this.config.margins.top + fontSize * 2) + ')').text(this.responseData.Axes.X.Title);

            if (this.config.colorScale) {
                this.svg.selectAll('.heat-axis').remove();
                this.heatAxis = d3.svg.axis();
                this.heatAxis.scale(this.heatScale).orient('right').ticks(7);

                var rightMargin1 = 1.2 * fontSize;
                var rightMargin2 = Math.max(0, this.config.margins.right - rightMargin1);
                var heatAxisGroup = this.svg.append('g').attr('class', 'heat-axis').attr('transform', 'translate(' + (this.config.width - rightMargin2) + ', 0)');

                heatAxisGroup.call(this.heatAxis);
                heatAxisGroup.selectAll('g.tick.major text').style('font-size', Math.max(6, fontSize - 2) + 'px');

                heatAxisGroup.insert('rect', ':first-child').attr('class', 'axis-reference-gradient').attr('x', -fontSize).attr('y', this.config.margins.top + this.heatAxisMargin).attr('width', 1.2 * fontSize).attr('height', this.heatAxisHeight).attr('fill', 'url(#' + this.referenceGradientId + ')');
                heatAxisGroup.append("text").attr("class", "heat-axis-label").attr("transform", "rotate(-90 0 0), translate(" + (-this.config.margins.top - (this.config.height - this.config.margins.top - this.config.margins.bottom) / 2) + ", " + (rightMargin2 - 10) + ")").attr("x", 0).attr("y", 0).style("text-anchor", "middle").style('font-size', fontSize + 'px').text(this.responseData.ColorScaleTitle == null ? "Color Scale" : this.responseData.ColorScaleTitle);

                this.svg.selectAll(".heat-axis path").style("fill", "none");
            }
        };

        HeatmapCellChart.prototype.updateLegend = function () {
        };

        HeatmapCellChart.prototype.addLegend = function () {
        };

        HeatmapCellChart.prototype.updateScales = function () {
            this.scaleX = d3.scale.linear();
            this.scaleX.domain([0, this.responseData.Axes.X.GroupLabels.length + 1]);
            this.scaleX.range([this.config.margins.left, this.config.width - this.config.margins.right]);
            this.rectWidth = Math.floor(this.scaleX(0.45) - this.scaleX(0)) * 2;

            this.scaleY = d3.scale.linear();
            this.scaleY.domain([0, this.responseData.Axes.Y.GroupLabels.length + 1]);
            this.scaleY.range([this.config.margins.top, this.config.height - this.config.margins.bottom]);
            this.rectHeight = Math.floor(this.scaleY(0.45) - this.scaleY(0)) * 2;

            this.colorScale = d3.scale.linear();
            this.colorScale.domain(this.responseData.ColorScaleDomain);
            this.colorScale.range(this.responseData.ColorScaleRange);

            var min = this.responseData.ColorScaleDomain.reduce(function (prev, curr, i, arr) {
                return prev < curr ? prev : curr;
            });
            var max = this.responseData.ColorScaleDomain.reduce(function (prev, curr, i, arr) {
                return prev > curr ? prev : curr;
            });
            this.heatScale = d3.scale.linear();
            this.heatScale.domain([min, max]);
            var fontSize = this.config.fontSize != null ? this.config.fontSize : HeatmapCellChart.DEFAULT_FONT_SIZE;
            this.maxHeatAxisHeight = 20 * fontSize;
            this.heatAxisHeight = Math.min(this.maxHeatAxisHeight, this.config.height - this.config.margins.top - this.config.margins.bottom);
            this.heatAxisMargin = (this.config.height - this.config.margins.top - this.config.margins.bottom - this.heatAxisHeight) / 2;
            this.heatScale.range([this.config.height - this.config.margins.bottom - this.heatAxisMargin, this.config.margins.top + this.heatAxisMargin]);

            if (this.config.valueAsHeight) {
                this.heightScale = d3.scale.linear();
                this.heightScale.domain([min, max]);
                this.heightScale.range([0, this.rectHeight]);
            }
        };

        HeatmapCellChart.prototype.updateDefs = function () {
            var _this = this;
            this.referenceGradientId = this.config.chartContainerSelector.substring(1) + '-reference-gradient';
            if (this.defs == null) {
                this.defs = this.svg.append('defs');
                var referenceGradient = this.defs.append('linearGradient').attr('id', this.referenceGradientId).attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

                var colorScaleDomain = this.responseData.ColorScaleDomain;
                var domainSize = Math.abs(colorScaleDomain[0] - colorScaleDomain[colorScaleDomain.length - 1]);
                colorScaleDomain = colorScaleDomain.map(function (value, i, arr) {
                    return Math.abs(value - arr[0]) / domainSize;
                });
                referenceGradient.selectAll('stop').data(colorScaleDomain).enter().append('stop').attr('offset', function (d, i) {
                    return d * 100 + "%";
                }).attr('style', function (d, i) {
                    return "stop-color: " + _this.responseData.ColorScaleRange[i] + "; stop-opacity:1";
                });
            }
        };

        HeatmapCellChart.prototype.bindMouseEvents = function () {
            var that = this;
        };
        HeatmapCellChart.CHART_NAME = "HeatmapCellChart";

        HeatmapCellChart.DEFAULT_FONT_SIZE = 10;
        HeatmapCellChart.DEFAULT_RECT_FILL_OPACITY = 1;
        HeatmapCellChart.DEFAULT_RECT_STROKE = "grey";
        HeatmapCellChart.DEFAULT_RECT_STROKE_WIDTH = 1;
        HeatmapCellChart.DEFAULT_RECT_HOVER_STROKE_WIDTH = 2;
        HeatmapCellChart.DEFAULT_RECT_CLICK_STROKE = "black";
        HeatmapCellChart.DEFAULT_RECT_CLICK_STROKE_WIDTH = 2;
        return HeatmapCellChart;
    })(bs.Chart);
    bs.HeatmapCellChart = HeatmapCellChart;
})(bs || (bs = {}));
