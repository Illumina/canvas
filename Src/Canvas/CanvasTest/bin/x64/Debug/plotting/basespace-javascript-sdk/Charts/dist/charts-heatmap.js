var __extends = this.__extends || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    __.prototype = b.prototype;
    d.prototype = new __();
};
var bs;
(function (bs) {
    var HeatmapChart = (function (_super) {
        __extends(HeatmapChart, _super);
        function HeatmapChart(responseData, customConfig) {
            _super.call(this, responseData, $.extend(true, {
                margins: { right: 75 }
            }, customConfig));
            this.showDefaultAxes = true;
            this.laneHorizontalSpacing = 7;

            this.update(responseData);
        }
        HeatmapChart.prototype.update = function (updatedData, updatedSettings) {
            var _this = this;
            _super.prototype.update.call(this, updatedData, updatedSettings);

            if (this.defs == null) {
                this.defs = this.svg.append('defs');
                var referenceGradient = this.defs.append('linearGradient').attr('id', 'reference-gradient').attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

                referenceGradient.selectAll('stop').data(this.colorScaleDomain).enter().append('stop').attr('offset', function (d, i) {
                    return d * 100 + "%";
                }).attr('style', function (d, i) {
                    return "stop-color: " + _this.colorScaleRange[i] + "; stop-opacity:1";
                });
            }

            var gradients = this.defs.selectAll('.heatmap-gradient').data(this.responseData.Data);

            gradients.enter().append('linearGradient').attr('class', 'heatmap-gradient').attr('id', function (d, i) {
                return 'heatmap-gradient-' + i;
            }).attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

            gradients.attr('class', 'heatmap-gradient').attr('id', function (d, i) {
                return 'heatmap-gradient-' + i;
            }).transition().duration(this.config.transitionSpeed).attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

            gradients.exit().remove();

            var getColorStops = function (data, index) {
                var lastVal = -1;
                var colorStopArr = [];
                for (var i = 0; i < data.length; i++) {
                    var val = data[i];

                    if (val != lastVal) {
                        colorStopArr.push([i, val]);
                    }

                    lastVal = val;
                }
                return colorStopArr;
            };

            var setColorStyle = function (data, index) {
                var color = _this.colorScale(data);
                return "stop-color:" + color + ";stop-opacity: 1;";
            };

            var setColorOffset = function (data, index) {
                var offset = _this.gradientOffsetScale(index);
                return offset + "%";
            };

            var stops = gradients.selectAll('stop').data(getColorStops);

            stops.enter().append('stop').attr('offset', function (d, i) {
                return setColorOffset(d[1], d[0]);
            }).attr('style', function (d, i) {
                return setColorStyle(d[1], d[0]);
            });
            stops.exit().remove();
            stops.attr('offset', function (d, i) {
                return setColorOffset(d[1], d[0]);
            }).attr('style', function (d, i) {
                return setColorStyle(d[1], d[0]);
            });

            var colWidth = (this.config.width - (this.config.margins.left + this.config.margins.right)) / this.responseData.Data.length;
            var colHeight = Math.floor(this.config.height - (this.config.margins.top + this.config.margins.bottom));

            var g = this.svg.selectAll('g.heatmap').data([this.responseData.Data]);
            g.enter().insert('g', ':first-child').attr('class', 'heatmap');
            g.exit().remove();

            var rects = g.selectAll('rect').data(function (d) {
                return d;
            });
            rects.enter().append('rect').attr('fill', function (d, i) {
                return "url(#heatmap-gradient-" + i + ")";
            }).attr('stroke', '0').attr('stroke-width', '0').attr('x', function (d, i) {
                return _this.scaleX(i);
            }).attr('width', colWidth).attr('y', this.config.margins.top).attr('height', colHeight);
            rects.exit().remove();
            rects.attr('fill', function (d, i) {
                return "url(#heatmap-gradient-" + i + ")";
            }).attr('x', function (d, i) {
                return _this.scaleX(i);
            }).attr('width', colWidth).attr('y', this.config.margins.top).attr('height', colHeight);
        };

        HeatmapChart.prototype.updateAxes = function () {
            this.svg.selectAll('.axis-x, .axis-y, .heat-axis').remove();

            this.yAxis = d3.svg.axis();
            this.yAxis.scale(this.scaleY);
            this.yAxis.orient('left');
            this.yAxis.tickSize(-(this.config.width - (this.config.margins.left + this.config.margins.right)), 0, -(this.config.width - (this.config.margins.left + this.config.margins.right)));
            this.yAxis.tickSubdivide(1);
            this.svg.append('g').attr('transform', 'translate(' + (this.config.margins.left) + ', 0)').attr('class', 'axis-y').call(this.yAxis);

            this.xAxis = d3.svg.axis();
            this.xAxis.scale(this.scaleX);
            this.xAxis.orient('bottom');
            this.xAxis.tickSubdivide(1);

            this.xAxis.tickSize(-(this.config.height - (this.config.margins.top + this.config.margins.bottom)), 0, -(this.config.height - (this.config.margins.top + this.config.margins.bottom)));
            this.svg.append('g').attr('transform', 'translate(0, ' + (this.config.height - this.config.margins.bottom) + ')').attr('class', 'axis-x').call(this.xAxis);

            this.heatAxis = d3.svg.axis();
            this.heatAxis.scale(this.heatScale);
            this.heatAxis.orient('right');
            this.heatAxis.ticks(4);

            var heatAxisGroup = this.svg.append('g').attr('class', 'heat-axis').attr('transform', 'translate(' + (this.config.width - (this.config.margins.right * .75)) + ', 0)');

            heatAxisGroup.call(this.heatAxis);

            $('.heat-axis g:last-of-type text').text('% of Max');

            heatAxisGroup.insert('rect', ':first-child').attr('class', 'axis-reference-gradient').attr('x', -10).attr('y', this.config.margins.top).attr('width', 15).attr('height', (this.config.height - (this.config.margins.top + this.config.margins.bottom))).attr('fill', 'url(#reference-gradient)');
            this.svg.selectAll('.axis-x .domain, .axis-y .domain, .heat-axis .domain').remove();
        };

        HeatmapChart.prototype.updateScales = function () {
            this.colorScaleRange = ["white", "green", "yellow", "red"];
            this.colorScaleDomain = [0, 0.33, 0.66, 1];

            this.scaleX = d3.scale.linear();
            this.scaleX.domain([this.responseData.Axes.X.Min, this.responseData.Axes.X.Max]);
            this.scaleX.range([0 + this.config.margins.left, this.config.width - this.config.margins.right]);

            this.scaleY = d3.scale.linear();
            this.scaleY.domain([this.responseData.Axes.Y.Min, this.responseData.Axes.Y.Max]);
            this.scaleY.range([this.config.height - this.config.margins.bottom, this.config.margins.top]);

            this.heatScale = d3.scale.linear();
            this.heatScale.domain([0, 100]);
            this.heatScale.range([this.config.height - this.config.margins.bottom, this.config.margins.top]);

            this.gradientOffsetScale = d3.scale.linear();
            this.gradientOffsetScale.domain([0, this.responseData.Data[0].length]);
            this.gradientOffsetScale.range([0, 100]);

            this.colorScale = d3.scale.linear().domain([0, 100]);
            this.colorScale.domain(this.colorScaleDomain.map(this.colorScale.invert));
            this.colorScale.range(this.colorScaleRange);
        };
        HeatmapChart.CHART_NAME = "Heatmap";
        return HeatmapChart;
    })(bs.Chart);
    bs.HeatmapChart = HeatmapChart;
})(bs || (bs = {}));
