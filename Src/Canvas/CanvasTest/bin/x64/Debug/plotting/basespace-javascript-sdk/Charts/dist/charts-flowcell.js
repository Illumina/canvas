var __extends = this.__extends || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    __.prototype = b.prototype;
    d.prototype = new __();
};
var bs;
(function (bs) {
    var FlowcellChart = (function (_super) {
        __extends(FlowcellChart, _super);
        function FlowcellChart(responseData, customConfig) {
            _super.call(this, responseData, $.extend(true, {
                margins: { right: 60, left: 15, top: 15, bottom: 15 }
            }, customConfig));
            this.showDefaultAxes = false;

            this.update(responseData);
        }
        FlowcellChart.prototype.update = function (updatedData, updatedSettings) {
            var _this = this;
            _super.prototype.update.call(this, updatedData, updatedSettings);

            if (this.defs == null) {
                this.defs = this.svg.append('defs');
                var gradient = this.defs.append('linearGradient').attr('id', 'flowcell-color-scale-gradient').attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

                gradient.selectAll('stop').data(this.colorScaleDomain, function (d) {
                    return d;
                }).enter().append('stop').attr('offset', function (d, i) {
                    var percent = "";
                    percent = (d * 100) + "%";
                    return percent;
                }).attr('style', function (d, i) {
                    var color = _this.colorScaleRange[i];
                    return "stop-color:" + color + ";stop-opacity: 1;";
                });
            }

            var g = this.svg.selectAll('g.flowcell').data([this.responseData]);
            g.enter().append('g').attr('class', 'flowcell');
            g.exit().remove();

            var lanes = g.selectAll('g.lane').data(function (d) {
                return d.TileValuesByLane;
            });
            lanes.enter().append('g').attr('class', 'lane');
            lanes.exit().remove();

            var setRectX = function (laneIndex, index) {
                var positionWithinLane = Math.floor(index / _this.responseData.RowCount) * _this.laneColumnWidth;
                var rectX = (positionWithinLane + ((laneIndex) * _this.laneWidth + (laneIndex * _this.laneHorizontalSpacing)));
                return rectX + _this.config.margins.left;
            };

            var setRectY = function (index) {
                var rectY = (index == 0) ? 0 : (index % _this.responseData.RowCount) * _this.laneHeight;
                return rectY + _this.config.margins.top;
            };

            var setRectFill = function (d) {
                var fillColor;

                if ((d[1] == -1) || !_this.axisExists(_this.responseData.Axis))
                    fillColor = "#666";
                else
                    fillColor = _this.colorScale(d[1]);

                return fillColor;
            };

            var setRectTitle = function (d, i) {
                var data = d[1] != -1 ? d[1] : "No Data";
                return _this.responseData.Title + ": " + data + "\n" + "Lane: " + (d[0] + 1);
            };

            var setRectClass = function (d) {
                var className = "lane-cell";
                if (d[1] == -1 || !_this.axisExists(_this.responseData.Axis))
                    className = className + " no-data";
                return className;
            };

            var cells = lanes.selectAll('rect.lane-cell').data(function (d, index) {
                var arr = [];
                for (var i = 0; i < d.length; i++) {
                    arr.push([index, d[i] != null ? d[i] : -1]);
                }
                return arr;
            });
            cells.enter().append('rect').attr('class', setRectClass).attr('rel', 'tooltip').attr('width', this.laneColumnWidth).attr('height', this.laneHeight).attr('fill', setRectFill).attr('title', setRectTitle).attr('rel', 'tooltip').attr('x', function (d, i) {
                return setRectX(d[0], i);
            }).attr('y', function (d, i) {
                return setRectY(i);
            });
            cells.exit().remove();
            cells.attr('class', setRectClass).transition().duration(this.config.transitionSpeed).attr('width', this.laneColumnWidth).attr('height', this.laneHeight).attr('fill', setRectFill).attr('title', setRectTitle).attr('x', function (d, i) {
                return setRectX(d[0], i);
            }).attr('y', function (d, i) {
                return setRectY(i);
            });
        };

        FlowcellChart.prototype.addGradients = function () {
            var _this = this;
            var gradient = this.svg.append('defs').append('linearGradient').attr('id', 'flowcell-color-scale-gradient').attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');

            gradient.selectAll('stop').data(this.colorScaleDomain, function (d) {
                return d;
            }).enter().append('stop').attr('offset', function (d, i) {
                var percent = "";
                percent = (d * 100) + "%";
                return percent;
            }).attr('style', function (d, i) {
                var color = _this.colorScaleRange[i];
                return "stop-color:" + color + ";stop-opacity: 1;";
            });
        };

        FlowcellChart.prototype.updateAxes = function () {
            var _this = this;
            this.svg.selectAll('.axis').remove();
            this.yAxis = d3.svg.axis();
            this.yAxis.scale(this.scaleY);
            this.yAxis.orient('right');
            this.yAxis.tickSize(0, 0, 0);
            this.yAxis.tickSubdivide(1);

            var axisWidth = Math.min(30, this.config.width / 15);

            var axisGroup = this.svg.append('g').attr('transform', 'translate(' + (this.config.width - 30) + ', 0)').attr('class', 'axis');

            var setAxisFill = function (d) {
                if (!_this.axisExists(_this.responseData.Axis))
                    return '#666';
                else
                    return 'url(#flowcell-color-scale-gradient)';
            };

            var setAxisFillOpacity = function (d) {
                if (!_this.axisExists(_this.responseData.Axis))
                    return '0.25';
                else
                    return '1';
            };

            axisGroup.insert('rect', ':first-child').attr('x', -axisWidth).attr('y', this.config.margins.top).attr('width', axisWidth).attr('height', (this.config.height - (this.config.margins.top + this.config.margins.bottom))).attr('fill', setAxisFill).attr('fill-opacity', setAxisFillOpacity);

            axisGroup.call(this.yAxis);
        };

        FlowcellChart.prototype.axisExists = function (axis) {
            if (!axis || (axis.Min == 0 && axis.Max == 0))
                return false;
            else
                return true;
        };

        FlowcellChart.prototype.updateScales = function () {
            var that = this;

            this.colorScaleDomain = [0, 0.33, 0.66, 1];
            this.colorScaleRange = ["blue", "cyan", "yellow", "orange"];
            this.laneHorizontalSpacing = this.config.width * .015;
            console.log(this.laneHorizontalSpacing);

            this.laneWidth = ((this.config.width - (this.config.margins.left + this.config.margins.right)) / this.responseData.LaneCount) - this.laneHorizontalSpacing;
            this.laneColumnWidth = this.laneWidth / this.responseData.ColumnCount;
            this.laneHeight = ((this.config.height - (this.config.margins.top + this.config.margins.bottom)) / this.responseData.RowCount);

            var domainMin = this.responseData.Axis.Min;
            var domainMax = this.responseData.Axis.Max;

            if (this.axisExists(this.responseData.Axis) && (domainMin == domainMax)) {
                ++domainMax;
                --domainMin;
            }

            this.colorScale = d3.scale.linear().domain([domainMin, domainMax]);
            this.colorScale.domain(this.colorScaleDomain.map(this.colorScale.invert));
            this.colorScale.range(this.colorScaleRange);

            this.scaleY = d3.scale.linear();
            this.scaleY.domain([domainMin, domainMax]);
            this.scaleY.range([this.config.height - this.config.margins.bottom, 0 + this.config.margins.top]);
        };
        FlowcellChart.CHART_NAME = "Flowcell";
        return FlowcellChart;
    })(bs.Chart);
    bs.FlowcellChart = FlowcellChart;
})(bs || (bs = {}));
