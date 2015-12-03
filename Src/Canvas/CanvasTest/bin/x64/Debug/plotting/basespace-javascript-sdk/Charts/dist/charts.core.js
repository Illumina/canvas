var bs;
(function (bs) {
    var Chart = (function () {
        function Chart(responseData, customConfig) {
            this.legendItemHeight = 10;
            this.TOOLTIP_CLASS = "chart-tooltip";
            this.baseGradientDefinitions = [];
            this.config = {
                chartContainerSelector: null,
                fontSize: 10,
                height: 320,
                width: 480,
                margins: {
                    top: 30,
                    left: 80,
                    right: 30,
                    bottom: 50
                },
                strokeWidth: 1,
                preserveAspectRatio: 'xMidYMax',
                transitionSpeed: 250,
                contentType: '',
                labelBarValues: false,
                labelBarPostfix: '',
                color: '#98cc78',
                legend: true,
                svgLegend: false,
                xAxis: {
                    orient: null,
                    tickFormat: function (x) {
                        return bs.Util.trimAndAddUnits(x);
                    },
                    tickSize: null,
                    tickValues: null,
                    ticks: 10,
                    scale: null,
                    tickSubdivide: 1
                },
                yAxis: {
                    orient: 'left',
                    tickFormat: function (y) {
                        return bs.Util.trimAndAddUnits(y);
                    },
                    tickSize: null,
                    tickValues: null,
                    ticks: 8,
                    scale: null,
                    tickSubdivide: 1
                },
                scaleXMaxPadding: 0,
                scaleYMaxPadding: 0,
                scaleXMinPadding: 0,
                scaleYMinPadding: 0,
                scaleY: d3.scale.linear,
                scaleX: d3.scale.linear,
                gradientDefinitions: [],
                stepSize: 20
            };

            $.extend(true, this.config, customConfig);

            this.responseData = responseData;

            if (this.responseData.Series != null) {
                this.allSeries = this.responseData.Series;
                this.allData = this.flattenData();
                this.dataLookupIndex = this.allSeries[0].Fields;
            }

            this.svg = d3.select(this.config.chartContainerSelector).append('svg');
        }
        Chart.prototype.flattenData = function () {
            var flat = [];
            this.allSeries.forEach(function (series, index, allSeries) {
                series.Data.forEach(function (d, j, data) {
                    flat.push(d);
                });
            });
            return flat;
        };

        Chart.prototype.update = function (updatedData, updatedSettings) {
            if (updatedData != null)
                this.responseData = updatedData;
            if (updatedSettings != null)
                $.extend(true, this.config, updatedSettings);
            this.svg.attr('viewBox', '0 0 ' + this.config.width + ' ' + this.config.height).attr('width', '100%').attr('height', '100%').attr('preserveAspectRatio', this.config.preserveAspectRatio);
            this.updateScales();
            this.updateAxes();
            this.updateGradients();
        };

        Chart.prototype.updateScales = function () {
        };

        Chart.prototype.configureAxis = function (axis, params) {
            for (var key in params) {
                var paramValue = params[key];

                if (paramValue === undefined || paramValue === null)
                    continue;

                axis[key](paramValue);
            }
        };

        Chart.prototype.updateAxes = function () {
            this.svg.selectAll('.axis-x, .axis-y').remove();

            this.yAxis = d3.svg.axis();
            this.yAxis.scale(this.scaleY);
            this.yAxis.tickSize(-(this.config.width - ((this.config.margins.left) + this.config.margins.right)), 3, 0);

            this.configureAxis(this.yAxis, this.config.yAxis);

            this.svg.append('g').attr('class', 'axis-y').attr('transform', 'translate(' + this.config.margins.left + ', 0)').call(this.yAxis);

            this.xAxis = d3.svg.axis();
            this.xAxis.scale(this.scaleX);
            this.xAxis.orient('bottom');
            this.xAxis.tickSize(-(this.config.height - ((this.config.margins.top + this.config.margins.bottom))), 3, 0);
            var xRange = d3.max(this.scaleX.domain()) - d3.min(this.scaleX.domain());
            if (xRange < 10) {
                this.xAxis.ticks(xRange);
            }

            this.configureAxis(this.xAxis, this.config.xAxis);

            this.svg.append('g').attr('class', 'axis-x').attr('transform', 'translate(0, ' + (this.config.height - this.config.margins.bottom) + ')').call(this.xAxis);

            this.svg.select('.axis-y').append('text').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'rotate(-90 0 0), translate(' + (-(this.config.height - this.config.margins.top - this.config.margins.bottom) / 2 - this.config.margins.top) + ', ' + (-this.config.margins.left + 20) + ')').text(this.responseData.Axes.Y.Title);

            this.svg.select('.axis-x').append('text').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'translate(' + ((this.config.width - this.config.margins.left - this.config.margins.right) / 2 + this.config.margins.left) + ', ' + (this.config.margins.bottom - 20) + ')').text(this.responseData.Axes.X.Title);
        };

        Chart.prototype.updateZoom = function () {
            var _this = this;
            this.zoom = d3.behavior.zoom().x(this.scaleX).y(this.scaleY).scaleExtent([1, 10]).on("zoom", function (scale, translate) {
                _this.updateAxes();
            });
            this.svg.call(this.zoom);
        };

        Chart.prototype.updateGradients = function () {
            var allGradients = this.baseGradientDefinitions;

            if ("gradientDefinitions" in this.config) {
                allGradients = allGradients.concat(this.config.gradientDefinitions);
            }

            if (allGradients.length == 0)
                return;

            var defs = this.svg.selectAll('defs.gradients').data([allGradients]);
            defs.enter().append('defs').attr('class', 'gradients');
            defs.exit().remove();

            var gradients = defs.selectAll('linearGradient').data(function (gradientDef) {
                return gradientDef;
            });
            gradients.enter().append('linearGradient');
            gradients.exit().remove();
            gradients.attr('id', function (gradientDef) {
                return gradientDef.id;
            }).attr('x1', function (gradientDef) {
                return gradientDef.direction.toLowerCase() === Chart.GRADIENT_DIRECTION_HORIZONTAL ? '100%' : '0%';
            }).attr('x2', '0%').attr('y1', function (gradientDef) {
                return gradientDef.direction.toLowerCase() === Chart.GRADIENT_DIRECTION_VERTICAL ? '100%' : '0%';
            }).attr('y2', '0%');

            var stops = gradients.selectAll('stop').data(function (gradientOption) {
                return gradientOption.stops;
            });
            stops.enter().append('stop');
            stops.exit().remove();
            stops.attr('offset', function (stop) {
                return stop.offset;
            }).attr('style', function (stop) {
                return 'stop-color: ' + stop.color + '; stop-opacity: ' + stop.opacity + ';';
            });
        };

        Chart.prototype.addAxisLabels = function () {
        };

        Chart.prototype.addLegend = function () {
            var that = this;

            var legendWidth = 50;
            var legendHeight = 50;
            this.legend = d3.select(this.config.chartContainerSelector).append('div').attr('class', 'legend');

            this.legend.append('ul');

            this.svgLegend = this.svg.append('g').attr('class', 'legend').attr('transform', 'translate(' + (this.config.width - this.config.margins.right + 5) + ')');
        };

        Chart.prototype.postRender = function () {
            if (this.config.fontSize != null) {
            }
        };

        Chart.prototype.callTooltips = function (tooltipSelector) {
            if (typeof tooltipSelector === "undefined") { tooltipSelector = "." + this.TOOLTIP_CLASS; }
            if ("powerTip" in $.fn) {
                $(this.config.chartContainerSelector + " " + tooltipSelector).powerTip({
                    smartPlacement: true,
                    mouseOnToPopup: true,
                    followMouse: true,
                    placement: "n"
                });
            }
        };
        Chart.TOOLTIP_TEMPLATE = "<div class='bs-tooltip-content dark'><%= content %></div>";
        Chart.GRADIENT_DIRECTION_HORIZONTAL = "h";
        Chart.GRADIENT_DIRECTION_VERTICAL = "v";

        Chart.MIN_LOG_INPUT = 1;
        return Chart;
    })();
    bs.Chart = Chart;
})(bs || (bs = {}));
