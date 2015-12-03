var __extends = this.__extends || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    __.prototype = b.prototype;
    d.prototype = new __();
};
var bs;
(function (bs) {
    var CoverageBAFChart = (function (_super) {
        __extends(CoverageBAFChart, _super);
        function CoverageBAFChart(responseData, customConfig) {
            var _this = this;
            _super.call(this, responseData, customConfig);
            this.scatterPlotFillOpacity = .75;
            this.defs = null;
            this.config = $.extend(true, {
                yAxes: responseData.StackedAxes.Y.map(function (x) {
                    return _this.config.yAxis;
                }),
                scaleYs: responseData.StackedAxes.Y.map(function (x) {
                    return d3.scale.linear;
                }),
                interChartYPadding: 10
            }, this.config);

            this.update(responseData);
        }
        CoverageBAFChart.prototype.update = function (updatedData, updatedSettings) {
            this.updateDefs();

            _super.prototype.update.call(this, updatedData, updatedSettings);

            var responseData = this.responseData;

            var coverageBAFChart = this.svg.selectAll('g.coveragebafchart').data([responseData]);
            coverageBAFChart.enter().append('g').attr('class', 'coveragebafchart');
            coverageBAFChart.exit().remove();

            for (var i = 0; i < responseData.StackedSeries.length; i++) {
                // render series differently depending on type: line, scatter, heatmap.
                var lineSeries = [];
                var scatterSeries = [];
                var heatSeries = [];

                responseData.StackedSeries[i].forEach(function (series) {
                    switch (series.Type.toLowerCase()) {
                        case CoverageBAFChart.SERIES_TYPE_LINE:
                            lineSeries.push(series);
                            break;
                        case CoverageBAFChart.SERIES_TYPE_SCATTER:
                            scatterSeries.push(series);
                            break;
                        case CoverageBAFChart.SERIES_TYPE_HEATMAP:
                            heatSeries.push(series);
                            break;
                        default:
                            break;
                    }
                });

                if (heatSeries.length > 0 && (lineSeries.length > 0 || scatterSeries.length > 0)) {
                    throw "Cannot mix heat series with line or scatter series.";
                }

                if (heatSeries.length > 0) {
                    this.updateHeatSeries(coverageBAFChart, heatSeries, i);
                } else {
                    this.updateLineSeries(coverageBAFChart, lineSeries, i);
                    this.updateScatterSeries(coverageBAFChart, scatterSeries, i);
                }
            }

            this.updateLegend(i);

            this.callTooltips("[title]");
        };

        CoverageBAFChart.prototype.getYRangeStep = function (config) {
            return (config.height - config.margins.top - config.margins.bottom - (config.scaleYs.length - 1) * config.interChartYPadding) / config.scaleYs.length;
        };

        CoverageBAFChart.prototype.updateScaleX = function () {
            var responseData = this.responseData;
            var config = this.config;

            this.scaleX = config.scaleX();
            this.scaleX1 = config.scaleX();
            if (config.scaleX === d3.scale.log) {
                this.scaleX.clamp(true);
                this.scaleX1.clamp(true);
                var domainXMin = responseData.StackedAxes.X.Min == 0 ? bs.Chart.MIN_LOG_INPUT : responseData.StackedAxes.X.Min;
                this.scaleX.domain([domainXMin, responseData.StackedAxes.X.Max]);
                this.scaleX1.domain([0, responseData.StackedAxes.X.Max - domainXMin]);
            } else {
                this.scaleX.domain([responseData.StackedAxes.X.Min, responseData.StackedAxes.X.Max]);
                this.scaleX1.domain([0, responseData.StackedAxes.X.Max - responseData.StackedAxes.X.Min]);
            }
            this.scaleX.range([config.margins.left, config.width - config.margins.right]);
            this.scaleX1.range([0, config.width - config.margins.left - config.margins.right]);
            if (config.scaleX && config.scaleX.nice === true) {
                this.scaleX.nice();
                this.scaleX1.nice();
            }
        };

        CoverageBAFChart.prototype.updateScaleYs = function () {
            var responseData = this.responseData;
            var config = this.config;

            var yRangeStep = this.getYRangeStep(config);
            this.scaleYs = [];
            for (var i = 0; i < responseData.StackedAxes.Y.length; i++) {
                var configScaleY = config.scaleYs[i];
                var scaleY = configScaleY();
                this.scaleYs.push(scaleY);
                if (configScaleY == d3.scale.log) {
                    scaleY.clamp(true);
                    var domainYMin = responseData.StackedAxes.Y[i].Min == 0 ? bs.Chart.MIN_LOG_INPUT : responseData.StackedAxes.Y[i].Min;
                    scaleY.domain([domainYMin, responseData.StackedAxes.Y[i].Max]);
                } else {
                    scaleY.domain([responseData.StackedAxes.Y[i].Min, responseData.StackedAxes.Y[i].Max]);
                }
                scaleY.range([
                    config.height - config.margins.bottom - i * yRangeStep - i * config.interChartYPadding,
                    config.height - config.margins.bottom - (i + 1) * yRangeStep - i * config.interChartYPadding]);
                if (configScaleY && configScaleY.nice === true) {
                    scaleY.nice();
                }
            }
        };

        CoverageBAFChart.prototype.updateScales = function () {
            var responseData = this.responseData;
            var config = this.config;
            var yRangeStep = this.getYRangeStep(config);

            // x scale
            this.updateScaleX();

            // y scale(s)
            this.updateScaleYs();

            // color/heat scale if any
            var fontSize = config.fontSize ? config.fontSize : CoverageBAFChart.DEFAULT_FONT_SIZE;
            this.maxHeatAxisHeight = 20 * fontSize;
            this.heatAxisHeight = Math.min(this.maxHeatAxisHeight, yRangeStep);
            this.heatAxisMargin = (yRangeStep - this.heatAxisHeight) / 2;

            this.colorScales = [];
            this.heatScales = [];
            var plotYStart = config.height - config.margins.bottom;
            var plotYEnd;
            for (var i = 0; i < responseData.StackedAxes.Y.length; i++) {
                plotYEnd = plotYStart - yRangeStep;
                var series = responseData.StackedSeries[i];
                if (series.length > 0 && series[0].Type.toLowerCase() == CoverageBAFChart.SERIES_TYPE_HEATMAP) {
                    var min = Infinity;
                    var max = -Infinity;
                    series.forEach(function (s) {
                        s.Data.forEach(function (d) {
                            max = Math.max(max, _.max(d[s.Fields.y]));
                            min = Math.min(min, _.min(d[s.Fields.y]));
                        });
                    });

                    // color scale
                    this.colorScales.push(d3.scale.linear());
                    if (responseData.ColorScaleDomain && responseData.ColorScaleDomain[i]) {
                        this.colorScales[i].domain(responseData.ColorScaleDomain[i]);
                        this.colorScales[i].range(responseData.ColorScaleRange[i]);
                    } else {
                        this.colorScales[i].range([
                            responseData.ColorScaleRange[i][0],
                            responseData.ColorScaleRange[i][responseData.ColorScaleRange[i].length - 1]]);
                    }

                    // heat scale
                    this.heatScales.push(d3.scale.linear());
                    this.heatScales[i].domain([min, max]);
                    this.heatScales[i].range([plotYStart - this.heatAxisMargin, plotYEnd + this.heatAxisMargin]);
                } else {
                    this.colorScales.push(null);
                    this.heatScales.push(null);
                }
                plotYStart -= (yRangeStep + config.interChartYPadding);
            }
        };

        CoverageBAFChart.prototype.updateAxes = function () {
            var _this = this;
            var responseData = this.responseData;
            var config = this.config;
            var fontSize = config.fontSize ? config.fontSize : CoverageBAFChart.DEFAULT_FONT_SIZE;

            // x axis
            this.svg.selectAll('.axis-x').remove();

            this.xAxis = d3.svg.axis();
            this.xAxis.scale(this.scaleX);
            this.xAxis.orient('bottom');
            this.xAxis.tickSize(-(config.height - ((config.margins.top + config.margins.bottom))), 3, 0);
            var xRange = d3.max(this.scaleX.domain()) - d3.min(this.scaleX.domain());
            if (xRange < 10) {
                this.xAxis.ticks(xRange);
            }

            this.configureAxis(this.xAxis, config.xAxis);

            var tickValues = this.xAxis.tickValues();

            this.svg.append('g').attr('class', 'axis-x').attr('transform', 'translate(0, ' + (config.height - config.margins.bottom + 1) + ')').call(this.xAxis).selectAll('text').style('text-anchor', 'end').style('font-size', Math.round(fontSize * 0.8) + 'px').attr("dx", "-.8em").attr("dy", function (d, i) {
                if (responseData.StackedAxes.X.GroupLabelPositions) {
                    var diff = responseData.StackedAxes.X.GroupLabelPositions[i] - tickValues[i];
                    return _this.scaleX1(diff);
                } else {
                    return ".1em";
                }
            }).attr('transform', 'rotate(-90)').text(function (d, i) {
                if (responseData.StackedAxes.X.GroupLabels) {
                    return responseData.StackedAxes.X.GroupLabels[i];
                } else {
                    return d;
                }
            });

            this.svg.select('.axis-x path').style('fill', 'none').style('stroke', 'gray').style('opacity', '1').style('shape-rendering', 'crispEdges');

            this.svg.select('.axis-x').append('text').style('font-size', fontSize + 'px').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'translate(' + ((config.width - config.margins.left - config.margins.right) / 2 + config.margins.left) + ', ' + (config.margins.bottom - 20) + ')').text(responseData.StackedAxes.X.Title);

            // y axes and heat axes
            var yRangeStep = this.getYRangeStep(config);
            var plotYStart = config.height - config.margins.bottom;
            var plotYEnd;
            this.yAxes = [];
            this.heatAxes = [];
            for (var i = 0; i < this.scaleYs.length; i++) {
                var axisClass = 'axis-y chart-' + i;
                var axisClassSelector = '.' + axisClass.replace(/ /g, '.');
                this.svg.selectAll(axisClassSelector).remove();

                this.yAxes.push(d3.svg.axis());
                this.yAxes[i].scale(this.scaleYs[i]);
                this.yAxes[i].tickSize(-(config.width - ((config.margins.left) + config.margins.right)), 3, 0);

                this.configureAxis(this.yAxes[i], config.yAxes[i]);

                this.svg.append('g').attr('class', axisClass).attr('transform', 'translate(' + (config.margins.left - 1) + ', 0)').call(this.yAxes[i]).selectAll('text').style('font-size', Math.round(fontSize * 0.8) + 'px');

                this.svg.select(axisClassSelector + ' path').style('fill', 'none').style('stroke', 'gray').style('opacity', '1').style('shape-rendering', 'crispEdges');

                plotYEnd = plotYStart - yRangeStep;
                this.svg.select(axisClassSelector).append('text').style('font-size', fontSize + 'px').attr('class', 'label').attr('y', 0).attr('x', 0).attr('text-anchor', 'middle').attr('transform', 'rotate(-90 0 0), translate(' + (-(plotYStart + plotYEnd) / 2) + ', ' + (-config.margins.left + 2 * config.fontSize) + ')').text(responseData.StackedAxes.Y[i].Title);

                // heat axis
                if (this.heatScales[i] !== null) {
                    var heatAxisClass = 'heat-axis-' + i;
                    var heatAxisClassSelector = '.' + heatAxisClass.replace(/ /g, '.');
                    this.svg.selectAll(heatAxisClassSelector).remove();

                    this.heatAxes.push(d3.svg.axis());
                    this.heatAxes[i].scale(this.heatScales[i]).orient('right').ticks(7);

                    var rightMargin1 = 1.5 * fontSize;
                    var rightMargin2 = Math.max(0, this.config.margins.right - rightMargin1);
                    var heatAxisGroup = this.svg.append('g').attr('class', heatAxisClass).attr('transform', 'translate(' + (this.config.width - rightMargin2) + ', 0)');

                    heatAxisGroup.call(this.heatAxes[i]);
                    heatAxisGroup.selectAll('g.tick.major text').style('font-size', Math.max(6, fontSize * 0.8) + 'px');

                    heatAxisGroup.insert('rect', ':first-child').attr('class', 'axis-reference-gradient').attr('x', -fontSize).attr('y', plotYEnd + this.heatAxisMargin).attr('width', 1.2 * fontSize).attr('height', this.heatAxisHeight).attr('fill', 'url(#' + this.referenceGradientIds[i] + ')');

                    heatAxisGroup.append('text').attr('class', 'heat-axis-label-' + i).attr('transform', 'rotate(-90 0 0), translate(' + (-(plotYStart + plotYEnd) / 2) + ', ' + Math.min(3 * fontSize, rightMargin2 - fontSize) + ')').attr('x', 0).attr('y', 0).style('text-anchor', 'middle').style('font-size', fontSize + 'px').text(responseData.ColorScaleTitle[i] == null ? 'Color Scale' : responseData.ColorScaleTitle[i]);

                    this.svg.selectAll(heatAxisClassSelector + ' path').style("fill", "none");
                } else {
                    this.heatAxes.push(null);
                }

                plotYStart -= (yRangeStep + config.interChartYPadding);
            }

            // tick style
            this.svg.selectAll('g line.tick.minor').style('stroke', '#cccccc').style('stroke-width', '0.25');
            this.svg.selectAll('g.tick.major line').style('stroke', '#cccccc').style('stroke-width', '0.5');
        };

        CoverageBAFChart.prototype.updateLineSeries = function (coverageBAFChart, lineSeries, chartIndex) {
            var chart = this;
            var seriesClass = 'series line chart-' + chartIndex;
            var seriesClassSelector = '.' + seriesClass.replace(/ /g, '.');

            var series = coverageBAFChart.selectAll('g' + seriesClassSelector).data(lineSeries);
            series.enter().append('g');
            series.exit().remove();
            series.attr('class', function (d) {
                return seriesClass + ' ' + d.Color.toLowerCase().replace('#', '_');
            });

            // fields is for scaline the x and y coordinates, it is initialized before lineFn is called
            var fields;
            var lineFn = d3.svg.line().x(function (d) {
                return chart.scaleX(d[fields.x]);
            }).y(function (d) {
                return chart.scaleYs[chartIndex](d[fields.y]);
            });

            var line = series.selectAll('path').data(function (d) {
                return [d];
            });
            line.enter().append('path');
            line.exit().remove();
            line.attr("data-powertip", function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "powertip" in parentData.Accessors) {
                    return parentData.Accessors.powertip(d);
                } else {
                    return d.Title;
                }
            }).attr('title', function (d) {
                return d.Title;
            }).transition().duration(this.config.transitionSpeed).attr('d', function (d) {
                fields = d.Fields;
                return lineFn(d.Data);
            }).style('fill', 'none').style('stroke', function (d) {
                return d.Color.toLowerCase();
            }).style('stroke-width', this.config.strokeWidth);
        };

        CoverageBAFChart.prototype.updateScatterSeries = function (coverageBAFChart, scatterSeries, chartIndex) {
            var chart = this;

            var seriesClass = 'series scatter chart-' + chartIndex;
            var seriesClassSelector = '.' + seriesClass.replace(/ /g, '.');

            this.svg.selectAll('g' + seriesClassSelector).remove();
            var series = coverageBAFChart.selectAll('g' + seriesClassSelector).data(scatterSeries);
            series.enter().append('g');
            series.exit().remove();
            series.attr('class', seriesClass);

            var circle = series.selectAll('circle').data(function (d) {
                return d.Data;
            });
            circle.enter().append('circle');
            circle.exit().remove();
            circle.attr("data-powertip", function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "powertip" in parentData.Accessors) {
                    return parentData.Accessors.powertip(d);
                } else {
                    return d.Title;
                }
            }).attr('title', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                return d[parentData.Fields.tooltipTitle] + "\n" + d[parentData.Fields.y];
            }).transition().duration(this.config.transitionSpeed).attr('cx', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                return chart.scaleX(d[parentData.Fields.x]);
            }).attr('cy', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                return chart.scaleYs[chartIndex](d[parentData.Fields.y]);
            }).style('fill', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "fill" in parentData.Accessors) {
                    return parentData.Accessors.fill(d);
                } else {
                    return parentData.Color.toLowerCase();
                }
            }).style('fill-opacity', this.scatterPlotFillOpacity).style('stroke', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "stroke" in parentData.Accessors) {
                    return parentData.Accessors.fill(d);
                } else if (parentData.Stroke) {
                    return parentData.Stroke.toLowerCase();
                } else {
                    return parentData.Color.toLowerCase();
                }
            }).style('stroke-width', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ('CircleStrokeWidth' in parentData) {
                    return parentData.CircleStrokeWidth;
                } else {
                    return CoverageBAFChart.DEFAULT_SCATTER_CIRCLE_STROKE_WIDTH;
                }
            }).attr('r', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ('CircleRadius' in parentData) {
                    return parentData.CircleRadius;
                } else {
                    return CoverageBAFChart.DEFAULT_SCATTER_CIRCLE_RADIUS;
                }
            });
        };

        CoverageBAFChart.prototype.updateHeatSeries = function (coverageBAFChart, heatSeries, chartIndex) {
            var _this = this;
            var chart = this;

            var seriesClass = 'series heat chart-' + chartIndex;
            var seriesClassSelector = '.' + seriesClass.replace(/ /g, '.');

            var series = coverageBAFChart.selectAll('g' + seriesClassSelector).data(heatSeries);
            series.enter().append('g');
            series.exit().remove();
            series.attr('class', function (d) {
                return seriesClass;
            });

            var column = series.selectAll('g.column').data(function (d) {
                return d.Data;
            });
            column.enter().append('g');
            column.exit().remove();
            column.attr('class', 'column');

            var yRange = Math.abs(chart.scaleYs[chartIndex].domain()[1] - chart.scaleYs[chartIndex].domain()[0]);
            var rect = column.selectAll('rect').data(function (d, i) {
                var series = d3.select(this.parentNode.parentNode).datum();
                var rectHeight = yRange / d[series.Fields.y].length;
                var dataArr = [];
                for (var j = 0; j < d[series.Fields.y].length; j++) {
                    var value = d[series.Fields.y][j];
                    var color = '' + chart.colorScales[chartIndex](value);
                    if (color == '#ffffff') {
                        continue;
                    }
                    var width = Math.ceil(chart.scaleX(d[series.Fields.x][1]) - chart.scaleX(d[series.Fields.x][0]));
                    var height = Math.ceil(Math.abs(chart.scaleYs[chartIndex](0) - chart.scaleYs[chartIndex](rectHeight)));
                    var rectData = {
                        x: chart.scaleX(d[series.Fields.x][0]),
                        y: chart.scaleYs[chartIndex]((j + 1) * rectHeight),
                        width: Math.max(1, width),
                        height: Math.max(1, height),
                        value: value
                    };
                    dataArr.push(rectData);
                }
                return dataArr;
            });
            rect.enter().append('rect');
            rect.exit().remove();

            rect.transition().duration(this.config.transitionSpeed).attr('width', function (d) {
                return d.width + 'px';
            }).attr('height', function (d) {
                return d.height + 'px';
            }).attr('x', function (d) {
                return d.x;
            }).attr('y', function (d) {
                return d.y;
            });

            rect.style('fill', function (d) {
                return _this.colorScales[chartIndex](d.value);
            }).style('fill-opacity', CoverageBAFChart.DEFAULT_RECT_FILL_OPACITY).style('stroke', function (d) {
                return _this.colorScales[chartIndex](d.value);
            }).style('stroke-width', CoverageBAFChart.DEFAULT_RECT_STROKE_WIDTH);

            rect.append('title').text(function (d) {
                return d.value.toFixed(2);
            });
        };

        /*
        updateGradients(): void {
        }
        */
        CoverageBAFChart.prototype.updateDefs = function () {
            var responseData = this.responseData;
            if (!('ColorScaleDomain' in responseData)) {
                return;
            }
            if (!('ColorScaleRange' in responseData)) {
                return;
            }

            if (this.defs === null) {
                this.referenceGradientIds = [];
                this.defs = this.svg.append('defs'); // for reference heat scales
                for (var i = 0; i < responseData.StackedSeries.length; i++) {
                    if (responseData.ColorScaleDomain[i] && responseData.ColorScaleRange[i] && responseData.ColorScaleDomain[i].length > 0 && responseData.ColorScaleRange[i].length > 0) {
                        // make sure this id is unique
                        var referenceGradientId = this.config.chartContainerSelector.substring(1) + '-reference-gradient-' + i;
                        this.referenceGradientIds.push(referenceGradientId);
                        var referenceGradient = this.defs.append('linearGradient').attr('id', referenceGradientId).attr('x1', '0%').attr('x2', '0%').attr('y1', '100%').attr('y2', '0%');
                        var colorScaleDomain = responseData.ColorScaleDomain[i];
                        var domainSize = Math.abs(colorScaleDomain[0] - colorScaleDomain[colorScaleDomain.length - 1]);

                        // for the color scale axis
                        colorScaleDomain = colorScaleDomain.map(function (value, i, arr) {
                            return Math.abs(value - arr[0]) / domainSize;
                        });
                        referenceGradient.selectAll('stop').data(colorScaleDomain).enter().append('stop').attr('offset', function (d) {
                            return d * 100 + "%";
                        }).attr('style', function (d, j) {
                            return "stop-color: " + responseData.ColorScaleRange[i][j] + "; stop-opacity:1";
                        });
                    } else {
                        this.referenceGradientIds.push(null);
                    }
                }
            }
        };

        /*
        getColorToLegend(series: IChartResponseSeries[]): any {
        var colorToLegend = {};
        for (var i = 0; i < series.length; i++) {
        colorToLegend[series[i].Color.toLowerCase()] = series[i].Title;
        }
        
        return colorToLegend;
        }
        */
        CoverageBAFChart.prototype.getMeanLegendItemPixelsPerChar = function (svgLegendList, fontSize) {
            var nChars = 0;
            var svgItems = svgLegendList.selectAll('text').data(function (d) {
                return d.map(function (s) {
                    nChars += s.Title.length;
                    return { title: s.Title, color: s.Color.toLowerCase() };
                });
            });
            svgItems.enter().append('text');
            svgItems.exit().remove();
            svgItems.text(function (d) {
                return d.title;
            }).attr('x', fontSize * 1.5).attr('y', function (d, i) {
                return i * fontSize;
            }).attr('alignment-baseline', 'middle').style('font-size', fontSize + 'px').style('fill', function (d) {
                return d.color;
            });

            var nPixels = 0;
            svgItems[0].forEach(function (text) {
                nPixels += text.getBBox().width;
            });

            return nChars > 0 ? nPixels / nChars : 0;
        };

        CoverageBAFChart.prototype.updateLegend = function (chartIndex) {
            var responseData = this.responseData;
            var config = this.config;
            var fontSize = config.fontSize ? config.fontSize : CoverageBAFChart.DEFAULT_FONT_SIZE;
            var yRangeStep = this.getYRangeStep(config);

            this.svgLegends = [];
            var plotYStart = config.height - config.margins.bottom;
            var plotYEnd;
            for (var i = 0; i < responseData.StackedSeries.length; i++) {
                plotYEnd = plotYStart - yRangeStep;

                var legendClass = 'legend-' + i;
                var legendClassSelector = '.' + legendClass.replace(/ /g, '.');
                this.svg.select(legendClassSelector).remove();

                var series = responseData.StackedSeries[i];

                // no legend for heatmap (for now)
                if (series[0].Type.toLowerCase() == CoverageBAFChart.SERIES_TYPE_HEATMAP) {
                    this.svgLegends.push(null);
                } else {
                    if (series[0].Type.toLocaleLowerCase() == CoverageBAFChart.SERIES_TYPE_LINE) {
                        var newSeries = [];
                        var seen = [];
                        for (var j = 0; j < series.length; j++) {
                            if (!(seen.indexOf(series[j].Title) >= 0)) {
                                seen.push(series[j].Title);
                                newSeries.push(series[j]);
                            }
                        }
                        series = newSeries;
                    }
                    this.svgLegends.push(this.svg.append('g').attr('class', legendClass));
                    this.svgLegends[i].attr('transform', 'translate(' + (config.width - config.margins.right + fontSize / 2) + ')');

                    var legendHeight = (2 + series.length) * fontSize;
                    var margin = Math.max(0, yRangeStep - legendHeight) / 2;

                    // legend items
                    var svgLegendList = this.svgLegends[i].selectAll('g').data([series]);
                    svgLegendList.enter().append('g');
                    svgLegendList.exit().remove();
                    svgLegendList.attr('transform', 'translate(' + fontSize / 2 + ', ' + (plotYEnd + margin + fontSize) + ')');

                    var pixelsPerChar = this.getMeanLegendItemPixelsPerChar(svgLegendList, fontSize);
                    var maxLegendItemLength = Math.floor((config.margins.right - 2.5 * fontSize) / pixelsPerChar);

                    var title2Tooltip = {};
                    var svgItems = svgLegendList.selectAll('text').data(function (d) {
                        return d.map(function (s) {
                            var displayTitle = s.Title;
                            if (displayTitle.length > maxLegendItemLength) {
                                displayTitle = displayTitle.substr(0, maxLegendItemLength - 3) + '...';
                            }
                            while (displayTitle in title2Tooltip) {
                                displayTitle += ' ';
                            }
                            title2Tooltip[displayTitle] = s.Title;
                            return { title: displayTitle, color: s.Color.toLowerCase() };
                        });
                    });
                    svgItems.enter().append('text');
                    svgItems.exit().remove();
                    svgItems.text(function (d) {
                        return d.title;
                    }).attr('x', fontSize * 1.5).attr('y', function (d, i) {
                        return i * fontSize;
                    }).attr('alignment-baseline', 'middle').style('font-size', fontSize + 'px').style('fill', function (d) {
                        return d.color;
                    });
                    svgItems.append('title').text(function (d) {
                        return title2Tooltip[d.title];
                    });

                    var svgItemRects = svgLegendList.selectAll('rect').data(function (d) {
                        return d.map(function (s) {
                            return s.Color.toLowerCase();
                        });
                    });
                    svgItemRects.enter().append('rect');
                    svgItemRects.exit().remove();
                    svgItemRects.attr('width', fontSize + 'px').attr('height', (fontSize * 0.6) + 'px').attr('x', 0).attr('y', function (d, i) {
                        return (i - 0.4) * fontSize;
                    });
                    svgItemRects.style('fill', function (d) {
                        return d;
                    }).style('fill-opacity', CoverageBAFChart.DEFAULT_RECT_FILL_OPACITY).style('stroke-width', 0);
                }

                plotYStart -= (yRangeStep + config.interChartYPadding);
            }
        };
        CoverageBAFChart.CHART_NAME = "CoverageBAFChart";
        CoverageBAFChart.SERIES_TYPE_LINE = "line";
        CoverageBAFChart.SERIES_TYPE_SCATTER = "scatter";
        CoverageBAFChart.SERIES_TYPE_HEATMAP = "heat";

        CoverageBAFChart.DEFAULT_SCATTER_CIRCLE_RADIUS = 3.5;
        CoverageBAFChart.DEFAULT_SCATTER_CIRCLE_STROKE_WIDTH = 0.75;
        CoverageBAFChart.DEFAULT_FONT_SIZE = 10;
        CoverageBAFChart.DEFAULT_RECT_FILL_OPACITY = 1;
        CoverageBAFChart.DEFAULT_RECT_STROKE_WIDTH = 0;
        return CoverageBAFChart;
    })(bs.Chart);
    bs.CoverageBAFChart = CoverageBAFChart;
})(bs || (bs = {}));
