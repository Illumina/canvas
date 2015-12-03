var __extends = this.__extends || function (d, b) {
    for (var p in b) if (b.hasOwnProperty(p)) d[p] = b[p];
    function __() { this.constructor = d; }
    __.prototype = b.prototype;
    d.prototype = new __();
};
var bs;
(function (bs) {
    var LineChart = (function (_super) {
        __extends(LineChart, _super);
        function LineChart(responseData, customConfig) {
            _super.call(this, responseData, customConfig);
            this.scatterPlotFillOpacity = .75;
            this.scatterPlotStrokeWidth = 0.75;

            this.baseGradientDefinitions = [LineChart.DEFAULT_SCATTER_FILL_GRADIENT];

            this.update(responseData, customConfig);

            this.bindMouseEvents();
        }
        LineChart.prototype.update = function (updatedData, updatedSettings) {
            _super.prototype.update.call(this, updatedData, updatedSettings);

            var lineSeries = [], scatterSeries = [];

            _.each(this.responseData.Series, function (series) {
                switch (series.Type.toLowerCase()) {
                    case LineChart.SERIES_TYPE_LINE:
                        lineSeries.push(series);
                        break;
                    case LineChart.SERIES_TYPE_SCATTER:
                        scatterSeries.push(series);
                        break;
                    default:
                        break;
                }
            });

            var linechart = this.svg.selectAll('g.linechart').data([this.responseData]);
            linechart.enter().append('g').attr('class', 'linechart');
            linechart.exit().remove();

            this.updateLineSeries(linechart, lineSeries);
            this.updateScatterSeries(linechart, scatterSeries);
            this.updateLegend();

            this.callTooltips("[title]");
        };

        LineChart.prototype.updateLineSeries = function (linechart, lineSeries) {
            var _this = this;
            var lookupIndex;

            var series = linechart.selectAll('g.series.line').data(lineSeries);
            series.enter().append('g');
            series.exit().remove();
            series.datum(function (d) {
                lookupIndex = d.Fields;
                return d;
            }).attr('class', function (d) {
                return 'series line ' + d.Color.toLowerCase().replace('#', '_');
            });

            var lineFn = d3.svg.line().x(function (d) {
                return _this.scaleX(d[lookupIndex.x]);
            }).y(function (d) {
                return _this.scaleY(d[lookupIndex.y]);
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
                return lineFn(d.Data);
            }).style('fill', 'none').style('stroke', function (d) {
                return d.Color.toLowerCase();
            }).style('stroke-width', this.config.strokeWidth);
        };

        LineChart.prototype.updateScatterSeries = function (linechart, scatterSeries) {
            var _this = this;
            var lookupIndex, powertipAccessor, fillAccessor;

            var scatterPlotCircleRadius;

            var series = linechart.selectAll('g.series.scatter').data(scatterSeries);
            series.enter().append('g');
            series.exit().remove();
            series.datum(function (d) {
                lookupIndex = d.Fields;
                scatterPlotCircleRadius = (d.PointSize === undefined ? LineChart.DEFAULT_SCATTER_CIRCLE_RADIUS : Math.round(d.PointSize) / 2);
                return d;
            }).attr('class', function (d) {
                return 'series scatter';
            });

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
                return d[lookupIndex.tooltipTitle] + "\n" + d[lookupIndex.y];
            }).transition().duration(this.config.transitionSpeed).attr('cx', function (d) {
                return _this.scaleX(d[lookupIndex.x]);
            }).attr('cy', function (d) {
                return _this.scaleY(d[lookupIndex.y]);
            }).style('fill', function (d) {
                var parentData = d3.select(this.parentNode).datum();
                if ("Accessors" in parentData && "fill" in parentData.Accessors) {
                    return parentData.Accessors.fill(d);
                } else {
                    return parentData.Color.toLowerCase();
                }
            }).style('fill-opacity', this.scatterPlotFillOpacity).style('stroke-width', this.scatterPlotStrokeWidth).attr('r', scatterPlotCircleRadius);
        };

        LineChart.prototype.updateLegend = function () {
            if (this.legend == null)
                this.legend = d3.select(this.config.chartContainerSelector).append('div').attr('class', 'legend');
            if (this.svgLegend == null) {
                this.svgLegend = this.svg.append('g').attr('class', 'legend');
            }

            if (!this.config.legend) {
                this.legend.style('display', 'none');
                this.svgLegend.style('display', 'none');
            } else {
                this.legend.style('display', this.config.svgLegend ? 'none' : 'block');
                this.svgLegend.style('display', this.config.svgLegend ? 'block' : 'none');
            }

            var filterBases = function (seriesArr, includeOrExclude) {
                var arr = seriesArr.filter(function (series) {
                    if (!("Title" in series))
                        return !includeOrExclude;
                    else if ((series.Title.toLowerCase() == "a") || (series.Title.toLowerCase() == "c") || (series.Title.toLowerCase() == "g") || (series.Title.toLowerCase() == "t") || (series.Title.toLowerCase() == "u"))
                        return includeOrExclude;
                    else
                        return !includeOrExclude;
                });

                return arr;
            };

            var legendList = this.legend.selectAll('ul.bases').data([this.responseData]);
            legendList.enter().append('ul').attr('class', 'bases');
            legendList.exit().remove();

            var baseLi = legendList.selectAll('li').data(function (d) {
                return filterBases(d.Series, true);
            });
            baseLi.enter().append('li');
            baseLi.exit().remove();
            baseLi.text(function (d) {
                return d.Title;
            }).attr('class', function (d) {
                return d.Color.toLowerCase().replace('#', '_');
            }).style('color', function (d) {
                return d.Color.toLowerCase();
            });

            var extraList = this.legend.selectAll('ul.extra').data([this.responseData]);
            extraList.enter().append('ul');
            extraList.exit().remove();
            extraList.attr('class', 'extra');

            var extraItems = extraList.selectAll('li').data(function (d) {
                return filterBases(d.Series, false);
            });
            extraItems.enter().append('li');
            extraItems.exit().remove();
            extraItems.text(function (d) {
                return d.Title;
            }).attr('class', function (d) {
                return d.Color.toLowerCase().replace('#', '_');
            }).style('color', function (d) {
                return d.Color.toLowerCase();
            });

            var that = this;
            var baseSeries = filterBases(this.responseData.Series, true);
            var extraSeries = filterBases(this.responseData.Series, false);

            this.svgLegend.attr('transform', 'translate(' + (this.config.width - this.config.margins.right + 5) + ')');

            var fontSize = this.config.fontSize != null ? this.config.fontSize : 10;
            var maxTitleLength = Math.max.apply(null, this.responseData.Series.map(function (series) {
                return series.Title.length;
            }));
            maxTitleLength = Math.max(maxTitleLength, baseSeries.length * 2);

            var svgLegendData = {
                x: 0, y: 0, rx: 10, ry: 10,
                width: Math.min(this.config.margins.right - 5, maxTitleLength * fontSize + fontSize / 2),
                height: Math.min(this.config.height, (extraSeries.length + (baseSeries.length > 0 ? 1 : 0) + 1) * fontSize)
            };

            var svgLegendBox = this.svgLegend.selectAll('rect').data([svgLegendData]);
            svgLegendBox.enter().append('rect');
            svgLegendBox.exit().remove();
            svgLegendBox.attr('x', function (d) {
                return d.x;
            }).attr('y', function (d) {
                return d.y;
            }).attr('rx', function (d) {
                return d.rx;
            }).attr('ry', function (d) {
                return d.ry;
            }).attr('width', function (d) {
                return d.width;
            }).attr('height', function (d) {
                return d.height;
            }).style('stroke', 'black').style('stroke-width', 1).style('opacity', 0.5).style('fill', 'none');

            var svgBaseLegendList = this.svgLegend.selectAll('g.bases').data([baseSeries]);
            svgBaseLegendList.enter().append('g');
            svgBaseLegendList.exit().remove();
            svgBaseLegendList.attr('class', 'bases').attr('transform', 'translate(' + fontSize + ', ' + fontSize + ')');

            var svgBaseItems = svgBaseLegendList.selectAll('text').data(function (d) {
                return d.map(function (series) {
                    return { title: series.Title, color: series.Color, fontSize: fontSize };
                });
            });

            svgBaseItems.enter().append('text');
            svgBaseItems.exit().remove();
            svgBaseItems.text(function (d) {
                return d.title;
            }).attr('x', function (d, i) {
                return i * d.fontSize * 2;
            }).attr('y', 0).attr('alignment-baseline', 'middle').style('font-size', function (d) {
                return d.fontSize + 'px';
            }).style('fill', function (d) {
                return d.color;
            }).on('mouseover', function (d1, i) {
                that.svg.selectAll('.linechart path').attr('class', function (d2, i) {
                    return d1.title == d2.Title ? 'highlighted' : '';
                });
                that.svgLegend.selectAll('text').style('text-decoration', function (d2, i) {
                    return d1.title == d2.title ? 'underline' : '';
                });
                $(this).css('cursor', 'pointer');
            }).on('mouseout', function (d1, i) {
                that.svg.selectAll('.linechart path').attr('class', '');
                that.svgLegend.selectAll('text').style('text-decoration', '');
                $(this).css('cursor', 'auto');
            });
            svgBaseItems.append("title").text(function (d) {
                return d.title;
            });

            var svgExtraLegendList = this.svgLegend.selectAll('g.extra').data([extraSeries]);
            svgExtraLegendList.enter().append('g');
            svgExtraLegendList.exit().remove();
            svgExtraLegendList.attr('class', 'extra').attr('transform', 'translate(' + fontSize / 2 + ', ' + fontSize * (baseSeries.length > 0 ? 2 : 1) + ')');

            var svgExtraItems = svgExtraLegendList.selectAll('text').data(function (d) {
                return d.map(function (series) {
                    return { title: series.Title, color: series.Color, fontSize: fontSize };
                });
            });

            svgExtraItems.enter().append('text');
            svgExtraItems.exit().remove();
            svgExtraItems.text(function (d) {
                return d.title;
            }).attr('x', 0).attr('y', function (d, i) {
                return i * d.fontSize;
            }).attr('alignment-baseline', 'middle').style('font-size', function (d) {
                return d.fontSize + 'px';
            }).style('fill', function (d) {
                return d.color;
            }).on('mouseover', function (d1, i) {
                that.svg.selectAll('.linechart path').attr('class', function (d2, i) {
                    return d1.title == d2.Title ? 'highlighted' : '';
                });
                that.svgLegend.selectAll('text').style('text-decoration', function (d2, i) {
                    return d1.title == d2.title ? 'underline' : '';
                });
                $(this).css('cursor', 'pointer');
            }).on('mouseout', function (d1, i) {
                that.svg.selectAll('.linechart path').attr('class', '');
                that.svgLegend.selectAll('text').style('text-decoration', '');
                $(this).css('cursor', 'auto');
            });
            svgExtraItems.append("title").text(function (d) {
                return d.title;
            });
        };

        LineChart.prototype.addLegend = function () {
            _super.prototype.addLegend.call(this);

            var legendList = this.legend.select('ul');

            for (var i = 0; i < this.responseData.Series.length; i++) {
                var currentSeries = this.responseData.Series[i];
                var title = currentSeries.Title;
                var colorClass = currentSeries.Color.toLowerCase().replace('#', '_');

                legendList.append('li').text(title + ' (' + currentSeries.Color + ')').attr('style', 'color: ' + currentSeries.Color.toLowerCase());
            }
        };

        LineChart.prototype.updateScales = function () {
            this.scaleX = this.config.scaleX();
            if (this.config.scaleX === d3.scale.log) {
                this.scaleX.clamp(true);
                var domainXMin = this.responseData.Axes.X.Min == 0 ? bs.Chart.MIN_LOG_INPUT : this.responseData.Axes.X.Min;
                this.scaleX.domain([domainXMin, this.responseData.Axes.X.Max]);
            } else {
                this.scaleX.domain([this.responseData.Axes.X.Min, this.responseData.Axes.X.Max]);
            }
            this.scaleX.range([this.config.margins.left, this.config.width - this.config.margins.right]);
            if (this.config.scaleX && this.config.scaleX.nice === true) {
                this.scaleX.nice();
            }

            this.scaleY = this.config.scaleY();
            if (this.config.scaleY === d3.scale.log) {
                this.scaleY.clamp(true);
                var domainYMin = this.responseData.Axes.Y.Min == 0 ? bs.Chart.MIN_LOG_INPUT : this.responseData.Axes.Y.Min;
                this.scaleY.domain([domainYMin, this.responseData.Axes.Y.Max]);
            } else {
                this.scaleY.domain([this.responseData.Axes.Y.Min, this.responseData.Axes.Y.Max]);
            }
            this.scaleY.range([this.config.height - this.config.margins.bottom, this.config.margins.top]);
            if (this.config.scaleY && this.config.scaleY.nice === true) {
                this.scaleY.nice();
            }
        };

        LineChart.prototype.bindMouseEvents = function () {
            var that = this;

            $(document).on('mouseenter', this.config.chartContainerSelector + ' path', function () {
                var title = $(this).attr('title');

                $(that.config.chartContainerSelector + ' div.legend ul li').each(function () {
                    var color = $(this).css('color');
                    if (color == 'white' || color == 'rgb(255, 255, 255)') {
                        color = $(this).css('background-color');
                    }
                    if (title == $(this).text()) {
                        $(this).css('background-color', color);
                        $(this).css('color', 'white');
                    } else {
                        $(this).css('background-color', 'white');
                        $(this).css('color', color);
                    }
                });

                $(that.config.chartContainerSelector + ' g.legend text').each(function () {
                    var legendText = $(this).clone().children().remove().end().text();
                    if (title == legendText) {
                        $(this).css('text-decoration', 'underline');
                    } else {
                        $(this).css('text-decoration', '');
                    }
                });
            });

            $(document).on('mouseleave', this.config.chartContainerSelector + ' path', function () {
                $(that.config.chartContainerSelector + ' div.legend ul li').each(function () {
                    var color = $(this).css('color');
                    if (color == 'white' || color == 'rgb(255, 255, 255)') {
                        color = $(this).css('background-color');
                    }
                    $(this).css('background-color', 'white');
                    $(this).css('color', color);
                });

                $(that.config.chartContainerSelector + ' g.legend text').css('text-decoration', '');
            });

            $(document).on('mouseenter', this.config.chartContainerSelector + ' div.legend ul li', function () {
                var title = $(this).text();
                var color = $(this).css('color');
                if (color == 'white' || color == 'rgb(255, 255, 255)') {
                    color = $(this).css('background-color');
                }
                $(this).css('background-color', color);
                $(this).css('color', 'white');
                $(that.config.chartContainerSelector + ' path').each(function () {
                    if (title == $(this).attr('title')) {
                        $(this).attr('class', 'highlighted');
                    } else {
                        $(this).attr('class', '');
                    }
                });
                $(this).css('cursor', 'pointer');
            });

            $(document).on('mouseleave', this.config.chartContainerSelector + ' div.legend ul li', function () {
                var color = $(this).css('color');
                if (color == 'white' || color == 'rgb(255, 255, 255)') {
                    color = $(this).css('background-color');
                }
                $(this).css('color', color);
                $(this).css('background-color', 'white');
                $(that.config.chartContainerSelector + ' path').attr('class', '');
                $(this).css('cursor', 'auto');
            });
        };
        LineChart.CHART_NAME = "LineChart";
        LineChart.SERIES_TYPE_LINE = "line";
        LineChart.SERIES_TYPE_SCATTER = "scatter";

        LineChart.DEFAULT_SCATTER_CIRCLE_RADIUS = 3.5;

        LineChart.DEFAULT_SCATTER_CIRCLE_COLOR = "#98cc78";

        LineChart.DEFAULT_LINE_VERTEX_FIELD_INDICES = {
            x: 0,
            y: 1
        };

        LineChart.DEFAULT_SCATTER_CIRCLE_FIELD_INDICES = {
            x: 0,
            tooltipTitle: 1,
            y: 2
        };

        LineChart.DEFAULT_SCATTER_FILL_GRADIENT = {
            id: "scatter-main",
            stops: [
                {
                    color: "#2c53a1",
                    offset: "0%",
                    opacity: .5
                },
                {
                    color: "#2c53a1",
                    offset: "50%",
                    opacity: 1
                }
            ],
            direction: bs.Chart.GRADIENT_DIRECTION_VERTICAL
        };
        return LineChart;
    })(bs.Chart);
    bs.LineChart = LineChart;
})(bs || (bs = {}));
