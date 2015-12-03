var bs;
(function (bs) {
    (function (Charts) {
        (function (_Data) {
            _Data.AMPLICON_COVERAGE_Y_COLUMN = "coverage";
            _Data.AMPLICON_COVERAGE_TOOLTIP_TITLE = "region";
            _Data.LOW_COVERAGE_MULTIPLIER = 0.2;
            _Data.MOVING_AVG_PERIOD = 5;

            _Data.BELOW_THRESHOLD_GRADIENT = {
                direction: bs.Chart.GRADIENT_DIRECTION_VERTICAL,
                id: "below-threshold",
                stops: [
                    {
                        color: "red",
                        offset: "0%",
                        opacity: .5
                    },
                    {
                        color: "red",
                        offset: "50%",
                        opacity: 1
                    }
                ]
            };

            function parseAmpliconCoverageRegionCsv(csvString, title, xLabel, yLabel, lookupIndices, addLowCoverageSeries, addMovingAverage) {
                if (typeof lookupIndices === "undefined") { lookupIndices = bs.LineChart.DEFAULT_SCATTER_CIRCLE_FIELD_INDICES; }
                if (typeof addLowCoverageSeries === "undefined") { addLowCoverageSeries = true; }
                if (typeof addMovingAverage === "undefined") { addMovingAverage = true; }
                var arr = d3.csv.parseRows(csvString, function (row, i) {
                    return [i].concat(row);
                });

                arr.shift();

                var chartResponse = {
                    Axes: {
                        X: {
                            Min: 1,
                            Max: arr.length,
                            Title: xLabel
                        },
                        Y: {
                            Min: 0,
                            Max: d3.max(arr, function (d) {
                                return parseInt(d[lookupIndices.y]);
                            }),
                            Title: yLabel
                        }
                    },
                    Title: title,
                    Series: [
                        {
                            Type: "Scatter",
                            Title: "Title",
                            Color: "green",
                            Fields: lookupIndices,
                            Data: arr,
                            Accessors: {
                                powertip: function (d) {
                                    return _.template(bs.Chart.TOOLTIP_TEMPLATE, { content: "<p><strong>Region</strong>: " + d[lookupIndices.tooltipTitle] + "</p>" + "<p><strong>Coverage</strong>: " + d[lookupIndices.y] + "</p>" });
                                }
                            }
                        }
                    ]
                };

                if (addLowCoverageSeries) {
                    var lowCoverageThreshold = getLowCoverageThreshold(chartResponse.Series[0]);
                    var lowCoverageThresholdSeries = getLowCoverageThresholdSeries(chartResponse);
                    chartResponse.Series[0].Accessors["fill"] = function (d) {
                        if (parseInt(d[lookupIndices.y]) < lowCoverageThreshold) {
                            return "url(#" + _Data.BELOW_THRESHOLD_GRADIENT.id + ")";
                        } else {
                            return "url(#" + bs.LineChart.DEFAULT_SCATTER_FILL_GRADIENT.id + ")";
                        }
                    };
                    chartResponse.Series.push(lowCoverageThresholdSeries);
                }

                if (addMovingAverage) {
                    var movingAverageSeries = getMovingAverageSeries(chartResponse);
                    chartResponse.Series.push(movingAverageSeries);
                }

                return chartResponse;
            }
            _Data.parseAmpliconCoverageRegionCsv = parseAmpliconCoverageRegionCsv;

            function getLowCoverageThreshold(series) {
                var meanCoverage = d3.mean(series.Data, function (d) {
                    return parseInt(d[series.Fields.y]);
                }), multiplier = _Data.LOW_COVERAGE_MULTIPLIER;

                return meanCoverage * multiplier;
            }
            _Data.getLowCoverageThreshold = getLowCoverageThreshold;

            function getLowCoverageThresholdSeries(chartResponse) {
                var lowCoverageThreshold = Math.max(bs.Chart.MIN_LOG_INPUT, getLowCoverageThreshold(chartResponse.Series[0])), lookupIndices = bs.LineChart.DEFAULT_LINE_VERTEX_FIELD_INDICES, lowCoverageSeriesDataFirst = [], lowCoverageSeriesDataLast = [];

                lowCoverageSeriesDataFirst[lookupIndices.x] = chartResponse.Axes.X.Min;
                lowCoverageSeriesDataFirst[lookupIndices.y] = lowCoverageThreshold;
                lowCoverageSeriesDataLast[lookupIndices.x] = chartResponse.Axes.X.Max;
                lowCoverageSeriesDataLast[lookupIndices.y] = lowCoverageThreshold;

                var lowCoverageSeriesData = [lowCoverageSeriesDataFirst, lowCoverageSeriesDataLast];

                return {
                    Type: "Line",
                    Title: "Low Coverage Threshold",
                    Color: "red",
                    Fields: lookupIndices,
                    Data: lowCoverageSeriesData,
                    Accessors: {
                        powertip: function (d) {
                            return _.template(bs.Chart.TOOLTIP_TEMPLATE, { content: "<p>Low Coverage Threshold (" + lowCoverageThreshold.toFixed(2) + ")</p>" });
                        }
                    }
                };
            }
            _Data.getLowCoverageThresholdSeries = getLowCoverageThresholdSeries;

            function getMovingAverageSeries(chartResponse) {
                var period = _Data.MOVING_AVG_PERIOD, averager = new bs.MathUtils.MovingAverager(period), averageSeriesData = [], scatterLookupIndices = chartResponse.Series[0].Fields, lineLookupIndices = bs.LineChart.DEFAULT_LINE_VERTEX_FIELD_INDICES;

                _.each(chartResponse.Series[0].Data, function (d, index) {
                    var averagePoint = [];
                    averagePoint[lineLookupIndices.x] = index;
                    averagePoint[lineLookupIndices.y] = Math.max(bs.Chart.MIN_LOG_INPUT, averager.add(parseInt(d[scatterLookupIndices.y])));
                    averageSeriesData.push(averagePoint);
                });

                return {
                    Type: "Line",
                    Title: "Moving Average",
                    Color: "orange",
                    Fields: lineLookupIndices,
                    Data: averageSeriesData,
                    Accessors: {
                        powertip: function (d) {
                            return _.template(bs.Chart.TOOLTIP_TEMPLATE, { content: "<p>Moving Average</p>" });
                        }
                    }
                };
            }
            _Data.getMovingAverageSeries = getMovingAverageSeries;
        })(Charts.Data || (Charts.Data = {}));
        var Data = Charts.Data;
    })(bs.Charts || (bs.Charts = {}));
    var Charts = bs.Charts;
})(bs || (bs = {}));
