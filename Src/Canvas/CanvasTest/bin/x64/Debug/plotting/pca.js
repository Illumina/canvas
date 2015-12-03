// chartSelector: selector of the div container
// data: string containing content of pca.csv
pca = function (chartSelector, data) {
    var rows = d3.csv.parseRows(data);
    if (rows.length < 4) { return; }
    var sampleIDs = rows[0].slice(1);
    var sampleGroups = rows[1].slice(1);
    var pc1 = rows[2].slice(1).map(function (x) { return +x; });
    var pc2 = rows[3].slice(1).map(function (x) { return +x; });

    var scatterData = {
        Axes: {
            X: {
                Min: Infinity,
                Max: -Infinity,
                Title: rows[2][0]
            },
            Y: {
                Min: Infinity,
                Max: -Infinity,
                Title: rows[3][0]
            }
        },
        Title: "PCA",
        Series: [
        ]
    };

    var maxSampleGroupLength = 0;
    var sampleGroup2Series = {};
    sampleGroups.forEach(function (sampleGroup) {
        if (!(sampleGroup in sampleGroup2Series)) {
            sampleGroup2Series[sampleGroup] = {
                Type: "Scatter",
                Title: sampleGroup,
                Color: "",
                PointSize: 5,
                Fields: { x: 0, y: 1, tooltipTitle: 2 },
                Accessors: {
                    powertip: function (d) {
                        return "(" + d[0] + ", " + d[1] + "): " + d[2];
                    }
                },
                Data: []
            };
            scatterData.Series.push(sampleGroup2Series[sampleGroup]);
            maxSampleGroupLength = Math.max(maxSampleGroupLength, sampleGroup.length);
        }
    });

    var nSeries = scatterData.Series.length;
    var colors = []; (new KolorWheel("#0000FF")).abs(0, -1, -1, nSeries).each(function () { colors.push(this.getHex()); });
    scatterData.Series.forEach(function (series, i) {
        series.Color = colors[i];
    });

    sampleIDs.forEach(function (sampleID, i) {
        var sampleGroup = sampleGroups[i];
        var x = pc1[i];
        var y = pc2[i];
        scatterData.Axes.X.Min = Math.min(scatterData.Axes.X.Min, x);
        scatterData.Axes.X.Max = Math.max(scatterData.Axes.X.Max, x);
        scatterData.Axes.Y.Min = Math.min(scatterData.Axes.Y.Min, y);
        scatterData.Axes.Y.Max = Math.max(scatterData.Axes.Y.Max, y);
        sampleGroup2Series[sampleGroup].Data.push([x, y, sampleID]);
    });

    var log2Ceil = function (x){
        var sign = x >= 0 ? 1 : -1;
        x = Math.abs(x);
        var log2X = Math.log(x) / Math.LN2;
        if (sign > 0) {
            log2X = Math.ceil(log2X);
        } else {
            log2X = Math.floor(log2X);
        }
        return sign * Math.pow(2, log2X);
    }

    var log2Floor = function (x) {
        var sign = x >= 0 ? 1 : -1;
        x = Math.abs(x);
        var log2X = Math.log(x) / Math.LN2;
        if (sign > 0) {
            log2X = Math.floor(log2X);
        } else {
            log2X = Math.ceil(log2X);
        }
        return sign * Math.pow(2, log2X);
    }

    scatterData.Axes.X.Min = log2Floor(scatterData.Axes.X.Min);
    scatterData.Axes.X.Max = log2Ceil(scatterData.Axes.X.Max);
    scatterData.Axes.Y.Min = log2Floor(scatterData.Axes.Y.Min);
    scatterData.Axes.Y.Max = log2Ceil(scatterData.Axes.Y.Max);

    $(chartSelector).width($(chartSelector).width() + 10 * maxSampleGroupLength - 30);
    var config = {
        chartContainerSelector: chartSelector,
        width: $(chartSelector).width(),
        height: $(chartSelector).height(),
        margins: {
            top: 30,
            left: 50,
            right: 10 * maxSampleGroupLength,
            bottom: 50
        },
        svgLegend: true,
        scaleX: d3.scale.linear,
        scaleY: d3.scale.linear
    };
    //config.scaleX.nice = true;
    //config.scaleY.nice = true;
    var chart = new bs.LineChart(scatterData, config);
    chart.svg.selectAll(".tick.major text").style("font-size", "10px");
    return chart;
}