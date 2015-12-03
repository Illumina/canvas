drawMAPlot = function (chartSelector, data, enlarge) {
	var rows = d3.csv.parse(data);
	
	// Create MA plot
	var maData = {
		Response: {
			Axes: {
				X: {
					Min: 0.1,
					Max: 100,
					Title: "Normalized Mean Count"
				},
				Y: {
					Min: -0.5,
					Max: 0.5,
					Title: "Log2 Fold Change"
				}
			},
			Title: "MA Plot",
			Series: [
				{
					Type: "Scatter",
					Title: "Significant",
					Color: "Red",
					PointSize: 3 * enlarge,
					Fields: {
						tooltipTitle:0,
						x: 1,
						y: 2
					},
					Data: [],
					Accessors: {
						"powertip": function (d){ return d[0] + ": (" + d[1].toPrecision(3) + ", " + d[2].toPrecision(3) + ")"; }
					}
				},
				{
					Type: "Scatter",
					Title: "Insignificant",
					Color: "Green",
					PointSize: 3 * enlarge,
					Fields: {
						tooltipTitle:0,
						x: 1,
						y: 2
					},
					Data: [],
					Accessors: {
						"powertip": function (d){ return d[0] + ": (" + d[1].toPrecision(3) + ", " + d[2].toPrecision(3) + ")"; }
					}
				}
			]
		}
	};
	
	for (var i = 0; i < rows.length; i++){
		var row = rows[i];
		var marker = row[""];
		var x = +row["baseMean"], y = +row["log2FoldChange"];
		if (isNaN(x) || isNaN(y)){ continue; }
		if (x < maData.Response.Axes.X.Min){ continue; }
		var pAdj = +row["padj"];
		maData.Response.Axes.X.Max = Math.max(maData.Response.Axes.X.Max, x);
		maData.Response.Axes.Y.Min = Math.min(maData.Response.Axes.Y.Min, y);
		maData.Response.Axes.Y.Max = Math.max(maData.Response.Axes.Y.Max, y);
		if (!isNaN(x) && pAdj <= 0.05){
			maData.Response.Series[0].Data.push([marker, x, y]);
		}else{
			maData.Response.Series[1].Data.push([marker, x, y]);
		}
	}
	
	var yMax = Math.max(Math.abs(maData.Response.Axes.Y.Min), maData.Response.Axes.Y.Max);
	maData.Response.Axes.Y.Min = -yMax;
	maData.Response.Axes.Y.Max = yMax;
	
	var tickValues = [];
	var tickValue;
	for (tickValue = 0.1; tickValue < maData.Response.Axes.X.Max; tickValue *= 10){
		tickValues.push(tickValue);
	}
	tickValues.push(tickValue);
	maData.Response.Axes.X.Max = tickValue;
	
	var maxSeriesTitleLength = _.max(maData.Response.Series.map(function (s){ return s.Title.length; }));
	var fontSize = 10 * enlarge;
	var width = 550 * enlarge;
	var height = 350 * enlarge;
	var margins = {
					top: 30 * enlarge,
					left: 30 * enlarge,
					right: fontSize * maxSeriesTitleLength,
					bottom: 50 * enlarge
	};
	$(chartSelector).width(width);
	$(chartSelector).height(height);
	bs.Chart.MIN_LOG_INPUT = 0.1;
	var maChart = new bs.LineChart(maData.Response, {
		chartContainerSelector: chartSelector,
		width: $(chartSelector).width(),
		height: $(chartSelector).height(),
		preserveAspectRatio: 'none',
		margins: margins,
		scaleX: (function () {
			var scaleX = d3.scale.log;
			scaleX.nice = false;
			return scaleX;
		})(),
		xAxis:{
			tickValues: tickValues
		},
		legend: true,
		svgLegend: true,
		fontSize: fontSize
	});
	maChart.svg.selectAll("text").style("font-size", fontSize + "px");
	var legendRect = maChart.svgLegend.select("rect");
	var legendItems = maChart.svgLegend.selectAll("text");
	var legendWidth = Math.min(margins.right - 5, _.max(legendItems[0].map(function (r) { return r.getBBox().width; })) + fontSize);
	legendRect.attr("width", legendWidth);
	
	return maChart;
}