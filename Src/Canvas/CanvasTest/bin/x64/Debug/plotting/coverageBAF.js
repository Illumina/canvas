// Functions for SV plotting
function removeVCFHeader(vcfFileAsString) {
    return vcfFileAsString.substr((vcfFileAsString.indexOf("#CHROM")));
}

function parseVCFRow(row) {
    if (!(row[0] + row[1] == "##")) {
        for (var key in row) {
            if (key === "INFO") {
                row[key] = row[key].split(";");
                s = "{ "
                for (var i = 0; i < row[key].length; i++) {
                    if (i > 0) {
                        s = s + ','
                    }
                    if (row[key][i].indexOf("=") != -1) {
                        s = s + "\"" + row[key][i].split("=")[0] + "\" :  \"" + row[key][i].split("=")[1] + "\"";
                    }
                    else {
                        s = s + "\"" + row[key][i] + "\" :  \"" + row[key][i] + "\"";
                    }
                }
                row[key] = JSON.parse(s + " }");
            }
        }
    }
}

function isPassingSV(parsedRow) {
    if (parsedRow["FILTER"] == "PASS") {
        return true;
    }
    return false;
}

function isSV(parsedRow) {
    switch (parsedRow["INFO"].SVTYPE) {
        case "BND":
            return true;
            break;
        case "DEL":
            return true;
            break;
        case "DUP":
            return true;
            break;
        case "INS":
            return true;
            break;
        case "INV":
            return true;
            break;
        default:
            return false;
            break;
    }
}

function getSVType(parsedRow) {
    return parsedRow["INFO"].SVTYPE;
}

function findMate(variant, SVlist) {
    for (var i = 0; i < SVlist.length; i++) {
        if (SVlist[i].INFO.MATEID === undefined) {
            continue;
        }
        if (variant.ID.slice(0, -1) === SVlist[i].INFO.MATEID.slice(0, -1) && variant.ID != SVlist[i].ID) {
            return SVlist[i];
        }
    }
}

function getParabolaPoints(root1, root2, max, numPoints) {
    var k = -max / (Math.pow((root2 - root1) / 2, 2));
    var points = []
    for (var i = 0; i <= numPoints; i++) {
        var x = Math.min(root1, root2) + i * Math.abs(root2 - root1) / numPoints;
        var y = k * (x - root1) * (x - root2);
        points.push([x, y]);
    }
    return points;
}

// parsing logic taken from TumorNormalReportGenerator.AddCNVChart
function parseRow(row) {
    for (var key in row) {
        if (key == "Chromosome") { continue; }
        row[key] = +row[key];
    }
    row["End"] -= 1; // 0-based, inclusive

    var sum = 0;
    for (var i = 0; true; i++){
        var key = "VariantFrequencyBin" + i;
        if (!(key in row)) { break; }
        sum += row[key];
    }

    row["VariantFrequencyBinSum"] = sum;
}

var CNVType = {
    Reference: "Reference",
    LOH: "LOH",
    Loss: "Loss",
    Gain: "Gain"
};
function getCNVType(row) {
    var type = CNVType.Reference;

    var LOH = row["CopyNumber"] == row["MajorChromosomeCount"];
    var refCN = isNaN(row["ReferencePloidy"]) ? 2 : row["ReferencePloidy"];
    if (row["CopyNumber"] < refCN) {
        type = CNVType.Loss;
    } else if (row["CopyNumber"] > refCN) {
        type = CNVType.Gain;
    } else if (LOH && row["CopyNumber"] >= 2) {
        type = CNVType.LOH;
    }
    return type;
}

function getChromosomeStartEnd(rows) {
    var chroms = [];
    var chromosomeLength = {};
    for (var i = 0; i < rows.length; i++) {
        var row = rows[i];
        var chrom = row["Chromosome"];
        if (!(chrom in chromosomeLength)) {
            chromosomeLength[chrom] = 0;
            chroms.push(chrom);
        }
        chromosomeLength[chrom] = Math.max(row["End"] + 1, chromosomeLength[chrom]);
    }

    var chromosomeStartEnd = {};
    var start = 0;
    for (var i = 0; i < chroms.length; i++) {
        var chrom = chroms[i];
        chromosomeStartEnd[chrom] = { Start: start, End: 0 }; // both inclusive
        start += chromosomeLength[chrom];
        chromosomeStartEnd[chrom].End = start - 1;
    }

    return chromosomeStartEnd;
}

function reduceLongChromosomes(chromosomeStartEnd){
	var maxWidthRatio = 2.5;
	var maxLength = 0;
	var minLength = Infinity;
	var chromStart = [];
	for (chrom in chromosomeStartEnd){
		var length = chromosomeStartEnd[chrom].End - chromosomeStartEnd[chrom].Start + 1;
		maxLength = Math.max(maxLength, length);
		minLength = Math.min(minLength, length);
		chromStart.push([chrom, chromosomeStartEnd[chrom].Start]);
	}
	
	var reductionFactor = 0;
	if (maxLength / minLength > maxWidthRatio){
		var newMaxLength = Math.round(maxWidthRatio * minLength);
		reductionFactor = (maxLength - newMaxLength) / (maxLength - minLength);
	}
	
	chromStart.sort(function (a, b) { return a[1] - b[1]; });
	var chromosomeCorrectionFactors = {};
	var prevEnd = -1;
	chromStart.forEach(function (x){
		var chrom = x[0];
		var oldLength = chromosomeStartEnd[chrom].End - chromosomeStartEnd[chrom].Start + 1;
		var newLength = oldLength - (oldLength - minLength) * reductionFactor;
		chromosomeCorrectionFactors[chrom] = newLength / oldLength;
		chromosomeStartEnd[chrom].Start = prevEnd + 1; // inclusive
		chromosomeStartEnd[chrom].End = prevEnd + newLength; // inclusive; = chromosomeStartEnd[chrom].Start + newLength - 1
		prevEnd = chromosomeStartEnd[chrom].End;
	});
	
	return chromosomeCorrectionFactors;
}

drawCoverageBAF = function (chartSelector, dataStr, enlarge, somaticSVVcfStr) {
    var vcfRows;
    var SVlist;
    var lineSeries;
    if (somaticSVVcfStr != "") {
        vcfRows = d3.tsv.parse(removeVCFHeader(somaticSVVcfStr).substr(1));
        SVlist = [];
        vcfRows.forEach(function (row) {
            parseVCFRow(row);
            if (isSV(row) && isPassingSV(row)) {
                SVlist.push(row)
            }
        });

        lineSeries = []
    }

    var haveSVs = SVlist.length !== 0;
    var haveCanvasData = dataStr !== "";
    if (!haveSVs && !haveCanvasData) { return; }

	var rows = d3.tsv.parse(dataStr.substr(1));
	rows.forEach(function (row) { parseRow(row); });
	var chromStartEnd = getChromosomeStartEnd(rows);
	
	var chromCorrectionFactors = reduceLongChromosomes(chromStartEnd);
	
	var scatterSeries = [
		{
		    "Title": CNVType.Gain,
		    "Color": "rgb(213,94,0)" // vermillion
		},
		{
			"Title": CNVType.Loss,
			"Color": "rgb(0,158,115)" //bluish green
		},
		{
			"Title": CNVType.LOH,
			"Color": "rgb(0,114,178)" // blue
		},
		{
			"Title": CNVType.Reference,
			"Color": "Black"
		}
	];
	for (var i = 0; i < scatterSeries.length; i++) {
		scatterSeries[i].Type = "Scatter";
		scatterSeries[i].CircleRadius = 1.5 * enlarge;
		scatterSeries[i].CircleStrokeWidth = 0;
		scatterSeries[i].Fields = { "x": 0, "y": 1, "tooltipTitle": 1 };
		scatterSeries[i].Accessors = {
			"powertip": function (d) { return d[1]; }
		};
		scatterSeries[i].Data = [];
	}

	var scatterSeriesByType = {
		Gain: scatterSeries[0],
		Loss: scatterSeries[1],
		LOH: scatterSeries[2],
		Reference: scatterSeries[3]
	};

	var heatSeries = [
		{
			"Type": "Heat",
			"Title": "",
			"Color": "Blue",
			"Fields": { "x": 0, "y": 1 },
			"Data": []
		}
	];

	var data = {
		Response: {
			"Title": "Coverage and Minor Allele Frequency",
			"StackedAxes": {
				"X": {
					"Min": 0,
					"Max": 0,
					"Title": "Genomic Location"
				},
				"Y": [
					
				]
			},
			"StackedSeries": [],
			"ColorScaleTitle": ["Percentage", null],
			"ColorScaleDomain": [[0, 30], null],
			"ColorScaleRange": [["white", "blue"], null],
			"Series": null // dummy
		}
	};

	if (haveCanvasData) {
	    data.Response.StackedAxes.Y.push({
	        "Min": 0,
	        "Max": 1,
	        "Title": "B Allele Frequency"
	    });
	    data.Response.StackedAxes.Y.push({
	        "Min": 0,
	        "Max": 0,
	        "Title": "Normalized Coverage"
	    });
	    data.Response.StackedSeries.push(heatSeries);
	    data.Response.StackedSeries.push(scatterSeries);
	}
    var SVPlotYAxisTitle;

    if (SVlist[0]) {
        if (SVlist[0].INFO.SOMATICSCORE) {
        SVPlotYAxisTitle = "SV Somatic Score";
		}
		else {
            SVPlotYAxisTitle = "SV Qual Score";
        }
    }
    

	if (haveSVs) {
	    data.Response.StackedAxes.Y.push({
	        "Min": 0,
	        "Max": 0,
	        "Title": SVPlotYAxisTitle
	    });
        
	    data.Response.StackedSeries.push(lineSeries);
	    for (var i = 0; i < SVlist.length; i++) {
	        variant = SVlist[i];
	        var SVTYPE = getSVType(variant);
	        var svLabel;
	        if (variant.CHROM.toLowerCase() === "mt" || variant.CHROM.toLowerCase() === "chrm") {
	            continue;
	        }
	        var start = Number(chromStartEnd[variant.CHROM].Start) + Number(variant.POS);
	        var end;
	        var svScore;
	        if (variant.INFO.SOMATICSCORE) {
	            svScore = Number(variant.INFO.SOMATICSCORE);
	        }
	        else {
	            svScore = Number(variant.QUAL);
	        }
	        
	        if (SVTYPE === "INS") {
	            continue;
	        }
	        else if (SVTYPE === "BND") {
	            mate = findMate(variant, SVlist);
                // A mate might not exist if it was on a decoy contig and filtered out before coverageBAF.html/js got it
	            if (!mate) { continue; }

	            if (mate.CHROM.toLowerCase() === "mt" || mate.CHROM.toLowerCase() === "chrm") {
	                continue;
	            }
	            end = Number(chromStartEnd[mate.CHROM].Start) + Number(mate.POS);
	            color = "rgb(0,114,178)"; // blue
	            svLabel = "Translocation";
	        }
	        else {
	            end = Number(start) + Number(variant.INFO.SVLEN);
	            if (SVTYPE === "DUP") {
	                color = "red";
	                svLabel = "Duplication";
	            }
	            else if (SVTYPE === "INV") {
	                color = "rgb(213,94,0)"; // vermillion
	                svLabel = "Inversion";
	            }
	            else if (SVTYPE === "DEL") {
	                color = "rgb(0,158,115)"; //bluish green
	                svLabel = "Deletion";
	            }
	        }
	        // only work with SVs greater than 100kbp in length per 
	        if (Math.abs(start - end) < 100000) {
	            continue;
	        }
	        parabolaPoints = getParabolaPoints(start, end, svScore, 500);
	        // lineSeriesByType[SVTYPE].Data = lineSeriesByType[SVTYPE].Data.concat(parabolaPoints);
	        data.Response.StackedAxes.Y[2].Max = Math.max(data.Response.StackedAxes.Y[2].Max, svScore);

	        data.Response.StackedSeries[2].push({
	            "Title": svLabel,
	            "Color": color,
	            "Type": "Line",
	            "Data": parabolaPoints,
	            "Fields": { "x": 0, "y": 1, "tooltipTitle": 1 },
	        });
	    }
	}

	var prevRow = null;
	var k = 0; // row count for the current chromosome
	var prevVariantFrequencyBinSum = 0;
	var maxPercentage = 20;
	for (var i = 0; i < rows.length; i++) {
		var row = rows[i];
		if (isNaN(row["NormalizedCoverage"])) { continue; }
		if (prevRow === null || row["Chromosome"] != prevRow["Chromosome"]) { k = 0; }
		var cnvType = getCNVType(row);
		var start = chromStartEnd[row["Chromosome"]].Start + row["Start"] * chromCorrectionFactors[row["Chromosome"]];
		var end = chromStartEnd[row["Chromosome"]].Start + row["End"] * chromCorrectionFactors[row["Chromosome"]];
		var mid = (start + end) / 2.0;
		var coverage = row["NormalizedCoverage"];
		if (coverage >= 14) { coverage = 14; }
		scatterSeriesByType[cnvType].Data.push([mid, coverage]);
		data.Response.StackedAxes.Y[1].Max = Math.max(data.Response.StackedAxes.Y[1].Max, coverage);

		if (isNaN(row["VariantFrequencyBinSum"]) || row["VariantFrequencyBinSum"] == 0) { continue; }
		var percentage = 0;
		var percentages = [];
		for (var j = 0; true; j++){
			var key = "VariantFrequencyBin" + j;
			if (!(key in row)) { break; }
			percentage += row[key] / row["VariantFrequencyBinSum"] * 100;
			if (j % 4 == 3) { // merge 100 bins into 25 bins
				percentages.push(percentage);
				percentage = 0;
			}
		}

		if (heatSeries[0].Data.length > 0) {
		    var prevData = heatSeries[0].Data[heatSeries[0].Data.length - 1];
		    if (start - prevData[0][1] > 1) { // gap: reset k
		        k = 0;
		        prevVariantFrequencyBinSum = 0;
		    }
		}
		if (k % 4 > 0) { // merge upto 4 segments
		    var prevData = heatSeries[0].Data.pop();
		    start = prevData[0][0];
		    percentages = percentages.map(function (p, j) {
		        return (p * row["VariantFrequencyBinSum"] + prevData[1][j] * prevVariantFrequencyBinSum) / (row["VariantFrequencyBinSum"] + prevVariantFrequencyBinSum);
		    });
		} else {
		    prevVariantFrequencyBinSum = 0;
		}
		heatSeries[0].Data.push([[start, end], percentages]);
		maxPercentage = Math.max(maxPercentage, _.max(percentages));

		prevRow = row;
		prevVariantFrequencyBinSum += row["VariantFrequencyBinSum"];
		k++;
	}

	data.Response.ColorScaleDomain[0][1] = maxPercentage;
	var chromMids = [];
	var chromEnds = [];
	for (var chrom in chromStartEnd) {
		chromMids.push([chrom, (chromStartEnd[chrom].Start + chromStartEnd[chrom].End) / 2]);
		chromEnds.push(chromStartEnd[chrom].End);
	}
	chromMids.sort(function (a, b) { return a[1] - b[1]; });
	chromEnds.sort(function (a, b) { return a - b; });
	data.Response.StackedAxes.X.Max = chromEnds[chromEnds.length - 1];
	var xTickLabelPositions = chromMids.map(function (x) { return x[1]; });
	var xTickLabels = chromMids.map(function (x) { return x[0]; });
	data.Response.StackedAxes.X.GroupLabels = xTickLabels;
	data.Response.StackedAxes.X.GroupLabelPositions = xTickLabelPositions;


    var format = d3.format();
    function formatCoverageAxisWithPlus(d) {
        s = format(d);
        return d === 14 ? s + "+" : s;
    }

    var coverageFormat = { tickFormat: formatCoverageAxisWithPlus}



	var chart = new bs.CoverageBAFChart(data.Response, {
	    chartContainerSelector: chartSelector,
	    width: $(chartSelector).width(),
	    height: $(chartSelector).height(),
	    preserveAspectRatio: 'none',
	    fontSize: 14 * enlarge,
	    margins: {
	        top: 30 * enlarge,
	        left: 60 * enlarge,
	        right: 120 * enlarge,
	        bottom: 60 * enlarge
	    },
	    xAxis: {
	        tickValues: chromEnds
	    },
        yAxes: [{ }, coverageFormat, { }],
		interChartYPadding: 10 * enlarge
	});
	return chart;
}
