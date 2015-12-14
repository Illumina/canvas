args = commandArgs(TRUE)

args[1] = "CanvasSomaticTestData\\BoringTumor\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\BoringTumor\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\T8N0\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\T8N0\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\T2N6\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\T2N6\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\Colo829\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\Colo829\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\Colo829Pure\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\Colo829Pure\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\Gel9FFPE\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\Gel9FFPE\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\Gel9FF\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\Gel9FF\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\Gel35FFPE\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\Gel35FFPE\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\Gel35FF\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\Gel35FF\\GaussianMixtureModel.txt"

args[1] = "CanvasSomaticTestData\\HCC2218\\CNVModeling.txt"
args[2] = "CanvasSomaticTestData\\HCC2218\\GaussianMixtureModel.txt"


modelPoints = read.table(args[1], header=FALSE, sep="\t", nrows=36)
segments2 = read.table(args[1], header=FALSE, sep="\t", skip=37)
gmm = read.table(args[2], header=FALSE, sep="\t", skip=1, nrows=36)
segments = read.table(args[2], header=FALSE, sep="\t", skip=39)

nSegments = dim(segments)[1]
nCols = dim(segments)[2]
postProbs = as.matrix(segments[, 3:nCols])
assignment = apply(postProbs, 1, which.max)

nModels = dim(gmm)[1]
col = rainbow(nModels)
plot(modelPoints$V1, modelPoints$V2, col="blue", bg="blue", pch=23, xlab="MAF", ylab="Coverage")

if (all(segments2$V8 == -1)){
	pch = assignment
}else{
	pch = as.character(segments2$V8)
}

points(segments$V1, segments$V2, col=col[assignment], pch=pch)

library(car)
for (i in 1:dim(gmm)[1]){
	mu = c(gmm[i, ]$V6, gmm[i, ]$V7)
	sigma = matrix(as.matrix(gmm[i, 8:11]), nrow=2, byrow=TRUE)
	omega = gmm[i, ]$V5
	if (omega > 0){
		ellipse(mu, sigma, 1, col="red")
		points(mu[1], mu[2], col="red", bg="red", pch=19)
	}else{
		ellipse(mu, sigma, 1, col="red", center.pch=15)
		points(mu[1], mu[2], col="green", bg="green", pch=19)
	}
}

points(modelPoints$V1, modelPoints$V2, col="blue", bg="blue", pch=23)

