## Visualizes values of the best Canvas purity model predictions ## 

#First read in the arguments listed at the command line
args=(commandArgs(TRUE))

if (length(args)==0)
{
    print("No arguments supplied.")
    print("Usage: plotSomaticModelMatrix.r [Canvas TempCNV* directory] [output directory]")
} else if (length(args)!=2)
{
    print("Incorrect number of arguments.")
    print("Usage: plotSomaticModelMatrix.r [Canvas TempCNV* directory] [output directory]")

} else 
{
    basedir = args[1]
    outdir = args[2]
    path = paste(basedir,"/CNVModeling.txt",sep="")
    if (file.exists(path))
    {
        # main plotting function
        print("Ceating Canvas somatic model plot.")
        no_col = count.fields(path, sep = "\t")
        transition = which(diff(no_col)!=0)
        observed_model = read.table(path, skip=transition+1)
        expected_model = read.table(path, skip=1, nrow=transition)
        outfile = paste(outdir, "CanvasSomaticModel.png", sep="")
        png(outfile, width = 680, height = 680)
        plot(observed_model[,1], observed_model[,2], xlim=c(0,0.5), xlab="MAF", ylab="Coverage", pch=19, col = "gray", cex=1.3,main="Canvas somatic model", cex.axis=1.4)
        points(expected_model[,1], expected_model[,2], xlim=c(0,0.5), xlab="MAF", ylab="Coverage", pch=19, col="red", cex=1.7, cex.axis=1.4)
        text(expected_model[,1], expected_model[,2], labels=expected_model[,3], cex= 1.4, offset=10)
        legend("bottomleft", legend=c("Observed Coverage/MAFs", "Expected + CN GTs"),col=c("gray","red"),cex=1.2,bty="n",pch=19, title="Legend")
        dev.off()
    }
    else{
        print("File CNVModeling does not exist. Check that Canvas TempCNV* directory is specified correctly.")
	}
}
