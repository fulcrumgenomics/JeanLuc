# Script to generate a chart of quality score distribution in a file
# @author nhomer

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
outputFile   <- args[2]
bamFile  <- ifelse(length(args) < 3, NA, args[3])
subtitle <- ifelse(length(args) < 4, "", args[4])

# Figure out where the metrics and the histogram are in the file and parse them out
startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

firstBlankLine=0

for (i in 1:length(startFinder))
{
        if (startFinder[i] == "") {
                if (firstBlankLine==0) {
                        firstBlankLine=i+1
                } else {
                        secondBlankLine=i+1
                        break
                }
        }
}

metrics <- read.table(metricsFile, header=T, nrows=1, sep="\t", skip=firstBlankLine)
histogram <- read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine)

histogram$coverage = histogram$coverage[!is.na(histogram$count)]
histogram$count = histogram$count[!is.na(histogram$count)]

# Then plot the histogram as a PDF
pdf(outputFile)

meanCoverage <- metrics$MEAN_COVERAGE;
percentOfMean <- histogram$coverage / meanCoverage; # x-axis
percentCovered <- rep(0, length(histogram$count)); # y-axis
# must do a cumulative sume of percentCovered
totalCount = sum(as.numeric(histogram$count));
for (i in 1:length(percentCovered)) {
    percentCovered[i] = 100.0 * sum(as.numeric(histogram$count[i:length(percentCovered)])) / totalCount;
}

ymin = percentCovered[round(meanCoverage+1)];
ymax = min(100,max(percentCovered))

plot(x=percentOfMean,
     y=percentCovered,
     xlim=c(0, 1.0),
     ylim=c(ymin, ymax),
     type="l",
     main=paste("WGS Base Coverage Plot", ifelse(is.na(bamFile),"",paste("\nin file",bamFile))," ",ifelse(subtitle == "","",paste("(",subtitle,")",sep="")),sep=""),
     xlab="Fold Coverage of Mean",
     ylab="% of Bases Covered",
     col="blue",
     lwd=5);
dev.off()

