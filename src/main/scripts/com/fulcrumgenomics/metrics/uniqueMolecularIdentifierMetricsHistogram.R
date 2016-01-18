# Script to generate a chart of unique molecular identifiers in a SAM or BAM file.  The input data
# should be created by com.fulcrumgenomics.metrics.CollectUniqueMolecularIdentifierMetrics
#
# @author Nils Homer

# Parse the arguments
args        = commandArgs(trailing=T);
metricsFile = args[1];
outputFile  = args[2];
bamFile     = ifelse(length(args) < 3, NA, args[3]);

# Figure out where the metrics and the histogram are in the file and parse them out
startFinder = scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE);

firstBlankLine=0;

for (i in 1:length(startFinder))
{
	if (startFinder[i] == "") {
		if (firstBlankLine==0) {
			firstBlankLine = i+1;
		} else {
			secondBlankLine = i+1;
			break;
		}
	}
}

metrics   = read.table(metricsFile, header=T, nrows=1, sep="\t", skip=firstBlankLine);
histogram = read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine);

inputs = c(
		list(histogram$OBSERVED),
		list(histogram$FROM_DUPLICATE_SET),
		list(histogram$FROM_BACKGROUND),
		list(histogram$UNIFORM_RANDOM)
		);

x = list();
y = list();
labels = c("OBSERVED", "FROM_DUPLICATE_SET", "FROM_BACKGROUND", "UNIFORM_RANDOM");
stopifnot(length(labels) == length(inputs));
colors = c("black", "red", "green3", "blue", "cyan")[1:length(labels)];
ltys = 1:length(colors);
x_upper_limit = max(metrics$MEAN_READS_PER_OBSERVED_UMI + (3 * metrics$SD_READS_PER_OBSERVED_UMI), 25);

for (i in 1:length(inputs)) {
	breaks=seq(from=-1, to=max(histogram[,2:ncol(histogram)]), by=1)+0.5;
	h = hist(inputs[[i]], breaks=breaks, plot=FALSE);
	breaks = 1:length(h$counts);
	for (i in 1:length(breaks)) {
		breaks[i] = (h$breaks[i] + h$breaks[i+1]) / 2;
	}
    # ignore zero-count, so start at 2
	x = append(x, list(breaks[2:min(x_upper_limit, length(breaks))]));
	y = append(y, list(h$counts[2:min(x_upper_limit, length(breaks))]));
}

x_range = range(x);
y_range = range(y);

pdf(outputFile);

plot(x = x_range,
		y = y_range,
		xlab="Reads per UMI",
		ylab="Counts",
		main=paste("Histogram of Reads per UMI for:", bamFile, sep="\n"), 
		xlim=x_range,
		ylim=y_range,
		type="n");

for (i in 1:length(x)) {
	lines(x = x[[i]],
			y = y[[i]],
			col = colors[i],
			lwd=1.5,
			lty=ltys[i]);
	points(x = x[[i]],
			y = y[[i]],
			col = colors[i],
			pch=20);
}


legend(x="topright", legend=labels, col=colors, lwd=1.5, lty=ltys);

dev.off();
