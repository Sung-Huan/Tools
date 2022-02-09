args = commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Should provide two arguments: paramFile inFile outFile", call. = FALSE)
}

paramFile <- args[1]
inFile <- args[2]
outFile <- args[3]
library(PerseusR)
parameters <- parseParameters(paramFile)
mdata <- read.perseus(inFile)
counts <- main(mdata)
mName <- colnames(counts)
counts <- as.matrix(counts)
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[, 1])
if (!is.installed('preprocessCore')) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
    BiocManager::install("preprocessCore")
}
library(preprocessCore)
qdata <- normalize.quantiles(counts)
qdata <- data.frame(qdata)
colnames(qdata) <- mName
outdata <- matrixData(main = qdata, annotCols = annotCols(mdata),
                      annotRows = annotRows(mdata))
print(paste('writing to', outFile))
write.perseus(outdata, outFile)