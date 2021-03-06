#!/usr/bin/env Rscript
# **Authors** Marie Lefebvre  INRA - UMR 1332 BFP - marie.lefebvre|at|bordeaux.inra.fr

"Parser of Rnmr1D  - Galaxy parser

Usage:
  Rnmr1D_parser.R [options]

Options:
  -h, --help           	    Show this screen.
  -z, --zip <file>     	    the full path name of the ZIP file (raw.zip)
  -s, --samples <file>      the full path name of the Sample file (tabular format)
  -p, --proccmd <file>      the full path name of the Macro-commands file for processing (text format)
  -b, --bucfile <file>      the full path name of the file of bucket's zones (tabular format)
  -c, --cpu <n>             the number of cores [default: 4]
  -o, --outdir <path>       the full path name of the directory to output the resulting files
  -l, --logfile <file>      the full path name of the LOG file [default: stderr]
" -> doc

options(warn=-1)
Write.LOG <- function(..., logfile=stderr()) cat(sprintf(...), sep='', file=logfile, append=TRUE)

library(docopt)

# Retrieve the command-line arguments
opts <- docopt(doc)


# -------------------------------------------- INPUT -------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# ----- Get & Check Inputs for Rnmr1D

# Zip file
if ( is.null(opts$zip)) stop("Parser : The full path name of the ZIP file (-z option) must be specified\n")
if ( ! file.exists(opts$zip)) stop("Parser : File ", opts$zip, " does not exist\n")
zipFile <- opts$zip

# Sample file
if ( is.null(opts$samples)) stop("Parser : The full path name of the samples file (-s option) must be specified\n")
if ( ! file.exists(opts$samples)) stop("Parser : File ", opts$samples, " does not exist\n")
samplesFile <- opts$samples

# Bucket file
if ( is.null(opts$bucfile)) stop("Parser : The full path name of the bucket file file (-b option) must be specified\n")
if ( ! file.exists(opts$bucfile)) stop("Parser : File ", opts$bucfile, " does not exist\n")
bucketFile <- opts$bucfile

# Macro command file
if ( is.null(opts$proccmd)) stop("Parser : The full path name of the macro-commands file (-p option) must be specified\n")
if ( ! file.exists(opts$proccmd)) stop("Parser : File ", opts$proccmd, " does not exist\n")
macroFile <- opts$proccmd

# Output directory
if ( is.null(opts$outdir)) stop("Parser : The directory to output the resulting files (-o option) must be specified\n")
if ( ! file.exists(opts$outdir)) stop("Parser : Directory ", opts$outdir, " does not exist\n")
outdir <- opts$outdir

# The Supevisor Log file (optional)
LOGFILE <- stderr()
if (!is.null(opts$logfile) && file.exists(opts$logfile)) LOGFILE <- opts$logfile


# -------------------------------------------- MAIN --------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

# ----- Run Rnmr1D
cmd <- paste("Rnmr1D -z ",zipFile," -s ",samplesFile," -p ",macroFile," -b ",bucketFile," -o ",outdir," -l ",LOGFILE)
system(cmd, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
Write.LOG("Parser: ------------------------------------\n", logfile=LOGFILE)
Write.LOG("Parser: Process parser\n", logfile=LOGFILE)
Write.LOG("Parser: ------------------------------------\n", logfile=LOGFILE)

# ----- Get & Check Inputs for Wrapper

# Data matrix
dataMatrix <- file.path(outdir,"data_matrix.txt") # output of Rnmr1D
if ( is.null(dataMatrix)) stop("Parser: The full path name of the data_matrix file must be specified\n")
if ( ! file.exists(dataMatrix)) stop("Parser: File ", dataMatrix, " does not exist\n")
dataMatrixFile <- dataMatrix

# Samples
samplesFile <- file.path(outdir,"samples.csv") # output of Rnmr1D
if ( is.null(samplesFile)) stop("Parser: The full path name of the samples.csv file must be specified\n")
if ( ! file.exists(samplesFile)) stop("Parser: File ", samplesFile, " does not exist\n")

# Variables
#variablesFile <- file.path(outdir,"") # output of Rnmr1D
#if ( is.null(variablesFile)) stop("Parser: The full path name of the variable file must be specified\n")
#if ( ! file.exists(variablesFile)) stop("Parser: File ", variablesFile, " does not exist\n")

# Factors
factorsFile <- file.path(outdir,"factors") # output of Rnmr1D
if ( is.null(factorsFile)) stop("Parser: The full path name of the factors file must be specified\n")
if ( ! file.exists(factorsFile)) stop("Parser: File ", factorsFile, " does not exist\n")

# ----- Initialize output file

dataMatrix.outname <- file.path(outdir,"dataMatrix.tsv")
samplesFile.outname <- file.path(outdir,"sampleMetadata.tsv")
variablesFile.outname <- file.path(outdir,"variableMetadata.tsv")

# ----- Core

factors <- read.table(file=factorsFile, header=FALSE, stringsAsFactors=FALSE, dec=".", na.strings = "NA", sep=";")
dataMatrix <- read.table(file=dataMatrixFile, header=TRUE, stringsAsFactors=FALSE, dec=".", na.strings = "NA", sep="\t")
# Select factors
samplesMetadata <- dataMatrix[1:dim(factors)[1]]
# Remove factors from matrix
Write.LOG("Parser: Transform data matrix\n", logfile=LOGFILE)
dataMatrix <- cbind( dataMatrix[1], dataMatrix[(dim(factors)[1]+1):length(dataMatrix)] )
interm <- dataMatrix[,-1]
dataMatrix <- dataMatrix[,1]
dataMatrix <- data.frame(dataMatrix, interm)
# Transposition of matrix without factors
dataMatrix <- t(dataMatrix)
# Buckets order
buckets.name <- row.names(dataMatrix)[-1]
order <- seq(1, length(buckets.name))
buckets <- data.frame(buckets.name, order)
if ( ! is.null(dataMatrix) ) {
     write.table(dataMatrix, file=dataMatrix.outname, sep="\t", dec=".", row.names=TRUE, col.names=FALSE, quote=F)
     Write.LOG("Parser: Write new dataMatrix\n", logfile=LOGFILE)
}else {
     stop("Parser: Data matrix is null. The format might be wrong. Contact administrator for help.\n")
}
# Samples x factors
if ( ! is.null(samplesMetadata) ) {
     write.table(samplesMetadata, file=samplesFile.outname, sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=F)
     Write.LOG("Parser: Write sample meta data\n", logfile=LOGFILE)
}else {
     stop("Parser: Sample meta data is null. The format might be wrong. Contact administrator for help.\n")
}
# Variables x order
if ( ! is.null(buckets) ) {
     write.table(buckets, file=variablesFile.outname, sep="\t", dec=".", col.names=TRUE, row.names=FALSE, quote=F)
     Write.LOG("Parser: Write variable meta data\n", logfile=LOGFILE)
}else {
     stop("Parser: Variable meta data is null. The format might be wrong. Contact administrator for help.\n")
}
Write.LOG("\nParser: DONE\n", logfile=LOGFILE)












