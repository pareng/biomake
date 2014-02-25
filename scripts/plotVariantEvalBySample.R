#! /usr/bin/Rscript --vanilla

## This scripts plots the results from GATK variant call evaluation with sample stratification.

## Usage: plotVariantEvalBySample.R in.grp out.pdf

## The input file should be a GATK report created by running GATK VariantEval with options:
##  -ST Filter -ST Sample --doNotUseAllStandardModules --evalModule CompOverlap \
##  --evalModule CountVariants --evalModule TiTvVariantEvaluator --dbsnp <dbSNP>

## This script requires the gsalib package
## It can bee installed by calling: install.packages("gsalib")

## Useful GATK pages for interpreting the VariantEval report:
## http://www.broadinstitute.org/gatk/guide/article?id=2361
## http://www.broadinstitute.org/gatk/guide/tagged?tag=varianteval

library(gsalib)  
library(lattice)

## Get input filename from commandline argument
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  message("Usage: plotVariantEvalBySample.R in.grp out.pdf\nSee comments in script for further usage information.\n")
  stop("incorrect number of arguments")
}
in.fn <- args[1]
out.fn <- args[2]

## Read input file
ev <- gsa.read.gatkreport(in.fn)

## Open plot file
pdf(out.fn)

## Plot results from CompOverlap

x <- ev$CompOverlap
x <- x[x$Sample != "all", ]


print(stripplot(nEvalVariants ~ Sample | Filter + Novelty, data=x,
                subset=(Filter != "raw" & Novelty != "all"),
                scales=list(x=list(draw=FALSE)),
                main="Number of variants per sample\n(variants are classified as known/novel and called/filtered)",
                ylab="Number of variants", xlab="Sample"
                ))

print(stripplot(compRate ~ Sample | Filter, data=x,
                subset=(Novelty == "all"), index.cond=list(c(3,1,2)),
                scales=list(x=list(draw=FALSE)), layout=c(3,1),
                main="Frequency of known sites",
                ylab="Known sites (%)", xlab="Sample",
                ))

print(stripplot(concordantRate ~ Sample | Filter, data=x,
                subset=(Novelty == "known"), index.cond=list(c(3,1,2)),
                scales=list(x=list(draw=FALSE)), layout=c(3,1),
                main="Concordance with dbSNP for variants at known sites",
                ylab="Concordant variants (%)", xlab="Sample"
                ))


## Plot results from CountVariants

x <- ev$CountVariants
x <- x[x$Sample != "all", ]

count.metrics <- c("Number of SNPs" = "nSNPs",
                   "Number of insertions" = "nInsertions",
                   "Number of deletions" = "nDeletions",
                   "Number of complex variants" = "nComplex",
                   "Number of other variants (MNP, symbolic and mixed)" = "I(nMNPs + nSymbolic + nMixed)",
                   "Number of heterozygous sites" = "nHets",
                   "Number of homozygous-reference sites" = "nHomRef",
                   "Number of homozygous-alternative sites" = "nHomVar",
                   "Number of \"no calls\"" = "nNoCalls",
                   "Number of private variants" = "nSingletons"
                   )

for(metric.name in names(count.metrics)) {
  print(stripplot(as.formula(paste(count.metrics[metric.name], "~ Sample | Filter + Novelty")),
                  data=x, subset=(Filter != "raw" & Novelty != "all"),
                  scales=list(x=list(draw=FALSE)),
                  main=paste(metric.name, "per sample"), ylab=metric.name, xlab="Sample"
                  ))
}


## Plot results from TiTvVariantEvaluator

x <- ev$TiTvVariantEvaluator
x <- x[x$Sample != "all", ]

print(stripplot(tiTvRatio ~ Sample | Filter + Novelty, data=x,
                subset=(Filter != "raw" & Novelty != "all"),
                scales=list(x=list(draw=FALSE)),
                main="Transition / transversion ratio for each sample",
                ylab="Transition / transversion ratio", xlab="Sample"
                ))

## Close plot file

dev.off()

## Done.

sessionInfo()

