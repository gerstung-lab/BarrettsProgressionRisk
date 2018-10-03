# Run this script from within the directory containing bam files.

if (length(list.files('.', 'bam')) <= 0)
  stop("No bam files found in current directory. Exiting.")


suppressPackageStartupMessages( library(Biobase) )
suppressPackageStartupMessages( library(QDNAseq) )
suppressPackageStartupMessages( library(tidyverse) )

# get/read bin annotations
bins <- getBinAnnotations(binSize = 15)

saveRDS(bins, "15kbp.rds")
#bins <- readRDS("15kbp.rds")

# write bin annotations
pData(bins) %>%
  mutate_at(vars(start, end), funs(as.integer)) %>%
  write_tsv("15kbp.txt")

# process BAM files obtaining read counts within bins
readCounts <- binReadCounts(bins, path = ".")

pData(readCounts) %>%
  rownames_to_column(var = "sample") %>%
  select(id = name, everything()) %>%
  write_tsv("readCountSummary.txt")

exportBins(readCounts, "readCounts.txt", logTransform = FALSE)

saveRDS(readCounts, "readCounts.rds")
#readCounts <- readRDS("readCounts.rds")

# apply filters for which bins are used, including loess residuals of calibration set
# and blacklisted regions both taken from bin annotations
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)

# fix for issue with copy number segmentation arising from zero count bins
#fData(readCountsFiltered)[which(assayData(readCountsFiltered)$counts == 0),]$use <- FALSE

# plot median read counts as a function of GC content and mappability as an isobar plot
pdf("isobar.pdf")
isobarPlot(readCountsFiltered)
dev.off()

# estimate the correction for GC content and mappability
readCountsCorrected <- estimateCorrection(readCountsFiltered)

# output raw and fitted read counts
features <- fData(readCountsCorrected) %>%
  as.data.frame %>%
  rownames_to_column(var = "location") %>%
  transmute(location, chrom = chromosome, start = as.integer(start), end = as.integer(end))

rawReadCounts <- assayData(readCountsCorrected)$counts %>%
  as.data.frame %>%
  rownames_to_column(var = "location")
features %>%
  left_join(rawReadCounts, by = "location") %>%
  write_tsv("rawReadCounts.txt")

fittedReadCounts <- assayData(readCountsCorrected)$fit %>%
  as.data.frame %>%
  rownames_to_column(var = "location") %>%
  mutate_if(is.numeric, funs(round(., digits = 3)))
features %>%
  left_join(fittedReadCounts, by = "location") %>%
  write_tsv("fittedReadCounts.txt")

# noise plot showing relationship between the observed standard deviation in the data
# and its read depth
pdf("noise.pdf")
noisePlot(readCountsCorrected)
dev.off()

save.image("qdnaseq.RData")
