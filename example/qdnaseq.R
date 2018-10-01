library(Biobase)
library(QDNAseq)
library(tidyverse)

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

# apply correction for GC content and mappability
copyNumber <- correctBins(readCountsCorrected)

# normalize bins based on median corrected read counts
copyNumberNormalized <- normalizeBins(copyNumber)

# smooth outlier bins
copyNumberSmoothed <- smoothOutlierBins(copyNumberNormalized)

# export corrected copy number profiles as tab-delimited file and IGV track
#exportBins(copyNumberSmoothed, "copyNumberSmoothed.txt", logTransform = TRUE)
exportBins(copyNumberSmoothed, "copyNumberSmoothed.igv", format = "igv", logTransform = TRUE)

# apply circular binary segmentation (CBS) algorithm from DNAcopy package
# (by default this uses a log2 transformation)
copyNumberSegmented <- segmentBins(copyNumberSmoothed, transformFun = "sqrt")

# normalize segmented data using recursive search procedure
copyNumberSegmented <- normalizeSegmentedBins(copyNumberSegmented)

# export segmented copy number profile
exportBins(copyNumberSegmented, "copyNumber.txt", logTransform = TRUE)
exportBins(copyNumberSegmented, "copyNumberSegmented.txt", type = "segments", logTransform = TRUE)

pdf("copyNumberSegmented.pdf", width = 8, height = 6)
plot(copyNumberSegmented)
dev.off()

save.image("qdnaseq.RData")

samples <- read_tsv("samples.txt")

ids <- samples$ID
samples <- samples$Sample
names(samples) <- ids

for (id in sampleNames(copyNumberSegmented))
{
  sample <- samples[id]
  title <- id
  if (!is.na(sample) && id != sample) title <- paste(id, sample, sep = " ")
  png(paste(id, "copyNumberSegmented.png", sep = "."), width = 1200, height = 800)
  plot(copyNumberSegmented[,id], main = title)
  dev.off()
}

sessionInfo()

