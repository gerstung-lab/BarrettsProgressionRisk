
library(BarrettsProgressionRisk)


pr = predictRisk(path='example/')

# Plot raw data
plotSegmentData(pr)

# output the risk table per sample
predictions(pr)

# Get recommendations per sample pair, assuming these were at different timepoint
rx(pr, 'example/demo_file.txt')

