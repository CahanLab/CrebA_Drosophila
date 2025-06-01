library(ggplot2)

raw_data = read.csv("../output/range_peak_regions/fkh_sage_intersect_raw.narrowPeak", sep = '\t', header = FALSE)
rownames(raw_data) = raw_data$V4

reformat_data = read.csv('../output/range_peak_regions/fkh_sage_intersect_untrimmed.narrowPeak', sep = '\t', header = FALSE)
rownames(reformat_data) = raw_data$V4

reformat_data = reformat_data[rownames(raw_data), ]

dist_diff = raw_data$V10 - reformat_data$V10

data=data.frame(value=dist_diff)

# basic histogram
p <- ggplot(data, aes(x=value)) + 
  geom_histogram()
