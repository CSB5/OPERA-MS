
args = commandArgs(trailingOnly=TRUE)

dist.mat <- args[1]
thres <- args[2]
outfile <- args[3]
dist.dat <- read.table(dist.mat)#0.01

## clustering
cluster.full <- hclust(as.dist(dist.dat), method = 'single' )

clusters <- cutree(cluster.full, h=thres)
write.table(data.frame(clusters), outfile, sep='\t')
