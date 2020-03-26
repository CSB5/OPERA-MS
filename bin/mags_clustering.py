import sys
import pandas as pd
#from scipy.cluster.hierarchy import dendrogram, linkage
from scipy import cluster


dist_mat = sys.argv[1]
thres = sys.argv[2]
outfile = sys.argv[3]


#with open(dist_mat)
dist_dat = pd.read_csv(dist_mat, sep="\t")#0.01

cluster_full = cluster.hierarchy.single(dist_dat)
clusters = cluster.hierarchy.cut_tree(cluster_full, height = thres)

print(clusters)

#linked = linkage(dist_dat, 'single')

#dendrogram(linked,
#            orientation='top',
#            labels=labelList,
#            distance_sort='descending',
#            show_leaf_counts=True)


## clustering
#cluster.full <- hclust(as.dist(dist.dat), method = 'single' )

#clusters <- cutree(cluster.full, h=thres)
#write.table(data.frame(clusters), outfile, sep='\t')
