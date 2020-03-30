import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import sys


dist_mat = sys.argv[1]
thres = float(sys.argv[2])
outfile = sys.argv[3]

dist_dat = pd.read_csv(dist_mat, sep="\t", index_col=0, header = 0)#0.01
column_names = dist_dat.columns
dist_dat.dropna(axis=1, how='all', inplace = True)

clustering = AgglomerativeClustering(linkage = "single", n_clusters = None, compute_full_tree = True, distance_threshold = thres, affinity = "precomputed").fit(dist_dat)
with open(outfile, "w") as fp:
    fp.write("clusters\n")
    for clust in zip(column_names, clustering.labels_):
        fp.write("{}\t{}\n".format(clust[0], clust[1]))
