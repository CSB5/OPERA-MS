mkdir OPERA-MS-LOG/
cp --parent  $1/intermediate_files/*log $1/intermediate_files/*/*out $1/intermediate_files/*/*err $1/intermediate_files/*/*log $1/intermediate_files/strain_analysis/*/*/*log* $1/intermediate_files/opera_long_read/GAPFILLING/*log $1/intermediate_files/opera_long_read/GAPFILLING/*err $1/intermediate_files/opera_long_read/GAPFILLING/TILLING/*out $1/intermediate_files/opera_long_read/GAPFILLING/TILLING/*err $1/intermediate_files/reference_clustering/NUCMER_OUT/*LOG* $1/intermediate_files/reference_clustering/NUCMER_OUT/*log* $1/intermediate_files/reference_clustering/MASH/*.err $1/intermediate_files/reference_clustering/MASH/*.out OPERA-MS-LOG/
rm  OPERA-MS-LOG.tar.gz
tar -cvzf OPERA-MS-LOG.tar.gz OPERA-MS-LOG/
rm -r OPERA-MS-LOG/
