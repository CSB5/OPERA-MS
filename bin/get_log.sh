mkdir OPERA-MS-LOG/

argpath=`readlink -f $1`

cp --parent  $argpath/intermediate_files/*log $argpath/intermediate_files/*/*out $argpath/intermediate_files/*/*err $argpath/intermediate_files/*/*log $argpath/intermediate_files/strain_analysis/*/*/*log* $argpath/intermediate_files/opera_long_read/GAPFILLING/*log $argpath/intermediate_files/opera_long_read/GAPFILLING/*err $argpath/intermediate_files/opera_long_read/GAPFILLING/TILLING/*out $argpath/intermediate_files/opera_long_read/GAPFILLING/TILLING/*err $argpath/intermediate_files/reference_clustering/NUCMER_OUT/*LOG* $argpath/intermediate_files/reference_clustering/NUCMER_OUT/*log* $argpath/intermediate_files/reference_clustering/MASH/*.err $argpath/intermediate_files/reference_clustering/MASH/*.out OPERA-MS-LOG/
rm  OPERA-MS-LOG.tar.gz
tar -cvzf OPERA-MS-LOG.tar.gz OPERA-MS-LOG/
rm -r OPERA-MS-LOG/
