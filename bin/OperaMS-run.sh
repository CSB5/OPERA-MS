#!/bin/sh

perl bin/runOperaMS.pl $1 >> $2/run_log.txt 2>&1
#perl /mnt/projects/bertrandd/sigma/OperaMS2.0/runOperaMS.pl $OPERAMSconfig >> $OPERAMSlogdir/run_log.txt 2>&1
