
tools_dir = tools_opera_ms

ln -s $(shell pwd)/src_utils/OPERA-MS-UTILS.py $(shell pwd)/

install: sigma opera short_read_analysis mummer perlmodule 

clean : sigmaclean operaclean

sigma: 
	cd SIGMA&&make;

opera: 
	cd OPERA-LG&&make install;

short_read_analysis:
	cd SHORT_READ_ANALYSIS/src/&&make;

mummer:
	cd $(tools_dir)&&sh install_mummer3.23.sh

perlmodule:
	cd $(tools_dir)&&perl install_perl_module.pl;

sigmaclean:
	cd SIGMA&&make clean;

operaclean:
	cd OPERA-LG&&make clean;
