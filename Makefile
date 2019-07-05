
install: sigma opera short_read_analysis mummer perlmodule 

clean : sigmaclean operaclean

sigma: 
	cd SIGMA&&make;

opera: 
	cd OPERA-LG&&make install;

short_read_analysis:
	cd SHORT_READ_ANALYSIS/src/&&make;

mummer:
	cd utils&&sh install_mummer3.23.sh

perlmodule:
	cd utils&&perl install_perl_module.pl;

sigmaclean:
	cd SIGMA&&make clean;

operaclean:
	cd OPERA-LG&&make clean;
