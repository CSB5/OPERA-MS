
install: sigma opera perlmodule mummer

clean : sigmaclean operaclean

sigma: 
	cd SIGMA&&make;

opera: 
	cd OPERA-LG&&make install;

mummer:
	cd utils&&sh install_mummer3.23.sh

perlmodule:
	cd utils&&perl install_perl_module.pl;

sigmaclean:
	cd SIGMA&&make clean;

operaclean:
	cd OPERA-LG&&make clean;
