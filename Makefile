
install: sigma opera

clean : sigmaclean operaclean

sigma: 
	cd SIGMA&&make;

opera: 
	cd OPERA-LG&&make install;

sigmaclean:
	cd SIGMA&&make clean;

operaclean:
	cd OPERA-LG&&make clean;
