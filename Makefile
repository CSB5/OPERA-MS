
install: sigma opera

clean : sigmaclean operaclean

sigma: 
	cd sigma_testing_jim&&make;

opera: 
	cd OPERA-LG_v2.1.0&&make install;

sigmaclean:
	cd sigma_testing_jim&&make clean;

operaclean:
	cd OPERA-LG_v2.1.0&&make clean;
