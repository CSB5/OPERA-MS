install:
	cd src; make clean;sed -i "s/int Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS = 10000000;/int Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS = 100000;/" Configure.cpp; sed -i "s/int Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS = 1;/int Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS = 5;/" Configure.cpp;make install; mv ../bin/OPERA-LG ../bin/OPERA-LG-fast; make clean;sed -i "s/int Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS = 100000;/int Configure::MAX_NUMBER_OF_PARTIAL_SCAFFOLDS = 10000000;/" Configure.cpp; sed -i "s/int Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS = 5;/int Configure::PERCENTAGE_OF_INCREASING_THRESHOLD_GRAPHS = 1;/" Configure.cpp;make install;

clean:
	rm bin/OPERA-LG; cd src; make clean;
