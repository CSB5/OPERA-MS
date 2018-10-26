install:
	cd src; sed -i "s/Configure::ONLY_BUNDLE = false;/Configure::ONLY_BUNDLE = true;/" Configure.cpp; make install; mv ../bin/OPERA-LG ../../bin/bundler; make clean;sed -i "s/Configure::ONLY_BUNDLE = true;/Configure::ONLY_BUNDLE = false;/" Configure.cpp; make install

clean:
	rm bin/OPERA-LG; cd src; make clean;
