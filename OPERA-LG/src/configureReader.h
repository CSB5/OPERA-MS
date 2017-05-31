#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include "Configure.h"
#include "LibInfo.h"
#include "CommonFunction.h"

using namespace std;

class configureReader
{
public:
	configureReader(void);
	~configureReader(void);

 private:
	bool firstLib;

	// Methods
public:
	// Read the configuration file
	int ReadConfigFile( string fileName );

private:
	// anylize the parameters
	int AnalyzeParameters( string line );
};
