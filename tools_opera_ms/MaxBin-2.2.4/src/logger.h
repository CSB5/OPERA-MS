#ifndef __LOGGER_H__
#define __LOGGER_H__

#include "global_inc.h"
#include <fstream>
#include <time.h>

class logger
{
	public:
		logger(char *input_filename);
		~logger();
		void writelog(char *c, bool is_verbose);
		void writelog(char const *c, bool is_verbose);
		void writetime(time_t t);
		void setVerbose(bool is_verbose);
	private:
		// Variables
		char *file;
		fstream *log;
		bool verbosetype;
		char str[1024];
};


#endif
