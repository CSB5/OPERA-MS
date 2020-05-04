#include "logger.h"

logger::logger(char *input_filename)
{
	file = input_filename;
	log = new fstream(file, ios::out);
	verbosetype = false;
}

logger::~logger()
{
	log->close();
	delete(log);
}

void logger::setVerbose(bool is_verbose)
{
	verbosetype = is_verbose;
}

void logger::writelog(char *c, bool is_verbose)
{
	if (verbosetype == true || is_verbose == true)
	{
		printf("%s", c);
		log->write(c, strlen(c));
	}
}

void logger::writelog(char const *c, bool is_verbose)
{
	if (verbosetype == true || is_verbose == true)
	{
		printf("%s", c);
		log->write(c, strlen(c));
	}
}

void logger::writetime(time_t t)
{
	int day, hour, min, sec;
	min = (int)t / 60;
	sec = (int)t % 60;
	hour = min / 60;
	min = min % 60;
	day = hour / 24;
	hour = hour % 24;
	sprintf(str, "Elapsed time: %2d days %02d:%02d:%02d\n", day, hour, min, sec);
	printf("%s", str);
	log->write(str, strlen(str));
}

