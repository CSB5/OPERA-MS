#include "AbundanceLoader.h"

#define VERY_SMALL_NUM 0.0001

AbundanceLoader::AbundanceLoader(char *file)
{
	num = 0;
	max = 0;
	open_flag = false;
	parse_flag = false;
	ADD_SIZE = 1024;
	parse(file);
}

AbundanceLoader::~AbundanceLoader()
{
	int i;
	for (i = 0; i < num; i++)
	{
		free(header[i]);
	}
	free(header);
	free(abund);
}

int AbundanceLoader::getNum()
{
	return num;
}

bool AbundanceLoader::is_open()
{
	return open_flag;
}

bool AbundanceLoader::is_parse_success()
{
	return parse_flag;
}

char* AbundanceLoader::get_parse_error()
{
	return str_err;
}

double AbundanceLoader::getAbundance(int index)
{
	if (index < num)
	{
		return abund[index];
	}
	else
	{
		return (double)-1;
	}
}

double AbundanceLoader::getAbundance(char *header)
{
	abund_it = abund_hash.find(header);
	if (abund_it != abund_hash.end())
	{
		return(abund_hash[header]);
	}

	return (double)-1;
}

void AbundanceLoader::parse(char *file)
{
	char str[1024];
	char *c;
	fstream *fs = new fstream(file, ios::in);

	if (fs->is_open() == false)
	{
		return;
	}
	else
	{
		open_flag = true;
	}

	while (!fs->eof())
	{
		memset(str, '\0', 1024);
		memset(str_err, '\0', 1024);
		fs->getline(str, 1023);
		memcpy(str_err, str, 1023);
		while (str[strlen(str) - 1] == '\n' || str[strlen(str) - 1] == '\r' || str[strlen(str) - 1] == ' ' || str[strlen(str) - 1] == '\t')
		{
			str[strlen(str) - 1] = '\0';
		}
		if (str[0] == '\0')
		{
			continue;
		}

		c = str;
		while (*c != '\t' && *c != ' ' && *c != ',' && *c != ';')
		{
			c++;
			if (*c == '\0')
			{
				return;
			}
		}
		*c = '\0';
		c++;
		while (*c == '\t' && *c == ' ' && *c == ',' && *c == ';')
		{
			c++;
		}
		if (*c == '\0')
		{
			return;
		}
		else
		{
			_insertElement(str, atof(c));
		}
		if (fs->eof())
		{
			break;
		}
		else if (fs->fail())
		{
			fs->clear();
		}
	}
	parse_flag = true;
	
	fs->close();
	delete(fs);
}

void AbundanceLoader::_insertElement(char *h, double ab)
{
	char **tempheader;
	double *tempabund;
	if (num >= max)
	{
		// Increase array size
		max = max + ADD_SIZE;
		tempheader = (char**)malloc(sizeof(char*) * max);
		memset(tempheader, '\0', sizeof(char*) * max);
		tempabund = (double*)malloc(sizeof(double) * max);
		memset(tempabund, '\0', sizeof(double) * max);
		if (num > 0)
		{
			memcpy(tempheader, header, sizeof(char*) * num);
			memcpy(tempabund, abund, sizeof(double) * num);
			free(header);
			free(abund);
		}
		header = tempheader;
		abund = tempabund;
	}
	// Check trailing spaces
	while (h[strlen(h) - 1] == ' ' || h[strlen(h) - 1] == '\t')
	{
		h[strlen(h) - 1] = '\0';
	}

	header[num] = (char*)malloc(sizeof(char) * (strlen(h) + 1));
	memset(header[num], '\0', sizeof(char) * (strlen(h) + 1));
	memcpy(header[num], h, sizeof(char) * strlen(h));
	if (ab < VERY_SMALL_NUM)
	{
		ab = VERY_SMALL_NUM;
	}
	abund[num] = ab;

	abund_hash[header[num]] = ab;
	num++;
}

