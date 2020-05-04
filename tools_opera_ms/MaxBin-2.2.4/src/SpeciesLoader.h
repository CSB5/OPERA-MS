#ifndef __SPECIESLOADER_H__
#define __SPECIESLOADER_H__

#include "global_inc.h"
#include <fstream>

typedef struct species_info_str
{
	char *id;
	char *name;
	char *NCBIid;
	char *domain;
	char *phylum;
	char *classs;
	char *order;
	char *family;
	char *genus;
	char *species;
	int length;
	string *seq;
	bool complete;
} species_info_str;

class SpeciesLoader
{
	public:
		SpeciesLoader(char *species_list, char *dir);
		~SpeciesLoader();
		int getSpeciesNum();
		char* getID(int index);
		char* getName(int index);
		char* getNCBIID(int index);
		char* getDomain(int index);
		char* getPhylum(int index);
		char* getClass(int index);
		char* getOrder(int index);
		char* getFamily(int index);
		char* getGenus(int index);
		char* getSpecies(int index);
		string* getSeq(int index);
		bool getComplete(int index);
		bool isError();

	private:
		// Variables
		int SIZE_INCREMENT;
		int species_num;
		int maxnum;
		species_info_str *info;
		bool err;

		// Functions
		void init();
		void getSeqs(char *list, char *dir);
		char* getToken(char **input, char sep);
		void copyString(char **dest, char *source);
};




#endif