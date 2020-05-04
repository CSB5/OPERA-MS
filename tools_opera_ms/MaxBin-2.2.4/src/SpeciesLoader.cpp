#include "SpeciesLoader.h"

SpeciesLoader::SpeciesLoader(char *species_list, char *dir)
{
	init();
	getSeqs(species_list, dir);
}

SpeciesLoader::~SpeciesLoader()
{
	int i;
	if (info != NULL)
	{
		for (i = 0; i < species_num; i++)
		{
			free(info[i].id);
			free(info[i].name);
			free(info[i].NCBIid);
			free(info[i].domain);
			free(info[i].phylum);
			free(info[i].classs);
			free(info[i].order);
			free(info[i].family);
			free(info[i].genus);
			free(info[i].species);
			delete(info[i].seq);
		}
		free(info);
	}
}

void SpeciesLoader::init()
{
	SIZE_INCREMENT = 128;
	species_num = 0;
	maxnum = 0;
	info = NULL;
	err = false;
}

bool SpeciesLoader::isError()
{
	return err;
}

void SpeciesLoader::getSeqs(char *list, char *dir)
{
	fstream *fs = new fstream(list);
	char str[1024], *pch, *token;
	species_info_str *tempinfo;
	int i, curr_num;

	if (!(fs->is_open()))
	{
		delete(fs);
		err = true;
		return;
	}

	while(!fs->eof())
	{
		memset(str, '\0', 1024);
		fs->getline(str, 1023);
		pch = str;
		i = 0;
		while(!(*pch == '\0' || *pch == '\r' || *pch == '\n'))
		{
			token = getToken(&pch, '\t');
			switch (i)
			{
				case 0: //taxon oid
					if (species_num >= maxnum)
					{
						maxnum = maxnum + SIZE_INCREMENT;
						tempinfo = (species_info_str*)malloc(sizeof(species_info_str) * maxnum);
						memset(tempinfo, '\0', sizeof(species_info_str) * maxnum);
						if (species_num > 0)
						{
							memcpy(tempinfo, info, sizeof(species_info_str) * species_num);
							free(info);
						}
						info = tempinfo;
					}
					curr_num = species_num;
					species_num++;
					copyString(&(info[curr_num].id), token);
				break;

				case 1: // Domain
					copyString(&(info[curr_num].domain), token);
				break;

				case 2: // Complete or draft
					if (strcmp(token, "Finished") == 0)
					{
						info[curr_num].complete = true;
					}
				break;

				case 3: // Name
					copyString(&(info[curr_num].name), token);
				break;

				case 4: // Object ID - ignore for now
				break;

				case 5: // NCBI Taxon ID
					copyString(&(info[curr_num].NCBIid), token);
				break;

				case 6: // Phylum
					copyString(&(info[curr_num].phylum), token);
				break;

				case 7: // Class
					copyString(&(info[curr_num].classs), token);
				break;

				case 8: // Order
					copyString(&(info[curr_num].order), token);
				break;

				case 9: // Family
					copyString(&(info[curr_num].family), token);
				break;

				case 10: // Genus
					copyString(&(info[curr_num].genus), token);
				break;

				case 11: // Species
					copyString(&(info[curr_num].species), token);
				break;

				default:
				break;
			}
			i++;
			if (i > 11)
			{
				break;
			}
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
	fs->close();
	delete(fs);

	for (i = 0; i < species_num; i++)
	{
#ifdef WIN32
		sprintf(str, "%s\\%s\\%s.fna", dir, info[i].id, info[i].id);
#else
		sprintf(str, "%s/%s/%s.fna", dir, info[i].id, info[i].id);
#endif //WIN32
		fs = new fstream(str);
		if (!(fs->is_open()))
		{
			delete(fs);
			err = true;
			return;
		}
		info[i].seq = new string();
		while(!(fs->eof()))
		{
			memset(str, '\0', 1024);
			fs->getline(str, 1023);
			if (str[0]!= '\0' && str[0] != '\r' && str[0] != '\n')
			{
				if (str[0] == '>')
				{
					if (info[i].length > 0)
					{
						info[i].seq->append("N");
						info[i].length = info[i].length + 1;
					}
				}
				else
				{
					info[i].seq->append(str);
					info[i].length = info[i].length + strlen(str);
				}
			}
		}
		fs->close();
		delete(fs);
	}
}

char* SpeciesLoader::getToken(char **input, char sep)
{
	char *token;
	token = *input;
	while(**input != sep && **input != '\n' && **input != '\r' && **input != '\0')
	{
		(*input)++;
	}
	**input = '\0';
	(*input)++;
	return token;
}

void SpeciesLoader::copyString(char **dest, char *source)
{
	int len = strlen(source);
	*dest = (char*)malloc(sizeof(char) * (len + 1));
	memset(*dest, '\0', sizeof(char) * (len + 1));
	memcpy(*dest, source, sizeof(char) * len);
}

int SpeciesLoader::getSpeciesNum()
{
	return species_num;
}

char* SpeciesLoader::getID(int index)
{
	return info[index].id;
}

char* SpeciesLoader::getName(int index)
{
	return info[index].name;
}

char* SpeciesLoader::getNCBIID(int index)
{
	return info[index].NCBIid;
}

char* SpeciesLoader::getDomain(int index)
{
	return info[index].domain;
}

char* SpeciesLoader::getPhylum(int index)
{
	return info[index].phylum;
}

char* SpeciesLoader::getClass(int index)
{
	return info[index].classs;
}

char* SpeciesLoader::getOrder(int index)
{
	return info[index].order;
}

char* SpeciesLoader::getFamily(int index)
{
	return info[index].family;
}

char* SpeciesLoader::getGenus(int index)
{
	return info[index].genus;
}

char* SpeciesLoader::getSpecies(int index)
{
	return info[index].species;
}

bool SpeciesLoader::getComplete(int index)
{
	return info[index].complete;
}


string* SpeciesLoader::getSeq(int index)
{
	return info[index].seq;
}
