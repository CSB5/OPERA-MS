#include "kmerMap.h"

kmerMap::kmerMap(int inputlen, bool is_symmetric)
{
	len_mer = inputlen;
	entry_num = (int)pow((float)4, (float)len_mer);
	symmetric = is_symmetric;
	buildKmerMappingTable();
}

kmerMap::~kmerMap()
{
	if (kmer_small != NULL)
	{
		free(kmer_small);
	}
	if (kmer_large != NULL)
	{
		free(kmer_large);
	}
	free(kmer_map);
}

int kmerMap::getMapping(char *inputKmer)
{
	unsigned int i;
	i = calNum(inputKmer);
	if (err != false)
	{
		return(kmer_map[i]);
	}
	else
	{
		return -1;
	}
}

int kmerMap::getReverseMapping(char *inputKmer)
{
	unsigned int i;
	char *tmprev;
	tmprev = (char*)malloc(len_mer + 1);
	memset(tmprev, '\0', len_mer + 1);
	revSeq(inputKmer, tmprev);
	i = calNum(tmprev);
	free(tmprev);
	if (err != false)
	{
		return(kmer_map[i]);
	}
	else
	{
		return -1;
	}
}

int kmerMap::getReverseMapping(int index)
{
	unsigned int i;
	char *tmpchr, *tmprev;
	tmpchr = (char*)malloc(len_mer + 1);
	memset(tmpchr, '\0', len_mer + 1);
	tmprev = (char*)malloc(len_mer + 1);
	memset(tmprev, '\0', len_mer + 1);
	numToString(index, tmpchr);
	revSeq(tmpchr, tmprev);
	i = calNum(tmprev);
	free(tmpchr);
	free(tmprev);
	if (err != false)
	{
		return(kmer_map[i]);
	}
	else
	{
		return -1;
	}
}

int kmerMap::getEntryNum()
{
	if (symmetric == true)
	{
		if (len_mer % 2 == 0)
		{
			return(((int)pow((float)4, (float)len_mer) + (int)pow((float)4, (float)(len_mer / 2)))/ 2);
		}
		else
		{
			return((int)pow((float)4, (float)len_mer) / 2);
		}
	}
	else
	{
		return entry_num;
	}
}

void kmerMap::buildKmerMappingTable()
{
	int i, j, k, p, match;
	char *tmpchr, *tmprev;
	if (symmetric == true)
	{
		if (len_mer % 2 == 0)
		{
			i = ((int)pow((float)4, (float)len_mer) + (int)pow((float)4, (float)(len_mer / 2)))/ 2;
		}
		else
		{
			i = (int)pow((float)4, (float)len_mer) / 2;
		}
		kmer_small = (int*)malloc(i * sizeof(int));
		kmer_large = (int*)malloc(i * sizeof(int));
		memset(kmer_small, '\0', i * sizeof(int));
		memset(kmer_large, '\0', i * sizeof(int));
		tmpchr = (char*)malloc(len_mer + 1);
		memset(tmpchr, '\0', len_mer + 1);
		tmprev = (char*)malloc(len_mer + 1);
		memset(tmprev, '\0', len_mer + 1);
		kmer_num = 0;
		for (i = 0; i < entry_num; i++)
		{
			numToString(i, tmpchr);
			revSeq(tmpchr, tmprev);
			j = calNum(tmprev);
			if (i <= j)
			{
				k = i;
				p = j;
			}
			else
			{
				p = i;
				k = j;
			}
			match = 0;
			for (j = 0; j < kmer_num; j++)
			{
				if (kmer_small[j] == k)
				{
					match = 1;
				}
			}
			if (match == 0)
			{
				kmer_small[kmer_num] = k;
				kmer_large[kmer_num] = p;
				kmer_num++;
			}
		}
		// Build mapping table
		kmer_map = (int*)malloc(sizeof(int) * entry_num);
		memset(kmer_map, '\0', sizeof(int) * entry_num);
		for (i = 0; i < kmer_num; i++)
		{
			kmer_map[kmer_small[i]] = i;
			kmer_map[kmer_large[i]] = i;
		}
		free(tmpchr);
		free(tmprev);
	}
	else
	{
		kmer_map = (int*)malloc(sizeof(int) * entry_num);
		memset(kmer_map, '\0', sizeof(int) * entry_num);
		for (i = 0; i < entry_num; i++)
		{
			kmer_map[i] = i;
		}
		kmer_small = NULL;
		kmer_large = NULL;
	}
}

void kmerMap::revSeq(char *seq, char *result)
{
	int i, len;
	char ch;
	len = strlen(seq);
	for (i = 0; i < len; i++)
	{
		switch (*(seq + i))
		{
			case 'A':
			case 'a':
				ch = 'T';
			break;

			case 'T':
			case 't':
				ch = 'A';
			break;

			case 'C':
			case 'c':
				ch = 'G';
			break;

			case 'G':
			case 'g':
				ch = 'C';
			break;
		}
		*(result + len - i - 1) = ch;
	}
}

unsigned int kmerMap::calNum(char *str)
{
	int i;
	unsigned int ret;
	unsigned _int64 num;
	num = 0;
	err = true;
	for (i = 0; i < len_mer; i++)
	{
		num = num * 4;
		switch (str[i])
		{
			case 'A':
			case 'a':
				num = num + 0;
			break;

			case 'T':
			case 't':
				num = num + 1;
			break;

			case 'C':
			case 'c':
				num = num + 2;
			break;

			case 'G':
			case 'g':
				num = num + 3;
			break;

			default:
				num = 0;
				err = false;
			break;
		}
		if (err == false)
		{
			break;
		}
	}
	ret = (unsigned int)num;
	return ret;
}

void kmerMap::numToString(int num, char *str)
{
	int i, j, k;
	i = num;
	for (j = len_mer - 1; j >= 0; j--)
	{
		k = i % 4;
		i = i / 4;
		switch (k)
		{
			case 0:
				str[j] = 'A';
			break;

			case 1:
				str[j] = 'T';
			break;

			case 2:
				str[j] = 'C';
			break;

			case 3:
				str[j] = 'G';
			break;

			default:
			break;
		}
	}
}
