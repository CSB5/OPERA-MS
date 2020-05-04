#include "fastaReader.h"

const int fastaReader::ADD_SIZE = 1024;

fastaReader::fastaReader(char *inputFile)
{
	fs = new fstream(inputFile, ios::in);
	if (fs->is_open() == false)
	{
		is_fasta = false;
		return;
	}
	is_fasta = true;
	init();
	parse();
}

fastaReader::fastaReader()
{
	is_fasta = true;
	fs = NULL;
	init();
}

fastaReader::~fastaReader()
{
	unsigned int i;
	if (header != NULL)
	{
		for (i = 0; i < seqNum; i++)
		{
			free(header[i]);
		}
		free(header);
	}
	if (seq != NULL)
	{
		for (i = 0; i < seqNum; i++)
		{
			free(seq[i]);
		}
		free(seq);
	}
	if (seqLen != NULL)
	{
		free(seqLen);
	}
	if (coverage != NULL)
	{
		free(coverage);
	}
	if (mark != NULL)
	{
		free(mark);
	}
	if (offset != NULL)
	{
		free(offset);
	}
	if (fs != NULL)
	{
		delete(fs);
	}
}

bool fastaReader::isFasta()
{
	return is_fasta;
}

unsigned int fastaReader::getNum()
{
	return seqNum;
}

char* fastaReader::getNext()
{
	if (currentPos >= seqNum)
	{
		return NULL;
	}
	else
	{
		currentPos++;
		return seq[currentPos - 1] + offset[currentPos - 1];
	}
}

char* fastaReader::getNextMarked()
{
	while (mark[currentPos] == false && currentPos < seqNum)
	{
		currentPos++;
	}
	if (currentPos >= seqNum)
	{
		return NULL;
	}
	currentPos++;
	return seq[currentPos - 1] + offset[currentPos - 1];
}

char* fastaReader::getSeqByNum(unsigned int num)
{
	if (num >= seqNum)
	{
		return NULL;
	}
	else
	{
		currentPos = num;
		return seq[num] + offset[num];
	}
}

char* fastaReader::getFullSeqByNum(unsigned int num)
{
	if (num >= seqNum)
	{
		return NULL;
	}
	else
	{
		currentPos = num;
		return seq[num];
	}
}

char* fastaReader::getHeaderByNum(unsigned int num)
{
	if (num >= seqNum)
	{
		return NULL;
	}
	else
	{
		return header[num];
	}
}

float fastaReader::getCoverageByNum(unsigned int num)
{
	if (num >= seqNum)
	{
		return (float)0;
	}
	else
	{
		return coverage[num];
	}
}

unsigned int fastaReader::getSeqLenByNum(unsigned int num)
{
	if (num >= seqNum)
	{
		return 0;
	}
	else
	{
		return seqLen[num] - offset[num];
	}
}

unsigned int fastaReader::getFullSeqLenByNum(unsigned int num)
{
	if (num >= seqNum)
	{
		return 0;
	}
	else
	{
		return seqLen[num];
	}
}

unsigned int fastaReader::getMaxLen()
{
	unsigned int maxLen, i;

	maxLen = 0;
	for (i = 0; i < seqNum; i++)
	{
		if (maxLen < (seqLen[i] - offset[i]))
		{
			maxLen = seqLen[i] - offset[i];
		}
	}
	return maxLen;
}

char* fastaReader::getCurrentHeader()
{
	return header[currentPos - 1];
}

void fastaReader::resetCurrent()
{
	currentPos = 0;
}

void fastaReader::setMark(unsigned int num)
{
	if (num >= seqNum)
	{
		return;
	}
	mark[num] = true;
}

void fastaReader::setMark(char *id)
{
	int i;
	i = getHeaderIndex(id);
	if (i != -1)
	{
		mark[i] = true;
	}
}

void fastaReader::resetMark()
{
	memset(mark, '\0', sizeof(bool) * maxNum);
}

bool fastaReader::isMark(int i)
{
	return(mark[i]);
}

void fastaReader::setOffset(unsigned int num, unsigned int input_offset)
{
	offset[num] = input_offset;
}

char* fastaReader::getPartialSeq(unsigned int num, unsigned int pos_start, unsigned int pos_end)
{
	char *retch;
	if (pos_end > seqLen[num])
	{
		pos_end = seqLen[num];
	}
	pos_start--;
	pos_end--;
	retch = (char*)malloc(pos_end - pos_start + 2);
	memset(retch, '\0', pos_end - pos_start + 2);
	memcpy(retch, seq[num] + pos_start, pos_end - pos_start + 1);
	return retch;
}

void fastaReader::separateSeq(unsigned int num, unsigned int cut1, unsigned int cut2)
{
	char *ch, *ch2, *tmpheader2;
	//char *tmpheader1;

	// sequences
	ch = (char*)malloc(cut1);
	memset(ch, '\0', cut1);
	memcpy(ch, seq[num], cut1 - 1);
	ch2 = (char*)malloc(seqLen[num] - cut2 + 1);
	memset(ch2, '\0', seqLen[num] - cut2 + 1);
	memcpy(ch2, seq[num] + cut2, seqLen[num] - cut2);

	// headers
	//tmpheader1 = (char*)malloc(strlen(header[num]) + 4);
	//sprintf(tmpheader1, "%s_c1", header[num]);
	tmpheader2 = (char*)malloc(strlen(header[num]) + 4);
	sprintf(tmpheader2, "%s_c2", header[num]);

	// Modify original sequence
	free(seq[num]);
	//free(header[num]);
	seq[num] = ch;
	//header[num] = tmpheader1;
	charSum = charSum - seqLen[num];
	seqLen[num] = cut1 - 1;
	charSum = charSum + seqLen[num];
	// Add a new sequence
	addSeq(ch2, tmpheader2);
	free(ch2);
	free(tmpheader2);
}

void fastaReader::addSeq(const char *input_seq, const char *input_header)
{
	char **tmpseq, **tmpheader;
	bool *tmpmark;
	float *tmpcov;
	unsigned int *tmpseqlen, *tmpoffset;

	if (seqNum > 0 && seqNum % 100000 == 0)
	{
		printf("Loaded %d sequences\n", seqNum);
	}

	if (seqNum >= maxNum)
	{
		// Append arrays
		maxNum = maxNum + ADD_SIZE;
		tmpseq = (char**)malloc(sizeof(char*) * maxNum);
		tmpheader = (char**)malloc(sizeof(char*) * maxNum);
		tmpseqlen = (unsigned int*)malloc(sizeof(unsigned int) * maxNum);
		tmpmark = (bool*)malloc(sizeof(bool) * maxNum);
		tmpoffset = (unsigned int*)malloc(sizeof(unsigned int) * maxNum);
		tmpcov = (float*)malloc(sizeof(float) * maxNum);
		memset(tmpseq, '\0', sizeof(char*) * maxNum);
		memset(tmpheader, '\0', sizeof(char*) * maxNum);
		memset(tmpmark, '\0', sizeof(bool) * maxNum);
		memset(tmpseqlen, '\0', sizeof(unsigned int) * maxNum);
		memset(tmpoffset, '\0', sizeof(unsigned int) * maxNum);
		memset(tmpcov, '\0', sizeof(float) * maxNum);
		if (seq != NULL)
		{
			memcpy(tmpseq, seq, sizeof(char*) * seqNum);
			memcpy(tmpheader, header, sizeof(char*) * seqNum);
			memcpy(tmpmark, mark, sizeof(bool) * seqNum);
			memcpy(tmpseqlen, seqLen, sizeof(unsigned int) * seqNum);
			memcpy(tmpoffset, offset, sizeof(unsigned int) * seqNum);
			memcpy(tmpcov, coverage, sizeof(float) * seqNum);
			free(seq);
			free(header);
			free(mark);
			free(seqLen);
			free(offset);
			free(coverage);
		}
		seq = tmpseq;
		header = tmpheader;
		mark = tmpmark;
		seqLen = tmpseqlen;
		offset = tmpoffset;
		coverage = tmpcov;
	}
	seqLen[seqNum] = strlen(input_seq);
	charSum = charSum + seqLen[seqNum];
	seq[seqNum] = (char*)malloc(seqLen[seqNum] + 1);
	header[seqNum] = (char*)malloc(strlen(input_header) + 1);
	memset(seq[seqNum], '\0', seqLen[seqNum] + 1);
	memset(header[seqNum], '\0', strlen(input_header) + 1);
	memcpy(seq[seqNum], input_seq, seqLen[seqNum]);
	memcpy(header[seqNum], input_header, strlen(input_header));
	convertToUpper(seq[seqNum]);
	//coverage[seqNum] = getCoverage(header[seqNum]);
	//trimHeader(header[seqNum], ' ');
	headerHash[header[seqNum]] = seqNum;
	seqNum++;
}

void fastaReader::removeSeq(unsigned int num)
{
	if (num >= seqNum)
	{
		return;
	}
	if (num < seqNum - 1)
	{
		memcpy(seq + num, seq + num + 1, sizeof(char*) * (seqNum - num - 1));
		memcpy(header + num, header + num + 1, sizeof(char*) * (seqNum - num - 1));
		memcpy(mark + num, mark + num + 1, sizeof(char*) * (seqNum - num - 1));
		charSum = charSum - seqLen[num];
		memcpy(seqLen + num, seqLen + num + 1, sizeof(char*) * (seqNum - num - 1));
		memcpy(offset + num, offset + num + 1, sizeof(char*) * (seqNum - num - 1));
		memcpy(coverage + num, coverage + num + 1, sizeof(char*) * (seqNum - num - 1));
	}
	seqNum--;
}

void fastaReader::splitMarkedSeq(fastaReader **marked, fastaReader **unmarked)
{
	unsigned int i;
	fastaReader *frmark, *frunmark;
	frmark = new fastaReader();
	frunmark = new fastaReader();

	for (i = 0; i < seqNum; i++)
	{
		if (mark[i] == true)
		{
			frmark->addSeq(seq[i], header[i]);
		}
		else
		{
			frunmark->addSeq(seq[i], header[i]);
		}
	}

	*marked = frmark;
	*unmarked = frunmark;
}

fastaReader* fastaReader::getMarkedSeq()
{
	unsigned int i;
	fastaReader *frmark;
	frmark = new fastaReader();

	for (i = 0; i < seqNum; i++)
	{
		if (mark[i] == true)
		{
			frmark->addSeq(seq[i], header[i]);
		}
	}

	return frmark;
}

fastaReader* fastaReader::getUnmarkedSeq()
{
	unsigned int i;
	fastaReader *frunmark;
	frunmark = new fastaReader();

	for (i = 0; i < seqNum; i++)
	{
		if (mark[i] == false)
		{
			frunmark->addSeq(seq[i], header[i]);
		}
	}

	return frunmark;
}

void fastaReader::init()
{
	header = NULL;
	seq = NULL;
	seqLen = NULL;
	mark = NULL;
	coverage = NULL;
	offset = NULL;
	currentPos = 0;
	seqNum = 0;
	maxNum = 0;
	charSum = 0;
}

void fastaReader::reverseSeq(unsigned int num)
{
	int i, len;
	char *temp;
	if (num >= seqNum)
	{
		return;
	}
	len = strlen(seq[num]);
	temp = (char*)malloc((sizeof(char) * len) + 1);
	memset(temp, '\0', (sizeof(char) * len) + 1);
	for (i = len - 1; i >= 0; i--)
	{
		temp[len - i - 1] = revChar(seq[num][i]);
	}
	free(seq[num]);
	seq[num] = temp;
}

char fastaReader::revChar(char ch)
{
	switch (ch)
	{
		case 'A':
		case 'a':
			return 'T';
		break;

		case 'T':
		case 't':
			return 'A';
		break;

		case 'C':
		case 'c':
			return 'G';
		break;

		case 'G':
		case 'g':
			return 'C';
		break;

		default:
			return 'N';
		break;
	}
}

void fastaReader::convertToUpper(char *seq)
{
	unsigned int len, i;
	len = strlen(seq);
	for (i = 0; i < len; i++)
	{
		if (seq[i] >= 97 && seq[i] <= 122)
		{
			seq[i] = seq[i] - 32;
		}
	}
}

void fastaReader::trimHeader(char *input_header, char div)
{
	char *c;
	c = input_header;
	while (*c != div && *c != '\0')
	{
		c++;
	}
	if (*c == div)
	{
		*c = '\0';
	}
}

void fastaReader::parse()
{
	string tempstr;
	unsigned int i;
	char *curr_header;
	bool start_fasta;
	start_fasta = false;
	curr_header = NULL;
	while (!fs->eof())
	{
		memset(str, '\0', 1024);
		fs->getline(str, 1023);
		while (str[strlen(str) - 1] == '\n' || str[strlen(str) - 1] == '\r')
		{
			str[strlen(str) - 1] = '\0';
		}
		if (str[0] == '>')
		{
			if (start_fasta == false)
			{
				start_fasta = true;
			}
			else
			{
				addSeq(tempstr.data(), curr_header);
			}
			if (curr_header != NULL)
			{
				free(curr_header);
			}
			curr_header = (char*)malloc(strlen(str));
			memset(curr_header, '\0', strlen(str));
			memcpy(curr_header, str + 1, strlen(str) - 1);
			for (i = 0; i < strlen(curr_header); i++)
			{
				if (curr_header[i] == ' ' || curr_header[i] == '\t')
				{
					curr_header[i] = '\0';
					break;
				}
			}
			tempstr.clear();
		}
		else if (start_fasta == false)
		{
			continue;
		}
		else if (str[0] != '\0')
		{
			if (tempstr.empty() == true)
			{
				tempstr.assign(str);
			}
			else
			{
				tempstr.append(str);
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
	addSeq(tempstr.data(), curr_header);
	if (curr_header != NULL)
	{
		free(curr_header);
	}
	fs->close();
}

int fastaReader::getHeaderIndex(char *header)
{
	map<char*, int, str_cmp>::iterator iter, iter2, iter3;
	iter = headerHash.begin();
	iter2 = headerHash.find(header);
	iter3 = headerHash.end();
	if (headerHash.find(header) != headerHash.end())
	{
		return(headerHash[header]);
	}
	return -1;
}

void fastaReader::duplicateSeqByNum(unsigned int num)
{
	char *c;
	c = (char*)malloc(strlen(header[num]) + 3);
	sprintf(c, "%s_r", header[num]);
	addSeq(seq[num], c);
	setCoverage(getHeaderIndex(c), coverage[num]);
	free(c);
}

unsigned int fastaReader::findDuplicateOrigin(unsigned int num)
{
	char *c;
	int len, ret;
	c = header[num];
	len = strlen(c);
	c[len - 2] = '\0';
	ret = getHeaderIndex(c);
	c[len - 2] = '_';
	return ret;
}

float fastaReader::getCoverage(char *input_header)
{
	char *c, *c2;
	float ret = (float)0;
	c = input_header;
	while (!(*c == 'c' && *(c + 1) == 'v' && *(c + 2) == 'g' && *(c + 3) == '_') && *c != '\0')
	{
		c++;
	}
	if (*c == '\0')
	{
		return ret;
	}
	c = c + 4;
	c2 = c + 1;
	while(*c2 != '_')
	{
		c2++;
	}
	*c2 = '\0';
	ret = (float)atof(c);
	*c2 = '_';

	return ret;
}

void fastaReader::setCoverage(unsigned int num, float input_cov)
{
	coverage[num] = input_cov;
}

unsigned int fastaReader::getCharSum()
{
	return charSum;
}
