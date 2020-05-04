#ifndef __FASTAREADER_H__
#define __FASTAREADER_H__

#include "global_inc.h"
#include <map>

class fastaReader
{
	public:
		fastaReader(char *inputFile);
		fastaReader();
		~fastaReader();
		bool isFasta();
		bool isMark(int i);
		unsigned int getNum();
		char* getNext();
		char *getNextMarked();
		char* getSeqByNum(unsigned int num);
		char* getFullSeqByNum(unsigned int num);
		char* getHeaderByNum(unsigned int num);
		float getCoverageByNum(unsigned int num);
		char* getCurrentHeader();
		unsigned int getSeqLenByNum(unsigned int num);
		unsigned int getFullSeqLenByNum(unsigned int num);
		unsigned int getMaxLen();
		char* getPartialSeq(unsigned int num, unsigned int pos_start, unsigned int pos_end);
		void resetCurrent();
		void reverseSeq(unsigned int num);
		void setMark(unsigned int num);
		void setMark(char *id);
		void resetMark();
		void setOffset(unsigned int num, unsigned int input_offset);
		void splitMarkedSeq(fastaReader **marked, fastaReader **unmarked);
		void addSeq(const char *seq, const char *header);
		void removeSeq(unsigned int num);
		void separateSeq(unsigned int num, unsigned int cut1, unsigned int cut2);
		fastaReader* getMarkedSeq();
		fastaReader* getUnmarkedSeq();
		int getHeaderIndex(char *header);
		unsigned int getCharSum();
		void duplicateSeqByNum(unsigned int num);
		unsigned int findDuplicateOrigin(unsigned int num);

	private:
		// Structure
		struct str_cmp
		{
			bool operator() (char* const a, char* const b) const
			{
				return (strcmp(a, b) < 0) ? true : false;
			}
		};

		// Variables
		static const int ADD_SIZE;
		bool is_fasta;
		unsigned int currentPos;
		fstream *fs;
		char **header;
		char **seq;
		bool *mark;
		float *coverage;
		unsigned int *offset;
		unsigned int *seqLen;
		unsigned int seqNum;
		unsigned int maxNum;
		unsigned int charSum;
		char str[10240];
		map<char*, int, str_cmp> headerHash;
		// Functions
		void init();
		void parse();
		char revChar(char ch);
		void convertToUpper(char *seq);
		void trimHeader(char *input_header, char div);
		float getCoverage(char *input_header);
		void setCoverage(unsigned int num, float input_cov);
};


#endif

