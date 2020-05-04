/* Take a nucleotide sequence and spit out the 6 translated open reading frames as protein sequences.
 *
 * WMA March 29, 2013  
 *  
 * Contents:
 *   1.  esl_trans_6frames - turn a nucleotide sequence into 6 protein sequences, one for each reading frame
 *   2.  esl_trans_s2p - turn a nucleotide sequence into its corresponding protein sequence
 *   3.  esl_trans_seq_stop_split - split a protein sequence which contains stop codons
 *   4.  practice driver
 */
   
#include "esl_config.h"

#include <stdio.h>
#include <string.h>
    
#include "easel.h"
#include "esl_translate.h"
#include "esl_sq.h"
#include "esl_sqio.h"

/*****************************************************************
 * 1. generating the 6 reading frames for a nucleo sequence
 *****************************************************************/

/* Function: esl_trans_6frames()
 * Synopsis: Get the 6 translated open reading frames from a nucleotide sequence.
 *
 * Args:     ESL_SQ *in     the input dna/rna sequence
 *           ESL_SQ **out   a pointer to a caller created array with 6 items of type ESL_SQ*
 *                          where the resulting protein sequences will be placed
 *
 * Returns:  <eslOK> if it worked, <eslEMEM> if some allocation or copy failed
 *
 */

 int esl_trans_6frame(ESL_SQ *in, ESL_SQ **out)
 {
   int x;
   
   for(x = 0; x < 3; x++)
   {
     if(eslOK != esl_trans_s2p(in, &out[x], x, 0)) return eslFAIL;
   }
     
   for(x = 3; x < 6; x++)
   {
     if(eslOK != esl_trans_s2p(in, &out[x], x-3, 1)) return eslFAIL;
   }
   
   return eslOK;
 }

 
 /* Function: esl_trans_orf()
  * Synopsis: Get the open reading frames from a nucleotide sequence.
  *
  * Args:     ESL_SQ *in      the input dna/rna sequence
  *           ESL_SQ ***out   a pointer to a new array holding a sequence for each open
  *                           reading frame of sufficient size
  *           int* retSeq     location where size of returned array will be placed
  *
  * Returns:  <eslOK> if it worked, <eslEMEM> if some allocation or copy failed
  *
  * Note:   This function is tightly coupled with esl_trans_seq_stop_split()
  *         Specifically, the format expected for the names is parsed and
  *         modified to convert from residue index numbers to nucleo
  *         index from the original .fa file.
  */
 
 int esl_trans_orf(ESL_SQ *in, ESL_SQ ***out, int *retSeq, int minimumResidues)
 {
   int x, y;           //loop counters
   int status;         //ESL_ALLOC required
   
   ESL_SQ *frame;    //hold the current reading frames
   ESL_SQ **splits[6];   //hold the open reading frame shards when split by stop codons
   int stops[6];         //hold the number of shards per each reading frame
   int totalStops = 0;   //hold the total number of shards among all frames to be output
   int currentWrite = 0; //the next location to place a new reading frame shard into the output
   
   char workbench[256];  //used to build new sequence names
   
   char *toStart, *fromStart;  //pointer to beginning of residue index in sequence name
   int toNum, fromNum;         //residue index extracted from sequence name
   
   if(minimumResidues <= 0) minimumResidues = 15;  //default minimum shard size
   
   if(out == NULL) return eslFAIL;
   
   for(x = 0; x < 3; x++)
   {
     //get one of the three forward reading frames
     if((status = esl_trans_s2p(in, &frame, x, 0)) != eslOK) return status; 
     //split the frame into shards between stop codons
     if((status = esl_trans_seq_stop_split(frame, &(splits[x]), &stops[x])) != eslOK) return status;
     esl_sq_Destroy(frame);
     for(y = 0; y < stops[x]; y++)
     {
       //extract the to/from out of the back of the sq name
       //mess hack gross, but we want nucleotide indexes not residue
       //expect "NAME_XXXtoYYY" where x is the first residue index in the shard and y is the last
       toStart = splits[x][y]->name; 
       while(*toStart != '\0') toStart+=1;
       while(*toStart != 'o') toStart-=1;
       toStart+=1;
       toNum = atoi(toStart);
       fromStart = toStart;
       while(*fromStart != '_') fromStart-=1;
       fromStart+=1;
       toStart[-2]='\0';
       fromNum = atoi(fromStart);
       toStart[-2]='t';
            
       *fromStart = '\0';
       strcpy(workbench, splits[x][y]->name);
       //substitute residue index numbers for nucleotide ones
       sprintf(splits[x][y]->name, "%s%dto%d", workbench, x+1+3*(fromNum-1), x+3*(toNum));
       }
   }
   for(x = 0; x < 3; x++)
   {
     //same as the last for loop, except reverse compliment the inputs
     if((status = esl_trans_s2p(in, &frame, x, 1)) != eslOK) return status;
     if((status = esl_trans_seq_stop_split(frame, &(splits[x+3]), &stops[x+3])) != eslOK) return status;
     esl_sq_Destroy(frame);
     for(y = 0; y < stops[x+3]; y++)
     {
       //extract the to/from out of the back of the sq name
       //mess hack gross, but we want nucleotide indexes not residue
       toStart = splits[x+3][y]->name; 
       while(*toStart != '\0') toStart+=1;
       while(*toStart != 'o') toStart-=1;
       toStart+=1;
       toNum = atoi(toStart);
       fromStart = toStart;
       while(*fromStart != '_') fromStart-=1;
       fromStart+=1;
       toStart[-2]='\0';
       fromNum = atoi(fromStart);
       toStart[-2]='t';
              
       *fromStart = '\0';
       strcpy(workbench, splits[x+3][y]->name);
       sprintf(splits[x+3][y]->name, "%s%dto%d", workbench, (int)in->n-(x+1+3*(fromNum-1)), (int)in->n-(x+3*(toNum)));
       }
   }
   
   //count the number of shards which are at least the minimum size
   for(x = 0; x < 6; x++)
   {
     for(y = 0; y < stops[x]; y++)
     {
       if(splits[x][y]->n >= minimumResidues) totalStops++;
     }
   }
   
   ESL_ALLOC(*out, totalStops * sizeof(ESL_SQ*));
      
   //copy shards of minimum size into newly allocated output array
   //or erase the ones that are too small
   for(x = 0; x < 6; x++)
   {
     for(y = 0; y < stops[x]; y++)
     {
       if(splits[x][y]->n >= minimumResidues)
       {
         (*out)[currentWrite++] = splits[x][y];
       }
       else
       {
         esl_sq_Destroy(splits[x][y]);
       }
     }
     free(splits[x]);
   }
   (*retSeq) = totalStops;

   if(*retSeq) return eslOK;

   ERROR:

   if(frame) esl_sq_Destroy(frame);

   for(x = 0; x < 6; x++)
   {
     for(y = 0; y < stops[x]; y++)
     {
       if(splits[x][y]) esl_sq_Destroy(splits[x][y]);
     }
     free(splits[x]);
   }
   
   *retSeq = 0;

   return eslFAIL;
 }

/*****************************************************************
 * 2. generating the reading frames for a nucleo sequence
 *****************************************************************/

/* Function: esl_trans_s2p()
 * Synopsis: Turn a nucleotide sequence into a protein sequence.
 *
 * Purpose:  Input is a dna/rna sequence.  Read codons and output a new sequence in the protein alphabet.
 *
 * Args:     ESL_SQ *in      a dna/rna sequence
 *           ESL_SQ **out    the location where the new protein sequence will be output
 *           int frameshift  The number of bases shifted off the front to cause a shift in the reading frame
 *                           Will fail if you shift further than sequence available
 *           int rcFlag      True if you want a reverse compliment reading frame, false if you do not
 *
 * Returns:  <eslOK>   on success, out is valid
 *           <eslEMEM> an allocation failed
 *           <eslFAIL> frameshifted too much
 *   
 * Note:     If this function throws a bad return when rc is true, the input sequence
 *           will be left in an altered state.
 */

int esl_trans_s2p(ESL_SQ *in, ESL_SQ **out, int frameshift, int rcFlag)
{
  // The encoding for this is taken from squid:  A=0, C=1, G=2, U/T=3, 
  // code[0] corresponds to AAA, code[1] is AAC... code[4] is ACA... 
  // and so on up to 63 being UUU. 64 is a sentinel. Regular 20 amino codes and '*' for stop
  // the nucleotide indices match well with the easel alphabet index
  // but the actual translation still needs to be hard coded
  char code[] = {'K','N','K','N','T','T','T','T','R','S','R','S',
                 'I','I','M','I','Q','H','Q','H','P','P','P','P',
                 'R','R','R','R','L','L','L','L','E','D','E','D',
                 'A','A','A','A','G','G','G','G','V','V','V','V',
                 '*','Y','*','Y','L','F','L','F','*','C','W','C',
                 'L','F','L','F'};

  int status;

  int codon;     //progress in counting current codon
  char *aaseq;   //hold the protein sequence to be output
  char *aaptr;   //pointer records progress in writing to output
  char *readseq; //pointer records progress in reading nucleotide sequence
  int read_dg;   //index into digital sequence
  
  ESL_ALPHABET *abc = esl_alphabet_Create(eslDNA);
  char errbuf[256]; //validateseq demands this
  
  char namestring[256];
  
  (*out) = NULL;

  if(frameshift >= in->n) return eslFAIL;
  if(!abc) goto ERROR;
  
  //make sure we have a nucleotide sequence; could use esl_abc_ValidateSeq but that wants too
  //much boilerplate for the simple bit I need done. doesn't help that i don't care if there are U or T
  //characters but that would test against two alphabets
  if(in->seq)
  {
    if(eslOK != esl_abc_ValidateSeq(abc, in->seq, in->n, errbuf)) goto ERROR;
  }
  else if(in->dsq)
  {
    if(in->abc->type != eslRNA && in->abc->type != eslDNA) goto ERROR;
  }
  else
  {
    goto ERROR;
  }

  
  //apply the reverse compliment
  if(rcFlag) {if(esl_sq_ReverseComplement(in) != eslOK) goto ERROR;}
  
  
  ESL_ALLOC(aaseq, (in->n+1) * sizeof(char));
  aaptr = aaseq;
  
  if(in->seq) //text sequence
  { 
    //get an alphabet to do the lookup with.
    //an ordinary text sequence doesn't have in->abc
    //if it has one that is not a standard dna/rna alphabet
    //then this code won't work. I wanted to use an alphabet if available, could save some allocating time that way
    //if we're calling this repeatedly
    //but the compiler complains about "pointer qualifiers" so nevermind
    
    readseq = in->seq+frameshift;
      
    //as long as there are at least 3 nucleotides left, pull and translate another codon
    for (; *readseq != '\0' && *(readseq+1) != '\0' && *(readseq+2) != '\0'; readseq += 3)
    {
      codon = abc->inmap[(int)*(readseq)] * 16 + abc->inmap[(int)*(readseq+1)] * 4 + abc->inmap[(int)*(readseq+2)];
      if(codon > 63 || codon < 0) break;

      *aaptr = code[codon];
      aaptr += 1;
    }
    *aaptr = '\0';
  }
  else if(in->dsq)  //do it digitally
  { 
    if(in->dsq == NULL) goto ERROR;
    
    read_dg = 1+frameshift; //add one here because digital index 0 is a sentinel
    for(;in->dsq[read_dg] != 255 && in->dsq[read_dg+1] != 255 && in->dsq[read_dg+2] != 255; read_dg += 3)
    {
      codon = in->dsq[read_dg] * 16 + in->dsq[read_dg+1] * 4 + in->dsq[read_dg+2];
      if(codon > 63 || codon < 0) break;
      *aaptr = code[codon];
      aaptr += 1;
    }
    *aaptr = '\0';
  }
  else
  {
    goto ERROR;
  }
  
  //modify name to record any reading frame adjustments
  sprintf(namestring, "%s_s%d", in->name, frameshift);
  if(rcFlag) strcat(namestring, "_rc");
  *out = esl_sq_CreateFrom(namestring, aaseq, in->desc, in->acc, in->ss);
        
  if(aaseq != NULL) free(aaseq);
  
  //return the input to its original state
  if(rcFlag) {if(esl_sq_ReverseComplement(in) != eslOK) goto ERROR;}
  
  if(abc) esl_alphabet_Destroy(abc);
  if(*out) return eslOK;
  
  ERROR:
    
  if(abc) esl_alphabet_Destroy(abc);
  if(aaseq != NULL) free(aaseq);
  (*out) = NULL;
  
  return eslEMEM;
}

/*************************************************************************
 * 3. Divide a protein sequence with stop codons into multiple sequences
 *************************************************************************/

/* Function: esl_trans_seq_stop_split()
 * Synopsis: Split a protein sequence with stop codons in it.
 *
 * Purpose:  Input is a protein sequence.  Look for stop codon marks and output an array listing all
 *           fragments separated by a stop codon.
 *
 * Args:     ESL_SQ *in      a protein sequence
 *           ESL_SQ ***out   the location where an array of new protein sequences will be placed.
 *                           caller is responsible for freeing this array when done with it.
 *           int *outcount   Address of an integer.  This is where the size of the returned array will be placed.
 *                          
 * Returns:  <eslOK>   on success, out is valid
 *           <eslFAIL> something wrong. don't trust out
 *   
 * Note:     This function is tightly coupled with esl_trans_orf().
 *           If the modification to in->name is changed then the corresponding
 *           name parse in trans_orf must be modified as well.
 */

int esl_trans_seq_stop_split(ESL_SQ *in, ESL_SQ ***out, int *outCount)
{
  int status;
  
  int x, y;          //loop counters
  int nextSeqOut;    //index of the next open location in the output sequence array
  int front;         //front of the segment of sequence currently being read
  
  char* buff;        //temporary home of output sequence before calling createFrom
  char name[256];    //workbench for building the name of each output sequence
  
  ESL_ALLOC(buff, (in->n+1) * sizeof(char));
  
  *outCount = 1;
  
  if(in->seq) //text mode
  {
    //count how many sequences are present. minimum size is one non-stop residue
    x = 1;
    while(in->seq[x] != '\0')
    {
      if(in->seq[x] == '*' && in->seq[x-1] != '*') (*outCount)++;
      x++;
    }
    
    ESL_ALLOC(*out, sizeof(ESL_SQ*) * *outCount);
    
    x = front = 0;
    nextSeqOut = 0;
    
    //continue until the sequence front steps past the end of the list
    while(front < in->n)
    {
      //x is the location currently being read, current segment is from front to x
      x++;
      if(in->seq[x] == '\0' || in->seq[x] == '*') //if we see something that ends a segment
      {
        if(x - front > 0) //if there is at least one residue
        {
          //build name
          sprintf(name, "%s_%dto%d", in->name, front+1, x);
          
          //build temporary sequence string
          strncpy(buff, in->seq+front, x-front);
          buff[x-front] = '\0';
          
          //load output array
          (*out)[nextSeqOut++] = esl_sq_CreateFrom(name, buff, in->desc, in->acc, in->ss);
        }
        //step the front to the beginning of the next sequence
        front = x+1;
      }
    }
  }
  else if(in->dsq) //digital mode
  {
    //start a little different because dsq has a sentinel in position 0
    x = 2;
    while(in->dsq[x] != 255) //until the end sentinal, count sequences with at least one residue
    {
      if(in->abc->inmap[(int)'*'] == in->dsq[x] && in->abc->inmap[(int)'*'] != in->dsq[x-1]) (*outCount)++;
      x++;
    }

    ESL_ALLOC(*out, sizeof(ESL_SQ*) * *outCount);
    
    x = front = 1;    
    nextSeqOut = 0;
    
    while(front < in->n+2) //as long as we have residues left
    {
      x++;
      if(in->dsq[x] == 255 || in->abc->inmap[(int)'*'] == in->dsq[x]) //if we see something that finishes a sequence
      {
        if(x - front > 0) //have at least one residue in the sequence
        {
          //build name
          sprintf(name, "%s_%dto%d", in->name, front, x-1);
          
          //build temporary sequence
          for(y = 0; y < x-front; y++) buff[y] = in->abc->sym[in->dsq[front+y]];
          buff[x-front] = '\0';
          
          //load output
          (*out)[nextSeqOut++] = esl_sq_CreateFrom(name, buff, in->desc, in->acc, in->ss);
        }
        front = x+1;
      }
    }
  }
  else
  {
    goto ERROR;
  }
  
  free(buff);
  
  return eslOK;
  
  ERROR:

  if(buff) free(buff);

  return eslFAIL;
}

/*------------ end, esl_translate --------------------*/

/*****************************************************************
 * 4. practice driver.
 *****************************************************************/ 
#ifdef eslTRANSLATE_TESTDRIVE
/* gcc -g -Wall -o esl_translate_utest -L. -I. -D eslTRANSLATE_TESTDRIVE esl_translate.c -leasel -lm
 */
   
#include "esl_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_translate.h"
#include "esl_sq.h"
#include "esl_sqio.h"

int main(int argc, char **argv)
{
  ESL_SQFILE        *sqfp   = NULL;      
  ESL_SQ            *sq   = NULL;
  ESL_SQ            *dsq  = NULL;
  ESL_SQ            **prot;
  int               c;
  ESL_ALPHABET      *abc, *prot_abc;
  
  ESL_SQ            *prot6[6];
    
  int x;
  
  abc = esl_alphabet_Create(eslDNA);
  prot_abc = esl_alphabet_Create(eslAMINO);
  
  if(argc != 2)
  {
    printf("You need to pass an argument for a filepath to a dna/rna fasta file\n");
    exit(0);
  }
    
  if(eslOK != esl_sqfile_Open(argv[1], eslSQFILE_FASTA, NULL, &sqfp)) 
  {
    printf("Invalid filepath: %s\n", argv[1]);
    exit(0);
  }
  
  sq = esl_sq_Create();
  if(sq == NULL)
  {
    printf("could not allocate new sequence\n");
    exit(0);
  }
  
  if(esl_sqio_Read(sqfp, sq) != eslOK)
  {
    printf("Not a valid fasta file %s\n", argv[1]);
    exit(0);
  }

  dsq = esl_sq_Create();
  if(dsq == NULL)
  {
    printf("could not allocate digital sequence\n");
    exit(0);
  }
  
  if(esl_sq_Copy(sq, dsq) != eslOK)
  {
    printf("could not copy sequence\n");
    exit(0);
  }
  
  if(esl_sq_Digitize(abc, dsq) != eslOK)
  {
    printf("could not digitize sequence\n");
    exit(0);
  }

  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, 0);
  
  if(esl_trans_6frame(sq, prot6) != eslOK)
  {
    printf("could not generate six frame translation\n");
    exit(0);
  }
   
  for(x = 0; x < 6; x++)
  {
    esl_sqio_Write(stdout, prot6[x], eslSQFILE_FASTA, 0);
  }
    
  if(esl_trans_orf(dsq, &prot, &c, 10) != eslOK)
  {
    printf("could not translate open reading frames\n");
    exit(0);
  }
  
  for(x = 0; x < c; x++)
  {
    esl_sqio_Write(stdout, prot[x], eslSQFILE_FASTA, 0);
  }

  return 0;
}

#endif /*eslTRANSLATE_TESTDRIVE*/
/*----------------- end, practice driver ----------------------------*/
