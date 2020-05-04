/* Take a nucleotide sequence and spit out the 6 translated open reading frames as protein sequences.
*/

#ifndef eslTRANSLATE_INCLUDED
#define eslTRANSLATE_INCLUDED

#include "esl_sq.h"

extern int esl_trans_orf(ESL_SQ *in, ESL_SQ ***out, int *retSeq, int minimumResidues);
extern int esl_trans_6frame         (ESL_SQ *in, ESL_SQ **out);
extern int esl_trans_s2p            (ESL_SQ *in, ESL_SQ **out, int frameshift, int rcFlag);
extern int esl_trans_seq_stop_split (ESL_SQ *in, ESL_SQ ***out, int *outCount);

#endif
