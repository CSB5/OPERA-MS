/* Easel's portable, threadsafe random number generator.
 */
#ifndef eslRANDOM_INCLUDED
#define eslRANDOM_INCLUDED

#define eslRND_FAST     0
#define eslRND_MERSENNE 1

typedef struct {
  int      type;		/* eslRND_FAST | eslRND_MERSENNE               */
  int      mti;			/* current position in mt[] table              */
  uint32_t mt[624];		/* state of the Mersenne Twister               */
  uint32_t x;			/* state of the Knuth generator                */
  uint32_t seed;		/* seed used to init the RNG                   */
} ESL_RANDOMNESS;

/* esl_rnd_Roll(a) chooses a uniformly distributed integer
 * in the range 0..a-1, given an initialized ESL_RANDOMNESS r.
 */
#define esl_rnd_Roll(r, a)    ((int) (esl_random(r) * (a)))

/* 1. The ESL_RANDOMNESS object.
 */
extern ESL_RANDOMNESS *esl_randomness_Create(uint32_t seed);
extern ESL_RANDOMNESS *esl_randomness_CreateFast(uint32_t seed);
extern ESL_RANDOMNESS *esl_randomness_CreateTimeseeded(void); /* DEPRECATED */
extern void            esl_randomness_Destroy(ESL_RANDOMNESS *r);
extern int             esl_randomness_Init(ESL_RANDOMNESS *r, uint32_t seed);
extern uint32_t        esl_randomness_GetSeed(const ESL_RANDOMNESS *r);

/* 2. The generator, esl_random().
 */
extern double esl_random(ESL_RANDOMNESS *r);

/* 3. Debugging/development tools.
 */
extern int esl_randomness_Dump(FILE *fp, ESL_RANDOMNESS *r);

/* 4. Other fundamental sampling (including Gaussian, gamma).
 */
extern double esl_rnd_UniformPositive(ESL_RANDOMNESS *r);
extern double esl_rnd_Gaussian(ESL_RANDOMNESS *r, double mean, double stddev);
extern double esl_rnd_Gamma(ESL_RANDOMNESS *r, double a);

/* 5. Multinomial sampling from discrete probability n-vectors.
 */
extern int    esl_rnd_DChoose   (ESL_RANDOMNESS *r, const double *p,   int N);
extern int    esl_rnd_FChoose   (ESL_RANDOMNESS *r, const float  *p,   int N);
extern int    esl_rnd_DChooseCDF(ESL_RANDOMNESS *r, const double *cdf, int N);
extern int    esl_rnd_FChooseCDF(ESL_RANDOMNESS *r, const float  *cdf, int N);


#endif /*eslRANDOM_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.1b1; May 2013
 * Copyright (C) 2013 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
