#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
double RandomReal(double LowerLimit, double UpperLimit)
{
 const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default; //T is the type of random number generator.
  r = gsl_rng_alloc (T);

    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
 

  gsl_rng_set(r, mySeed); // Seed with time
  double u = gsl_rng_uniform (r);
  double RandomRealResult = (UpperLimit - LowerLimit)*u + LowerLimit;
  gsl_rng_free (r);
  return(RandomRealResult);
}