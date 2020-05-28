#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>


double RandomReal(double LowerLimit, double UpperLimit){
    const gsl_rng * r;
    gsl_rng_env_steup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    struct timeval tv;
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    gsl_rng_set(r, mySeed);

    double u = gsl_rng_uniform(r);
    double result = (UpperLimit - LowerLimit) * u + LowerLimit
    gsl_rng_free(r);

    return(result);
}
