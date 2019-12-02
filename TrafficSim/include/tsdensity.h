#ifndef TS_DENSITY_H
#define TS_DENSITY_H

#include <petscsys.h> /* primarily for PetscMath.h */

typedef enum {
	      TS_NO_ARRIVAL,
	      TS_POISSON_STATIC,/* Poisson dist w/ time-independent rate */
	      TS_POISSON_DYNAMIC, /* Poisson dist w/ time-dependent rate (e.g. a Hidden Markov Model where the transition densities given the hidden state are Poisson) */
	      TS_POISSON_STATIC_MIXTURE, /* static mixture of independent Poisson distributions */
	      TS_POISSON_DYNAMIC_MIXTURE /* HMM with conditional densities being a Poisson mixture */
} TSArrivalDistributionType;

typedef struct {
  PetscReal lambda;
  PetscReal *samples;
  PetscInt  nsamples=0;
} TSExponentialDistribution;

typedef struct {
  PetscReal lambda;
  PetscInt  *samples;
  PetscInt  nsamples=0;
} TSPoissonDistribution;


extern PetscErrorCode TSExponentialSample(TSExponentialDistribution*, PetscInt);

extern PetscErrorCode TSPoissonSample(TSPoissonDistribution*, PetscInt);


struct _ts_DynamicExponentialDistribution {

  PetscInt                                       num_timepoints=0;

  PetscInt*                                      timepoints;

  /* array of exponential distributions, one for each timepoint */
  TSExponentialDistribution*                     dists;
  /* pointer to function whose first parameter is lambda, second parameter
     is time, and has possibly some other parameters, and computes the desired
     value of lambda at the given time. */
  PetscErrorCode (*lambda_gen)(PetscReal*, PetscReal, ...);
};


typedef struct _ts_DynamicExponentialDistribution TSDynamicExponentialDistribution;


typedef enum {TS_LINEAR, TS_GENERALIZED_LOGISTIC, TS_DEFAULT,TS_DECIDE} TSSpeedDensityModel;

extern PetscErrorCode TSSpeedFromLinearDensity(PetscReal*, PetscReal, PetscReal);

extern PetscErrorCode TSSpeedFromLogisticDensity(PetscReal*, PetscReal, PetscReal);

#endif
