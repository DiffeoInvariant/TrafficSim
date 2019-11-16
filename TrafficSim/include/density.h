#ifndef TS_DENSITY_H
#define TS_DENSITY_H

#include "highway.h"
#include <petscsys.h> /* primarily for PetscMath.h */



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

  TSExponentialDistribution*                     dist;
  /* pointer to function whose first parameter is lambda, second parameter
     is time, and has possibly some other parameters, and computes the desired
     value of lambda at the given time. */
  (PetscErrorCode)(PetscReal*, PetscReal, ...)*  lambda_generator;
};
  

typedef struct _ts_DynamicExponentialDistribution TSDynamicExponentialDistribution;

#endif
