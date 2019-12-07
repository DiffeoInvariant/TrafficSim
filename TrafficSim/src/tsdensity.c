#include "../include/tsdensity.h"
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <petscsys.h>


PetscErrorCode TSUniformSample(PetscReal *x, PetscReal low, PetscReal high)
{
  PetscFunctionBeginUser;
        *x = (double)rand() / nextafter((double)RAND_MAX, DBL_MAX);

	if(PetscUnlikely(low != 0.0)){
		*x += low;
	}
	if(PetscUnlikely(high - low != 1.0)){
		*x *= (high - low);
	}

  PetscFunctionReturn(0);
}

PetscErrorCode TSPoissonSample(PetscReal lambda, PetscInt *sample)
{
  PetscErrorCode ierr;
  PetscReal      lim, tot, uspl;
  PetscInt       spl;

  PetscFunctionBegin;

  spl = 0;
  lim = PetscExpReal(-lambda);
  ierr = TSUniformSample(&tot, 0.0, 1.0);CHKERRQ(ierr);
  while(tot > lim){
    ++spl;
    ierr = TSUniformSample(&uspl, 0.0, 1.0);CHKERRQ(ierr);
    tot *= uspl;
  }

  *sample = spl;

  PetscFunctionReturn(0);
}

PetscErrorCode TSSpeedFromLinearDensity(PetscReal* speed, PetscReal rho, PetscReal speed_limit, PetscReal rho_limit)
{
  PetscFunctionBegin;
  *speed = speed_limit * (1 - rho/rho_limit);
  PetscFunctionReturn(0);
}

  
  
  
