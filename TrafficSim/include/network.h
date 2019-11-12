#ifndef TNET_H
#define TNET_H

#include <petscsnes.h>
#include <petscdmnetwork.h>
#include <mpi.h>
#include "highway.h"

typedef enum {
	      TS_POISSON_STATIC,/* Poisson dist w/ time-independent rate */
	      TS_POISSON_DYNAMIC, /* Poisson dist w/ time-dependent rate (e.g. a Hidden Markov Model where the transition densities given the hidden state are Poisson) */
	      TS_POISSON_STATIC_MIXTURE, /* static mixture of independent Poisson distributions */
	      TS_POISSON_DYNAMIC_MIXTURE /* HMM with conditional densities being a Poisson mixture */
} TSArrivalDistributionType;

extern PetscErrorCode TrafficNetworkCreate(MPI_Comm comm, DM* network_dm);


#endif
