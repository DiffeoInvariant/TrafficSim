#ifndef TNET_H
#define TNET_H

#include <petscsnes.h>
#include <petscdmnetwork.h>
#include <mpi.h>
#include "highway.h"

typedef enum {
	      TS_NO_ARRIVAL,
	      TS_POISSON_STATIC,/* Poisson dist w/ time-independent rate */
	      TS_POISSON_DYNAMIC, /* Poisson dist w/ time-dependent rate (e.g. a Hidden Markov Model where the transition densities given the hidden state are Poisson) */
	      TS_POISSON_STATIC_MIXTURE, /* static mixture of independent Poisson distributions */
	      TS_POISSON_DYNAMIC_MIXTURE /* HMM with conditional densities being a Poisson mixture */
} TSArrivalDistributionType;

struct _ts_HighwayVertex {
  /* set these to NULL if not applicable for this vertex */
  TSInterchangeCtx          *intc_ctx;
  TSHighwayEntryCtx         *entr_ctx;
  TSHighwayExitCtx          *exit_ctx;
  Mat                       *jac;
  /* Can people enter or exit the highway system here? 
     If so, what distribution does that follow? */
  TSArrivalDistributionType entry_dist=TS_NO_ARRIVAL;
  TSExitType                exit_dist=TS_NO_EXIT;

  PetscInt                  is_boundary=0; /* 0 = not a boundary vertex,
					      1 = boundary vertex, cars can
					      enter here, but no travel
					      in the other direction,
					      2 = boundary vertex, cars can
					      exit but not enter.
					      3 = boundary vertex, cars can
					      both enter and exit */
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayVertex *TSHighwayVertex;

struct _ts_Network {
  DM              network;
  MPI_comm        comm;
  PetscInt        l_nedges, l_nvertices; /* local num of edges and vertices */
  PetscInt        g_nedges, g_nvertices; /* global num of edges and vertices */

  PetscInt*       edgelist;/*local edge list */
  PetscInt        g_discrete_dimension, l_discrete_dimension; /* number of global and local
							   nodes used in DMDA discretization */
  Vec             l_X, l_dXdt; /* vectors used for local function evaluation at nodes (X is (rho, v)) */

  TSHighway       highways;
  TSHighwayVertex vertices;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_Network *TSNetwork;

extern PetscErrorCode TrafficNetworkCreate(TSNetwork*);


#endif
