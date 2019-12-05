#ifndef TNET_H
#define TNET_H

#include <petscsnes.h>
#include <petscdmnetwork.h>
#include <petscts.h>
#include <mpi.h>
#include "tshighway.h"

typedef enum {TS_LINEAR_STEADY_STATE, TS_NONLINEAR_STEADY_STATE, TS_TIME_STEP} TSProblemType;

struct _ts_HighwayVertex {
  /* set these to NULL if not applicable for this vertex */
  TSInterchangeCtx          *intc_ctx=NULL;
  TSInterchangeCtx          *intc2_ctx=NULL; /* if there's a split/join */
  TSHighwayEntryCtx         *entr_ctx=NULL;
  TSHighwayExitCtx          *exit_ctx=NULL;
  Mat                       *jac=NULL;
  /* Can people enter or exit the highway system here? 
     If so, what distribution does that follow? */
  TSArrivalDistributionType entry_dist=TS_NO_ARRIVAL;
  TSExitType                exit_dist=TS_NO_EXIT;
  
  PetscInt                  id;
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

  PetscInt        highway_key, vertex_key;

  PetscInt*       edgelist;/*local edge list */
  PetscInt        g_discrete_dimension, l_discrete_dimension; /* number of global and local
							   nodes used in DMDA discretization */
  Vec             g_X, l_X, l_dXdt; /* vectors used for local function evaluation at nodes (X is (rho, v)) */

  TSHighway       highways;
  TSHighwayVertex vertices;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_Network *TSNetwork;

extern PetscErrorCode TSNetworkCreate(MPI_Comm,TSNetwork*);

extern PetscErrorCode TSNetworkCreateWithStructure(MPI_Comm, TSNetwork*, DM, PetscInt, const char*);

extern PetscErrorCode TSNetworkDistribute(MPI_Comm, TSNetwork);

extern PetscErrorCode TSNetworkDestroy(TSNetwork);

extern PetscErrorCode TSNetworkSetHighways(TSNetwork,PetscInt,TSHighway,
						PetscInt,TSHighwayVertex,
						PetscInt[]);

extern PetscErrorCode TSNetworkSetUp(TSNetwork);

extern PetscErrorCode TSNetworkSetSNES(TSNetwork, SNES);

extern PetscErrorCode TSNetworkGetSNES(TSNetwork, SNES);

extern PetscErrorCode TSNetworkSetTS(TSNetwork, TS);

extern PetscErrorCode TSNetworkGetTS(TSNetwork, TS);

extern PetscErrorCode TSNetworkGetDM(TSNetwork, DM);

extern PetscErrorCode TSNetworkSetDM(TSNetwork, DM);

extern PetscErrorCode TSNetworkSetMaxTime(TSNetwork, PetscReal);

extern PetscErrorCode TSNetworkGetMaxTime(TSNetwork, PetscReal*);

extern PetscErrorCode TSNetworkSetMaxTimeSteps(TSNetwork, PetscInt);

extern PetscErrorCode TSNetworkGetMaxTimeSteps(TSNetwork, PetscInt*);

extern PetscErrorCode TSNetworkSetProblemType(TSNetwork, TSProblemType);

extern PetscErrorCode TSNetworkSolve(TSNetwork, Vec);

extern PetscErrorCode TSNetworkSetFromOptions(TSNetwork);



#endif
