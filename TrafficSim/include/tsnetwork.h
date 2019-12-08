#ifndef TNET_H
#define TNET_H

#include <petscsnes.h>
#include <petscdmnetwork.h>
#include <petscts.h>
#include <mpi.h>
#include "tshighway.h"

/*typedef enum {TS_LINEAR_STEADY_STATE, TS_NONLINEAR_STEADY_STATE, TS_TIME_STEP} TSProblemType;*/

struct _ts_HighwayVertex {
  /* set these to NULL if not applicable for this vertex */
  TSInterchangeCtx          *intc_ctx;/*=NULL;*/
  TSInterchangeCtx          *intc2_ctx;/*=NULL;  if there's a split/join */
  TSHighwayEntryCtx         *entr_ctx;/*=NULL;*/
  TSHighwayExitCtx          *exit_ctx;/*=NULL;*/

  PetscReal                 rho, v, speed_limit, rho_limit;
  Mat                       *jac;/*=NULL;*/
  /* Can people enter or exit the highway system here? 
     If so, what distribution does that follow? */
  TSArrivalDistributionType arrival_dist;/*=TS_NO_ARRIVAL;*/
  TSExitType                exit_dist;/*=TS_NO_EXIT;*/
  
  PetscInt                  id;
  PetscInt                  is_boundary;/*=0; 0 = not a boundary vertex,
					      1 = boundary vertex, cars can
					      enter here, but no travel
					      in the other direction,
					      2 = boundary vertex, cars can
					      exit but not enter.
					      3 = boundary vertex, cars can
					      both enter and exit */
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayVertex *TSHighwayVertex;


extern PetscErrorCode VertexCreateJacobian(DM dm, TSHighwayVertex vertex, PetscInt v, Mat* Jpre, Mat *J[]);

extern PetscErrorCode VertexGetArrivalRate(TSHighwayEntryCtx* vertex_ctx, PetscReal t, PetscReal* lambda);

extern PetscErrorCode VertexGetExitRate(TSHighwayExitCtx* vertex_ctx, PetscReal t, PetscReal* erate);

extern PetscErrorCode VertexGetNumArrivals(TSHighwayEntryCtx* vertex_ctx, PetscReal t, PetscInt* narrival);

extern PetscErrorCode VertexGetTrafficFlowRate(PetscReal rho, PetscReal v, PetscReal* tflowrate);

extern PetscErrorCode VertexGetCurrentTrafficCount(PetscReal rho, PetscReal v, PetscReal deltat, PetscInt* count);

extern PetscErrorCode VertexGetDensityFactor(TSHighwayEntryCtx* vertex_ctx, TSHighwayExitCtx* vertex_exctx, PetscReal rho, PetscReal v, PetscReal t, PetscReal deltat, PetscReal* density_fact);
  
struct _ts_Network {
  DM              network;
  MPI_Comm        comm;
  PetscInt        l_nedges, l_nvertices; /* local num of edges and vertices */
  PetscInt        g_nedges, g_nvertices; /* global num of edges and vertices */

  PetscInt        highway_key, vertex_key;
  PetscBool       manual_jacobian;/*=PETSC_FALSE; does the network have a manually-provided Jacobian? */

  PetscInt*       edgelist;/*local edge list */
  PetscInt        g_discrete_dimension, l_discrete_dimension; /* number of global and local
							   nodes used in DMDA discretization */
  Vec             g_X, l_X, l_dXdt; /* vectors used for local function evaluation at nodes (X is (rho, v)) */

  TSHighway       highways;
  TSHighwayVertex vertices;

  TSArrivalDistributionType arrival_dist_t;/*=TS_POISSON_DYNAMIC;*/
  TSSpeedDensityModel       speed_model;/*=TS_LINEAR;*/
  TSProblemType             problem_type;/*=TS_TIME_STEP;*/
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_Network *TSNetwork;

extern PetscErrorCode TSNetworkCreate(MPI_Comm,TSNetwork*);

extern PetscErrorCode TSNetworkCreateWithStructure(TSNetwork* network, DM* netdm, PetscInt network_case, const char* filename);

extern PetscErrorCode TSNetworkDistribute(MPI_Comm, TSNetwork, PetscBool);

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

extern PetscErrorCode TSHighwayNetIFunction(TS, PetscReal, Vec, Vec, Vec, void*);

extern PetscErrorCode TSNetworkSolve(TSNetwork, Vec);

extern PetscErrorCode TSNetworkSetFromOptions(TSNetwork);



#endif
