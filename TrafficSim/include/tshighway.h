#ifndef TS_HIGHWAY_H
#define TS_HIGHWAY_H

#include <petsc.h>
#include <mpi.h>
#include "tsdensity.h"

#define TS_COUNTY_NAME_LEN 3
#define TS_MAX_CITY_NAME_LEN 50
#define TS_MAX_ROAD_NAME_LEN 100


typedef enum {
	      TS_NS,
	      TS_EW,
	      TS_UNKNOWN_DIRECTION
} TSRoadDirection;


typedef enum {
	      TS_INTERSTATE,
	      TS_HIGHWAY,
	      TS_SURFACE_STREET
} TSRoadType;


struct _ts_HighwayEntryCtx{
  PetscReal                 postmile;
  TSArrivalDistributionType arrival_dist;
  PetscInt                  num_arrival_params;
  PetscReal*                arrival_params;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayEntryCtx TSHighwayEntryCtx;

extern PetscErrorCode HighwayEntryCreate(TSHighwayEntryCtx** entry, const PetscReal postmile,
				         const TSArrivalDistributionType arrival_dist_t,
				         PetscInt* n_params, PetscReal* params);

extern PetscErrorCode StaticPoissonHighwayEntryCreate(TSHighwayEntryCtx**, const PetscReal, const PetscReal*);

extern PetscErrorCode HighwayEntryDestroy(TSHighwayEntryCtx*);

typedef enum {
	      TS_NO_EXIT,
	      TS_STATIC_EXIT,/* time-independent exit probability for each passing car */
	      TS_DYNAMIC_EXIT /* time-dependent exit probability for each passing car */
} TSExitType;

/* struct to hold params for multi-parameter distribution; TODO: maybe put this in the .c file? */
struct _ts_multiparams{
  PetscInt   n;
  PetscReal* params=NULL;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_multiparams TSExitParams;
  
struct _ts_HighwayExitCtx{
  PetscReal    postmile;
  TSExitType   type;
  TSExitParams prob_params;
}PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayExitCtx TSHighwayExitCtx;

extern PetscErrorCode HighwayExitCreate(TSHighwayExitCtx** exit, PetscReal postmile,
				        TSExitType exit_t, TSExitParams* params);

extern PetscErrorCode HighwayExitDestroy(TSHighwayExitCtx*);

extern PetscErrorCode HighwayExitGetParams(TSHighwayExitCtx*, PetscInt*, PetscReal*);

extern PetscErrorCode HighwayExitSetParams(TSHighwayExitCtx* exit, PetscInt nparam, PetscReal* params);

struct _ts_InterchangeCtx{

  PetscInt  exit_road_id=-1; /* id for the road you're leaving */
  PetscInt  entry_road_id=-1; /* id for the road you're entering */
  PetscReal postmile=-1.0; /* posted mile marker location*/
  PetscReal continue_p = 1.0; /* probability that a car will go through this interchange,
				 only relevant if this interchange is part of a split vertex*/
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_InterchangeCtx TSInterchangeCtx;


typedef struct {
  PetscReal   rho; /* traffic density */
  PetscReal   v; /* traffic speed */
  PetscReal   q; /* traffix flux, q = rho * v */
} TSHighwayTrafficField;

typedef struct {
  PetscReal rho_ahead, v_ahead; /* boundary values ahead */
  PetscReal rho_behind, v_behind; /* boundary values behind */
} TSHighwayTrafficBoundary;

struct _ts_HighwayLocalSolverCtx {
  TSHighwayTrafficBoundary bc;
  Vec                      x;
  TSHighwayTrafficField    *old_rho_v_q;
  PetscReal                dt;
  DM                       da;
  PetscInt                 discrete_dimension; /* number of nodes used in DMDA discretization */
  Mat                      jac; /* Jacobian */
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayLocalSolverCtx TSHighwayLocalSolverCtx;


struct _ts_HighwayIdenfication {
  
  PetscInt           district=-1; /* district ID (e.g. in CA, comes from Caltrans)*/
  TSRoadDirection    direction;
  PetscReal          postmile=-1.0;
  PetscInt           city_name_length = 0;
  char               city_name[TS_MAX_CITY_NAME_LEN];
  char               county[TS_COUNTY_NAME_LEN];
  /*if two letters, add a trailing X (e.g. LA -> LAX)*/
  PetscInt           name_length = 0;
  char               name[TS_MAX_ROAD_NAME_LEN];
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayIdentification TSHighwayIdentification;


struct _ts_HighwayTrafficData {
  PetscReal          back_peak_hour=-1;
  PetscReal          back_peak_month=-1;
  PetscReal          back_aadt=-1; /* aadt = annual average daily traffic*/
  PetscReal          ahead_peak_hour=-1;
  PetscReal          ahead_peak_month=-1;
  PetscReal          ahead_aadt=-1;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayTrafficData TSHighwayTrafficData;


struct _ts_HighwayCtx{

  PetscInt                 id;/*internal identification for easy indexing*/
  
  TSHighwayIdentification* id_info;

  /* size and location data */
  PetscReal                length=-1; /* total length of this highway */
  PetscInt                 num_entries=0;
  TSHighwayEntryCtx*       entries=NULL;
  PetscInt                 num_exits=0;
  TSHighwayExitCtx*        exits=NULL;
  PetscInt                 num_interchanges=1;/* num interchanges per pointer below (if allocated)*/
  TSInterchangeCtx*        interchanges=NULL;
  TSInterchangeCtx*        interchanges2=NULL;
  TSInterchangeCtx*        interchanges3=NULL;
  PetscReal                speed_limit=-1;
  PetscInt                 num_lanes=-1;
  /* recorded traffic data */
  TSHighwayTrafficData*    traffic_data;

  /* simulation data */
  TSHighwayLocalSolverCtx* solver_ctx;

} PETSC_ATTRIBUTEALIGNED(sizeof(PetscScalar));

typedef struct _ts_HighwayCtx* TSHighway;

extern PetscErrorCode HighwayCreate(MPI_Comm, TSHighway*);

extern PetscErrorCode HighwayDestroy(MPI_Comm, TSHighway*);

/*extern PetscErrorCode HighwayCreateJacobian(TSHighway, Mat*, Mat*[]);

  extern PetscErrorCode HighwayDestroyJacobian(TSHighway);*/

extern PetscErrorCode HighwayComputeSpeedFromDensity(TSHighway, TSSpeedDensityModel);

extern PetscErrorCode HighwayGetCurrentDensity(TSHighway, PetscInt, PetscReal*, TSSpeedDensityModel);

extern PetscErrorCode HighwayGetCurrentSpeed(TSHighway, PetscInt, PetscReal*, TSSpeedDensityModel);

/* flux is speed times density, and can be integrated and compared to counts */
extern PetscErrorCode HighwayGetCurrentFlux(TSHighway, PetscInt, PetscReal*, TSSpeedDensityModel);

extern PetscErrorCode HighwaySetIdentification(TSHighway, const PetscInt*, const TSRoadDirection*,
					       const PetscReal*, const PetscReal*, const PetscInt*,
					       char*, char*, const PetscInt*, char*);

extern PetscErrorCode HighwaySetTrafficData(TSHighway, const PetscInt*, TSHighwayEntryCtx*,
					    const PetscInt*, TSHighwayExitCtx*, const PetscInt*,
					    TSInterchangeCtx*, const PetscReal*, const PetscReal*,
					    const PetscReal*, const PetscReal*, const PetscReal*,
					    const PetscReal*, const PetscInt*);


extern PetscErrorCode HighwaySetInfo(TSHighway highway, const PetscInt* district,
				    const TSRoadDirection* direction, const PetscReal* postmile,
				    const PetscInt* city_name_len, const char* city_name,
				    const char* county_name, const PetscInt* highway_name_len,
				    const char* highway_name, const PetscInt* num_entries,
				    TSHighwayEntryCtx* entries, const PetscInt* num_exits,
				    TSHighwayExitCtx* exits, const PetscInt* num_interchanges,
				    TSInterchangeCtx* interchanges, const PetscReal* back_peak_hr,
				    const PetscReal* back_peak_month, const PetscReal* ahead_peak_hr,
				     const PetscReal* ahead_peak_month, const PetscReal* ahead_aadt,
				     const PetscReal* speed_limit, const PetscInt* num_lanes);

extern PetscErrorCode HighwaySetUp(TSHighway);

extern PetscErrorCode HighwaySetCounty(TSHighway, const char*);

extern PetscErrorCode HighwaySetCity(TSHighway, const char*, PetscInt);

extern PetscErrorCode HighwaySetName(TSHighway const char*, PetscInt);

extern PetscErrorCode HighwaySetBackPeakHour(TSHighway, PetscReal);

extern PetscErrorCode HighwaySetAheadPeakHour(TSHighway, PetscReal);

extern PetscErrorCode HighwaySetPeakHour(TSHighway ctx, PetscReal ahead_peak, PetscReal behind_peak);

extern PetscErrorCode HighwaySetBackPeakMonth(TSHighway, PetscReal);

extern PetscErrorCode HighwaySetAheadPeakMonth(TSHighway, PetscReal);

extern PetscErrorCode HighwaySetPeakMonth(TSHighway ctx, PetscReal ahead_peak, PetscReal behind_peak);

extern PetscErrorCode HighwaySetBackAADT(TSHighway, PetscReal);

extern PetscErrorCode HighwaySetAheadAADT(TSHighway, PetscReal);

extern PetscErrorCode HighwaySetAADT(TSHighway ctx, PetscReal ahead_aadt, PetscReal back_aadt);



#endif
