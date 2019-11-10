#ifndef TNET_H
#define TNET_H

#include <petscsnes.h>
#include <petscdmnetwork.h>
#include <mpi.h>

#define COUNTY_NAME_LEN 3
#define MAX_CITY_NAME_LEN 50
#define MAX_ROAD_NAME_LEN 100


typedef enum {
	      TS_NS,
	      TS_EW
} TSRoadDirection;


typedef enum {
	      TS_INTERSTATE,
	      TS_HIGHWAY,
	      TS_SURFACE_STREET
} TSRoadType;


typedef enum {
	      TS_POISSON_STATIC,/* Poisson dist w/ time-independent rate */
	      TS_POISSON_DYNAMIC, /* Poisson dist w/ time-dependent rate (e.g. a Hidden Markov Model where the transition densities given the hidden state are Poisson) */
	      TS_POISSON_STATIC_MIXTURE, /* static mixture of independent Poisson distributions */
	      TS_POISSON_DYNAMIC_MIXTURE /* HMM with conditional densities being a Poisson mixture */
} TSArrivalDistributionType;

#define TS_SIZEOF_BIGGEST(a, b, c) \
  #if sizeof(a) > sizeof(b) \
    #if sizeof(a) > sizeof(c) \
      sizeof(a) \
    #else \
      sizeof(c) \
    #endif \
  #else \
    #if sizeof(b) > sizeof(c) \			\
      sizeof(b) \
    #else \
      sizeof(c) \
    #endif \
  #endif 


struct _ts_HighwayEntryCtx{
  PetscReal                 postmile;
  TSArrivalDistributionType arrival_dist;
  PetscInt                  num_arrival_params;
  PetscReal*                arrival_params;
} PETSC_ATTRIBUTEALIGNED(TS_SIZEOF_BIGGEST(TSArrivalDistributionType, PetscReal, PetscReal*));

typedef struct _ts_HighwayEntryCtx TSHighwayEntryCtx;

extern PetscErrorCode CreateHighwayEntry(TSHighwayEntryCtx* entry, const PetscReal postmile,
				         const TSArrivalDistributionType arrival_dist_t,
				         const PetscInt* n_params, const PetscReal* params);

extern PetscErrorCode CreateStaticPoissonHighwayEntry(TSHighwayEntryCtx*, const PetscReal, const PetscReal*);


typedef enum {
	      TS_STATIC_EXIT,/* time-independent exit probability for each passing car */
	      TS_DYNAMIC_EXIT /* time-dependent exit probability for each passing car */
} TSExitType;

/* struct to hold params for multi-parameter distribution; TODO: maybe put this in the .c file? */
struct _ts_multiparams{
  PetscInt   n;
  PetscReal* params;
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscReal*));

typedef union {
  PetscReal p;
  struct _ts_multiparams params;
} TSExitParams;
  
struct _ts_HighwayExitCtx{
  PetscReal    postmile;
  TSExitType   type;
  TSExitParams prob_params;
}PETSC_ATTRIBUTEALIGNED(TS_SIZEOF_BIGGEST(TSExitParams, TSExitType, PetscReal));

typedef struct _ts_HighwayExitCtx TSHighwayExitCtx;

extern PetscErrorCode CreateHighwayExit(TSHighwayExitCtx* exit, PetscReal postmile,
				        TSExitType exit_t, TSExitParams* params);

extern PetscErrorCode CreateHighwayExit(TSHighwayExitCtx* exit, PetscReal postmile,
				        TSExitType exit_t, PetscInt n_params,
				        PetscReal* params);



struct _ts_InterchangeCtx{

  PetscInt  exit_road_id; /* id for the road you're leaving */
  PetscInt  entry_road_id; /* id for the road you're entering */
  PetscReal postmile; /* posted mile marker location*/
} PETSC_ATTRIBUTEALIGNED(sizeof(PetscReal));

typedef struct _ts_InterchangeCtx TSInterchangeCtx;


struct _ts_HighwayCtx{

  PetscInt           id;/*internal identification for easy indexing*/
  PetscInt           district; /* district ID from caltrans*/
  TSRoadDirection    direction;
  PetscReal          postmile;
  PetscInt           num_entries;
  TSHighwayEntryCtx* entries;
  PetscInt           num_exits;
  TSHighwayExitCtx*  exits;
  PetscInt           num_interchanges;
  TSInterchangeCtx*  interchanges;
  PetscInt           city_name_length = 0;
  char               city_name[MAX_CITY_NAME_LEN];
  char               county[COUNTY_NAME_LEN];
  /*if two letters, add a trailing X (e.g. LA -> LAX)*/
  PetscInt           name_length = 0;
  char               name[MAX_ROAD_NAME_LEN];
  PetscReal          back_peak_hour;
  PetscReal          back_peak_month;
  PetscReal          back_aadt; /* aadt = annual average daily traffic*/
  PetscReal          ahead_peak_hour;
  PetscReal          ahead_peak_month;
  PetscReal          ahead_aadt; 
} PETSC_ATTRIBUTEALIGNED(TS_SIZEOF_BIGGEST(TSHighwayEntryCtx*, TSRoadDirection, PetscReal));

typedef struct _ts_HighwayCtx* TSHighway;

extern PetscErrorCode TrafficNetworkCreate(MPI_Comm comm, DM* network_dm);


extern PetscErrorCode CreateHighway(TSHighway* highway, const PetscInt* district,
				    const TSRoadDirection* direction, const PetscReal* postmile,
				    const PetscInt* city_name_len, const char* city_name,
				    const char* county_name, const PetscInt* highway_name_len,
				    const char* highway_name, const PetscInt* num_entries,
				    const TSHighwayEntryCtx* entries, const PetscInt* num_exits,
				    const TSHighwayExitCtx* exits, const PetscInt* num_interchanges,
				    const TSInterchangeCtx* interchanges, const PetscReal* back_peak_hr,
				    const PetscReal* back_peak_month, const PetscReal* ahead_peak_hr,
				    const PetscReal* ahead_peak_month, const PetscReal* ahead_aadt);				    

extern PetscErrorCode SetHighwayCounty(TSHighway, const char*);

extern PetscErrorCode SetHighwayCity(TSHighway, const char*, PetscInt);

extern PetscErrorCode SetHighwayName(TSHighway const char*, PetscInt);

extern PetscErrorCode SetHighwayBackPeakHour(TSHighway, PetscReal);

extern PetscErrorCode SetHighwayAheadPeakHour(TSHighway, PetscReal);

extern PetscErrorCode SetHighwayPeakHour(TSHighway ctx, PetscReal ahead_peak, PetscReal behind_peak);

extern PetscErrorCode SetHighwayBackPeakMonth(TSHighway, PetscReal);

extern PetscErrorCode SetHighwayAheadPeakMonth(TSHighway, PetscReal);

extern PetscErrorCode SetHighwayPeakMonth(TSHighway ctx, PetscReal ahead_peak, PetscReal behind_peak);

extern PetscErrorCode SetHighwayBackAADT(TSHighway, PetscReal);

extern PetscErrorCode SetHighwayAheadAADT(TSHighway, PetscReal);

extern PetscErrorCode SetHighwayAADT(TSHighway ctx, PetscReal ahead_aadt, PetscReal back_aadt);


#endif
