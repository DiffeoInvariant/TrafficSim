#ifndef TNET_H
#define TNET_H

#include <petscsnes.h>
#include <petscdmnetwork.h>

#define COUNTY_NAME_LEN 3
#define MAX_CITY_NAME_LEN 50
#define MAX_ROAD_NAME_LEN 100


typedef enum {
	      NS,
	      EW
} RoadDirection;

struct _p_RoadCtx{

  PetscInt      id;/*internal identification for easy indexing*/
  PetscInt      district; /* district ID from caltrans*/
  RoadDirection direction;
  PetscReal     postmile;
  PetscInt      city_name_length = 0;
  char          city_name[MAX_CITY_NAME_LEN];
  char          county[COUNTY_NAME_LEN];
  /*if two letters, add a trailing X (e.g. LA -> LAX)*/
  PetscInt      name_length = 0;
  char          name[MAX_ROAD_NAME_LEN];
  PetscReal     back_peak_hour;
  PetscReal     back_peak_month;
  PetscReal     back_aadt; /* aadt = annual average daily traffic*/
  PetscReal     ahead_peak_hour;
  PetscReal     ahead_peak_month;
  PetscReal     ahead_aadt; 
};

typedef struct _p_RoadCtx RoadCtx;


struct _p_InterchangeCtx{

  PetscInt  exit_road_id; /* id for the road you're leaving */
  PetscInt  entry_road_id; /* id for the road you're entering */
  PetscReal postmile; /* posted mile marker location*/
};

typedef struct _p_InterchangeCtx InterchangeCtx;



extern PetscErrorCode SetRoadCounty(RoadCtx*, const char*);

extern PetscErrorCode SetRoadCity(RoadCtx*, const char*, PetscInt);

extern PetscErrorCode SetRoadName(RoadCtx*, const char*, PetscInt);

extern PetscErrorCode SetRoadBackPeakHour(RoadCtx*, PetscReal);

extern PetscErrorCode SetRoadAheadPeakHour(RoadCtx*, PetscReal);

extern PetscErrorCode SetRoadPeakHour(RoadCtx* ctx, PetscReal ahead_peak, PetscReal behind_peak);

extern PetscErrorCode SetRoadBackPeakMonth(RoadCtx*, PetscReal);

extern PetscErrorCode SetRoadAheadPeakMonth(RoadCtx*, PetscReal);

extern PetscErrorCode SetRoadPeakMonth(RoadCtx* ctx, PetscReal ahead_peak, PetscReal behind_peak);

extern PetscErrorCode SetRoadBackAADT(RoadCtx*, PetscReal);

extern PetscErrorCode SetRoadAheadAADT(RoadCtx*, PetscReal);

extern PetscErrorCode SetRoadAADT(RoadCtx* ctx, PetscReal ahead_aadt, PetscReal back_aadt);






#endif
