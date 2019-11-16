#include "../include/highway.h"

/* subroutines for TSHighway handline */

/*
  HighwayCreate - Create TSHighway object

  Input parameters:
  comm - MPI communicator

  Output Parameters:
  highway - pointer to memory location we want to allocate for the highway
 */

PetscErrorCode HighwayCreate(MPI_Comm comm, TSHighway* highway)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscNew(highway);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* 
   HighwaySetIdentification - Set identifying information for a TSHighway
   object. 

   Input parameters:
   district - district identification
   direction - which direction does the highway run? (TS_NS, TS_EW)
   length - how long is the road (in desired distance units)
   postmile - posted mile at data collection point
   highway_name_len - length of highway name in characters
   highway_name - name of highway
   county_name - name of county -- MUST BE THREE UPPER-CASE CHARACTERS
   city_name_len - length of name of city in characters
   city_name - name of city data is collected in
   
   Output parameters:
   highway - highway whose info we want to set

   Notes:
   If you don't want to set a particular parameter, pass NULL. If you pass NULL as a name length and do not pass NULL for the corresponding name parameter, 
   TSim will try to guess the length of the name, but this is UNSAFE; DO NOT DO IT IF YOU HAVE ANY CHOICE
 */

PetscErrorCode HighwaySetIdentification(TSHighway highway, const PetscInt* district,
					const TSRoadDirection* direction, const PetscReal* length,
					const PetscReal* postmile, const PetscInt* highway_name_len,
					char* highway_name, char* county_name,
					const PetscInt* city_name_len, char* city_name)
{
  PetscFunctionBegin;
  if PetscLikely(district){
    highway->district = *district;
  } else {
    highway->district = -1;
  }
  
  if PetscLikely(direction){
    highway->direction = *direction;
  }
  else {
    highway->direction = TS_UNKNOWN_DIRECTION;
  }

  if(length){
    if(*length > 0){
      highway->length = *length;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Highway length must be positive.\n");
    }
  }
  
  if PetscLikely(postmile){
    highway->postmile = *postmile;
  } else {
    highway->postmile = -1.0;
  }

  if(highway_name){
    if PetscLikely(highway_name_len){
      if PetscUnlikely(*highway_name_len > TS_MAX_ROAD_NAME_LEN){
	  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "highway_name_len (including null-terminator) must be less than or equal to TS_MAX_ROAD_NAME_LEN.\n");
      } else {
	highway->name_length = *highway_name_len;
	highway->name = highway_name;
      }
    } else {
      /* DANGER WARNING */
      PetcsInt namelen = (PetscInt)(sizeof(highway_name)/sizeof(highway_name[0]));
      if PetscUnlikely(namelen > TS_MAX_ROAD_NAME_LEN){
	  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "highway_name_len (including null-terminator) must be less than or equal to TS_MAX_ROAD_NAME_LEN.\n");
      } else {
	highway->name_length = namelen;
        highway->name = highway_name;
      }
    }
  }
  else {
    highway->name_length = -1;
  }

  if(county_name){
    if PetscUnlikely(sizeof(county_name)/sizeof(county_name[0]) != TS_COUNTY_NAME_LEN){
	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "county_name must have exactly TS_COUNTY_NAME_LEN characters.\n");
    } else {
      highway->county = county_name;
    }
  }

  if(city_name){
    if PetscLikely(city_name_len){
	if PetscUnlikely(*city_name_len > TS_MAX_CITY_NAME_LEN){
	    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ , "city_name_len must be less than or equal to TS_MAX_CITY_NAME_LEN.\n");
        } else {
	  highway->city_name_length = *city_name_len;
	  highway->city_name = city_name;
	}
      }
    else {
      PetscInt cnamelen = (PetscInt)(sizeof(city_name)/sizeof(city_name[0]));
      if(PetscUnlikely(cnamelen > TS_MAX_CITY_NAME_LEN)){
	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ , "city_name_len must be less than or equal to TS_MAX_CITY_NAME_LEN.\n");
      } else {
	highway->city_name_length = cnamelen;
	highway->city_name = city_name;
      }
    }
  }

  PetscFunctionReturn(0);
}

/*
  HighwaySetTrafficData - Set traffic count data for a TSHighway

  Input parameters:
  num_entries - how many places can cars enter the highway?
  entries - array of TSHighwayEntryCtx containing relevant info to the simulation
  num_exits, exits, num_interchanges, interchanges - same as above, but exits and interchanges instead of entries

  ahead/back input parameters: "ahead" means the count of cars passing by the count station BEFORE they get to the station (and generally intersection, if the station is at intersection), and "back" means the count AFTER the cars have passed the count station. For North-South highways, "ahead" is north of the location, and "back" is south of the location. For East-West highways, "ahead" is east of the station, and "back" is west of the location. (Notation based on Caltrans traffic data reporting notation).

  peak_month - estimate of the avg. daily traffic (ADT) in the largest single month (in terms of ADT) of the year
  peak_hr - highest one-hour count of the year, excluding some (Caltrans uses 30-50) "outlier" hours with significantly higher-than-normal counts.
  aadt - annual average daily traffic (AADT) is the average daily traffic count over one year
  speed_limit - speed limit on this highway
  num_lanes - number of lanes in each direction
 */
PetscErrorCode HighwaySetTrafficData(TSHighway highway, const PetscInt* num_entries,
				     TSHighwayEntryCtx* entries, const PetscInt* num_exits,
				     TSHighwayExitCtx* exits, const PetscInt* num_interchanges,
				     TSInterchangeCtx* interchanges, const PerscReal* back_peak_hr,
				     const PetscReal* back_peak_month, const PetscReal* back_aadt,
				     const PetscReal* ahead_peak_hr, const PetscReal* ahead_peak_month,
				     const PetscReal* ahead_aadt, const PetscReal* speed_limit,
				     const PetscInt* num_lanes)
{
  PetscFunctionBegin;
  if(entries){
    if(num_entries){
      if(*num_entries >= 0){
	highway->num_entries = *num_entries;
	highway->entries = entries;
      } else {
	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Number of entries must be positive or zero.\n");
      }
    }
  }/*entries*/

  if(exits){
    if(num_exits){
      if(*num_exits >= 0){
	highway->num_exits = *num_exits;
	highway->exits = exits;
      } else {
	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Number of exits must be positive or zero.\n");
      }
    } 
  }/*exits*/


  if(interchanges){
    if(num_interchanges){
      if(*num_interchanges >= 0){
	highway->num_interchanges = *num_interchanges;
	highway->interchanges = interchanges;
      } else {
	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Number of entries must be positive or zero.\n");
      }
    }
  }/*interchanges*/

  if(back_peak_hr){
    if(*back_peak_hr >=0){
      highway->back_peak_hour = *back_peak_hour;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*back_peak_hr must be positive or zero.\n");
    }
  }

  if(ahead_peak_month){
    if(*ahead_peak_month >=0){
      highway->ahead_peak_month = *ahead_peak_month;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*ahead_peak_month must be positive or zero.\n");
    }
  }

    if(back_peak_month){
    if(*back_peak_month >=0){
      highway->back_peak_month = *back_peak_month;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*back_peak_month must be positive or zero.\n");
    }
  }

  if(ahead_peak_hr){
    if(*ahead_peak_hr >=0){
      highway->ahead_peak_hour = *ahead_peak_hour;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*ahead_peak_hr must be positive or zero.\n");
    }
  }


  if(back_aadt){
    if(*back_aadt >=0){
      highway->back_aadt = *back_aadt;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*back_aadt must be positive or zero.\n");
    }
  }

  if(ahead_aadt){
    if(*ahead_aadt >=0){
      highway->ahead_aadt = *ahead_aadt;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*ahead_aadt must be positive or zero.\n");
    }
  }

  if(speed_limit){
    if(*speed_limit >=0){
      highway->speed_limit = *speed_limit;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*speed_limit must be positive or zero.\n");
    }
  }

  if(num_lanes){
    if(*num_lanes >=0){
      highway->num_lanes = *num_lanes;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "*num_lanes must be positive or zero.\n");
    }
  }
  
  PetscFunctionReturn(0);
}

#if 0
/*
  HighwayCreateJacobian - create (preallocate) Jacobian matrices for a highway section (at the start, in the middle, and at the end of the segment).

  Input Parameters:
  highway - a TSHighway object
  J_pre - array of precomputed Jacobian matrices (set to NULL if not available) to reuse. This is only to avoid unnecessary calls to malloc, the values need not be correct, or even set at all if the nonzero structure is correctly set.

  Output Parameter:
  J_out - array of three preallocated (but not computed) Jacobian matrices 

 */
/*
	  
PetscErrorCode HighwayCreateJacobian(TSHighway highway, Mat* J_pre, Mat* J_out[])
{
  PetscErrorCode ierr;
  Mat*           J_highway;
  PetscInt       M, rows[2], cols[2], *nz;
  PetscScalar    *aa;


  PetscFunctionBegin;
  if(J_pre){
    *J_out = J_pre;
    highway->jac = J_pre;
    ierr = PetscObjectReference((PetscObject)(J_pre[0])); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  ierr = PetscMalloc1(3, &J_highway);CHKERRQ(ierr);
*/
  /* this highway segment's Jacobian */
				     /*
  ierr = DMSetMatrixStructureOnly(highway->da, PETSC_TRUE);CHKERRQ(ierr);
  ierr = DMCreateMatrix(highway->da, &J_highway[0]);CHKERRQ(ierr);
  ierr = DMSetMatrixStructureOnly(highway->da, PETSC_FALSE);CHKERRQ(ierr);
				     */
  /* Jacobian for ahead highway-highway junction (ahead of monitoring station).

     The variables to be solved for are rho and v, so the Jacobian structure is
     [ [1  dv/drho]
       [drho/dv r ] ]
   */
				     /*
  ierr = MatGetSize(J_highway[0], &M, NULL);CHKERRQ(ierr);
  ierr = PetscCalloc2(M, &nz, 4, &aa);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_SELF, &J_highway[1]);CHKERRQ(ierr);
  ierr = MatSetSizes(J_highway[1], PETSC_DECIDE, PETSC_DECIDE, M, 2);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J_highway[1]);CHKERRQ(ierr);
  ierr = MatSetOption(J_highway[1], MAT_STRUCTURE_ONLY, PETSC_TRUE);CHKERRQ(ierr);
				     */
  /* nonzero structure by row */
				     /*
  nz[0] = 2;
  nz[1] = 2;


  PetscFunctionReturn(0);
}
				     */

#endif

  

  
					
					
