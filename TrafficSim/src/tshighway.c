#include "../include/tshighway.h"
#include <cstdlib.h>

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

PetscErrorCode HighwayDestroy(TSHighway highway){
  PetscErrorCode ierr;
  PetscInt       id;
  TSExitParams   *epars;
  PetscFunctionBegin;

  ierr = PetscFree(highway->id_info);CHKERRQ(ierr);

  if(highway->entries){
    for(id = 0; id < highway->num_entries; ++id){
      ierr = PetscFree((highway->entries[id]).arrival_params);CHKERRQ(ierr);
    }
    ierr = PetscFree(highway->entries);CHKERRQ(ierr);
  }

  if(highway->exits){
    for(id = 0; id < highway->num_exits; ++id){
      ierr = PetscFree(highway->exits[id].prob_params.params.params);CHKERRQ(ierr);
    }
    ierr = PetscFree(highway->exits);CHKERRQ(ierr);
  }
  
  ierr = PetscFree(highway->interchanges);CHKERRQ(ierr);
  ierr = PetscFree(highway->traffic_data);CHKERRQ(ierr);

  if(highway->solver_ctx){
    ierr = PetscFree(highway->solver_ctx->old_rho_v_q);CHKERRQ(ierr);
    ierr = MatDestroy(highway->solver_ctx->jac);CHKERRQ(ierr);
  }
  
  ierr = PetscFree(highway->solver_ctx);CHKERRQ(ierr);

  ierr = PetscFree(highway);CHKERRQ(ierr);

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
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* allocate memory for info */
  ierr = PetscNew(&(highway->id_info));
  
  
  if PetscLikely(district){
    highway->id_info->district = *district;
  } else {
    highway->id_info->district = -1;
  }
  
  if PetscLikely(direction){
    highway->id_info->direction = *direction;
  }
  else {
    highway->id_info->direction = TS_UNKNOWN_DIRECTION;
  }

  if(length){
    if(*length > 0){
      highway->id_info->length = *length;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Highway length must be positive.\n");
    }
  }
  
  if PetscLikely(postmile){
    highway->id_info->postmile = *postmile;
  } else {
    highway->id_info->postmile = -1.0;
  }

  if(highway_name){
    if PetscLikely(highway_name_len){
      if PetscUnlikely(*highway_name_len > TS_MAX_ROAD_NAME_LEN){
	  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "highway_name_len (including null-terminator) must be less than or equal to TS_MAX_ROAD_NAME_LEN.\n");
      } else {
	highway->id_info->name_length = *highway_name_len;
	highway->id_info->name = highway_name;
      }
    } else {
      /* DANGER WARNING */
      PetcsInt namelen = (PetscInt)(sizeof(highway_name)/sizeof(highway_name[0]));
      if PetscUnlikely(namelen > TS_MAX_ROAD_NAME_LEN){
	  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "highway_name_len (including null-terminator) must be less than or equal to TS_MAX_ROAD_NAME_LEN.\n");
      } else {
	highway->id_info->name_length = namelen;
        highway->id_info->name = highway_name;
      }
    }
  }
  else {
    highway->id_info->name_length = -1;
  }

  if(county_name){
    if PetscUnlikely(sizeof(county_name)/sizeof(county_name[0]) != TS_COUNTY_NAME_LEN){
	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "county_name must have exactly TS_COUNTY_NAME_LEN characters.\n");
    } else {
      highway->id_info->county = county_name;
    }
  }

  if(city_name){
    if PetscLikely(city_name_len){
	if PetscUnlikely(*city_name_len > TS_MAX_CITY_NAME_LEN){
	    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ , "city_name_len must be less than or equal to TS_MAX_CITY_NAME_LEN.\n");
        } else {
	  highway->id_info->city_name_length = *city_name_len;
	  highway->id_info->city_name = city_name;
	}
      }
    else {
      PetscInt cnamelen = (PetscInt)(sizeof(city_name)/sizeof(city_name[0]));
      if(PetscUnlikely(cnamelen > TS_MAX_CITY_NAME_LEN)){
	SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ , "city_name_len must be less than or equal to TS_MAX_CITY_NAME_LEN.\n");
      } else {
	highway->id_info->city_name_length = cnamelen;
	highway->if_info->city_name = city_name;
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
  PetscErrorCode ierr;
  size_t       nentry, nexit, ninter;
  PetscFunctionBegin;
  
  nentry = num_entries ? (*num_entries >= 0 ? *num_entries : 0) : 0;
  nexit = num_exits ? (*num_exits >=0 ? *num_exits : 0): 0;
  ninter = num_interchanges ? (*num_interchanges >= 0 ? *num_interchanges : 0) : 0;

  ierr = PetscCalloc3(nentry, &(highway->entries), nexit, &(highway->exits), ninter, &(highway->interchanges));

  /* maybe could delete a lot of the below code, but I think I should still
     return an error if an argument is out of range (so that the program doesn't fail 
     or misbehave silently if the end user made a typo). */
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


/*
  HighwayCreateJacobian - create (preallocate) Jacobian matrices for a highway section (at the start, in the middle, and at the end of the segment).

  Input Parameters:
  highway - a TSHighway object
  J_pre - array of precomputed Jacobian matrices (set to NULL if not available) to reuse. This is only to avoid unnecessary calls to malloc, the values need not be correct, or even set at all if the nonzero structure is correctly set.

  Output Parameter:
  J_out - array of three preallocated (but not computed) Jacobian matrices 

 */
	  
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
  /* this highway segment's Jacobian */
				     
  ierr = DMSetMatrixStructureOnly(highway->da, PETSC_TRUE);CHKERRQ(ierr);
  ierr = DMCreateMatrix(highway->da, &J_highway[0]);CHKERRQ(ierr);
  ierr = DMSetMatrixStructureOnly(highway->da, PETSC_FALSE);CHKERRQ(ierr);
				     
  /* Jacobian for ahead highway-highway junction (ahead of monitoring station).

     The variables to be solved for are rho and v, so the Jacobian structure is
     [ [1  dv/drho]
       [drho/dv 1 ] ]
   */
				     
  ierr = MatGetSize(J_highway[0], &M, NULL);CHKERRQ(ierr);
  ierr = PetscCalloc2(M, &nz, 4, &aa);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_SELF, &J_highway[1]);CHKERRQ(ierr);
  ierr = MatSetSizes(J_highway[1], PETSC_DECIDE, PETSC_DECIDE, M, 2);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J_highway[1]);CHKERRQ(ierr);
  ierr = MatSetOption(J_highway[1], MAT_STRUCTURE_ONLY, PETSC_TRUE);CHKERRQ(ierr);
				     
  /* nonzero structure by row for behind vertex*/
  				     
  nz[0] = 2;
  nz[1] = 2;
  rows[0] = 0; rows[1] = 1;
  cols[0] = 0; cols[1] = 1;
  ierr = MatAIJSetPreallocation(J_highway[1], 0, nz);CHKERRQ(ierr);
  ierr = MatSetValues(J_highway[1], 2, rows, 2, cols, aa, INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(J_highway[1], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J_highway[1], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* same for ahead vertex */
  ierr = MatCreate(PETSC_COMM_SELF, &J_highway[2]);CHKERRQ(ierr);
  ierr = MatSetSizes(J_highway[2], PETSC_DECIDE, PETSC_DECIDE, M, 2);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J_highway[2]);CHKERRQ(ierr);
  ierr = MatSetOption(J_highway[2], MAT_STRUCTURE_ONLY, PETSC_TRUE);CHKERRQ(ierr);
				      				     
  nz[0] = 0;
  nz[1] = 0;
  nz[M-2] = 2;
  nz[M-1] = 2;
 
  rows[0] = M-2; rows[1] = M-1;
  ierr = MatAIJSetPreallocation(J_highway[2], 0, nz);CHKERRQ(ierr);
  ierr = MatSetValues(J_highway[2], 2, rows, 2, cols, aa, INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(J_highway[2], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J_highway[2], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = PetscFree2(nz, aa);

  *J_out = J_highway;
  highway->jac = J_highway;


  PetscFunctionReturn(0);
}

PetscErrorCode HighwayDestroyJacobian(TSHighway highway)
{
  PetscErrorCode ierr;
  PetscInt       i;
  Mat            *J = highway->jac;

  PetscFunctionBegin;
  if(J){
    for(i=0; i < 3; ++i){
      ierr = MatDestroy(&J[i]);CHKERRQ(ierr);
    }
  }
  ierr = PetscFree(J);
  PetscFunctionReturn(0);
}
				     


PetscErrorCode HighwaySetUp(TSHighway highway)
{
  PetscErrorCode ierr;
  DMDALocalInfo  info;

  PetscFunctionBegin;
  ierr = DMDACreate1d(PETSC_COMM_SELF, DM_BOUNDARY_GHOSTED, highway->discrete_dimension, 2, 1, NULL, &(highway->da));CHKERRQ(ierr);
  ierr = DMSetFromOptions(highway->da);CHKERRQ(ierr);
  ierr = DMDASetFieldName(highway->da, 0, "rho");CHKERRQ(ierr);
  ierr = DMDASetFieldName(highway->da, 1, "v");CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(highway->da, 0, highway->length, 0, 0, 0, 0);CHKERRQ(ierr);
  ierr = DMDACreateGlobalVector(highway->da, &(highway->X));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode HighwayExitCreate(TSHighwayExitCtx** exit, PetscReal postmile, TSExitType exit_t, TSExitParams* params)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscNew(exit);CHKERRQ(ierr);
  (*exit)->postmile = postmile;
  (*exit)->type = exit_t;
  if(params){
    (*exit)->params = *params;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode HighwayExitDestroy(TSHighwayExitCtx* exit)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(exit);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode HighwayExitGetParams(TSHighwayExitCtx* exit, PetscInt* nparam, PetscReal* params)
{
  PetscFunctionBegin;
  *nparam = exit->prob_params.n;
  params = exit->prob_params.params;
  PetscFunctionReturn(0);
}

PetscErrorCode HighwayExitSetParams(TSHighwayExitCtx* exit, PetscInt nparam, PetscReal* params)
{
  PetscFunctionBegin;
  exit->prob_params.n = nparam;
  exit->prob_params.params = params;
  PetscFunctionReturn(0);
}
  
  
  
PetscErrorCode HighwayEntryCreate(TSHighwayEntryCtx** entry, const PetscReal postmile,
					  const TSArrivalDistributionType arrival_dist_t, PetscInt* n_params,
					  PetscReal* params)
 {
   PetscErrorCode ierr;
   PetscFunctionBegin;
   ierr = PetscNew(entry);CHKERRQ(ierr);
   (*entry)->postmile = postmile;
   (*entry)->arrival_dist = arrival_dist_t;
   if(n_params && params){
     (*entry)->num_arrival_params = *n_params;
     (*entry)->arrival_params = params;
   }

   PetscFunctionReturn(0);
 }

 PetscErrorCode HighwayEntryDestroy(TSHighwayEntryCtx* entry)
 {
   PetscErrorCode ierr;
   PetscFunctionBegin;
   ierr = PetscFree(entry);CHKERRQ(ierr);
   PetscFunctionReturn(0);
 }


PETSC_STATIC_INLINE PetscScalar drhodx(TSHighwayTrafficField *x, PetscInt i, PetscInt end, PetscReal dx)
{
  if(i == 0){
     return (x[i+1].rho - x[i-1].rho)/dx;
  } else if(i == end){
    return (x[i].rho - x[i-1].rho)/dx;
  } else {
    return (x[i+1].rho - x[i-1].rho)/(2*dx);
  }
}
     
    
/* highway local IFunction using an explicit first-order Lax-Friedrichs scheme.
   Note: this is a very bad scheme for this problem and its initial conditions 
   (for ex., it has strong dissipation and dispersion (Chu, C. K. (1978), Numerical Methods in Fluid Mechanics, Advances in Applied Mechanics, 18, Page 304), particularly for initial/boundary conditions with discontinuities (Thomas, J. W. (1995), Numerical Partial Differential Equations: Finite Difference Methods, Texts in Applied Mathematics, 22, Sec. 7.8)). 

   However, it can be extended and turned into a higher-order method relatively easily, and if the boundary conditions are sufficiently smooth, we know that the method will be stable iff
   |u(rho) * dt/dx| <= 1 (which is to say, iff |dt| <= dx/max(u(rho)) (usually = lim_{r->0}(dx/u(r)))) */
PetscErrorCode HighwayLocalIFunction_LaxFriedrichs(TSHighway highway, DMDALocalInfo *info,
				     PetscReal t, TSHighwayTrafficField *x,
				     TSHighwayTrafficField *xdot, PetscScalar *f)
{
  PetscErrorCode ierr;
  PetscInt       i, start, n, end;
  PetscReal      dx=highway->length/(info->mx - 1), dt=highway->dt;
  PetscScalar    u_avg, xold_avg, dudrho,drho, rho0, rho1, u0, u1;
  TSHighwayTrafficField *xold=highway->old_rho_v;

  PetscFunctionBegin;

  ierr = DMDAGetCorners(highway->da, &start, 0, 0, &n, 0, 0);CHKERRQ(ierr);
  end = start + n - 1;
  /* interior points*/
  for(i = start+1; i < end; ++i){
    old_u_avg = (xold[i+1].v + xold[i-1].v) * 0.5;
    rho0 = xold[i-1].rho;
    rho1 = xold[i+1].rho;
    u0 = xold[i-1].v;
    u1 = xold[i+1].v;

    xold_avg = 0.5*(rho0 + rho1);
    drho = (x[i].rho - xold_avg);
    /*drho + u(rho)*drho/dx*dt = 0 */
    f[2*(i-1) + 2] = drho + dt * old_u_avg * drhodx(xold, i, end, dx);
    
    ierr = TSSpeedDerivativeFromLinearDensity(&dudrho, highway->speed_limit, highway->rho_limit);CHKERRQ(ierr);
    /* du - du/drho * drho = 0*/
    f[2*(i-1)+3] = (x[i].v - old_u_avg) - dudrho * drho;
  }

  /* boundary points, apply characteristic equation here */
  
    
    
    
  
				     



   
#endif

  

  
					
					
