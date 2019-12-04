#include "../include/tsnetwork.h"

PetscErrorCode TSNetworkCreate(MPI_Comm comm, TSNetwork* network){
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscNew(network);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode TSNetworkDestroy(TSNetwork network){
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscFree(network->edgelist);CHKERRQ(ierr);
  ierr = PetscFree(network->highways);CHKERRQ(ierr);
  ierr = PetscFree(network->vertices);CHKERRQ(ierr);
  ierr = VecDestroy(network->l_X);CHKERRQ(ierr);
  ierr = VecDestroy(network->l_dXdt);CHKERRQ(ierr);
  ierr = DMDestroy(network->network);CHKERRQ(ierr);
  ierr = PetscFree(network);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TSNetworkCreateWithStructure(MPI_Comm comm, TSNetwork* network, PetscInt network_case, const char* filename)
{
  PetscErrorCode     ierr;
  PetscInt           nhighway, nentry, nexit, ninterchange,nvertices, nlane;
  PetscInt           district, glob_id, i, discrete_dimension, *edgelist;
  PetscMPIInt        size, rank;
  TSNetwork          net=NULL;
  TSHighway          highways=NULL;
  TSHighwayVertex    vertices=NULL;
  TSInterchangeCtx*  interchanges=NULL;
  TSHighwayEntryCtx* entries=NULL;
  TSHighwayExitCtx*  exits=NULL;
  TSRoadDirection    direction;
  PetscReal          speed_lim, *arrival_params=NULL, *exit_params=NULL;
  TSExitParams       epars1, epars2, epars3, epars4, epars5;

  PetscFunctionBegin;

  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);
  ierr = PetscCalloc1(1, &net);CHKERRQ(ierr);
  net->comm = comm;
  *network = net;

  /* create a TSNetwork on process 0 */
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Setting up TSNetwork in pre-set case %D.\n", network_case);CHKERRQ(ierr);

  switch(network_case) {

  case 0:
    /* network_case = 0: simple ring of 7 end-to-end connected highways with an 8th (slightly narrower) highway (H4)
       cutting across the middle, entries (later exits) at V0 and V3, exits (later entries) at V6, V4, V5: 

       
       V0 --H0--> V1 --H1--> V2 --H2--> V3
       ^                     |          |
       |                     |          |
       H7                    H4         H3
       |                     |          |
       |                     v          v
       V6       <---H6----   V5 <--H5-- V4

    */
    net-> g_nedges = 8;
    net-> l_nedges = 8;
    net-> g_nvertices = 7;
    net-> l_nvertices = 7;
    nentry = 7;
    nexit = 7;

    if(!rank){
      nedges = net->g_ndeges;
      nvertices = net->g_nvertices;
      ierr = PetscCalloc1(2 * g_nedges, &edgelist);CHKERRQ(ierr);
      edgelist[0] = 0; edgelist[1] = 1; /* H0 */
      edgelist[2] = 1; edgelist[3] = 2; /*H1 */
      edgelist[4] = 2; edgelist[5] = 3; /* H2 */
      edgelist[6] = 3; edgelist[7] = 4; /* H3 */
      edgelist[8] = 2; edgelist[9] = 6; /*H4*/
      edgelist[10] = 4; edgelist[11] = 5; /*H5*/
      edgelist[12] = 5; edgelist[13] = 6; /*H6*/
      edgelist[14] = 6; edgelist[15] = 0; /* H7 */

      /* allocate highways and vertices, incl entries, exits, and interchanges */
      ninterchange = nvertex;
      ierr = PetscCalloc7(nvertices, &vertices,
			  ninterchange, &interchanges,
			  nentry, &entries,
			  nentry * 8, &arrival_params,
			  nexit * 8, &exit_params,
			  nexit, &exits,
			  nedges, &highways);CHKERRQ(ierr);
      
      for(i = 0; i < nvertices; ++i){
	vertices[i].id=i;
	switch(i){
	case 0:
	  /* V0, one entry, one exit, one interchange */
	  interchanges[i].exit_road_id = 7;
	  interchanges[i].entry_road_id = 0;
	  interchanges[i].postmile = 0.0;
	  vertices[i].intc_ctx = interchanges + i;/* pointer arithmetic */

	  entries[i].postmile = 0.0;
	  entries[i].arrival_dist = TS_POISSON_DYNAMIC;
	  /* the N arrival parameters for TS_POISSON_DYNAMIC are given as
	     an array of length 2N with entries
	     [t0,lambda0, t1, lambda1, ..., tN-1, lambdaN-1],
	     where the ti are the times at which lambdai becomes the arrival
	     distribution's parameter. */
	  entries[i].num_arrival_params = 4;
	  arrival_params[0] = 0.0; arrival_params[1] = 100.0;
	  arrival_params[2] = 2.5; arrival_params[3] = 15.0; /* 2.5 hour morning rush hour */
	  arrival_params[4] = 5.0; arrival_params[5] = 30.0;
	  arrival_params[6] = 9.0; arrival_params[7] = 8.0;
	  entries[i].arrival_params = arrival_params;/*giving entries the whole pointer
						       containing all the arrival params,
						       which is stupid and unsafe (since I now
						       have to ensure that nobody changes anybody
						       else's arrival parameters, but eh, whatever for now).*/
	  vertices[i].entr_ctx = entries[i];
	  vertices[i].arrival_dist = TS_POISSON_DYNAMIC;

	  exits[i].postmile = 0.0;
	  exits[i].type = TS_DYNAMIC_EXIT;
	  /* same deal as with arrival params, but the odd-numbered-indexed values are probabilities
	     between 0 and 1 (inclusive) that any given car will leave the highway on this exit 
	     in the given timeframe.*/
	  vertices[i].exit_dist = TS_DYNAMIC_EXIT;
	  exit_params[0] = 0.0; exit_params[1] = 0.001;
	  exit_params[2] = 2.5; exit_params[3] = 0.005;
	  exit_params[4] = 5.0; exit_params[5] = 0.21;
	  exit_params[6] = 9.0; exit_params[7] = 0.02;
	  epars1.n = 4;
	  epars1.params = exit_params;
	  exits[i].prob_params = epars1;
	  
	  vertices[i].exit_ctx = exits[i];
	  break;
	case 1:
	  /* V1, one interchange, one exit, one entry */
	  interchanges[i].exit_road_id = 0;
	  interchanges[i].entry_road_id = 1;
	  interchanges[i].postmile = 10.0;
	  vertices[i].intc_ctx = interchanges + i;/* pointer arithmetic */

	  entries[i].postmile = 0.0;
	  entries[i].arrival_dist = TS_POISSON_DYNAMIC;
	  /* the N arrival parameters for TS_POISSON_DYNAMIC are given as
	     an array of length 2N with entries
	     [t0,lambda0, t1, lambda1, ..., tN-1, lambdaN-1],
	     where the ti are the times at which lambdai becomes the arrival
	     distribution's parameter. */
	  entries[i].num_arrival_params = 4;
	  arrival_params[8 * i + 0] = 0.0; arrival_params[8 * i + 1] = 140.0;
	  arrival_params[8 * i + 2] = 2.5; arrival_params[8 * i + 3] = 51.0; /* 2.5 hour morning rush hour */
	  arrival_params[8 * i + 4] = 5.0; arrival_params[8 * i + 5] = 20.0;
	  arrival_params[8 * i + 6] = 9.0; arrival_params[8 * i + 7] = 43.0;
	  entries[i].arrival_params = arrival_params + 8 * i;/* hey, we're adding some pointer arithmetic to this unholy assignment protocol we're using for the arrival/exit parameters. Why? Because the universe hates you, that's why.*/
	  vertices[i].entr_ctx = entries[i];
	  vertices[i].arrival_dist = TS_POISSON_DYNAMIC;

	  exits[i].postmile = 10.0;
	  exits[i].type = TS_DYNAMIC_EXIT;
	  /* same deal as with arrival params, but the odd-numbered-indexed values are probabilities
	     between 0 and 1 (inclusive) that any given car will leave the highway on this exit 
	     in the given timeframe.*/
	  vertices[i].exit_dist = TS_DYNAMIC_EXIT;
	  exit_params[8 * i + 0] = 0.0; exit_params[8 * i + 1] = 0.08;
	  exit_params[8 * i + 2] = 2.5; exit_params[8 * i + 3] = 0.005;
	  exit_params[8 * i + 4] = 5.0; exit_params[8 * i + 5] = 0.24;
	  exit_params[8 * i + 6] = 9.0; exit_params[8 * i + 7] = 0.04;
	  epars2.n = 4;
	  epars2.params = exit_params + 8 * i;
	  exits[i].prob_params = epars2;
	  
	  vertices[i].exit_ctx = exits[i];
	  
	  
	  
	  

	}
      }
	
    
    




  } /* end switch */
  
  
