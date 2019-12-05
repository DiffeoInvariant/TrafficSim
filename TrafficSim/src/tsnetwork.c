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


PetscErrorCode TSNetworkDistribute(MPI_Comm comm, TSNetwork network)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank, size, tag=0;
  PetscInt       i, e, v, l_nedges, l_nvertices, g_nedges, g_nvertices, estart, eend, tag;
  PetscInt       *eowners, *edgelist = network->edgelist, *nvtx=NULL, *vtx_done=NULL;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);

  g_ndeges = network->l_nedges;
  g_nvertices = network->l_nvertices;
  /* everybody knows the number of local and global edges (and we store extra edges on rank 0) */
  ierr = MPI_Bcast(&g_nedges, 1, MPIU_INT, 0, comm);CHKERRQ(ierr);
  l_nedges = g_nedges/size;
  if(!rank){
    l_nedges += g_nedges - size * (g_nedges/size);
  }
  network->g_nedges = g_nedges;
  network->l_nedges = l_ndeges;

  ierr = PetscCalloc3(size+1, &eowners, size, &nvtx, g_nvertices, &vtx_done);CHKERRQ(ierr);
  /*gather the number of local edges into eowners*/
  ierr = MPI_Allgather(&l_nedges, 1, MPIU_INT, eowners+1, 1, MPIU_INT, PETSC_COMM_WORLD);CHKERRQ(ierr);

  eowners[0] = 0;
  for(i=2; i <= size; ++i){
    eowners[i] += eowners[i-1];
  }
  /*get local edge range */
  estart = eowners[rank];
  eend = eowners[rank+1];
  /* distribute edgelist */
  if(!rank){
    for(i=1; i < size; ++i){
      ierr = MPI_Send(edgelist + 2 * eowners[i], 2 * (eowners[i+1] - eowners[i]),
		      MPIU_INT, i, tag, comm);CHKERRQ(ierr);
    }
  } else {
    MPI_Status status;
    ierr = PetscMalloc1(2*(eend-estart), &edgelist);CHKERRQ(ierr);
    ierr = MPI_Recv(edgelist, 2*(eend-estart), MPIU_INT, 0, tag, comm, &status);CHKERRQ(ierr);
  }

  network->edgelist = edgelist;
  /* compute global and local number of non-ghost vertices for each process (and send to everybody)*/
  if(!rank){
    for(i=0; i < size; ++i){
      for(e=eowners[i]; e < eowners[i+1]; ++e){
	v = edgelist[2*e];
	if(!vtx_done[v]){
	  ++(nvtx[i]);
	  vtx_done[v] = 1;
	}
	v = edgelist[2*e+1];
	if(!vtx_done[v]){
	  ++(nvtx[i]);
	  vtx_done[v] = 1;
	}
      }/*end for e*/
    }/*end for i*/
  }
  /*send it out*/
  ierr = MPI_Bcast(&g_nvertices, 1, MPIU_INT, 0, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Scatter(nvtx, 1, MPIU_INT, &l_nvertices, 1, MPIU_INT, 0, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = PetscFree3(eowners, nvtx, vtx_done);CHKERRQ(ierr);

  net->g_nvertices = g_nvertices;
  net->l_nvertices = l_nvertices;

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
  PetscBool          ts_distribute;

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
    nentry = 5;
    nexit = 5;
    nedges = 8;
    nvertices = 7;
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
      ninterchange = nvertex + 2;
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
	  vertices[i].intc_ctx = interchanges;
	  /* set for the edge too */
	  highways[0].interchanges = interchanges;
	  highways[7].interchanges = interchanges;
	  highways[0].length = 10.0;
	  
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
	  highways[i].entries = entries[i];
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
	  highways[7].exits = exits;
	  vertices[i].exit_ctx = exits[i];
	  highways[i].num_lanes = 4;
	  highways[i].speed_limit=70.0;
	  break;
	case 1:
	  /* V1, one interchange, one exit, one entry */
	  interchanges[i].exit_road_id = 0;
	  interchanges[i].entry_road_id = 1;
	  interchanges[i].postmile = 10.0;
	  
	  vertices[i].intc_ctx = interchanges + i;/* pointer arithmetic */
	  highways[1].interchanges = interchanges + i;
	  highways[1].length = 10.0;
	  highways[0].interchanges2 = interchanges + i;
	  
	  entries[i].postmile = 10.0;
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
	  highways[i].entries = entries[i];
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
	  highways[0].exits = exits + i;
	  highways[1].num_lanes = 3;
	  highways[1].speed_limit=70.0;
	  break;
	case 2:
	  /* V2, two interchanges, no entries or exits */
	  interchanges[i].exit_road_id = 1;
	  interchanges[i].entry_road_id = 2;
	  interchanges[i].postmile = 20.0;
	  interchanges[i].continue_p = 0.65;
	  
	  interchanges[3].exit_road_id = 1;
	  interchanges[3].entry_road_id = 4;
	  interchanges[3].postmile = 20.0;
	  interchanges[3].continue_p = 0.35;
	  
	  vertices[i].intc_ctx = interchanges + i;/* pointer arithmetic */
	  vertices[i].intc2_ctx = interchanges + 3;
	  
	  highways[2].interchanges = interchanges + 2;
	  highways[4].interchanges = interchanges + 3;
	  highways[1].interchanges2 = interchanges + 2;
	  highways[1].interchanges3 = interchanges + 3;
	  highways[2].num_lanes = 3;
	  highways[2].speed_limit = 80.0;
	  highways[4].num_lanes = 2;
	  highways[4].speed_limit = 55.0;
	  highways[2].length = 8.0;
	  highways[4].length = 12.0;
	  break;
	case 3:
	  /* V3, one interchange, one exit, one entry */
	  interchanges[4].exit_road_id = 2;
	  interchanges[4].entry_road_id = 3;
	  interchanges[4].postmile = 28.0;
	  vertices[i].intc_ctx = interchanges + 4;/* pointer arithmetic */
	  highways[2].interchanges2 = interchanges2 = interchanges + 4;
	  highways[3].interchanges = interchanges + 4;
	  highways[3].length = 12.0;
	  highways[3].num_lanes = 3;
	  highways[3].speed_limit = 65.0;
	  
	  entries[2].postmile = 28.0;
	  entries[2].arrival_dist = TS_POISSON_DYNAMIC;
	  /* the N arrival parameters for TS_POISSON_DYNAMIC are given as
	     an array of length 2N with entries
	     [t0,lambda0, t1, lambda1, ..., tN-1, lambdaN-1],
	     where the ti are the times at which lambdai becomes the arrival
	     distribution's parameter. */
	  entries[2].num_arrival_params = 4;
	  arrival_params[8 * 2 + 0] = 0.0; arrival_params[8 * 2 + 1] = 10.0;
	  arrival_params[8 * 2 + 2] = 2.5; arrival_params[8 * 2 + 3] = 1.01; /* 2.5 hour morning rush hour */
	  arrival_params[8 * 2 + 4] = 5.0; arrival_params[8 * 2 + 5] = 9.0;
	  arrival_params[8 * 2 + 6] = 9.0; arrival_params[8 * 2 + 7] = 3.0;
	  entries[2].arrival_params = arrival_params + 8 * 2;/* hey, we're adding some pointer arithmetic to this unholy assignment protocol we're using for the arrival/exit parameters. Why? Because the universe hates you, that's why.*/
	  vertices[i].entr_ctx = entries[2];
	  vertices[i].arrival_dist = TS_POISSON_DYNAMIC;
	  highways[3].entries = entries + 2;
	  
	  exits[2].postmile = 28.0;
	  exits[2].type = TS_DYNAMIC_EXIT;
	  /* same deal as with arrival params, but the odd-numbered-indexed values are probabilities
	     between 0 and 1 (inclusive) that any given car will leave the highway on this exit 
	     in the given timeframe.*/
	  vertices[i].exit_dist = TS_DYNAMIC_EXIT;
	  exit_params[8 * 2 + 0] = 0.0; exit_params[8 * 2 + 1] = 0.08;
	  exit_params[8 * 2 + 2] = 2.5; exit_params[8 * 2 + 3] = 0.005;
	  exit_params[8 * 2 + 4] = 5.0; exit_params[8 * 2 + 5] = 0.24;
	  exit_params[8 * 2 + 6] = 9.0; exit_params[8 * 2 + 7] = 0.04;
	  epars3.n = 4;
	  epars3.params = exit_params + 8 * 2;
	  exits[2].prob_params = epars3;
	  highways[2].exits = exits + 2;
	  vertices[i].exit_ctx = exits[2];
	  break;
	 case 4:
	  /* V4, one interchange, one exit, one entry */
	  interchanges[5].exit_road_id = 3;
	  interchanges[5].entry_road_id = 5;
	  interchanges[5].postmile = 40.0;
	  vertices[i].intc_ctx = interchanges + 5;/* pointer arithmetic */
	  highways[3].interchanges2 = interchanges + 5;
	  highways[5].interchanges = interchanges;
	  highways[5].length = 8.0;
	  highways[5].num_lanes = 3.0;
	  highways[5].speed_limit = 65.0;
	  
	  entries[3].postmile = 40.0;
	  entries[3].arrival_dist = TS_POISSON_DYNAMIC;
	  /* the N arrival parameters for TS_POISSON_DYNAMIC are given as
	     an array of length 2N with entries
	     [t0,lambda0, t1, lambda1, ..., tN-1, lambdaN-1],
	     where the ti are the times at which lambdai becomes the arrival
	     distribution's parameter. */
	  entries[3].num_arrival_params = 4;
	  arrival_params[8 * 3 + 0] = 0.0; arrival_params[8 * 3 + 1] = 10.0;
	  arrival_params[8 * 3 + 2] = 2.5; arrival_params[8 * 3 + 3] = 1.01; /* 2.5 hour morning rush hour */
	  arrival_params[8 * 3 + 4] = 5.0; arrival_params[8 * 3 + 5] = 9.0;
	  arrival_params[8 * 3 + 6] = 9.0; arrival_params[8 * 3 + 7] = 3.0;
	  entries[3].arrival_params = arrival_params + 8 * 3;/* hey, we're adding some pointer arithmetic to this unholy assignment protocol we're using for the arrival/exit parameters. Why? Because the universe hates you, that's why.*/
	  vertices[i].entr_ctx = entries[3];
	  highways[5].entries = entries + 3;
	  vertices[i].arrival_dist = TS_POISSON_DYNAMIC;

	  exits[3].postmile = 40.0;
	  exits[3].type = TS_DYNAMIC_EXIT;
	  /* same deal as with arrival params, but the odd-numbered-indexed values are probabilities
	     between 0 and 1 (inclusive) that any given car will leave the highway on this exit 
	     in the given timeframe.*/
	  vertices[i].exit_dist = TS_DYNAMIC_EXIT;
	  exit_params[8 * 3 + 0] = 0.0; exit_params[8 * 3 + 1] = 0.08;
	  exit_params[8 * 3 + 2] = 2.5; exit_params[8 * 3 + 3] = 0.005;
	  exit_params[8 * 3 + 4] = 5.0; exit_params[8 * 3 + 5] = 0.24;
	  exit_params[8 * 3 + 6] = 9.0; exit_params[8 * 3 + 7] = 0.04;
	  epars4.n = 4;
	  epars4.params = exit_params + 8 * 3;
	  exits[3].prob_params = epars4;
	  highways[3].exits = exits + 3;
	  vertices[i].exit_ctx = exits[3];
	  break;
	case 5:
	/* V5, two interchanges, no exits or entries */
	  interchanges[6].exit_road_id = 6;
	  interchanges[6].entry_road_id = 5;
	  interchanges[6].postmile = 48.0;
	  vertices[i].intc_ctx = interchanges + 6;/* pointer arithmetic */
	  highways[5].interchanges2 = interchanges + 6;

	  interchanges[7].exit_road_id = 6;
	  interchanges[7].entry_road_id = 4;
	  interchanges[7].postmile = 48.0;
	  vertices[i].intc2_ctx = interchanges + 7;/* pointer arithmetic */
	  highways[4].interchanges2 = interchanges + 7;
	  highways[6].interchanges = interchanges + 6;
	  highways[6].interchanges2 = interchanges + 7;
	  highways[6].num_lanes = 3;
	  highways[6].length = 20.0;
	  highways[6].speed_limit = 65.0
	  break;
	 case 6:
	  /* V6, one interchange, one exit, one entry */
	  interchanges[8].exit_road_id = 3;
	  interchanges[8].entry_road_id = 5;
	  interchanges[8].postmile = 68.0;
	  vertices[i].intc_ctx = interchanges + 8;/* pointer arithmetic */
	  highways[7].interchanges = interchanges + 8;
	  highways[7].num_lanes = 4;
	  highways[7].length = 12.0;
	  
	  entries[4].postmile = 68.0;
	  entries[4].arrival_dist = TS_POISSON_DYNAMIC;
	  /* the N arrival parameters for TS_POISSON_DYNAMIC are given as
	     an array of length 2N with entries
	     [t0,lambda0, t1, lambda1, ..., tN-1, lambdaN-1],
	     where the ti are the times at which lambdai becomes the arrival
	     distribution's parameter. */
	  entries[4].num_arrival_params = 4;
	  arrival_params[8 * 4 + 0] = 0.0; arrival_params[8 * 4 + 1] = 10.0;
	  arrival_params[8 * 4 + 2] = 2.5; arrival_params[8 * 4 + 3] = 1.01; /* 2.5 hour morning rush hour */
	  arrival_params[8 * 4 + 4] = 5.0; arrival_params[8 * 4 + 5] = 9.0;
	  arrival_params[8 * 4 + 6] = 9.0; arrival_params[8 * 4 + 7] = 3.0;
	  entries[4].arrival_params = arrival_params + 8 * 3;/* hey, we're adding some pointer arithmetic to this unholy assignment protocol we're using for the arrival/exit parameters. Why? Because the universe hates you, that's why.*/
	  vertices[i].entr_ctx = entries[4];
	  vertices[i].arrival_dist = TS_POISSON_DYNAMIC;
	  highways[7].entries = entries + 4;
	  highways[7].exits = exits;
	  exits[4].postmile = 68.0;
	  exits[4].type = TS_DYNAMIC_EXIT;
	  /* same deal as with arrival params, but the odd-numbered-indexed values are probabilities
	     between 0 and 1 (inclusive) that any given car will leave the highway on this exit 
	     in the given timeframe.*/
	  vertices[i].exit_dist = TS_DYNAMIC_EXIT;
	  exit_params[8 * 4 + 0] = 0.0; exit_params[8 * 4 + 1] = 0.08;
	  exit_params[8 * 4 + 2] = 2.5; exit_params[8 * 4 + 3] = 0.005;
	  exit_params[8 * 4 + 4] = 5.0; exit_params[8 * 4 + 5] = 0.24;
	  exit_params[8 * 4 + 6] = 9.0; exit_params[8 * 4 + 7] = 0.04;
	  epars5.n = 4;
	  epars5.params = exit_params + 8 * 4;
	  exits[4].prob_params = epars5;

	  highways[6].exits = exits + 4;
	  vertices[i].exit_ctx = exits[4];
	  break;

        default:
	  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong value of i in TSNetworkCreateWithStructure (from line 385).");

	}
      }/*end vertex for*/
      net->highways = highways;
      net->vertices = vertices;
      
     }/*if(!rank)*/
    /* global highway ids */


  default:
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong value of i in TSNetworkCreateWithStructure (from line 394).");
    } /* end network case switch */
 
  for(i=0; i < nedges; ++i){
    highways[i].id = i;
  }

  *network = net;
   net->g_nedges = nedges;
   if(!rank){
     net->l_nedges = nedges;
     net->l_nvertices = nvertices;
   }
   net->g_nvertices = nvertices;

   ierr = PetscOptionsGetBool(NULL, NULL, "-ts_distribute", &ts_distribute, NULL);CHKERRQ(ierr);
   if(ts_distribute){
     ierr = TSNetworkDistribute(comm, net);CHKERRQ(ierr);
   }
     


  PetscFunctionReturn(0);
}
  
  
