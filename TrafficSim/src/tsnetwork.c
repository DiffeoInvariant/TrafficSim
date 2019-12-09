#include "../include/tsnetwork.h"
#include "../include/tsdensity.h"
#include "../include/tshighway.h"
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
  ierr = VecDestroy(&network->l_X);CHKERRQ(ierr);
  ierr = VecDestroy(&network->l_dXdt);CHKERRQ(ierr);
  ierr = DMDestroy(&network->network);CHKERRQ(ierr);
  ierr = PetscFree(network);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode TSNetworkCleanUp(TSNetwork network){
  PetscErrorCode ierr;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(network->comm, &rank);CHKERRQ(ierr);
  ierr = PetscFree(network->edgelist);CHKERRQ(ierr);
  
  if(!rank){
    ierr = PetscFree2(network->vertices,network->highways);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


PetscErrorCode VertexGetArrivalRate(TSHighwayEntryCtx *vertex_ctx, PetscReal t, PetscReal* lambda)
{
  TSHighwayEntryCtx       *entr;
  PetscInt                n, i;

  PetscFunctionBegin;
  if(vertex_ctx){
    entr = vertex_ctx;
    if(entr->arrival_dist == TS_POISSON_DYNAMIC){
      n = entr->num_arrival_params;
      for(i=0; i < 2*n - 3; i += 2){
	/* find the time */
	if(entr->arrival_params[i] <= t && entr->arrival_params[i+2] > t){
	  *lambda = entr->arrival_params[i+1];
	} else if(entr->arrival_params[i+2] == t){
	  *lambda = entr->arrival_params[i+3];
	}
      }
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Arrival distributions other than TS_POISSON_DYNAMIC are not currently supported.");
    }

  } else {
    *lambda = 0.0;
  }

  PetscFunctionReturn(0);
}


PetscErrorCode VertexGetExitRate(TSHighwayExitCtx* vertex_ctx, PetscReal t, PetscReal* erate)
{
  TSExitParams     epars;
  PetscInt         n, i;

  PetscFunctionBegin;

  if(vertex_ctx){
    if(vertex_ctx->type == TS_DYNAMIC_EXIT){
      epars = vertex_ctx->prob_params;
      n = epars.n;
      for(i=0; i < 2*n - 3; i += 2){
	/* find the time */
	if(epars.params[i] <= t && epars.params[i+2] > t){
	  *erate = epars.params[i+1];
	} else if(epars.params[i+2] == t || (i == 2*n-4 && epars.params[i+2] < t)){
	  *erate = epars.params[i+3];
	}
      }
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Exit distributions other than TS_DYNAMIC_EXIT are not currently supported.");
    }
  } else {
    *erate = 0.0;
  }

  PetscFunctionReturn(0);
}
      
      
PetscErrorCode VertexGetTrafficFlowRate(PetscReal rho, PetscReal v, PetscReal* tflowrate)
{
  PetscFunctionBegin;
  if(!tflowrate) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "tflowrate cannot be NULL.");
  *tflowrate = rho * v;
  PetscFunctionReturn(0);
}
  
    
PetscErrorCode VertexGetCurrentTrafficCount(PetscReal rho, PetscReal v, PetscReal deltat, PetscInt *count)
{
  PetscErrorCode ierr;
  PetscReal      flowrate;

  PetscFunctionBegin;
  ierr = VertexGetTrafficFlowRate(rho, v, &flowrate);CHKERRQ(ierr);

  *count = (PetscInt)(flowrate * deltat);

  PetscFunctionReturn(0);
}
  
	    
PetscErrorCode VertexGetNumArrivals(TSHighwayEntryCtx* vertex_ctx, PetscReal t, PetscInt *narrival)
{
  PetscErrorCode ierr;
  PetscReal      lambda;

  PetscFunctionBegin;

  ierr = VertexGetArrivalRate(vertex_ctx, t, &lambda);CHKERRQ(ierr);
  ierr = TSPoissonSample(lambda, narrival);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode VertexGetDensityFactor(TSHighwayEntryCtx* vertex_ctx, TSHighwayExitCtx* vertex_exctx, PetscReal rho, PetscReal v, PetscReal t, PetscReal deltat, PetscReal* density_fact)
{
  PetscErrorCode ierr;
  PetscInt       ncurr, nnew, ntot;
  PetscReal      new_density_factor, erate;

  PetscFunctionBegin;

  ierr = VertexGetCurrentTrafficCount(rho, v, deltat, &ncurr);CHKERRQ(ierr);
  ierr = VertexGetNumArrivals(vertex_ctx, t, &nnew);CHKERRQ(ierr);
  ierr = VertexGetExitRate(vertex_exctx, t, &erate);CHKERRQ(ierr);
  ntot = (PetscInt)((1.0 - erate)*(PetscReal)ncurr) + nnew;

  new_density_factor = ((PetscReal)ntot)/((PetscReal)ncurr);

  *density_fact = new_density_factor;

  PetscFunctionReturn(0);
}
    
    
    

PetscErrorCode TSNetworkDistribute(MPI_Comm comm, TSNetwork network, PetscBool manual_jacobian)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank, size, tag=0;
  PetscInt       i, e, v, l_nedges, l_nvertices, g_nedges, g_nvertices, estart, eend, vstart, vend;
  PetscInt       *eowners, *edgelists[1]={NULL}, *edgelist = network->edgelist, *nvtx=NULL, *vtx_done=NULL;
  TSHighway       highway=NULL,highways=NULL;
  TSHighwayVertex vertex=NULL,vertices=NULL;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);

  g_nedges = network->l_nedges;
  g_nvertices = network->l_nvertices;
  /* everybody knows the number of local and global edges (and we store extra edges on rank 0) */
  ierr = MPI_Bcast(&g_nedges, 1, MPIU_INT, 0, comm);CHKERRQ(ierr);
  l_nedges = g_nedges/size;
  if(!rank){
    l_nedges += g_nedges - size * (g_nedges/size);
  }
  network->g_nedges = g_nedges;
  network->l_nedges = l_nedges;

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

  network->g_nvertices = g_nvertices;
  network->l_nvertices = l_nvertices;

  edgelists[0] = edgelist;
  ierr = DMNetworkSetSizes(network->network, 1, &l_nvertices, &l_nedges, 0, NULL);CHKERRQ(ierr);
  ierr = DMNetworkSetEdgeList(network->network, edgelists, NULL);CHKERRQ(ierr);
  ierr = DMNetworkLayoutSetUp(network->network);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(network->network, &vstart, &vend);CHKERRQ(ierr);
  /*if NOT root, need to allocate memory for vertices*/
  if(rank){
    ierr = PetscCalloc2(vend-vstart, &vertices, l_nedges, &highways);CHKERRQ(ierr);
  }

  highways = network->highways;
  vertices = network->vertices;
  
  ierr = DMNetworkGetEdgeRange(network->network, &estart, &eend);CHKERRQ(ierr);
  for(e = estart; e < eend; ++e){
    ierr = DMNetworkAddComponent(network->network, e, network->highway_key, &highways[e - estart]);CHKERRQ(ierr);
    /* discrete_dimension is the number of nodes in one highway's discretization */
    ierr = DMNetworkAddNumVariables(network->network, e, 2 * highways[e-estart].discrete_dimension);CHKERRQ(ierr);
  }

  for(v = vstart; v < vend; ++v){
    ierr = DMNetworkAddComponent(network->network, v, network->vertex_key, &vertices[v-vstart]);CHKERRQ(ierr);

    ierr = DMNetworkAddNumVariables(network->network, v, 2);CHKERRQ(ierr);
  }

  if(size > 1){
    DM               plex;
    PetscPartitioner part;
    ierr = DMNetworkGetPlex(network->network, &plex);CHKERRQ(ierr);
    ierr = DMPlexGetPartitioner(plex, &part);CHKERRQ(ierr);
    ierr = PetscPartitionerSetType(part, PETSCPARTITIONERSIMPLE);CHKERRQ(ierr);
  }

  ierr = DMSetUp(network->network);CHKERRQ(ierr);
  ierr = TSNetworkCleanUp(network);CHKERRQ(ierr);

  /* distribute data and create vectors*/
  ierr = DMNetworkDistribute(&(network->network), 0);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(network->network, &network->g_X);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(network->network, &network->l_X);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(network->network, &network->l_dXdt);CHKERRQ(ierr);

  /* set up highways after distribution*/
  ierr = DMNetworkGetVertexRange(network->network, &vstart, &vend);CHKERRQ(ierr);

  ierr = DMNetworkHasJacobian(network->network, manual_jacobian, manual_jacobian);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(network->network, &estart, &eend);CHKERRQ(ierr);
  for(e=estart; e < eend; ++e){
    ierr = DMNetworkGetComponent(network->network, e, 0, &network->highway_key, (void**)&highway);CHKERRQ(ierr);
    ierr = HighwaySetUp(highway);CHKERRQ(ierr);

    if(manual_jacobian){
      Mat *jac;
      ierr = HighwayCreateJacobian(highway, NULL, &jac);CHKERRQ(ierr);
      ierr = DMNetworkEdgeSetMatrix(network->network, e, jac);CHKERRQ(ierr);
    }
  }/* end loop over edges */

  if(manual_jacobian){
    ierr = DMNetworkGetVertexRange(network->network, &vstart, &vend);CHKERRQ(ierr);
    Mat *jac;
    for(v = vstart; v < vend; ++v){
      ierr = DMNetworkGetComponent(network->network, v, 0, &network->vertex_key, (void**)&vertex);CHKERRQ(ierr);
      ierr = VertexCreateJacobian(network->network, vertex, v, NULL, &jac);CHKERRQ(ierr);
      ierr = DMNetworkVertexSetMatrix(network->network, v, jac);CHKERRQ(ierr);
      
      vertex->jac = jac;
    }/* end loop over vertices */
  }

  network->manual_jacobian = manual_jacobian;
      
  PetscFunctionReturn(0);
}


PetscErrorCode TSHighwayNetIFunction(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void* ctx)
{
  /* IFunction for d(rho)/dt + u(rho) * d(rho)/dx = (sum over vertices(cars entering - cars leaving)) */
  PetscErrorCode ierr;
  TSNetwork      net = (TSNetwork)ctx;
  DM             netdm;
  Vec            lX, lXdot, lF, lXold;
  const PetscInt *conn_comp;
  PetscInt       vbehind, vahead, offset_behind, offset_ahead, var_offset;
  PetscInt       v, vstart, vend, e, estart, eend, nend;
  PetscBool      ghost;
  PetscScalar    *farr, *vtxf, *hwyf;
  TSHighway       highway;
  TSHighwayVertex vertex;
  DMDALocalInfo  info;
  TSHighwayTrafficField *hwyx, *hwyxdot, *vtxx;
  const PetscScalar *xarr, *xdotarr, *xoldarr;
  PetscReal      dt, density_factor;

  PetscFunctionBegin;

  lX = net->l_X;
  lXdot = net->l_dXdt;

  ierr = TSGetSolution(ts, &lXold);CHKERRQ(ierr);
  ierr = TSGetDM(ts, &netdm);CHKERRQ(ierr);
  ierr = TSGetTimeStep(ts, &dt);CHKERRQ(ierr);
  ierr = DMGetLocalVector(netdm, &lF);CHKERRQ(ierr);

  /* zero rhs result before computing it */
  ierr = VecSet(F, 0.0);CHKERRQ(ierr);
  ierr = VecSet(lF, 0.0);CHKERRQ(ierr);

  
  /* update ghost values for lX and lXdot */
  ierr = DMGlobalToLocalBegin(netdm, X, INSERT_VALUES, lX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(netdm, X, INSERT_VALUES, lX);CHKERRQ(ierr);
  
  ierr = DMGlobalToLocalBegin(netdm, Xdot, INSERT_VALUES, lXdot);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(netdm, Xdot, INSERT_VALUES, lXdot);CHKERRQ(ierr);

  ierr = VecGetArrayRead(lX, &xarr);CHKERRQ(ierr);
  ierr = VecGetArrayRead(lXdot, &xdotarr);CHKERRQ(ierr);
  ierr = VecGetArrayRead(lXold, &xoldarr);CHKERRQ(ierr);
  ierr = VecGetArray(lF, &farr);CHKERRQ(ierr);

  /* get vertex information */
  ierr = DMNetworkGetVertexRange(netdm, &vstart, &vend);CHKERRQ(ierr);
  for(v = vstart; v < vend; ++v){
    ierr = DMNetworkIsGhostVertex(netdm, v, &ghost);CHKERRQ(ierr);
    if(ghost) continue;

    ierr = DMNetworkGetComponent(netdm, v, 0, &net->vertex_key, (void**)&vertex);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(netdm, v, &var_offset);CHKERRQ(ierr);

    /* value of x = (rho, v) at the vertex */
    vtxx = (TSHighwayTrafficField*)(xarr + var_offset);
    /* value of F at the vertex */
    vtxf = (PetscScalar*)(farr + var_offset);
    /* drho/dt: NOTE: THIS FUNCTION (VertexGetDensityFactor) IS WHERE THE RANDOMNESS
     HAPPENS.*/
    ierr = VertexGetDensityFactor(vertex->entr_ctx, vertex->exit_ctx, vtxx[0].rho, vtxx[0].v, t, dt, &density_factor);CHKERRQ(ierr);

    vtxf[0] = (density_factor - 1.0) * vtxx[0].rho;
    if(net->speed_model == TS_LINEAR){
      /* dv/dt*/
      vtxf[1] = vtxf[0] * -1 * vertex->speed_limit / vertex->rho_limit;
    }

  }/*end loop over vertices */

  ierr = DMNetworkGetEdgeRange(netdm, &estart, &eend);CHKERRQ(ierr);
  for(e = estart; e < eend; ++e){
    ierr = DMNetworkGetComponent(netdm, e, 0, &net->highway_key, (void**)&highway);CHKERRQ(ierr);

    ierr = DMNetworkGetVariableOffset(netdm, e, &var_offset);CHKERRQ(ierr);

    hwyx = (TSHighwayTrafficField*)(xarr + var_offset);
    hwyxdot = (TSHighwayTrafficField*)(xdotarr + var_offset);
    hwyf = (PetscScalar*)(farr + var_offset);

    highway->dt = dt;
    highway->old_rho_v = (TSHighwayTrafficField*)(xoldarr + var_offset);

    /*boundary values from connected vertices */
    ierr = DMNetworkGetConnectedVertices(netdm, e, &conn_comp);CHKERRQ(ierr);
    vbehind = conn_comp[0];
    vahead = conn_comp[1];
    ierr = DMNetworkGetVariableOffset(netdm, vbehind, &offset_behind);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(netdm, vahead, &offset_ahead);CHKERRQ(ierr);

    ierr = DMNetworkGetComponent(netdm, vbehind, 0, &net->vertex_key, (void**)&vertex);CHKERRQ(ierr);
    vtxx = (TSHighwayTrafficField*)(xarr + offset_behind);
    vtxf = (PetscScalar*)(farr + offset_behind);

    ierr = DMDAGetLocalInfo(highway->da, &info);CHKERRQ(ierr);
    ierr = HighwayLocalIFunction_LaxFriedrichs(highway, &info, t, hwyx, hwyxdot, hwyf,
					       vertex->rho, vertex->v);CHKERRQ(ierr);
    
    /* evaluate behind boundary */
    hwyf[0] = hwyx[0].rho - vtxx[0].rho;
    vtxf[0] -= hwyx[0].rho * hwyx[0].v;

    ierr = DMNetworkGetComponent(netdm, vahead, 0, &net->vertex_key, (void**)&vertex);CHKERRQ(ierr);
    vtxx = (TSHighwayTrafficField*)(xarr + offset_ahead);
    vtxf = (PetscScalar*)(farr + offset_ahead);
    nend = highway->discrete_dimension - 1;
    /*evaluate ahead boundary */
    hwyf[2*nend] = hwyx[nend].rho - vtxx[0].rho;
    vtxf[0] += hwyx[2*nend].rho * hwyx[2*nend].v;
  }/*end loop over edges */

  ierr = VecRestoreArrayRead(lX, &xarr);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(lXdot, &xdotarr);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(lXold, &xoldarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(lF, &farr);CHKERRQ(ierr);
  
  ierr = DMLocalToGlobalBegin(netdm, lF, ADD_VALUES, F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(netdm, lF, ADD_VALUES, F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(netdm, &lF);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
    

    
    
				  
  


PetscErrorCode TSNetworkCreateWithStructure(TSNetwork* network, DM* netdm, PetscInt network_case, PetscInt nodes_per_highway, const char* filename)
{
  PetscErrorCode     ierr;
  PetscInt           nentry, nexit, ninterchange,nvertices, nedges;
  PetscInt           i, *edgelist;
  PetscMPIInt        size, rank;
  TSNetwork          net=NULL;
  TSHighway          highways=NULL;
  TSHighwayVertex    vertices=NULL;
  TSInterchangeCtx*  interchanges=NULL;
  TSHighwayEntryCtx* entries=NULL;
  TSHighwayExitCtx*  exits=NULL;
  /*TSRoadDirection    direction;*/
  PetscReal          *arrival_params=NULL, *exit_params=NULL;
  TSExitParams       epars1, epars2, epars3, epars4, epars5;
  PetscBool          ts_distribute, manual_jacobian;
  MPI_Comm           comm;

  PetscFunctionBegin;

  ierr = PetscObjectGetComm((PetscObject)(*netdm), &comm);CHKERRQ(ierr);

  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);
  ierr = PetscCalloc1(1, &net);CHKERRQ(ierr);

  ierr = DMNetworkCreate(PETSC_COMM_WORLD, netdm);CHKERRQ(ierr);
  

  ierr = DMNetworkRegisterComponent((*netdm), "HighwayStruct", sizeof(struct _ts_HighwayCtx), &(net->highway_key));CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent((*netdm), "HighwayVertexStruct", sizeof(struct _ts_HighwayVertex), &(net->vertex_key));CHKERRQ(ierr);
  
  net->comm = comm;
  net->network = *netdm;
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
      nedges = net->g_nedges;
      nvertices = net->g_nvertices;
      ierr = PetscCalloc1(2 * net->g_nedges, &edgelist);CHKERRQ(ierr);
      edgelist[0] = 0; edgelist[1] = 1; /* H0 */
      edgelist[2] = 1; edgelist[3] = 2; /*H1 */
      edgelist[4] = 2; edgelist[5] = 3; /* H2 */
      edgelist[6] = 3; edgelist[7] = 4; /* H3 */
      edgelist[8] = 2; edgelist[9] = 6; /*H4*/
      edgelist[10] = 4; edgelist[11] = 5; /*H5*/
      edgelist[12] = 5; edgelist[13] = 6; /*H6*/
      edgelist[14] = 6; edgelist[15] = 0; /* H7 */

      /* allocate highways and vertices, incl entries, exits, and interchanges */
      ninterchange = nvertices + 2;
      ierr = PetscCalloc2(nvertices, &vertices,
			  nedges, &highways);CHKERRQ(ierr);
      
      ierr = PetscCalloc5(ninterchange, &interchanges,
			  nentry, &entries,
			  nentry * 8, &arrival_params,
			  nexit * 8, &exit_params,
			  nexit, &exits);CHKERRQ(ierr);
      
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
          highways[0].discrete_dimension = nodes_per_highway;
	  
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
	  vertices[i].entr_ctx = &entries[i];
	  highways[i].entries = &entries[i];
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
	  vertices[i].exit_ctx = &exits[i];
	  highways[i].num_lanes = 4;
	  highways[i].speed_limit=70.0;
	  vertices[0].speed_limit = 70.0;
	  vertices[0].rho_limit = 10.0;
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
	  
	  vertices[i].entr_ctx = &entries[i];
	  highways[i].entries = &entries[i];
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
	  
	  vertices[i].exit_ctx = &exits[i];
	  highways[0].exits = exits + i;
	  highways[1].num_lanes = 3;
          highways[1].discrete_dimension = nodes_per_highway;
	  highways[1].speed_limit=70.0;
	  vertices[1].speed_limit=70.0;
	  vertices[1].rho_limit = 10.0;
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
	  vertices[i].speed_limit = 70.0;
	  vertices[i].rho_limit = 10.0;
	  highways[4].num_lanes = 2;
	  highways[4].speed_limit = 55.0;
	  highways[2].length = 8.0;
	  highways[4].length = 12.0;
          highways[2].discrete_dimension = nodes_per_highway;
          highways[4].discrete_dimension = nodes_per_highway;
	  break;
	case 3:
	  /* V3, one interchange, one exit, one entry */
	  interchanges[4].exit_road_id = 2;
	  interchanges[4].entry_road_id = 3;
	  interchanges[4].postmile = 28.0;
	  vertices[i].intc_ctx = interchanges + 4;/* pointer arithmetic */
	  highways[2].interchanges2 = interchanges + 4;
	  highways[3].interchanges = interchanges + 4;
	  highways[3].length = 12.0;
	  highways[3].num_lanes = 3;
	  highways[3].speed_limit = 65.0;
	  vertices[i].speed_limit = 65.0;
	  vertices[i].rho_limit = 8.0;
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
	  vertices[i].entr_ctx = entries + 2;
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
	  vertices[i].exit_ctx = exits +2;
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
          highways[5].discrete_dimension = nodes_per_highway;
	  vertices[i].speed_limit = 65.0;
	  vertices[i].rho_limit = 10.0;
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
	  vertices[i].entr_ctx = entries + 3;
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
	  vertices[i].exit_ctx = exits + 3;
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
	  highways[6].speed_limit = 65.0;
          highways[6].discrete_dimension = nodes_per_highway;
	  vertices[i].speed_limit = 65.0;
	  vertices[i].rho_limit = 10.0;
	  break;
	 case 6:
	  /* V6, one interchange, one exit, one entry */
	  interchanges[8].exit_road_id = 3;
	  interchanges[8].entry_road_id = 5;
	  interchanges[8].postmile = 68.0;
	  vertices[i].intc_ctx = interchanges + 8;/* pointer arithmetic */
	  highways[7].interchanges = interchanges + 8;
	  highways[7].num_lanes = 4;
	  highways[7].speed_limit = 70.0;
	  vertices[i].speed_limit = 70.0;
	  vertices[i].rho_limit = 10.0;
	  highways[7].length = 12.0;
          highways[7].discrete_dimension = nodes_per_highway;
	  
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
	  vertices[i].entr_ctx = entries + 4;
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
	  vertices[i].exit_ctx = exits + 4;
	  break;

        default:
	  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong value of i in TSNetworkCreateWithStructure (from line 826).");

	}
      }/*end vertex for*/
      net->highways = highways;
      net->vertices = vertices;
      
     }/*if(!rank)*/
    /* global highway ids */


  default:
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Wrong value of i in TSNetworkCreateWithStructure (from line 838).");
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
     ierr = PetscOptionsGetBool(NULL, NULL, "-manual_jacobian", &manual_jacobian, NULL);CHKERRQ(ierr);
     ierr = TSNetworkDistribute(comm, net, manual_jacobian);CHKERRQ(ierr);
   }
     
  PetscFunctionReturn(0);
}
  

PetscErrorCode VertexCreateJacobian(DM dm, TSHighwayVertex vertex, PetscInt v, Mat* Jpre, Mat *J[])
{
  PetscErrorCode ierr;
  Mat            *vjac;
  PetscInt       nedges, e, i, M, N, *rows, *cols;
  PetscBool      isself;
  const PetscInt *edges, *conn_comps;
  PetscScalar    *zeros;

  PetscFunctionBegin;
  
  ierr = DMNetworkGetSupportingEdges(dm, v, &nedges, &edges);CHKERRQ(ierr);
  if (nedges <= 0){
    SETERRQ2(PETSC_COMM_SELF,1,"%d vertex, nedges %d\n",v,nedges);
  }

  /* each connected edge gets two Jacobians, one vertex-edge and one vertex-vertex */
  ierr = PetscCalloc1(2 * nedges + 1, &vjac);CHKERRQ(ierr);

  /* this vertex gets a dense block of zeros */
  ierr = DMNetworkGetNumVariables(dm, v, &M);CHKERRQ(ierr);
  if(M != 2){
    SETERRQ1(PETSC_COMM_SELF, 1, "Number of vertex variables M != 2.", M);
  }
  ierr = PetscMalloc3(M, &rows, M, &cols, M * M, &zeros);CHKERRQ(ierr);
  ierr = PetscArrayzero(zeros, M * M);CHKERRQ(ierr);
  for(i=0; i < M; ++i){
    rows[i] = i;
  }

  for(e=0; e < nedges; ++e){
    /* create Jacobian for vertex v to edge e, J(v, e) */
    ierr = DMNetworkGetConnectedVertices(dm, edges[e], &conn_comps);CHKERRQ(ierr);
    isself = (v == conn_comps[0]) ? PETSC_TRUE : PETSC_FALSE;

    if(Jpre){
      if(isself){
	vjac[2 * e + 1] = Jpre[0];
      } else {
	vjac[2 * e + 1] = Jpre[1];
      }
      vjac[2 * e + 2] = Jpre[2];
      ierr = PetscObjectReference((PetscObject)(vjac[2 * e + 1]));CHKERRQ(ierr);
      ierr = PetscObjectReference((PetscObject)(vjac[2 * e + 2]));CHKERRQ(ierr);
    } else {
      ierr = MatCreate(PETSC_COMM_SELF, &vjac[2 * e + 1]);CHKERRQ(ierr);
      ierr = DMNetworkGetNumVariables(dm, edges[e], &N);CHKERRQ(ierr);
      ierr = MatSetSizes(vjac[2 * e + 1], PETSC_DECIDE, PETSC_DECIDE, M, N);CHKERRQ(ierr);
      ierr = MatSetFromOptions(vjac[2 * e + 1]);CHKERRQ(ierr);
      ierr = MatSetOption(vjac[2 * e + 1], MAT_STRUCTURE_ONLY, PETSC_TRUE);CHKERRQ(ierr);
      ierr = MatSeqAIJSetPreallocation(vjac[2 * e + 1], 2, NULL);CHKERRQ(ierr);
      /* J(e, v) */
      if(N){
	if(isself){
	  for(i=0; i < 2; ++i){
	    /* coupling to incoming edge */
	    cols[i] = i;
	  }
	} else {
	  /* coupling to outgoing edge */
	  cols[0] = N-2;
	  cols[1] = N-1;
	}
	ierr = MatSetValues(vjac[2 * e + 1], 2, rows, 2, cols, zeros, INSERT_VALUES);CHKERRQ(ierr);
      }
      ierr = MatAssemblyBegin(vjac[2 * e + 1], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(vjac[2 * e + 1], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      /* Jacobian between this vertex and connected vertices, which is
	 just zero, since all relevant information as to what's happening
	 at this vertex is contained solely in the incoming/outgoing edges
	 and the vertex. */
      ierr = MatCreate(PETSC_COMM_SELF, &vjac[2 * e + 2]);CHKERRQ(ierr);
      ierr = MatSetSizes(vjac[2 * e + 2], PETSC_DECIDE, PETSC_DECIDE, M, M);CHKERRQ(ierr);
      ierr = MatSetFromOptions(vjac[2 * e + 2]);CHKERRQ(ierr);
      ierr = MatSetOption(vjac[2 * e + 2], MAT_STRUCTURE_ONLY, PETSC_TRUE);CHKERRQ(ierr);
      ierr = MatSeqAIJSetPreallocation(vjac[2 * e + 2], 1, NULL);CHKERRQ(ierr);
      ierr = MatAssemblyBegin(vjac[2 * e + 2], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(vjac[2 * e + 2], MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      }/* end if(N) */
    }/*end loop over e */
    ierr = PetscFree3(rows, cols, zeros);CHKERRQ(ierr);
    
    vertex->jac = vjac;

    PetscFunctionReturn(0);
  
}

static PetscErrorCode TSNetworkHighwaySetSolution(TSHighway highway, PetscReal rho)
{
  PetscErrorCode ierr;
  TSHighwayTrafficField *x;
  PetscInt       i, start, n, end;

  PetscFunctionBegin;
  ierr = DMDAVecGetArray(highway->da, highway->X, &x);CHKERRQ(ierr);
  ierr = DMDAGetCorners(highway->da, &start, 0, 0, &n, 0, 0);CHKERRQ(ierr);

  end = start + n;
  for(i = start; i < end; ++i){
    x[i].rho = 0.0;
    x[i].v = highway->speed_limit;
  }

  ierr = DMDAVecRestoreArray(highway->da, highway->X, &x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
  
  

PetscErrorCode TSNetworkSetSolution(TSNetwork network, DM netdm, Vec X)
{
  PetscErrorCode ierr;
  PetscInt       i, n, vstart, vend, offsetstart, offsetend, var_offset;
  PetscInt       e, estart, eend;
  Vec            lX;
  PetscScalar    *xarr;
  TSHighway      highway;
  const PetscInt *conn_comp;
  const PetscScalar *xarrview;

  PetscFunctionBegin;

  ierr = VecSet(X, 0.0);CHKERRQ(ierr);
  ierr = DMGetLocalVector(netdm, &lX);CHKERRQ(ierr);
  ierr = VecGetArray(lX, &xarr);CHKERRQ(ierr);

  ierr = DMNetworkGetEdgeRange(netdm, &estart, &eend);CHKERRQ(ierr);
  for(e = estart; e < eend; ++e){
    ierr = DMNetworkGetVariableOffset(netdm, e, &var_offset);CHKERRQ(ierr);
    ierr = DMNetworkGetComponent(netdm, e, 0, &network->highway_key, (void**)&highway);CHKERRQ(ierr);
    ierr = TSNetworkHighwaySetSolution(highway, 0.0);CHKERRQ(ierr);

    ierr = VecGetSize(highway->X, &n);CHKERRQ(ierr);
    ierr = VecGetArrayRead(highway->X, &xarrview);CHKERRQ(ierr);

    for(i=0; i<n; ++i){
      (xarr+var_offset)[i] = xarrview[i];
    }

    ierr = DMNetworkGetConnectedVertices(netdm, e, &conn_comp);CHKERRQ(ierr);
    vstart = conn_comp[0];
    vend = conn_comp[1];

    ierr = DMNetworkGetVariableOffset(netdm, vstart, &offsetstart);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(netdm, vend, &offsetend);CHKERRQ(ierr);

    (xarr + offsetstart)[0] = 0.0;
    (xarr + offsetend)[0] = 0.0;

    ierr = VecRestoreArrayRead(highway->X, &xarrview);CHKERRQ(ierr);
  }/*end edge loop */
  
  ierr = DMLocalToGlobalBegin(netdm, lX, ADD_VALUES, X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(netdm, lX, ADD_VALUES, X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(netdm, &lX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
    
