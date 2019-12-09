#include <tsnetwork.h>
#include <tshighway.h>
#include <petscdmnetwork.h>





int main(int argc, char** argv)
{
  PetscErrorCode   ierr;
  TSNetwork        network;
  DMNetworkMonitor monitor;
  DM               netdm;
  TS               ts;
  PetscInt         max_steps, nodes_per_highway, estart, eend, i, v, vstart, vend;
  PetscReal        ht;
  PetscBool        flag, manual_jacobian=PETSC_TRUE;
  TSHighway        highway;
  TSHighwayVertex  vertex;


  ierr = PetscInitialize(&argc, &argv, NULL, NULL);if(ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL, NULL, "-nodesperhighway", &nodes_per_highway, &flag);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL, NULL, "-max_steps", &max_steps, &flag);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &ht, &flag);CHKERRQ(ierr);
  ierr = TSNetworkCreateWithStructure(&network, &netdm, 0, nodes_per_highway, NULL);CHKERRQ(ierr);
  
  ierr = TSCreate(PETSC_COMM_WORLD, &ts);CHKERRQ(ierr);
  ierr = TSSetDM(ts, (DM)netdm);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts, NULL, TSHighwayNetIFunction, network);CHKERRQ(ierr);
  ierr = TSSetMaxSteps(ts, max_steps);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts, ht);CHKERRQ(ierr);
  ierr = TSSetType(ts, TSBEULER);CHKERRQ(ierr);

  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ierr = TSNetworkSetSolution(network, netdm, network->g_X);CHKERRQ(ierr);

  ierr = TSSolve(ts, network->g_X);CHKERRQ(ierr);

  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = VecDestroy(&network->g_X);CHKERRQ(ierr);
  ierr = VecDestroy(&network->l_X);CHKERRQ(ierr);
  ierr = VecDestroy(&network->l_dXdt);CHKERRQ(ierr);

  ierr = DMNetworkGetEdgeRange(netdm, &estart, &eend);CHKERRQ(ierr);
  for(i=estart; i < eend; ++i){
    ierr = DMNetworkGetComponent(netdm, i, 0, &network->highway_key, (void**)&highway);CHKERRQ(ierr);
    ierr = HighwayDestroyJacobian(highway);CHKERRQ(ierr);
    ierr = PetscFree(highway->old_rho_v);CHKERRQ(ierr);
    ierr = VecDestroy(&highway->X);CHKERRQ(ierr);
  }
  ierr = DMNetworkGetVertexRange(netdm, &vstart, &vend);CHKERRQ(ierr);
  for(v=vstart; v<vend; ++v){
    ierr = DMNetworkGetComponent(netdm, v, 0, &network->vertex_key, (void**)&vertex);CHKERRQ(ierr);
    for(i=0; i < 3; ++i){
      ierr = MatDestroy(&vertex->jac[i]);CHKERRQ(ierr);
    }
  }
  ierr = DMDestroy(&netdm);CHKERRQ(ierr);
  ierr = PetscFree(network);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
    
  
