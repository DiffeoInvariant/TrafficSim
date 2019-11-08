#include "../include/network.h"

extern PetscErrorCode SetRoadCounty(RoadCtx* ctx, const char* county_name)
{
  /* county_name MUST be three characters only */
  PetscFunctionBegin;
  if(county_name == NULL){
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_NULL, "county_name must be a non-null, three-character char array.");
  }
  PetscInt namelen = sizeof(county_name) / sizeof(county_name[0]);
  if(namelen != COUNTY_NAME_LEN){
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "county_name must contain precisely three characters.");
  }

  ctx->county = county_name;

  PetscFunctionReturn(0);
}
