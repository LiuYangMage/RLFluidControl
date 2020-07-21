/* ------------------------------------------------------------------------- *
 * StokesBC() - Calculate the high-order boundary conditions for Stokes flow *
 *                                                                           *
 * This routine simply sets the non-linear terms to zero.                    *
 *                                                                           *
 * RCS Information                                  
 * ---------------
 * $Author: ssherw $
 * $Date: 2006/05/13 09:32:45 $
 * $Source: /homedir/cvs/Nektar/Fourier/src/stokes.C,v $
 * $Revision: 1.2 $
 * ------------------------------------------------------------------------- */

#include "nektarF.h"

void StokesBC(Domain *omega){
  int      nq;
  
  /* Set the non-linear terms to zero */

  nq = omega->U->htot*omega->U->nz;
  dzero (nq, omega->Uf->base_h, 1);
  dzero (nq, omega->Vf->base_h, 1);
  dzero (nq, omega->Wf->base_h, 1);

  omega->Uf->Set_state('p');
  omega->Vf->Set_state('p');
  omega->Wf->Set_state('p');

  return;
}

