#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "nektarF.h"
#include "map.h"
#include <rfftw.h>
extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;

/* ------------------------------------------------------------------------ *
 *
 * INITIALIZATION
 * ==============
 *
 * SaveBCs()
 * saveBC()
 * ResetICs()
 * 
 * ------------------------------------------------------------------------ */

static int CountBCtype (Bndry *Xbc, char type);
static void BndryJbwd(Bndry *Bc, double *d);

#ifdef TO_DELETE
/* ------------------------------------------------------------------------ *
 * SaveBCs() - save all the velocity boundary conditions
 * ------------------------------------------------------------------------ */

static void SaveBCs (Domain *omega)
{
  omega->ubc1 = saveBC (omega->Ubc);
  omega->vbc1 = saveBC (omega->Vbc);
#if DIM == 3
  omega->wbc1 = saveBC (omega->Wbc);
#endif

  /* IFVERBOSE showsavebc (omega); */
  
  return;
}

/* ------------------------------------------------------------------------ *
 * saveBC() - save a velocity boundary condition of type 'V'
 * ------------------------------------------------------------------------ */

double *saveBC (Bedge *Xbc)
{
  int     np     = Xbc->bedge->l,
          nz     = Xbc->elmt->nz,
          ntot   = np * nz,                /* number of points per processor */
          NZ     = option("NZ"),
          ntot0  = np * NZ,                /* number of points on processor 0 */
          nbcs   = CountBCtype (Xbc, 'V'); /* number of saved BCs */
  double *data0, *tmp0,
         *data   = dvector (0, nbcs*ntot-1),
         *tmp    = dvector (0, NZ-1);
  register int i=0, j;

  /* copy vaules into data */
  for ( ; Xbc; Xbc = Xbc->next)
    if (Xbc->type == 'V') {
      dcopy (ntot, Xbc->bc.value, 1, data + i*ntot, 1);      
      i++;
    }

#ifndef SPMD
#if DIM == 3 /* Move to Physical space */
  for (j = 0; j < nbcs; j++)
    for (i = 0; i < np; i++) {
      dcopy  (NZ,   data + j*ntot + i, np, tmp, 1);
      realft (NZ/2, tmp, 1);
      dcopy  (NZ,   tmp, 1, data + j*ntot + i, np);
    }
#endif /* DIM == 3 */
  free (tmp);
  return data;
#else  /* ifndef SPMD */
  ROOTONLY {
    data0 = dvector (0, nbcs*ntot0-1); /* BC data is saved only on processor 0 */
    tmp0  = dvector (0, nbcs*ntot0-1);
  }
  dcopy_all0 (nbcs*ntot, data, data0);
  ROOTONLY {
    for (i = 0; i < nprocs; i++)
      for (j = 0; j < nbcs; j++)
	dcopy (ntot, data0 + i*ntot*nbcs + j*ntot,   1,
                      tmp0 + i*ntot + j*ntot*nprocs, 1);
    for (j = 0; j < nbcs; j++)
      for (i = 0; i < np; i++) {
	dcopy  (NZ,   tmp0 + j*ntot0 + i, np, tmp, 1);
	realft (NZ/2, tmp, 1);
	dcopy  (NZ,   tmp, 1, data0 + j*ntot0 + i, np);
      }
    free (tmp0);
  }
  free (data);
  free (tmp);
  return data0;
#endif /* ifndef SPMD */
}
#endif

/* ------------------------------------------------------------------------ *
 * ResetICs() -- reset velocity initial conditions                          *
 * ------------------------------------------------------------------------ */
void ResetICs (Domain *omega)
{
  Element_List *U = omega->U,
               *V = omega->V,
               *W = omega->W,
               *P = omega->P;
  #ifdef THERMO
  Element_List *T = omega->T;
  #endif

  // just to be safe
  U->Set_state('p');
  V->Set_state('p');
  W->Set_state('p');
  P->Set_state('p');
 #ifdef THERMO
  T->Set_state('p');
 #endif

  if (option("dealias")) { // Need to use 3/2 rule again
    U->Trans(U, P_to_F);
    U->Trans(U, F_to_P32);
    V->Trans(V, P_to_F);
    V->Trans(V, F_to_P32);
    W->Trans(W, P_to_F);
    W->Trans(W, F_to_P32);
   #ifdef THERMO
    T->Trans(T, P_to_F);
    T->Trans(T, F_to_P32);
   #endif
  }
  
  Mapfield (omega, 1);

  return;
}

/* ------------------------------------------------------------------------ *
 * 
 * UPDATING 
 * ========
 * 
 * gather_baseP() 
 * UpdateVbc()
 * MapField()
 * 
 * ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ *
 * gather_baseP() - gather base pressure: 
 * ------------------------------------------------------------------------ */

void gather_baseP (Domain *omega, double *basep)
{
  int idbasep = iparam("IDEBASEP") - 1,
    vertbasep = iparam("IDVBASEP") - 1;

  int nztot = omega->P->nztot,
         nz = omega->P->nz;
  double *tbasep = (double *) malloc (nz*sizeof(double));
  register int i;
  for(i = 0; i < nz; ++i)
    tbasep[i] = omega->P->flevels[i]->flist[idbasep]->vert[vertbasep].hj[0];
  MPI_Gather (tbasep, nz, MPI_DOUBLE, basep, nz, MPI_DOUBLE, 0, 
	      MPI_COMM_WORLD);
  ROOT rfftw(rplan_inv, 1, (FFTW_COMPLEX *) basep, 1, 0, 0, 0, 0);
  free (tbasep);
  return;
}

static void BndryJbwd(Bndry *Bc, double *d){
  Element *E = Bc->elmt;
  Basis   *b = E->getbasis();
  int      i;

  dsmul(E->qa, Bc->bvert[0], b->vert[0].a, 1, d, 1);
  daxpy(E->qa, Bc->bvert[1], b->vert[1].a, 1, d, 1);
  
  for(i = 0; i < E->edge[Bc->face].l;++i)
    daxpy(E->qa, Bc->bedge[0][i], b->edge[0][i].a, 1, d, 1);
}

/* ------------------------------------------------------------------------ *
 * UpdateVbc() -- update velocity boundary conditions
 * ------------------------------------------------------------------------ */

void UpdateVbc (Domain *omega)
{
#ifndef TMAP // for stationary geometries speeds need to be adjusted only once.
  static int updated = 0;
  if (updated) 
    return;
  else
    updated++; // next time
#endif

  // code uses lots of statics that don't need to be there for stationary
  // geometries but it's not worth cleaning up for the extra storage.
  Mapping     *mapx = omega->mapx,
          *mapy = omega->mapy;
  int nz    = omega->U->nz,
      nztot = omega->U->nztot;
  register int k;

#ifndef CYL_SHORTCUT // general treatment of the BCs
#

#ifndef VAR_BNDRY
  // Prepare for 'V' type BCs
  static int checked = 0, ell = 0, bcisV = 0;
  static double ***VlBC;
  if (!checked) {
    bcisV = CountBCtype (omega->Ubc[0], 'V'); // count them
    if (bcisV) {
      VlBC = (double ***) calloc(nz, sizeof(double **));
      for (k = 0; k < nz; k++)
	VlBC[k] = dmatrix(0, bcisV, 0, 5); // storage for the vertex values
    }
    checked++; // not to repeat check
  }
#else
  // To handle 'v' type BCs
  static double  *tmpu, *tmpv, *tmpw;
#endif

  Element_List *U = omega->U,
               *V = omega->V,
               *W = omega->W;
  dparam_set("t", mapy->time);

  // invFFT the maps to physical space
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
#ifdef TMAP // Need the mesh velocities
  // invFFT the maps to physical space
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
#endif  

  // Get the BCs in Physical Space
  TransformBCs(omega, F_to_P);
  register int Vl;

  for (k = 0; k < nz; k++) {
    Bndry   *Ubc  = omega->Ubc[k],
            *Vbc  = omega->Vbc[k],
            *Wbc  = omega->Wbc[k];
    int kloc = parid(k);
    Vl = 0;
    dparam_set("z", zmesh(k)); // Do 3 dparam_set() calls nz times at the most
    dparam_set("ORIGINX", mapx->d[kloc]);
    dparam_set("ORIGINY", mapy->d[kloc]);
    
    for(Bndry *uBc=Ubc, *vBc=Vbc, *wBc=Wbc; uBc; 
	uBc=uBc->next, vBc=vBc->next, wBc=wBc->next){
      if (uBc->type == 'V') {
	Element *UE, *VE, *WE;
	register int i;
	
#ifdef VAR_BNDRY
	if (!tmpu) {
	  tmpu = dvector(0, QGmax-1);
	  tmpv = dvector(0, QGmax-1);
	  tmpw = dvector(0, QGmax-1); // no de-aliasing here
	}
	// get the boundary elements
	UE = U->flevels[k]->flist[uBc->elmt->id];
	VE = V->flevels[k]->flist[vBc->elmt->id];
	WE = W->flevels[k]->flist[wBc->elmt->id];
	// update the BCs in modal space
	uBc->elmt->update_bndry(uBc,0);
	vBc->elmt->update_bndry(vBc,0);
	wBc->elmt->update_bndry(wBc,0);
	// get the quadrature space representation
	BndryJbwd(uBc, tmpu);
	BndryJbwd(vBc, tmpv);
	BndryJbwd(wBc, tmpw);
	
	for (i = 0; i < UE->qa; i++) {
#ifdef TMAP
	  tmpu[i] -= mapx->t[kloc] + mapx->z[kloc]*tmpw[i];
	  tmpv[i] -= mapy->t[kloc] + mapy->z[kloc]*tmpw[i];
#else
	  tmpu[i] -= mapx->z[kloc]*tmpw[i];
	  tmpv[i] -= mapy->z[kloc]*tmpw[i];
#endif
	}
	// set the vertex values
	uBc->bvert[0] = tmpu[0];
	uBc->bvert[1] = tmpu[UE->qa-1];
	
	vBc->bvert[0] = tmpv[0];
	vBc->bvert[1] = tmpv[UE->qa-1];
	
	uBc->elmt->JtransEdge(uBc,0,0,tmpu);
	vBc->elmt->JtransEdge(vBc,0,0,tmpv);
#else
	if (ell < bcisV*nz) {
	  VlBC[k][ell%bcisV][0] = uBc->bvert[0];
	  VlBC[k][ell%bcisV][1] = uBc->bvert[1];
	  VlBC[k][ell%bcisV][2] = vBc->bvert[0];
	  VlBC[k][ell%bcisV][3] = vBc->bvert[1];
	  VlBC[k][ell%bcisV][4] = wBc->bvert[0];
	  VlBC[k][ell%bcisV][5] = wBc->bvert[1];
	  ell++;
	}
	uBc->bvert[0] = VlBC[k][Vl][0] -mapx->z[kloc]*VlBC[k][Vl][4];
	uBc->bvert[1] = VlBC[k][Vl][1] -mapx->z[kloc]*VlBC[k][Vl][5];
	vBc->bvert[0] = VlBC[k][Vl][2] -mapy->z[kloc]*VlBC[k][Vl][4];
	vBc->bvert[1] = VlBC[k][Vl][3] -mapy->z[kloc]*VlBC[k][Vl][5];
#ifdef TMAP
	uBc->bvert[0] -= mapx->t[kloc];
	uBc->bvert[1] -= mapx->t[kloc];
	vBc->bvert[0] -= mapy->t[kloc];
	vBc->bvert[1] -= mapy->t[kloc];
#endif
	Vl++;
#endif
      }
    }
  }

  // Get the BCs in Fourier Space
  TransformBCs(omega, P_to_F);

  // FFT the maps back to fourier space
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
#ifdef TMAP
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
#endif
#
#else // ifdef CYL_SHORTCUT 
#
  // Case of cylinder flow with uniform inflow, U=(a,b,0)
  // and side BCs
  static double *prev_mapxt, *prev_mapyt;
  
  // no need to FFT anything - previous values of maps are
  // stored in the matrices below.
  if (!prev_mapxt) {
    prev_mapxt = dvector(0, nz-1); dzero(nz, prev_mapxt, 1);
    prev_mapyt = dvector(0, nz-1); dzero(nz, prev_mapyt, 1);
  }

  for (k = 0; k < nz; k++) {
    Bndry   *Ubc  = omega->Ubc[k],
            *Vbc  = omega->Vbc[k];
    int kloc = parid(k);

    for(Bndry *uBc=Ubc, *vBc=Vbc; uBc; uBc=uBc->next, vBc=vBc->next){
      if (uBc->type == 'V') {
	uBc->bvert[0] += prev_mapxt[k] - mapx->t[kloc];
	uBc->bvert[1] += prev_mapxt[k] - mapx->t[kloc];
	vBc->bvert[0] += prev_mapyt[k] - mapy->t[kloc];
	vBc->bvert[1] += prev_mapyt[k] - mapy->t[kloc];
      }
    }
    prev_mapxt[k] = mapx->t[kloc];
    prev_mapyt[k] = mapy->t[kloc];
  }
#endif
  return;
}

/* ------------------------------------------------------------------------ *
 * MapField() - computes U = U' - xt and V = V' - yt 
 *              (argument and result in Fourier -Modal or Quadrature- space)
 * ------------------------------------------------------------------------ */

void MapField (Domain *omega, int dir)
{
  Element_List *U = omega->U,
               *V = omega->V,
               *W = omega->W;
  Nek_Trans_Type f_to_p = F_to_P,
                 p_to_f = P_to_F;
  int revert = 0, mdir = dir;

  if (option("dealias")) {
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
  }

  // make sure velocity fields are in quadrature space
  if (U->fhead->state == 't') { // all 3 velocity fields should be same type
    U->Trans(U, J_to_Q);
    V->Trans(V, J_to_Q);
    W->Trans(W, J_to_Q);
    revert++;
  }

  U->Trans(U, f_to_p);
  V->Trans(V, f_to_p);
  W->Trans(W, f_to_p);

  if (dir == 0) mdir = -1;

  Mapfield (omega, mdir);

  U->Trans(U, p_to_f);
  V->Trans(V, p_to_f);
  W->Trans(W, p_to_f);

  if (revert && (dir != 0)) {
    U->Trans(U, Q_to_J);
    V->Trans(V, Q_to_J);
    W->Trans(W, Q_to_J);
  }

  return;
}

/* ------------------------------------------------------------------------ *
 * Mapfield() - computes U = U' - xt - xz W' and V = V' - yt - yz W'
 *              (argument and result in Physical Quadrature space)
 * ------------------------------------------------------------------------ */

void Mapfield (Domain *omega, int dir)
{
  Element_List *U = omega->U,
               *V = omega->V,
               *W = omega->W;

  Mapping       *mapx = omega->mapx,
            *mapy = omega->mapy;
		
  register int k, j, i, plane;
  rfftw_plan MPlan     = rplan, 
             MPlan_inv = rplan_inv;


  if (abs(dir) == 2) dir /= 2; // handle the extra cases due to MapField

  int NZ = U->nz;
  int NZTOT = U->nztot;
  int nZtot = NZTOT; // keep this value as it may change because of dealiasing

  if (option("dealias")) {
    NZ    = 3*NZ/2;
    NZTOT = 3*NZTOT/2;
    MPlan     = rplan32;    
    MPlan_inv = rplan_inv32;

    // zero out for dealiasing
    dzero (NZTOT-nZtot, mapx->z+nZtot, 1);
    dzero (NZTOT-nZtot, mapy->z+nZtot, 1);
#ifdef TMAP
    dzero (NZTOT-nZtot, mapx->t+nZtot, 1);
    dzero (NZTOT-nZtot, mapy->t+nZtot, 1);
#endif
  }

  // invFFT the maps to physical space
  rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
#ifdef TMAP
  rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
#endif
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*W->htot;
    for (i = 0; i < W->htot; i++) {
#ifdef TMAP // time dependent mesh
      U->base_h[j+i] -= dir * (mapx->t[plane] + mapx->z[plane] 
			       * W->base_h[j+i]);
      V->base_h[j+i] -= dir * (mapy->t[plane] + mapy->z[plane] 
			       * W->base_h[j+i]);				   
#else // only Z-deformation
      U->base_h[j+i] -= dir * mapx->z[plane] 
	                    * W->base_h[j+i];
      V->base_h[j+i] -= dir * mapy->z[plane] 
		            * W->base_h[j+i];
#endif
    }
  }

  // FFT the maps back to fourier space
  rfftw(MPlan, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(MPlan, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
#ifdef TMAP
  rfftw(MPlan, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(MPlan, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
#endif
  return;
}

#ifdef AXIAL
void Mapfield_w (Domain *omega, int dir)
{
  Element_List    *W = omega->W;
  Mapping       *mapz = omega->mapz;		
  register int k, j, i, plane;
  rfftw_plan MPlan     = rplan, 
             MPlan_inv = rplan_inv;


  if (abs(dir) == 2) dir /= 2; // handle the extra cases due to MapField

  int NZ = W->nz;
  int NZTOT = W->nztot;
  int nZtot = NZTOT; // keep this value as it may change because of dealiasing

  if (option("dealias")) {
    NZ    = 3*NZ/2;
    NZTOT = 3*NZTOT/2;
    MPlan     = rplan32;    
    MPlan_inv = rplan_inv32;

    // zero out for dealiasing
    dzero (NZTOT-nZtot, mapz->z+nZtot, 1);
    dzero (NZTOT-nZtot, mapz->t+nZtot, 1);
  }

  // invFFT the maps to physical space

  rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) mapz->z, 1, 0, 0, 0, 0);
  rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) mapz->t, 1, 0, 0, 0, 0);

  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*W->htot;
    for (i = 0; i < W->htot; i++) {

      W->base_h[j+i] -= dir * (mapz->t[plane] + mapz->z[plane] 
			       * W->base_h[j+i]);			   
    }
  }

  // FFT the maps back to fourier space
  rfftw(MPlan, 1, (FFTW_COMPLEX *) mapz->z, 1, 0, 0, 0, 0);
  rfftw(MPlan, 1, (FFTW_COMPLEX *) mapz->t, 1, 0, 0, 0, 0);
  return;
}
#endif
/* ------------------------------------------------------------------------ *
 * Compare() -- check answer against exact solution or zero if undefined    *
 * ------------------------------------------------------------------------ */

void Compare (Domain *omega, ACTION space)
{
  Element_List *U = omega->U,
               *V = omega->V,
               *W = omega->W,
               *P = omega->P;
#ifdef AXIAL
  int           nq = U->htot*U->nz;
  double *Wstore = dvector(0,nq);
  dcopy (nq, W->base_h, 1, Wstore, 1);	  
#endif   
  Nek_Trans_Type f_to_p = F_to_P,
                 p_to_f = P_to_F;
  int revert = 0, revertP = 0, NZ = U->nz;

#ifdef MAP 
  Mapping     *mapx = omega->mapx,
          *mapy = omega->mapy;

  if (option("dealias")) { // need dealiasing for correct Mapfield 
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
      NZ = 3*NZ/2;
  }
#endif

  // make sure velocity fields are in quadrature space
  if (U->fhead->state == 't') { // all 3 velocity fields should be same type
    U->Trans(U, J_to_Q);
    V->Trans(V, J_to_Q);
    W->Trans(W, J_to_Q);
    revert++;
  }

  if (P->fhead->state == 't') {
    P->Trans(P, J_to_Q);
    revertP++;
  }

  if (space == Fourier) {
    U->Trans(U, f_to_p);
    V->Trans(V, f_to_p);
    W->Trans(W, f_to_p);
    P->Trans(P, F_to_P); // P doesn't need dealiasing
  }

#ifdef MAP
  Mapfield (omega, -1);
#ifdef AXIAL
  Mapfield_w(omega, -1);
#endif
  if (option("dealias")) { // revert back to nz planes
    U->Trans(U, P_to_F32);
    U->Trans(U, F_to_P);
    V->Trans(V, P_to_F32);
    V->Trans(V, F_to_P);
    W->Trans(W, P_to_F32);
    W->Trans(W, F_to_P);
  }

  dparam_set("t", mapx->time);

  // invFFT the maps
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
#endif

#ifndef PLANAR_ERROR
#ifdef MAP
  GTerror(omega->U, omega->mapx, omega->mapy, omega->soln[0]);
  GTerror(omega->V, omega->mapx, omega->mapy, omega->soln[1]);
  GTerror(omega->W, omega->mapx, omega->mapy, omega->soln[2]);
  GTerror(omega->P, omega->mapx, omega->mapy, omega->soln[3]);
#else
  GTerror(omega->U, omega->soln[0]);
  GTerror(omega->V, omega->soln[1]);
  GTerror(omega->W, omega->soln[2]);
  GTerror(omega->P, omega->soln[3]);
#endif
#else
  for (int i = 0; i < U->nz; i++) {
    dparam_set("z", zmesh(i));
#ifdef MAP
    int iloc = parid(i);
    dparam_set("ORIGINX", mapx->d[iloc]);
    dparam_set("ORIGINY", mapy->d[iloc]);
#endif
    printf("Plane %d - ", iloc); U->flevels[i]->PTerror(omega->soln[0]); 
    printf("Plane %d - ", iloc); V->flevels[i]->PTerror(omega->soln[1]);
    printf("Plane %d - ", iloc); W->flevels[i]->PTerror(omega->soln[2]);
    printf("Plane %d - ", iloc); P->flevels[i]->PTerror(omega->soln[3]);
  }
#endif

#ifdef MAP
  // FFT the maps back
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);

  if (option("dealias")) { // Need to use 3/2 rule again
    U->Trans(U, P_to_F);
    U->Trans(U, F_to_P32);
    V->Trans(V, P_to_F);
    V->Trans(V, F_to_P32);
    W->Trans(W, P_to_F);
    W->Trans(W, F_to_P32);
  }
  //here w is changed into w'.

  Mapfield (omega, 1);
#ifdef AXIAL
  dcopy (nq, Wstore, 1, W->base_h, 1);	  
#endif	
#endif

  if (space == Fourier) {
    U->Trans(U, p_to_f);
    V->Trans(V, p_to_f);
    W->Trans(W, p_to_f);
    P->Trans(P, P_to_F); // P doesn't need dealiasing
  }

  if (revert) { // return velocity fields to modal storage if necessary
    U->Set_state('t');
    V->Set_state('t');
    W->Set_state('t');
    //    U->Trans(U, Q_to_J);
    //    V->Trans(V, Q_to_J);
    //    W->Trans(W, Q_to_J);
  }

  if (revertP) // return pressure to modal storage if necessary
    P->Set_state('t');
    //    P->Trans(P, Q_to_J);
  #ifdef AXIAL
  free(Wstore);
  #endif
  return;
}

/* ------------------------------------------------------------------------ *
 * 
 * TEST
 * ====
 * 
 * UpdateFFs()
 * 
 * ------------------------------------------------------------------------ */

#ifdef TEST

/* ------------------------------------------------------------------------ *
 * UpdateFFs() -- update time dependent forcing                             *
 * ------------------------------------------------------------------------ */

void UpdateFFs (Domain *omega)
{
#ifdef DOESNT_WORK
  const int     nz   = omega->U->nz,
                nel  = omega->U->nel;
  Element      *eU;
  Mapping          *mapx = omega->mapx,
               *mapy = omega->mapy;
  double       dummy;
  register int i, j, k, skip;

  if (omega->ForceFuncs) {
    Coord X; 
    X.x = dvector(0,QGmax*QGmax-1);
    X.y = dvector(0,QGmax*QGmax-1);
    
    const double dispx = mapx->d[0];
    const double dispy = mapy->d[0];
    const double time  = mapx->time;
    
    for (i = 0; i < nel; i++) {
      eU = omega->U->flist[i];
      int qtot = eU->qtot;
      eU->coord(&X);
      
      for (j = skip = 0; j < qtot; j++, skip += qtot) {
	drive_force (X.x[j]+dispx, X.y[j]+dispy, 0.0, time,
		     omega->ForceFuncs[0]+skip+j, omega->ForceFuncs[1]+skip+j,
		     &dummy);
      }
    }
  free_dvector (X.x,0);
  free_dvector (X.y,0);
  }

#if DIM == 3
  for (i = 0; i < ntot; i++)
    //#else
  for (k = 0; k < nz; k++)
    for (i = 0; i < ntot; i++)
      drive_force (x[i]+omega->mapx->d[parid(k)],
		   y[i]+omega->mapy->d[parid(k)],
		   z[parid(k)], omega->mapx->time,
		   (*FFx->base)+k*ntot+i, (*FFy->base)+k*ntot+i, (*FFz->base)+k*ntot+i);
  
  Transform (FFx, *FFx->base, Fourier);
  Transform (FFy, *FFy->base, Fourier);
  Transform (FFz, *FFz->base, Fourier);
#endif
#endif // if 0  
  return;
}

#endif

/* ------------------------------------------------------------------------ *
 * CountBCtype() - count bc of given type
 * ------------------------------------------------------------------------ */

static int CountBCtype (Bndry *Xbc, char type)
{
  int count = 0;
  for (Bndry *Vbc = Xbc; Vbc; Vbc = Vbc->next)
    if (Vbc->type == type) count++;
  return count;
}
