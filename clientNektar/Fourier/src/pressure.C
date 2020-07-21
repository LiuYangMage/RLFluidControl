/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: /users/tcew/Hybrid/Fourier/src/RCS/pressure_new.C,v $
 * $Revision: 1.1 $
 * $Date: 1996/07/01 15:26:38 $ 
 * $Author: tcew $ 
 * $State: Exp $ 
 *---------------------------------------------------------------------------*/
#include <mpi.h>
#include "nektarF.h"
// ce107
#include <rfftw.h>
extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;
#ifdef MAP
#include "map.h"
static void time_derivative_pbc (Bndry *B, double *w, 
				 double *ux, double *uy, 
				 double *vx, double *vy, 
				 double *wx, double *wy, 
				 Mapping *mapx, Mapping *mapy, int j);
static void Addaccel(Bndry *B, Domain *omega, int k);
#endif

/*
 * Build boundary conditions for the pressure
 */

Bndry *BuildPBCs(Element_List *P, Bndry *temp)
{
  int      nbcs = 0, Je = iparam("INTYPE");
  Bndry    *Pbc;
  Element  *E;
  register  int i,j;
  
  /* Count the number of boundaries to create */
  if(!temp) 
    return (Bndry*)0;

  while (temp[nbcs++].next);
  Pbc = (Bndry*)calloc(nbcs,sizeof(Bndry));

  for(i = 0; i < nbcs; ++i) {
    Pbc[i].id          = i;
    Pbc[i].type        = temp[i].type;
    Pbc[i].elmt        = P->flist[temp[i].elmt->id];
    Pbc[i].face        = temp[i].face;
    Pbc[i].next        = Pbc + i + 1;
  }
  Pbc[nbcs-1].next = (Bndry*) NULL;
  
  /* clear vertex solve mask */
  for(E = P->fhead; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      E->vert[i].solve = 1;
  
  /* Translate the boundary conditions */
  
  for(i = 0; i < nbcs; ++i) {
    switch (Pbc[i].type) {
    case 'O':
      Pbc[i].type = 'o';                       /* fall through */
      E = Pbc[i].elmt;
      
      for(j = 0; j < DIM; ++j)
	E->vert[E->vnum(Pbc[i].face,j)].solve = 0;
      
      E->MemBndry(Pbc+i,Pbc[i].face,(iparam("EQTYPE") == Rotational)? Je:1);
      
      break;     	  
      
    case 'W': // Set usrtype for wall pressure BC
      Pbc[i].usrtype = 'w';

    case 'V': case 'v':          /* allocate additional multi-step storage */
//zwang 20170108      
    case 'S': case 's':          /* allocate additional multi-step storage */
      Pbc[i].type = 'F';
      Pbc[i].elmt->Surface_geofac(Pbc+i);
      Pbc[i].elmt->MemBndry(Pbc+i,Pbc[i].face,Je);
      break;
    default:
      Pbc[i].elmt->MemBndry(Pbc+i,Pbc[i].face,1);
      break;
    }
  }
  
  return  bsort(Pbc,nbcs);
}

/* ------------------------------------------------------------------------- *
 * SetPBCs() -  Set boundary conditions for the pressure                     *
 *                                                                           *
 *                   Je                     1                                *
 *        dP/dn = SUM    beta[q] * [ N(u) - - curl ( curl u ) ] * n          *
 *                   q=0                    R                                *
 *                                                                           *
 * where n is the unit outward normal along the edge, u is the velocity      *
 * field, and N(u) are the non-linear terms in the momentum equation.  This  *
 * routine computes the RHS (already integrated) of the above equation at    *
 * the current time level and saves it in the corresponding Bedge for Press. *
 * ------------------------------------------------------------------------- */

static void CalcPbc (Domain *omega, Bndry *B, 
		     Element_List *V[3], Element_List *N[3], double nu);
static void IntPbc  (Domain *omega, Bndry *B, int Je);
static void outflow (Domain *omega, Bndry *B);

void SetPBCs(Domain *omega){
  Bndry    *Pbc = omega->Pbc[0];
  Element_List  *V[3], *N[3];
  int       Je   =  iparam("INTYPE");
  double    nu   =  dparam("KINVIS");

  V[0] =  omega->U;
  V[1] =  omega->V;
  V[2] =  omega->W;
  
  N[0] =  omega->Uf;
  N[1] =  omega->Vf;
  N[2] =  omega->Wf;
  
  /* Get the integration coefficients */  
  while (Pbc) {
    switch (Pbc->type) {
    case 'D': case 'N': case 'P':
      break;
      
    case 'o':
      if (iparam("EQTYPE") == Rotational){
	outflow (omega, Pbc);
	IntPbc  (omega, Pbc, Je);
      }
      break;
      
    case 'F': case 'R': {
      CalcPbc(omega, Pbc, V,N,nu);
      IntPbc (omega, Pbc, Je);
      break;
    }
    default:
      error_msg(SetPBCs -- unknown pressure b.c.)
	break;
    }
    
    Pbc = Pbc->next;
  }

  return;
}

static void IntPbc(Domain *omega, Bndry *B_base, int Je){
  register int i,j,k;
  const    int face = B_base->face;
  int      ll;
  double   beta[3];
  Bndry    *B;
  double   *tmp = dvector(0,LGmax-1);

  getbeta (beta);        
  
  for(k=0;k<omega->U->nz;++k){
    B = omega->Pbc[k]+B_base->id;
    /* integrate vertices */
    dzero(DIM,tmp,1);
    for(i = 0; i < DIM; ++i)
      for(j = 0; j < Je; ++j)
	tmp[i] += beta[j]*B->bvert[i+j*DIM];
    
    dcopy(DIM * --j, B->bvert, 1, B->bvert+DIM,1);
    dcopy(DIM      ,   tmp,      1, B->bvert,1);
    
    ll = B->elmt->edge[face].l;
    dzero(ll,tmp,1);
    for(j = 0; j < Je; ++j)
      daxpy(ll, beta[j], B->bedge[0] + j*ll, 1, tmp, 1);
    
    dcopy(ll * --j, B->bedge[0], 1, B->bedge[0] + ll, 1);
    dcopy(ll      , tmp        , 1, B->bedge[0]     , 1);
  }
  free(tmp);
}
 
#if 0
static void CalcPbc(Domain *omega, Bndry *B_base, Element_List *VL[3],
		    Element_List *NL[3], double nu)
{
  int     id    = B_base->elmt->id;
  int     i,j;

  Element *V[3], *N[3];
  Bndry   *B;
 
  int     qa = VL[0]->flist[id]->qa, tot = qa*VL[0]->flist[id]->qb;
  double  bet;  
  double  *Q    = dvector(0, tot*12-1);
  double  *QU   = Q+tot;
  double  *Uy   = QU+tot;
  double  *Vx   = Uy+tot;
  double  *QV   = Vx+tot;
  double  *Wz   = QV+tot;
  double  *Uzz  = Wz+tot;
  double  *Vzz  = Uzz+tot;
  double  *Wzx  = Vzz+tot;
  double  *Wzy  = Wzx+tot;
  double  *GxWx = Wzy+tot;
  double  *GxWy = GxWx+tot;

  for(j = 0; j < VL[0]->nz; ++j){
    B = omega->Pbc[j]+B_base->id;
    
    for(i =0; i < 3; ++i){
      V[i] = VL[i]->flevels[j]->flist[id];
      N[i] = NL[i]->flevels[j]->flist[id];
    }

    dzero(tot*12,Q,1);
    
  //ce107 - do this here instead of convective.C & analyzer.C
 #if defined (MAP)
    if (B->usrtype != 'w') { // do this for all modes, non-wall boundaries

      if (option("farfield")) { // U,V=/=constant in the farfield
	double *workspace = dvector(0,11*tot-1),
	  *Ux = workspace,
	  *Vy = Ux + tot,
	  *w  = Vy + tot,
	  *wx = w  + tot,
	  *wy = wx + tot,
	  *Wx = wy + tot,
	  *Wy = Wx + tot,
	  *ux = Wy + tot,
	  *uy = ux + tot,
	  *vx = uy + tot,
	  *vy = vx + tot;

	V[2]->GetFace(*V[2]->h, B->face, w);

	V[0]->Grad_d( Ux, Uy, 0, 'a');
	V[1]->Grad_d( Vx, Vy, 0, 'a');
	 V[2]->Grad_d( Wx, Wy, 0, 'a');

	// this needs to be done here before the GetFace()...
	/*******   Q = Vx-Uy   ********/
	dvsub  (tot, Vx,  1, Uy,  1, Q, 1);      
	
	V[0]->GetFace(Ux, B->face, Ux);
	V[0]->GetFace(Uy, B->face, Uy);
	V[1]->GetFace(Vx, B->face, Vx);
	V[1]->GetFace(Vy, B->face, Vy);
	V[2]->GetFace(Wx, B->face, Wx);
	V[2]->GetFace(Wy, B->face, Wy);
	
	V[0]->InterpToFace1(B->face,Ux,ux);
	V[0]->InterpToFace1(B->face,Uy,uy);
	V[1]->InterpToFace1(B->face,Vx,vx);
	V[1]->InterpToFace1(B->face,Vy,vy);
	V[2]->InterpToFace1(B->face,Wx,wx);
	V[2]->InterpToFace1(B->face,Wy,wy);
	
	time_derivative_pbc (B, w, ux, uy, vx, vy, wx, wy,
			     omega->mapx, omega->mapy, parid(j));
	free (workspace);
      } else { // U,V=constant in the farfield, W can vary with z
	// NOTE: Needs fixing
	// Get W on the face
	// V[2]->GetFace(*V[2]->h, face, w);
 
	time_derivative_pbc (B, (double *)NULL, 
			     (double *)NULL, (double *)NULL, 
			     (double *)NULL, (double *)NULL, 
			     (double *)NULL, (double *)NULL, 
			     omega->mapx, omega->mapy, parid(j));
	/***** quick note ... (Wx-Uz)z = (Wx)z - Uzz = (Wz)x-Uzz ********/
	V[0]->Grad_d( 0, Uy, 0, 'y');
	V[1]->Grad_d( Vx, 0, 0, 'x');
	
	/*******   Q = Vx-Uy   ********/
	dvsub  (tot, Vx,  1, Uy,  1, Q, 1);      
      }
    } else
#endif // ifdef MAP
      { // wall boundaries
	/***** quick note ... (Wx-Uz)z = (Wx)z - Uzz = (Wz)x-Uzz ********/
	V[0]->Grad_d( 0, Uy, 0, 'y');
	V[1]->Grad_d( Vx, 0, 0, 'x');
	
	/*******   Q = Vx-Uy   ********/
	 dvsub  (tot, Vx,  1, Uy,  1, Q, 1);
      }

    /*******   GxWy = (Vx-Uy)x  GxWx = (Vx-Uy)y ******/
    V[0]->Grad_h(Q,  GxWy, GxWx, 0, 'a');
    
    if(parid(j) < 2)
      dzero(tot, Wz, 1);
    else{
      if(parid(j)%2)
	dsmul(tot,Beta(j), VL[2]->flevels[j-1]->flist[id]->h[0], 1, Wz, 1);
      else
	dsmul(tot,-Beta(j),VL[2]->flevels[j+1]->flist[id]->h[0], 1, Wz, 1);
    }
    
    V[2]->Grad_h(Wz,  Wzx,  Wzy, 0, 'a'); 
    
    bet = -Beta(j)*Beta(j);
    dsmul(tot, bet, *V[0]->h,1, Uzz,1);
    dsmul(tot, bet, *V[1]->h,1, Vzz,1);
    
    /*******   QU= Wzx-Uzz   ********/
    dvsub  (tot, Wzx, 1, Uzz, 1, QU, 1);
    /*******   QV= Wzy-Vzz   ********/
    dvsub  (tot, Wzy, 1, Vzz, 1, QV, 1);
    /*******  GxWx = GxWx+QU = (Vx-Uy)y + (Wx-Uz)z  *******/
    dvadd (tot, GxWx, 1, QU, 1, GxWx, 1);
     /*******  GxWy = QV-GxWy = (Wy-Vz)z - (Vx-Uy)x   *******/
    dvsub (tot, QV, 1, GxWy, 1, GxWy, 1);  
    
    //    nu *=-1.0;
    dsvtvp (tot, -nu, GxWx, 1, *N[0]->h, 1, GxWx, 1);
    dsvtvp (tot, -nu, GxWy, 1, *N[1]->h, 1, GxWy, 1);
    
    V[0]->GetFace(GxWx, B->face, GxWx);
    V[0]->GetFace(GxWy, B->face, GxWy);
    
    V[0]->InterpToFace1(B->face,GxWx,Q);
    V[0]->InterpToFace1(B->face,GxWy,GxWx);
    
    if(V[0]->curvX){
      dvmul  (qa,B->nx.p,1,Q   ,1,Q,1);
      dvvtvp (qa,B->ny.p,1,GxWx,1,Q,1,Q,1);
    }
    else {
      dsmul  (qa,B->nx.d,Q   ,1,Q,1);
      dsvtvp (qa,B->ny.d,GxWx,1,Q,1,Q,1);
    }
 #ifdef MAP
    // ce107
    // store old values
    if (B->usrtype != 'w') {
      Uy[0] = B->bvert[0];
      Uy[1] = B->bvert[1];
      dcopy(B->elmt->edge[B->face].l,B->bedge[0],1,Uy+2,1);
    }
#endif
    B->elmt->MakeFlux(B,0,Q);
#ifdef MAP
    // add stored values
    if (B->usrtype != 'w') {
      B->bvert[0] += Uy[0];
      B->bvert[1] += Uy[1];
      //    dvadd(B->elmt->edge[B->face].l,Uy+2,1,B->bedge[0],1,B->bedge[0],1);
      daxpy(B->elmt->edge[B->face].l,1.0,Uy+2,1,B->bedge[0],1);
    }
#endif
    //    B->bvert[0] = 0.0;
    //    B->bvert[1] = 0.0;
    //    dzero(B->elmt->edge[B->face].l,B->bedge[0],1);
  }
  free(Q);
}


#else
 
static void CalcPbc(Domain *omega, Bndry *B_base, Element_List *VL[3],
		    Element_List *NL[3], double nu)
{
  int     id    = B_base->elmt->id;
  int     i,j;

  Element *V[3], *N[3];
  Bndry   *B;

  int     qa = VL[0]->flist[id]->qa, tot = qa*VL[0]->flist[id]->qb;
  double  bet;  
  double  *Q    = dvector(0, tot*12-1);
  double  *QU   = Q+tot;
  double  *Uy   = QU+tot;
  double  *Vx   = Uy+tot;
  double  *QV   = Vx+tot;
  double  *Wz   = QV+tot;
  double  *Uzz  = Wz+tot;
  double  *Vzz  = Uzz+tot;
  double  *Wzx  = Vzz+tot;
  double  *Wzy  = Wzx+tot;
  double  *GxWx = Wzy+tot;
  double  *GxWy = GxWx+tot;

  for(j = 0; j < VL[0]->nz; ++j){
    B = omega->Pbc[j]+B_base->id;
    
    for(i =0; i < 3; ++i){
      V[i] = VL[i]->flevels[j]->flist[id];
      N[i] = NL[i]->flevels[j]->flist[id];
    }

    dzero(tot*12,Q,1);
    
    /***** quick note ... (Wx-Uz)z = (Wx)z - Uzz = (Wz)x-Uzz ********/
    V[0]->Grad_d( 0, Uy, 0, 'y');
    V[1]->Grad_d( Vx, 0, 0, 'x');
    
    /*******   Q = Vx-Uy   ********/
    dvsub  (tot, Vx,  1, Uy,  1, Q, 1);
    /*******   GxWy = (Vx-Uy)x  GxWx = (Vx-Uy)y ******/
    V[0]->Grad_h(Q,  GxWy, GxWx, 0, 'a');
    
#ifdef FLOK
    if(VL[0]->nz == 2){
      if(j%2)
	dsmul(tot, Beta(j), VL[2]->flevels[j-1]->flist[id]->h[0], 1, Wz, 1);
      else
	dsmul(tot,-Beta(j), VL[2]->flevels[j+1]->flist[id]->h[0], 1, Wz, 1);
    }
    else
      dsmul(tot,-Beta(j), VL[2]->flevels[j]->flist[id]->h[0], 1, Wz, 1);
#else
    if(parid(j) < 2)
      dzero(tot, Wz, 1);
    else{
      if(parid(j)%2)
	dsmul(tot,Beta(j), VL[2]->flevels[j-1]->flist[id]->h[0], 1, Wz, 1);
      else
	dsmul(tot,-Beta(j),VL[2]->flevels[j+1]->flist[id]->h[0], 1, Wz, 1);
    }
#endif
    
    V[2]->Grad_h(Wz,  Wzx,  Wzy, 0, 'a'); 
    
    bet = -Beta(j)*Beta(j);
    dsmul(tot, bet, *V[0]->h,1, Uzz,1);
    dsmul(tot, bet, *V[1]->h,1, Vzz,1);
    
    /*******   QU= Wzx-Uzz   ********/
    dvsub  (tot, Wzx, 1, Uzz, 1, QU, 1);
    /*******   QV= Wzy-Vzz   ********/
    dvsub  (tot, Wzy, 1, Vzz, 1, QV, 1);
    /*******  GxWx = GxWx+QU = (Vx-Uy)y + (Wx-Uz)z  *******/
    dvadd (tot, GxWx, 1, QU, 1, GxWx, 1);
    /*******  GxWy = QV-GxWy = (Wy-Vz)z - (Vx-Uy)x  *******/
    dvsub (tot, QV, 1, GxWy, 1, GxWy, 1);  
    
    dsvtvp (tot, -nu, GxWx, 1, *N[0]->h, 1, GxWx, 1);
    dsvtvp (tot, -nu, GxWy, 1, *N[1]->h, 1, GxWy, 1);
    
    V[0]->GetFace(GxWx, B->face, GxWx);
    V[0]->GetFace(GxWy, B->face, GxWy);
    
    V[0]->InterpToFace1(B->face,GxWx,Q);
    V[0]->InterpToFace1(B->face,GxWy,GxWx);
    
    if(V[0]->curvX){
      dvmul  (qa,B->nx.p,1,Q   ,1,Q,1);
      dvvtvp (qa,B->ny.p,1,GxWx,1,Q,1,Q,1);
    }
    else {
      dsmul  (qa,B->nx.d,Q   ,1,Q,1);
      dsvtvp (qa,B->ny.d,GxWx,1,Q,1,Q,1);
    }
    B->elmt->MakeFlux(B,0,Q);
  }
  free(Q);
}

#endif
/*
 * PI = p + 1/2 U.U   (outflow boundary conditions)
 */

static void outflow (Domain *omega, Bndry *B){

  // Have to fix this to make sure that u.u is done in physical

  const    int id = B->elmt->id, face = B->face;
  Element  *U,*V,*W;
  int       q, vid;
  int       Bid = B->id;
  int       i,j,m;
#ifdef MAP /* for the flexible code this needs to change - keep in mind 
              that for dealising we need 3*nz/2 values for mapvel. */
  double mapvel[3];
  mapvel[0] = omega->mapx->t[0];
  mapvel[1] = omega->mapy->t[0];
  mapvel[2] = 0.;
#endif
  U = omega->U->flevels[0]->flist[id];
  V = omega->V->flevels[0]->flist[id];
  W = omega->W->flevels[0]->flist[id];

  if(U->identify() == Nek_Tri)
    if(face == 0){
      q = U->qa;
      vid = 1;
    }
    else{
      q = U->qb;
      vid = 2;
    }    
  else
    switch(face){
    case 0:
      q = U->qa;
      vid = 1;
      break;
    case 1:
      q = U->qb;
      vid = 2;
      break;
    case 2:
      q = U->qa;
      vid = 2;
    case 3:
      q = U->qb;
      vid = 3;
      break;
    }

  int np = q+1;
  int nz = omega->U->nz; 
  int NZ = omega->U->nztot; 
  int nZ = NZ; // keep this value as NZ may change because of dealiasing
  int procid = option("PROCID");

  // set the RFFT structs
  rfftw_plan Plan     = rplan, 
             Plan_inv = rplan_inv;
  if (option("dealias")) {
    NZ = 3*NZ/2;
    Plan     = rplan32;    
    Plan_inv = rplan_inv32;
  }

  double *tmp  = dvector(0, 3*np*nz-1);
  double *uout = dvector(0, 3*np*nZ-1); // here we need nZ
  double **u   = dmatrix(0, np-1, 0, NZ-1);
  double *fvec = dvector(0, NZ-1);

  dzero(np*NZ, u[0], 1);

  for(i = 0; i < nz; ++i){
    U = omega->U->flevels[i]->flist[id];
    V = omega->V->flevels[i]->flist[id];
    W = omega->W->flevels[i]->flist[id];
    
    tmp[    3*i*np] = U->vert[vid].hj[0];
    tmp[(3*i+1)*np] = V->vert[vid].hj[0];
    tmp[(3*i+2)*np] = W->vert[vid].hj[0];
    
    U->GetFace(*U->h, face, tmp+1+     3*i*np);
    V->GetFace(*V->h, face, tmp+1+ (3*i+1)*np);
    W->GetFace(*W->h, face, tmp+1+ (3*i+2)*np);
  }
  
  MPI_Allgather (tmp, 3*nz*np, MPI_DOUBLE, uout, 3*nz*np, MPI_DOUBLE, 
		 MPI_COMM_WORLD);
  
  for (i = 0; i < np; ++i) {
    for (j = 0; j < 3; ++j) {
      dcopy  (nZ, uout + np*j + i, 3*np, fvec, 1);
      if (option("dealias")) dzero (NZ-nZ, fvec+nZ, 1); // zero out...
      rfftw(Plan_inv, 1, (FFTW_COMPLEX *) fvec, 1, 0, 0, 0, 0);
#ifdef MAP
      for (m = 0; m < NZ; ++m)
	u[i][m] += (fvec[m]+mapvel[j])*(fvec[m] + mapvel[j]);
#else
      for (m = 0; m < NZ; ++m)
	u[i][m] += fvec[m]*fvec[m];
#endif
    }
    dscal (NZ, .5, u[i], 1);
    rfftw(Plan, 1, (FFTW_COMPLEX *) u[i], 1, 0, 0, 0, 0);
  }


  for(i = 0; i < nz; ++i){
    B = omega->Pbc[i]+Bid;
    B->bvert[0] = u[0][i+procid*nz];
    dcopy(q, u[1]+i+procid*nz, NZ, tmp, 1);
    B->elmt->JtransEdge(B,face,0,tmp);
  }

  free(tmp); free(uout); free_dmatrix(u, 0, 0); free(fvec);

  return;
}

#ifdef MAP
/*
 * Add mapping acceleration term to the pressure BCs for the case of 
 * implicit treatment of the extra mapping terms. Only needs to be done
 * for the wall as the cancellation doesn't happen there.
 */
void AddAccelPBCs(Domain *omega)
{
  register int k;
  const int NZ = omega->U->nz;
  Bndry *Pbc;

  for (k = 0; k < NZ; k++) {
    Pbc = omega->Pbc[k];
    while (Pbc) {
      if (Pbc->type == 'F') // do this for all Newman BCs.
	Addaccel(Pbc, omega, parid(k));
      Pbc = Pbc->next;
    }
  }
  
  return;
}

/*
 * Add mapping acceleration term to the pressure BCs for the case of 
 * implicit treatment of the extra mapping terms. Only correct for 2D/Fourier.
 */

static void Addaccel(Bndry *B, Domain *omega, int k){
    const int    qa = B->elmt->qa;
    double       *Q = dvector(0,qa-1),
                 *store = dvector(0,LGmax-1),
                 Ax, Ay,
                 nu = dparam("KINVIS");

    if (B->usrtype == 'w') { // need to have the full "acceleration" term
      Ax = omega->mapx->tzz[k] * nu - omega->mapx->tt[k];
      Ay = omega->mapy->tzz[k] * nu - omega->mapy->tt[k];
    } else { // the other term cancels out
      Ax = omega->mapx->tzz[k] * nu;
      Ay = omega->mapy->tzz[k] * nu;
    }

    if(B->elmt->curvX) {
      dsmul (qa,Ax,B->nx.p,1,Q,1);
      daxpy (qa,Ay,B->ny.p,1,Q,1);
    }
    else {
      dfill (qa,Ax*B->nx.d,Q,1);
      dsadd (qa,Ay*B->ny.d,Q,1,Q,1);
    }

    /* store old values */
    store[0] = B->bvert[0];
    store[1] = B->bvert[1];
    dcopy(B->elmt->edge[B->face].l,B->bedge[0],1,store+2,1);
    
    B->elmt->MakeFlux(B,0,Q);
    
    /* add stored values */
    B->bvert[0] += store[0];
    B->bvert[1] += store[1];
//    dvadd(B->elmt->edge[B->face].l,store+2,1,B->bedge[0],1,B->bedge[0],1);
    daxpy(B->elmt->edge[B->face].l,1.0,store+2,1,B->bedge[0],1);
    free(Q); free(store);
    return;
}

/* ------------------------------------------------------------------------- *
 * time_derivative_pbc () -- Compute the du/dt term for the Pressure BC's    *
 *                                                                           *
 * This function computes the -n.dU/dt term for the pressure BC's. In the    *
 * transformed frame:                                                        *
 *                                                                           *
 * du/dt = d(u' - mapx_t - w' mapx_z)/dt' + mapx_t du/dx + mapy_t du/dy      *
 * dv/dt = d(v' - mapy_t - w' mapy_z)/dt' + mapx_t dv/dx + mapy_t du/dy      *
 * dw/dt = dw'/dt' + mapx_t dw/dx + mapy_t dw/dy                             *
 * that is:                                                                  *
 * du/dt = du'/dt' - mapx_tt - dw'/dt' mapx_z - w' mapx_zt +                 *
 *         mapx_t du/dx + mapy_t du/dy                                       *
 * dv/dt = dv'/dt' - mapy_tt - dw'/dt' mapy_z - w' mapy_zt +                 *
 *         mapx_t dv/dx + mapy_t dv/dy                                       *
 * dw/dt = dw'/dt' + mapx_t dw/dx + mapy_t dw/dy                             *
 *                                                                           *
 * ------------------------------------------------------------------------- */

static void time_derivative_pbc (Bndry *B, double *w, 
				 double *ux, double *uy, 
				 double *vx, double *vy, 
				 double *wx, double *wy, 
				 Mapping *mapx, Mapping *mapy, int j)
{
  double   mvx, mvy, max, may;
  const int qa = B->elmt->qa;
  double    *Q = dvector(0,qa-1);
  // farfield = 0, acc_impl = 0 -> Pbc_case = 0 (default)
  // farfield = 0, acc_impl = 1 -> Pbc_case = 1 
  // farfield = 1, acc_impl = 0 -> Pbc_case = 2 
  // farfield = 1, acc_impl = 1 -> Pbc_case = 3 
  int Pbc_case = 2*option("farfield") + option("acc_impl");

  switch (Pbc_case) {
  case 0:
#if defined (TMAP) && !defined(DAVE)
    // U,V,W=constant in the farfield, ux = uy = vx = vy = wx = wy = 0,
    // the default case for cylinder calculations
    // explicit treatment - add (mapx_tt,mapy_tt)
    max = mapx->tt[j];
    may = mapy->tt[j];
    // Assuming that du'/dt' = dv'/dt' = w' = dw'/dt' = 0.0 at the boundary:
    if(B->elmt->curvX) {
      dsmul (qa,max,B->nx.p,1,Q,1);
      daxpy (qa,may,B->ny.p,1,Q,1);
    }
    else {
      dfill (qa,max*B->nx.d,Q,1);
      dsadd (qa,may*B->ny.d,Q,1,Q,1);
    }
    B->elmt->MakeFlux(B,0,Q);
#else // otherwise zero out everything
    B->bvert[0] = 0.0;
    B->bvert[1] = 0.0;
    dzero(B->elmt->edge[B->face].l,B->bedge[0],1);
#endif
    break;
  case 1:
    // U,V=constant in the farfield, ux = uy = vx = vy = 0,
    // implicit treatment of the mapping terms
    // no need to add (mapx_tt,mapy_tt) 
    B->bvert[0] = 0.0;
    B->bvert[1] = 0.0;
    dzero(B->elmt->edge[B->face].l,B->bedge[0],1);
    break;
#ifdef DOESNT_WORK
  case 2:
    // U,V=/=constant in the farfield 
    mvx = mapx->t[j];
    mvy = mapy->t[j];
    max = mapx->tt[j];
    may = mapy->tt[j];
    // Assuming that du'/dt' = dv'/dt' = 0.0 at the boundary:
    // explicit treatment of the mapping terms - add (mapx_tt,mapy_tt)
    dscal (qa, -mvx, ux, 1);
    daxpy (qa, -mvy, uy, 1, ux, 1);
    dsadd (qa, max, ux, 1, ux, 1);
    dscal (qa, -mvx, vx, 1);
    daxpy (qa, -mvy, vy, 1, vx, 1);
    dsadd (qa, may, vx, 1, vx, 1);
    if(B->elmt->curvX) {
      dvmul (qa,B->nx.p,1,ux,1,Q,1);
      dvvtvp (qa,B->ny.p,1,vx,1,Q,1,Q,1);
    }
    else {
      dsmul (qa,B->nx.d,ux,1,Q,1);
      dsvtvp (qa,B->ny.d,vx,1,Q,1,Q,1);
    }
    B->elmt->MakeFlux(B,0,Q);
    break;
  case 3:
    // U,V=/=constant in the farfield 
    mvx = mapx->t[j];
    mvy = mapy->t[j];
    // Assuming that du'/dt' = dv'/dt' = 0.0 at the boundary:
    // implicit treatment of the mapping terms
    // no need to add the (mapx_tt,mapy_tt) as they cancel out */
    dscal (qa, -mvx, ux, 1);
    daxpy (qa, -mvy, uy, 1, ux, 1);
    dscal (qa, -mvx, vx, 1);
    daxpy (qa, -mvy, vy, 1, vx, 1);
    if(B->elmt->curvX) {
      dvmul (qa,B->nx.p,1,ux,1,Q,1);
      dvvtvp (qa,B->ny.p,1,vx,1,Q,1,Q,1);
    }
    else {
      dsmul (qa,B->nx.d,ux,1,Q,1);
      dsvtvp (qa,B->ny.d,vx,1,Q,1,Q,1);
    }
    B->elmt->MakeFlux(B,0,Q);
    break;
#endif
  default:
    error_msg(time_derivative_pbc -- invalid case.)  
    break;
  }
  free(Q); 
  return;
}

#endif /* ifdef MAP */


#ifdef SPM
void SPM_SetPBCs(Domain *omega);

static void SPM_CalcPbc (Domain *omega, Bndry *B, 
		                     Element_List *V[3]);

void SPM_SetPBCs(Domain *omega){
  Bndry    *Pbc = omega->Pbc[0];
  Element_List  *V[3], *N[3];
  int       Je   =  iparam("INTYPE");
  double    nu   =  dparam("KINVIS");

  V[0] =  omega->Uf;
  V[1] =  omega->Vf;
  V[2] =  omega->Wf;
  
  /* Get the integration coefficients */  
  while (Pbc) {
    switch (Pbc->type) {
    case 'D': case 'N': case 'P':
      break;
      
    case 'o':
      if (iparam("EQTYPE") == Rotational){
  fprintf(stderr,"SPM_SetPBcs not implemeted for rotional case yet ! \n");
  exit(-1);
//	outflow (omega, Pbc);
//	IntPbc  (omega, Pbc, Je);
      }
      break;
      
    case 'F': case 'R': {
      SPM_CalcPbc(omega, Pbc, V);
      IntPbc (omega, Pbc, Je);
      break;
    }
    default:
      error_msg(SetPBCs -- unknown pressure b.c.)
	break;
    }
    
    Pbc = Pbc->next;
  }

  return;
}

  
static void SPM_CalcPbc(Domain *omega, Bndry *B_base, Element_List *VL[3])
{
  int     id    = B_base->elmt->id;
  int     i,j;

  Element *V[3];
  Bndry   *B;

  int     qa = VL[0]->flist[id]->qa, tot = qa*VL[0]->flist[id]->qb;
  double  bet;  
  double  *Q    = dvector(0, tot*5-1);
  double  *Qx   = Q+tot;
  double  *Qy   = Qx+tot;
  double  *uf   = Qy+tot;
  double  *vf   = uf+tot;

  for(j = 0; j < VL[0]->nz; ++j){
    B = omega->Pbc[j]+B_base->id;
    
    for(i =0; i < 3; ++i){
      V[i] = VL[i]->flevels[j]->flist[id];
    }

    dzero(tot*5,Q,1);
    
    V[0]->GetFace(*V[0]->h, B->face, uf);
    V[1]->GetFace(*V[1]->h, B->face, vf);
    
    V[0]->InterpToFace1(B->face,uf,Qx);
    V[1]->InterpToFace1(B->face,vf,Qy);
    
    if(V[0]->curvX){
      dvmul  (qa,B->nx.p,1,Qx   ,1,Q,1);
      dvvtvp (qa,B->ny.p,1,Qy,1,Q,1,Q,1);
    }
    else {
      dsmul  (qa,B->nx.d,Qx   ,1,Q,1);
      dsvtvp (qa,B->ny.d,Qy,1,Q,1,Q,1);
    }

//    dscal(tot, 1./dparam("DELT"), Q, 1);
    dscal(tot, getgamma()/dparam("DELT"), Q, 1);

    B->elmt->MakeFlux(B,0,Q);
  }
  free(Q);
}

#endif
