/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $RCSfile:
 * $Revision:
 * $Author: 
 * $Date: 
 * $State:
 * ------------------------------------------------------------------------- */
#ifdef FORCES
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <veclib.h>
#include "nektarF.h"
#ifdef MAP
#include "map.h"
#endif

#include <rfftw.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) > (b) ? (b) : (a))

extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;

static int init=1; /* externals */

void set_SurGeofac(Bndry *Pbc, Bndry *Ubc);

/*-------------------------------------------------------------------*
 * This is a function to calculate the forces of the fluid acting on *
 * a body pressumed to be represented by the flag 'W'.               *
 *                                                                   *
 *   F_i = P n_i - RHO*KINVIS*T_ij n_j  (RHO=1)                     *
 *                                                                   *
 * Note: This works out force for previous time step                 *
 *-------------------------------------------------------------------*/

void forces(Domain *omega, double time){
  static int step = 0, hisstep = option("hisstep");

  register int i,j,k;
  int      nz = option("NZ"), nztot = option("NZTOT");
  int      nprocs = option("NPROCS");
  int      eid,qa,face,lmax=0, trip = 0;
  Bndry   *B = omega->Ubc[0];
  double  *workspace, *fpk[3], *fvk[3], *send, *recv, *w;
  double   kinvis = dparam("KINVIS");
  double  *wk = dvector(0,max(LGmax+1,QGmax)*QGmax-1);
  double  **D = dmatrix(0,6,0,QGmax*QGmax-1);
  double  *ux = D[0], *uy = D[1];
  double  *vx = D[2], *vy = D[3];
  double  *wx = D[4], *wy = D[5], *p = D[6];
  Element_List *U;
  Element *eU, *eV, *eW, *eP;
  FILE    *fout = omega->fce_file;
  Basis   *b;
    
  workspace = dvector(0,12*nztot-1);
  dzero(12*nztot,workspace,1);

  send   = workspace;
  recv   = workspace+6*nztot-1;

  fpk[0] = workspace;
  fpk[1] = workspace +   nz;
  fpk[2] = workspace + 2*nz;
  fvk[0] = workspace + 3*nz;
  fvk[1] = workspace + 4*nz;
  fvk[2] = workspace + 5*nz;

  /* write out header */
  if(init){
    ROOT{
      fprintf(fout,"# Force acting on body\n");
      fprintf(fout,"# \n");
      fprintf(fout,"# Time  (Fx-press, Fx-visc)Fx "
	      " (Fy-press, Fy-visc) Fy (Fz-press, Fz-visc) Fz\n");
    }
    /* need to set up surface geometric factors in U from P */
    for(k = 0; k < nz; ++k)
      set_SurGeofac(omega->Pbc[k],omega->Ubc[k]);
      
    init = 0;
  }
    
  for(;B; B = B->next)
//    if(B->type == 'W'){
    if(B->usrtype == 'W'){
      trip = 1;
      U    = omega->U;
      face = B->face;
      eid  = B->elmt->id;
      qa   = U->flist[eid]->qa;
	
      if(U->flist[eid]->lmax != lmax){
	b    = U->flist[eid]->getbasis();
	lmax = U->flist[eid]->lmax;
	getzw(qa,&w,&w,'a');
      }
	
      for(k = 0; k < nz; ++k){
	eU = omega->U->flevels[k]->flist[eid];
	eV = omega->V->flevels[k]->flist[eid];
	eW = omega->W->flevels[k]->flist[eid];
	eP = omega->P->flevels[k]->flist[eid];
	  
	eP->Jbwd(eP,b); eP->state = 't';
	eU->Jbwd(eU,b); eU->state = 't';
	eV->Jbwd(eV,b); eV->state = 't';
	eW->Jbwd(eW,b); eW->state = 't';
	  
	eU->Grad_d(ux, uy, 0, 'a');
	eV->Grad_d(vx, vy, 0, 'a');
	eW->Grad_d(wx, wy, 0, 'a');
	  
	/* get appropriate face values from grad(V) */
	for(i = 0; i < 6; ++i){
	  eU->GetFace(D[i] ,face,wk);
	  eU->InterpToFace1(face,wk,D[i]);
	}

	eP->GetFace(eP->h[0],face,wk);
	eP->InterpToFace1(face,wk,p);

	// calculate modal forces - nz,uz,vz and wz are all zero in four. exp

	if(eU->curvX)
	  for(i = 0; i < qa; ++i){
	    fpk[0][k] += p[i]*B->nx.p[i]*B->sjac.p[i]*w[i];
	    fpk[1][k] += p[i]*B->ny.p[i]*B->sjac.p[i]*w[i];
	    fpk[2][k] += 0.0;

	    fvk[0][k] -= (2.0*ux[i]*B->nx.p[i]+(uy[i]+vx[i])*
			  B->ny.p[i])*B->sjac.p[i]*w[i];
	    fvk[1][k] -= (2.0*vy[i]*B->ny.p[i]+(uy[i]+vx[i])*
			  B->nx.p[i])*B->sjac.p[i]*w[i];
	    fvk[2][k] -= (wx[i]*B->nx.p[i] + wy[i]*
			  B->ny.p[i])*B->sjac.p[i]*w[i];
	  }
	else
	  for(i = 0; i < qa; ++i){
	    fpk[0][k] += p[i]*B->nx.d*B->sjac.d*w[i];
	    fpk[1][k] += p[i]*B->ny.d*B->sjac.d*w[i];
	    fpk[2][k] += 0.0;

	    fvk[0][k] -= (2.0*ux[i]*B->nx.d+(uy[i] + vx[i])*
			  B->ny.d)*B->sjac.d*w[i];
	    fvk[1][k] -= (2.0*vy[i]*B->ny.d+(uy[i] + vx[i])*
			  B->nx.d)*B->sjac.d*w[i];
	    fvk[2][k] -= (wx[i]*B->nx.d + wy[i]*
			  B->ny.d)*B->sjac.d*w[i];
	  }
      }    
    }

  if(trip){
    MPI_Allgather (send, 6*nz, MPI_DOUBLE, recv, 6*nz, MPI_DOUBLE, 
		   MPI_COMM_WORLD);

    fpk[0] = workspace;
    fpk[1] = workspace +   nztot;
    fpk[2] = workspace + 2*nztot;
    fvk[0] = workspace + 3*nztot;
    fvk[1] = workspace + 4*nztot;
    fvk[2] = workspace + 5*nztot;

    for(i = 0; i < nprocs; ++i)
      for(k = 0; k < 3; ++k) {
#ifdef DDC
	dcopy(nz,recv+i*6*nz + k*nz       ,1,fpk[k]+i*nz,1);
	dcopy(nz,recv+i*6*nz + 3*nz + k*nz,1,fvk[k]+i*nz,1);
#else
	for (j = 0; j < nz; j++)
	  fpk[k][i*nz+j] = recv[(6*i+k)*nz +j];
	for (j = 0; j < nz; j++)
	  fvk[k][i*nz+j] = recv[(6*i+3+k)*nz +j];
#endif
      }
     
//  double ffx = 0.;
//  double ffy = 0.;
//  double ffz = 0.;

double Cf,Cd;
    ROOT {
      if (!(step%hisstep)) 
	fprintf(fout,"%lf ( %lf %lf ) %lf ( %lf %lf ) %lf ( %lf %lf ) %lf \n",
		time,
		fpk[0][0],kinvis*fvk[0][0],fpk[0][0]+kinvis*fvk[0][0],
		fpk[1][0],kinvis*fvk[1][0],fpk[1][0]+kinvis*fvk[1][0],
		fpk[2][0],kinvis*fvk[2][0],fpk[2][0]+kinvis*fvk[2][0]);
      step++;

//      if(dparam("PIPE"))
//        ffz = 4.*(fpk[2][0]+kinvis*fvk[2][0])/4./atan(1.);
//     MPI_Bcast(&ffz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    Cf = fpk[1][0]+kinvis*fvk[1][0];
    Cd = fpk[0][0]+kinvis*fvk[0][0];
    }

    MPI_Bcast(&Cf,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Cd,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier  (MPI_COMM_WORLD);                  /* sync before work */
//   omega->fz = ffz;
    dparam_set("D_DRAGCOEFF", 2.0*Cd);
    dparam_set("D_LIFTCOEFF", 2.0*Cf);

#ifdef MAP // we are in Fourier space
        
#ifdef SPM

 #ifdef OLDFFTS
    realft (nztot/2, omega->mapx->f, 1);
    realft (nztot/2, omega->mapy->f, 1);
 #else
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) omega->mapx->f, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) omega->mapy->f, 1, 0, 0, 0, 0);
 #endif

 //  ROOTONLY
//   for (i = 0; i < nztot; i += 1)
//     fprintf(stderr,"m = %d  fx = %g   fy = %g \n",i,omega->mapx->f[i], omega->mapy->f[i]);

   if(strcmp(omega->vSPM[0].shape[0].type,"wall") != 0) //for case 'wall', force from forces() is used
    for (i = 0; i < nztot; i++) {
     int offset = i*10;
     omega->mapx->f[i] /= omega->vSPM[0].shape[0].data[offset+7];
     omega->mapy->f[i] /= omega->vSPM[0].shape[0].data[offset+7];
    }

  #ifdef OLDFFTS
    realft (nztot/2, omega->mapx->f, -1);
    realft (nztot/2, omega->mapy->f, -1);
#else
    rfftw(rplan, 1, (FFTW_COMPLEX *) omega->mapx->f, 1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) omega->mapy->f, 1, 0, 0, 0, 0);
#endif

#endif

    int n = iparam ("NFILTER");

    for (i = 0; i < nztot; i++) {
      omega->mapx->f[i] = (fpk[0][i]+kinvis*fvk[0][i]);
      omega->mapy->f[i] = (fpk[1][i]+kinvis*fvk[1][i]);
//      ROOT fprintf(stderr,"mapx = %g   mapy = %g  \n",omega->mapx->f[i],omega->mapy->f[i]);
#ifdef AXIAL
      omega->mapz->f[i] = (fpk[2][i]+kinvis*fvk[2][i]);
#endif	  
    }
    
    if(n>0)  
    for (i = nztot-n; i < nztot; i++) {
      omega->mapx->f[i] = 0.;
      omega->mapy->f[i] = 0.;
#ifdef AXIAL
      omega->mapz->f[i] = 0.;
#endif	  
    }
#endif
  }
 
//force on additional surfaces
  double   fp1[3] = {0., 0., 0.},fv1[3] = {0., 0., 0.};
  for(B=omega->Ubc[0]; B; B = B->next)
    if(B->usrtype == 'T'){
	 trip = 1;
 	 U    = omega->U;
	 face = B->face;
	 eid  = B->elmt->id;
	 qa   = U->flist[eid]->qa;

	 if(U->flist[eid]->lmax != lmax){
	  b    = U->flist[eid]->getbasis();
	  lmax = U->flist[eid]->lmax;
	  getzw(qa,&w,&w,'a');
	}

	eU = omega->U->flevels[0]->flist[eid];
	eV = omega->V->flevels[0]->flist[eid];
	eW = omega->W->flevels[0]->flist[eid];
	eP = omega->P->flevels[0]->flist[eid];

	eP->Jbwd(eP,b); eP->state = 't';
	eU->Jbwd(eU,b); eU->state = 't';
	eV->Jbwd(eV,b); eV->state = 't';
	eW->Jbwd(eW,b); eW->state = 't';

	eU->Grad_d(ux, uy, 0, 'a');
	eV->Grad_d(vx, vy, 0, 'a');
	eW->Grad_d(wx, wy, 0, 'a');

	/* get appropriate face values from grad(V) */
	for(i = 0; i < 6; ++i){
	  eU->GetFace(D[i] ,face,wk);
	  eU->InterpToFace1(face,wk,D[i]);
	}

	eP->GetFace(eP->h[0],face,wk);
	eP->InterpToFace1(face,wk,p);

	// calculate modal forces - nz,uz,vz and wz are all zero in four. exp

	if(eU->curvX)
	  for(i = 0; i < qa; ++i){
	    fp1[0] += p[i]*B->nx.p[i]*B->sjac.p[i]*w[i];
	    fp1[1] += p[i]*B->ny.p[i]*B->sjac.p[i]*w[i];
	    fp1[2] += 0.0;

	    fv1[0] -= (2.0*ux[i]*B->nx.p[i]+(uy[i]+vx[i])*
		      B->ny.p[i])*B->sjac.p[i]*w[i];
	    fv1[1] -= (2.0*vy[i]*B->ny.p[i]+(uy[i]+vx[i])*
		      B->nx.p[i])*B->sjac.p[i]*w[i];
	    fv1[2] -= (wx[i]*B->nx.p[i] + wy[i]*
		      B->ny.p[i])*B->sjac.p[i]*w[i];

	  }
	else
	  for(i = 0; i < qa; ++i){
	    fp1[0] += p[i]*B->nx.d*B->sjac.d*w[i];
	    fp1[1] += p[i]*B->ny.d*B->sjac.d*w[i];
	    fp1[2] += 0.0;

	    fv1[0] -= (2.0*ux[i]*B->nx.d+(uy[i] + vx[i])*
		      B->ny.d)*B->sjac.d*w[i];
	    fv1[1] -= (2.0*vy[i]*B->ny.d+(uy[i] + vx[i])*
		      B->nx.d)*B->sjac.d*w[i];
	    fv1[2] -= (wx[i]*B->nx.d + wy[i]*
		      B->ny.d)*B->sjac.d*w[i];
	  }
      }

    for(i = 0; i < 3; ++i)
      fv1[i] *= kinvis;

  //  if(trip)
    {
    ROOT {
      if (!((step-1)%hisstep))
     {
       FILE *fout1;
	     char *buf1 = (char*) calloc(BUFSIZ,sizeof(char));
	     sprintf(buf1, "HydroForce_T.fce");
       fout1 = fopen(buf1,"aw");
       fprintf(fout1,"%lf ( %lf %lf ) %lf ( %lf %lf ) %lf ( %lf %lf ) %lf \n",
	      time,fp1[0],fv1[0],fp1[0]+fv1[0],fp1[1],fv1[1],fp1[1]+fv1[1],fp1[2],
	      fv1[2],fp1[2]+fv1[2]);
        
        fflush(fout1);
        fclose(fout1);
        free(buf1);
      }
     }
    }

  double   fp2[3] = {0., 0., 0.},fv2[3] = {0., 0., 0.};
  for(B=omega->Ubc[0]; B; B = B->next)
    if(B->usrtype == 'X'){
	 trip = 1;
 	 U    = omega->U;
	 face = B->face;
	 eid  = B->elmt->id;
	 qa   = U->flist[eid]->qa;

	 if(U->flist[eid]->lmax != lmax){
	  b    = U->flist[eid]->getbasis();
	  lmax = U->flist[eid]->lmax;
	  getzw(qa,&w,&w,'a');
	}

	eU = omega->U->flevels[0]->flist[eid];
	eV = omega->V->flevels[0]->flist[eid];
	eW = omega->W->flevels[0]->flist[eid];
	eP = omega->P->flevels[0]->flist[eid];

	eP->Jbwd(eP,b); eP->state = 't';
	eU->Jbwd(eU,b); eU->state = 't';
	eV->Jbwd(eV,b); eV->state = 't';
	eW->Jbwd(eW,b); eW->state = 't';

	eU->Grad_d(ux, uy, 0, 'a');
	eV->Grad_d(vx, vy, 0, 'a');
	eW->Grad_d(wx, wy, 0, 'a');

	/* get appropriate face values from grad(V) */
	for(i = 0; i < 6; ++i){
	  eU->GetFace(D[i] ,face,wk);
	  eU->InterpToFace1(face,wk,D[i]);
	}

	eP->GetFace(eP->h[0],face,wk);
	eP->InterpToFace1(face,wk,p);

	// calculate modal forces - nz,uz,vz and wz are all zero in four. exp

	if(eU->curvX)
	  for(i = 0; i < qa; ++i){
	    fp2[0] += p[i]*B->nx.p[i]*B->sjac.p[i]*w[i];
	    fp2[1] += p[i]*B->ny.p[i]*B->sjac.p[i]*w[i];
	    fp2[2] += 0.0;

	    fv2[0] -= (2.0*ux[i]*B->nx.p[i]+(uy[i]+vx[i])*
		      B->ny.p[i])*B->sjac.p[i]*w[i];
	    fv2[1] -= (2.0*vy[i]*B->ny.p[i]+(uy[i]+vx[i])*
		      B->nx.p[i])*B->sjac.p[i]*w[i];
	    fv2[2] -= (wx[i]*B->nx.p[i] + wy[i]*
		      B->ny.p[i])*B->sjac.p[i]*w[i];
	  }
	else
	  for(i = 0; i < qa; ++i){
	    fp2[0] += p[i]*B->nx.d*B->sjac.d*w[i];
	    fp2[1] += p[i]*B->ny.d*B->sjac.d*w[i];
	    fp2[2] += 0.0;

	    fv2[0] -= (2.0*ux[i]*B->nx.d+(uy[i] + vx[i])*
		      B->ny.d)*B->sjac.d*w[i];
	    fv2[1] -= (2.0*vy[i]*B->ny.d+(uy[i] + vx[i])*
		      B->nx.d)*B->sjac.d*w[i];
	    fv2[2] -= (wx[i]*B->nx.d + wy[i]*
		      B->ny.d)*B->sjac.d*w[i];
	  }
      }

    for(i = 0; i < 3; ++i)
      fv2[i] *= kinvis;

//    if(trip)
    {
    ROOT {
      if (!((step-1)%hisstep))
     {
       FILE *fout2;
	     char *buf2 = (char*) calloc(BUFSIZ,sizeof(char));
	     sprintf(buf2, "HydroForce_X.fce");
       fout2 = fopen(buf2,"aw");
      fprintf(fout2,"%lf ( %lf %lf ) %lf ( %lf %lf ) %lf ( %lf %lf ) %lf \n",
	      time,fp2[0],fv2[0],fp2[0]+fv2[0],fp2[1],fv2[1],fp2[1]+fv2[1],fp2[2],
	      fv2[2],fp2[2]+fv2[2]);
        
        fflush(fout2);
        fclose(fout2);
        free(buf2);
      }
     }
    }

  free(wk); free_dmatrix(D,0,0);
  free(workspace);
  return;
}

void set_SurGeofac(Bndry *Pbc, Bndry *Ubc){
  Bndry *Ebc;

  for(;Pbc; Pbc = Pbc->next)
    if(Pbc->type == 'F')
      for(Ebc = Ubc; Ebc; Ebc = Ebc->next)
	if((Ebc->elmt->id == Pbc->elmt->id)&&(Ebc->face == Pbc->face)){
	  Ebc->sjac = Pbc->sjac;
	  Ebc->nx   = Pbc->nx;
	  Ebc->ny   = Pbc->ny;
	} 
}
#endif // ifdef FORCES
