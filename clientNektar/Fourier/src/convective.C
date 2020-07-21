/* ------------------------------------------------------------------------- *
 * VdgradV() - Calculate convective form of the nonlinear terms              *
 *                                                                           *
 * The following function calculates the nonlinear portion of the Navier-    *
 * Stokes equation in convective form in order to minimize aliasing errors.  *
 * In this form, the non-linear part looks like:                             *
 *                                                                           *
 *                                      dV_i                                 *
 *                         N_i = -( V_j ---- )                               *
 *                                      dx_j                                 *
 *                                                                           *
 * RCS Information                                                           *
 * ---------------                                                           
 * $Author:
 * $Date: 
 * $Source: 
 * $Revision:
 * ------------------------------------------------------------------------- */

#include "nektarF.h"
#include <math.h>
#include <rfftw.h>

#ifdef SPM
#include "C_SPM.h"
#endif

#ifdef THERMO
void compute_temperature_range(Domain *omega);
void compute_temperature_variation (Domain *omega);

#ifdef NATURAL_CONVECTION
void add_buoyancy_force(Domain *omega);
#endif

#ifdef CONJUGATE_HEAT

#ifdef SPM
void add_heat_source(Domain *omega, C_SPM *SPMv);
#else
void add_heat_source(Domain *omega);
#endif

#endif

#endif

extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;

static double *grad_wk=(double*)0, *workspace_end;
extern double **u2nd;
void   addfields     (Element_List *V[]);

void VdgradV(Domain *omega)
{
  Element_List *U  =  omega->U,  *V  =  omega->V,  *W  = omega->W,
               *Ux =  omega->P,  *Uy = omega->Uf, *Uz = omega->Uf,
               *Vx =  omega->P,  *Vy = omega->Vf, *Vz = omega->Vf,
               *Wx =  omega->P,  *Wy = omega->Wf, *Wz = omega->Wf,
#ifdef THERMO
               *T  =  omega->T,  *Tf = omega->Tf,
               *Tx =  omega->P,  *Ty = omega->Tf, *Tz = omega->Tf,
#endif
               *nul = (Element_List*)0, *Div = omega->P;
#ifdef MAP
  Element_List *Wlast = omega->Wlast;
#endif
  int           nq = U->htot*U->nz;
  int           nq32 = 3*U->htot*U->nz/2;
  char          trip = 'a';

  static Nek_Trans_Type f_to_p = F_to_P,
                        p_to_f = P_to_F;

  Element_List *P  =  omega->P;
  Element_List *Ut  =  omega->Ut,  *Vt  =  omega->Vt,  *Wt  = omega->Wt;
  double      *afx0 = omega->afx[0], *afx1 = omega->afx[1];
  double      *afy0 = omega->afy[0], *afy1 = omega->afy[1];
  double      *afz0 = omega->afz[0], *afz1 = omega->afz[1];
#ifdef THERMO
  double      *aft0 = omega->aft[0], *aft1 = omega->aft[1];
#endif

//#if defined(THERMO) && defined(MAP)
//              fprintf(stderr," Thermo && MAP doesn't work right now ! \n");
//              exit(-1);
//#endif

  static int step = 1, hisstep = option("hisstep");
  static int timestep = 0;
  static double time = dparam("STARTIME");
  static double evm_time = dparam("STARTIME");

#ifdef MAP
  static int not_first = 1;
  static double *Fx, *Fy, *Fz, *tmpstore;
#ifdef THERMO
  static double *Ft; //for temperature
#endif
  register int i, j, k, plane;
#endif

  int worksize = nq;
  if (option("dealias"))
    worksize = 3*nq/2;

  if(!grad_wk){
    if (option("dealias")) {
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
 //     worksize = 3*nq/2;
    }
 #ifndef THERMO
    grad_wk = dvector(0,worksize+4*nq-1);
 #else
    grad_wk = dvector(0,worksize+5*nq-1);
 #endif
    workspace_end = grad_wk + worksize;

    //    dzero(worksize,grad_wk,1);
#ifdef MAP
 #ifndef THERMO
    Fx       = dvector(0,4*worksize-1);
 #else
    Fx       = dvector(0,5*worksize-1);
 #endif
    Fy       = Fx + worksize;
    Fz       = Fy + worksize;
 #ifndef THERMO
    tmpstore = Fz + worksize;
 #else
    Ft       = Fz + worksize;
    tmpstore = Ft + worksize;
 #endif
    if (Wlast->fhead->state == 't') {// This only needs to be done after IOSTEP
      Wlast->Trans (Wlast, J_to_Q);
      Wlast->Trans(Wlast, f_to_p);
    }
#endif
  }

  dcopy(worksize, afx0, 1, afx1, 1);// afx[1] stores the force at previous time step 
  dcopy(worksize, afy0, 1, afy1, 1);// afx[1] stores the force at previous time step
  dcopy(worksize, afz0, 1, afz1, 1);// afx[1] stores the force at previous time step

  dzero(worksize, afx0, 1);
  dzero(worksize, afy0, 1);
  dzero(worksize, afz0, 1);

#ifdef THERMO
  dcopy(worksize, aft0, 1, aft1, 1);// aft[1] stores the force at previous time step
  dzero(worksize, aft0, 1);
#endif

  dcopy(nq, U->base_h, 1, workspace_end     , 1); //back value in (Q,F)-space
  dcopy(nq, V->base_h, 1, workspace_end+  nq, 1);
  dcopy(nq, W->base_h, 1, workspace_end+2*nq, 1);
  
/********the following part is for computing adjust force and statistics*******/  
//  if (P->fhead->state == 't')
  char save_state_P = P->fhead->state;
  P->Set_state('p');
  P->Trans (P, J_to_Q);
  dcopy(nq, P->base_h, 1, workspace_end+3*nq, 1);
 
 #ifdef THERMO
  dcopy(nq, T->base_h, 1, workspace_end+4*nq, 1);
 #endif

   double adjust_force = 0.;
  if( dparam("ADJUSTFORCE") )
    {
      adjust_force = compute_adjustment_forces(omega);
      omega->adjust_force = adjust_force;

      ROOT fprintf(stdout,"Adjust force = %lf \n", adjust_force);
       
//      if(dparam("PIPE"))
//         dsadd(nq, adjust_force, afz0, 1, afz0, 1);
//       if(dparam("CHANNEL"))
//         dsadd(nq, adjust_force, afx0, 1, afx0, 1);
    }

/*for dealiasing mode, only nq of worksize is used */
/*
  U->Trans(U, f_to_p); 
  V->Trans(V, f_to_p); 
  W->Trans(W, f_to_p); // here we are in phsyical/quadrature space 
#ifdef THERMO
  T->Trans(T, f_to_p); // here we are in phsyical/quadrature space
#endif
  P->Trans(P, f_to_p); // here we are in phsyical/quadrature space 

  
#if 0
  double fx = dparam("FFX");
  double fy = dparam("FFY");
  double fz = dparam("FFZ");
#ifdef THERMO
  double ft = dparam("FFT");
#endif
  if(fx) dsadd(nq, fx, afx0, 1,afx0, 1);
  if(fy) dsadd(nq, fy, afy0, 1,afy0, 1);
  if(fz) dsadd(nq, fz, afz0, 1,afz0, 1);
#ifdef THERMO
  if(ft) dsadd(nq, ft, aft0, 1,aft0, 1);

#ifdef NATURAL_CONVECTION
   dvadd(nq, T->base_h, 1, afx0, 1, afx0, 1);
   dvadd(nq, T->base_h, 1, afy0, 1, afy0, 1);
   dvadd(nq, T->base_h, 1, afz0, 1, afz0, 1);
#endif

#endif

#endif
*/

  dcopy(worksize,omega->uk[1],1,omega->uk[0],1); // u values at n-1 step
	dcopy(worksize,omega->vk[1],1,omega->vk[0],1);
	dcopy(worksize,omega->wk[1],1,omega->wk[0],1);
#ifdef THERMO
	dcopy(worksize,omega->tk[1],1,omega->tk[0],1);
#endif
	dcopy(worksize,omega->pk[1],1,omega->pk[0],1);

//  dcopy(worksize,U->base_h,1,omega->uk[1],1);
//	dcopy(worksize,V->base_h,1,omega->vk[1],1);
//	dcopy(worksize,W->base_h,1,omega->wk[1],1);
//#ifdef THERMO
//	dcopy(worksize,T->base_h,1,omega->tk[1],1);
//#endif
//	dcopy(worksize,P->base_h,1,omega->pk[1],1);

  U->Trans(U->htot, U->base_h, omega->uk[1], f_to_p); // in 3/2(Q,P)-space
  V->Trans(V->htot, V->base_h, omega->vk[1], f_to_p); // in 3/2(Q,P)-space
  W->Trans(W->htot, W->base_h, omega->wk[1], f_to_p); // in 3/2(Q,P)-space
#ifdef THERMO
  T->Trans(T->htot, T->base_h, omega->tk[1], f_to_p); // in 3/2(Q,P)-space
#endif
  P->Trans(P->htot, P->base_h, omega->pk[1], f_to_p); // in 3/2(Q,P)-space
 // dcopy(nq, P->base_h, 1, Pqstore, 1);

  compute_velocity_range(omega);
#ifdef THERMO

  compute_temperature_range(omega);
#endif

#ifdef ENTROPYVISCOSITY
  compute_velocity_entropy_variation(omega);
#ifdef THERMO
//  ROOTONLY fprintf(stderr, "warning: dealiasing mode needs debugging following for heat transfer ...\n");
//  exit(-1);
  compute_temperature_variation(omega);
#endif

#endif
 
  int ev_step = (dparam("EV_STEP") ? dparam("EV_STEP") : 2);
#ifdef ENTROPYVISCOSITY
   if(timestep < ev_step)
#endif
  ROOTONLY fprintf(stdout,"Maximum velocity ...... ............. %g \n", omega->maximal_velocity);

//  if( dparam("STATISTICS") )
  if( (dparam("STATISTICS")) &&(timestep%hisstep == 0) )
//  if( dparam("STATISTICS") )
  {
    if(dparam("CHANNEL"))
    statistics_output(omega,timestep,evm_time);
    if(dparam("CYLINDER") || dparam("PIPE") || dparam("SPHERE") )
    output_point_value(omega,timestep,evm_time);
  }
      
//copy back the Fourier coefficients
//  dcopy(nq, workspace_end,       1,  U->base_h, 1);
//  dcopy(nq, workspace_end+   nq, 1,  V->base_h, 1);
//  dcopy(nq, workspace_end+ 2*nq, 1,  W->base_h, 1);
// dcopy(nq, workspace_end+ 3*nq, 1,  P->base_h, 1);
//#ifdef THERMO
//  dcopy(nq, workspace_end+ 4*nq, 1,  T->base_h, 1);
//#endif
  
  dzero(worksize, Ut->base_h, 1);
  dzero(worksize, Vt->base_h, 1);
  dzero(worksize, Wt->base_h, 1);
//#endif
/****************end of statistics and adjust force*************************/



  if (step == 1 && option("verbose")) {
    // Take Ux, Uz, Vy, Vz, Wz derivatives in fourier/quadrature space
    // Can get Ux & Uz at the same time (P isn't used anywhere else)
    U->Grad(Ux,nul,Uz,'a'); Uz->Set_state('p');
    // need to do get Vy separately due to lack of space
    V->Grad(nul,Vy,nul,'y'); 
    daxpy  (nq, 1.0, Vy->base_h,1,Ux->base_h,1); 
    // Get the z-derivative now
    V->Grad(nul,nul,Vz,'z'); Vz->Set_state('p');
    W->Grad(nul,nul,Wz,'z'); Wz->Set_state('p');
    // Add up the Wz component of divergence
    daxpy  (nq, 1.0, Wz->base_h,1,Ux->base_h,1); 
  
    Div->Set_state('p');
    double GDiv = GL2(Div);
#ifdef MAP
    ROOT printf ("Time %lf Momentum Mx = %#lf My = %#lf Mz = %#lf Mass flux = %#lf L2 Divergence = %#lf\n", time, VolInt (U, omega->mapx->t[0]), VolInt (V, omega->mapy->t[0]), VolInt (W, 0.), VolInt (Div, 0.), GDiv);
#else
    ROOT printf ("Time %lf Momentum Mx = %#lf My = %#lf Mz = %#lf Mass flux = %#lf L2 Divergence = %#lf\n", time, VolInt (U), VolInt (V), VolInt (W), VolInt (Div), GDiv);
#endif
    //ROOT printf ("Time %lf Momentum Mx = %#lf My = %#lf Mz = %#lf Mass flux = %#lf L2 Divergence = %#lf\n", time, VolInt (U), VolInt (V), VolInt (W), VolInt (Div), GDiv);
    Div->Set_state('t');
    time += hisstep*dparam("DELT");
  } else {
    // Take z derivatives in fourier/quadrature space
    U->Grad(nul,nul,Uz,'z'); Uz->Set_state('p');
    V->Grad(nul,nul,Vz,'z'); Vz->Set_state('p');
    W->Grad(nul,nul,Wz,'z'); Wz->Set_state('p');
#ifdef THERMO
    T->Grad(nul,nul,Tz,'z'); Tz->Set_state('p');
#endif
  }

  // Transform velocity to physical/quadrature space
  U->Trans(U, f_to_p); 
  V->Trans(V, f_to_p); 
  W->Trans(W, f_to_p); 
#ifdef THERMO
  T->Trans(T, f_to_p); 
#endif

//#if defined(MAP) && defined(SPM) //concentration should also be in 3/2 space
//  omega->CONCENTR->Trans(omega->CONCENTR, P_to_F); //go to (F,P) 
//  omega->CONCENTR->Trans(omega->CONCENTR, f_to_p); //go to 3/2*nq space
//#endif

  // Transform the pressure
  P->Trans(P, f_to_p);
  // Transform z derivative of velocity to physical/quadrature space
  Uz->Trans(Uz, f_to_p);
  Vz->Trans(Vz, f_to_p);
  Wz->Trans(Wz, f_to_p);
#ifdef THERMO
  Tz->Trans(Tz, f_to_p); 
#endif

  if (option("dealias")) {
    trip = 'A';
    nq   = 3*nq/2;
  }

//  if(option("STATAVG"))
  if(option("timeavg"))
   {

   if( dparam("PIPE") || dparam("CYLINDER") )
    {
#ifndef THERMO
     Element_List  **VF = (Element_List**) malloc((4)*sizeof(Element*));
#else
     Element_List  **VF = (Element_List**) malloc((5)*sizeof(Element*));
#endif

     VF[0]   = omega->U;  VF[1]   = omega->V;  VF[2]   = omega->W;
     VF[3]   = omega->P;
#ifdef THERMO
     VF[4]   = omega->T;
#endif

     addfields(VF);
    }
    else
     {
     fprintf(stderr,"need to debug for dealiasing .... \n");
     exit(-1);
     dvvtvp(nq, U->base_h, 1, U->base_h, 1, u2nd[0], 1, u2nd[0], 1);
     dvvtvp(nq, V->base_h, 1, V->base_h, 1, u2nd[1], 1, u2nd[1], 1);
     dvvtvp(nq, W->base_h, 1, W->base_h, 1, u2nd[2], 1, u2nd[2], 1);
     dvvtvp(nq, U->base_h, 1, V->base_h, 1, u2nd[3], 1, u2nd[3], 1);
#ifdef THERMO
     dvvtvp(nq, U->base_h, 1, T->base_h, 1, u2nd[4], 1, u2nd[4], 1);
#endif
     }
  }

#ifdef SMAGORINSKY
  if (iparam("ISO")!= 0)
    IsoStats(omega);
#endif

#ifdef MAP
  // Get dW/dt
  Mapping *mapx = omega->mapx,
      *mapy = omega->mapy;

  double *u = U->base_h,
         *v = V->base_h,
         *w = W->base_h;
  double nu = dparam("KINVIS");
#ifdef THERMO
  double kappa = dparam("KINVIS")/dparam("PRANDTL");
#endif
  double *dwdt   = tmpstore,
         *Uz_sav = tmpstore,
         *Vz_sav = tmpstore,
         *Wz_sav = tmpstore;
 #ifdef THERMO
  double *Tz_sav = tmpstore;
 #endif

  dvsub (nq, W->base_h, 1, Wlast->base_h, 1, dwdt, 1); 
  dscal (nq, 1./dparam("DELT"), dwdt, 1);
  if (!not_first) {
    // Fix in the absence of Wlast in the restart file:
    // for stability zero-out the value of dwdt on the first timestep
    if (dsum(nq, Wlast->base_h, 1) == 0.0) 
      dzero (nq, dwdt, 1);
    not_first++;
  }
  dcopy (nq, W->base_h, 1, Wlast->base_h, 1);

  int NZ    = U->nz;
  int NZTOT = U->nztot;
  int nZtot = NZTOT; // keep this value as it may change because of dealiasing
  // Set the RFFT structs
  rfftw_plan CPlan     = rplan, 
             CPlan_inv = rplan_inv;
  if (option("dealias")) {
    NZ    = 3*NZ/2;
    NZTOT = 3*NZTOT/2;
    CPlan     = rplan32; 
    CPlan_inv = rplan_inv32;
  }

  // Take the maps to physical space

  if (option("dealias")) {
    // zero out for dealiasing
    dzero (NZTOT-nZtot,   mapx->z+nZtot, 1);
    dzero (NZTOT-nZtot,   mapy->z+nZtot, 1);
    dzero (NZTOT-nZtot,  mapx->zz+nZtot, 1);
    dzero (NZTOT-nZtot,  mapy->zz+nZtot, 1);
#ifdef TMAP
    dzero (NZTOT-nZtot,  mapx->tz+nZtot, 1);
    dzero (NZTOT-nZtot,  mapy->tz+nZtot, 1);
    dzero (NZTOT-nZtot, mapx->tzz+nZtot, 1);
    dzero (NZTOT-nZtot, mapy->tzz+nZtot, 1);
    dzero (NZTOT-nZtot,  mapx->tt+nZtot, 1);
    dzero (NZTOT-nZtot,  mapy->tt+nZtot, 1);
#endif
  }
  // invFFT the maps to physical space
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)   mapx->z, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)   mapy->z, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)  mapx->zz, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)  mapy->zz, 1, 0, 0, 0, 0);
#ifdef TMAP
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)  mapx->tz, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)  mapy->tz, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *) mapx->tzz, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *) mapy->tzz, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)  mapx->tt, 1, 0, 0, 0, 0);
  rfftw(CPlan_inv, 1, (FFTW_COMPLEX *)  mapy->tt, 1, 0, 0, 0, 0);
#endif

  double *Px = Fx,
         *Py = Fy;
 
  P->Grad_h(P->base_h, Px, Py, (double *)NULL, trip);
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++)
    for (i = 0, j = k*U->htot; i < U->htot; i++, j++)
      Fz[j] = mapx->z[plane] * Px[j] + mapy->z[plane] * Py[j];
#ifdef TMAP
  if (!option("acc_impl")) {
    for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++)
      for (i = 0, j = k*U->htot; i < U->htot; i++, j++) {
	Fx[j] = mapx->tzz[plane] * nu - mapx->tt[plane] - 
	        mapx->z[plane] * dwdt[j] - mapx->zz[plane] * w[j]*w[j] - 
                2. * mapx->tz[plane] * w[j];
	Fy[j] = mapy->tzz[plane] * nu - mapy->tt[plane] - 
	        mapy->z[plane] * dwdt[j] - mapy->zz[plane] * w[j]*w[j] - 
	        2. * mapy->tz[plane] * w[j];	        
      }
  } else {
    for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++)
      for (i = 0, j = k*U->htot; i < U->htot; i++, j++) {
	Fx[j] = -(mapx->z[plane] * dwdt[j] + mapx->zz[plane] * w[j]*w[j] +
		  2. * mapx->tz[plane] * w[j]);
	Fy[j] = -(mapy->z[plane] * dwdt[j] + mapy->zz[plane] * w[j]*w[j] +
		  2. * mapy->tz[plane] * w[j]);
      }
  }
#else // ifdef TMAP
    for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++)
      for (i = 0, j = k*U->htot; i < U->htot; i++, j++) {
	Fx[j] = -(mapx->z[plane] * dwdt[j] + mapx->zz[plane] * w[j]*w[j]);
	Fy[j] = -(mapy->z[plane] * dwdt[j] + mapy->zz[plane] * w[j]*w[j]);
      }
#endif // ifdef TMAP

  // Uf = U.Grad u

  // save Uz
  dcopy  (nq, Uz->base_h, 1, Uz_sav, 1);
  // form (u.nabla) u
  dvmul  (nq, W->base_h,1,Uz->base_h,1,grad_wk,1);
  U->Grad_h(U->base_h, Ux->base_h,Uy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Ux->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Uy->base_h,1,grad_wk,1,grad_wk,1);
  
//#ifdef SPM //zwang 07192017, no acceleration inside solid not sure about it right now !
//  for(i=0; i<nq; ++i)
//    Fx[i] = (1.-omega->CONCENTR->base_h[i])*Fx[i];
//#endif

  // subtract the (u.nabla) u term from Fx
  daxpy  (nq, -1.0, grad_wk, 1, Fx, 1);

//#ifdef ENTROPYVISCOSITY
// dcopy  (nq, grad_wk, 1, Ut->base_h, 1); //(u, \nabla) \cdot u
//#endif
  // form Uz - mapx_z Ux - mapy_z Uy
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (U->htot, -mapx->z[plane], Ux->base_h+j, 1, Uz_sav+j, 1);
    daxpy (U->htot, -mapy->z[plane], Uy->base_h+j, 1, Uz_sav+j, 1);
  }
  // Store it in Uz
  dcopy  (nq, Uz_sav, 1, Uz->base_h, 1); 

  // Vf = U.Grad v

  // save Vz
  dcopy  (nq, Vz->base_h, 1, Vz_sav, 1);
  // form (u.nabla) v
  dvmul  (nq, W->base_h,1,Vz->base_h,1,grad_wk,1);
  V->Grad_h(V->base_h,Vx->base_h,Vy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Vx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Vy->base_h,1,grad_wk,1,grad_wk,1);

//#ifdef SPM //zwang 07192017, no acceleration inside solid not sure about it right now !
//  for(i=0; i<nq; ++i)
//    Fy[i] = (1.-omega->CONCENTR->base_h[i])*Fy[i];
//#endif
  // subtract the (u.nabla) v term from Fy
  daxpy  (nq, -1.0, grad_wk, 1, Fy, 1);

//#ifdef ENTROPYVISCOSITY
//  dcopy  (nq, grad_wk, 1, Vt->base_h, 1);
//#endif

  // form Vz - mapx_z Vx - mapy_z Vy
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (V->htot, -mapx->z[plane], Vx->base_h+j, 1, Vz_sav+j, 1);
    daxpy (V->htot, -mapy->z[plane], Vy->base_h+j, 1, Vz_sav+j, 1);
  }
  // Store it in Vz
  dcopy  (nq, Vz_sav, 1, Vz->base_h, 1);

//#ifdef SPM //zwang 07192017, no acceleration inside solid not sure about it right now !
//  for(i=0; i<nq; ++i)
//    Vz->base_h[i] = (1.-omega->CONCENTR->base_h[i])*Vz->base_h[i];
//#endif

  // Wf = U.Grad w
  dcopy  (nq, Wz->base_h, 1, Wz_sav, 1);
  // form (u.nabla) w
  dvmul  (nq, W->base_h,1,Wz->base_h,1,grad_wk,1);
  W->Grad_h(W->base_h, Wx->base_h,Wy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Wx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Wy->base_h,1,grad_wk,1,grad_wk,1);

//#ifdef SPM //zwang 07192017, no acceleration inside solid not sure about it right now !
//  for(i=0; i<nq; ++i)
//   Fz[i] = (1.-omega->CONCENTR->base_h[i])*Fz[i];
//#endif

  // subtract the (u.nabla) w term from Fz
  daxpy  (nq, -1.0, grad_wk, 1, Fz, 1);

//#ifdef ENTROPYVISCOSITY
//  dcopy  (nq, grad_wk, 1, Wt->base_h, 1);
//#endif
  // I am not sure how these two terms come to the Fx,Fy. 
  // I add the term related to z by following Fx and Fy
  // subtract the (map[xy]_z (u.nabla) w) term from F[xy]
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (W->htot, -mapx->z[plane], grad_wk+j, 1, Fx+j, 1);
    daxpy (W->htot, -mapy->z[plane], grad_wk+j, 1, Fy+j, 1);
  }
  // we can now reuse grad_wk

#ifdef THERMO

  dcopy  (nq, Tz->base_h, 1, Tz_sav, 1);
  // form (u.nabla) T
  dvmul  (nq, W->base_h,1,Tz->base_h,1,grad_wk,1);
  T->Grad_h(T->base_h, Tx->base_h,Ty->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Tx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Ty->base_h,1,grad_wk,1,grad_wk,1); //grad_wk stores the convective term of temperature

  dsmul (nq, -1.0, grad_wk, 1, Ft, 1); //Ft stores the convective term

  // form Tz - mapx_z Tx - mapy_z Ty
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (T->htot, -mapx->z[plane], Tx->base_h+j, 1, Tz_sav+j, 1);
    daxpy (T->htot, -mapy->z[plane], Ty->base_h+j, 1, Tz_sav+j, 1);
  }	

  dcopy  (nq, Tz_sav, 1, Tz->base_h, 1);
#endif

  // form Wz - mapx_z Wx - mapy_z Wy
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (W->htot, -mapx->z[plane], Wx->base_h+j, 1, Wz_sav+j, 1);
    daxpy (W->htot, -mapy->z[plane], Wy->base_h+j, 1, Wz_sav+j, 1);
  }	  

  // Get Wxx
  Wx->Grad_h(Wx->base_h,       grad_wk,(double *)NULL,(double *)NULL,trip);
  // Get Wyy
  Wy->Grad_h(Wy->base_h,(double *)NULL,    Wx->base_h,(double *)NULL,trip);
  // Store Wxx+Wyy in grad_wk
  daxpy  (nq, 1.0, Wx->base_h, 1, grad_wk, 1);
  // add the (map[xy]_z nu (Wxx+Wyy)) term to F[xy]
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (W->htot, mapx->z[plane]*nu, grad_wk+j, 1, Fx+j, 1);
    daxpy (W->htot, mapy->z[plane]*nu, grad_wk+j, 1, Fy+j, 1);
  }

  // Store (Wz - mapx_z Wx - mapy_z Wy) in Wz
  dcopy  (nq, Wz_sav, 1, Wz->base_h, 1);


  // Add (map[xy]_z (Wz - mapx_z Wx - mapy_z Wy) + map[xy]_zz W) to [UV]z
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (W->htot, mapx->z[plane], Wz->base_h+j, 1, Uz->base_h+j, 1);
    daxpy (W->htot, mapx->zz[plane], W->base_h+j, 1, Uz->base_h+j, 1);
    daxpy (W->htot, mapy->z[plane], Wz->base_h+j, 1, Vz->base_h+j, 1);
    daxpy (W->htot, mapy->zz[plane], W->base_h+j, 1, Vz->base_h+j, 1);
  }

  Uz->Grad_h(Uz->base_h, U->base_h,Ux->base_h,(double *)NULL,trip);

//#ifdef SPM //zwang 07192017, no acceleration inside solid not sure about it right now !
//  for(i=0; i<nq; ++i) {
//    U->base_h[i]  = (1.-omega->CONCENTR->base_h[i])*U->base_h[i];
//    Ux->base_h[i] = (1.-omega->CONCENTR->base_h[i])*Ux->base_h[i];
//  }
//#endif
  // Do the final subtraction from Fx
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (U->htot, -mapx->z[plane]*nu, U->base_h+j, 1, Fx+j, 1);
    daxpy (U->htot, -mapy->z[plane]*nu, Ux->base_h+j, 1, Fx+j, 1);
  }

  Vz->Grad_h(Vz->base_h, V->base_h,Vx->base_h,(double *)NULL,trip);

//#ifdef SPM //zwang 07192017, no acceleration inside solid not sure about it right now !
//  for(i=0; i<nq; ++i) {
//    V->base_h[i]  = (1.-omega->CONCENTR->base_h[i])*V->base_h[i];
//    Vx->base_h[i] = (1.-omega->CONCENTR->base_h[i])*Vx->base_h[i];
//  }
//#endif

  // Do the final subtraction from Fy
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (V->htot, -mapx->z[plane]*nu, V->base_h+j, 1, Fy+j, 1);
    daxpy (V->htot, -mapy->z[plane]*nu, Vx->base_h+j, 1, Fy+j, 1);
  }
  Wz->Grad_h(Wz->base_h, W->base_h,Wx->base_h,(double *)NULL,trip);

//#ifdef SPM //zwang 07192017, no acceleration inside solid not sure about it right now !
//  for(i=0; i<nq; ++i) {
//    W->base_h[i]  = (1.-omega->CONCENTR->base_h[i])*W->base_h[i];
//    Wx->base_h[i] = (1.-omega->CONCENTR->base_h[i])*Wx->base_h[i];
//  }
//#endif

  // Do the final subtraction from Fz
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (W->htot, -mapx->z[plane]*nu, W->base_h+j, 1, Fz+j, 1);
    daxpy (W->htot, -mapy->z[plane]*nu, Wx->base_h+j, 1, Fz+j, 1);
  }

#ifdef THERMO

  Tz->Grad_h(Tz->base_h, T->base_h,Tx->base_h,(double *)NULL,trip);

  // Do the final subtraction from Tz
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    j = k*U->htot;
    daxpy (T->htot, -mapx->z[plane]*kappa, T->base_h+j, 1, Ft+j, 1);
    daxpy (T->htot, -mapy->z[plane]*kappa, Tx->base_h+j, 1, Ft+j, 1);
  }
#endif
  
//#ifdef ENTROPYVISCOSITY
//  dvadd (nq, Fx, 1, Ut->base_h, 1, Ut->base_h,1); //Ut->base_h has force only
//  dvadd (nq, Fy, 1, Vt->base_h, 1, Vt->base_h,1);
//  dvadd (nq, Fz, 1, Wt->base_h, 1, Wt->base_h,1);
//#endif
  // Copy the non-linear terms sofar to V and clean up F[xyz]
  dcopy  (nq, Fx, 1, U->base_h, 1); dzero (nq, Fx, 1);
  dcopy  (nq, Fy, 1, V->base_h, 1); dzero (nq, Fy, 1);
  dcopy  (nq, Fz, 1, W->base_h, 1); dzero (nq, Fz, 1);

#ifdef THERMO

  dcopy  (nq, Ft, 1, T->base_h, 1); dzero (nq, Ft, 1);
#endif

  // Transform them to Fourier Space
  U->Trans(U, p_to_f);
  V->Trans(V, p_to_f);
  W->Trans(W, p_to_f);

  // Transform the term needing further derivatives to Fourier/Quadrature space
  Uz->Trans(Uz, p_to_f);
  Vz->Trans(Vz, p_to_f);
  Wz->Trans(Wz, p_to_f);

#ifdef THERMO
  T->Trans(T, p_to_f);
  Tz->Trans(Tz, p_to_f);
#endif

  // Reset nq
  if (option("dealias"))
    nq = U->htot*U->nz;

  // Copy to F[xyz] -> replace with Trans() that works with vectors

  dcopy  (nq, U->base_h, 1, Fx, 1);
  dcopy  (nq, V->base_h, 1, Fy, 1);
  dcopy  (nq, W->base_h, 1, Fz, 1);
#ifdef THERMO
  dcopy  (nq, T->base_h, 1, Ft, 1);
#endif

  // Get the extra z-derivatives
  Uz->Grad_z(U); U->Set_state('p');
  Vz->Grad_z(V); V->Set_state('p');
  Wz->Grad_z(W); W->Set_state('p');
#ifdef THERMO
  Tz->Grad_z(T); T->Set_state('p');
#endif

  // Add to F[xyz]
  daxpy   (nq, nu, U->base_h, 1, Fx, 1);
  daxpy   (nq, nu, V->base_h, 1, Fy, 1);
  daxpy   (nq, nu, W->base_h, 1, Fz, 1);
#ifdef THERMO
  daxpy   (nq, kappa, T->base_h, 1, Ft, 1);
#endif
  // Restore Fourier fields
  dcopy(nq, workspace_end,       1,  U->base_h, 1);
  dcopy(nq, workspace_end+ nq,   1,  V->base_h, 1);
  dcopy(nq, workspace_end+ 2*nq, 1,  W->base_h, 1);
  dcopy(nq, workspace_end+ 3*nq, 1,  P->base_h, 1);
#ifdef THERMO
  dcopy(nq, workspace_end+ 4*nq, 1,  T->base_h, 1);
#endif

  // Get d^2/dz^2 terms
  Grad_zz(U, Uz);
  Grad_zz(V, Vz);
  Grad_zz(W, Wz);
#ifdef THERMO
  Grad_zz(T, Tz);
#endif

#if 0
// afx,afy,afz stored the nonlinear force at current time step.
#ifdef ENTROPYVISCOSITY
  dcopy(nq, Uz->base_h, 1, grad_wk, 1);//backup Fourier/quadrature values
  Uz->Trans(Uz, f_to_p); // go to physical

//#ifdef SPM
//  for(i=0; i<3*nq/2; ++i) //nq has been changed
//    Uz->base_h[i] = (1.- omega->CONCENTR->base_h[i])*Uz->base_h[i];
//  dcopy(3*nq/2, Uz->base_h, 1, grad_wk, 1);//backup Physical/quadrature values
//#endif
//  dsvtvp (nq, -nu, Uz->base_h, 1,Ut->base_h, 1, Ut->base_h, 1);//add to afx0
  dsvtvp (3*nq/2, -nu, Uz->base_h, 1,Ut->base_h, 1, Uz->base_h, 1);//add to afx0
  if (option("dealias")) {
  Uz->Trans(Uz, p_to_f); // from 3/2 nq to Fourier
  Uz->Trans(Uz, F_to_P); // to nq points
  }
  dcopy(nq, Uz->base_h, 1, afx0, 1);

//#ifdef SPM
//  dcopy(3*nq/2, grad_wk, 1, Uz->base_h, 1);//SPM copy back physical/quadrature values, 
//  Uz->Trans(Uz, f_to_p); // from 3/2 nq to Fourier
//#else
  dcopy(nq, grad_wk, 1, Uz->base_h, 1);//NO_SPM copy back Fourier/quadrature values, 
//#endif

  dcopy(nq, Vz->base_h, 1, grad_wk, 1);//backup Fourier/quadrature values
  Vz->Trans(Vz, f_to_p);
 //#ifdef SPM
//  for(i=0; i<3*nq/2; ++i) //nq has been changed
//    Vz->base_h[i] = (1.- omega->CONCENTR->base_h[i])*Vz->base_h[i];
//  dcopy(3*nq/2, Vz->base_h, 1, grad_wk, 1);//backup Physical/quadrature values
//#endif
//  dsvtvp (nq, -nu, Vz->base_h, 1,Vt->base_h, 1, Vt->base_h, 1);//add to afx0
  dsvtvp (3*nq/2, -nu, Vz->base_h, 1,Vt->base_h, 1, Vz->base_h, 1);//add to afx0
  if (option("dealias")) {
  Vz->Trans(Vz, p_to_f); // from 3/2 nq to Fourier
  Vz->Trans(Vz, F_to_P); // to nq points
  }
  dcopy(nq, Vz->base_h, 1, afy0, 1);

//#ifdef SPM
//  dcopy(3*nq/2, grad_wk, 1, Vz->base_h, 1);//SPM copy back physical/quadrature values, 
//  Vz->Trans(Vz, f_to_p); // from 3/2 nq to Fourier
//#else
  dcopy(nq, grad_wk, 1, Vz->base_h, 1);//copy back Fourier/quadrature values
//#endif

  dcopy(nq, Wz->base_h, 1, grad_wk, 1);//backup Fourier/quadrature values
  Wz->Trans(Wz, f_to_p);
//#ifdef SPM
//  for(i=0; i<3*nq/2; ++i) //nq has been changed
//    Wz->base_h[i] = (1.- omega->CONCENTR->base_h[i])*Wz->base_h[i];
//  dcopy(3*nq/2, Wz->base_h, 1, grad_wk, 1);//backup Physical/quadrature values
//#endif
//  dsvtvp (nq, -nu, Wz->base_h, 1,Wt->base_h, 1, Wt->base_h, 1);//add to afx0
  dsvtvp (3*nq/2, -nu, Wz->base_h, 1,Wt->base_h, 1, Wz->base_h, 1);//add to afx0
  if (option("dealias")) {
  Wz->Trans(Wz, p_to_f); // from 3/2 nq to Fourier
  Wz->Trans(Wz, F_to_P); // to nq points
  }
  dcopy(nq, Wz->base_h, 1, afz0, 1);

  dcopy(nq, grad_wk, 1, Wz->base_h, 1);//copy back Fourier/quadrature values

#endif
#endif

  // Add F[xyz] to get full non-linear terms in [UVW]z
  dsvtvp (nq, -nu, Uz->base_h, 1, Fx, 1, Uz->base_h, 1);
  dsvtvp (nq, -nu, Vz->base_h, 1, Fy, 1, Vz->base_h, 1);
  dsvtvp (nq, -nu, Wz->base_h, 1, Fz, 1, Wz->base_h, 1);
#ifdef THERMO
  dsvtvp (nq, -kappa, Tz->base_h, 1, Ft, 1, Tz->base_h, 1);
#endif


  if (step == 1 && option("verbose")) {
    double magf[3], gmagf[3], invq;
    invq = max(ddot(nq, U->base_h, 1, U->base_h, 1), 1.);
    magf[0] = ddot(nq, Uz->base_h, 1, Uz->base_h, 1)/invq;
    invq = max(ddot(nq, V->base_h, 1, V->base_h, 1), 1.);
    magf[1] = ddot(nq, Vz->base_h, 1, Vz->base_h, 1)/invq;
    invq = max(ddot(nq, W->base_h, 1, W->base_h, 1), 1.);
    magf[2] = ddot(nq, Wz->base_h, 1, Wz->base_h, 1)/invq;
    MPI_Reduce(magf, gmagf, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    ROOT {
      double dt = dparam("DELT");
      printf("\t dt|Fx| (map_z) = %#9.6e\n", dt*gmagf[0]);
      printf("\t dt|Fy| (map_z) = %#9.6e\n", dt*gmagf[1]);
      printf("\t dt|Fz| (map_z) = %#9.6e\n", dt*gmagf[2]);
    }
  }

  // FFT the maps back to fourier space
  rfftw(CPlan, 1, (FFTW_COMPLEX *)   mapx->z, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *)   mapy->z, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *)  mapx->zz, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *)  mapy->zz, 1, 0, 0, 0, 0);
#ifdef TMAP
  rfftw(CPlan, 1, (FFTW_COMPLEX *)  mapx->tz, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *)  mapy->tz, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *) mapx->tzz, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *) mapy->tzz, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *)  mapx->tt, 1, 0, 0, 0, 0);
  rfftw(CPlan, 1, (FFTW_COMPLEX *)  mapy->tt, 1, 0, 0, 0, 0);
#endif

#else //No MAP
  // Uf = U.Grad u
  dvmul  (nq, W->base_h,1,Uz->base_h,1,grad_wk,1);
  U->Grad_h(U->base_h, Ux->base_h,Uy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Ux->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Uy->base_h,1,grad_wk,1,Uz->base_h,1);
  dscal  (nq,-1.0,Uz->base_h,1);

  // Vf = U.Grad v
  dvmul  (nq, W->base_h,1,Vz->base_h,1,grad_wk,1);
  V->Grad_h(V->base_h,Vx->base_h,Vy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Vx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Vy->base_h,1,grad_wk,1,Vz->base_h,1);
  dscal  (nq,-1.0,Vz->base_h,1);

  // Wf = U.Grad w
  dvmul  (nq, W->base_h,1,Wz->base_h,1,grad_wk,1);
  W->Grad_h(W->base_h, Wx->base_h,Wy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Wx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Wy->base_h,1,grad_wk,1,Wz->base_h,1);
  dscal  (nq,-1.0,Wz->base_h,1);
#ifdef THERMO
  // Tf = U.Grad T
  dvmul  (nq, W->base_h,1,Tz->base_h,1,grad_wk,1);
  T->Grad_h(T->base_h, Tx->base_h,Ty->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Tx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Ty->base_h,1,grad_wk,1,Tz->base_h,1);
  dscal  (nq,-1.0,Tz->base_h,1);
//  for(int i=0; i<nq; ++i)
//    fprintf(stderr," T = %g  Tz = %g \n",T->base_h[i],Tz->base_h[i]);
#ifdef NATURAL_CONVECTION
  add_buoyancy_force(omega);
#endif

#ifdef CONJUGATE_HEAT
 #ifdef SPM
  add_heat_source(omega,omega->vSPM);
 #else
  add_heat_source(omega);
 #endif
#endif

#endif


//#ifdef ENTROPYVISCOSITY //save for EVM
//  dcopy(nq, Uz->base_h, 1, afx0, 1);
//  dcopy(nq, Vz->base_h, 1, afy0, 1);
//  dcopy(nq, Wz->base_h, 1, afz0, 1);
//#ifdef THERMO
//  dcopy(nq, Tz->base_h, 1, aft0, 1);
//#endif
//#endif //EVM

  // Transform U.Grad U to Fourier/Quadrature space
  Uz->Trans(Uz, p_to_f);
  Vz->Trans(Vz, p_to_f);
  Wz->Trans(Wz, p_to_f);
#ifdef THERMO
  Tz->Trans(Tz, p_to_f);
#endif

  if (option("dealias"))
    nq = U->htot*U->nz;

#endif // end of dinfine MAP
  dcopy(nq, workspace_end     , 1, U->base_h, 1);
  dcopy(nq, workspace_end+  nq, 1, V->base_h, 1);
  dcopy(nq, workspace_end+2*nq, 1, W->base_h, 1);
  dcopy(nq, workspace_end+3*nq, 1, P->base_h, 1);
#ifdef THERMO
  dcopy(nq, workspace_end+4*nq, 1, T->base_h, 1);
#endif

#ifdef ENTROPYVISCOSITY
  dcopy(nq, Uz->base_h, 1, Ut->base_h, 1);
  dcopy(nq, Vz->base_h, 1, Vt->base_h, 1);
  dcopy(nq, Wz->base_h, 1, Wt->base_h, 1);

  Ut->Trans(Ut, f_to_p);
  Vt->Trans(Vt, f_to_p);
  Wt->Trans(Wt, f_to_p);


  dcopy(nq32, Ut->base_h, 1, afx0, 1);
  dcopy(nq32, Vt->base_h, 1, afy0, 1);
  dcopy(nq32, Wt->base_h, 1, afz0, 1);

 #ifdef THERMO
  dcopy(nq, Tz->base_h, 1, Ut->base_h, 1);
  Ut->Trans(Ut, f_to_p);
  dcopy(nq32, Ut->base_h, 1, aft0, 1);
#endif

#endif

  P->Set_state(save_state_P);
  step == hisstep ? step = 1 : step++;
  timestep ++;
  evm_time += dparam("DELT");

  return;
}
