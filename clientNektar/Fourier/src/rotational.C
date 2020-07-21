/* ------------------------------------------------------------------------- *
 * VxOmega() - Calculate the nonlinear terms in rotational form              *
 *                                                                           *
 * The following function calculates the nonlinear portion of the Navier-    *
 * Stokes equation in rotational form to minimize work.  It may be expressed *
 * as:                                                                       *
 *                                                                           *
 *     N(V) = [ v.wz - w.wy ] i + [ w.wx - u.wz ] j + [ u.wy - v.wx ] k      *
 *                                                                           *
 * where V = (u,v,w) are the three velocity components and (wx,wy,wz) are    *
 * the corresponding components of the vorticty vector.                      *
 *                                                                           *
 * RCS Information
 * ---------------
 * $Author: ssherw $
 * $Date: 2004/01/28 00:04:21 $
 * $Source: /homedir/cvs/Nektar/Fourier/src/rotational.C,v $
 * $Revision: 1.1 $
 * ------------------------------------------------------------------------- */

#include "nektarF.h"

extern double **u2nd;

/* -----------------  2D Nonlinear Terms  ------------------- */
void Vorticity(Domain *omega){
  Element_List *U  = omega->U,     *V  = omega->V,     *W  = omega->W;
  Element_List *W1 = omega->Uf,    *W2 = omega->Vf,    *W3 = omega->Wf;
  Element_List *nul = (Element_List*)NULL;
  double *wk = dvector(0, U->nz*U->htot-1);

  static int step = 1, hisstep = option("hisstep");
  static double time = dparam("STARTIME");

  U->Set_state('p');
  V->Set_state('p');
  W->Set_state('p');

  if (step == 1 && option("verbose")) {
    // Take Ux, Vy, Wz derivatives in fourier/quadrature space
    U->Grad(W1,nul,nul,'x'); 
    V->Grad(nul,W2,nul,'y');  
    daxpy  (U->htot*U->nz, 1.0, W2->base_h,1,W1->base_h,1); 
    V->Grad(nul,nul,W2,'z'); 
    daxpy  (U->htot*U->nz, 1.0, W2->base_h,1,W1->base_h,1); 
  
    W1->Set_state('p');
    double GDiv = GL2(W1);
#ifdef MAP
    ROOT printf ("Time %lf Momentum Mx = %#lf My = %#lf Mz = %#lf Mass flux = %#lf L2 Divergence = %#lf\n", time, VolInt (U, omega->mapx->t[0]), VolInt (V, omega->mapy->t[0]), VolInt (W, 0.), VolInt (W1, 0.), GDiv);
#else
    ROOT printf ("Time %lf Momentum Mx = %#lf My = %#lf Mz = %#lf Mass flux = %#lf L2 Divergence = %#lf\n", time, VolInt (U), VolInt (V), VolInt (W), VolInt (W1), GDiv);
#endif
    W1->Set_state('t');
    time += hisstep*dparam("DELT");
  }

  W->Grad(nul,  W1, nul, 'y');
  V->Grad(nul, nul,  W2, 'z');
  dvsub(U->htot*U->nz, W1->base_h, 1, W2->base_h, 1, W1->base_h, 1);

  U->Grad(nul, nul,  W2, 'z');
  W->Grad( W3, nul, nul, 'x');
  dvsub(U->htot*U->nz, W2->base_h, 1, W3->base_h, 1, W2->base_h, 1);

  U->Grad(nul,  W3, nul, 'y');
  dsmul(W->htot*W->nz, -1.0,  W3->base_h, 1, wk, 1);
  V->Grad( W3, nul, nul, 'x');
  dvadd(W->htot*W->nz, W3->base_h, 1, wk, 1, W3->base_h, 1);

  free(wk);

  step == hisstep ? step = 1 : step++;

  return;
}

void FFT_AtoA(Element_List *EL, double *base_from);
void AtoA_FFT(Element_List *EL, double *base_to);
void VxOmega(Domain *omega)
{
  int    nprocs = option("NPROCS");
  int     ntot  = omega->U->htot;
  int    ntotz  = (ntot+nprocs-1)/nprocs,
          nwork = ntotz*omega->U->nztot;
	  if (option("dealias")) nwork=3*nwork/2;

  double  **u   = dmatrix(0, 2, 0, nwork-1),
          **q   = dmatrix(0, 2, 0, nwork-1),
          **res = dmatrix(0, 2, 0, nwork-1);
  register int k;
  
  Element_List *U  = omega->U, *V = omega->V, *W = omega->W;
  Element_List *Qx = omega->Uf, *Qy = omega->Vf, *Qz = omega->Wf;
  
  // Vorticity -> omega->Uf
  Vorticity(omega);
  
  AtoA_FFT (U, u[0]);    /* Transform to physical space */
  AtoA_FFT (V, u[1]);   
  AtoA_FFT (W, u[2]);  

  AtoA_FFT (Qx, q[0]);
  AtoA_FFT (Qy, q[1]);
  AtoA_FFT (Qz, q[2]);
  
  for (k = 0; k < nwork; ++k) {
    res[0][k] = u[1][k] * q[2][k] - u[2][k] * q[1][k];
    res[1][k] = u[2][k] * q[0][k] - u[0][k] * q[2][k];
    res[2][k] = u[0][k] * q[1][k] - u[1][k] * q[0][k];
  }

  if(option("STATAVG")){
#ifdef MAP
    if (omega->mstat->n) {
      double tolstat = dparam("TOLSTAT");
      for (int i = 0; i < omega->mstat->n; i++)
	if (   (fabs(omega->mapx->d[0] - omega->mstat->x[i]) < tolstat)
	    && (fabs(omega->mapy->d[0] - omega->mstat->y[i]) < tolstat)
	    && (omega->mapx->t[0]*omega->mstat->sgnx[i] >= 0.)
	    && (omega->mapy->t[0]*omega->mstat->sgny[i] >= 0.)) {
	  // Get the unmapped velocity
	  if (omega->mapx->t[0] != 0.)
	    dsadd (nwork, omega->mapx->t[0], u[0], 1, u[0], 1);
	  if (omega->mapy->t[0] != 0.)
	    dsadd (nwork, omega->mapy->t[0], u[1], 1, u[1], 1);
	  dvvtvp(nwork, u[0], 1, u[0], 1, u2nd[0], 1, u2nd[0], 1);
	  dvvtvp(nwork, u[1], 1, u[1], 1, u2nd[1], 1, u2nd[1], 1);
	  dvvtvp(nwork, u[2], 1, u[2], 1, u2nd[2], 1, u2nd[2], 1);
	  dvvtvp(nwork, u[0], 1, u[1], 1, u2nd[3], 1, u2nd[3], 1);
	  omega->mstat->nstat[i]++;
	}
    } else
#endif
    dvvtvp(nwork, u[0], 1, u[0], 1, u2nd[0], 1, u2nd[0], 1);
    dvvtvp(nwork, u[1], 1, u[1], 1, u2nd[1], 1, u2nd[1], 1);
    dvvtvp(nwork, u[2], 1, u[2], 1, u2nd[2], 1, u2nd[2], 1);
    dvvtvp(nwork, u[0], 1, u[1], 1, u2nd[3], 1, u2nd[3], 1);
  }

  FFT_AtoA (Qx, res[0]);   /* Transform back */
  FFT_AtoA (Qy, res[1]);
  FFT_AtoA (Qz, res[2]);

  Qx->Set_state('p');
  Qy->Set_state('p');
  Qz->Set_state('p');
  
  free_dmatrix(u,0,0);
  free_dmatrix(q,0,0);
  free_dmatrix(res,0,0);
  return;
}
