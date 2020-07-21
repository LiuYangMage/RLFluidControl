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
#include "varv.h"

static double *grad_wk=(double*)0, *workspace_end;
extern double **u2nd;
extern double **u2nde;

void VdgradV_VarV(Domain *omega)
{
  Element_List *U  =  omega->U,  *V  =  omega->V,  *W  = omega->W,
               *Ux =  omega->P,  *Uy = omega->Uf, *Uz = omega->Uf,
               *Vx =  omega->P,  *Vy = omega->Vf, *Vz = omega->Vf,
               *Wx =  omega->P,  *Wy = omega->Wf, *Wz = omega->Wf,
               *nul= (Element_List*)0, *Div = omega->P;
  int           nq = U->htot*U->nz;
  char        trip = 'a';
  varv_variables *vvexp;
  static Nek_Trans_Type f_to_p = F_to_P,
                        p_to_f = P_to_F;

  static int step = 1, hisstep = option("hisstep");
  static double time = dparam("STARTIME");

  if(!grad_wk){
    int worksize = nq;
    if (option("dealias")) {
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
      worksize = 3*nq/2;
    }
    grad_wk = dvector(0,worksize+3*nq-1);
    workspace_end = grad_wk + worksize;
    //    dzero(worksize,grad_wk,1);
  }

  dcopy(nq, U->base_h, 1, workspace_end     , 1);
  dcopy(nq, V->base_h, 1, workspace_end+  nq, 1);
  dcopy(nq, W->base_h, 1, workspace_end+2*nq, 1);
  
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
//    ROOT printf ("Time %lf Momentum Mx = %#lf My = %#lf Mz = %#lf Mass flux = %#lf L2 Divergence = %#lf\n", time, VolInt (U), VolInt (V), VolInt (W), VolInt (Div), GDiv);
    Div->Set_state('t');
    time += hisstep*dparam("DELT");
  } else {
    // Take z derivatives in fourier/quadrature space
    U->Grad(nul,nul,Uz,'z'); Uz->Set_state('p');
    V->Grad(nul,nul,Vz,'z'); Vz->Set_state('p');
    W->Grad(nul,nul,Wz,'z'); Wz->Set_state('p');
  }

  // Transform velocity to physical/quadrature space
  U->Trans(U, f_to_p); 
  V->Trans(V, f_to_p); 
  W->Trans(W, f_to_p); 

  // Transform z derivative of velocity to physical/quadrature space
  Uz->Trans(Uz, f_to_p);
  Vz->Trans(Vz, f_to_p);
  Wz->Trans(Wz, f_to_p);

  if (option("dealias")) {
    trip = 'A';
    nq   = 3*nq/2;
  }

  if(option("STATAVG")){
    dvvtvp(nq, U->base_h, 1, U->base_h, 1, u2nd[0], 1, u2nd[0], 1);
    dvvtvp(nq, V->base_h, 1, V->base_h, 1, u2nd[1], 1, u2nd[1], 1);
    dvvtvp(nq, W->base_h, 1, W->base_h, 1, u2nd[2], 1, u2nd[2], 1);
    dvvtvp(nq, U->base_h, 1, V->base_h, 1, u2nd[3], 1, u2nd[3], 1);

  }
  if (option("VARV"))
    vvexp = CsCalc(omega);
  if (iparam("ISO")!= 0)
    IsoStats(omega);
  // Uf = U.Grad u
  dvmul  (nq, W->base_h,1,Uz->base_h,1,grad_wk,1);
  U->Grad_h(U->base_h, Ux->base_h,Uy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Ux->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Uy->base_h,1,grad_wk,1,Uz->base_h,1);
  dscal  (nq,-1.0,Uz->base_h,1);
  if (option("VARV")) dvadd(nq,Uz->base_h,1,vvexp->NLx,1,Uz->base_h,1);

  // Vf = U.Grad v
  dvmul  (nq, W->base_h,1,Vz->base_h,1,grad_wk,1);
  V->Grad_h(V->base_h,Vx->base_h,Vy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Vx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Vy->base_h,1,grad_wk,1,Vz->base_h,1);
  dscal  (nq,-1.0,Vz->base_h,1);
  if (option("VARV")) dvadd(nq,Vz->base_h,1,vvexp->NLy,1,Vz->base_h,1);

  // Wf = U.Grad w
  dvmul  (nq, W->base_h,1,Wz->base_h,1,grad_wk,1);
  W->Grad_h(W->base_h, Wx->base_h,Wy->base_h,(double *)NULL,trip);
  dvvtvp (nq, U->base_h,1,Wx->base_h,1,grad_wk,1,grad_wk,1);
  dvvtvp (nq, V->base_h,1,Wy->base_h,1,grad_wk,1,Wz->base_h,1);
  dscal  (nq,-1.0,Wz->base_h,1);
  if (option("VARV")) dvadd(nq,Wz->base_h,1,vvexp->NLz,1,Wz->base_h,1);
  
  // Transform U.Grad U to Fourier/Quadrature space
  Uz->Trans(Uz, p_to_f);
  Vz->Trans(Vz, p_to_f);
  Wz->Trans(Wz, p_to_f);

  if (option("dealias"))
    nq = U->htot*U->nz;

  dcopy(nq, workspace_end     , 1, U->base_h, 1);
  dcopy(nq, workspace_end+  nq, 1, V->base_h, 1);
  dcopy(nq, workspace_end+2*nq, 1, W->base_h, 1);

  if (option("VARV"))
  {
	  free(vvexp->NLx);free(vvexp->NLy);free(vvexp->NLz);
  }

  step == hisstep ? step = 1 : step++;

  return;
}


void IsoStats(Domain *omega){
  register int n,j,k,i, nel,nq,nz;
  double  *dudx , *dudy , *dvdx , *dvdy,*dudz,*dvdz,*dwdz,*dwdx,*dwdy;
  double **tmp,*wk;
  Element_List  *U  = omega->U,  *V  =  omega->V,  *W  =  omega->W,
    *Uz = omega->Uf, *Vz = omega->Vf, *Wz = omega->Wf;
  double divave=0.0, uave=0.0,vave=0.0,wave=0.0,tkeave=0.0,dudx2=0.0;
  double dudx3=0.0,dudy2=0.0,dudy3=0.0,dudz2=0.0,dudz3=0.0;
  double dvdx2=0.0,dvdx3=0.0,dvdy2=0.0,dvdy3=0.0,dvdz2=0.0,dvdz3=0.0;
  double dwdx2=0.0,dwdx3=0.0,dwdy2=0.0,dwdy3=0.0,dwdz2=0.0,dwdz3=0.0;
  double con,rsmall = 1E-15,u2ave=0.0,v2ave=0.0,w2ave=0.0,
  Relx=0.0,Rely=0.0,Relz=0.0;
  double w;
  double skewux=0.0,skewuy=0.0,skewuz=0.0,skewvx=0.0,skewvy=0.0;
  double skewvz=0.0,skewwx=0.0,skewwy=0.0,skewwz=0.0;
  FILE *fp=omega->his_file;

  nz       = option("NZ");

  nq = U->htot*U->nz;
  if (option("dealias")) 
    nq = 3*U->htot*U->nz/2;
  tmp = dmatrix(0,8,0,nq);
  wk  = dvector(0,nq);
  dudx= tmp[0]; 
  dudy= tmp[1];
  dvdx= tmp[2];
  dvdy= tmp[3];
  dudz= tmp[4];
  dvdz= tmp[5];
  dwdz= tmp[6];
  dwdy= tmp[7];
  dwdx= tmp[8];

  /* Zeroing the values */
  divave=0.0;uave=0.0;vave=0.0;wave=0.0;tkeave=0.0;dudx2=0.0;
  dudx3=0.0;dudy2=0.0;dudy3=0.0;dudz2=0.0;dudz3=0.0;
  dvdx2=0.0;dvdx3=0.0;dvdy2=0.0;dvdy3=0.0;dvdz2=0.0;dvdz3=0.0;
  dwdx2=0.0;dwdx3=0.0;dwdy2=0.0;dwdy3=0.0;dwdz2=0.0;dwdz3=0.0;
  con=0.0;rsmall = 1E-15;
  skewux=0.0,skewuy=0.0,skewuz=0.0,skewvx=0.0,skewvy=0.0;
  skewvz=0.0,skewwx=0.0,skewwy=0.0,skewwz=0.0;

 
  U->Grad_h(U->base_h,dudx,dudy,(double *)NULL,'a');
  V->Grad_h(V->base_h,dvdx,dvdy,(double *)NULL,'a');
  W->Grad_h(W->base_h,dwdx,dwdy,(double *)NULL,'a');
  for (i=0;i<nq;i++)
    {
      
      uave = uave + U->base_h[i];
      vave = vave + V->base_h[i];
      wave = wave + W->base_h[i];
      u2ave = u2ave + U->base_h[i]*U->base_h[i];
      v2ave = v2ave + V->base_h[i]*V->base_h[i];
      w2ave = w2ave + W->base_h[i]*W->base_h[i];
      
      tkeave = tkeave +  u2ave + v2ave + w2ave ;
      
      dudx2 = dudx2 + dudx[i]*dudx[i];
      dudx3 = dudx3 + dudx[i]*dudx[i]*dudx[i];
      dudy2 = dudy2 + dudy[i]*dudy[i];
      dudy3 = dudy3 + dudy[i]*dudy[i]*dudy[i];
      dudz2 = dudz2 + Uz->base_h[i] * Uz->base_h[i] ;
      dudz3 = dudz3 +  Uz->base_h[i] * Uz->base_h[i] *Uz->base_h[i];

      dvdx2 = dvdx2 + dvdx[i]*dvdx[i];
      dvdx3 = dvdx3 + dvdx[i]*dvdx[i]*dvdx[i];
      dvdy2 = dvdy2 + dvdy[i]*dvdy[i];
      dvdy3 = dvdy3 + dvdy[i]*dvdy[i]*dvdy[i];
      dvdz2 = dvdz2 + Vz->base_h[i] * Vz->base_h[i];
      dvdz3 = dvdz3 + Vz->base_h[i] * Vz->base_h[i] *Vz->base_h[i];
      
      dwdx2 = dwdx2 + dwdx[i]*dwdx[i];
      dwdx3 = dwdx3 + dwdx[i]*dwdx[i]*dwdx[i];
      dwdy2 = dwdy2 + dwdy[i]*dwdy[i];
      dwdy3 = dwdy3 + dwdy[i]*dwdy[i]*dwdy[i];
      dwdz2 = dwdz2 + Wz->base_h[i] * Wz->base_h[i];
      dwdz3 = dwdz3 + Wz->base_h[i] * Wz->base_h[i] *Wz->base_h[i];
      
      divave = divave + dudx[i] + dvdy[i] + Wz->base_h[i];
      
    }
  
  con = 1.0/(nq);
  divave = divave * con;
  
  uave = uave * con;
  vave = vave * con;
  wave = wave * con;
  
  u2ave = u2ave * con;
  v2ave = v2ave * con;
  w2ave = w2ave * con;

  tkeave = 0.5 * tkeave * con;

  dudx2 = dudx2 * con;
  dudx3 = dudx3 * con;
  dudy2 = dudy2 * con;
  dudy3 = dudy3 * con;
  dudz2 = dudz2 * con;
  dudz3 = dudz3 * con;

  dvdx2 = dvdx2 * con;
  dvdx3 = dvdx3 * con;
  dvdy2 = dvdy2 * con;
  dvdy3 = dvdy3 * con;
  dvdz2 = dvdz2 * con;
  dvdz3 = dvdz3 * con;

  dwdx2 = dwdx2 * con;
  dwdx3 = dwdx3 * con;
  dwdy2 = dwdy2 * con;
  dwdy3 = dwdy3 * con;
  dwdz2 = dwdz2 * con;
  dwdz3 = dwdz3 * con;

  Relx = (sqrt(u2ave) /dparam("KINVIS-ORIG")) * sqrt(u2ave/dudx2);
  Rely = (sqrt(v2ave) /dparam("KINVIS-ORIG")) * sqrt(v2ave/dvdy2);
  Relz = (sqrt(w2ave) /dparam("KINVIS-ORIG")) * sqrt(w2ave/dwdz2);

  skewux = dudx3 / sqrt((pow(dudx2,3)));
  skewuy = dudy3 / sqrt((pow(dudy2,3)));
  skewuz = dudz3 / sqrt((pow(dudz2,3)));

  skewvx = dvdx3 / sqrt((pow(dvdx2,3)));
  skewvy = dvdy3 / sqrt((pow(dvdy2,3)));
  skewvz = dvdz3 / sqrt((pow(dvdz2,3)));

  skewwx = dwdx3 / sqrt((pow(dwdx2,3)));
  skewwy = dwdy3 / sqrt((pow(dwdy2,3)));
  skewwz = dwdz3 / sqrt((pow(dwdz2,3)));

  gdsum(&Relx,1,&w);
  gdsum(&tkeave,1,&w);
  gdsum(&Rely,1,&w);
  gdsum(&Relz,1,&w);
  Relx   = Relx   / option("NPROCS");
  Rely   = Rely   / option("NPROCS");
  Relz   = Relz   / option("NPROCS");
  tkeave = tkeave / option("NPROCS");
  
  ROOTONLY printf("%e  ",Relx);
  ROOTONLY printf("%e  ",Rely);
  ROOTONLY printf("%e  ",Relz);
  ROOTONLY printf("%e  :Rex Rey Rez tke\n",tkeave);
  
  free_dmatrix(tmp,0,0); free(wk);
}
