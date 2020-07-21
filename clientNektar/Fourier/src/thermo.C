/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <veclib.h>
#include <time.h>
#include "nektarF.h"
#include "thermo2d.h"
int alternating_switch = 0;

void MakeF_T(Domain *omega, ACTION act)
{
#ifdef THERMO
  int j;
  switch(act)
    {
    case Prep:{ /* put field in physical space for waveprop */
      Jtransbwd(omega->T,omega->Ts);
      break;
    }
    case Rotational:{
      break;
    }
    case Convective:{
      VdgradT (omega);
      goto AddForcing_T;
    }
    case Stokes:{
      break;
    }
    case Alternating:{
      switch (alternating_switch){
      case 0:{
	VdgradT (omega);
	alternating_switch = 1;
	break;
      }
      case 1:{
	divTV (omega);
	alternating_switch = 0;
	break;
      }
      }
	goto AddForcing_T;
    }
      
    AddForcing_T:{
      double f;
      Element *Tf   = *omega->Tf;
      int ntot = Tf[0].qa*Tf[0].qb;
      int k;
      int nel = countelements(Tf);
      /* Standard Forcing on constant plane */
      
      if ((f = dparam("FFTX")) != 0.0) 
	for(k = 0; k < nel; ++k){
	  ntot = Tf[k].qa*Tf[k].qb;
	  dsadd(ntot, f, *Tf[k].h, 1, *Tf[k].h, 1);
	}
      break;
    }
      
    case Viscous:{
      double dt = dparam("DELT");      
      double Pe = dparam("PECLET");
      double   PeDt  = Pe / dt;
      register int i,j,k;
      Element *T, *Tf, *Ts;
      int nel = countelements(omega->T);
      T    =  omega->T ;
      Ts   =  omega->Ts;
      Tf   = *omega->Tf;
      
      InnerProduct(Tf,Tf);
      
      for(k = 0; k < nel; ++k){
	for(i = 0; i < Tf[k].Nmodes; ++i)
	  Tf[k].vert->hj[i]= -PeDt*Tf[k].vert->hj[i];
	
	/* Put initial guess as solution from last time step */
	dcopy(T[k].Nmodes, Ts[k].vert->hj, 1, T[k].vert->hj, 1);
      }
      break;
    }

    case Post:{
      int    Nm,k;
      double theta = dparam("THETA");
      int    nel = countelements(omega->T);
      double fac   = 1/(1.0 - theta);
      for(k = 0; k < nel; ++k){
	Nm = omega->Ts[k].Nmodes;
	dsvtvp(Nm, -theta, omega->Ts[k].vert->hj, 1,
	       omega->T[k].vert->hj, 1, omega->Ts[k].vert->hj, 1);
	
	dsmul(Nm, fac, omega->Ts[k].vert->hj, 1, omega->Ts[k].vert->hj, 1);
      }	
      set_state(omega->Ts,'t');
      break;
    }
    default:
      error_msg(MakeF--unknown type of action);
      break;
    }
  
#endif
  return;
}

void VdgradT(Domain *omega){
#ifdef THERMO
  Element_List *U  =  omega->U,  *V  =  omega->V,  *W  = omega->W,
               *T  =  omega->T,  *Tf =  omega->Tf;
  register int i,j,k;
  
  int          nq = U->htot*U->nz;
  char         trip = 'a';
  static Nek_Trans_Type f_to_p = F_to_P,
                        p_to_f = P_to_F;
  
  if (option("dealias")) {
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
      worksize = 3*nq/2;
   }
  
  
#endif
  return;
}


void divTV(Domain *omega){
#ifdef THERMO
  Element *U, *V, *T, *Tz;
  register int i,j,k;
  double  *dx, *dy, *tu, *wk, **tmp  = dmatrix(0,3,0,QGmax*QGmax-1);
  int     nel = countelements(omega->T), nq;
  
  dx = tmp[0];  dy = tmp[1];  tu = tmp[2]; wk = tmp[3];

  
  U  =  omega->U,  V  =  omega->V, 
  T  =  omega->T,  Tz = *omega->Tf;
  
  /* Calculate Div(TV) */
  for(k = 0; k < nel; ++k){
    nq = T[k].qa*T[k].qb;
    if(T[k].qa != U[k].qa || T[k].qb != U[k].qb){
      EJtransbwd(Tz+k, U+k);
      dvmul(nq, Tz[k].h[0], 1, T[k].h[0], 1, tu, 1);
      Egradv(T+k, tu, dx, (double*)NULL, 'x',wk);
      
      EJtransbwd(Tz+k, V+k);
      dvmul(nq, Tz[k].h[0], 1, T[k].h[0], 1, tu, 1);
      Egradv(T+k, tu, (double*)NULL, dy, 'y',wk);
      
      dsmul(nq,  -1.0, dx, 1, Tz[k].h[0], 1);
      dsvtvp(nq, -1.0, dy, 1, Tz[k].h[0], 1, Tz[k].h[0], 1);
    }
    else{
      dvmul(nq, U[k].h[0], 1, T[k].h[0], 1, tu, 1);
      Egradv(T+k, tu, dx, (double*)NULL, 'x',wk);
      dsmul(nq,  -1.0, dx, 1, Tz[k].h[0], 1);
      
      dvmul(nq, V[k].h[0], 1, T[k].h[0], 1, tu, 1);
      Egradv(T+k, tu, (double*)NULL, dx, 'y',wk);
      dsvtvp(nq, -1.0, dx, 1, Tz[k].h[0], 1, Tz[k].h[0], 1);
    }
  }
  set_state(Tz,'p'); 
  
  free_dmatrix(tmp,0,0);
#endif
  return;
}









