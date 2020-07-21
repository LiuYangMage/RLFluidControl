
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include "nektarF.h"

#include "nse_util.h"


static double *escales1 = (double*)0;
static double *escales2 = (double*)0;
static double q2_variation = 0.;

void compute_velocity_range(Domain *omega){
    register int i, j, k, plane;
    int          qa,qb,qc,m;
    Element_List *Ut,*Vt,*Wt;
    Element      *E,*F,*G;
    int          nz = option("NZ"), nztot = option("NZTOT");
    int          nprocs = option("NPROCS");
    int          eDIM = 3;

    Ut = omega->Ut;
    Vt = omega->Vt;
    Wt = omega->Wt;

    int NZ       = Ut->nz;
    int htot     = Ut->nz*Ut->htot;
    if (option("dealias")) {
      NZ   = 3*NZ/2;
      htot = 3*htot/2;
    }

    double umax = -1e6,umax_tmp = -1e6;
    double uu_max = -1e6,uu_min = 1e6;
    double vv_max = -1e6,vv_min = 1e6;
    double ww_max = -1e6,ww_min = 1e6;
    double uu,vv,ww;

    dcopy(htot,omega->uk[1],1,Ut->base_h,1);
    dcopy(htot,omega->vk[1],1,Vt->base_h,1);
    dcopy(htot,omega->wk[1],1,Wt->base_h,1);

    for (k = 0; k < NZ; k++)
//     for(E=Ut->flevels[k]->fhead,F=Vt->flevels[k]->fhead,G=Wt->flevels[k]->fhead;E;E=E->next,F=F->next,G=G->next)
    for(int eid =0; eid<Ut->nel; ++eid)
      {
	      E   = Ut->flevels[0]->flist[eid];
        qa = E->qa;   qb = E->qb;

       int offset = Ut->htot*k+eid*E->qtot; 

        for ( i=0; i<qa; i++ )
        for ( j=0; j<qb; j++ )
        {    
          m = i*qb + j;

          uu = Ut->base_h[offset+m];
          vv = Vt->base_h[offset+m];
          ww = Wt->base_h[offset+m];

          umax_tmp  = uu*uu+vv*vv+ww*ww;

          umax = max(umax,umax_tmp);

          uu_max = max(uu,uu_max);
          uu_min = min(uu,uu_min);

          vv_max = max(vv,vv_max);
          vv_min = min(vv,vv_min);

          ww_max = max(ww,ww_max);
          ww_min = min(ww,ww_min);

        }
    }

//    gdmax(&umax,  1,  &umax_tmp);
    MPI_Allreduce (&umax, &umax_tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&umax_tmp,1,&umax,1);

    MPI_Allreduce (&uu_max, &umax_tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&umax_tmp,1,&uu_max,1);
//    gdmax(&uu_max,  1,  &umax_tmp);
    uu_min = -uu_min; 
    MPI_Allreduce (&uu_min, &umax_tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&umax_tmp,1,&uu_min,1);
//    gdmax(&uu_min,  1,  &umax_tmp);

    MPI_Allreduce (&vv_max, &umax_tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&umax_tmp,1,&vv_max,1);
//    gdmax(&vv_max,  1,  &umax_tmp);
    vv_min = -vv_min; 
    MPI_Allreduce (&vv_min, &umax_tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&umax_tmp,1,&vv_min,1);
//    gdmax(&vv_min,  1,  &umax_tmp);

    MPI_Allreduce (&ww_max, &umax_tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&umax_tmp,1,&ww_max,1);
//    gdmax(&ww_max,  1,  &umax_tmp);
    ww_min = -ww_min; 
    MPI_Allreduce (&ww_min, &umax_tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&umax_tmp,1,&ww_min,1);
//    gdmax(&ww_min,  1,  &umax_tmp);
    

    omega->velocity_range[0][0] = uu_max;
    omega->velocity_range[0][1] = -uu_min;

    omega->velocity_range[1][0] = vv_max;
    omega->velocity_range[1][1] = -vv_min;

    omega->velocity_range[2][0] = ww_max;
    omega->velocity_range[2][1] = -ww_min;

    omega->maximal_velocity = sqrt(umax);

    if(umax != umax) {
      fprintf(stderr,"wrong value ...\n");
      exit(-1);
     }


}

#ifdef ENTROPYVISCOSITY

void compute_velocity_entropy_variation (Domain *omega){
    register int i,j,k;
    int          qa,qb,qc,m;
    Element_List *Ut,*Vt,*Wt;
    Element      *E,*F,*G;
    int          eDIM = 3;
  	double       *za,*wa,*zb,*wb,*zc,*wc;
    double       weight;
    int          nprocs = option("NPROCS");
    int          nz = option("NZ"), nztot = option("NZTOT");
    int          htot = nz*omega->Ut->htot;
    int          htot_nz = nz*omega->Ut->htot;

    if (option("dealias")) {
      nz   = 3*nz/2;
      htot = 3*htot/2;
    }

    double       **u_range = omega->velocity_range;
    double       max_velocity = omega->maximal_velocity;
    double       uu,vv,ww;
    double       velocity_mag;
    double       entropy;

    double umean = 0.5*(u_range[0][0]+u_range[0][1]);
    double vmean = 0.5*(u_range[1][0]+u_range[1][1]);
    double wmean = 0.5*(u_range[2][0]+u_range[2][1]);

//    double       *up1  = omega->up[1],
//                 *vp1  = omega->vp[1],
//                 *wp1  = omega->wp[1];

    Ut = omega->Ut;
    Vt = omega->Vt;
    Wt = omega->Wt;

    dcopy(htot,omega->uk[1],1,Ut->base_h,1);
    dcopy(htot,omega->vk[1],1,Vt->base_h,1);
    dcopy(htot,omega->wk[1],1,Wt->base_h,1);
//    dzero(htot,omega->up[1],1);
//    dzero(htot,omega->vp[1],1);
//    dzero(htot,omega->wp[1],1);

    double volume = 0.,integ_entropy = 0.;
    double max_entropy = -1e6, min_entropy = 1e6;
    for (k = 0; k < nz; k++)
//    for(E=Ut->flevels[k]->fhead,F=Vt->flevels[k]->fhead,G=Wt->flevels[k]->fhead;E;E=E->next,F=F->next,G=G->next)
    for(int eid =0; eid<Ut->nel; ++eid)
    {
	      E   = Ut->flevels[0]->flist[eid];
	      F   = Vt->flevels[0]->flist[eid];
	      G   = Wt->flevels[0]->flist[eid];
        qa = E->qa;   qb = E->qb;
        E->GetZW(&za, &wa, &zb, &wb, &zc, &wc);

        int offset = k*Ut->htot+eid*E->qtot; 

        for ( i=0; i<qa; i++ )
        for ( j=0; j<qb; j++ )
        {    
          m = i*qb + j;
          
          if( E->curvX )
            weight = wb[j]*wa[i]*E->geom->jac.p[m];
          else
            weight = wb[j]*wa[i]*E->geom->jac.d;

//          uu = E->h[0][m];
//          vv = F->h[0][m];
//          ww = G->h[0][m];
          uu = Ut->base_h[offset+m];
          vv = Vt->base_h[offset+m];
          ww = Wt->base_h[offset+m];

          velocity_mag = sqrt(uu*uu+vv*vv+ww*ww);
          
//          entropy = (uu-umean)*(uu-umean)+
//                    (vv-vmean)*(vv-vmean)+
//                    (ww-wmean)*(ww-wmean);
          entropy = (velocity_mag-0.5*max_velocity)
                   *(velocity_mag-0.5*max_velocity);

          volume += weight;
          integ_entropy += entropy*weight;

          max_entropy = max(max_entropy,entropy);
          min_entropy = min(min_entropy,entropy);

//          up1[offset+m] += uu;
//          vp1[offset+m] += vv;
//          wp1[offset+m] += ww;
        }
    }

    double tmp = 0.;
    
    MPI_Allreduce (&volume, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(1,&tmp,1,&volume,1);
    
    MPI_Allreduce (&integ_entropy, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(1,&tmp,1,&integ_entropy,1);
//    gdsum(&volume,1,&tmp);
//    gdsum(&integ_entropy, 1, &tmp);

    double mean_entropy = integ_entropy/volume;

    MPI_Allreduce (&max_entropy, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&tmp,1,&max_entropy,1);
//    gdmax(&max_entropy,1,&tmp);

    min_entropy = -min_entropy;
    MPI_Allreduce (&min_entropy, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&tmp,1,&min_entropy,1);
//    gdmax(&min_entropy,1,&tmp);
    
    omega->entropy_variation = max(max_entropy-mean_entropy,
                                   mean_entropy-(-min_entropy)); 
    
//    MPI_Allreduce (up1, &up1[Ut->htot], Ut->htot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    dcopy(Ut->htot,&up1[Ut->htot],1,up1,1);

//    MPI_Allreduce (vp1, &vp1[Ut->htot], Ut->htot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    dcopy(Ut->htot,&vp1[Ut->htot],1,vp1,1);

//    MPI_Allreduce (wp1, &wp1[Ut->htot], Ut->htot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    dcopy(Vt->htot,&wp1[Ut->htot],1,wp1,1);

//    int int_tmp;
//
//    MPI_Allreduce (&nz, &int_tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

//    for (i=0; i<Ut->htot; ++i) {
//      up1[i] /=(double) nztot;
//       vp1[i] /=(double) nztot;
//       wp1[i] /=(double) nztot;
//    }
}


static double setQfactor(int ind, int limit, int low, int high);
void SetReducedPhys(Element_List *In, Element_List *Out, Element_List *Work);

void compute_entropy_viscosity(Domain *omega){
    register int i,j,q;
    int          qa,qb,qc,m,n,L;
    double       ns_residual[3];
  	double       *null= (double*)0;
	  double       *za,*wa,*zb,*wb,*zc,*wc;
	  double       dt    = dparam("DELT");
	  double       dtinv = 1./dt;
    int          id0,id1;
    double        medge,elen,mean_edge;
    Element      *E,*Et,*Ef,*F,*Ft,*Ff,*G,*Gt,*Gf,*St,*Rt;
    double       dxa, dxb, sca, scb, scc, scd;
    
    int          nz = option("NZ"), nztot = option("NZTOT");
    int          nprocs  = option("NPROCS");
    int          myid    = option("PROCID");
    int          hjtot   = omega->U->nz*omega->U->hjtot;
    int          htot    = omega->U->nz*omega->U->htot;
    int          htot_nz = omega->U->nz*omega->U->htot;

    static       Nek_Trans_Type f_to_p = F_to_P,
                                p_to_f = P_to_F;
    char         trip = 'a';
    
    if (option("dealias")) {
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
      nz   = 3*nz/2;
      htot = 3*htot/2;
      trip = 'A';
    }

	  double       Re     = 1./dparam("KINVIS");
    double       niu    = dparam("KINVIS");
    double       CA    = dparam("STABILIZATION_ALPHA");
    double       CB    = dparam("STABILIZATION_BETA");
    double       WALE_A = dparam("DWALE_ALPHA")?dparam("DWALE_ALPHA"):0.325;
    double       DSIJ   = dparam("DSIJ_SCALE")? dparam("DSIJ_SCALE"):1.0;
    double       my_radius    = dparam("RADIUS")? dparam("RADIUS"):0.5;
    double       cutoff_accel    = dparam("CUT_OFF_ACCEL")? dparam("CUT_OFF_ACCEL"):0.2;
    double       cutoff_dist    = dparam("CUT_OFF_DIST")? dparam("CUT_OFF_DIST"):2.0;
    double       alpha_reduction = dparam("EVM_ALPHA_REDUCTION")? dparam("EVM_ALPHA_REDUCTION"):0.01;
    int          use_reduction = iparam("IBL_VISC_REDUCTION");

    int          iostep = iparam("IOSTEP");

    double       **u_range = omega->velocity_range;
    double       u_mean = 0.5*(u_range[0][0]+u_range[0][1]);
    double       v_mean = 0.5*(u_range[1][0]+u_range[1][1]);
    double       w_mean = 0.5*(u_range[2][0]+u_range[2][1]);

    double       min_alpha;
    double       alpha_scale;

    double       qscale = 1.;
    
    double entropy_residual;
    double entropy_viscosity,max_viscosity = -1e6,element_max_viscosity;
    double element_diameter,qll_diameter;
    double entropy_variation        = omega->entropy_variation;
    double first_order_viscosity    = 0.;
    double max_velocity             = omega->maximal_velocity;
    double A                        = dparam("VAN_DRIEST_A")? dparam("VAN_DRIEST_A"):26.;
    double boundary_layer_thickness = 0.37*pow(Re,-0.2);

    double max_sensor = 0.0;

#ifdef MAP
     Mapping *mapx = omega->mapx,
             *mapy = omega->mapy;

//     realft (NZ/2, mapx->t, 1); //Transform to Physical space
//     realft (NZ/2, mapx->tt,1);
//     realft (NZ/2, mapy->t, 1);
//     realft (NZ/2, mapy->tt,1);

     double map_vx = mapx->t[0];
     double map_vy = mapy->t[0];

     double map_ax = mapx->tt[0];
     double map_ay = mapy->tt[0];
     
//     realft (NZ/2, mapx->t, -1); //Transform to Fourier space
//     realft (NZ/2, mapx->tt,-1);
//     realft (NZ/2, mapy->t, -1);
//     realft (NZ/2, mapy->tt,-1);
//     double effect_Re = sqrt((-1.0+map_vx)*(-1.0+map_vx)+map_vy*map_vy)*Re;
     
//     ROOTONLY fprintf(stderr,"Effective Reynolds number: %lf  \n",effect_Re);
#endif
    if(dparam("VAN_DRIEST")&&(!dparam("RE_TAU"))&&(!(dparam("CYLINDER"))))
     {
       fprintf(stderr,"param RE_TAU is ZERO \n");
       exit(1);
     }

    double damping_length;
    double AA;
    static int step = 1;
    double Re_tau = dparam("RE_TAU");
    
//    double actived_mode = (dparam("ACTIVE_MODE"))? dparam("ACTIVE_MODE")/2. : 5.*sqrt((double) (nztot/2));
//    double *smoothing_factor = dvector(0, nz-1);
//    dzero(nz,smoothing_factor,1);
//
//
//   for(int k = 0, j = 0; k < nz; k=k+2, ++j)
//     {
//       if(parid(k)  > actived_mode)
//          {
//             double fac = (double)((double)(parid(k) - nztot)*(parid(k) - nztot))/
//                        ((double)((parid(k) - actived_mode)*(parid(k) - actived_mode)));
//             smoothing_factor[k]   = exp(-fac);
//             smoothing_factor[k+1] = exp(-fac);
//          }
//
//    //   fprintf(stdout," k = %d  fac = %lf \n", parid(k),smoothing_factor[k]);
//    }

//    double  **d = dmatrix(0,8,0,QGmax*QGmax*QGmax-1);
//    double  *Ux=d[0], *Uy=d[1], *Uz=d[2];
//    double  *Vx=d[3], *Vy=d[4], *Vz=d[5];
//    double  *Wx=d[6], *Wy=d[7], *Wz=d[8];

    if((CA==0)||(CB==0.))
     {
      ROOTONLY fprintf(stderr,"Forget to set STABILIZATION_ALPHA or STABILIZATION_BETA parameter ! \n");
      exit(1);
     }


    Element_List *U    =  omega->U,  *V    =  omega->V,   *W  = omega->W,
                 *Uf   =  omega->Uf, *Vf   =  omega->Vf,  *Wf = omega->Wf,
                 *Ut   =  omega->Ut, *Vt   =  omega->Vt,  *Wt = omega->Wt,
                 *P    =  omega->P,  *Pf   =  omega->Pf,  *Pt = omega->Pt;
#ifdef SPM
    Element_List *CONCENTR = omega->CONCENTR; 
#endif

    double       *u0   =  omega->uk[0],     *v0    =  omega->vk[0],     *w0   = omega->wk[0], // velocity value at time step n-1
                 *u1   =  omega->uk[1],     *v1    =  omega->vk[1],     *w1   = omega->wk[1], // velocity value at time step n
                 *p0   =  omega->pk[0],     *p1    =  omega->pk[1];

//    double       *uf0   =  omega->uf[0],    *vf0    =  omega->vf[0],    *wf0  = omega->wf[0], // convective term at time step n-1
//                 *uf1   =  omega->uf[1],    *vf1    =  omega->vf[1],    *wf1  = omega->wf[1]; // convective term at time step n

    double       *ut0  =  omega->ut[0],    *vt0   =  omega->vt[0],    *wt0 = omega->wt[0], // used as temps arrary for gradients
                 *ut1  =  omega->ut[1],    *vt1   =  omega->vt[1],    *wt1 = omega->wt[1]; // used for calculating laplacians

    double       *up0   =  omega->up[0],    *vp0   =  omega->vp[0],     *wp0  = omega->wp[0], // stored pressure graidients term
                 *up1   =  omega->up[1],    *vp1   =  omega->vp[1],     *wp1  = omega->wp[1]; // stored pressure graidients term

    double       *us0   =  omega->uss[0],   *vs0   =  omega->vss[0],   *ws0  = omega->wss[0], // backup arrary
                 *us1   =  omega->uss[1],   *vs1   =  omega->vss[1],   *ws1  = omega->wss[1];

//#ifdef WALE
   double       **st   = omega->st;
//#endif

    double       *afx0   = omega->afx[0],    *afy0   =  omega->afy[0], *afz0  = omega->afz[0]; //force at time step n
    double       *afx1   = omega->afx[1],    *afy1   =  omega->afy[1], *afz1  = omega->afz[1];

    double       *visc_tmp = omega->visc_tmp, *visc_ave = omega->visc_ave, *visc_bak = omega->visc_bak;
    double       *visc_max_indices = omega->visc_max_indices;      
#ifdef NATURAL_CONVECTION
    double       *tk0   =  omega->tk[1],     *tk1    =  omega->tk[0];
#endif

    double BBeta[3];
    getbeta(BBeta);

    Ut->Set_state('p');
    Vt->Set_state('p');
    Wt->Set_state('p');
    Pt->Set_state('p');

#ifdef NS_RESIDUAL
    dvsub  (htot, u1, 1, u0,  1, Ut->base_h, 1);
    dvsub  (htot, v1, 1, v0,  1, Vt->base_h, 1);
    dvsub  (htot, w1, 1, w0,  1, Wt->base_h, 1);

    dsmul (htot, dtinv, Ut->base_h,  1, Ut->base_h, 1); 
    dsmul (htot, dtinv, Vt->base_h,  1, Vt->base_h, 1); 
    dsmul (htot, dtinv, Wt->base_h,  1, Wt->base_h, 1); 
    
 //   memset(ut1, '\0', sizeof(double)*htot);
 //   memset(vt1, '\0', sizeof(double)*htot);
 //   memset(wt1, '\0', sizeof(double)*htot);

    dsmul (htot, 0.5, afx0, 1, ut1, 1);
    daxpy (htot, 0.5, afx1, 1, ut1, 1);
    dsmul (htot, 0.5, afy0, 1, vt1, 1);
    daxpy (htot, 0.5, afy1, 1, vt1, 1);
    dsmul (htot, 0.5, afz0, 1, wt1, 1);
    daxpy (htot, 0.5, afz1, 1, wt1, 1);

    dvsub  (htot, Ut->base_h, 1, ut1,  1, ut1, 1); //force term
    dvsub  (htot, Vt->base_h, 1, vt1,  1, vt1, 1);
    dvsub  (htot, Wt->base_h, 1, wt1,  1, wt1, 1);
//
    dsmul (htot, 0.5, p0, 1, Pt->base_h, 1);
    daxpy (htot, 0.5, p1, 1, Pt->base_h, 1);
//    Pt->Grad      (Ut,Vt,Wt,'a');
    Pt ->Grad_h(Pt->base_h,Ut->base_h,Vt->base_h,null,trip);
    dvadd(htot, ut1, 1, Ut->base_h,  1, ut1, 1);
    dvadd(htot, vt1, 1, Vt->base_h,  1, vt1, 1);

    Pt->Trans(Ut, p_to_f); //go to (Q,F)-space 
    Ut->Grad_z(Wt); Wt->Set_state('p'); //in (Q,F)-space
    Wt->Trans(Wt, f_to_p); //go to (Q,P)-space 

    dvadd(htot, wt1, 1, Wt->base_h,  1, wt1, 1);

// laplacians u
    dsmul (htot, 0.5, u0, 1, Ut->base_h, 1);
    daxpy (htot, 0.5, u1, 1, Ut->base_h, 1);
    dsmul (htot, 0.5, v0, 1, Vt->base_h, 1);
    daxpy (htot, 0.5, v1, 1, Vt->base_h, 1);
    dsmul (htot, 0.5, w0, 1, Wt->base_h, 1);
    daxpy (htot, 0.5, w1, 1, Wt->base_h, 1);

    Ut ->Grad_h(Ut->base_h,ut0,vt0,null,trip);

    dcopy(htot, ut0, 1, up0, 1); //for divergence
    
    dcopy(htot, ut0, 1, Pt->base_h, 1);
    Pt->Grad_h(Pt->base_h,ut0,null,null,trip);

    dcopy(htot, vt0, 1, Pt->base_h, 1);
    Pt->Grad_h(Pt->base_h,null,vt0,null,trip);

    Ut->Trans(Ut, p_to_f); //go to (Q,F)-space 
    Ut->Grad_z(Pt); Pt->Set_state('p'); //in (Q,F)-space
    Pt->Grad_z(U); U->Set_state('p');// in (Q,F)-space
    U->Trans(U, f_to_p); //go to (Q,P)-space 

    dvadd(htot, ut0, 1, Uf->base_h,  1, Uf->base_h, 1);
    dvadd(htot, vt0, 1, Uf->base_h,  1, Uf->base_h, 1);
    dsmul(htot, -niu, U->base_h,  1, Uf->base_h, 1); // laplacians on x direction

    Ut->Trans(Ut, f_to_p); //go to (Q,P)-space 
//
//
//
    Vt ->Grad_h(Vt->base_h, ut0,vt0,null,trip);

    dvadd(htot, vt0, 1, up0,  1, up0, 1); //divergence


    dcopy(htot, ut0, 1, Pt->base_h, 1);
    Pt->Grad_h(Pt->base_h,ut0,null,null,trip);

    dcopy(htot, vt0, 1, Pt->base_h, 1);
    Pt->Grad_h(Pt->base_h, null,vt0,null,trip);

    Vt->Trans(Vt, p_to_f); //go to (Q,F)-space 
    Vt->Grad_z(Pt);Pt->Set_state('p'); //in (Q,F)-space
    Pt->Grad_z(V); V->Set_state('p');// in (Q,F)-space
    V->Trans(V, f_to_p); //go to (Q,P)-space 

    dvadd(htot, ut0, 1, Vf->base_h,  1, Vf->base_h, 1);
    dvadd(htot, vt0, 1, Vf->base_h,  1, Vf->base_h, 1);
    dsmul(htot, -niu, V->base_h,  1, Vf->base_h, 1); // laplacians on y direction

    Vt->Trans(Vt, f_to_p); //go to (Q,P)-space 
//
//
//
    Wt ->Grad_h(Wt->base_h, ut0,vt0,null,trip);
    dcopy(htot, ut0, 1, Pt->base_h, 1);
    Pt->Grad_h(Pt->base_h, ut0,null,null,trip);

    dcopy(htot, vt0, 1, Pt->base_h, 1);
    Pt->Grad_h(Pt->base_h, null,vt0,null,trip);

    Wt->Trans(Wt, p_to_f); //go to (Q,F)-space 
    Wt->Grad_z(Pt); Pt->Set_state('p'); //in (Q,F)-space
    Pt->Grad_z(W); W->Set_state('p');// in (Q,F)-space
    W->Trans(W, f_to_p); //go to (Q,P)-space 

    Pt->Trans(Pt, f_to_p); //go to (Q,P)-space 
    dvadd(htot, Pt->base_h, 1, up0,  1, up0, 1); // now up0 has the divergence in (Q,P)-space

    dvadd(htot, ut0, 1, Wf->base_h,  1, Wf->base_h, 1);
    dvadd(htot, vt0, 1, Wf->base_h,  1, Wf->base_h, 1);
    dsmul(htot, -niu, W->base_h,  1, Wf->base_h, 1); // laplacians

    Wt->Trans(Wt, f_to_p); //go to (Q,P)-space 
//
//
//
    dvadd  (htot, ut1, 1, Uf->base_h,  1, ut1, 1);
    dvadd  (htot, vt1, 1, Vf->base_h,  1, vt1, 1);
    dvadd  (htot, wt1, 1, Wf->base_h,  1, wt1, 1);
    
#ifdef NATURAL_CONVECTION
    double     Pr = dparam("PRANDTL");
    double     Ra = dparam("DRAYLEIGH");
  switch (iparam("DIRGRAVITY"))
   {
      case 0:
//         dvsub (htot, ut1, 1, tk0, 1, ut1,1);
         for(i=0; i<htot; ++i) ut1[i] -= Ra*Pr*tk0[i];
      break;
      case 1:
//         dvadd (htot, vt1, 1, tk0, 1, vt1,1);
         for(i=0; i<htot; ++i) vt1[i] -= Ra*Pr*tk0[i];
      break;
      case 2:
//         dvadd (htot, wt1, 1, tk0, 1, wt1,1);
         for(i=0; i<htot; ++i) wt1[i] -= Ra*Pr*tk0[i];
      break;
   }
#endif
    dvmul  (htot, Ut->base_h, 1, ut1, 1, Pt->base_h, 1);
    dvvtvp (htot, Vt->base_h, 1, vt1, 1, Pt->base_h, 1, Pt->base_h,1);
    dvvtvp (htot, Wt->base_h, 1, wt1, 1, Pt->base_h, 1, Pt->base_h,1);// final residual
//    for (int k = 0; k < nz; k++)
//      for(i=0; i<Ut->htot; ++i) {
//         Pt->base_h[k*Ut->htot+i] = (Ut->base_h[k*Ut->htot+i]-up1[i])*ut1[k*Ut->htot+i]
//                                   +(Vt->base_h[k*Ut->htot+i]-vp1[i])*vt1[k*Ut->htot+i]
//                                   +(Wt->base_h[k*Ut->htot+i]-wp1[i])*wt1[k*Ut->htot+i];
//    }
#ifdef SPM
    for (int k = 0; k < nz; k++)
      for(i=0; i<Ut->htot; ++i) 
        Pt->base_h[k*Ut->htot+i] *=(1.0-CONCENTR->base_h[k*Ut->htot+i]); 
#endif

#else  // the following is residual shown by Equ.2 in Guermond et al J Sci Comput 49:35-50 2011

   #ifdef MAP
    fprintf(stderr, "Not work ....afx have convective term ...\n");
    exit(-1);
   #endif
   
//    fprintf(stderr, "Need to debug for 3/2 dealiasing rule for this branch ....\n");
//    exit(-1);

    dvmul  (htot, u0, 1, u0, 1, ut0, 1);
    dvvtvp (htot, v0, 1, v0, 1, ut0, 1, ut0,1);
    dvvtvp (htot, w0, 1, w0, 1, ut0, 1, ut0,1);// u^{n-1}*u^{n-1} in (Q,P)-space
    dsmul  (htot, 0.5, ut0,  1, ut0, 1); 

    dvmul  (htot, u1, 1, u1, 1, ut1, 1);
    dvvtvp (htot, v1, 1, v1, 1, ut1, 1, ut1,1);
    dvvtvp (htot, w1, 1, w1, 1, ut1, 1, ut1,1);// u^{n}*u^{n} in (Q,P)-space
    dsmul  (htot, 0.5, ut1,  1, ut1, 1); 

    dvsub  (htot, ut1, 1, ut0,  1, Ut->base_h, 1);
    dsmul  (htot, dtinv,  Ut->base_h,  1, wt1, 1); //wt1 has the \frac{\partial 0.5(u^2)}{\partial t} 

    dsmul  (htot, 0.5, ut0,  1, Ut->base_h, 1); 
    daxpy  (htot, 0.5, ut1,  1, Ut->base_h, 1); 

    Ut ->Grad_h(Ut->base_h,vt0,wt0,null,trip);
//    Ut ->Grad_d(vt0,wt0,NULL,'a');
    Ut->Trans(Ut, p_to_f); //go to (Q,F)-space 
    Ut->Grad_z(Pt); Pt->Set_state('p'); //in (Q,F)-space
    dcopy(htot, Pt->base_h, 1, vt1, 1); //in (Q,F)-space


    dcopy(htot, vt0, 1, Pt->base_h, 1);
//    Pt->Grad_d(vt0,null,null,'x'); //vt0 has \frac{\partial^2{0.5(u^2)}}{\partial x^2} in (Q,P)-space
    Pt->Grad_h(Pt->base_h,vt0,null,null,trip);

    dcopy(htot, wt0, 1, Pt->base_h, 1);
//    Pt->Grad_d(null,wt0,null,'y');//wt0 has \frac{\partial^2{0.5(u^2)}}{\partial y^2} in (Q,P)-space
    Pt->Grad_h(Pt->base_h,null,wt0,null,trip);

    dcopy(htot, vt1, 1, Pt->base_h, 1);
    Pt->Grad_z(W); W->Set_state('p');// in (Q,F)-space
    W->Trans(W, f_to_p); //go to (Q,P)-space 


    daxpy  (htot, -niu, vt0,  1, wt1, 1); 
    daxpy  (htot, -niu, wt0,  1, wt1, 1); 
    daxpy  (htot, -niu, W->base_h,  1, wt1, 1);  //subtract Laplacian 

    dvadd  (htot, p0, 1, ut0,  1, Pt->base_h, 1);
    dvmul  (htot, u0, 1, Pt->base_h, 1, Ut->base_h, 1); //u*(0.5 u^2+p)
    dvmul  (htot, v0, 1, Pt->base_h, 1, Vt->base_h, 1); //v*(0.5 u^2+p)
    dvmul  (htot, w0, 1, Pt->base_h, 1, Wt->base_h, 1); //w*(0.5 u^2+p)

    dvadd  (htot, p1, 1, ut1,  1, Pt->base_h, 1);
    dvmul  (htot, u1, 1, Pt->base_h, 1, vt0, 1); //u*(0.5 u^2+p)
    dvmul  (htot, v1, 1, Pt->base_h, 1, wt0, 1); //v*(0.5 u^2+p)
    dvmul  (htot, w1, 1, Pt->base_h, 1, vt1, 1); //w*(0.5 u^2+p)

    dsmul  (htot, 0.5, Ut->base_h,  1, Ut->base_h, 1); 
    daxpy  (htot, 0.5, vt0,  1, Ut->base_h, 1); 
//    Ut->Grad_d(vt0,null,null,'x'); //
    Ut->Grad_h(Ut->base_h,vt0,null,null,trip);
    dvadd  (htot, wt1, 1, vt0,  1, wt1, 1);

    dsmul  (htot, 0.5, Vt->base_h,  1, Vt->base_h, 1); 
    daxpy  (htot, 0.5, wt0,  1, Vt->base_h, 1); 
//    Vt->Grad_d(null,wt0,null,'y'); //
    Vt->Grad_h(Vt->base_h,null,wt0,null,trip);
    dvadd  (htot, wt1, 1, wt0,  1, wt1, 1);

    dsmul  (htot, 0.5, Wt->base_h,  1, Wt->base_h, 1); 
    daxpy  (htot, 0.5, vt1,  1, Wt->base_h, 1); 
    Wt->Trans(Wt, p_to_f); //go to (Q,F)-space 
    Wt->Grad_z(W); W->Set_state('p');// in (Q,F)-space
    W->Trans(W, f_to_p); //go to (Q,P)-space 
    dvadd  (htot, wt1, 1, W->base_h,  1, U->base_h, 1);  //now U->base_h has the residual except the  gradient^2

    dsmul (htot, 0.5, u0, 1, Ut->base_h, 1);
    daxpy (htot, 0.5, u1, 1, Ut->base_h, 1);
    dsmul (htot, 0.5, v0, 1, Vt->base_h, 1);
    daxpy (htot, 0.5, v1, 1, Vt->base_h, 1);
    dsmul (htot, 0.5, w0, 1, Wt->base_h, 1);
    daxpy (htot, 0.5, w1, 1, Wt->base_h, 1);

   #ifdef MAP
    dsmul (htot, 0.5, afx0, 1, ut0, 1);
    daxpy (htot, 0.5, afx1, 1, ut0, 1);
    dsmul (htot, 0.5, afy0, 1, vt0, 1);
    daxpy (htot, 0.5, afy1, 1, vt0, 1);
    dsmul (htot, 0.5, afz0, 1, wt0, 1);
    daxpy (htot, 0.5, afz1, 1, wt0, 1);
    
    dvmul  (htot, Ut->base_h, 1, ut0, 1, Pt->base_h, 1); //u*fx
    dvvtvp (htot, Vt->base_h, 1, vt0, 1, Pt->base_h, 1, Pt->base_h,1); //+v*fy
    dvvtvp (htot, Wt->base_h, 1, wt0, 1, Pt->base_h, 1, Pt->base_h,1); //+v*fy

    dvsub  (htot, U->base_h, 1, Pt->base_h,  1, U->base_h, 1);
   #endif


    Ut->Grad_h(Ut->base_h,ut0,vt0,null,trip);
  //  Ut ->Grad_d(ut0,vt0,NULL,'a');
    Ut->Trans(Ut, p_to_f); //go to (Q,F)-space 
    Ut->Grad_z(Pt); Pt->Set_state('p'); //in (Q,F)-space
    Pt->Trans(Pt, f_to_p); //go to (Q,P)-space 

    dcopy(htot, ut0, 1, up0, 1); //for divergence

    dvmul  (htot, Pt->base_h, 1, Pt->base_h, 1, Ut->base_h, 1); //
    dvvtvp (htot, ut0, 1, ut0, 1, Ut->base_h, 1, Ut->base_h,1);
    dvvtvp (htot, vt0, 1, vt0, 1, Ut->base_h, 1, Ut->base_h,1);//


    Vt->Grad_h(Vt->base_h,ut0,vt0,null,trip);
//    Vt ->Grad_d(ut0,vt0,NULL,'a');
    Vt->Trans(Vt, p_to_f); //go to (Q,F)-space 
    Vt->Grad_z(Pt);Pt->Set_state('p'); //in (Q,F)-space
    Pt->Trans(Pt, f_to_p); //go to (Q,P)-space 

    dvadd(htot, vt0, 1, up0,  1, up0, 1);

    dvvtvp (htot, Pt->base_h, 1, Pt->base_h, 1, Ut->base_h, 1, Ut->base_h,1);
    dvvtvp (htot, ut0, 1, ut0, 1, Ut->base_h, 1, Ut->base_h,1);
    dvvtvp (htot, vt0, 1, vt0, 1, Ut->base_h, 1, Ut->base_h,1);//

    Wt->Grad_h(Wt->base_h,ut0,vt0,null,trip);
//    Wt ->Grad_d(ut0,vt0,NULL,'a');
    Wt->Trans(Wt, p_to_f); //go to (Q,F)-space 
    Wt->Grad_z(Pt); Pt->Set_state('p'); //in (Q,F)-space
    Pt->Trans(Pt, f_to_p); //go to (Q,P)-space 
    
    dvadd(htot, Pt->base_h, 1, up0,  1, up0, 1); // now up0 has the divergence in (Q,P)-space

    dvvtvp (htot, Pt->base_h, 1, Pt->base_h, 1, Ut->base_h, 1, Ut->base_h,1);
    dvvtvp (htot, ut0, 1, ut0, 1, Ut->base_h, 1, Ut->base_h,1);
    dvvtvp (htot, vt0, 1, vt0, 1, Ut->base_h, 1, Ut->base_h,1);//

    dsmul  (htot, niu, Ut->base_h,  1, Ut->base_h, 1);
    dvadd  (htot, Ut->base_h, 1, U->base_h,  1, Pt->base_h, 1);   //now Pt->base_h has all the residual

#endif


    dcopy(hjtot, U->base_hj, 1, us1, 1); // backup ?
    dcopy(hjtot, V->base_hj, 1, vs1, 1);
    dcopy(hjtot, W->base_hj, 1, ws1, 1);
  
    dsmul(htot_nz, BBeta[0], omega->u[0],  1, U->base_h, 1); //in (Q,F)-sapce
    daxpy(htot_nz, BBeta[1], omega->u[1],  1, U->base_h, 1); // 
    dsmul(htot_nz, BBeta[0], omega->v[0],  1, V->base_h, 1);
    daxpy(htot_nz, BBeta[1], omega->v[1],  1, V->base_h, 1); // 
    dsmul(htot_nz, BBeta[0], omega->w[0],  1, W->base_h, 1);
    daxpy(htot_nz, BBeta[1], omega->w[1],  1, W->base_h, 1); // 

    U->Grad_z_d(st[2]);
    U->Grad_d(st[0], st[1], 0, 'a'); // NOTE: 'a' instead of 'A'
    // Now (st[0],st[1],st[2])[htot_nz] contains -(d/dx,d/dy,d/dz) of U^{*,n+1} in (Q,F)-space
    V->Grad_z_d(st[5]);
    V->Grad_d(st[3], st[4], 0, 'a'); // NOTE: 'a' instead of 'A'

    W->Grad_z_d(st[8]);
    W->Grad_d(st[6], st[7], 0, 'a'); // NOTE: 'a' instead of 'A'

    U->Trans(U, f_to_p); //go to (Q,P)-space 
    V->Trans(V, f_to_p); //go to (Q,P)-space 
    W->Trans(W, f_to_p); //go to (Q,P)-space 

#ifdef DYNAMIC_ALPHA
//    dvmul  (htot, U->base_h, 1, U->base_h,  1, Ut->base_h, 1);   //
//    dvvtvp (htot, V->base_h, 1, V->base_h, 1,  Ut->base_h, 1, Ut->base_h,1);
//    dvvtvp (htot, W->base_h, 1, W->base_h, 1,  Ut->base_h, 1, Ut->base_h,1);
//    dsmul  (htot, 0.5, Ut->base_h,  1, Ut->base_h, 1);  //Ut->base_h has 0.5*(u*u+v*v+w*w) in (Q,P)-space
    
    //
    SetReducedPhys(U, Ut, Wt); //Ut input Vt output ---base_h has the reduced velocity in (Q,P)-space  //Wt workspace
    dcopy(htot, Ut->base_h, 1, ut0, 1);// ut0 has the reduced velocity u

    SetReducedPhys(V, Ut, Wt); //Ut input Vt output ---base_h has the reduced velocity in (Q,P)-space  //Wt workspace
    dcopy(htot, Ut->base_h, 1, vt0, 1);

    SetReducedPhys(W, Ut, Wt); //Ut input Vt output ---base_h has the reduced velocity in (Q,P)-space  //Wt workspace
    dcopy(htot, Ut->base_h, 1, wt0, 1);
#endif

    if (dparam("DAMPING_LENGTH"))
        damping_length = dparam("DAMPING_LENGTH");
    else
        damping_length = 50.;

    double y_plus,yy;

    L = iparam("MODES");;

    qa = U->flevels[0]->fhead->qa,
    qb = U->flevels[0]->fhead->qb;
  	Coord X;
	  X.x  = dvector(0, qa*qb-1);
	  X.y  = dvector(0, qa*qb-1);
    
    Metric *lambda    = (Metric*)calloc(1, sizeof(Metric));
	  lambda->p = dvector(0, qa*qb-1);

    double vis_plan;
    double vis_mean = 0.;
    double num_elem = 0.;

    double mean_vis = 0.;
    double dist,Rex,Cf,u_star;
    double weight =0., volume = 0.;
    double damping_factor = 1.;
    double dz = dparam("LZ")/option("NZTOT");
    
    double a,b,c,s,AC,AC1,Radius,D;


  if(step == 1)
   {
    escales1 = dvector(0, U->nel-1);
    escales2 = dvector(0, U->nel-1);
    dzero(U->nel, escales1, 1);
    dzero(U->nel, escales2, 1);
  
    for(E=U->fhead;E;E=E->next){

 #if 1
     if(E->identify() == Nek_Tri){
      getzw(E->qa, &za, &wa, 'a');
      getzw(E->qb+1, &zb, &wb, 'b');
      
      dxa = 0.5*(za[1]-za[0]);
      dxb = 0.5*(zb[1]-zb[0]);

      sca = (E->vert[1].x-E->vert[0].x)*(E->vert[1].x-E->vert[0].x)+
	    (E->vert[1].y-E->vert[0].y)*(E->vert[1].y-E->vert[0].y);
      
      scb = (E->vert[2].x-E->vert[1].x)*(E->vert[2].x-E->vert[1].x)+
	    (E->vert[2].y-E->vert[1].y)*(E->vert[2].y-E->vert[1].y);
      
      scc = (E->vert[2].x-E->vert[0].x)*(E->vert[2].x-E->vert[0].x)+
	    (E->vert[2].y-E->vert[0].y)*(E->vert[2].y-E->vert[0].y);
      
      escales1[E->id] = min(sqrt(sca), sqrt(min(scb, scc)));
      escales2[E->id] = min(dxa*sqrt(sca), dxb*sqrt(min(scb, scc)));
    }
    else{
      getzw(E->qa, &za, &wa, 'a');
      getzw(E->qb, &zb, &wb, 'a');
    //  getzw(L+1, &za, &wa, 'a');
    //  getzw(L+1, &zb, &wb, 'a');
      
      dxa = 0.5*(za[1]-za[0]);
      dxb = 0.5*(zb[1]-zb[0]);

/*      
      sca = (E->vert[1].x-E->vert[0].x)*(E->vert[1].x-E->vert[0].x)+
	    (E->vert[1].y-E->vert[0].y)*(E->vert[1].y-E->vert[0].y);
      
      scb = (E->vert[2].x-E->vert[1].x)*(E->vert[2].x-E->vert[1].x)+
	    (E->vert[2].y-E->vert[1].y)*(E->vert[2].y-E->vert[1].y);
      
      scc = (E->vert[3].x-E->vert[2].x)*(E->vert[3].x-E->vert[2].x)+
            (E->vert[3].y-E->vert[2].y)*(E->vert[3].y-E->vert[2].y);
      
      scd = (E->vert[3].x-E->vert[0].x)*(E->vert[3].x-E->vert[0].x)+
	    (E->vert[3].y-E->vert[0].y)*(E->vert[3].y-E->vert[0].y);
      
      escales1[E->id] = min(sqrt(min(sca,scc)), sqrt(min(scb, scd)));
      escales2[E->id] = min(dxa*sqrt(min(sca,scc)), dxb*sqrt(min(scb, scd)));
*/  

      sca = (E->vert[1].x-E->vert[0].x)*(E->vert[1].x-E->vert[0].x)+
	    (E->vert[1].y-E->vert[0].y)*(E->vert[1].y-E->vert[0].y);

      scb = (E->vert[2].x-E->vert[3].x)*(E->vert[2].x-E->vert[3].x)+
	    (E->vert[2].y-E->vert[3].y)*(E->vert[2].y-E->vert[3].y);

      scc = (E->vert[2].x-E->vert[0].x)*(E->vert[2].x-E->vert[0].x)+
	    (E->vert[2].y-E->vert[0].y)*(E->vert[2].y-E->vert[0].y);

      scd = (E->vert[3].x-E->vert[1].x)*(E->vert[3].x-E->vert[1].x)+
	    (E->vert[3].y-E->vert[1].y)*(E->vert[3].y-E->vert[1].y);

//      double sx = sqrt(0.5*(sca+scb));
//      double sy = sqrt(0.5*(scc+scd));

//      escales1[E->id] = min(sx, sy);
//      escales2[E->id] = min(dxa*sx, dxb*sy);

      escales1[E->id] = min(sqrt(min(sca,scb)), sqrt(min(scc, scd)));
      escales2[E->id] = min(dxa*sqrt(min(sca,scb)), dxb*sqrt(min(scc, scd)));
//
//      escales1[E->id]=pow(sqrt(0.5*(sca+scb))*dxa*sqrt(0.5*(scc+scd))*dxb,0.5);
//      escales2[E->id]=pow(sqrt(0.5*(sca+scb))*dxa*sqrt(0.5*(scc+scd))*dxb,0.5);

//      escales1[E->id]=pow(sqrt(0.5*(sca+scb))*sqrt(0.5*(scc+scd))*dz,1.0/3.0)/iparam("MODES");
//      escales2[E->id]=pow(sqrt(0.5*(sca+scb))*sqrt(0.5*(scc+scd))*dz,1.0/3.0)/iparam("MODES");

//      fprintf(stdout,"escales[%d] = %lf dxa = %lf dxb = %lf sca = %lf scb = %lf scc = %lf  scd = %lf \n",E->id,escales[E->id],dxa,dxb, sca,scb,scc,scd);
     }
 #else
    register double lx, ly;
    if(E->identify() == Nek_Tri)
      {
      
	    lx = E->vert[1].x-E->vert[0].x; ly = E->vert[1].y-E->vert[0].y;
	    a = sqrt(lx*lx+ly*ly);
	
	    lx = E->vert[2].x-E->vert[1].x; ly = E->vert[2].y-E->vert[1].y;
	    b = sqrt(lx*lx+ly*ly);
	
	    lx = E->vert[0].x-E->vert[2].x; ly = E->vert[0].y-E->vert[2].y;
	    c = sqrt(lx*lx+ly*ly);
	
	    s = 0.5 * ( a + b + c );
	
	    AC  = asin(sqrt((E->vert[1].y-E->vert[0].y)*
			  (E->vert[1].y-E->vert[0].y))/a);
	    AC1 = asin(sqrt((E->vert[2].y-E->vert[0].y)*
			  (E->vert[2].y-E->vert[0].y))/c);
	
	     D = 0.5 * a * c * sin(AC-AC1) * M_PI * M_PI /(iparam("MODES")*   
						      iparam("MODES"));
      }

    if(E->identify() == Nek_Quad)
     {
	     E->GetZW(&za, &wa, &zb, &wb, &zc, &wc);
       
      double elmt_area = 0.0;
      for (int i=0; i<qa; i++ )
      for (int j=0; j<qb; j++ )
       {    
        
         if( E->curvX )
        weight = wb[j]*wa[i]*E->geom->jac.p[m];
       else
        weight = wb[j]*wa[i]*E->geom->jac.d;
        
        elmt_area += weight;
      }

      D = sqrt(elmt_area)/iparam("MODES");
     }

	    D = pow(sqrt(D*D) * dparam("LZ") /(nztot),0.333);  

      escales1[E->id] = D;
      escales2[E->id] = D;

#endif

    }

   } //step ==1
    step++;
    double max_q2_tmp = -1e9;
    double min_q2_tmp = 1e9;
    for (int k = 0; k < nz; k++)
    {
     int elem_num = 0;

     int plane = myid*nz+k;
     for(int eid =0; eid<U->nel; ++eid)
      {

	     E   = U->flevels[0]->flist[eid];
	     F   = V->flevels[0]->flist[eid];
	     G   = W->flevels[0]->flist[eid];
	     Rt  = Pt->flevels[0]->flist[eid];

	     E->GetZW(&za, &wa, &zb, &wb, &zc, &wc);
       qa = E->qa;   qb = E->qb;
	     E->coord(&X);

        element_diameter = escales1[eid];
//        qll_diameter = element_diameter/iparam("MODES");
        qll_diameter = escales2[eid];

       int offset = U->htot*k+eid*E->qtot; 

        double local_max_velocity = -10e6;
        double local_max_k = -1e7;
#ifdef MAP
        double local_max_vel_map = -1e7;
#endif
        for ( int i=0; i<qa; i++ )
         for ( int j=0; j<qb; j++ )
          {
            m = i*qb + j;

            local_max_velocity =  max( sqrt(U->base_h[offset+m]*U->base_h[offset+m]
                                           +V->base_h[offset+m]*V->base_h[offset+m]
                                           +W->base_h[offset+m]*W->base_h[offset+m]),
                                              local_max_velocity);
#ifdef MAP
            local_max_vel_map =  max( sqrt((U->base_h[offset+m]+map_vx)*(U->base_h[offset+m]+map_vx)
                                           +(V->base_h[offset+m]+map_vy)*(V->base_h[offset+m]+map_vy)
                                           +W->base_h[offset+m]*W->base_h[offset+m]),
                                           local_max_vel_map);
#endif
          }

//        first_order_viscosity =  CB* local_max_velocity* element_diameter;
        first_order_viscosity =  CB* local_max_velocity* qll_diameter;

        if(local_max_k<1e-10) local_max_k = 1e-3*entropy_variation;

        double alpha = CA;

        volume = 0.;
        element_max_viscosity = 0.;

#ifdef DYNAMIC_ALPHA
        double numerator = 0.;
        double denominator = 0.;
        for ( int i=0; i<qa; i++ )
         for ( int j=0; j<qb; j++ )
          {
            m = i*qb + j;

//            numerator += (Ut->base_h[offset+m]-Vt->base_h[offset+m])
//                        *(Ut->base_h[offset+m]-Vt->base_h[offset+m]);
//            denominator += Ut->base_h[offset+m]*Ut->base_h[offset+m];

            numerator += (U->base_h[offset+m]-ut0[offset+m])*(U->base_h[offset+m]-ut0[offset+m])
                        +(V->base_h[offset+m]-vt0[offset+m])*(V->base_h[offset+m]-vt0[offset+m])
                        +(W->base_h[offset+m]-wt0[offset+m])*(W->base_h[offset+m]-wt0[offset+m]);

            denominator += U->base_h[offset+m]*U->base_h[offset+m]
                          +V->base_h[offset+m]*V->base_h[offset+m]
                          +W->base_h[offset+m]*W->base_h[offset+m];

          }

        double elmtSensor = sqrt(numerator/denominator);


        max_sensor = max(max_sensor,elmtSensor);
#endif


//        double element_pres_grad = 0.0;
//        double element_volume = 0.0;
        double elmt_max_x = 0.0;
        double elmt_max_y = 0.0;
        for (int i=0; i<qa; i++ )
        for (int j=0; j<qb; j++ )
         {    
          m = i*qb + j;
          
          entropy_residual = fabs( Pt->base_h[offset+m]);

          alpha_scale = 1.0;

#ifdef DYNAMIC_ALPHA
          alpha *=elmtSensor;
#endif

	        entropy_viscosity = alpha * qll_diameter * qll_diameter *
					                        			    entropy_residual / entropy_variation;
         
    
//    #ifdef OUTPUT_VISCOSITY
//           omega->visc1->base_h[offset+m] = alpha;
//    #endif

          if(element_max_viscosity<entropy_viscosity)
           {
             element_max_viscosity = entropy_viscosity;

             elmt_max_x= X.x[m];
             elmt_max_y= X.y[m];
           }

        }
        double element_entropy_viscosity = min(element_max_viscosity,first_order_viscosity);
          
#ifdef DYNAMIC_ALPHA
    #ifdef OUTPUT_VISCOSITY
        dfill(qa*qb, elmtSensor, omega->visc1->base_h+offset, 1);
    #endif
#endif

        damping_factor  = 1.;
        AA = A;
        
        if(dparam("DWALL_ADAPTING"))
         {

           double dl = damping_length;
           if(dparam("CHANNEL"))
//            y_plus = Re_tau*(1.-fabs(elmt_cy));
            y_plus = Re_tau*(1.-fabs(elmt_max_y));

           if(dparam("PIPE"))
//            y_plus = Re_tau*(0.5-sqrt(elmt_cx*elmt_cx+elmt_cy*elmt_cy));
            y_plus = Re_tau*(0.5-sqrt(elmt_max_x*elmt_max_x+elmt_max_y*elmt_max_y));

           if(dparam("CYLINDER"))
            {
             int spm_offset = plane*9;

                 my_radius = dparam("RADIUS");
        

            double dr = sqrt(pow(elmt_max_x,2.)+pow(elmt_max_y,2));
            dist = dr-my_radius;
             Rex = dist*Re;
           //  Cf = pow(2.*log10(Rex)-0.65,-2.3);
             Cf = 0.058*pow(Re,-0.2);
             u_star = sqrt(0.5*Cf);
             y_plus = Re*u_star*dist;

            }
           
             double A_plus = A;
             double dl_plus = dl;



            if( (y_plus<dl_plus) )
             { 
               damping_factor = pow(1.-exp(-y_plus/A_plus),2.);
             

               element_entropy_viscosity *= damping_factor;
          //    if(k == 0 )
          //    ROOTONLY fprintf(stderr,"A = %lf dl = %lf  factor = %lf  y_plus = %lf cx = %lf  cy = %lf \n",
          //         A_plus,dl_plus,damping_factor,y_plus,elmt_max_x, elmt_max_y);
             }

          } //end of wall adapting

        double scaling_factor = 1.0;
        double dr = sqrt(pow(elmt_max_x,2.)+pow(elmt_max_y,2));
        double elmt_angle = atan2(elmt_max_y,elmt_max_x);
        double dist = dr-my_radius;


        double ex = elmt_max_x/dr; 
        double ey = elmt_max_y/dr;
        if( use_reduction && (dist<cutoff_dist*boundary_layer_thickness))
         {
//#ifdef MAP
//     //     double effect_Re = sqrt((-1.0+map_vx)*(-1.0+map_vx)+map_vy*map_vy)*Re;
//    //      if((effect_Re>=1.0e5)&&(effect_Re<=2.0e5))
//           double effect_Re = local_max_vel_map*Re;
//          if(effect_Re<=3.5e5)
//#endif
//        if( use_reduction && (dist<cutoff_dist*boundary_layer_thickness) && (local_mean_vel<= 0.99))
//        if( use_reduction && (dist<boundary_layer_thickness) && ((map_ax*ex+map_ay*ey)>0.1) )
         {
//           scaling_factor = 1.0+tanh(-Am/cutoff_accel);
           scaling_factor = pow(1.-exp(-dist/boundary_layer_thickness/cutoff_dist),2.);
           element_entropy_viscosity *= scaling_factor;
#if 0
         double WALE_viscosity = -10e6; 
        for ( int i=0; i<qa; i++ )
         for ( int j=0; j<qb; j++ )
          {
            m = i*qb + j;
           
            double gij_gij_11 = st[0][offset+m]*st[0][offset+m]+st[1][offset+m]*st[3][offset+m]+st[2][offset+m]*st[6][offset+m];
            double gij_gij_12 = st[0][offset+m]*st[1][offset+m]+st[1][offset+m]*st[4][offset+m]+st[2][offset+m]*st[7][offset+m];
            double gij_gij_13 = st[0][offset+m]*st[2][offset+m]+st[1][offset+m]*st[5][offset+m]+st[2][offset+m]*st[8][offset+m];

            double gij_gij_21 = st[3][offset+m]*st[0][offset+m]+st[4][offset+m]*st[3][offset+m]+st[5][offset+m]*st[6][offset+m];
            double gij_gij_22 = st[3][offset+m]*st[1][offset+m]+st[4][offset+m]*st[4][offset+m]+st[5][offset+m]*st[7][offset+m];
            double gij_gij_23 = st[3][offset+m]*st[2][offset+m]+st[4][offset+m]*st[5][offset+m]+st[5][offset+m]*st[8][offset+m];

            double gij_gij_31 = st[6][offset+m]*st[0][offset+m]+st[7][offset+m]*st[3][offset+m]+st[8][offset+m]*st[6][offset+m];
            double gij_gij_32 = st[6][offset+m]*st[1][offset+m]+st[7][offset+m]*st[4][offset+m]+st[8][offset+m]*st[7][offset+m];
            double gij_gij_33 = st[6][offset+m]*st[2][offset+m]+st[7][offset+m]*st[5][offset+m]+st[8][offset+m]*st[8][offset+m];
                             
            double gji_gji_11 = st[0][offset+m]*st[0][offset+m]+st[3][offset+m]*st[1][offset+m]+st[6][offset+m]*st[2][offset+m];
            double gji_gji_12 = st[0][offset+m]*st[3][offset+m]+st[3][offset+m]*st[4][offset+m]+st[6][offset+m]*st[5][offset+m];
            double gji_gji_13 = st[0][offset+m]*st[6][offset+m]+st[3][offset+m]*st[7][offset+m]+st[6][offset+m]*st[8][offset+m];

            double gji_gji_21 = st[1][offset+m]*st[0][offset+m]+st[4][offset+m]*st[1][offset+m]+st[7][offset+m]*st[2][offset+m];
            double gji_gji_22 = st[1][offset+m]*st[3][offset+m]+st[4][offset+m]*st[4][offset+m]+st[7][offset+m]*st[5][offset+m];
            double gji_gji_23 = st[1][offset+m]*st[6][offset+m]+st[4][offset+m]*st[7][offset+m]+st[7][offset+m]*st[8][offset+m];

            double gji_gji_31 = st[2][offset+m]*st[0][offset+m]+st[5][offset+m]*st[1][offset+m]+st[8][offset+m]*st[2][offset+m];
            double gji_gji_32 = st[2][offset+m]*st[3][offset+m]+st[5][offset+m]*st[4][offset+m]+st[8][offset+m]*st[5][offset+m];
            double gji_gji_33 = st[2][offset+m]*st[6][offset+m]+st[5][offset+m]*st[7][offset+m]+st[8][offset+m]*st[8][offset+m];


            double gkk_gkk_11 = (st[0][offset+m]+st[4][offset+m]+st[8][offset+m])*(st[0][offset+m]+st[4][offset+m]+st[8][offset+m]);
            double gkk_gkk_22 = gkk_gkk_11;
            double gkk_gkk_33 = gkk_gkk_11;
           
            double vij_11 = 0.5*(gij_gij_11+gji_gji_11)-1./3.*gkk_gkk_11;
            double vij_12 = 0.5*(gij_gij_12+gji_gji_12);
            double vij_13 = 0.5*(gij_gij_13+gji_gji_13);

            double vij_21 = 0.5*(gij_gij_21+gji_gji_21);
            double vij_22 = 0.5*(gij_gij_22+gji_gji_22)-1./3.*gkk_gkk_22;
            double vij_23 = 0.5*(gij_gij_23+gji_gji_23);

            double vij_31 = 0.5*(gij_gij_31+gji_gji_31);
            double vij_32 = 0.5*(gij_gij_32+gji_gji_32);
            double vij_33 = 0.5*(gij_gij_33+gji_gji_33)-1./3.*gkk_gkk_33;

            double vkk = vij_11*vij_11+vij_12*vij_12+vij_13*vij_13
                        +vij_21*vij_21+vij_22*vij_22+vij_23*vij_23
                        +vij_31*vij_31+vij_32*vij_32+vij_33*vij_33;

            double sij_11 = 2.0*st[0][offset+m];
            double sij_12 = st[1][offset+m]+st[3][offset+m];
            double sij_13 = st[2][offset+m]+st[6][offset+m];

            double sij_21 = st[3][offset+m]+st[1][offset+m];
            double sij_22 = 2.0*st[4][offset+m];
            double sij_23 = st[5][offset+m]+st[7][offset+m];

            double sij_31 = st[6][offset+m]+st[2][offset+m];
            double sij_32 = st[7][offset+m]+st[5][offset+m];
            double sij_33 = 2.0*st[8][offset+m];
            
            double skk = sij_11*sij_11+sij_12*sij_12+sij_13*sij_13
                        +sij_21*sij_21+sij_22*sij_22+sij_23*sij_23
                        +sij_31*sij_31+sij_32*sij_32+sij_33*sij_33;

            double  denominator = pow(skk,2.5)+pow(vkk,1.25);
            double local_WALE_viscosity = 0.0;
            if(fabs(denominator)>1e-6) 
            local_WALE_viscosity = 
            WALE_A*WALE_A*qll_diameter*qll_diameter
                                          *pow(vkk,1.5)/denominator;

//            if(local_WALE_viscosity<0.0) local_WALE_viscosity = 0.0;

            WALE_viscosity = max(WALE_viscosity,local_WALE_viscosity);

//        if(denominator<1e-6)
//        fprintf(stderr, "WALE = %g  denominator = %g \n",local_WALE_viscosity,denominator);
          }

//        fprintf(stderr, "WALE = %g entropy = %g \n",WALE_viscosity, element_entropy_viscosity);
           scaling_factor = WALE_viscosity/element_entropy_viscosity;
           element_entropy_viscosity = min(WALE_viscosity,element_entropy_viscosity);
#endif
         }
        }

#ifdef OUTPUT_VISCOSITY
  #ifdef DYNAMIC_ALPHA
         dsmul(qa*qb, scaling_factor, omega->visc1->base_h+offset, 1, omega->visc1->base_h+offset, 1);
  #else
         dfill(qa*qb, scaling_factor, omega->visc1->base_h+offset, 1);
  #endif
#endif

           for (int i=0; i<qa; i++ )
           for (int j=0; j<qb; j++ )
           {    
             m = i*qb + j;
             visc_tmp[k*U->htot+eid*E->qtot+m] = element_entropy_viscosity;
           }
          
         }

       }

#ifndef EVM_NOT_APPLIED
      dcopy(U->htot,visc_ave,1,visc_bak,1);
      dzero(U->htot, visc_ave, 1);

     for(int j=0; j<Ut->htot; j++)
      {
        visc_ave[j] = visc_tmp[j];
        for(int i=1; i<nz; ++i)
           visc_ave[j] = max(visc_ave[j],visc_tmp[i*Ut->htot+j]);
      }

//       MPI_Allreduce (visc_ave, visc_tmp, Ut->htot, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
       MPI_Allreduce (visc_ave, &Pt->base_h[0], Ut->htot, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//       dcopy(U->htot,visc_tmp,1,visc_ave,1);
       dcopy(U->htot,Pt->base_h,1,visc_ave,1);


       if(step > 2)
       for(int i=0; i<U->htot; ++i)
         visc_ave[i] =2./3.*visc_ave[i]+1./3.*visc_bak[i];

    #ifdef OUTPUT_VISCOSITY
     for (int k = 0; k < U->nz; k++)
         dcopy(U->htot, visc_ave, 1, omega->visc2->base_h+k*U->htot, 1);
    #endif

       for(int i=0; i<U->htot; ++i)
         max_viscosity = max(max_viscosity,visc_ave[i]);

       

      Ut->Set_state('p');
      Vt->Set_state('p');
      Wt->Set_state('p');

//In (Q,F)-space
     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
        ut0[k*U->htot+j] = visc_ave[j]*st[0][k*U->htot+j];
//        ut0[k*U->htot+j] = visc_ave[j]*(st[0][k*U->htot+j]+st[0][k*U->htot+j]);
//        ut0[k*U->htot+j] = visc_tmp[k*U->htot+j]*(st[0][k*U->htot+j]+st[0][k*U->htot+j]);
//     

     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
        vt0[k*U->htot+j] = visc_ave[j]*st[1][k*U->htot+j];
//        vt0[k*U->htot+j] = visc_ave[j]*(st[1][k*U->htot+j]+st[3][k*U->htot+j]);
//        vt0[k*U->htot+j] = visc_tmp[k*U->htot+j]*(st[1][k*U->htot+j]+st[3][k*U->htot+j]);
//
     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
        wt0[k*U->htot+j] = visc_ave[j]*st[2][k*U->htot+j];
//        wt0[k*U->htot+j] = visc_ave[j]*(st[2][k*U->htot+j]+st[6][k*U->htot+j]);
//        wt0[k*U->htot+j] = visc_tmp[k*U->htot+j]*(st[2][k*U->htot+j]+st[6][k*U->htot+j]);

//      Ut->Trans(U->htot, ut0, ut0, p_to_f); // in(Q,F)-space
//      Vt->Trans(U->htot, vt0, vt0, p_to_f); // in(Q,F)-space
//      Wt->Trans(U->htot, wt0, wt0, p_to_f); // in(Q,F)-space
      FourierList_Iprod_div(Ut,	ut0, vt0, wt0, Pt->base_hj);
      daxpy(hjtot, Re, Ut->base_hj, 1, Uf->base_hj, 1);

      Ut->Set_state('p');
      Pt->Set_state('p');
//
//
//     

     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
        ut0[k*U->htot+j] = visc_ave[j]*st[3][k*U->htot+j];
//        ut0[k*U->htot+j] = visc_ave[j]*(st[3][k*U->htot+j]+st[1][k*U->htot+j]);
//        ut0[k*U->htot+j] = visc_tmp[k*U->htot+j]*(st[3][k*U->htot+j]+st[1][k*U->htot+j]);

     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
        vt0[k*U->htot+j] = visc_ave[j]*st[4][k*U->htot+j];
//        vt0[k*U->htot+j] = visc_ave[j]*(st[4][k*U->htot+j]+st[4][k*U->htot+j]);
//        vt0[k*U->htot+j] = visc_tmp[k*V->htot+j]*(st[4][k*U->htot+j]+st[4][k*U->htot+j]);
//

     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
        wt0[k*U->htot+j] = visc_ave[j]*st[5][k*U->htot+j];
//        wt0[k*U->htot+j] = visc_ave[j]*(st[5][k*U->htot+j]+st[7][k*U->htot+j]);
//        wt0[k*U->htot+j] = visc_tmp[k*V->htot+j]*(st[5][k*U->htot+j]+st[7][k*U->htot+j]);


//      Ut->Trans(U->htot, ut0, ut0, p_to_f); // in(Q,F)-space
//      Vt->Trans(U->htot, vt0, vt0, p_to_f); // in(Q,F)-space
//      Wt->Trans(U->htot, wt0, wt0, p_to_f); // in(Q,F)-space
      FourierList_Iprod_div(Vt,	ut0, vt0, wt0, Pt->base_hj);
      daxpy(hjtot, Re, Vt->base_hj, 1, Vf->base_hj, 1);

      Vt->Set_state('p');
      Pt->Set_state('p');
//
//

     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
       ut0[k*U->htot+j] = visc_ave[j]*st[6][k*U->htot+j];
//       ut0[k*U->htot+j] = visc_ave[j]*(st[6][k*U->htot+j]+st[2][k*U->htot+j]);
//       ut0[k*U->htot+j] = visc_tmp[k*U->htot+j]*(st[6][k*U->htot+j]+st[2][k*U->htot+j]);


     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
       vt0[k*U->htot+j] = visc_ave[j]*st[7][k*U->htot+j];
//       vt0[k*U->htot+j] = visc_ave[j]*(st[7][k*U->htot+j]+st[5][k*U->htot+j]);
//       vt0[k*U->htot+j] = visc_tmp[k*U->htot+j]*(st[7][k*U->htot+j]+st[5][k*U->htot+j]);

     for (int k = 0; k < U->nz; k++)
      for(int j=0; j<U->htot; ++j)
       wt0[k*U->htot+j] = visc_ave[j]*st[8][k*U->htot+j];
//       wt0[k*U->htot+j] = visc_ave[j]*(st[8][k*U->htot+j]+st[8][k*U->htot+j]);
//       wt0[k*U->htot+j] = visc_tmp[k*U->htot+j]*(st[8][k*U->htot+j]+st[8][k*U->htot+j]);

//      Ut->Trans(U->htot, ut0, ut0, p_to_f); // in(Q,F)-space
//      Vt->Trans(U->htot, vt0, vt0, p_to_f); // in(Q,F)-space
//      Wt->Trans(U->htot, wt0, wt0, p_to_f); // in(Q,F)-space
      FourierList_Iprod_div(Wt,	ut0, vt0, wt0, Pt->base_hj);
      daxpy(hjtot, Re, Wt->base_hj, 1, Wf->base_hj, 1);
 
      Wt->Set_state('p');
      Pt->Set_state('p');
#endif
//

    
//    dcopy(htot, ut1, 1, U->base_h, 1); // copy back
//    dcopy(htot, vt1, 1, V->base_h, 1);
//    dcopy(htot, wt1, 1, W->base_h, 1);

    dcopy(hjtot, us1, 1, U->base_hj, 1); // copy back
    dcopy(hjtot, vs1, 1, V->base_hj, 1);
    dcopy(hjtot, ws1, 1, W->base_hj, 1);


    MPI_Allreduce (&max_viscosity, &entropy_viscosity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&entropy_viscosity,1,&max_viscosity,1);

#ifdef DYNAMIC_ALPHA
    MPI_Allreduce (&max_sensor, &entropy_viscosity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dcopy(1,&entropy_viscosity,1,&max_sensor,1);
#endif

    double ht= (double) htot;
    double total_div = 0.;
    for(int i=0; i<htot; ++i)
        total_div += up0[i]*up0[i] ;

    MPI_Allreduce (&total_div, &entropy_viscosity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(1,&entropy_viscosity,1,&total_div,1);
    MPI_Allreduce (&ht, &entropy_viscosity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(1,&entropy_viscosity,1,&ht,1);
    total_div /= ht;
    total_div = sqrt(total_div);

    #ifdef OUTPUT_VISCOSITY
      omega->visc1->Set_state('p');
      omega->visc2->Set_state('p');
      omega->visc1->Trans(omega->visc1, p_to_f); // Get us back to Fourier Space
      omega->visc1->Trans(omega->visc1, Q_to_J); // Get us back to modal space
      omega->visc2->Trans(omega->visc2, p_to_f); // Get us back to Fourier Space
      omega->visc2->Trans(omega->visc2, Q_to_J); // Get us back to modal space
      omega->visc1->Set_state('t');
      omega->visc2->Set_state('t');
    #endif 

//    min_alpha *= -1.;
//    gdmax(&total_div, 1, &entropy_viscosity);
//    gdmax(&min_alpha, 1, &entropy_viscosity);
//    min_alpha *= -1.;

//    ROOTONLY fprintf(stdout,"Minimum alpha ....................... %lf \n", min_alpha);
    ROOTONLY fprintf(stdout,"Maximum entropy viscosity ........... %g \n", max_viscosity);
#ifdef DYNAMIC_ALPHA
    ROOTONLY fprintf(stdout,"Maximum  Sensor ......... ........... %g \n", max_sensor);
#endif
    ROOTONLY fprintf(stdout,"Maximum velocity ...... ............. %g \n", max_velocity);
    ROOTONLY fprintf(stdout,"RMS   divergence ...... ............. %g \n", total_div);
//#ifdef DYNAMIC_ALPHA
//#endif

//	  free_dvector(lambda->p,qa*qb-1);
    free(lambda);

  	free(X.x);
	  free(X.y);
}

   
static void OrthoInnerProduct(Element_List *EL, Element_List *ELf);
static void OrthoJTransBwd(Element_List *EL, Element_List *ELf);

void SetReducedPhys(Element_List *In, Element_List *Out, Element_List *Work) {
     // In->base_h velocity in (Q,P)-space
     // Out->base_h reduced-order velocity in (Q,P)-space

    Element *E, *F;
    static int INIT = 1;

    int     nz = option("NZ"), nztot = option("NZTOT");
    int     nprocs  = option("NPROCS");
    int     myid    = option("PROCID");
    int     hjtot   = In->nz*In->hjtot;
    int     htot    = In->nz*In->htot;
    int     htot_nz = In->nz*In->htot;

    static  Nek_Trans_Type f_to_p = F_to_P,
                           p_to_f = P_to_F;
    char    trip = 'a';
    int     actived_mode = iparam("NACTIVE_MODE");

    if (option("dealias")) {
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
      trip = 'A';
      nz   = 3*nz/2;
      htot = 3*htot/2;
    }

    int Nm=In->flevels[0]->flist[0]->Nmodes;
    int cnt;
    double *filter = dvector(0,Nm-1);
    double *coeff_tmp = dvector(0,Nm-1);
    double evm_theta = (dparam("EVM_THETA"))? dparam("EVM_THETA"): 1.0;
    double evm_eps = 1.0/pow(Nm,evm_theta)/log(Nm);
    int NCutM = 1;
    double Qa,Qb;
    if (iparam("N_CUT_EVM"))
    NCutM = iparam("N_CUT_EVM");
      else
    NCutM = 1; //

    if(NCutM<1) NCutM=1;

    In->Set_state('p');
    In->Trans(In->htot, In->base_h, Work->base_h, p_to_f); // to (Q,F)-space

   if(INIT)
    {
      reset_bases();
      init_ortho_basis();
      INIT = 0;
    }
    
 //   Work->Set_state('p');
 //  Work->Trans(Work, Q_to_J); // Get us back to modal space
 //   Work->Set_state('t');

    OrthoInnerProduct(Work, Out);

    for (int k = 0; k < Out->nz; k++)
     {

     for(int eid =0; eid<Out->nel; ++eid)
      {
	    
       E   = Out->flevels[k]->flist[eid];
       F   = In->flevels[k]->flist[eid];
       int lmax = E->lmax;
    
//       int offset_j = k*Out->hjtot+eid*E->Nmodes;
  /* 
       dfill(Nm,0.0, filter,1);
      
       cnt = 0;
       for (int i = 0; i < NCutM; i++){
        for (int j = NCutM; j < lmax; j++){
         Qb = setQfactor (j, lmax, NCutM, 0);
         filter[cnt+j] = exp(-Qb); 
         }
          cnt += lmax;
        }

        for (int i = NCutM; i < lmax; i++){
          Qa = setQfactor (i, lmax, NCutM, 0);

        for (int j = 0; j < NCutM; j++)
           filter[cnt+j] = exp(-Qa); 

        for (int j = NCutM; j < lmax; j++){
          Qb = setQfactor (j, lmax, NCutM, 0);
          filter[cnt+j] = exp(-(Qa+Qb)); 
          }
          cnt += lmax;
         }
          dvmul(Nm, filter,1, E->vert->hj, 1, E->vert->hj,1);
 */   
//          int offset_j = k*Out->hjtot+eid*E->Nmodes;

       #if 0
          dcopy(Nm, E->vert->hj, 1, coeff_tmp, 1);
          dzero(Nm, E->vert->hj, 1);
          cnt = 0;
          for (int i = 0; i < NCutM+1; ++i)
            {
               dcopy(NCutM, coeff_tmp+cnt,1, E->vert->hj+cnt,1);
               cnt = i*lmax;
            }
       #endif
          for(int i=0; i<lmax; ++i)
            for(int j=0; j<lmax; ++j)
              if( (i+j) >= (NCutM+1))
                E->vert->hj[i*lmax+j] = 0.0;
     
 //              dcopy(Nm, coeff_tmp,1, F->vert->hj,1);

      }
   }


//    Out->Set_state('t');
//    Out->Trans(Out, J_to_Q); // Get us back to physical space
//    Out->Set_state('P');
//    Out->Trans(Out, f_to_p); // to (Q,P)-space

    OrthoJTransBwd(Out, Work);
    Work->Set_state('p');
    Work->Trans(Work->htot, Work->base_h, Out->base_h, f_to_p); // to (Q,P)-space

 //   OrthoJTransBwd(In, Work);
 //   Work->Set_state('p');
 //   Work->Trans(Work->htot, Work->base_h, In->base_h, f_to_p); // to (Q,P)-space

    free(filter);
    free(coeff_tmp);

}

static double setQfactor(int ind, int limit, int low, int high){ //SVV
 double num   = ind - (limit - high),
        denom = ind - low + 1 ,
        Q;

 if(ind>limit)
   num = 0.0;

 Q = num*num/(denom*denom);
 return Q;
}

/* Inner project w.r.t. orthogonal modes */

static void OrthoInnerProduct(Element_List *EL, Element_List *ELf){
  Element *U, *Uf;

  int nz = option("NZ");
//  double   **ba,**bb,*wa,*wb;
//  int qa = EL->flevels[0]->flist[0]->qa;
//  int qb = EL->flevels[0]->flist[0]->qb;
//
//  get_moda_GL (qa, &ba);
//  get_moda_GL (qb, &bb); 
//  if (option("dealias"))  nz   = 3*nz/2;
//  int L = EL->flevels[0]->flist[0]->lmax;
//  double *work = dvector(0,L*L-1);

  for (int k = 0; k < nz; k++)
   {
     for(int eid =0; eid<EL->nel; ++eid)
      {
        U  = EL->flevels[k]->flist[eid];
        Uf = ELf->flevels[k]->flist[eid];

        Uf->Ofwd(U->h[0], Uf->vert->hj, Uf->lmax);
//       dgemm('N','T', L, L, L, 1.0,   bb[0], L,    U->vert->hj, L, 0.0, work, L);
//       dgemm('N','N', L, L, L, 1.0,  work, L, ba[0], L, 0.0,     Uf->vert->hj, L);
      }
   }

//  free(work);

  return;
}

static void OrthoJTransBwd(Element_List *EL, Element_List *ELf){
  Element *U, *Uf;

  int nz = option("NZ");
 // if (option("dealias"))  nz   = 3*nz/2;

  for (int k = 0; k < nz; k++)
   {
     for(int eid =0; eid<EL->nel; ++eid)
      {
        U  = EL->flevels[k]->flist[eid];
        Uf = ELf->flevels[k]->flist[eid];

        Uf->Obwd(U->vert->hj, Uf->h[0], Uf->lmax);
        
      }
   }

  return;
}


#endif
