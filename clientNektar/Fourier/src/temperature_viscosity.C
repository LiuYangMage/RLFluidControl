#include <time.h>
#include "nektarF.h"
#include <stdio.h>

#include "nse_util.h"

#define max(a,b) ( (b) < (a) ? (a) : (b) )
#define min(a,b) ( (b) > (a) ? (a) : (b) )

static double *escales1 = (double*)0;
static double *escales2 = (double*)0;
static double *temp_visc_tmp = (double*)0;
static double *temp_visc_ave = (double*)0;
static double *temp_visc_bak = (double*)0;

void NSE_Iprod_div(Element *E,double *u, double *v);
void NSE_Iprod_div_3d(Element_List *UL, double *u, double *v,double *w);
void NSE_Iprod_div(Element_List *UL, double *u, double *v);
void NSE_Iprod_div_3dd(Element *U, double *u, double *v,double *w);

static int evm_step = 1;
static double temperature_range[2];
static double temperature_entropy_variation;
static double *temperature_entropy_viscosities; 
static double *ts;//for backup

void compute_temperature_range(Domain *omega){
    register int i,j,k;
    int          qa,qb,qc,m;
    double       *Tk0 = omega->tk[1];
    int          nz     = omega->T->nz;
    int          htot   = omega->T->htot;
  	
    double       tmp = 0.;
    double       T_max = -1e8;
    double       T_min = 1e8;

    int NZ       = omega->T->nz;
    int htot_nz  = omega->T->nz*omega->T->htot;
    if (option("dealias")) {
      NZ   = 3*NZ/2;
      htot_nz = 3*htot_nz/2;
    }

     for(i=0; i<htot_nz; ++i)
      {
       T_max = max(Tk0[i],T_max);
       T_min = min(Tk0[i],T_min);
      }

    gdmax(&T_max,  1,  &tmp);

    T_min = -T_min;
    gdmax(&T_min,  1,  &tmp);
    T_min = -T_min;

    temperature_range[0] = T_max;
    temperature_range[1] = T_min;

#ifndef ENTROPYVISCOSITY
    ROOTONLY fprintf(stdout,"temperature_range = %g  %g \n", temperature_range[1],temperature_range[0]);
#endif
}


#ifdef ENTROPYVISCOSITY
void compute_temperature_variation (Domain *omega){
    register int i,j,k,l;
    int          qa,qb,qc,m;
    Element_List *T = omega->T,
                 *Tf = omega->Tf;
    double       **Tk = omega->tk;
    Element      *E,*F,*G;
	  double       *za,*wa,*zb,*wb,*zc,*wc;
    double       weight;
    int          nz     = T->nz;
    int          htot   = T->htot;
    
    int NZ      = nz;
    int htot_nz = htot*nz;
    if (option("dealias")) {
      NZ   = 3*NZ/2;
      htot_nz = 3*htot_nz/2;
    }

    double       temperature_mean  = 0.;
    double       entropy;
    
    temperature_mean = 0.5*(temperature_range[0]+temperature_range[1]);

//    temperature_entropy_variation = 0.;
//
    dcopy(htot_nz, Tk[1], 1, Tf->base_h, 1);

    double tt = 0.;
    double volume = 0.,integ_entropy = 0.;
    double max_entropy = -1e6, min_entropy = 1e6;

//    Coord X;
//	  X.x  = dvector(0, QGmax*QGmax-1);
//	  X.y  = dvector(0, QGmax*QGmax-1);

    for (k = 0; k < NZ; k++)
  //   for(E=Tf->flevels[0]->fhead;E;E=E->next)
     for(int eid =0; eid<Tf->nel; ++eid)
      {
	      E   = Tf->flevels[0]->flist[eid];
        qa = E->qa;   qb = E->qb;
        E->GetZW(&za, &wa, &zb, &wb, &zc, &wc);
        
//	      E->coord(&X);

        int offset = k*Tf->htot+eid*E->qtot; 
        for ( i=0; i<qa; i++ )
        for ( j=0; j<qb; j++ )
        {    
          m = i*qb + j;
          
          if( E->curvX )
            weight = wb[j]*wa[i]*E->geom->jac.p[m];
          else
            weight = wb[j]*wa[i]*E->geom->jac.d;

//          tt = E->h[0][m];
          tt = Tf->base_h[offset+m];

//          if((tt<-0.25) || (tt>1.25))
//           {
//            ROOT fprintf(stderr, " x = %g  y= %g r = %g  k = %d  t = %g \n",X.x[m],X.y[m],sqrt(X.x[m]*X.x[m]+X.y[m]*X.y[m]),k,tt);
//           }
          
          entropy = (tt - temperature_mean)
                   *(tt - temperature_mean);

          volume += weight;
          integ_entropy += entropy*weight;

          max_entropy = max(max_entropy,entropy);
          min_entropy = min(min_entropy,entropy);
        }
    }

    double tmp = 0.;
    
    gdsum(&volume,1,&tmp);
    gdsum(&integ_entropy, 1, &tmp);

    double mean_entropy = integ_entropy/volume;

     gdmax(&max_entropy,1,&tmp);

    min_entropy = -min_entropy;
    gdmax(&min_entropy,1,&tmp);
    
     temperature_entropy_variation = max(max_entropy-mean_entropy,
                                               mean_entropy-(-min_entropy)); 
    evm_step ++;

//    ROOTONLY fprintf(stdout,"temperature_variation[%d] = %lf  %lf  %lf %lf \n", l,temperature_entropy_variation[l],mean_entropy,min_entropy,max_entropy); 
}


void SetReducedPhys(Element_List *In, Element_List *Out, Element_List *Work);

void compute_temperature_viscosity(Domain *omega){
    register int i,j,k,l,q;

    int          qa,qb,qc,m,n;
    double       ns_residual[3];
  	double       *null= (double*)0;
	  double       *za,*wa,*zb,*wb,*zc,*wc;
	  double       dt    = dparam("DELT");
	  double       dtinv = 1./dt;
    int          id0,id1;
    double        medge,elen;
    Element      *E,*F,*G,*H,*R;

	  double       Re     = 1./dparam("KINVIS");
	  double       Pr     = dparam("PRANDTL");
 #ifdef SPM       //diffusion of the solid over fluid
    double       Kr     = dparam("PRANDTL_RATIO")? dparam("PRANDTL_RATIO"):1.0;
 #endif

    double       alpha  = dparam("TEM_STABILIZATION_ALPHA");
    double       beta   = dparam("TEM_STABILIZATION_BETA");
    int          L      = iparam("MODES");;

    double       max_viscosities     = 0.;
    double       max_sensor          = 0.0;


    static int   step = 1;

    double       dxa, dxb, sca, scb, scc, scd;

    double damping_length;
    if (dparam("DAMPING_LENGTH"))
        damping_length = dparam("DAMPING_LENGTH");
    else
        damping_length = 50.;

    double y_plus,yy;
    double Re_tau = dparam("RE_TAU");
    double dist,Rex,Cf,u_star;
    double A = 26.;
   	
    Coord X;
	  X.x  = dvector(0, QGmax*QGmax-1);
	  X.y  = dvector(0, QGmax*QGmax-1);

    if((alpha==0)||(beta==0.))
     {
       ROOTONLY fprintf(stderr,"Forget to set TEM_STABILIZATION_ALPHA or TEM_STABILIZATION_BETA parameter ! \n");
      exit(1);
     }

    Element_List *U    =  omega->U,   *V    =  omega->V,   *W  = omega->W,
                 *Uf   =  omega->Uf,  *Vf   =  omega->Vf,  *Wf = omega->Wf,
                 *Ut   =  omega->Ut,  *Vt   =  omega->Vt,  *Wt = omega->Wt,
                 *P    =  omega->P,   *Pf   =  omega->Pf,  *Pt = omega->Pt;
   #ifdef SPM 
    Element_List  *C    =  omega->CONCENTR;
   #endif

    double       *u0   =  omega->uk[0],    *v0    =  omega->vk[0],    *w0  = omega->wk[0], // velocity at time step n
                 *u1   =  omega->uk[1],    *v1    =  omega->vk[1],    *w1  = omega->wk[1]; // velocity at time step n-1
     double      *Tk1   =  omega->tk[1],     *Tk0   =  omega->tk[0]; // temperature value at time step n, n-1

    double       *ut0  =  omega->ut[0],    *vt0   =  omega->vt[0],    *wt0 = omega->wt[0], // used as temps arrary for gradients
                 *ut1  =  omega->ut[1],    *vt1   =  omega->vt[1],    *wt1 = omega->wt[1]; // used for calculating laplacians

    double       *us0   =  omega->uss[0],   *vs0   =  omega->vss[0],   *ws0  = omega->wss[0], // backup arrary
                 *us1   =  omega->uss[1],   *vs1   =  omega->vss[1],   *ws1  = omega->wss[1];


     Element_List  *Tf    =  omega->Tf,        *T     =  omega->T;

     Metric *lambda    = (Metric*)calloc(1, sizeof(Metric));
	   lambda->p = dvector(0, QGmax*QGmax-1);

		 double RePr = Re * Pr;
     double kappa = 1.0/RePr;
      
     double temperature_mean = 0.5*(temperature_range[0]+temperature_range[1]);
     
     int          nz     = T->nz;
     int          htot   = T->htot;
     int          hjtot  = T->hjtot;
    
     static       Nek_Trans_Type f_to_p = F_to_P,
                                p_to_f = P_to_F;
     char         trip = 'a';

     int NZ       = nz;
     int htot_nz  = htot*nz;
     int htot_nz_A  = htot*nz;
     int hjtot_nz = hjtot*nz;
     if (option("dealias")) {
      NZ   = 3*NZ/2;
      htot_nz_A = 3*htot_nz_A/2;
      trip = 'A';
     }
     double BBeta[3];
     getbeta(BBeta);

     if(!temp_visc_tmp) 
     temp_visc_tmp = dvector(0, htot_nz_A-1);

     if(!temp_visc_ave)
     {
      temp_visc_ave = dvector(0, Ut->htot-1);
      temp_visc_bak = dvector(0, Ut->htot-1);
     }

     Ut->Set_state('p');
     Vt->Set_state('p');
     Wt->Set_state('p');
     Pt->Set_state('p');

     dvsub(htot_nz_A, Tk1, 1, Tk0, 1, ut1, 1);
     dsmul(htot_nz_A, dtinv,  ut1, 1, ut1, 1); 

     dvadd(htot_nz_A, Tk0, 1, Tk1, 1, Pt->base_h, 1); 
     dsmul(htot_nz_A, 0.5, Pt->base_h, 1, Pt->base_h, 1);

//     Pt->Grad_d(Ut->base_h,Vt->base_h,Wt->base_h,'a');
 //    Pt->Grad_z(Wt);
 //    Pt->Grad_d(Ut->base_h,Vt->base_h,NULL,'a');
     Pt->Grad_h(Pt->base_h,Ut->base_h,Vt->base_h,null,trip);

     Pt->Trans(Wt, p_to_f); //go to (Q,F)-space 
     Wt->Grad_z(Pt); Pt->Set_state('p'); //in (Q,F)-space
     Pt->Trans(Wt, f_to_p); //go to (Q,P)-space 

     dsmul (htot_nz_A, 0.5, u0, 1, ut0, 1);
     daxpy (htot_nz_A, 0.5, u1, 1, ut0, 1);
     dsmul (htot_nz_A, 0.5, v0, 1, vt0, 1);
     daxpy (htot_nz_A, 0.5, v1, 1, vt0, 1);
     dsmul (htot_nz_A, 0.5, w0, 1, wt0, 1);
     daxpy (htot_nz_A, 0.5, w1, 1, wt0, 1);
      
  	 dvmul  (htot_nz_A, ut0,1,Ut->base_h,1,Pt->base_h,1);
  	 dvvtvp (htot_nz_A, vt0,1,Vt->base_h,1,Pt->base_h,1,Pt->base_h,1);
  	 dvvtvp (htot_nz_A, wt0,1,Wt->base_h,1,Pt->base_h,1,Pt->base_h,1);

     dvadd(htot_nz_A, ut1, 1, Pt->base_h, 1, Pt->base_h, 1);

      
     Ut->Grad_h(Ut->base_h,ut1,null,null,trip);
     Vt->Grad_h(Vt->base_h,null,vt1,null,trip);
     Wt->Trans(Ut, p_to_f); //go to (Q,F)-space 
     Ut->Grad_z(Wt); Wt->Set_state('p'); //in (Q,F)-space
     Wt->Trans(Wt, f_to_p); //go to (Q,P)-space 
//     Ut->Grad_d(ut1, NULL, NULL, 'x');
//     Vt->Grad_d(NULL,vt1, NULL, 'y');
//     Wt->Trans(Wt, P_to_F); //go to (Q,F)-space 
  //   Wt->Grad_z(Ut);

     dvadd(htot_nz_A, ut1, 1, vt1, 1, ut1, 1);
     dvadd(htot_nz_A, Wt->base_h, 1, ut1, 1, ut1, 1);
#ifdef SPM
     for(k = 0; k<NZ; ++k)
      for(i =0; i<T->nel; ++i)
       {
         E = C->flevels[0]->flist[i];
          for(j=0; j<E->qtot; ++j)
          {
            int offset = k*T->htot+i*E->qtot+j;
             if (Kr > 1.)
              {
//               if(C->base_h[offset]<0.5) //in fluid side
//                ut1[offset] *=-kappa/Kr;
//               else               //in solid side
//                ut1[offset] *=-kappa;
//
               ut1[offset] *=-(kappa/Kr*(1.0-C->base_h[offset])+kappa*C->base_h[offset]);
              }

             if (Kr < 1.)
              {
//               if(C->base_h[offset]>0.5) //in solid side
//                ut1[offset] *= -kappa*Kr;
//               else
//                ut1[offset] *= -kappa;
//
               ut1[offset] *=-(kappa*(1.0-C->base_h[offset])+kappa*Kr*C->base_h[offset]);
              }

         }
//         int group_id =  E->group_id; 
//
//         if(group_id > 1)
//           for(j=0; j<E->qtot; ++j)
//            E->h[0][j] -= qs; 
       }
#else
     dsmul(htot_nz_A, -1./RePr, ut1, 1, ut1, 1); 
#endif


     dvadd(htot_nz, Pt->base_h, 1,ut1, 1, Pt->base_h, 1); //now Pt->base_h stores all the residual
 
   #ifdef CONJUGATE_HEAT
     double       qs = dparam("DHEAT_SOURCE");
     for(k = 0; k<NZ; ++k)
      for(i =0; i<T->nel; ++i)
       {
         E = C->flevels[0]->flist[i];
         F = Pt->flevels[0]->flist[i];
//         int group_id =  E->group_id; 
//
//         if(group_id > 1)
//           for(j=0; j<E->qtot; ++j)
//            E->h[0][j] -= qs; 
           for(j=0; j<E->qtot; ++j)
           {
            int offset = k*T->htot+i*E->qtot+j;
//            if(E->h[0][j]>0.5)
//               F->h[0][j] -= qs; 
            if(C->base_h[offset]>0.5) //in solid side
              Pt->base_h[offset] -=qs;
           }
       }
    #endif


     dsadd(htot_nz_A, -temperature_mean, Tk1, 1, Ut->base_h, 1);

     dvmul(htot_nz_A, Ut->base_h, 1, Pt->base_h, 1, Pt->base_h, 1); // Pt->base_h has the entropy residual;

     double entropy_residual;
     double entropy_viscosity,max_viscosity = -1e6,element_max_viscosity;
     double element_diameter,qll_diameter;
     double first_order_viscosity = 0.;

     double entropy_variation = temperature_entropy_variation;
     double weight =0., volume = 0.;
      
#ifdef DYNAMIC_ALPHA
//      dvadd  (htot_nz_A, Tk1, 1, Tk1,  1, Ut->base_h, 1);   //
//      dsmul  (htot_nz_A, 0.5, Ut->base_h,  1, Ut->base_h, 1);  //Ut->base_h has 0.5*(T*T) in (Q,P)-space
      dcopy(htot_nz_A, Tk1, 1, Ut->base_h, 1);  //original T
      SetReducedPhys(Ut, Vt, Wt); //Ut input
      dcopy(htot_nz_A, Vt->base_h, 1, vt1, 1);  //reduced T

      dcopy(htot_nz_A, Tk1, 1, ut1, 1);  //original T*T

#endif
//     dcopy(htot, u0, 1, Ut->base_h, 1);
//     dcopy(htot, v0, 1, Vt->base_h, 1);
//     dcopy(htot, w0, 1, Wt->base_h, 1);
    dsmul(htot_nz_A, BBeta[0], u1,  1, Ut->base_h, 1);
    daxpy(htot_nz_A, BBeta[1], u0,  1, Ut->base_h, 1); // 
    dsmul(htot_nz_A, BBeta[0], v1,  1, Vt->base_h, 1);
    daxpy(htot_nz_A, BBeta[1], v0,  1, Vt->base_h, 1); // 
    dsmul(htot_nz_A, BBeta[0], w1,  1, Wt->base_h, 1);
    daxpy(htot_nz_A, BBeta[1], w0,  1, Wt->base_h, 1); // 

//     dcopy(htot, Tt0, 1, T->base_h, 1);
//     dcopy(hjtot, T->base_hj, 1, ts, 1);
  if(step == 1)
   {
    escales1 = dvector(0, Ut->nel-1);
    escales2 = dvector(0, Ut->nel-1);
    dzero(Ut->nel, escales1, 1);
     dzero(Ut->nel, escales2, 1);
  
    for(E=Ut->fhead;E;E=E->next){
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
      
      dxa = 0.5*(za[1]-za[0]);
      dxb = 0.5*(zb[1]-zb[0]);

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

//      fprintf(stdout,"escales[%d] = %lf dxa = %lf dxb = %lf sca = %lf scb = %lf scc = %lf  scd = %lf \n",E->id,escales[E->id],dxa,dxb, sca,scb,scc,scd);
     }
    }

     } //step ==1
    step++;

    for (k = 0; k < NZ; k++)
     for(int eid =0; eid<Pt->nel; ++eid)
       {
	       E   = Ut->flevels[0]->flist[eid];
	       F   = Vt->flevels[0]->flist[eid];
	       G   = Wt->flevels[0]->flist[eid];
	       R   = Pt->flevels[0]->flist[eid];

	       E->GetZW(&za, &wa, &zb, &wb, &zc, &wc);
         qa = E->qa;   qb = E->qb;

	       E->coord(&X);

         element_diameter = escales1[eid];
         qll_diameter = escales2[eid];
//         qll_diameter = escales1[eid];
  
          int offset = Ut->htot*k+eid*E->qtot; 
          double local_max_velocity = -1e8; 
         for ( i=0; i<qa; i++ )
         for ( j=0; j<qb; j++ )
            {
               m = i*qb + j;
              
            local_max_velocity =  max( sqrt(Ut->base_h[offset+m]*Ut->base_h[offset+m]
                                           +Vt->base_h[offset+m]*Vt->base_h[offset+m]
                                           +Wt->base_h[offset+m]*Wt->base_h[offset+m]),
                                              local_max_velocity);
            }

        first_order_viscosity =  beta* local_max_velocity* qll_diameter;
//        first_order_viscosity =  beta* local_max_velocity* element_diameter;
#ifdef DYNAMIC_ALPHA

        double numerator = 0.;
        double denominator = 0.;
        for ( int i=0; i<qa; i++ )
         for ( int j=0; j<qb; j++ )
          {
            m = i*qb + j;

            numerator += (ut1[offset+m]-vt1[offset+m])
                        *(ut1[offset+m]-vt1[offset+m]);
            denominator += ut1[offset+m]*ut1[offset+m];

          }
        double elmtSensor = sqrt(numerator/denominator);

        max_sensor = max(max_sensor,elmtSensor);
#endif

        volume = 0.;
        element_max_viscosity = 0.;
        for ( i=0; i<qa; i++ )
        for ( j=0; j<qb; j++ )
         {    
          m = i*qb + j;
          
          entropy_residual = fabs( Pt->base_h[offset+m]);

#ifdef DYNAMIC_ALPHA
          alpha *= elmtSensor;
#endif
	        entropy_viscosity = (alpha * qll_diameter * qll_diameter *
					                      entropy_residual / entropy_variation);

          element_max_viscosity = max(element_max_viscosity,entropy_viscosity);

//         if(entropy_residual>1e-5)
//          ROOTONLY fprintf(stdout, "residual = %lf  entropy_variation = %lf  entropy_viscosity = %lf qll_diameter = %lf \n", entropy_residual, entropy_variation, entropy_viscosity, qll_diameter);

         }
        
        double element_entropy_viscosity = min(element_max_viscosity,first_order_viscosity);
        
         double damping_factor  = 1.;

        if(dparam("DWALL_ADAPTING"))
            {
            double elmt_cx = 0.; //find the centriod of an element
            double elmt_cy = 0.;
            for(int v = 0; v < E->Nverts; ++v)
            {
             elmt_cx += E->vert[v].x; 
             elmt_cy += E->vert[v].y; 
            }
             elmt_cx /= (E->Nverts);
             elmt_cy /= (E->Nverts);

           if(dparam("CHANNEL"))
            y_plus = Re_tau*(1.-fabs(elmt_cy));

           if(dparam("PIPE"))
            y_plus = Re_tau*(0.5-sqrt(elmt_cx*elmt_cx+elmt_cy*elmt_cy));

           if(dparam("CYLINDER"))
            {
             dist = sqrt(pow(elmt_cx,2.)+pow(elmt_cy,2))-0.5;
             Rex = dist*Re;
             Cf = pow(2.*log10(Rex)-0.65,-2.3);
             u_star = sqrt(0.5*Cf);
             y_plus = Re*u_star*dist;
            }
          
            if(y_plus<damping_length)
             damping_factor = pow(1.-exp(-y_plus/A),2.);
            }

        element_entropy_viscosity *= damping_factor;

        max_viscosity = max(max_viscosity,element_entropy_viscosity);

        for ( i=0; i<qa; i++ )
        for ( j=0; j<qb; j++ )
         {    
           m = i*qb + j;
//           lambda->p[m] = element_entropy_viscosity;
           temp_visc_tmp[offset+m] = element_entropy_viscosity;
         }
    
     }

      dcopy(Ut->htot,temp_visc_ave,1,temp_visc_bak,1);
      dzero(Ut->htot, temp_visc_ave, 1);
#if 0
      for(int i=0; i<nz; ++i)
       for(int j=0; j<Ut->htot; ++j)
          temp_visc_ave[j] += temp_visc_tmp[i*Ut->htot+j];

      MPI_Allreduce (temp_visc_ave, temp_visc_tmp, Ut->htot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dcopy(Ut->htot,temp_visc_tmp,1,temp_visc_ave,1);

      int tmp = 0;
      int total_nz;
      MPI_Allreduce (&nz, &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      icopy(1,&tmp,1,&total_nz,1);

      for(int i=0; i<Ut->htot; ++i)
         temp_visc_ave[i] /= (double)total_nz;
#else  
     for(int j=0; j<Ut->htot; j++)
      {
        temp_visc_ave[j] = temp_visc_tmp[j];
        for(int i=1; i<NZ; ++i)
           temp_visc_ave[j] = max(temp_visc_ave[j],temp_visc_tmp[i*Ut->htot+j]);
      }

       MPI_Allreduce (temp_visc_ave, temp_visc_tmp, Ut->htot, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
       dcopy(Ut->htot,temp_visc_tmp,1,temp_visc_ave,1);
#endif
//       for(int i=0; i<Ut->htot; ++i)
//         temp_visc_ave[i] =2./3.*temp_visc_ave[i]+1./3.*temp_visc_bak[i];
       

       for(int i=0; i<Ut->htot; ++i)
         max_viscosity = max(max_viscosity,temp_visc_ave[i]);

#ifdef SPM
       for(int i=0; i<C->nel; ++i)
        {
             E = C->flevels[0]->flist[i];
            for(j=0; j<E->qtot; ++j)
             {
               int offset = i*E->qtot+j;
//              if(E->h[0][j]>0.5)
//                temp_visc_ave[offset] += 1./(RePr*kr)-1./RePr; 
        //      if(E->h[0][j]<=0.5)
        //        temp_visc_ave[offset] += 1./(RePr)-kr/RePr; 
             if (Kr > 1.)
              {
//               if(E->h[0][j]<0.5)
//                temp_visc_ave[offset] +=(-kappa+kappa/Kr);
                temp_visc_ave[offset] +=(-kappa+kappa/Kr)*(1.0-E->h[0][j]);
              }
             if (Kr < 1.)
              {
//               if(E->h[0][j]>0.5)
//                temp_visc_ave[offset] +=(-kappa+kappa*Kr);
                temp_visc_ave[offset] +=(-kappa+kappa*Kr)*E->h[0][j];
              }

             }
        }
#endif


#if 0       
       dcopy(htot, Tk0, 1, Ut->base_h, 1);
       
       Ut->Trans(Ut, P_to_F); // transformed to Fourier space
       
       dcopy(htot, Ut->base_h, 1, ut0, 1);
//
       
   for (int k = 0; k < nz; k++)
//   if(k != 1 || parid(k) != 0)//not sure about this !!!
//   if(parid(k) != 1)//not sure about this !!!
    {

     for(int eid =0; eid<Ut->nel; ++eid)
//     for (int k = 0; k < nz; k++)
     {
	       E = Ut->flevels[k]->flist[eid];

         dcopy(E->qtot, temp_visc_ave+eid*E->qtot, 1, lambda->p, 1);
//         dscal(E->qtot, RePr/kr, lambda->p, 1);
         dscal(E->qtot, 1./kappa, lambda->p, 1);


         E->HelmHoltz(lambda);
      }

    }

    dcopy(hjtot, Ut->base_hj, 1, us0, 1); // backup ?

    dcopy(htot, ut0, 1, Ut->base_h, 1);


   for (int k = 0; k < nz; k++)
//   for(int eid =0; eid<U->nel; ++eid)
//    if(k != 1 || parid(k) != 0)
//   if(parid(k) != 1)//not sure about this !!!
    {
//      for(Et=Ut->flevels[k]->fhead, Ft=Vt->flevels[k]->fhead, Gt=Wt->flevels[k]->fhead,eid=0; Et;Et=Et->next,Ft=Ft->next,Gt=Gt->next,eid++)
     for(int eid =0; eid<Ut->nel; ++eid)
//     for (int k = 0; k < nz; k++)
      {
	       E = Ut->flevels[k]->flist[eid];

         dvmul  (E->qtot, temp_visc_ave+eid*E->qtot, 1, E->h[0], 1, E->h[0], 1);

//         dscal  (E->qtot, RePr/kr, E->h[0], 1);
         dscal  (E->qtot, 1./kappa, E->h[0], 1);
      }
      
    }


      Ut->Iprod(Ut);

    for (int k = 0; k < nz; k=k+2)
//     if(k != 1 || parid(k) != 0)
//   if(parid(k) != 1)//not sure about this !!!
     {
	     double beta = Beta(k)*Beta(k);
       dsmul (2*Ut->hjtot,  beta, Ut->flevels[k]->base_hj,  1, Ut->flevels[k]->base_hj, 1); 
       
     }

    dvadd(hjtot, Ut->base_hj, 1, us0, 1, Ut->base_hj, 1); 

    dvadd(hjtot, Tf->base_hj, 1, Ut->base_hj,  1, Tf->base_hj, 1); // add the EVM term
#else
    dsmul(htot_nz_A, BBeta[0], Tk1,  1, Pt->base_h, 1);
    daxpy(htot_nz_A, BBeta[1], Tk0,  1, Pt->base_h, 1); // 
    
//    Pt->Trans(Wt, P_to_F); //go to (Q,F)-space 
    Pt->Trans(Wt, p_to_f); //go to (Q,F)-space 
    
    Wt->Grad_d(ut1, vt1, wt1, 'a'); // NOTE: 'a', not trip in this line, because in (Q,F)-space

    for (int k = 0; k < nz; k++)
      for(int j=0; j<Ut->htot; ++j)
        {
           ut1[k*Ut->htot+j] *= temp_visc_ave[j];
           vt1[k*Ut->htot+j] *= temp_visc_ave[j];
           wt1[k*Ut->htot+j] *= temp_visc_ave[j];
        }
       
     FourierList_Iprod_div(Pt,	ut1, vt1, wt1, Vt->base_hj);
     daxpy(hjtot_nz, RePr, Pt->base_hj, 1, Tf->base_hj, 1); 
     Pt->Set_state('p');
#endif
    

 //   dcopy(hjtot, ts, 1, T->base_hj, 1); // copy back

    gdmax(&max_viscosity,1,&entropy_viscosity);

#ifdef DYNAMIC_ALPHA
//    MPI_Allreduce (&max_sensor, &entropy_viscosity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
 //   dcopy(1,&entropy_viscosity,1,&max_sensor,1);
    gdmax(&max_sensor,1,&entropy_viscosity);
#endif

#ifdef DYNAMIC_ALPHA
    ROOTONLY fprintf(stdout,"Maximum  Temperature Sensor ......... ... %g \n", max_sensor);
#endif
    ROOTONLY fprintf(stdout,"Maximum temperature entropy viscosity ... %g \n", max_viscosity);
    ROOTONLY fprintf(stdout,"Temperature range ................[%g, %g] \n", temperature_range[1],temperature_range[0]);

  	free(X.x);
	  free(X.y);
    free(lambda->p);
    free(lambda);
//    free(entropy_viscosities);
}

#endif



