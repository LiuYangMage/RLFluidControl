
#include <time.h>
#include "nektarF.h"
#include <stdio.h>

#ifdef SPM
#include "C_SPM.h"
#endif

static int init=1; /* externals */
static int nqp=0;
static double *nusselt_vector,*coord_vector;
static double *coord_x;
static double *coord_y;
static int *map_vector;
void compute_local_Nusselt(Domain *omega);
void set_SurGeofac(Bndry *Pbc, Bndry *Ubc);
void output_temperature(Domain *omega, double time);
void setup_initial_temperature(Domain *omega);

static double * points_locations;
static double * data;
static int temp_init = 1;

#ifdef NATURAL_CONVECTION
void add_buoyancy_force(Domain *omega);

void add_buoyancy_force(Domain *omega)
 {
  Element_List *U   =  omega->U,   *V   =  omega->V,   *W   = omega->W,
               *Uf  =  omega->Uf,  *Vf  =  omega->Vf,  *Wf  = omega->Wf,
               *T   =  omega->T;
  #ifdef SPM
  Element_List *C = omega->CONCENTR;
  #endif

  register int i;
  int          nq = U->htot*U->nz;
  double       Ra = dparam("DRAYLEIGH");
  double       Pr = dparam("PRANDTL");

//  if(Ra == 0.)
//   {
//    fprintf(stderr," Rayleigh number is zero... \n");
//    exit(-1);
//   }
  
  if (option("dealias")) 
    nq *= 3/2;

  switch (iparam("DIRGRAVITY"))
   {
      case 0:
#ifdef SPM
       for(i=0; i<nq; ++i) Uf->base_h[i] += Ra*Pr*T->base_h[i]*(1.- C->base_h[i]);
#else
       for(i=0; i<nq; ++i) Uf->base_h[i] += Ra*Pr*T->base_h[i];
#endif
      //   dvadd (nq, T->base_h, 1, Uf->base_h, 1, Uf->base_h,1);
      break;
      case 1:
//         dvadd (nq, T->base_h, 1, Vf->base_h, 1, Vf->base_h,1);
#ifdef SPM
       for(i=0; i<nq; ++i) Vf->base_h[i] += Ra*Pr*T->base_h[i]*(1.- C->base_h[i]);
#else
       for(i=0; i<nq; ++i) Vf->base_h[i] += Ra*Pr*T->base_h[i];
#endif
      break;
      case 2:
//         dvadd (nq, T->base_h, 1, Wf->base_h, 1, Wf->base_h,1);
#ifdef SPM
       for(i=0; i<nq; ++i) Wf->base_h[i] += Ra*Pr*T->base_h[i]*(1.- C->base_h[i]);
#else
       for(i=0; i<nq; ++i) Wf->base_h[i] += Ra*Pr*T->base_h[i];
#endif
      break;
      default:
         fprintf(stderr,"parameter DIR_GRAVITY = %d Wrong value of direction of gravity \n",iparam("DIR_GRAVITY"));
         exit(-1);
      break;
   }


 }
#endif

#ifdef CONJUGATE_HEAT
#ifdef SPM
void add_heat_source(Domain *omega, C_SPM *SPMv);
#else
void add_heat_sorce(Domain *omega);
#endif

void regenerate_temperature_system(Domain *omega) ;

#ifdef SPM
void add_heat_source(Domain *omega, C_SPM *SPMv)
 {
  Element_List *C   =  omega->CONCENTR;
#else
void add_heat_source(Domain *omega)
 {
#endif
  Element_List  *T   =  omega->T, *Tf = omega->Tf;
  Element      *el,*fl;
  register int i,j,k;
  int          nz = T->nz;
  int          nq = T->htot*T->nz;
  double       qs = dparam("DHEAT_SOURCE");

  if(qs == 0.)
   {
    fprintf(stderr," Heat source is zero... \n");
    exit(-1);
   }
   
  //heat generation only exists in solid...
  //
//  for(i=0; i<nq; ++i) Tf->base_h[i] += C->base_h[i]*qs;
#ifdef SPM
  double *ind_func;
  int part_ID;   
  for (part_ID = 0; part_ID < SPMv[0].num_partls; ++part_ID){
     ind_func = SPMv[0].gaussian[part_ID];
#endif
   for(k = 0; k<nz; ++k)
    {
      int nid = 0,np=0;
     for(i =0; i<T->nel; ++i)
      {
       el = Tf->flevels[k]->flist[i];
//       int group_id =  el->group_id; 
//
//        if(group_id > 1)
//         {
//          nid++;
//          int qtot = el->qa*el->qb;
//          for(j=0; j<qtot; ++j)
//           {
//            np++;
//            el->h[0][j] += qs; 
//           }
//         }
#if 0
       fl = C->flevels[k]->flist[i];
          for(j=0; j<el->qtot; ++j)
            if(fl->h[0][j]>0.5)
              el->h[0][j] += qs;
#endif
         if(SPMv[0].mobile[part_ID]) //has heat source
          for(j=0; j<el->qtot; ++j)
            if(ind_func[k*Tf->htot+i*el->qtot+j]>0.5) 
              el->h[0][j] += qs;
              
       }
//    fprintf(stdout,"3 Number of solid elements = %d quad points = %d \n",nid,np);
     }
#ifdef SPM
  }
#endif

 }

void setup_initial_temperature(Domain *omega)
 {
    Element_List  *T   =  omega->T, *CONTRA = omega->CONCENTR;
    register int i,j,k;
    int          nz = T->nz;
    int          nq = T->htot*T->nz;
    double       T_init = dparam("DCABLE_TEMP_INIT");

    for(i=0; i<nq; ++i)
      T->base_h[i] = T_init*CONTRA->base_h[i];

    T->Set_state('p');
    T->Trans(T, P_to_F); 
    T->Trans(T, Q_to_J);
    T->Set_state('t');
 }

void regenerate_temperature_system(Domain *omega) {
  register int i,j,k;
  Element_List *C   =  omega->CONCENTR, *T   =  omega->T, *Tf = omega->Tf;
  double dt = dparam("DELT");
  double Re = 1.0 / dparam("KINVIS");
  double Pr = dparam("PRANDTL");
  double kr = dparam("RATIO_KAPPA") == 0? 1.:dparam("RATIO_KAPPA");

  for(k = 0; k < T->nz; ++k){
    omega->Tsys[k]->lambda = (Metric*) calloc(T->nel, sizeof(Metric));
  
    int nid = 0;
    for(i=0;i<T->nel;++i)  {
     double Pr_tmp = Pr;

     bool solid = true;
     Element      *el = C->flevels[k]->flist[i];
     Element      *fl = T->flevels[k]->flist[i];
     Element      *gl = Tf->flevels[k]->flist[i];
     
     int qtot = el->qa*el->qb;
	   for (int j = 0; j < qtot; j++)
      if( el->h[0][j] < 0.5)
       {
           solid = false;
           break;
       }

      if(solid)
      {
       nid++;

       Pr_tmp = Pr/kr;
       fl->group_id = 2; 
       gl->group_id = 2; 
      }

 	  omega->Tsys[k]->lambda[i].d = Beta(k)*Beta(k) + Re*Pr_tmp*getgamma(1)/dt;
   }
    
//    fprintf(stdout,"1 Number of solid elements = %d \n",nid);

    ROOT fprintf(stdout,"Level: %d \n", k);

    if(!(parid(k) % 2)) {
      ROOT fprintf(stdout,"SPM Generating temperature system [."); 
      ROOT fflush(stdout);
      GenMat (T->flevels[k],omega->Tbc[k],omega->Tsys[k],omega->Tsys[k]->lambda,Helm);
      ROOT fprintf(stdout,"]\n");
    }
    else{
      omega->Tsys[k]->Gmat = omega->Tsys[k-1]->Gmat; 
      omega->Tsys[k]->Pmat = omega->Tsys[k-1]->Pmat;
      omega->Tsys[k]->rslv = omega->Tsys[k-1]->rslv;
    }
  }

}

#endif

void compute_local_nusselt(Domain *omega, double time)
 {
  Element_List *U   =  omega->U,   *V   =  omega->V,   *W   = omega->W,
               *T   =  omega->T;
  Element *eT;
  register int i,j,k;
  int      nz = option("NZ"), nztot = option("NZTOT");
  int      nprocs = option("NPROCS");
  int      eid,face,lmax=0, trip = 0;
  int      count = 0;
  Basis    *b;

  double   radius = dparam("RADIUS");
  double   T_infty = dparam("DTEMP_INTY");

  if(radius <1e-10)
   {
     fprintf(stderr,"Need to prescribe a value to parameter RADIUS \n");
     exit(1);
   }

  Bndry   *UB = omega->Ubc[0],
          *PB = omega->Pbc[0],
          *TB = omega->Tbc[0];

  int      qa = T->flevels[0]->flist[0]->qa;
  int      qb = T->flevels[0]->flist[0]->qb;
  
//  T->Set_state('p');
//  T->Trans(T, J_to_Q);
//  T->Trans(T, F_to_P); // here we are in phsyical/quadrature space
  

  double  *w;
  double  *wk = dvector(0,max(LGmax+1,QGmax)*QGmax-1);
  double  **D = dmatrix(0,2,0,QGmax*QGmax-1);
  double  *tx = D[0], *ty = D[1], *tt = D[2];
  
  /* write out header */
  if(init){

  Coord X;
	X.x  = dvector(0, qa*qb-1);
	X.y  = dvector(0, qa*qb-1);
  
  /* need to set up surface geometric factors in T from P */
  for(k = 0; k < nz; ++k)
    set_SurGeofac(omega->Pbc[k],omega->Ubc[k]);
  // count how many points on the surface
  count = 0;
  nqp = 0;
//  for(UB=omega->Ubc[0],TB=omega->Tbc[0];UB; UB = UB->next,TB=TB->next)
  for(UB=omega->Ubc[0];UB; UB = UB->next)
    for(TB = omega->Tbc[0]; TB; TB = TB->next)
	    if((TB->elmt->id == UB->elmt->id)&&(TB->face == UB->face))
        if( ((TB->type == 'V') ||(TB->type == 'F')) && (UB->type == 'W') )
         {
           nqp += qa;
           count ++;
         }


    ROOT fprintf(stdout,"%d Surfaces are found \n",count);

    nusselt_vector = dvector(0,nqp-1);
    map_vector     = ivector(0,nqp-1);
    coord_vector   = dvector(0,nqp-1);
    coord_x        = dvector(0,nqp-1);
    coord_y        = dvector(0,nqp-1);

    count = 0;
//  for(UB=omega->Ubc[0];UB; UB = UB->next)
//    if(UB->type == 'W'){
  for(UB=omega->Ubc[0];UB; UB = UB->next)
   for(TB = omega->Tbc[0]; TB; TB = TB->next)
	  if((TB->elmt->id == UB->elmt->id)&&(TB->face == UB->face))
     if( ((TB->type == 'V') ||(TB->type == 'F')) && (UB->type == 'W') ) {
       face = UB->face;
       eid  = UB->elmt->id;

      for(k = 0; k < 1; ++k){
	     eT = T->flevels[k]->flist[eid];

	     eT->coord(&X);

	     eT->GetFace(X.x ,face,wk);
	     eT->InterpToFace1(face,wk,X.x);
	     eT->GetFace(X.y ,face,wk);
	     eT->InterpToFace1(face,wk,X.y);
	     
       for(i = 0; i < qa; ++i){
         coord_vector[count*qa+i] = atan2(X.y[i],X.x[i]);
         coord_x[count*qa+i] = X.x[i]/radius;//normal outward
         coord_y[count*qa+i] = X.y[i]/radius;//normal outward
         map_vector[count*qa+i] = count*qa+i;
//         fprintf(stderr," theta = %g r %g \n",coord_vector[count*qa+i],sqrt(X.x[i]*X.x[i]+X.y[i]*X.y[i]));
        }
      }
      count++;
    }
//simple bubble sorting map_vector

     double dswap = 0.;
     double xswap = 0.;
     double yswap = 0.;
     int  iswap = 0;
     for(i=0; i<nqp; ++i)
       for(j=i+1; j<nqp; ++j)
        if(coord_vector[i]>coord_vector[j])
        {
           dswap = coord_vector[i];
           iswap = map_vector[i];

           xswap = coord_x[i];
           yswap = coord_y[i];

           coord_vector[i] = coord_vector[j];
           coord_x[i]      = coord_x[j];
           coord_y[i]      = coord_y[j];
           map_vector[i]   = map_vector[j];

           coord_vector[j] = dswap;
           map_vector[j]   = iswap;
           coord_x[j]      = xswap;
           coord_y[j]      = yswap;

        }

    init = 0;
    free(X.x);free(X.y);
  }

  dzero(nqp,nusselt_vector,1);

  count = 0;
  double mean_nu = 0.,area = 0.;
//  for(UB=omega->Ubc[0];UB; UB = UB->next)
//    if(UB->type == 'W'){
  for(UB=omega->Ubc[0];UB; UB = UB->next)
   for(TB = omega->Tbc[0]; TB; TB = TB->next)
	  if((TB->elmt->id == UB->elmt->id)&&(TB->face == UB->face))
     if( ((TB->type == 'V') ||(TB->type == 'F')) && (UB->type == 'W') ) {
      face = UB->face;
      eid  = UB->elmt->id;

      if(T->flist[eid]->lmax != lmax){
	       b    = T->flist[eid]->getbasis();
	       lmax = T->flist[eid]->lmax;
	       getzw(qa,&w,&w,'a');
      }

      for(k = 0; k < 1; ++k){
	     eT = T->flevels[k]->flist[eid];

	     eT->Jbwd(eT,b);  eT->state = 't';
       eT->Grad_d(tx, ty, 0, 'a');
	     
	     for(i = 0; i < 2; ++i){
	        eT->GetFace(D[i] ,face,wk);
	        eT->InterpToFace1(face,wk,D[i]);
	      }

	     eT->GetFace(eT->h[0],face,wk);
	     eT->InterpToFace1(face,wk,tt);
	     
       if(eT->curvX)
	      for(i = 0; i < qa; ++i){
   //
            double local_nu = tx[i]*UB->nx.p[i]+UB->ny.p[i]*ty[i];
            double delta_T = tt[i]-T_infty;
//            double local_nu = fabs(tx[i]*coord_x[map_vector[count*qa+i]]+ty[i]*coord_y[map_vector[count*qa+i]]);
   //fprintf(stdout,"nu = %g tx = %g ty = %g nx = %g ny = %g count = %d \n",local_nu,tx[i],ty[i],TB->nx.p[i],TB->ny.p[i],count);
	          nusselt_vector[map_vector[count*qa+i]] = local_nu/delta_T;
            mean_nu += local_nu*UB->sjac.p[i]*w[i];
            area += UB->sjac.p[i]*w[i];
           }
	    else
	       for(i = 0; i < qa; ++i){
            double local_nu = tx[i]*UB->nx.d+UB->ny.d*ty[i];
	          nusselt_vector[map_vector[count*qa+i]] = local_nu;
            mean_nu += local_nu*UB->sjac.d*w[i];
            area += UB->sjac.p[i]*w[i];
           }
      }
        
       count++;
    }

  if(area < 1e-8)
   {
    fprintf(stderr,"Cannot find any cylinder wall ? \n");
    exit(-1);
   }
   
  mean_nu /= area;

    ROOT {

        FILE *fout1,*fout2;
	      char *buff1 = (char*) calloc(BUFSIZ,sizeof(char));
	      char *buff2 = (char*) calloc(BUFSIZ,sizeof(char));
	      sprintf(buff1, "local_nusselts.txt");
	      sprintf(buff2, "mean_nusselts.txt");
        fout1 = fopen(buff1,"aw");
        fout2 = fopen(buff2,"aw");

         for(i=0; i<nqp; ++i)
            fprintf(fout1,"%g %g  %g \n",time,nusselt_vector[i],coord_vector[i]);
         
       fprintf(fout2,"%g %g \n",time,mean_nu);

        fflush(fout1);
        fclose(fout1);
        free(buff1);

        fflush(fout2);
        fclose(fout2);
        free(buff2);

    }

  free(wk); free_dmatrix(D,0,0);
 }

void output_temperature(Domain *omega, double time) {
 #ifdef SPM
   Element  *E;
   int nz = omega->T->nz;
   int htot = omega->T->htot*omega->T->nz;
   int i,j,k,m,plane, nel, qtot,offset = 0;
   int qa,qb;
   int part_ID;   

   double   *ha,*hb,*hc,*za,*wa,*zb,*wb,*zc,*wc;
   qa = omega->T->flevels[0]->flist[0]->qa;
   qb = qa;
	 ha = dvector(0,qa-1);
	 hb = dvector(0,qb-1);

   Coord XX,*YY;
   XX.x = dvector(0,2-1); 
   XX.y = XX.x + 1;

   int data_size = omega->vSPM[0].num_partls*3;

   if(temp_init)
   {
     points_locations = dvector(0,data_size-1);
     data = dvector(0, omega->vSPM[0].num_partls-1); 

     for(k=0,plane=nz*option("PROCID"); k<1; ++k,++plane) {
       offset = plane*9; //circle particle 
      for (part_ID = 0, i=0; part_ID < omega->vSPM[0].num_partls; ++part_ID, i=i+3){
        double x0 = omega->vSPM[0].shape[part_ID].data[offset+0];
        double y0 = omega->vSPM[0].shape[part_ID].data[offset+1];

        XX.x[0] = x0;
        XX.y[0] = y0;

        int *eid ;
        Find_local_coords(omega->T,&XX,1,&eid,&YY);

        points_locations[i] = (double) eid[0];
        points_locations[i+1] = YY->x[0];
        points_locations[i+2] = YY->y[0];
      }
     }

     temp_init = 0;
   }

  omega->T->Set_state('p');
  omega->T->Trans(omega->T, J_to_Q);
  omega->T->Trans(omega->T, F_to_P); // here we are in phsyical/quadrature space

  for(int k=0;k < 1; k++)
    for(int i=0, j=0;i<data_size; i=i+3, j++)
     {
        int eid      =  (int) (points_locations[i]);

        XX.x[0] = points_locations[i+1]; 
        XX.y[0] = points_locations[i+2]; 
              
		    E = omega->T->flevels[k]->flist[eid];
	      E->GetZW(&za, &wa, &zb, &wb, &zc, &wc);

        get_point_shape_2d(E,XX.x[0],XX.y[0],ha,hb);
		    data[j] = eval_field_at_pt_2d(qa,qb, E->h[0], ha, hb);
     }

//  for(k = 0; k < omega->T->nz; ++k)
  double t_max  = -10e6;
  double area   = 0.;
  double t_mean = 0.;
  double weight = 0.;
  for(k = 0; k < 1; ++k){
     for(i =0; i< omega->T->nel; ++i)
      {
        E = omega->T->flevels[k]->flist[i];
        
        getzw(qa, &za, &wa, 'a');
        getzw(qb, &zb, &wb, 'a');

       for (int ii=0; ii<qa; ii++ )
        for (int jj=0; jj<qb; jj++ )
         {
            m = ii*qb + jj;

           if( E->curvX )
            weight = wa[ii]*wb[jj]*E->geom->jac.p[m];
           else
            weight = wa[ii]*wb[jj]*E->geom->jac.d;

           t_max = max(t_max,E->h[0][m]);

           area += weight;
           t_mean += E->h[0][m]*weight;
         }
      }
   }
     t_mean /= area;

     FILE * file_point;
	   char *buf_point = (char*) calloc(BUFSIZ,sizeof(char));
	   sprintf(buf_point, "temp.txt");
     file_point = fopen(buf_point,"aw");
      
     for(i=0; i<omega->vSPM[0].num_partls; ++i)
       fprintf(file_point," %lf  %lf  %d \n",time, data[i], i);

       fprintf(file_point," %lf  %lf  %d \n",time, t_max, omega->vSPM[0].num_partls);
       fprintf(file_point," %lf  %lf  %d \n",time, t_mean, omega->vSPM[0].num_partls+1);
     
     fflush(file_point);
     fclose(file_point);
     free(buf_point);
     
     free(XX.x);
     free(ha);
     free(hb);
 #else
     fprintf(stderr,"This function has to be working with SPM....\n");
     exit(-1);
 #endif

 }
