#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include "nektarF.h"

#include <rfftw.h>
extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;
#ifdef MAP
#include <map.h>
#endif

static double **mean_new;
static double **mean_new1;

void solve(Element_List *U, Element_List *Uf,
	   Bndry **Ubc,Bsystem **Ubsys,SolveType Stype,int step, double *scal);
void set_SurGeofac(Bndry *Pbc, Bndry *Ubc);

#ifdef MAP
void compute_rms_pressure(Domain *omega,Mapping *mapx, Mapping *mapy, int step);
#endif
void compute_cylindrical_vorticity(Domain *omega,int step);
void compute_local_vorticity_flux(Domain *omega, int step, double time);

void RMS_SetPBCs(Domain *omega);
static int pres_init=1; /* externals */
static int nqp=0;
static double *coord_vector;
static double *coord_x;
static double *coord_y;
static int *map_vector;

double compute_adjustment_forces(Domain *omega){
    register int i,j,k;
    int          qa,qb,qc,m;
    Element_List *U,*V,*W;
    Element      *E,*F,*G;
	  double       *za,*wa,*zb,*wb,*zc,*wc;
    double       weight;
    double       tmp;

    int          nz = option("NZ"), nztot = option("NZTOT");
    int          nprocs = option("NPROCS");
    int          htot = omega->U->nz*omega->U->htot;
	  
    double       dt    = dparam("DELT");

    const double u_bulk = dparam("U_BULK");
    if(u_bulk == 0.)
     {
      ROOTONLY fprintf(stderr,"Forget to set U_BULK parameter ! \n");
      exit(1);
     }

//  	Coord X;
//	  X.x  = dvector(0, QGmax*QGmax*QGmax-1);
//	  X.y  = dvector(0, QGmax*QGmax*QGmax-1);

    U = omega->U;
    V = omega->V;
    W = omega->W;
    double umean = 0.,umean_tmp =0.;
    double volume = 0.,volume_tmp =0.;
//    double *u    = dvector(0,qt-1),
//   ROOT
    for (k = 0; k < 1; k++)
    {    // ZERO plane only
   if(parid(k) == 0)
     for(int eid =0; eid<W->nel; ++eid)
      {
	      E = U->flevels[k]->flist[eid];
	      F = V->flevels[k]->flist[eid];
	      G = W->flevels[k]->flist[eid];
        
        qa = E->qa;   qb = E->qb;
        getzw(E->qa, &za, &wa, 'a');
        getzw(E->qb, &zb, &wb, 'a');

       for ( i=0; i<qa; i++ )
        for ( j=0; j<qb; j++ )
        {    
          m = i*qb + j;
          if( G->curvX )
            weight = wa[i]*wb[j]*E->geom->jac.p[m];
          else
            weight = wa[i]*wb[j]*E->geom->jac.d;

          volume += weight;
//          umean += sqrt(E->h[0][m]*E->h[0][m]+F->h[0][m]*F->h[0][m]+G->h[0][m]*G->h[0][m])*weight;
          if(dparam("CHANNEL"))
          umean += E->h[0][m]*weight;
          if(dparam("PIPE"))
          umean += G->h[0][m]*weight;

        }

      }
    }

    tmp = 0.;
//    MPI_Bcast(&umean,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Allreduce (&volume, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(1,&tmp,1,&volume,1);

    tmp = 0.;
    MPI_Allreduce (&umean, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(1,&tmp,1,&umean,1);

    int itmp = 0;
    int total_nz;
    MPI_Allreduce (&nz, &itmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    icopy(1,&itmp,1,&total_nz,1);

    umean /= volume;

//    ROOT fprintf(stdout," umean = %lf,  volume = %lf  flux = %lf  total_nz = %d \n",umean,volume,umean*volume, total_nz);
    
    double fx = (u_bulk-umean)/dt;

    const double alpha = 0.75;
    const double beta  = 0.25;

    fx = alpha * fx + beta *omega->adjust_force ;

    double ff =max( fabs(dparam("FFX")), max(fabs(dparam("FFY")),fabs(dparam("FFZ"))) );

    if(ff <1e-6)
      ff = 5;

    const double fx_min = -ff;
    const double fx_max = ff;

    fx = max(fx,fx_min);
    fx = min(fx,fx_max);

//    omega->fx = fx;
//	  free(X.x);
//	  free(X.y);
    
    return fx;
}


void  statistics_setup(Domain *omega){
	Element_List *U,*V,*W, *P;
  Element      *E;
  register int i,j,k;
  bool         FLAG;
  int          eDIM = omega->U->fhead->dim()+1;
  double      *send_yy = (double*)0;
  double      *yy    = (double*)0;
  double      **umean = (double**)0;
  double      **urms  = (double**)0;
  int          nz = option("NZ"), nztot = option("NZTOT");
  
  U = omega->U;
  V = omega->V;
  W = omega->W;
//  P = omega->P;

  int qa = omega->U->fhead->qa;
  int qtot = omega->U->fhead->qtot;
	Coord X;
	X.x  = dvector(0, qtot-1);
	X.y  = dvector(0, qtot-1);
  
  double Ny =  dparam("BOXNY");
  double Nc =  dparam("CIR_NC");
  double Nr =  dparam("CIR_NR");
  if( (Ny == 0)&&(dparam("CAHNNEL")))
   {
      ROOTONLY fprintf(stderr,"Forget to set BOX_NY parameter ! \n");
      exit(1);
   }
  if( ((Nc == 0)||(Nr == 0))&&(dparam("PIPE")))
   {
      ROOTONLY fprintf(stderr,"Forget to set CIR_NC or CIR_NR parameter ! \n");
      exit(1);
   }

 if(dparam("CHANNEL"))
  {
   int num_points = qa*(int)Ny;
   send_yy = dvector(0, num_points-1);
   dfill(num_points, 10e6, send_yy, 1);
    
   U = omega->U;

   ROOTONLY
   {
    k = 0;
    for(int eid =0; eid<U->nel; ++eid)
    {
	   E = U->flevels[0]->flist[eid];
	   E->coord(&X);
    
     for(i=0;i<E->qtot;++i)
     {
	     FLAG = false;  
//       for(j=0;j<num_points;++j)
     for(j=0;j<k;++j)
      if(fabs(X.y[i]-send_yy[j])<1e-10)
       {
          FLAG = true;
          break;
       }

       if(!FLAG)
       {
        send_yy[k] = X.y[i]; 
        k++;
       }
       }
      }

     double temp;
	   for(i=0; i<num_points; ++i)
	    for(j=i; j<num_points; ++j)
       if(send_yy[i] > send_yy[j])
       {
           temp   = send_yy[j];
           send_yy[j] = send_yy[i];
           send_yy[i] = temp;
       }
     }

   MPI_Bcast(send_yy,num_points,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
      
    k = 0;
    for(i=0; i<num_points; ++i)
      if(send_yy[i] != 10e6)
            k++;

    omega->ny = k;

      yy = dvector(0, k-1);
      umean = dmatrix(0,4,0,k-1);
      urms  = dmatrix(0,5,0,k-1);
      
      omega->yy = yy;
      omega->umean = umean;
      omega->urms  = urms;

    for(int i=0; i<k; ++i) {
        omega->yy[i] = send_yy[i];

    //    ROOTONLY fprintf(stderr,"yy[%d] = %lf  ny = %d \n",i,yy[i],k);
     }
   }
 

	free(X.x);
	free(X.y);
//  free(send_yy);
}

void statistics_output(Domain *omega, int step, double time){
  register int i,j,k;
	Element_List *U,*V,*W, *P;
  Element      *E,*F,*G,*H;
  int          eDIM = omega->U->fhead->dim()+1;
  int          ny = omega->ny;
  int          nz = option("NZ"), nztot = option("NZTOT");
  double       **umean = omega->umean;
  double       **urms  = omega->urms;
  double       *yy     = omega->yy;
  double       *visc   = omega->visc_ave;
  double       *xz,*xz_tmp;
  double       **umean_tmp;

  int qtot = omega->U->fhead->qtot;
	Coord X;
	X.x  = dvector(0, qtot-1);
	X.y  = dvector(0, qtot-1);

	xz_tmp  = dvector(0, ny-1);
	xz      = dvector(0, ny-1);
  dfill(ny, 0, xz_tmp, 1);
  dfill(ny, 0, xz,     1);

  umean_tmp = dmatrix(0,4,0,ny-1);

  for(int d=0; d<4; ++d)
     dzero(ny, umean[d], 1);

  for(int d=0; d<5; ++d)
   {
     dzero(ny, urms[d],  1);
     dzero(ny, umean_tmp[d], 1);
   }

//  U = omega->U;
//  V = omega->V;
//  W = omega->W;
//  P = omega->P;

  omega->Ut->Trans(omega->U->htot, omega->U->base_h, omega->Ut->base_h, F_to_P); // in (Q,P)-space
  omega->Vt->Trans(omega->V->htot, omega->V->base_h, omega->Vt->base_h, F_to_P); // in (Q,P)-space
  omega->Vt->Trans(omega->W->htot, omega->W->base_h, omega->Wt->base_h, F_to_P); // in (Q,P)-space
  omega->Pt->Trans(omega->P->htot, omega->P->base_h, omega->Pt->base_h, F_to_P); // in (Q,P)-space

  U = omega->Ut;
  V = omega->Vt;
  W = omega->Wt;
  P = omega->Pt;

// compute mean velocity
//ROOT
 for (k = 0; k < nz; k++)
 {    // ZERO plane only
//  if(parid(k) == 0)
    for(int eid =0; eid<W->nel; ++eid)
     {
	     E = U->flevels[k]->flist[eid];
	     F = V->flevels[k]->flist[eid];
	     G = W->flevels[k]->flist[eid];
	     H = P->flevels[k]->flist[eid];
	     
       E->coord(&X);
       for(i=0;i<E->qtot;++i)
        {
        for(j=0;j<ny;++j)
          if(fabs(X.y[i]-yy[j])<1e-10)
           {
            umean[0][j] += E->h[0][i];
            umean[1][j] += F->h[0][i];
            umean[2][j] += G->h[0][i];
//            umean[3][j] += visc[E->qtot*eid+i];
            umean[3][j] += H->h[0][i];

            urms[0][j] += E->h[0][i]*E->h[0][i];
            urms[1][j] += F->h[0][i]*F->h[0][i];
            urms[2][j] += G->h[0][i]*G->h[0][i];
            urms[3][j] += E->h[0][i]*F->h[0][i];
            urms[4][j] += H->h[0][i]*H->h[0][i];
            xz[j]        = xz[j]+1.;
           }
        }
     }
  }
  
 for(int d=0; d<4; ++d)
  {
    MPI_Allreduce (umean[d], umean_tmp[d], ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(ny,umean_tmp[d],1,umean[d],1);
  }


    MPI_Allreduce (xz, xz_tmp, ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(ny,xz_tmp,1,xz,1);

//    ROOT
//   for(j=0;j<ny;++j)
//      fprintf(stdout,"xz[%d] = %lf \n",j,xz[j]);

    for(j=0;j<ny;++j)
     {
         umean[0][j] /= xz[j];
         umean[1][j] /= xz[j];
         umean[2][j] /= xz[j];
         umean[3][j] /= xz[j];
     }


/*
    for(int d=0; d<eDIM+1; ++d)
        dzero(ny, umean_tmp[d], 1);

    dzero(ny, xz, 1);
// compute rms velocity
//ROOT
 for (k = 0; k < nz; k++)
 {    // ZERO plane only
//  if(parid(k) == 0)
    for(int eid =0; eid<W->nel; ++eid)
     {
	     E = U->flevels[k]->flist[eid];
	     F = V->flevels[k]->flist[eid];
	     G = W->flevels[k]->flist[eid];
	   
       E->coord(&X);
      for(i=0;i<E->qtot;++i)
      {
        for(j=0;j<ny;++j)
          if(fabs(X.y[i]-yy[j])<1e-10)
           {
            urms[0][j] += E->h[0][i]*E->h[0][i];
            urms[1][j] += F->h[0][i]*F->h[0][i];
            urms[2][j] += G->h[0][i]*G->h[0][i];
            urms[3][j] += E->h[0][i]*F->h[0][i];
            xz[j]++;
           }
       }
     }
 }
*/

  for(int d=0; d<5; ++d)
   {
    MPI_Allreduce (urms[d], umean_tmp[d], ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dcopy(ny,umean_tmp[d],1,urms[d],1);
   }

//    MPI_Allreduce (xz, xz_tmp, ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    dcopy(ny,xz_tmp,1,xz,1);

    for(j=0;j<ny;++j)
     {
         urms[0][j] /= xz[j];
         urms[1][j] /= xz[j];
         urms[2][j] /= xz[j];
         urms[3][j] /= xz[j];
         urms[4][j] /= xz[j];
     }

    for(j=0;j<ny;++j)
     {
         urms[0][j] -= umean[0][j]*umean[0][j];
         urms[1][j] -= umean[1][j]*umean[1][j];
         urms[2][j] -= umean[2][j]*umean[2][j];
         urms[3][j] -= umean[0][j]*umean[1][j];
         urms[4][j] -= umean[3][j]*umean[3][j];
     }

  ROOTONLY
  {
//      fprintf(stderr, "%d \n",ny);
   FILE *umean_out;
   char *buf = (char*) calloc(BUFSIZ,sizeof(char));
   sprintf(buf, "umean.txt",step);
   umean_out = fopen(buf,"aw");

   fprintf(umean_out,"y  umean  vmean  wmean  pmean  uu_rms  vv_rms  ww_rms  uv_rms  prms\n");
   for(int j=0; j<ceil(ny/2.); ++j)
   fprintf(umean_out,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
        1.-fabs(yy[j]), 0.5*(umean[0][j]+umean[0][ny-1-j]),
                        0.5*(umean[1][j]+umean[1][ny-1-j]),
                        0.5*(umean[2][j]+umean[2][ny-1-j]),
                        0.5*(umean[3][j]+umean[3][ny-1-j]),
                        0.5*(urms[0][j]+urms[0][ny-1-j]),
                        0.5*(urms[1][j]+urms[1][ny-1-j]),
                        0.5*(urms[2][j]+urms[2][ny-1-j]),
                        0.5*(urms[3][j]-urms[3][ny-1-j]),
                        0.5*(urms[4][j]+urms[4][ny-1-j]));

   fflush(umean_out);
   fclose(umean_out);
   free(buf);
  }

	free(X.x);
	free(X.y);

  free_dmatrix(umean_tmp,0,0);
  free(xz);
  free(xz_tmp);

}



void append_white_noise(Domain *omega){
  register int i,j,k;
	Element_List *U,*V,*W;
  Element      *E,*F,*G;
  int          nz = option("NZ"), nztot = option("NZTOT");

  double noise_amp =  dparam("NOISE_AMPLITUDE")? dparam("NOISE_AMPLITUDE"):0.001;  

  srand(time(NULL));

  U = omega->U;
  V = omega->V;
  W = omega->W;

  int htot = nz*U->htot; 

  for(i=0; i<htot; ++i)
    {
      U->base_h[i] += U->base_h[i]*noise_amp*((double)rand()/RAND_MAX-0.5);
      V->base_h[i] += V->base_h[i]*noise_amp*((double)rand()/RAND_MAX-0.5);
      W->base_h[i] += W->base_h[i]*noise_amp*((double)rand()/RAND_MAX-0.5);
    }

 }



void point_output_prepare(Domain *omega){
	  Element_List *U;
    Element      *E;
    double       *zz = (double*)0;
    int          eid;
    int          sum_tmp;
    char         status;
    double       PI = 4.*atan(1.);
    int          nz = option("NZ"), nztot = option("NZTOT");

    double Lx = dparam("BOX_LX");
    double Ly = dparam("BOX_LY");

    int    Nx =  (int) dparam("BOX_NX");
    int    Ny =  (int) dparam("BOX_NY");

    int    Nc =  (int) dparam("CIR_NC");
    int    Nr =  (int) dparam("CIR_NR");

    double x0 = dparam("X0")? dparam("X0"):0.5;
    double y0 = dparam("Y0");

    double radius = dparam("RADIUS")? dparam("RADIUS"):0.5;
    int    procid  = option("PROCID");

//    if(x0<0.5)
//     {
//      fprintf(stderr,"wrong value of X0, it should be large than 0.5 ! \n");
//      exit(1);
//     }

    int num_points      = 0;
    int num_x_points    = 0;
    int num_y_points    = 0;
    int num_hist_points = 0;

    int num_c_points    = 0;
    int num_r_points    = 0;
    int num_h_points    = 0;
    
    Coord X;

    if(dparam("CYLINDER"))
    {
      num_x_points    = Nx;
      num_y_points    = 3*Ny;
      num_hist_points = 7;
//    int num_points_c    = Nc;

      double len_x = Lx/(double)num_x_points;
      double len_y = Ly/(double)Ny;
    
      num_points   = num_x_points+num_y_points + num_hist_points;// +num_points_c;
	
	     X.x  = dvector(0, num_points-1);
	     X.y  = dvector(0, num_points-1);


       int point_offset = 0;

      for(int i=0; i<num_x_points; ++i)
      {
        int k = i;

        X.x[k] = x0+i*len_x;
        X.y[k] = 0.;
      }

        point_offset += num_x_points;

    for(int i=0; i<Ny; i++)
      {
        int k = point_offset+i;

        X.x[k] = 1.0;
        X.y[k] = -Ly/2.+len_y/2.+i*len_y;

        X.x[Ny+k] = 2.0;
        X.y[Ny+k] = -Ly/2.+len_y/2.+i*len_y;

        X.x[2*Ny+k] = 3.0;
        X.y[2*Ny+k] = -Ly/2.+len_y/2.+i*len_y;
      }
        
     if(procid==0)
      {
    #ifndef THERMO
       mean_new = dmatrix(0, num_x_points-1, 0, 3);
       dzero(num_x_points*4, mean_new[0], 1);

       mean_new1 = dmatrix(0, 3*Ny-1, 0, 8);
       dzero(3*Ny*9, mean_new1[0], 1);
    #else
       mean_new = dmatrix(0, num_x_points-1, 0, 4);
       dzero(num_x_points*5, mean_new[0], 1);

       mean_new1 = dmatrix(0, 3*Ny-1, 0, 11);
       dzero(3*Ny*12, mean_new1[0], 1);
    #endif
      }

        point_offset += num_y_points;

      {
        int k = point_offset;

        X.x[k] = 3.;
        X.y[k] = 0.;

        X.x[k+1] = 5.;
        X.y[k+1] = 0.;

        X.x[2+k] = 7.;
        X.y[2+k] = 0.;

        X.x[3+k] = 0.2;
        X.y[3+k] = 0.55;
        
        X.x[4+k] = 0.42;
        X.y[4+k] = 0.55;

        X.x[5+k] = 0.54;
        X.y[5+k] = 0.65;

        X.x[6+k] = 2.0;
        X.y[6+k] = 0.5;

       }
//    for(int i=0; i<Nc; ++i)
//      {
//        int k = point_offset+j;
//        double alpha = i*lc;
//
//        X.x[k] = radius*cos(PI-alpha);
//        X.y[k] = radius*sin(PI-alpha);
//      }
        
        point_offset += num_hist_points;
     }

   if(dparam("PIPE"))
    {

     double len_c = 2.*PI/(double)Nc;
     double len_r = radius/(double)Nr;

     double gamma = 2.25;
     
     int point_offset = 0;

     num_c_points    = Nc;
     num_r_points    = Nr;
     num_h_points    = 3;
     num_hist_points = num_h_points*Nc;//center line r=0.25 and r=0.495 3 positions

     num_points    = Nc*Nr+num_hist_points;
//     num_points    = Nc*Nr;

	   X.x  = dvector(0, num_points-1);
	   X.y  = dvector(0, num_points-1);
     
     omega->interp_r = dvector(0, Nr-1);
     
     int num_out = Nr*nztot;

     if(procid==0)
      {
       int nf = 16;
       if(iparam("IPIPE_RMS"))
         nf=nf+3*2;

       mean_new = dmatrix(0, num_out-1, 0, nf-1);
       dzero(num_out*nf, mean_new[0], 1);

       mean_new1 = dmatrix(0, Nr-1, 0, nf-1);
       dzero(Nr*nf, mean_new1[0], 1);
      }

     for(int i=0; i<Nc; ++i)
      for(int j=0; j<Nr; ++j)
        {
          int m = i*(Nr+num_h_points)+j;
          double theta = 0.5*len_c + i*len_c;
          double rr =  -radius*tanh(gamma*((double)(Nr-j)/(double)Nr-1.))/tanh(gamma);
//          double rr    = 0.5*len_r + j*len_r;
          omega->interp_r[j] = rr;
         
          X.x[m] = rr*cos(theta);
          X.y[m] = rr*sin(theta);

        }
        
//     point_offset += Nc*Nr;

     //hist points
      {

        double rr1[3] = {0,0.25,0.495};

       for(int i=0; i<Nc; ++i)
        for(int j=Nr,k=0; j<Nr+3; ++j,++k)
        {
          int m = i*(Nr+num_h_points)+j;
           double theta = 0.5*len_c + i*len_c;
           X.x[m] = rr1[k]*cos(theta);
           X.y[m] = rr1[k]*sin(theta);

        }

      }


    }

    Coord XX,*YY;
    XX.x = dvector(0,2-1); 
    XX.y = XX.x + 1;
//    YY.x = XX.x + 1;
//    YY.y = XX.x + 1;

 //   dfill(2,1e8,XX.x,1);

    int data_offset = 3;
    int data_size = data_offset*num_points;

    double * points_locations;
   // if(!points_locations)
     points_locations = dvector(0,data_size-1);
    
     dfill(data_size, 1e8, points_locations, 1);

    U = omega->U;
    int num = 0;
    for(int i=0; i<data_size; i=i+data_offset)
     {
         XX.x[0] = X.x[num]; 
         XX.y[0] = X.y[num]; 

         int *eid ;

//         Prt_find_local_coords(U->flevels[0],XX,&eid,&YY,status);
         Find_local_coords(U,&XX,1,&eid,&YY);
         
         points_locations[i] = (double) eid[0];
         points_locations[i+1] = YY->x[0];
         points_locations[i+2] = YY->y[0];
         
//         fprintf(stdout,"num = %d loc = %d x = %lf y = %lf \n",num, (int) points_locations[i], points_locations[i+1], points_locations[i+2]);
         num++;
    }

   gsync ();

   if(num != num_points )
      {
       fprintf(stderr," wrong statistics points ! total_points = %d num = %d \n",num_points,num);
       exit(-1);
     }


   omega->points_locations = points_locations;

   omega->num_x_points     = num_x_points;
   omega->num_y_points     = num_y_points;
   omega->num_hist_points  = num_hist_points;
   omega->num_c_points     = num_c_points;
   omega->num_r_points     = num_r_points;
   omega->num_points       = num_points;

   omega->data_size     = data_size;
   omega->data_offset   = data_offset;

//  ROOTONLY
//       for(int i=0,j=0; i<send_size; i=i+data_offset,++j)
//     fprintf(stdout," p = %d  index = %d  proc = %d   eid = %d  loc_x = %f  loc_y = %f  loc_z = %f \n",
//        j,(int)omega->points_location[i+0],(int)omega->points_location[i+1], (int)omega->points_location[i+2],omega->points_location[i+3],omega->points_location[i+4],omega->points_location[i+5]);

     free(X.x);
     free(X.y);

     free(XX.x);
}


void output_point_value(Domain *omega, int step, double time){
 Element  *E;
 Element_List  **V;
 int      eid, qa, qb, qc;
 double   *ha,*hb,*hc,*za,*wa,*zb,*wb,*zc,*wc;
 int      myid = pllinfo.procid;
 int      int_tmp;
 double   double_tmp;
 int      iostep = iparam("IOSTEP");
 int      hisstep   = option("hisstep");
 if(!hisstep) hisstep = iparam("HISSTEP");

 int      stastep   = iparam("STASTEP");
 if(!stastep)
  {
   fprintf(stderr,"forget to set STASTEP ! \n");
   exit(1);
  }

 int      nz = option("NZ"), nztot = option("NZTOT");
 int      procid  = option("PROCID");
 int      nprocs  = option("NPROCS");
 int      htot = omega->U->nz*omega->U->htot;

 double   Lx = dparam("BOX_LX");
 double   Ly = dparam("BOX_LY");
 int      Nx = (int)dparam("BOX_NX");
 int      Ny = (int)dparam("BOX_NY");
 int      Nc =  (int) dparam("CIR_NC");
 int      Nr =  (int) dparam("CIR_NR");
 double   radius = dparam("RADIUS");
 double   PI = 4.*atan(1.);
 double   x0 = dparam("X0")? dparam("X0"):0.5;

 int      num_points        = omega->num_points;
 int      num_x_points      = omega->num_x_points;
 int      num_y_points      = omega->num_y_points;
 int      num_c_points      = omega->num_c_points;
 int      num_r_points      = omega->num_r_points;
 int      num_hist_points   = omega->num_hist_points;
 int      data_size         = omega->data_size;
 int      data_offset       = omega->data_offset;
 double   *points_locations = omega->points_locations;
 double   *interp_r         = omega->interp_r;

 int      num_y_interp       = 3;

 double   len_x = Lx/(double)Nx;
 double   len_y = Ly/(double)Ny;
 double   len_c = 2.*PI/(double)Nc;
 double   len_r = radius/(double)Nr;

 static   int statistics_step = 1;

 double   fac = 1. /(double) statistics_step;

 #ifndef THERMO
 int      nfields = omega->U->fhead->dim()+2;
 #else
 int      nfields = omega->U->fhead->dim()+3;
 #endif
 int      nfields_plus = nfields;
 #ifndef THERMO
 if(nfields != 4)
 #else
 if(nfields != 5)
 #endif
   {
     fprintf(stderr," only works for 3D case \n");
     exit(-1);
   }
 

#ifdef OUTPUT_VISCOSITY
   nfields_plus += 2;
#endif

 if(iparam("IPIPE_RMS"))
   nfields_plus = nfields+3;

  V = (Element_List**) malloc((nfields_plus)*sizeof(Element_List*));
  
//now Ut,Vt,Wt have the vibration, slow part and fast part pressure rms, respectively. 
 if(iparam("IPIPE_RMS"))
  {
//#ifdef MAP
//    compute_rms_pressure(omega, omega->mapx, omega->mapy, step);
//#endif 
    compute_cylindrical_vorticity(omega,step);

    dcopy(htot, omega->Ut->base_h, 1, omega->Uf->base_h, 1);
    dcopy(htot, omega->Vt->base_h, 1, omega->Vf->base_h, 1);
    dcopy(htot, omega->Wt->base_h, 1, omega->Wf->base_h, 1);
  }
//  dcopy(htot, omega->P->base_h, 1, omega->Pt->base_h, 1);

//  dcopy(htot,omega->uk[1],1,omega->Ut->base_h,1);
//  dcopy(htot,omega->vk[1],1,omega->Vt->base_h,1);
//  dcopy(htot,omega->wk[1],1,omega->Wt->base_h,1);
//  dcopy(htot,omega->pk[1],1,omega->Pt->base_h,1);
  
   omega->Ut->Trans(omega->U->htot, omega->U->base_h, omega->Ut->base_h, F_to_P); // in (Q,P)-space
   omega->Vt->Trans(omega->V->htot, omega->V->base_h, omega->Vt->base_h, F_to_P); // in (Q,P)-space
   omega->Vt->Trans(omega->W->htot, omega->W->base_h, omega->Wt->base_h, F_to_P); // in (Q,P)-space
   omega->Pt->Trans(omega->P->htot, omega->P->base_h, omega->Pt->base_h, F_to_P); // in (Q,P)-space

 #ifdef THERMO
   omega->T->Trans(omega->T->htot, omega->T->base_h, omega->Tf->base_h, F_to_P); // in (Q,P)-space
 #endif
 
#ifdef OUTPUT_VISCOSITY
   omega->visc1->Set_state('p');
   omega->visc2->Set_state('p');
   omega->visc1->Trans(omega->visc1->htot, omega->visc1->base_h, omega->visc1->base_h, F_to_P); // in (Q,P)-space
   omega->visc2->Trans(omega->visc2->htot, omega->visc2->base_h, omega->visc2->base_h, F_to_P); // in (Q,P)-space
#endif

  V[0]   = omega->Ut;
  V[1]   = omega->Vt;
  V[2]   = omega->Wt;
  V[3]   = omega->Pt;

#if defined(THERMO) && !defined(OUTPUT_VISCOSITY)
  V[4] = omega->Tf;
#endif

#if defined(OUTPUT_VISCOSITY) && !defined(THERMO)
  V[4] = omega->visc1;
  V[5] = omega->visc2;
#endif

#if defined(OUTPUT_VISCOSITY) && defined(THERMO)
  V[4] = omega->Tf;
  V[5] = omega->visc1;
  V[6] = omega->visc2;
#endif


 if(iparam("IPIPE_RMS"))
 {
  V[6]   = omega->Uf;
  V[7]   = omega->Vf;
  V[8]   = omega->Wf;

//  V[7]   = omega->U;
//  V[8]   = omega->V;
//  V[9]   = omega->W;
//  V[10]  = omega->P;
 }
     

  Coord X;
  X.x = dvector(0,2-1); 
  X.y = X.x + 1;
  
  int send_length = num_points*nz*nfields_plus;
  double *data = dvector(0, send_length-1); 
  
#ifdef MAP
  rfftw_plan MPlan = rplan, 
             MPlan_inv = rplan_inv;
  int NZ = nz;
  int NZTOT =nztot;
   int nZtot = NZTOT; // keep this value as it may change because of dealiasing
  if (option("dealias")) {
      NZ    = 3*NZ/2;
      NZTOT = 3*NZTOT/2;
      MPlan     = rplan32;    
      MPlan_inv = rplan_inv32;

         // zero out for dealiasing
      dzero (NZTOT-nZtot, omega->mapx->z+nZtot, 1);
      dzero (NZTOT-nZtot, omega->mapy->z+nZtot, 1);
#ifdef TMAP
      dzero (NZTOT-nZtot, omega->mapx->t+nZtot, 1);
       dzero (NZTOT-nZtot, omega->mapy->t+nZtot, 1);
#endif
    }

        // invFFT the maps to physical space
   rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) omega->mapx->z, 1, 0, 0, 0, 0);
   rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) omega->mapy->z, 1, 0, 0, 0, 0);
#ifdef TMAP
   rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) omega->mapx->t, 1, 0, 0, 0, 0);
   rfftw(MPlan_inv, 1, (FFTW_COMPLEX *) omega->mapy->t, 1, 0, 0, 0, 0);
#endif

#endif
     
  E = V[0]->fhead;
	qa = E->qa;
	qb = E->qb;
	ha = dvector(0,qa-1);
	hb = dvector(0,qb-1);

  for(int k=0;k < nz; k++)
	  for(int n = 0; n < nfields_plus; ++n)
      for(int i=0, j=0;i<data_size; i=i+data_offset, j++)
       {
          int m = k*nfields_plus*num_points+n*num_points+j;

          eid      =  (int) (points_locations[i]);

          X.x[0] = points_locations[i+1]; 
          X.y[0] = points_locations[i+2]; 

		      E = V[n]->flevels[k]->flist[eid];
		 //     E = V[0]->flevels[k]->flist[eid];
     // 
     //     for (int ii=0; ii<qa*qb; ++ii)
     //     E->h[0][ii] = V[n]->flevels[k]->flist[eid]->h[0][ii]; //copy into E

	        E->GetZW(&za, &wa, &zb, &wb, &zc, &wc);

          get_point_shape_2d(E,X.x[0],X.y[0],ha,hb);
		      data[m] = eval_field_at_pt_2d(qa,qb, E->h[0], ha, hb);
        }
	
   gsync ();

   int comm_len = send_length*nprocs;  

   double *recv_buff = dvector(0, comm_len-1);
   
   MPI_Allgather(&data[0], send_length, MPI_DOUBLE, 
				 &recv_buff[0], send_length, MPI_DOUBLE, MPI_COMM_WORLD);
  
//  if( (step % (hisstep) == 0 ) && (procid == 0) )
   if( (procid == 0) )
   {
//      double max_pres = -1e6;
//      double min_pres = 1e6;

     if(dparam("CYLINDER"))
     {
      
      double **mean_x = dmatrix(0, num_x_points-1, 0, nfields-1);
      dzero(num_x_points*(nfields), mean_x[0], 1);
       
      for(int p=0; p<nprocs; ++p)
        for(int k=0;k < nz; k++)
//	       for(int n = 0; n < nfields; ++n)
           for(int j=0; j<num_x_points; ++j)
              {
//                 int m = p*nfields*num_points*nz 
//                   + k*nfields*num_points
//                   + n*num_points
//                   + j;
                 
                 int mu = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 0*num_points
                   + j;
                 
                 int mv = p*nfields*num_points*nz 
                   + k*nfields*num_points
                    + 1*num_points
                   + j;

                 int mw = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 2*num_points
                   + j;

                 int mp = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 3*num_points
                   + j;
      #ifdef THERMO
                 int mt = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 4*num_points
                   + j;
      #endif

                 double uu = recv_buff[mu];
                 double vv = recv_buff[mv];
                 double ww = recv_buff[mw];
                 double pp = recv_buff[mp];
      #ifdef THERMO
                 double tt = recv_buff[mt];
      #endif

#ifdef MAP
#ifdef TMAP // time dependent mesh
                 uu +=  (omega->mapx->t[k] + omega->mapx->z[k] * ww);
                 vv +=  (omega->mapy->t[k] + omega->mapy->z[k] * ww);
#else // only Z-deformation
                 uu  +=  omega->mapx->z[k] * ww;
                 vv  +=  omega->mapy->z[k] * ww;
#endif
#endif
//                 mean_x[j][n] += recv_buff[m]; 
                 mean_x[j][0] += uu; 
                 mean_x[j][1] += vv; 
                 mean_x[j][2] += ww; 
                 mean_x[j][3] += pp; 
      #ifdef THERMO
                 mean_x[j][4] += tt; 
      #endif

              }

     for(int j=0; j<num_x_points; ++j)
       for(int n=0; n<nfields; ++n)
        mean_x[j][n] /= (double)(nprocs*nz);

	        for(int n = 0; n < nfields; ++n)
           for(int i=0; i<num_x_points; ++i)
              mean_new[i][n] += mean_x[i][n];     //mean_x
       
       int offset = num_x_points;

    #ifndef THERMO 
      double **mean_y = dmatrix(0, num_y_points-1, 0, nfields+5-1);
      dzero(num_y_points*(nfields+5), mean_y[0], 1);
    #else
      double **mean_y = dmatrix(0, num_y_points-1, 0, nfields+7-1);
      dzero(num_y_points*(nfields+7), mean_y[0], 1);
    #endif

      for(int p=0; p<nprocs; ++p)
        for(int k=0;k < nz; k++)
//	       for(int n = 0; n < nfields; ++n)
           for(int j=offset,i=0; i<num_y_points; ++j,++i)
              {
                 int mu = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 0*num_points
                   + j;
                 
                 int mv = p*nfields*num_points*nz 
                   + k*nfields*num_points
                    + 1*num_points
                   + j;

                 int mw = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 2*num_points
                   + j;

                 int mp = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 3*num_points
                   + j;
      #ifdef THERMO
                 int mt = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 4*num_points
                   + j;
      #endif

                 double uu = recv_buff[mu];
                 double vv = recv_buff[mv];
                 double ww = recv_buff[mw];
                 double pp = recv_buff[mp];
      #ifdef THERMO
                 double tt = recv_buff[mt];
      #endif
#ifdef MAP
#ifdef TMAP // time dependent mesh
                 uu +=  (omega->mapx->t[k] + omega->mapx->z[k] * ww);
                 vv +=  (omega->mapy->t[k] + omega->mapy->z[k] * ww);
#else // only Z-deformation
                 uu  +=  omega->mapx->z[k] * ww;
                 vv  +=  omega->mapy->z[k] * ww;
#endif
#endif
                 mean_y[i][0] += uu;  // U
                 mean_y[i][1] += vv;  // V
                 mean_y[i][2] += ww;  // W
                 mean_y[i][3] += pp;  // P

                 mean_y[i][4] += uu*uu;
                 mean_y[i][5] += vv*vv;
                 mean_y[i][6] += ww*ww;
                 mean_y[i][7] += uu*vv;
                 mean_y[i][8] += pp*pp;
      #ifdef THERMO
                 mean_y[i][9] += tt;
                 mean_y[i][10] += tt*tt;
                 mean_y[i][11] += uu*tt;
      #endif
              }

      for(int i=0; i<num_y_points; ++i)
      #ifndef THERMO
        for(int n=0; n<nfields+5; ++n)
      #else
        for(int n=0; n<nfields+7; ++n)
      #endif
            mean_y[i][n] /= (double) (nprocs*nz);
      
      for(int i=0; i<num_y_points; ++i)
        {
          mean_y[i][4]   -= mean_y[i][0]*mean_y[i][0];
          mean_y[i][5]   -= mean_y[i][1]*mean_y[i][1];
          mean_y[i][6]   -= mean_y[i][2]*mean_y[i][2];
          mean_y[i][7]   -= mean_y[i][0]*mean_y[i][1];
          mean_y[i][8]   -= mean_y[i][3]*mean_y[i][3];
      #ifdef THERMO
          mean_y[i][10]   -= mean_y[i][9]*mean_y[i][9];
          mean_y[i][11]   -= mean_y[i][0]*mean_y[i][9];
      #endif
         }


      #ifndef THERMO
	        for(int n = 0; n < nfields+5; ++n)
      #else
	        for(int n = 0; n < nfields+7; ++n)
      #endif
           for(int i=0; i<num_y_points; ++i)
              mean_new1[i][n] += mean_y[i][n];     //mean_y

      if(statistics_step == stastep)
        {
          FILE *file_x;
	        char *buf_x = (char*) calloc(BUFSIZ,sizeof(char));
          FILE *file_y;
	        char *buf_y = (char*) calloc(BUFSIZ,sizeof(char));
	        sprintf(buf_x, "mean_x.txt");
	        sprintf(buf_y, "mean_y.txt");
 
          file_x = fopen(buf_x,"aw");
          file_y = fopen(buf_y,"aw");

         fac = 1./(double) statistics_step;
         // double gamma = 2.75;
	        for(int n = 0; n < nfields; ++n)
           for(int i=0; i<num_x_points; ++i)
              mean_new[i][n] *= fac;     //mean_x

      #ifndef THERMO
	         for(int n = 0; n < nfields+5; ++n)
      #else
	         for(int n = 0; n < nfields+7; ++n)
      #endif
           for(int i=0; i<num_y_points; ++i)
              mean_new1[i][n] *= fac;     //mean_y

           for(int j=0; j<num_x_points; ++j)
           {
              double xx = x0+len_x*j;
               fprintf(file_x," %lf  %lf",time, xx);
              for(int n=0; n<nfields; ++n)
                 fprintf(file_x," %lf  ", mean_new[j][n]);
                 fprintf(file_x,"\n");
            }

          for(int q=0; q<num_y_interp; ++q) 
            for(int i=0; i<num_y_points/num_y_interp; ++i)
           {
             int m = q*num_y_points/num_y_interp + i;

             double yy = -Ly/2.+len_y/2.+len_y*i;
             fprintf(file_y," %d  %lf ",q, yy);
                for(int n=0; n<nfields+5; ++n)
                 fprintf(file_y," %lf  ", mean_new1[m][n]);

              fprintf(file_y,"\n");
           }

            statistics_step = 1;
            fflush(file_x);
            fclose(file_x);
       
            fflush(file_y);
            fclose(file_y);
            free(buf_x);
            free(buf_y);
        }

          free_dmatrix(mean_x,0,0);
          free_dmatrix(mean_y,0,0);
      }
   
      if(dparam("PIPE"))
       {
#ifdef THERMO
         fprintf(stderr,"PointOutput not be debugged ...\n");
         exit(-1);
#endif
         int num_mean = nfields+12;

         if(iparam("IPIPE_RMS")) num_mean = num_mean+3*2;

//         double **mean_tmp = dmatrix(0, num_r_points-1, 0, nfields+5-1);
//         dzero(num_r_points*(nfields+5), mean_tmp[0], 1);
         int num_out = num_r_points*nz*nprocs;
         double **mean_tmp  = dmatrix(0, num_out-1, 0, num_mean-1);
         double **mean_tmp1 = dmatrix(0, num_r_points-1, 0, num_mean-1);
         dzero(num_out*(num_mean), mean_tmp[0], 1);
         dzero(num_r_points*(num_mean), mean_tmp1[0], 1);

         int pv = 0, ps = 0, pf =0,
             pv0 = 0, pv1 = 0, pv2 =0 ,pv3 = 0;
         for(int p=0; p<nprocs; ++p)
           for(int k=0;k < nz; k++)
//	         for(int n = 0; n < nfields; ++n)
             for(int j=0; j<num_c_points; ++j)
              for(int i=0; i<num_r_points; ++i)
               {
                 double phi = j*len_c;

                 int m = p*nz*num_r_points+k*num_r_points+i;

                 int mx = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 0*num_points
                    + j*(num_r_points+3) 
                   + i; 
                 int my = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 1*num_points
                   + j*(num_r_points+3)
                   + i;
                 int mz = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 2*num_points
                   + j*(num_r_points+3)
                   + i; 
                 int mp = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 3*num_points
                   + j*(num_r_points+3)
                   + i; 

            if(iparam("IPIPE_RMS")){
                   pv = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 4*num_points
                   + j*(num_r_points+3)
                   + i; 

                   ps = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 5*num_points
                   + j*(num_r_points+3)
                   + i; 

                   pf = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 6*num_points
                   + j*(num_r_points+3)
                   + i; 
/*
                   pv0 = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 7*num_points
                   + j*(num_r_points+3)
                   + i; 

                   pv1 = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 8*num_points
                   + j*(num_r_points+3)
                   + i; 

                   pv2 = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 9*num_points
                   + j*(num_r_points+3)
                   + i; 

                   pv3 = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 10*num_points
                   + j*(num_r_points+3)
                   + i;
*/
           }
 
                 double uu = recv_buff[mx];
                 double vv = recv_buff[my];
                 double ww = recv_buff[mz];
                 double pp = recv_buff[mp];

//                 max_pres = max(pp,max_pres);
//                 min_pres = min(pp,min_pres);

                 double uu1 = uu;
                 double vv1 = vv;
                 double ww1 = ww;
#ifdef MAP
#ifdef TMAP // time dependent mesh
                 uu1 +=  (omega->mapx->t[k] + omega->mapx->z[k] * ww1);
                 vv1 +=  (omega->mapy->t[k] + omega->mapy->z[k] * ww1);
#else // only Z-deformation
                 uu1  +=  omega->mapx->z[k] * ww1;
                 vv1  +=  omega->mapy->z[k] * ww1;
#endif
#endif
                 double ur   = cos(phi)*uu+sin(phi)*vv;
                 double uphi = -sin(phi)*uu+cos(phi)*vv;
                 double uz   = ww;

                 double ur1   = cos(phi)*uu1+sin(phi)*vv1;
                 double uphi1 = -sin(phi)*uu1+cos(phi)*vv1;
                 double uz1   = ww1;
                  
                  mean_tmp1[i][0] += ur1;     //U_r
                  mean_tmp1[i][1] += uphi1;   //U_phi
                  mean_tmp1[i][2] += uz1;    //U_z
                  mean_tmp1[i][3] += pp;   //P
                 
                  mean_tmp1[i][4] += ur1*ur1; //u_r*u_r
                  mean_tmp1[i][5] += uphi1*uphi1;
                  mean_tmp1[i][6] += uz1*uz1;
                  mean_tmp1[i][7] += uz1*ur1;
                  mean_tmp1[i][8] += pp*pp;

                  mean_tmp1[i][9]  += ur;     //U_r
                  mean_tmp1[i][10] += uphi;   //U_phi
                  mean_tmp1[i][11] += uz;    //U_z
                  mean_tmp1[i][12] += ur*ur; //u_r*u_r
                  mean_tmp1[i][13] += uphi*uphi;
                  mean_tmp1[i][14] += uz*uz;
                  mean_tmp1[i][15] += uz*ur;

                if(iparam("IPIPE_RMS")){
                  mean_tmp1[i][16] += recv_buff[pv];     // omega_r
                  mean_tmp1[i][17] += recv_buff[ps];   // omega_phi
                  mean_tmp1[i][18] += recv_buff[pf];     // omega_z
                  mean_tmp1[i][19] += recv_buff[pv]*recv_buff[pv];     // rms
                  mean_tmp1[i][20] += recv_buff[ps]*recv_buff[ps];     // rms
                  mean_tmp1[i][21] += recv_buff[pf]*recv_buff[pf];     // rms
//                  mean_tmp1[i][19] += recv_buff[pv0];  // rms omega_r 
//                  mean_tmp1[i][20] += recv_buff[pv1]; // rms omega_phi
//                  mean_tmp1[i][21] += recv_buff[pv2];   //rms omega_z
//                  mean_tmp1[i][22] += recv_buff[pv3];   //prms 3

//                  mean_tmp1[i][23] += recv_buff[pv]*recv_buff[pv];     // vibration presssure rms
//                 mean_tmp1[i][24] += recv_buff[ps]*recv_buff[ps];     // slow presussre rms
//                  mean_tmp1[i][25] += recv_buff[pf]*recv_buff[pf];     // fast pressure rms
//                  mean_tmp1[i][26] += recv_buff[pv0]*recv_buff[pv0];   // prms 0 
//                  mean_tmp1[i][27] += recv_buff[pv1]*recv_buff[pv1];   //prms 1
//                  mean_tmp1[i][28] += recv_buff[pv2]*recv_buff[pv2];   //prms 2
//                  mean_tmp1[i][29] += recv_buff[pv3]*recv_buff[pv3];   //prms 3

                  }

                  mean_tmp[m][0] += ur1;     //U_r
                  mean_tmp[m][1] += uphi1;   //U_phi
                  mean_tmp[m][2] += uz1;    //U_z
                  mean_tmp[m][3] += pp;   //P
                 
                  mean_tmp[m][4] += ur1*ur1; //u_r*u_r
                  mean_tmp[m][5] += uphi1*uphi1;
                  mean_tmp[m][6] += uz1*uz1;
                  mean_tmp[m][7] += uz1*ur1;
                  mean_tmp[m][8] += pp*pp;

                  mean_tmp[m][9]  += ur;     //U_r
                  mean_tmp[m][10] += uphi;   //U_phi
                  mean_tmp[m][11] += uz;    //U_z
                  mean_tmp[m][12] += ur*ur; //u_r*u_r
                  mean_tmp[m][13] += uphi*uphi;
                  mean_tmp[m][14] += uz*uz;
                  mean_tmp[m][15] += uz*ur;

                if(iparam("IPIPE_RMS")){
                  mean_tmp[m][16] += recv_buff[pv];     // mean vibration presssure rms
                  mean_tmp[m][17] += recv_buff[ps];   // mean slow presussre rms
                  mean_tmp[m][18] += recv_buff[pf];     // mean fast pressure rms
                  mean_tmp[m][19] += recv_buff[pv]*recv_buff[pv];     // vibration presssure rms
                  mean_tmp[m][20] += recv_buff[ps]*recv_buff[ps];     // slow presussre rms
                  mean_tmp[m][21] += recv_buff[pf]*recv_buff[pf];     // fast pressure rms
//                  mean_tmp[i][19] += recv_buff[pv0];  // prms 0 
//                  mean_tmp[i][20] += recv_buff[pv1]; //prms 1
//                  mean_tmp[i][21] += recv_buff[pv2];   //prms 2
//                  mean_tmp[i][22] += recv_buff[pv3];   //prms 3

//                  mean_tmp[i][23] += recv_buff[pv]*recv_buff[pv];     // vibration presssure rms
//                  mean_tmp[i][24] += recv_buff[ps]*recv_buff[ps];     // slow presussre rms
//                  mean_tmp[i][25] += recv_buff[pf]*recv_buff[pf];     // fast pressure rms
//                  mean_tmp[i][26] += recv_buff[pv0]*recv_buff[pv0];   // prms 0 
//                  mean_tmp[i][27] += recv_buff[pv1]*recv_buff[pv1];   //prms 1
//                  mean_tmp[i][28] += recv_buff[pv2]*recv_buff[pv2];   //prms 2
//                  mean_tmp[i][29] += recv_buff[pv3]*recv_buff[pv3];   //prms 3
                  }

               }

	       for(int n = 0; n < num_mean; ++n)
          for(int i=0; i<num_out; ++i)
            mean_tmp[i][n] /= (double) (num_c_points);
	       
          for(int i=0; i<num_out; ++i)
           {
            mean_tmp[i][4] -=mean_tmp[i][0]*mean_tmp[i][0];
            mean_tmp[i][5] -=mean_tmp[i][1]*mean_tmp[i][1];
            mean_tmp[i][6] -=mean_tmp[i][2]*mean_tmp[i][2];
            mean_tmp[i][7] -=mean_tmp[i][0]*mean_tmp[i][2];
            mean_tmp[i][8] -=mean_tmp[i][3]*mean_tmp[i][3];

            mean_tmp[i][12] -=mean_tmp[i][9]*mean_tmp[i][9];
            mean_tmp[i][13] -=mean_tmp[i][10]*mean_tmp[i][10];
            mean_tmp[i][14] -=mean_tmp[i][11]*mean_tmp[i][11];
            mean_tmp[i][15] -=mean_tmp[i][9]*mean_tmp[i][11];

           if(iparam("IPIPE_RMS")){
            mean_tmp[i][19] -=mean_tmp[i][16]*mean_tmp[i][16];
            mean_tmp[i][20] -=mean_tmp[i][17]*mean_tmp[i][17];
            mean_tmp[i][21] -=mean_tmp[i][18]*mean_tmp[i][18];
//            mean_tmp[i][26] -=mean_tmp[i][19]*mean_tmp[i][19];
//            mean_tmp[i][27] -=mean_tmp[i][20]*mean_tmp[i][20];
//            mean_tmp[i][28] -=mean_tmp[i][21]*mean_tmp[i][21];
//            mean_tmp[i][29] -=mean_tmp[i][22]*mean_tmp[i][22];
             }
           }

	        for(int n = 0; n < num_mean; ++n)
          for(int i=0; i<num_r_points; ++i)
            mean_tmp1[i][n] /= (double) (nprocs*num_c_points*nz);

          for(int i=0; i<num_r_points; ++i)
           {
            mean_tmp1[i][4] -=mean_tmp1[i][0]*mean_tmp1[i][0];
            mean_tmp1[i][5] -=mean_tmp1[i][1]*mean_tmp1[i][1];
            mean_tmp1[i][6] -=mean_tmp1[i][2]*mean_tmp1[i][2];
            mean_tmp1[i][7] -=mean_tmp1[i][0]*mean_tmp1[i][2];
            mean_tmp1[i][8] -=mean_tmp1[i][3]*mean_tmp1[i][3];

            mean_tmp1[i][12] -=mean_tmp1[i][9]*mean_tmp1[i][9];
            mean_tmp1[i][13] -=mean_tmp1[i][10]*mean_tmp1[i][10];
            mean_tmp1[i][14] -=mean_tmp1[i][11]*mean_tmp1[i][11];
            mean_tmp1[i][15] -=mean_tmp1[i][9]*mean_tmp1[i][11];

           if(iparam("IPIPE_RMS")){
            mean_tmp1[i][19] -=mean_tmp1[i][16]*mean_tmp1[i][16];
            mean_tmp1[i][20] -=mean_tmp1[i][17]*mean_tmp1[i][17];
            mean_tmp1[i][21] -=mean_tmp1[i][18]*mean_tmp1[i][18];
//            mean_tmp1[i][26] -=mean_tmp1[i][19]*mean_tmp1[i][19];
//            mean_tmp1[i][27] -=mean_tmp1[i][20]*mean_tmp1[i][20];
//            mean_tmp1[i][28] -=mean_tmp1[i][21]*mean_tmp1[i][21];
//            mean_tmp1[i][29] -=mean_tmp1[i][22]*mean_tmp1[i][22];
             }
           }

	        for(int n = 0; n < num_mean; ++n)
           for(int i=0; i<num_out; ++i)
              mean_new[i][n] += mean_tmp[i][n];     //U_r

	        for(int n = 0; n < num_mean; ++n)
           for(int i=0; i<num_r_points; ++i)
              mean_new1[i][n] += mean_tmp1[i][n];     //U_r

      if(statistics_step == stastep)
        {
         FILE *file_phase,*file_r,*file_new;
	       char *buf_phase = (char*) calloc(BUFSIZ,sizeof(char));
	       char *buf_r = (char*) calloc(BUFSIZ,sizeof(char));
	       char *buf_new = (char*) calloc(BUFSIZ,sizeof(char));
	       sprintf(buf_phase, "mean_phase.txt");
	       sprintf(buf_r, "mean_r.txt");
	       sprintf(buf_new, "mean_new.txt");
         file_phase = fopen(buf_phase,"aw");
         file_r = fopen(buf_r,"aw");
         file_new = fopen(buf_new,"aw");

//         fac = 1./(double) statistics_step;
         fac = 1./(double) stastep;
         // double gamma = 2.75;
	        for(int n = 0; n < num_mean; ++n)
           for(int i=0; i<num_out; ++i)
              mean_new[i][n] *= fac;     //U_r

	         for(int n = 0; n < num_mean; ++n)
           for(int i=0; i<num_r_points; ++i)
              mean_new1[i][n] *= fac;     //U_r

          for(int p=0; p<nprocs; ++p)
          for(int k=0;k < nz; k++)
           for(int i=num_r_points-1; i>=0; --i)
           {
              int m = p*nz*num_r_points+k*num_r_points+i;
              fprintf(file_phase," %lf  %lf ",time, radius-interp_r[i]);
                 for(int n=0; n<num_mean; ++n)
                   fprintf(file_phase," %lf  ", mean_tmp[m][n]);

              fprintf(file_phase," %d \n", i);

              fprintf(file_new," %lf  %lf ",time, radius-interp_r[i]);
                 for(int n=0; n<num_mean; ++n)
                   fprintf(file_new," %lf  ", mean_new[m][n]);

              fprintf(file_new," %d \n", i);
            }

           for(int i=num_r_points-1; i>=0; --i)
           {
              fprintf(file_r," %lf  %lf ",time, radius-interp_r[i]);
                 for(int n=0; n<num_mean; ++n)
                   fprintf(file_r," %lf  ", mean_new1[i][n]);
               fprintf(file_r," %d \n", i);
            }

          fflush(file_phase);
          fclose(file_phase);
          free(buf_phase);

          fflush(file_new);
          fclose(file_new);
          free(buf_new);
          
          dzero(num_out*num_mean, mean_new[0], 1);
          dzero(num_r_points*num_mean, mean_new1[0], 1);
          statistics_step = 0;

          fflush(file_r);
          fclose(file_r);
          free(buf_r);
       }
          free_dmatrix(mean_tmp,0,0);
          free_dmatrix(mean_tmp1,0,0);
          
//          fprintf(stdout,"Max Pressure = %g  Min pressure = %g \n",max_pres,min_pres);
     }


//     if( num_hist_points > 0 )
//      if(step == hisstep)
    if( (num_hist_points > 0) && (procid == 0) )
        {
         int offset = 0;

          if(dparam("CYLINDER"))
           offset =  num_x_points+num_y_points;
          if(dparam("PIPE"))
           offset =  num_c_points*num_r_points;
        
          int num_hist_data = nz*nprocs*num_hist_points;
          int nhp = nfields+3;
#ifdef  THERMO
              nhp ++;
#endif

           if(iparam("IPIPE_RMS"))
             nhp +=3;
#ifdef OUTPUT_VISCOSITY
             nhp +=2;
#endif

          double **hist_tmp = dmatrix(0, num_hist_data-1, 0, nhp-1);
          dzero(num_hist_data*(nhp), hist_tmp[0], 1);

         FILE * file_point;
	       char *buf_point = (char*) calloc(BUFSIZ,sizeof(char));
	       sprintf(buf_point, "hist.txt");
         file_point = fopen(buf_point,"aw");

        if(dparam("CYLINDER"))
         {
          for(int p=0; p<nprocs; ++p)
          for(int k=0;k < nz; k++)
	         for(int n = 0; n < nfields; ++n)
            for(int j=offset,i=0; j<num_points; ++j,++i)
            {
//              int m = p*nfields*num_points*nz 
//                   + k*nfields*num_points
//                   + n*num_points
//                   + j;
//

                 int mx = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 0*num_points
                    + j;

                 int my = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 1*num_points
                   + j;

                 int mz = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 2*num_points
                   + j;

                 int mp = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 3*num_points
                   + j;
      #ifdef THERMO
                 int mt = p*nfields*num_points*nz 
                   + k*nfields*num_points
                   + 4*num_points
                   + j;
      #endif
 
                 double uu = recv_buff[mx];
                 double vv = recv_buff[my];
                 double ww = recv_buff[mz];
                 double pp = recv_buff[mp];
      #ifdef THERMO
                 double tt = recv_buff[mt];
      #endif
#ifdef MAP
#ifdef TMAP // time dependent mesh
                 uu +=  (omega->mapx->t[k] + omega->mapx->z[k] * ww);
                 vv +=  (omega->mapy->t[k] + omega->mapy->z[k] * ww);
#else // only Z-deformation
                 uu  +=  omega->mapx->z[k] * ww;
                 vv  +=  omega->mapy->z[k] * ww;
#endif
#endif
                 int m = p*nz*num_hist_points+k*num_hist_points+i;

                 hist_tmp[m][0] = uu;
                 hist_tmp[m][1] = vv;
                 hist_tmp[m][2] = ww;
                 hist_tmp[m][3] = pp;
      #ifdef THERMO
                 hist_tmp[m][4] = tt;
      #endif
            }
             

          }

        if(dparam("PIPE"))
          {
//           int offset = num_c_points*num_r_points; 
//          int hists[3] = {4,48,63}; //approx r = 0.005,0.25,0.5
         int pv = 0, ps = 0, pf =0,
             pv0 = 0, pv1 = 0, pv2 =0 ,pv3 = 0;
#ifdef OUTPUT_VISCOSITY
         int ps1 = 0, ps2 = 0;
#endif
         for(int p=0; p<nprocs; ++p)
           for(int k=0;k < nz; k++)
//	         for(int n = 0; n < nfields; ++n)
             for(int j=0; j<num_c_points; ++j)
              for(int i=num_r_points,ii=0; i<num_r_points+3; ++i,++ii)
               {
                 int m = p*nz*num_c_points*3+k*num_c_points*3+j*3+ii;

                 double phi = j*len_c;

                 int mx = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 0*num_points
                    + j*(num_r_points+3) 
                   + i; 
                 int my = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 1*num_points
                   + j*(num_r_points+3)
                   + i;
                 int mz = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 2*num_points
                   + j*(num_r_points+3) 
                   + i; 
                 int mp = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 3*num_points
                   + j*(num_r_points+3)
                   + i; 

#ifdef OUTPUT_VISCOSITY
                   ps1 = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 4*num_points
                   + j*(num_r_points+3)
                   + i; 
                   ps2 = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 5*num_points
                   + j*(num_r_points+3)
                   + i; 
                
#else

            if(iparam("IPIPE_RMS")){
                   pv = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 4*num_points
                   + j*(num_r_points+3)
                   + i; 

                   ps = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 5*num_points
                   + j*(num_r_points+3)
                   + i; 

                   pf = p*nfields_plus*num_points*nz 
                   + k*nfields_plus*num_points
                   + 6*num_points
                   + j*(num_r_points+3)
                   + i; 
                }
#endif
 
                 double uu = recv_buff[mx];
                 double vv = recv_buff[my];
                 double ww = recv_buff[mz];
                 double pp = recv_buff[mp];

                 double uu1 = uu;
                 double vv1 = vv;
                 double ww1 = ww;
#ifdef MAP
#ifdef TMAP // time dependent mesh
                 uu1 +=  (omega->mapx->t[k] + omega->mapx->z[k] * ww1);
                 vv1 +=  (omega->mapy->t[k] + omega->mapy->z[k] * ww1);
#else // only Z-deformation
                 uu1  +=  omega->mapx->z[k] * ww1;
                 vv1  +=  omega->mapy->z[k] * ww1;
#endif
#endif
                 double ur   = cos(phi)*uu+sin(phi)*vv;
                 double uphi = -sin(phi)*uu+cos(phi)*vv;
                 double uz   = ww;

                 double ur1   = cos(phi)*uu1+sin(phi)*vv1;
                 double uphi1 = -sin(phi)*uu1+cos(phi)*vv1;
                 double uz1   = ww1;
#ifdef OUTPUT_VISCOSITY
                 double visc1 = recv_buff[ps1];
                 double visc2 = recv_buff[ps2];
#endif

                 hist_tmp[m][0] = ur1;
                 hist_tmp[m][1] = uphi1;
                 hist_tmp[m][2] = uz1;
                 hist_tmp[m][3] = pp;
                 hist_tmp[m][4] = ur;
                 hist_tmp[m][5] = uphi;
                 hist_tmp[m][6] = uz;
#ifdef OUTPUT_VISCOSITY
                 hist_tmp[m][7] = visc1;
                 hist_tmp[m][8] = visc2;
#else

               if(iparam("IPIPE_RMS")){
                 hist_tmp[m][7] = recv_buff[pv];
                 hist_tmp[m][8] = recv_buff[ps];
                 hist_tmp[m][9] = recv_buff[pf];
                }
#endif

          }

          }

          for(int i=0; i<num_hist_data; ++i)
            {
             fprintf(file_point," %lf  ",time);
              for(int n=0; n<nhp; ++n)
               fprintf(file_point," %lf  ", hist_tmp[i][n]);
               fprintf(file_point," \n");
            }

           free(hist_tmp);
           fflush(file_point);
           fclose(file_point);
           free(buf_point);
       }

   }

   gsync ();

   if(statistics_step == stastep)
   statistics_step=0;

   statistics_step++;
#ifdef MAP
  // FFT the maps back to fourier space
  rfftw(MPlan, 1, (FFTW_COMPLEX *) omega->mapx->z, 1, 0, 0, 0, 0);
  rfftw(MPlan, 1, (FFTW_COMPLEX *) omega->mapy->z, 1, 0, 0, 0, 0);
#ifdef TMAP
  rfftw(MPlan, 1, (FFTW_COMPLEX *) omega->mapx->t, 1, 0, 0, 0, 0);
  rfftw(MPlan, 1, (FFTW_COMPLEX *) omega->mapy->t, 1, 0, 0, 0, 0);
#endif
#endif

#ifdef OUTPUT_VISCOSITY
   omega->visc1->Set_state('t');
   omega->visc2->Set_state('t');
#endif

  free(V);
  free(X.x);
  free(data);
  free(recv_buff);
  free(ha);
  free(hb);
 //       free(hc);

}


#if 0
#ifdef MAP
void compute_rms_pressure(Domain *omega,Mapping *mapx, Mapping *mapy, int step)
{
    register int i,j,k;
    int          qa,qb,qc,m;
    int          nz = option("NZ"), nztot = option("NZTOT");
    int          plane;

    int          nprocs = option("NPROCS");
    int          procid = option("PROCID");
    int          hjtot  = omega->U->nz*omega->U->hjtot;
    int          htot   = omega->U->nz*omega->U->htot;
	  double       dt     = dparam("DELT");
	  double       dtinv  = 1./dt;

    Element      *eU,*eV,*eW,*eP,*eUt,*eVt,*eWt,*ePt,
                 *eUf, *eVf, *eWf;

    Element_List *U    =  omega->U,   *V    =  omega->V,    *W = omega->W,
                 *P    =  omega->P;

    Element_List *Ut   =  omega->Ut,  *Vt   =  omega->Vt,   *Wt = omega->Wt;
    Element_List *Uf   =  omega->Uf,  *Vf   =  omega->Vf,   *Wf = omega->Wf;
    Element_List *Pt   =  omega->Pt;


    double       *zerosv = (double *) NULL;

    double       *u0   =  omega->uk[1],     *v0    =  omega->vk[1],     *w0   = omega->wk[1], // velocity value at time step n
                 *u1   =  omega->uk[0],     *v1    =  omega->vk[0],     *w1   = omega->wk[0]; // velocity value at time step n-1

    double       *ut0  =  omega->ut[0],    *vt0   =  omega->vt[0],    *wt0 = omega->wt[0], // used as temps arrary for gradients
                 *ut1  =  omega->ut[1],    *vt1   =  omega->vt[1],    *wt1 = omega->wt[1]; // used for calculating laplacians

 //   double       **pipe_mean = omega->pipe_mean;

    double       *ux = omega->ut[0], *uy = omega->vt[0], *uz = omega->wt[0];
    double       *vx = omega->ut[1], *vy = omega->vt[1], *vz = omega->wt[1];
    double       *wx = omega->st[0], *wy = omega->st[1], *wz = omega->visc_tmp;


 //map stuff into physical space 
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->tz, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->tz, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->zz, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->zz, 1, 0, 0, 0, 0);

    Ut->Set_state('p');
    Vt->Set_state('p');
    Wt->Set_state('p');
    Pt->Set_state('p');


    char
    save_state_U = U->fhead->state,
    save_state_V = V->fhead->state,
    save_state_W = W->fhead->state,
    save_state_p = P->fhead->state,
    save_state_Uf = Uf->fhead->state,
    save_state_Vf = Vf->fhead->state,
    save_state_Wf = Wf->fhead->state;

    P->Set_state('p');
    U->Set_state('p');
    V->Set_state('p');
    W->Set_state('p');

    Uf->Set_state('p');
    Vf->Set_state('p');
    Wf->Set_state('p');


    dcopy(htot, u0, 1, Ut->base_h, 1);
    dcopy(htot, v0, 1, Vt->base_h, 1);
    dcopy(htot, w0, 1, Wt->base_h, 1);


//    Ut->Trans(Ut, F_to_P); // transformed to (Q,P)-space
//    Vt->Trans(Vt, F_to_P);
//    Wt->Trans(Wt, F_to_P);


    dvsub (htot, w0, 1, w1,  1, Pt->base_h, 1);
    dsmul (htot, dtinv, Pt->base_h,  1, Pt->base_h, 1); 
//    Pt->Trans(Pt, F_to_P);  //Pt has the (uz^(n)-uz^(n-1))/dt
    Pt ->Grad_d(vx,vy,NULL,'a'); //d/dx d/dy

    Wt->Grad_d(wx,wy,NULL,'a'); //duz/dx duz/dy
    Wt->Grad_z(P);     
    dcopy(htot, P->base_h, 1, wz, 1); //duz/dz

    dvmul  (htot, Ut->base_h, 1, wx, 1, P->base_h, 1);
    dvvtvp (htot, Vt->base_h, 1, wy, 1, P->base_h, 1, P->base_h,1);
    dvvtvp (htot, Wt->base_h, 1, wz, 1, P->base_h, 1, P->base_h,1); //P has the (u*\grad w)

    P->Grad_d(ux,uy,NULL,'a'); //dw/dx dw/dy P->base_h is free to use

    qa = Ut->flevels[0]->fhead->qa,
    qb = Ut->flevels[0]->fhead->qb;
  	Coord X;
	  X.x  = dvector(0, qa*qb-1);
	  X.y  = dvector(0, qa*qb-1);

    double  **D = dmatrix(0,1,0,QGmax*QGmax-1);
    double  *uz_x = D[0], *uz_y = D[1];

    int     nf = 3;
    double  **avgz   = dmatrix(0,nf-1,0,Ut->htot);
    double  *avg_tmp = dvector(0,Ut->htot-1);
     
    for(j=0; j<Ut->htot; ++j)
       for(int d = 0; d < nf; ++d)
        avgz[d][j] = 0.;
    
//    double       **pipe_mean = dmatrix(0,3,0, Ut->htot);

    for (k = 0, plane=nz*procid; k < nz; k++,plane++)
     {
     for(int eid =0; eid<Ut->nel; ++eid)
      {

	     eUt   = Ut->flevels[k]->flist[eid];
	     eVt   = Vt->flevels[k]->flist[eid];
	     eWt   = Wt->flevels[k]->flist[eid];
	     ePt   = Pt->flevels[k]->flist[eid];

       qa = eUt->qa;   qb = eUt->qb;
	     eUt->coord(&X);

       int offset = eid*qa*qb;
       int offset1 = Ut->htot*k+eid*qa*qb;
       for ( int i=0; i<qa; i++ )
        for ( int j=0; j<qb; j++ )
          {
            m = i*qb + j;
            double phi = atan2(X.x[m],X.y[m]);
            
            double ur   = cos(phi)*eUt->h[0][m]+sin(phi)*eVt->h[0][m];
            double uphi = -sin(phi)*eUt->h[0][m]+cos(phi)*eVt->h[0][m];
            double uz   = eWt->h[0][m];

            double a_rho_z = cos(phi)*mapx->z[plane]+sin(phi)*mapy->z[plane];
            double a_phi_z = cos(phi)*mapy->z[plane]-sin(phi)*mapx->z[plane];

            double a_rho_tz = cos(phi)*mapx->tz[plane]+sin(phi)*mapy->tz[plane];
            double a_phi_tz = cos(phi)*mapy->tz[plane]-sin(phi)*mapx->tz[plane];

            double a_rho_zz = cos(phi)*mapx->zz[plane]+sin(phi)*mapy->zz[plane];
            double a_phi_zz = cos(phi)*mapy->zz[plane]-sin(phi)*mapx->zz[plane];

            double vib_prms_0 = (cos(phi)*vx[offset1+m]+sin(phi)*vy[offset1+m])*a_rho_z
                                +(-sin(phi)*vx[offset1+m]+cos(phi)*vy[offset1+m])*a_phi_z;

            double  vib_prms_1 = 2.*(cos(phi)*wx[offset1+m]+sin(phi)*wy[offset1+m])*a_rho_tz
                               +2.*(-sin(phi)*wx[offset1+m]+cos(phi)*wy[offset1+m])*a_phi_tz;

            double  vib_prms_2 = (cos(phi)*ux[offset1+m]+sin(phi)*uy[offset1+m])*a_rho_z
                                +(-sin(phi)*ux[offset1+m]+cos(phi)*uy[offset1+m])*a_phi_z;

            double  vib_prms_3 = 2.*uz*(cos(phi)*wx[offset1+m]+sin(phi)*wy[offset1+m])*a_rho_zz
                                +2.*uz*(-sin(phi)*wx[offset1+m]+cos(phi)*wy[offset1+m])*a_phi_zz;

            U->base_h[offset1+m] = vib_prms_0;
            V->base_h[offset1+m] = vib_prms_1;
            W->base_h[offset1+m] = vib_prms_2;
            P->base_h[offset1+m] = vib_prms_3;

//            ROOTONLY fprintf(stderr,"vib0 = %lf  vib1 = %lf  vib2 = %lf  vib3 = %lf \n",vib_prms_0, vib_prms_1,vib_prms_2,vib_prms_3 );

            Uf->base_h[offset1+m] = vib_prms_0+vib_prms_1+vib_prms_2+vib_prms_3;  //Uf->hase has the vibration part  
//            Uf->base_h[offset1+m] = 0.;  //Uf->hase has the vibration part  

//            eUt->h[0][m] = ur - pipe_mean[0][offset+m];  //fluctuation
//            eVt->h[0][m] = uphi - pipe_mean[1][offset+m];
//            eWt->h[0][m] = uz - pipe_mean[2][offset+m];
//
//            Pt->base_h[offset1+m] = pipe_mean[0][offset+m]; //mean uz save for d(uz)/dz
            avgz[0][offset+m] += eUt->h[0][m];
            avgz[1][offset+m] += eVt->h[0][m];
            avgz[2][offset+m] += eWt->h[0][m];
          }
      }
     } //U,V,W have ur_prime. uphi_prime, uz_prime


    for(int d = 0; d < nf; ++d)
     {
       MPI_Allreduce (&avgz[d][0], avg_tmp, Ut->htot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
       dcopy(Ut->htot,avg_tmp,1,avgz[d],1);
     }

    int tmp = 0;
    int total_nz;
    MPI_Allreduce (&nz, &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    icopy(1,&tmp,1,&total_nz,1);

    for(i=0; i<Ut->htot; ++i)
    {
      avgz[0][i] /= (double)total_nz; //now avgz[d] have mean U_r,U_phi,U_z,P, in (Q,P)-space at each points
      avgz[1][i] /= (double)total_nz;
      avgz[2][i] /= (double)total_nz;
     for (k=0; k<nz; ++k)
      {
       Ut->base_h[k*Ut->htot+i] -= avgz[0][i];
       Vt->base_h[k*Ut->htot+i] -= avgz[1][i];
       Wt->base_h[k*Ut->htot+i] -= avgz[2][i];
      }
    }

    Ut->Grad_d(ux,uy,NULL,'a'); //du/dx du/dy
    Ut->Grad_z(Vf);     
    dcopy(htot, Vf->base_h, 1, uz, 1); //du/dz

    Vt->Grad_d(vx,vy,NULL,'a');  //dv/dx dv/dy
    Vt->Grad_z(Vf); 
    dcopy(htot, Vf->base_h, 1, vz, 1); //dv/dz

    Wt->Grad_d(wx,wy,NULL,'a');    //dw/dx dw/dy
    Wt->Grad_z(Vf);
    dcopy(htot, Vf->base_h, 1, wz, 1); //dw/dz

//Vf->base, Wf->base are free to use
    for (int k = 0; k < nz; k++)
     {
     for(int eid =0; eid<Ut->nel; ++eid)
      {

	     eUt   = Ut->flevels[k]->flist[eid];
	     eVt   = Vt->flevels[k]->flist[eid];
	     eWt   = Wt->flevels[k]->flist[eid];

	     ePt   = Pt->flevels[k]->flist[eid];
       
       qa = eUt->qa;   qb = eUt->qb;
	     eUt->coord(&X);
      
       for (m=0; m<qa*qb; ++m)
        ePt->h[0][m] = avgz[2][eid*qa*qb+m]; //copy mean uz

	     ePt->Grad_d(uz_x, uz_y, 0, 'a');

       int offset = k*Ut->htot+eid*qa*qb;
//	     
       for ( int i=0; i<qa; i++ )
        for ( int j=0; j<qb; j++ )
          {
            m = i*qb + j;
            double phi = atan2(X.x[m],X.y[m]);
            double rho = sqrt(X.x[m]*X.x[m]+X.y[m]*X.y[m]);
            
            double ur_prime   = cos(phi)*eUt->h[0][m]+sin(phi)*eVt->h[0][m];
            double uphi_prime = cos(phi)*eVt->h[0][m]-sin(phi)*eUt->h[0][m];
            double uz_prime   = eWt->h[0][m];
            
            double dur_dr = cos(phi)*cos(phi)*ux[offset+m]+cos(phi)*sin(phi)*vx[offset+m]
                           +sin(phi)*cos(phi)*uy[offset+m]+sin(phi)*sin(phi)*vy[offset+m];
//omit factor rho
            double duf_df = -sin(phi)*cos(phi)*vx[offset+m]+sin(phi)*sin(phi)*ux[offset+m]
                           +cos(phi)*cos(phi)*vy[offset+m]-sin(phi)*cos(phi)*uy[offset+m];
//omit factor rho
            double dur_df = -sin(phi)*cos(phi)*ux[offset+m]-sin(phi)*sin(phi)*vx[offset+m]
                           +cos(phi)*cos(phi)*uy[offset+m]+sin(phi)*cos(phi)*vy[offset+m];
             
            double duf_dr = cos(phi)*cos(phi)*vx[offset+m]-sin(phi)*cos(phi)*ux[offset+m]
                           +sin(phi)*cos(phi)*vy[offset+m]-sin(phi)*sin(phi)*uy[offset+m];

            double slow_prms1 = ux[offset+m]*ux[offset+m];
            double slow_prms2 = vy[offset+m]*vy[offset+m];
            double slow_prms3 = wz[offset+m]*wz[offset+m];

            double slow_prms4 = 2.*vx[offset+m]*uy[offset+m];
            double slow_prms5 = 2.*wx[offset+m]*uz[offset+m];
            double slow_prms6 = 2.*wy[offset+m]*vz[offset+m];
/*
            double slow_prms1 = dur_dr*dur_dr;
            double slow_prms2 = duf_df*duf_df;
            double slow_prms3 = wz[offset+m]*wz[offset+m];
            double slow_prms4 = 2.*dur_df*duf_dr;
            double slow_prms5 = 2.*(cos(phi)*uz[offset+m]+sin(phi)*vz[offset+m])*(cos(phi)*wx[offset+m]+sin(phi)*wy[offset+m]);
            double slow_prms6 = 2.*(cos(phi)*vz[offset+m]-sin(phi)*uz[offset+m])*(-sin(phi)*wx[offset+m]+cos(phi)*wy[offset+m]);
            double slow_prms7 = ur_prime*ur_prime/(rho*rho+1e-12);
            double slow_prms8 = 2.*ur_prime*duf_df/(rho+1e-12);
            double slow_prms9 = -2.*uphi_prime*duf_dr/(rho+1e-12);
*/

            double fast_prms = 2.*(cos(phi)*uz[offset+m]+sin(phi)*vz[offset+m])*(cos(phi)*uz_x[m]+sin(phi)*uz_y[m]);


//           ROOTONLY fprintf(stderr,"rho = %lf  slow1 = %lf  slow2 = %lf  slow3 = %lf  slow4 = %lf  slow5 = %lf  slow6 = %lf  offset = %d  \n"
//              ,rho ,slow_prms1, slow_prms2,slow_prms3,slow_prms4,slow_prms5,slow_prms6, offset+m);

           Vf->base_h[offset+m] = slow_prms1+slow_prms2+slow_prms3+slow_prms4+slow_prms5+slow_prms6;  //Vf->base_h
           Wf->base_h[offset+m] = fast_prms;  //Wf->base_h

          }
      }
     } //U,V,W have ur_prime. uphi_prime, uz_prime
    
    dcopy(htot, Uf->base_h, 1, Ut->base_h, 1);
    dcopy(htot, Vf->base_h, 1, Vt->base_h, 1);
    dcopy(htot, Wf->base_h, 1, Wt->base_h, 1);

    Ut->Set_state('p');
    Ut->Trans(Ut, P_to_F); // transformed to (Q,F)-space
    Ut->Iprod(Ut);
    Ut->Set_state('t');


    Vt->Set_state('p');
    Vt->Trans(Vt, P_to_F); // transformed to (Q,F)-space
    Vt->Iprod(Vt);
    Vt->Set_state('t');

    Wt->Set_state('p');
    Wt->Trans(Wt, P_to_F); // transformed to (Q,F)-space
    Wt->Iprod(Wt);
    Wt->Set_state('t');
//now we have all the RHS of pressure
 //setup pressure boundary   
    RMS_SetPBCs(omega);

    for(i=0; i<Ut->htot; ++i)
    ROOT fprintf(stderr,"Before Ut = %lf  Vt = %lf  Wt = %lf   i = %d \n",Ut->base_h[i],Vt->base_h[i],Wt->base_h[i],i);

     solve (Pt,Ut,omega->Pbc,omega->Pressure_sys,Helm,step,zerosv); //vibration part

     Pt->Trans(Pt, J_to_Q); // vibration pressure rms transformed to (Q,F)-space
     Pt->Set_state('p');
     Pt->Trans(Pt, F_to_P); // transformed to (Q,P)-space

     for(i=0; i<500; ++i)
     ROOT fprintf(stderr,"After Ptj = %lf  Pt = %lf i = %d \n",Pt->base_hj[i],Pt->base_h[i], i);

     dcopy(htot, Pt->base_h, 1, Uf->base_h, 1);

    solve (Pt,Vt,omega->Pbc,omega->Pressure_sys,Helm,step,zerosv); //vibration part
    Pt->Trans(Pt, J_to_Q); // vibration pressure rms transformed to (Q,F)-space
    Pt->Set_state('p');
    Pt->Trans(Pt, F_to_P); // transformed to (Q,P)-space
    dcopy(htot, Pt->base_h, 1, Vf->base_h, 1);

    solve (Pt,Wt,omega->Pbc,omega->Pressure_sys,Helm,step,zerosv); //fast part
    Pt->Trans(Pt, J_to_Q); // vibration pressure rms transformed to (Q,F)-space
    Pt->Set_state('p');
    Pt->Trans(Pt, F_to_P); // transformed to (Q,P)-space
    dcopy(htot, Pt->base_h, 1, Wf->base_h, 1);

    for(i=0; i<Ut->htot; ++i)
    ROOT fprintf(stderr,"After Ut = %lf  Vt = %lf  Wt = %lf   i = %d \n",Uf->base_h[i],Vf->base_h[i],Wf->base_h[i],i);


//maps go back to  Fourier space   
    rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->z,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->z,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->tz,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->tz,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->zz, 1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->zz, 1, 0, 0, 0, 0);

    U->Set_state(save_state_U);
    V->Set_state(save_state_V);
    W->Set_state(save_state_W);
    P->Set_state(save_state_W);

    Uf->Set_state(save_state_Uf);
    Vf->Set_state(save_state_Vf);
    Wf->Set_state(save_state_Wf);

    free_dmatrix(D,0,0);
    free_dmatrix(avgz,0,0);
    free(zerosv);
    free(avg_tmp);
}

#endif

#endif

void compute_cylindrical_vorticity(Domain *omega,int step)
{
    register int i,j,k;
    int          qa,qb,qc,m;
    int          nz = option("NZ"), nztot = option("NZTOT");
    int          plane;

    int          nprocs = option("NPROCS");
    int          procid = option("PROCID");
    int          hjtot  = omega->U->nz*omega->U->hjtot;
    int          htot   = omega->U->nz*omega->U->htot;
	  double       dt     = dparam("DELT");
	  double       dtinv  = 1./dt;

    Element      *eUt,*eVt,*eWt,*ePt;

    Element_List *Ut   =  omega->Ut,  *Vt   =  omega->Vt,   *Wt = omega->Wt;
    Element_List *Pt   =  omega->Pt;

    double       *u0   =  omega->uk[1],     *v0    =  omega->vk[1],     *w0   = omega->wk[1], // velocity value at time step n
                 *u1   =  omega->uk[0],     *v1    =  omega->vk[0],     *w1   = omega->wk[0]; // velocity value at time step n-1

    double       *ut0  =  omega->ut[0],    *vt0   =  omega->vt[0],    *wt0 = omega->wt[0], // used as temps arrary for gradients
                 *ut1  =  omega->ut[1],    *vt1   =  omega->vt[1],    *wt1 = omega->wt[1]; // used for calculating laplacians

 //   double       **pipe_mean = omega->pipe_mean;

    double       *ux = omega->ut[0], *uy = omega->vt[0], *uz = omega->wt[0];
    double       *vx = omega->ut[1], *vy = omega->vt[1], *vz = omega->wt[1];
    double       *wx = omega->st[0], *wy = omega->st[1], *wz = omega->visc_tmp;


    Ut->Set_state('p');
    Vt->Set_state('p');
    Wt->Set_state('p');
    Pt->Set_state('p');

    dcopy(htot, u0, 1, Ut->base_h, 1);
    dcopy(htot, v0, 1, Vt->base_h, 1);
    dcopy(htot, w0, 1, Wt->base_h, 1);


//    Ut->Trans(Ut, F_to_P); // transformed to (Q,P)-space
//    Vt->Trans(Vt, F_to_P);
//    Wt->Trans(Wt, F_to_P);


    Ut->Grad_d(ux,uy,NULL,'a'); //du/dx du/dy

    Ut->Trans(Ut, P_to_F); // transformed to Fourier space
    Ut->Grad_z(Pt);     
    Pt->Trans(Pt, F_to_P); // transformed to Fourier space
    dcopy(htot, Pt->base_h, 1, uz, 1); //du/dz

    Ut->Trans(Ut, F_to_P); // transformed to Fourier space

    Vt->Grad_d(vx,vy,NULL,'a');  //dv/dx dv/dy
    Vt->Trans(Vt, P_to_F); // transformed to Fourier space
    Vt->Grad_z(Pt); 
    Pt->Trans(Pt, F_to_P); // transformed to Fourier space
    dcopy(htot, Pt->base_h, 1, vz, 1); //dv/dz

    Vt->Trans(Vt, F_to_P); // transformed to Fourier space

    Wt->Grad_d(wx,wy,NULL,'a');  //dw/dx dw/dy
    Wt->Trans(Wt, P_to_F); // transformed to Fourier space
    Wt->Grad_z(Pt); 
    Pt->Trans(Pt, F_to_P); // transformed to Fourier space
    dcopy(htot, Pt->base_h, 1, wz, 1); //dw/dz

    Wt->Trans(Wt, F_to_P); // transformed to Fourier space

    qa = Ut->flevels[0]->fhead->qa,
    qb = Ut->flevels[0]->fhead->qb;
  	Coord X;
	  X.x  = dvector(0, qa*qb-1);
	  X.y  = dvector(0, qa*qb-1);

    for (int k = 0; k < nz; k++)
    {
     for(int eid =0; eid<Ut->nel; ++eid)
      {
	     eUt   = Ut->flevels[k]->flist[eid];
	     eVt   = Vt->flevels[k]->flist[eid];
	     eWt   = Wt->flevels[k]->flist[eid];
	     ePt   = Pt->flevels[k]->flist[eid];
       
       qa = eUt->qa;   qb = eUt->qb;
	     eUt->coord(&X);

       int offset = k*Ut->htot+eid*qa*qb;
//	     
       for ( int i=0; i<qa; i++ )
        for ( int j=0; j<qb; j++ )
          {
            m = i*qb + j;
            double phi = atan2(X.y[m],X.x[m]);
            double rho = sqrt(X.x[m]*X.x[m]+X.y[m]*X.y[m]);

            double ur_prime   = cos(phi)*eUt->h[0][m]+sin(phi)*eVt->h[0][m];
            double uphi_prime = cos(phi)*eVt->h[0][m]-sin(phi)*eUt->h[0][m];
            double uz_prime   = eWt->h[0][m];
            
            double dur_dr = cos(phi)*cos(phi)*ux[offset+m]+cos(phi)*sin(phi)*vx[offset+m]
                           +sin(phi)*cos(phi)*uy[offset+m]+sin(phi)*sin(phi)*vy[offset+m];
//omit factor rho
            double duf_df = -sin(phi)*cos(phi)*vx[offset+m]+sin(phi)*sin(phi)*ux[offset+m]
                           +cos(phi)*cos(phi)*vy[offset+m]-sin(phi)*cos(phi)*uy[offset+m];
//omit factor rho
            double dur_df = -sin(phi)*cos(phi)*ux[offset+m]-sin(phi)*sin(phi)*vx[offset+m]
                           +cos(phi)*cos(phi)*uy[offset+m]+sin(phi)*cos(phi)*vy[offset+m];
             
            double duf_dr = cos(phi)*cos(phi)*vx[offset+m]-sin(phi)*cos(phi)*ux[offset+m]
                           +sin(phi)*cos(phi)*vy[offset+m]-sin(phi)*sin(phi)*uy[offset+m];

            eUt->h[0][m] = -sin(phi)*wx[offset+m]+cos(phi)*wy[offset+m]
                           -cos(phi)*vz[offset+m]+sin(phi)*uz[offset+m];

            eVt->h[0][m] = cos(phi)*uz[offset+m]+sin(phi)*vz[offset+m]
                           -cos(phi)*wx[offset+m]-sin(phi)*wy[offset+m];

           if(rho>1e-8)
            eWt->h[0][m] = uphi_prime/rho+(duf_dr-dur_df);
           else
            eWt->h[0][m] = 0.;


          }
      }
    }

 }

void compute_local_vorticity_flux(Domain *omega, int step, double time)
 {
  Element_List *Ut   =  omega->Ut,   *Vt   =  omega->Vt,   *Wt   = omega->Wt;
  Element_List *Pt   =  omega->Pt,   *P    =  omega->P;
  Element *eUt, *eVt, *eWt,*ePt,*eP;
  register int i,j,k;
  int      nz = option("NZ"), nztot = option("NZTOT");
  int      nprocs = option("NPROCS");
  int      myid   = option("PROCID");
  int      hjtot  = Ut->nz*Ut->hjtot;
  int      htot   = Ut->nz*Ut->htot;
  int      eid,face,lmax=0, trip = 0;
  int      count = 0;
  Basis    *b;

  double   radius = dparam("RADIUS");

  if(radius <1e-10)
   {
     fprintf(stderr,"Need to prescribe a value to parameter RADIUS \n");
     exit(1);
   }

  double  *p0   =  omega->pk[1]; //pressure in (Q,P)-space at step n
  double  *ut0  =  omega->ut[0], *vt0  =  omega->vt[0],  *wt0   = omega->wt[0], // working space
          *ut1  =  omega->ut[1], *vt1  =  omega->vt[1],  *wt1   = omega->wt[1]; // working space

  Bndry   *UB = omega->Ubc[0],
          *PB = omega->Pbc[0];

  int      qa = Ut->flevels[0]->flist[0]->qa;
  int      qb = Ut->flevels[0]->flist[0]->qb;

  double  *w;
  double  *wk = dvector(0,max(LGmax+1,QGmax)*QGmax-1);
  double  **D = dmatrix(0,9,0,QGmax*QGmax-1);
  double  *ux = D[0], *uy = D[1];
  double  *vx = D[2], *vy = D[3];
  double  *wx = D[4], *wy = D[5];
  double  *px = D[6], *py = D[7];
  double  *pp = D[8], *pz = D[9];
  
  Coord X;
	X.x  = dvector(0, qa*qb-1);
	X.y  = dvector(0, qa*qb-1);
  
  /* write out header */
  if(pres_init){

  
  /* need to set up surface geometric factors in T from P */
  for(k = 0; k < nz; ++k)
    set_SurGeofac(omega->Pbc[k],omega->Ubc[k]);
  // count how many points on the surface
  count = 0;
  nqp = 0;
//count the number of faces 
  for(UB=omega->Ubc[0];UB; UB = UB->next)
     if(UB->type == 'W' )
      {
         nqp += qa;
         count ++;
      }


    ROOT fprintf(stdout,"%d Surfaces are found \n",count);

    map_vector     = ivector(0,nqp-1);
    coord_vector   = dvector(0,nqp-1);
    coord_x        = dvector(0,nqp-1);
    coord_y        = dvector(0,nqp-1);

    count = 0;
    for(UB=omega->Ubc[0];UB; UB = UB->next)
     if(UB->type == 'W' ) {
       face = UB->face;
       eid  = UB->elmt->id;

      for(k = 0; k < 1; ++k){
	     ePt = Pt->flevels[k]->flist[eid];

	     ePt->coord(&X);

	     ePt->GetFace(X.x ,face,wk);
	     ePt->InterpToFace1(face,wk,X.x);
	     ePt->GetFace(X.y ,face,wk);
	     ePt->InterpToFace1(face,wk,X.y);
	     
       for(i = 0; i < qa; ++i){
         coord_vector[count*qa+i] = atan2(X.y[i],X.x[i]);
         coord_x[count*qa+i] = X.x[i];//normal outward
         coord_y[count*qa+i] = X.y[i];//normal outward
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

    pres_init = 0;
  }

  double *workspace = dvector(0,10*nqp*nztot-1);
  dzero(10*nqp*nztot,workspace,1);

  double *send   = workspace;
  double *recv   = workspace+5*nqp*nztot-1;

  double *omega_theta_dr = workspace;
  double *omega_z_dr     = workspace+nqp*nz;
  double *pres_dtheta    = workspace+2*nqp*nz;
  double *pres           = workspace+3*nqp*nz;
  double *pres_dz        = workspace+4*nqp*nz;

  compute_cylindrical_vorticity(omega,step);

  char save_state_P = P->fhead->state;
   P->Set_state('p');

  Pt->Set_state('p');
  dcopy(htot, p0, 1, Pt->base_h, 1); //copy pressure to P->base_h
  Pt->Trans(Pt, P_to_F); // transformed to Fourier space
  Pt->Grad_z(P); //gradient on z direction
  P->Trans(P, F_to_P); // transformed to Fourier space

  for(k = 0; k < nz; ++k){
    int count = 0;
   for(UB=omega->Ubc[k];UB; UB = UB->next)
     if(UB->type == 'W'){
      face = UB->face;
      eid  = UB->elmt->id;

      if(Ut->flist[eid]->lmax != lmax){
	       b    = Ut->flist[eid]->getbasis();
	       lmax = Ut->flist[eid]->lmax;
	       getzw(qa,&w,&w,'a');
      }

	     eUt = Ut->flevels[k]->flist[eid];
	     eVt = Vt->flevels[k]->flist[eid];
	     eWt = Wt->flevels[k]->flist[eid];
	     ePt = Pt->flevels[k]->flist[eid];
	     eP  = P->flevels[k]->flist[eid];

//	     eUt->Jbwd(eUt,b);  eUt->state = 't'; //J to Q
       eUt->Grad_d(ux, uy, 0, 'a'); //omega_r
       eVt->Grad_d(vx, vy, 0, 'a'); //omega_theta
       eWt->Grad_d(wx, wy, 0, 'a'); //omega_z
       ePt->Grad_d(px, py, 0, 'a'); //p
	     
	     for(i = 0; i < 8; ++i){
	        eUt->GetFace(D[i] ,face,wk);
	        eUt->InterpToFace1(face,wk,D[i]);
	      }

	     ePt->GetFace(ePt->h[0],face,wk);
	     ePt->InterpToFace1(face,wk,pp);

	     eP->GetFace(eP->h[0],face,wk);
	     eP->InterpToFace1(face,wk,pz);

	      for(i = 0; i < qa; ++i){
            double xx = coord_x[map_vector[count*qa+i]];
            double yy = coord_y[map_vector[count*qa+i]];
            
            double theta = atan2(yy,xx);
            double rho = sqrt(xx*xx+yy*yy);

            double local_omega_c_dr = vx[i]*cos(theta)+sin(theta)*vy[i];
            double local_omega_z_dr = wx[i]*cos(theta)+sin(theta)*wy[i];
            double local_p_dc = rho*(-px[i]*sin(theta)+cos(theta)*py[i]);

//             fprintf(stderr,"xx = %g yy = %g  theta = %g  local_omega_c_dr = %g  local_omega_z_dr = %g  local_p_dc = %g  pres_z = %g   pres =  %g  count = %d \n",
//                          xx,yy,theta,local_omega_c_dr, local_omega_z_dr,local_p_dc, pz[i], pp[i],count);

	          omega_theta_dr[k*nqp+map_vector[count*qa+i]] = local_omega_c_dr;
	          omega_z_dr[k*nqp+map_vector[count*qa+i]] = local_omega_z_dr;
            pres_dtheta[k*nqp+map_vector[count*qa+i]]=local_p_dc;
            pres[k*nqp+map_vector[count*qa+i]]=pp[i];
            pres_dz[k*nqp+map_vector[count*qa+i]] = pz[i];

           }

       count++;
     }
        
    }

   
   MPI_Allgather(&send[0], 5*nqp*nz, MPI_DOUBLE, 
				 &recv[0], 5*nqp*nz, MPI_DOUBLE, MPI_COMM_WORLD);
   
   omega_theta_dr = workspace;
   omega_z_dr     = workspace+nqp*nztot;
   pres_dtheta    = workspace+2*nqp*nztot;
   pres           = workspace+3*nqp*nztot;
   pres_dz        = workspace+4*nqp*nztot;

  for(i = 0; i < nprocs; ++i)
	 for (k = 0; k < nz; k++)
	  for (j = 0; j < nqp; j++)
     {
	     omega_theta_dr[i*nz*nqp+k*nqp+j] = recv[i*5*nqp*nz +k*nqp+j];
	     omega_z_dr[i*nz*nqp+k*nqp+j]     = recv[i*5*nqp*nz +k*nqp+j+nz*nqp];
	     pres_dtheta[i*nz*nqp+k*nqp+j]    = recv[i*5*nqp*nz +k*nqp+j+2*nz*nqp];
	     pres[i*nz*nqp+k*nqp+j]           = recv[i*5*nqp*nz +k*nqp+j+3*nz*nqp];
	     pres_dz[i*nz*nqp+k*nqp+j]        = recv[i*5*nqp*nz +k*nqp+j+4*nz*nqp];
     }

    ROOT {
        
        int nqp_per_face = (iparam("LQUAD"))? iparam("LQUAD"):iparam("MODES")+1;
        int num_face = nqp/nqp_per_face;

        FILE *fout;
	      char *buff = (char*) calloc(BUFSIZ,sizeof(char));
	      sprintf(buff, "wall_vorticity.txt");
        fout = fopen(buff,"aw");

        for(k = 0; k < nztot; ++k)
         for(i=0; i<num_face; ++i)
          {
            int q1 = i*nqp_per_face+nqp_per_face-1;
            int q2 = i*nqp_per_face+nqp_per_face-2;

            double pz   =  0.5*(pres_dz[k*nqp+q1]+pres_dz[k*nqp+q2]);
            double omtr =  0.5*(omega_theta_dr[k*nqp+q1]+omega_theta_dr[k*nqp+q2]);
            double pt   =  0.5*(pres_dtheta[k*nqp+q1]+pres_dtheta[k*nqp+q2]);
            double omzr =  0.5*(omega_z_dr[k*nqp+q1]+omega_z_dr[k*nqp+q2]);
            double ps   =  0.5*(pres[k*nqp+q1]+pres[k*nqp+q2]);
            
            fprintf(fout,"%lf  %lf  %lf  %lf  %lf   %lf  %lf  %d \n",time, coord_vector[q1], pz,omtr,pt,omzr,ps,i);
          }

//        for(k = 0; k < nztot; ++k)
//         for(i=0; i<nqp; ++i)
//          fprintf(stderr," pres = %lf  pres_z = %lf   i= %d   k= %d\n",pres[k*nqp+i],pres_dz[k*nqp+i],i, k);

        fflush(fout);
        fclose(fout);
        free(buff);

    }

  P->Set_state(save_state_P);

  free(wk); free_dmatrix(D,0,0);
  free(X.x);free(X.y);
  free(workspace);
 }
