#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <polylib.h>
#include "veclib.h"
#include "hotel.h"
#include "nekstruct.h"
#include "Tri.h"
#include "Quad.h"
#include "nektarF.h"
#ifdef MAP
#include "map.h"
#endif
// ce107
#include <rfftw.h>

// global FFT structures
rfftw_plan rplan, rplan_inv;
rfftw_plan rplan32, rplan_inv32;

#ifdef FULL_FIELD
#define nrtot ntot
#endif

Element *ee;

#define WARN fprintf(stderr,"Call to void Element function\n");
#define Eloop(E) for(ee=E;ee;ee=ee->next)
#define Zloop(i) for(i = 0; i < nz; ++i)

Fourier_List::Fourier_List(){
  fhead   = (Element*)NULL;
  nel     = 0;
  hjtot   = 0;
  htot    = 0;
  base_h  = (double *)0;
  base_hj = (double *)0;
  
  nz      = 1;
  nztot   = 1;
  flevels = (Element_List**)NULL;
}

Fourier_List::Fourier_List(Element **hea, int n){
  flist  = hea;
  fhead  = *hea;
  nel    = countelements(*hea);
#ifdef DEBUG
  if(n!=nel)
    fprintf(stderr,"Fourier_List::Fourier_List error list incomplete\n");
#endif
  hjtot  = 0;
  htot   = 0;
  base_h = (double *)0;
  base_hj= (double *)0;

  nz      = 1;
  nztot   = 1;
  flevels = (Element_List**)NULL;
}

Element *Fourier_List::operator()(int i){
#ifdef DEBUG
  if(i>nel-1) 
    fprintf(stderr,"Fourier_List::operator() bounds error\n");
#endif
  return flist[i];
}

Element *Fourier_List::operator[](int i){
#ifdef DEBUG
  if(i>nel-1) 
    fprintf(stderr,"Fourier_List::operator() bounds error\n");
#endif
  return flist[i];
}

void Fourier_List::Cat_mem(){
  htot = 0;
  hjtot = 0;
  Eloop(fhead){
    htot  += ee->qa*ee->qb;
    hjtot += ee->Nmodes;
  }
  
  if (option("dealias") && (iparam("EQTYPE") == Convective)) nz = 3*nz/2;

  base_h  = dvector(0, nz*htot-1);
  dzero(nz*htot,   base_h, 1);

  nz = option("NZ"); // resets the value in case it's changed

  base_hj = dvector(0, nz*hjtot-1);
  dzero(nz*hjtot, base_hj, 1);

  Mem_shift(base_h, base_hj);
}


void Fourier_List::Mem_shift(double *new_h, double *new_hj){
  int i;
  base_h  = new_h;
  base_hj = new_hj;

  Zloop(i){
    flevels[i]->Mem_shift(new_h, new_hj);
    new_h  += htot;
    new_hj += hjtot;
    flevels[i]->htot = htot;
    flevels[i]->hjtot = hjtot;
  }
  //  new_hj -= hjtot;
  //  for(i = nz; i < 3*nz/2; ++i){
  //    flevels[i]->Mem_shift(new_h, new_hj);
  //    new_h  += htot;
  //    flevels[i]->htot = htot;
  //    flevels[i]->hjtot = hjtot;
  //  }    
}

void Fourier_List::Trans(Element_List *EL, Nek_Trans_Type ntt){
  // Treat FFT as a global operation
  switch (ntt) {
  case F_to_P:
  case P_to_F:
    FFT(EL, ntt);    
    break;
  case F_to_P32:
  case P_to_F32:
    FFT32(EL, ntt);
    break;
  default:
    int i;
    Zloop(i)
      flevels[i]->Trans(EL->flevels[i], ntt);
    break;
  }
}

void Fourier_List::Trans(int ntot, double *HL, double *EL, Nek_Trans_Type ntt)
{
  // Treat FFT as a global operation
  switch (ntt) {
  case F_to_P:
  case P_to_F:
    FFT(ntot, HL, EL, ntt);    
    break;
  case F_to_P32:
  case P_to_F32:
    FFT32(ntot, HL, EL, ntt);
    break;
  default:
    printf("Fourier_List::Trans() -- un-supported operations!\n");
    //error_msg(Fourier_List::Trans() -- unsupported operations!)
#if 0
    int i;
    Zloop(i)
      flevels[i]->Trans(EL->flevels[i], ntt);
#endif
    break;
  }
}

void Fourier_List::Iprod(Element_List *EL){
  int i;
  Zloop(i)
    flevels[i]->Iprod(EL->flevels[i]);
}

void Fourier_List::Grad(Element_List *AL, Element_List *BL, Element_List *CL,
			char Trip){
  if(nz > 1 && CL && (Trip == 'a' || Trip ==  'z'))
    Grad_z(CL);

  int i;
  if(AL && BL)
    Zloop(i)
      flevels[i]->Grad(AL->flevels[i], BL->flevels[i], 0, Trip);
  else if(AL)
    Zloop(i)
      flevels[i]->Grad(AL->flevels[i], 0, 0, Trip);
  else if(BL)
    Zloop(i)
      flevels[i]->Grad(0, BL->flevels[i], 0, Trip);
}

void Fourier_List::Grad_d(double *AL, double *BL, double *CL,
			char Trip){
  
  if(nz> 1 && CL && (Trip == 'a' || Trip ==  'z' || Trip == 'A' ||Trip == 'Z'))
    Grad_z_d(CL);

  int i;
  if (Trip == 'A' || Trip == 'X' || Trip == 'Y')
    nz = 3*nz/2;
  
  if(AL && BL)
    Zloop(i){
      flevels[i]->Grad_d(AL,BL,0,Trip);
      AL += htot; BL += htot;
    }
  else if(AL)
    Zloop(i){
      flevels[i]->Grad_d(AL, 0,0,Trip);
      AL += htot; 
    } 
  else if(BL)
    Zloop(i){
      flevels[i]->Grad_d(0, BL,0,Trip);
      BL += htot; 
    } 
  nz = option("NZ");
}

void Fourier_List::Grad_h(double *EL, double *AL, double *BL, double *CL,
			  char Trip){
  
  if(nz > 1 && CL && (Trip == 'a' || Trip ==  'z' || Trip == 'A' ||Trip == 'Z')) {
    fprintf(stderr,"Not implemented - not intended\n");
    exit (-1);
  }
  
  int i;
  if (Trip == 'A' || Trip == 'X' || Trip == 'Y') {
    nz = 3*nz/2;
    Trip = tolower(Trip);
  }

  if(AL && BL)
    Zloop(i){
      flevels[0]->Grad_h(EL,AL,BL,0,Trip);
      EL += htot; AL += htot; BL += htot;
    }
  else if(AL)
    Zloop(i){
      flevels[0]->Grad_h(EL, AL, 0,0,Trip);
      EL += htot; AL += htot; 
    } 
  else if(BL)
    Zloop(i){
      flevels[0]->Grad_h(EL, 0, BL,0,Trip);
      EL += htot; BL += htot; 
    } 
  nz = option("NZ");
}


void Fourier_List::HelmHoltz(Metric *lambda){
  fprintf(stderr,"Fourier_List::HelmHoltz not implemented\n");
}

void Fourier_List::Set_field(char *string){
  int i;
  dzero(nz*htot, base_h, 1);
  Zloop(i){
    dparam_set("z", zmesh(i));
    flevels[i]->Set_field(string);
  }
}

// Initialize RFFT structures
void init_rfft_struct()
{
  int nztot = option("NZTOT");

  rplan     = rfftw_create_plan(nztot, FFTW_FORWARD, 
				FFTW_MEASURE | FFTW_IN_PLACE, 
				REAL_TO_COMPLEX);
  rplan_inv = rfftw_create_plan(nztot, FFTW_BACKWARD, 
				FFTW_MEASURE | FFTW_IN_PLACE, 
				COMPLEX_TO_REAL);
  if (option("dealias")) {
    rplan32     = rfftw_create_plan(3*nztot/2, FFTW_FORWARD, 
				    FFTW_MEASURE | FFTW_IN_PLACE, 
				    REAL_TO_COMPLEX);
    rplan_inv32 = rfftw_create_plan(3*nztot/2, FFTW_BACKWARD, 
				    FFTW_MEASURE | FFTW_IN_PLACE, 
				    COMPLEX_TO_REAL);
  } else {
    rplan32     = rplan;
    rplan_inv32 = rplan_inv;
  }
  return;
}

#ifdef MAP
void Set_mapped_field(Element_List *EL, char *string, Mapping *mapx, Mapping *mapy){
  int i, nz = EL->nz, nztot = EL->nztot, htot = EL->htot;

  // invFFT the maps to physical space
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);

  dzero(nz*htot, EL->base_h, 1);
  Zloop(i){
    dparam_set("z", zmesh(i));
    dparam_set("ORIGINX", mapx->d[parid(i)]);
    dparam_set("ORIGINY", mapy->d[parid(i)]);
    dparam_set("t", mapy->time);
    EL->flevels[i]->Set_field(string);
  }
  // FFT the maps back to fourier space
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
}
#endif

void Fourier_List::Set_state(char st){
  int i;
  Zloop(i)
    flevels[i]->Set_state(st);
}

void Fourier_List::zerofield(){
  dzero(hjtot*nz, base_hj, 1);
  dzero(htot*nz,  base_h, 1);
}

Element_List* Fourier_List::gen_aux_field(char ty){
  Fourier_List *el = (Fourier_List*) new Fourier_List();

  el->nel         = nel;
  el->nz          = nz;
  el->nztot       = nztot;
  el->flevels     = (Element_List**) malloc(nz*sizeof(Element_List*));
  
  int i;
  Zloop(i){
    el->flevels[i] = flevels[i]->gen_aux_field(ty);
    el->htot = htot;
    el->hjtot = hjtot;
  }
    
  el->flist        = el->flevels[0]->flist;
  el->fhead        = el->flevels[0]->fhead;
  el->Cat_mem();

  return (Element_List*) el;
}

void Fourier_List::Terror(char *string){
  int i;
  Zloop(i)
    flevels[i]->Terror(string);
}

#ifdef MAP
void GTerror(Element_List *U, Mapping *mapx, Mapping *mapy, char *string)
#else
void GTerror(Element_List *U, char *string)
#endif
{
  int kloc,qt,trip=0, nz = U->nz, nztot = U->nztot;
  double li = 0.0, sob[2], gsob[2], areat = 0.0, l2a, h1a = 0.0,
         area,store;
  double *s,*utmp;
  Coord *X = (Coord *)malloc(sizeof(Coord));
  Element *E;
  register int k;

  if(string)
    vector_def("x y",string);
  
  X->x = dvector(0,QGmax*QGmax-1);
  X->y = dvector(0,QGmax*QGmax-1);
  utmp = dvector(0,QGmax*QGmax-1);

  sob[0]  = 0.;
  sob[1]  = 0.;
  gsob[0] = 0.;
  gsob[1] = 0.;

  for (k = 0; k < nz; k++) {
    areat = 0.;
    dparam_set("z", zmesh(k));
    kloc = parid(k);
    for(E=U->flevels[k]->fhead;E;E = E->next){
      qt = E->qtot;
      s  = E->h[0];
    
      if(E->state == 't'){ E->Trans(E,J_to_Q); trip = 1;}
      dcopy(qt,s,1,utmp,1);
    
      if(string){
	E->coord(X);
#ifdef MAP
	double orig;
	if ((orig = mapx->d[kloc]) != 0.0)
	    dsadd(qt, orig, X->x, 1, X->x, 1);
	if ((orig = mapy->d[kloc]) != 0.0)
	  dsadd(qt, orig, X->y, 1, X->y, 1);
#endif
	vector_set(qt,X->x,X->y,s);
	dvsub(qt,s,1,utmp,1,s,1);
      }

      li = ((store=E->Norm_li())>li)? store:li;
      E->Norm_l2m(&l2a,&area);
      E->Norm_h1m(&h1a,&area);
    
      sob[0] += l2a; 
      sob[1] += h1a;
      areat += area;
    
      if(trip){ E->state = 't'; trip = 0;}
      else     dcopy(qt,utmp,1,s,1);
    }
  }

  ROOT
    fprintf(stdout, "Field: %c ", U->fhead->type);

#ifdef PARALLEL    
  MPI_Reduce(&li, &store, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(sob, gsob, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  store = li;
  gsob[0] = sob[0];
  gsob[1] = sob[1];
#endif

  double dz = dparam("LZ")/nztot;
  ROOT
    fprintf(stdout,"%12.6lg %12.6lg %12.6lg (Linf  L2  H1)\n",
	    store,sqrt(dz*gsob[0]/areat),sqrt(dz*gsob[1]/areat));

  free(X->x); free(X->y); free(utmp); free((char *)X);

  return;
}

#ifndef MAP
void Compare (Domain *omega, ACTION space)
{
  Element_List *U = omega->U,
               *V = omega->V,
               *W = omega->W,
               *P = omega->P;
  Nek_Trans_Type f_to_p = F_to_P,
                 p_to_f = P_to_F;
  int revert = 0, revertP = 0, NZ = U->nz;

#ifdef MAP 
  Map     *mapx = omega->mapx,
          *mapy = omega->mapy;

  if (option("dealias")) { // need dealiasing for correct Mapfield 
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
      NZ = 3*NZ/2;
  }
#endif

  // make sure velocity fields are in quadrature space
  if (U->fhead->state == 't') { // all 3 velocity fields should be same type
    U->Trans(U, J_to_Q);
    V->Trans(V, J_to_Q);
    W->Trans(W, J_to_Q);
    revert++;
  }

  if (P->fhead->state == 't') {
    P->Trans(P, J_to_Q);
    revertP++;
  }

  if (space == Fourier) {
    U->Trans(U, f_to_p);
    V->Trans(V, f_to_p);
    W->Trans(W, f_to_p);
    P->Trans(P, F_to_P); // P doesn't need dealiasing
  }

#ifdef MAP
  Mapfield (omega, -1);

  if (option("dealias")) { // revert back to nz planes
    U->Trans(U, P_to_F32);
    U->Trans(U, F_to_P);
    V->Trans(V, P_to_F32);
    V->Trans(V, F_to_P);
    W->Trans(W, P_to_F32);
    W->Trans(W, F_to_P);
  }

  dparam_set("t", mapx->time);

  // invFFT the maps
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
#endif

#ifndef PLANAR_ERROR
#ifdef MAP
  GTerror(omega->U, omega->mapx, omega->mapy, omega->soln[0]);
  GTerror(omega->V, omega->mapx, omega->mapy, omega->soln[1]);
  GTerror(omega->W, omega->mapx, omega->mapy, omega->soln[2]);
  GTerror(omega->P, omega->mapx, omega->mapy, omega->soln[3]);
#else
  GTerror(omega->U, omega->soln[0]);
  GTerror(omega->V, omega->soln[1]);
  GTerror(omega->W, omega->soln[2]);
  GTerror(omega->P, omega->soln[3]);
#endif
#else
  for (int i = 0; i < U->nz; i++) {
    dparam_set("z", zmesh(i));
#ifdef MAP
    int iloc = parid(i);
    dparam_set("ORIGINX", mapx->d[iloc]);
    dparam_set("ORIGINY", mapy->d[iloc]);
#endif
    printf("Plane %d - ", iloc); U->flevels[i]->PTerror(omega->soln[0]); 
    printf("Plane %d - ", iloc); V->flevels[i]->PTerror(omega->soln[1]);
    printf("Plane %d - ", iloc); W->flevels[i]->PTerror(omega->soln[2]);
    printf("Plane %d - ", iloc); P->flevels[i]->PTerror(omega->soln[3]);
  }
#endif

#ifdef MAP
  // FFT the maps back
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);

  if (option("dealias")) { // Need to use 3/2 rule again
    U->Trans(U, P_to_F);
    U->Trans(U, F_to_P32);
    V->Trans(V, P_to_F);
    V->Trans(V, F_to_P32);
    W->Trans(W, P_to_F);
    W->Trans(W, F_to_P32);
  }

  Mapfield (omega, 1);
#endif

  if (space == Fourier) {
    U->Trans(U, p_to_f);
    V->Trans(V, p_to_f);
    W->Trans(W, p_to_f);
    P->Trans(P, P_to_F); // P doesn't need dealiasing
  }

  if (revert) { // return velocity fields to modal storage if necessary
    U->Set_state('t');
    V->Set_state('t');
    W->Set_state('t');
    //    U->Trans(U, Q_to_J);
    //    V->Trans(V, Q_to_J);
    //    W->Trans(W, Q_to_J);
  }

  if (revertP) // return pressure to modal storage if necessary
    P->Set_state('t');
    //    P->Trans(P, Q_to_J);

  return;
}
#endif

int parid(int i){
  int ret = option("NZ")*option("PROCID") + i;
  return ret;
}

double zmesh(int i){
  return  dparam("LZ")*(1.0*parid(i))/(option("NZTOT"));
}

double Beta2(int i){
#ifdef FLOK
  return  dparam("BETA")*dparam("BETA");
#else
#ifdef BETA0
  return 0;
#else
  return  (parid(i)==1) ? pow(2.0*M_PI*(option("NZTOT")/2)/dparam("LZ"),2):
    pow(2.0*M_PI*(parid(i)/2)/dparam("LZ"),2);
#endif    
#endif
}

double  Beta(int i){
#ifdef FLOK
  return  dparam("BETA");
#else
#ifdef BETA0
  return 0;
#else
  return  (parid(i)==1) ? 2.0*M_PI*(option("NZTOT")/2)/dparam("LZ"):
    2.0*M_PI*(parid(i)/2)/dparam("LZ"); 
#endif
#endif
}

int data_len(int nel,int *size, Element *E);
int count_facets(Element *E);
int countbcs(Domain *Omega);

double *Pget=(double*)0;
double *Pput=(double*)0;

void setup_transfer_space(Domain *omega)
{
  if(Pget) free(Pget);
  if(Pput) free(Pput);
  Element_List *U = omega->U;
  int nprocs = option("NPROCS");
  int fields = 3;
  int nz     = option("NZ");
  int NZ     = nz*nprocs;
  int nel    = countbcs(omega);
  int tot    = LGmax;/*Omega[0]->U->edge[0].l+2;*/
  int BC_ntot   = tot*nel*fields;
  int BC_nxy    = (BC_ntot+nprocs-1)/nprocs;
  int BC_ntotz  = BC_nxy*NZ;
  int nxy = (U->htot + nprocs - 1) / nprocs;
  int ntotz = nxy*NZ;
  //  int ntotz = QGmax*QGmax*U->nel*U->nz*nprocs;

  ntotz = max(ntotz, BC_ntotz);


  if (option("dealias")) ntotz = 3*ntotz/2;

  Pput = (double*) malloc(ntotz*sizeof(double));
  Pget = (double*) malloc(ntotz*sizeof(double));
  dzero(ntotz, Pput, 1);
  dzero(ntotz, Pget, 1);
}

void      packf(int nz, int ntot, double *from, double *data);
void    unpackf(int nz, int ntot, double *from, double *data);
#define packp FtPpackp
//void      packp(int nz, int   np, double *data, double *output);
void   PtFpackp(int nz, int   np, double *data, double *output);
void   FtPpackp(int nz, int   np, double *data, double *output);
#define unpackp PtFunpackp
//void    unpackp(int nz, int   np, double *data, double *output);
void PtFunpackp(int nz, int   np, double *data, double *output);
void FtPunpackp(int nz, int   np, double *data, double *output);
void   Exchange(int npts, double *a,    double *at);

void Fourier_List::FFT(Element_List *EL, Nek_Trans_Type ntt){
  if(nz>1){
    double *base_from;
    double *base_to;
    int     ntot;
    
    EL->Set_state(fhead->state);
    
    if(fhead->state == 'p'){
      base_from = base_h;
      base_to   = EL->base_h;
      ntot      = htot;
    }
    else{
      base_from = base_hj;
      base_to   = EL->base_hj;
      ntot      = hjtot;
    }

    int nprocs = option("NPROCS");
    int    nxy = (ntot+nprocs-1)/nprocs;
    int  ntotp = nxy*nprocs;
    int  ntotz = nxy*nztot;
    
    dzero(ntotz, Pput, 1);
    dzero(ntotz, Pget, 1);

    packf   ( nz,   ntot,  base_from, Pput);
    Exchange(ntotz, Pput,  Pget);
    unpackp ( nz,  ntotp,  Pget,      Pput);
    
#ifdef OLDFFTS
    int dir = (ntt==P_to_F) ? -1:1;
    for(int i = 0; i < nxy; ++i)
      realft(nztot/2, Pput+i*nztot, dir);
#else
    int mskip = nztot/2;
    if (ntt==F_to_P) {
#ifdef TESTFORZERO
      for(int i = 0; i < nxy; ++i)
	if (Pput[i*nztot+1] != 0.0) printf ("%d: FFT: %d - %g\n", option("PROCID"), i, Pput[i*nztot+1]);
#endif
      rfftw(rplan_inv, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
    } else {
      rfftw(rplan, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
    }
#endif
    packp   ( nz,   ntotp,  Pput, Pget);
    Exchange(ntotz, Pget,  Pput);
    unpackf ( nz,    ntot,  Pput, base_to);
  }
}

void Fourier_List::FFT(int ntot, double *base_from, double *base_to, Nek_Trans_Type ntt)
{
  // On input:
  //   base_from[this->nz*ntot] contains data to be transformed
  //
  // On exit:
  //   base_to[this->nz*ntot] contains data that has been transformed
  //
  // Note:
  //   base_from and base_to may or may not point to the same memory location
  //   if the same meory, then result will overwrite original data on exit
  //

  if(nz>1){
    //double *base_from;
    //double *base_to;
    //int     ntot;
#if 0    
    EL->Set_state(fhead->state);
    
    if(fhead->state == 'p'){
      base_from = base_h;
      base_to   = EL->base_h;
      ntot      = htot;
    }
    else{
      base_from = base_hj;
      base_to   = EL->base_hj;
      ntot      = hjtot;
    }
#endif

    int nprocs = option("NPROCS");
    int    nxy = (ntot+nprocs-1)/nprocs;
    int  ntotp = nxy*nprocs;
    int  ntotz = nxy*nztot;
    
    dzero(ntotz, Pput, 1);
    dzero(ntotz, Pget, 1);

    packf   ( nz,   ntot,  base_from, Pput);
    Exchange(ntotz, Pput,  Pget);
    unpackp ( nz,  ntotp,  Pget,      Pput);
    
#ifdef OLDFFTS
    int dir = (ntt==P_to_F) ? -1:1;
    for(int i = 0; i < nxy; ++i)
      realft(nztot/2, Pput+i*nztot, dir);
#else
    int mskip = nztot/2;
    if (ntt==F_to_P) {
#ifdef TESTFORZERO
      for(int i = 0; i < nxy; ++i)
	if (Pput[i*nztot+1] != 0.0) 
	  printf ("%d: FFT: %d - %g\n", 
		  mynode_1DZ() /*option("PROCID")*/, i, Pput[i*nztot+1]);
#endif
      rfftw(rplan_inv, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
    } else {
      rfftw(rplan, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
    }
#endif
    packp   ( nz,   ntotp,  Pput, Pget);
    Exchange(ntotz, Pget,  Pput);
    unpackf ( nz,    ntot,  Pput, base_to);
  }
}


void Fourier_List::FFT32(Element_List *EL, Nek_Trans_Type ntt){
  if(nz>1){
    double *base_from;
    double *base_to;
    int     ntot;
    
    EL->Set_state(fhead->state);
    
    if(fhead->state == 'p'){
      base_from = base_h;
      base_to   = EL->base_h;
      ntot      = htot;
    }
    else{
      base_from = base_hj;
      base_to   = EL->base_hj;
      ntot      = hjtot;
    }

    int nprocs = option("NPROCS");
    int    nxy = (ntot+nprocs-1)/nprocs;
    int  ntotp = nxy*nprocs;
    int  ntotz = nxy*nztot;
    
    dzero(3*ntotz/2, Pput, 1);
    dzero(3*ntotz/2, Pget, 1);

    if (ntt==F_to_P32) {
      packf   ( nz,   ntot,  base_from, Pput);
      Exchange(ntotz, Pput,  Pget);
      FtPunpackp ( nz,  ntotp,  Pget,      Pput);
    
      nz    = 3*nz/2;
      nztot = 3*nztot/2;
      ntotz = 3*ntotz/2;

#ifdef OLDFFTS
      for(int i = 0; i < nxy; ++i)
	realft(nztot/2, Pput+i*nztot, 1);
#else
#ifdef TESTFORZERO
      for(int i = 0; i < nxy; ++i)
	if (Pput[i*nztot+1] != 0.0) printf ("%d: FFT32: %d - %g\n", option("PROCID"), i, Pput[i*nztot+1]);
#endif
      int mskip = nztot/2;
      rfftw(rplan_inv32, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
#endif
      FtPpackp   ( nz,   ntotp,  Pput, Pget);
      Exchange(ntotz, Pget,  Pput);
      unpackf ( nz,    ntot,  Pput, base_to);
      
      nz = option("NZ");
      nztot = option("NZTOT");
    } else {
      nz    = 3*nz/2;
      nztot = 3*nztot/2;
      ntotz = 3*ntotz/2;

      packf   ( nz,   ntot,  base_from, Pput);
      Exchange(ntotz, Pput,  Pget);
      PtFunpackp ( nz,  ntotp,  Pget,      Pput);
    
#ifdef OLDFFTS
      for(int i = 0; i < nxy; ++i)
	realft(nztot/2, Pput+i*nztot, -1);
#else
      int mskip = nztot/2;
      rfftw(rplan32, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
#endif
      nz = option("NZ");
      nztot = option("NZTOT");
      ntotz = nxy*nztot;
      PtFpackp   ( nz,   ntotp,  Pput, Pget);
      Exchange(ntotz, Pget,  Pput);
      unpackf ( nz,    ntot,  Pput, base_to);
    }
  }
}

void Fourier_List::FFT32(int ntot, double *base_from, double *base_to, Nek_Trans_Type ntt)
{
  // On input:
  //   base_from[this->nz(3/2)*ntot] contains data to be transformed
  //
  // On exit:
  //   base_to[this->nz(3/2)*ntot] contains data that has been transformed
  //
  // Note:
  //   base_from and base_to may or may not point to the same memory location
  //   if the same meory, then result will overwrite original data on exit
  //

  if(nz>1){
    //double *base_from;
    //double *base_to;
    //int     ntot;
#if 0    
    EL->Set_state(fhead->state);
    
    if(fhead->state == 'p'){
      base_from = base_h;
      base_to   = EL->base_h;
      ntot      = htot;
    }
    else{
      base_from = base_hj;
      base_to   = EL->base_hj;
      ntot      = hjtot;
    }
#endif

    int nprocs = option("NPROCS");
    int    nxy = (ntot+nprocs-1)/nprocs;
    int  ntotp = nxy*nprocs;
    int  ntotz = nxy*nztot;
    
    dzero(3*ntotz/2, Pput, 1);
    dzero(3*ntotz/2, Pget, 1);

    if (ntt==F_to_P32) {
      packf   ( nz,   ntot,  base_from, Pput);
      Exchange(ntotz, Pput,  Pget);
      FtPunpackp ( nz,  ntotp,  Pget,      Pput);
    
      nz    = 3*nz/2;
      nztot = 3*nztot/2;
      ntotz = 3*ntotz/2;

#ifdef OLDFFTS
      for(int i = 0; i < nxy; ++i)
	realft(nztot/2, Pput+i*nztot, 1);
#else
#ifdef TESTFORZERO
      for(int i = 0; i < nxy; ++i)
	if (Pput[i*nztot+1] != 0.0) 
	  printf ("%d: FFT32: %d - %g\n", 
		  mynode_1DZ() /*option("PROCID")*/, i, Pput[i*nztot+1]);
#endif
      int mskip = nztot/2;
      rfftw(rplan_inv32, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
#endif
      FtPpackp   ( nz,   ntotp,  Pput, Pget);
      Exchange(ntotz, Pget,  Pput);
      unpackf ( nz,    ntot,  Pput, base_to);
      
      nz = option("NZ");
      nztot = option("NZTOT");
    } else {
      nz    = 3*nz/2;
      nztot = 3*nztot/2;
      ntotz = 3*ntotz/2;

      packf   ( nz,   ntot,  base_from, Pput);
      Exchange(ntotz, Pput,  Pget);
      PtFunpackp ( nz,  ntotp,  Pget,      Pput);
    
#ifdef OLDFFTS
      for(int i = 0; i < nxy; ++i)
	realft(nztot/2, Pput+i*nztot, -1);
#else
      int mskip = nztot/2;
      rfftw(rplan32, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
#endif
      nz = option("NZ");
      nztot = option("NZTOT");
      ntotz = nxy*nztot;
      PtFpackp   ( nz,   ntotp,  Pput, Pget);
      Exchange(ntotz, Pget,  Pput);
      unpackf ( nz,    ntot,  Pput, base_to);
    }
  }
}



void Fourier_List::Grad_z(Element_List *EL){
  if(fhead->state == 'p')
    Grad_z_d(EL->base_h);
  else
    Grad_z_d(EL->base_hj);  
}


void Fourier_List::Grad_z_d(double *EL){
  double *Ur;
  double *dUdz_r;
  int     ntot;
  int     i = 0;
  double *tmp, bet;
  
  if(fhead->state == 'p'){
    ntot   = htot;
    Ur     = base_h;     
  }
  else{
    ntot   = hjtot;
    Ur     = base_hj;    
  }
  
  dUdz_r = EL;

#ifdef FLOK
  if(nz == 2)
    {
      /* compute: dU/dz = i beta k ( Ur + i Ui) = beta k ( -Ui + i Ur ) */
      tmp = dvector(0,ntot-1);  // use tmp so that input and output can be same
      bet     = Beta(0);    dsmul(ntot,bet,Ur,1,tmp,1);
      bet     = -bet;       dsmul(ntot,bet,Ur+ntot,1,dUdz_r,1);
      dcopy(ntot,tmp,1,dUdz_r+ntot,1);
      free(tmp);
    }
  else
    {
      /* compute: dU/dz =  beta  Ur */
      bet       = Beta(0);
      dsmul(ntot,bet,Ur,1,dUdz_r,1);
    }
#else
  if(!parid(i)){
    dzero(2*ntot,dUdz_r,1);
    i += 2;
    if(i >= nz) return;
    dUdz_r += 2*ntot;
    Ur     += 2*ntot;
  }

  /* compute: dU/dz = i beta k ( Ur + i Ui) = beta k ( -Ui + i Ur ) */
  tmp = dvector(0,ntot-1);
  for(; i < nz; i += 2){
    bet     = Beta(i);    dsmul(ntot,bet,Ur,1,tmp,1);
    bet     = -bet;       dsmul(ntot,bet,Ur+ntot,1,dUdz_r,1);
    dcopy(ntot,tmp,1,dUdz_r+ntot,1);
    Ur     += 2*ntot;  
    dUdz_r += 2*ntot;  
  }
  free(tmp);
#endif
}

void Fourier_List::Grad_z_h(int nz, int ntot, double *HL, double *EL)
{
  // On input:
  //   HL[nz*ntot] contains velocity data in (*,F)-space
  //   nz : number of planes in z in HL
  //   ntot : dimension within each z plane
  //   EL[nz*ntot] to contain output results
  //
  // On exit:
  //   EL[nz*ntot] contains du/dz in (*,F)-space
  //   HL unchanged
  //

  double *Ur;
  double *dUdz_r;
  //int     ntot;
  int     i = 0;
  double *tmp, bet;

#if 0  
  if(fhead->state == 'p'){
    ntot   = htot;
    Ur     = base_h;     
  }
  else{
    ntot   = hjtot;
    Ur     = base_hj;    
  }
#endif

  Ur = HL;

  dUdz_r = EL;

  if(!parid(i)){
    dzero(2*ntot,dUdz_r,1);
    i += 2;
    if(i >= nz) return;
    dUdz_r += 2*ntot;
    Ur     += 2*ntot;
  }

  tmp = dvector(0,ntot-1);

  /* compute: dU/dz = i beta k ( Ur + i Ui) = beta k ( -Ui + i Ur ) */
  for(; i < nz; i += 2){
    bet     = Beta(i);    dsmul(ntot,bet,Ur,1,tmp,1);
    bet     = -bet;       dsmul(ntot,bet,Ur+ntot,1,dUdz_r,1);
    dcopy(ntot,tmp,1,dUdz_r+ntot,1);
    Ur     += 2*ntot;  
    dUdz_r += 2*ntot;  
  }

  free(tmp);
}

#ifdef MAP
void Grad_zz(Element_List *ELin, Element_List *ELout){
  if(ELin->fhead->state == 'p')
    Grad_zz_d(ELin, ELout->base_h);
  else
    Grad_zz_d(ELin, ELout->base_hj);
}
  
void Grad_zz_d(Element_List *ELin, double *ELout){
  double *Ur;
  double *dUdz_r;
  int     ntot;
  int     i = 0;
  double  bet2;
  
  if(ELin->fhead->state == 'p'){
    ntot   = ELin->htot;
    Ur     = ELin->base_h;     
  }
  else{
    ntot   = ELin->hjtot;
    Ur     = ELin->base_hj;    
  }
  
  dUdz_r = ELout;

  if(!parid(i)){
    dzero(2*ntot,dUdz_r,1);
    i += 2;
    if(i >= ELin->nz) return;
    dUdz_r += 2*ntot;
    Ur     += 2*ntot;
  }

  /* compute: d^2U/dz^2 = - (beta k)^2 ( Ur + i Ui) */
  for(; i < ELin->nz; i += 2){
    bet2    = -Beta2(i); dsmul(2*ntot,bet2,Ur,1,dUdz_r,1);
    Ur     += 2*ntot;  
    dUdz_r += 2*ntot;  
  }
}
#endif // ifdef MAP

/* pack from elements into local groups of nz */
/* direction of data flow ---->                 */
void packf(int nz, int ntot, double *from, double *data){
  //  dzero(nz*ntot, data, 1);
  for(int i = 0; i < nz; ++i)
    dcopy(ntot,from+i*ntot,1,data+i,nz);
}

/* unpack from elements into local groups of nz */
/* direction of data flow ---->                 */
void unpackf(int nz, int ntot, double *from, double *data){
  for(int i = 0; i < nz; ++i)
    dcopy(ntot,from+i, nz, data+i*ntot,1);
}

/* pack from groups of nztot into groups of nz */
void FtPpackp (int nz, int np, double *data, double *output){
  register int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  for(i = 0;i < nxy; ++i)
    for(j = 0; j < nprocs; ++j) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++) 
	output[i*nz+j*block+SLn] = data[i*nztot+j*nz+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nztot+j*nz, *SL_B = output+i*nz+j*block; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nztot+j*nz,output+i*nz+j*block);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nztot+j*nz,1,output+i*nz+j*block,1);
#endif // ifdef SL
    }
  return;
}

void PtFpackp (int nz, int np, double *data, double *output){
  register int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  if (option("dealias")) nztot = 3*nztot/2;
  for(i = 0;i < nxy; ++i)
    for(j = 0; j < nprocs; ++j) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++)
	output[i*nz+j*block+SLn] = data[i*nztot+j*nz+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nztot+j*nz, *SL_B = output+i*nz+j*block; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nztot+j*nz,output+i*nz+j*block);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nztot+j*nz,1,output+i*nz+j*block,1);
#endif // ifdef SL
    }
  return;
}

/* unpack vector of local groups of nz and put into groups of nztot */
void FtPunpackp (int nz, int np, double *data, double *output){
  int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  if (option("dealias")) {
    nztot = 3*nztot/2;
    dzero (nxy*nztot, output, 1);
  }
  for(j = 0; j < nprocs; ++j)
    for(i = 0;i < nxy; ++i) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++)
	output[i*nztot+j*nz+SLn] = data[i*nz+j*block+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nz+j*block, *SL_B = output+i*nztot+j*nz; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nz+j*block,output+i*nztot+j*nz);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nz+j*block,1,output+i*nztot+j*nz,1);
#endif // ifdef SL
    }
  return;
}

void PtFunpackp (int nz, int np, double *data, double *output){
  register int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  for(j = 0; j < nprocs; ++j)
    for(i = 0;i < nxy; ++i) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++)
	output[i*nztot+j*nz+SLn] = data[i*nz+j*block+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nz+j*block, *SL_B = output+i*nztot+j*nz; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nz+j*block,output+i*nztot+j*nz);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nz+j*block,1,output+i*nztot+j*nz,1);
#endif // ifdef SL
    }
  return;
}

void Exchange(int npts, double *a, double *at){
  int nprocs     = option("NPROCS");
  int block_size = npts / nprocs;
  int msg_size   = (int) block_size * sizeof(double);

  if (npts % nprocs) 
    fprintf(stderr,"Exchange: npts should be a multiple of nprocs!!!\n");

  /* copy the local block directly */
#ifdef PARALLEL
  MPI_Alltoall(a, msg_size, MPI_BYTE, at, msg_size, MPI_BYTE, MPI_COMM_WORLD);
#else
  dcopy(block_size,a,1,at,1);
#endif

  return;
}

#define DESCRIP 25
static int gettypes    (char *t, char *s),
           checkfmt    (char *format);

static char *fourier_hdr_fmt[] = { 
  "%-25s "            "Session\n",
  "%-25s "            "Created\n",
  "%-5c Fourier Hybrid      " "State 'p' = physical, 't' transformed\n",
  "%-5d %-5d %-5d %-5d   "    "Number of Elements; Dim of run; Lmax; NZ\n",
  "%-25d "            "Step\n",
  "%-25.6g "          "Time\n",
  "%-25.6g "          "Time step\n",
  "%-11.6g %-11.6g   " "Kinvis; LZ\n",
  "%-25s "            "Fields Written\n",
  "%-25s "            "Format\n"
  };
#define READLINE(fmt,arg) \
  if (fscanf(fp,fmt,arg)==EOF) return -1; fgets(buf,BUFSIZ,fp)

int readHeaderF(FILE* fp, Field *f){
  register int  nfields;
  char     buf[BUFSIZ];
  int      trip=0;
  if (fscanf (fp, "%s", buf) == EOF) return 0;            /* session name */
  f->name    = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               /* creation date */
  f->created = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  
  READLINE ("%c"  , &f->state);                  /* simulation parameters */

  /* check to see if it is standard 2d file */
  if(strstr(buf,"Fourier")){
    trip = 1;
    fscanf(fp,"%d%d%d%d",&f->nel,&f->dim,&f->lmax,&f->nz);fgets(buf,BUFSIZ,fp);
  }
  else{
    fscanf(fp,"%d%d%d",&f->nel,&f->dim,&f->lmax);fgets(buf,BUFSIZ,fp);
    f->nz = 1;
  }
  
  READLINE ("%d"  , &f->step);                  
  READLINE ("%lf" , &f->time);
  READLINE ("%lf" , &f->time_step);
  if(trip){
    fscanf(fp,"%lf%lf " ,&f->kinvis,&f->lz);
    fgets(buf,BUFSIZ,fp);
  }
  else{
    READLINE ("%lf" , &f->kinvis);
    f->lz = 1.0;
  }
  
  nfields = gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = (char *)strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("Warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);
  rewind(fp);
  return nfields;
}



/* ---------------------------------------------------------------------- *
 * readFieldF() -- Load a field file                                       *
 *                                                                        *
 * This function loads a field file into a Field structure.  The format   *
 * is a simple array, stored in 64-bit format.                            *
 *                                                                        *
 * Return value: number of variables loaded, 0 for EOF, -1 for error.     *
 * ---------------------------------------------------------------------- */


int readFieldF (FILE *fp, Field *f, Element_List *E)
{
  register int i, n, ntot, nfields;
  int      trip = 0;
  char     buf[BUFSIZ];
#ifdef _CRAY
  short *ssize;
#endif
  if (fscanf (fp, "%s", buf) == EOF) return 0;            /* session name */
  f->name    = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               /* creation date */
  f->created = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  
  READLINE ("%c"  , &f->state);                  /* simulation parameters */

  /* check to see if it is standard 2d file */
  if(strstr(buf,"Fourier")){
    trip = 1;
    fscanf(fp,"%d%d%d%d",&f->nel,&f->dim,&f->lmax,&f->nz);fgets(buf,BUFSIZ,fp);
  }
  else{
    fscanf(fp,"%d%d%d",&f->nel,&f->dim,&f->lmax);fgets(buf,BUFSIZ,fp);
    f->nz = 1;
  }

  READLINE ("%d"  , &f->step);                  
  READLINE ("%lf" , &f->time);
  READLINE ("%lf" , &f->time_step);
  if(trip){
    fscanf(fp,"%lf%lf " ,&f->kinvis,&f->lz);
    fgets(buf,BUFSIZ,fp);
  }
  else{
    READLINE ("%lf" , &f->kinvis);
    f->lz = 1.0;
  }
  nfields = gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = (char *)strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("Warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);

  /* allocate memory and load the data */
  
  int cnt = count_facets(E->fhead);
#ifndef FULL_FIELD
  int nltot, nrtot, skip1, skip2,
      proc   = option("PROCID"),
      NZ     = option("NZ"),
      NRZ    = min(NZ, f->nz - proc*NZ);
#endif

  switch (tolower(*f->format)) {
  case 'a':
    if(f->state == 't'){
      f->size = ivector(0,cnt-1);
      for(i = 0; i < cnt ; ++i) fscanf(fp,"%d",f->size+i);
    }
    else {
      fputs("Error: reading field file in physical space is "
	    "not implemented \n",stderr);
      exit(-1);
    }
    
    ntot = f->nz*data_len(f->nel,f->size, E->fhead);
#ifndef FULL_FIELD
    nltot = NZ*data_len(f->nel,f->size,E->fhead);
    nrtot = NRZ*data_len(f->nel,f->size,E->fhead);
    skip1 = nltot * proc;
    if (NRZ < 1) 
      return nfields;

    for (i = 0; i < skip1; i++)
      fgets(buf,BUFSIZ,fp);
#endif

    for(i = 0; i < nfields; ++i)
      if(!(f->data [i]))
	f->data[i] = dvector (0, nrtot-1);

    for (i = 0; i < nrtot; i++) {
      for (n = 0; n < nfields; n++)
	if (fscanf (fp, "%lf", f->data[n] + i) != 1) {
	  fprintf(stderr,"Error: reading field %c, point %d\n",f->type[n],i+1);
	  exit(-1);
	}
    }
    fgets (buf, BUFSIZ, fp);
    break;
    
  case 'b':
    if(f->state == 't'){
      f->size = ivector(0,cnt-1);
#ifdef _CRAY
      // int on Crays is 64b, short 32b - use short for portable fieldfiles
      ssize = (short *) calloc(cnt, sizeof(short));
      fread(ssize, sizeof(short), (size_t) cnt, fp);
      for (n = 0; n < cnt; n++)
	f->size[n] = (int) ssize[n];
      free (ssize);
#else
      fread(f->size, sizeof(int), (size_t) cnt,fp);
#endif
    }
    else {
      fputs("Error: reading field file in physical space is "
	    "not implemented \n",stderr);
      exit(-1);
    }

    ntot = f->nz*data_len(f->nel,f->size,E->fhead);

#ifndef FULL_FIELD
    nltot = NZ*data_len(f->nel,f->size,E->fhead);
    nrtot = NRZ*data_len(f->nel,f->size,E->fhead);
    skip1 = nltot * proc;
    skip2 = ntot - nrtot;
    if (NRZ < 1) 
      return nfields;

    fseek (fp, (long) (skip1*sizeof(double)), SEEK_CUR);
#endif

    for(i = 0; i < nfields; ++i)
      if(!(f->data [i]))
	f->data[i] = dvector (0, nrtot-1);
    
    for (n = 0; n < nfields; n++) {
      if (fread (f->data[n], sizeof(double), (size_t) nrtot, fp) != nrtot) {
	fprintf (stderr, "error reading field %c", f->type[n]);
	exit(-1);
      }
#ifndef FULL_FIELD
      fseek (fp, (long) (skip2*sizeof(double)), SEEK_CUR);
#endif
    }
    break;

  default:
    fprintf (stderr, "unknown field format -- %s\n", f->format);
    exit(-1);
    break;
  }

  return nfields;
}

/* ---------------------------------------------------------------------- *
 * writeField() -- Write a field file                                     *
 *                                                                        *
 * This function writes a field file from a Field structure.  The format  *
 * is a simple array, stored in 64-bit format.                            *
 *                                                                        *
 * Return value: number of variables written, -1 for error.               *
 * ---------------------------------------------------------------------- */

int writeFieldF (FILE *fp, Field *f, Element *E)
{
  int  nfields, ntot;
  char buf[BUFSIZ];
  register int i, n;
#ifdef _CRAY
  short *ssize;
#endif

  /* Write the header */
  fprintf (fp, fourier_hdr_fmt[0],  f->name);
  fprintf (fp, fourier_hdr_fmt[1],  f->created);
  fprintf (fp, fourier_hdr_fmt[2],  f->state);
  fprintf (fp, fourier_hdr_fmt[3],  f->nel, f->dim, f->lmax, f->nz);

  fprintf (fp, fourier_hdr_fmt[4],  f->step);
  fprintf (fp, fourier_hdr_fmt[5],  f->time);
  fprintf (fp, fourier_hdr_fmt[6],  f->time_step);

  fprintf (fp, fourier_hdr_fmt[7],  f->kinvis, f->lz);
  fprintf (fp, fourier_hdr_fmt[8],  f->type);
  fprintf (fp, fourier_hdr_fmt[9],  f->format);

  /* Write the field files  */
  if(f->state == 't'){
    ntot = f->nz*data_len(f->nel,f->size,E);
  }
  else{
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
      
  nfields = (int) strlen (f->type);

  int cnt = count_facets(E);
  switch (tolower(*f->format)) {
  case 'a':
    for(i = 0; i < cnt; ++i)
      fprintf(fp,"%d ",f->size[i]);
    fputc('\n',fp);
    for (i = 0; i < ntot; i++) {
      for (n = 0; n < nfields; n++)
	fprintf (fp, "%#16.10g ", f->data[n][i]);
      fputc ('\n', fp);
    }
    break;
    
  case 'b':
#ifdef _CRAY
    // int on Crays is 64b, short 32b - use short for portable fieldfiles
    ssize = (short *) calloc(cnt, sizeof(short));
    for (n = 0; n < cnt; n++)
      ssize[n] = (short) f->size[n];
    fwrite(ssize, sizeof(short),  cnt, fp);
    free (ssize);
#else
    fwrite(f->size, sizeof(int),  cnt, fp);
#endif
    for (n = 0; n < nfields; n++)
      if (fwrite (f->data[n], sizeof(double), ntot, fp) != ntot) {
	fprintf  (stderr, "error writing field %c", f->type[n]);
	exit(-1);
      }
    break;
    
  default:
    sprintf  (buf, "unknown format -- %s", f->format);
    fprintf(stderr,buf);
    exit(-1);
    break;
  }

  fflush (fp);
  return nfields;
}

#ifndef BTYPE
#if defined(i860) || defined (__alpha) || defined (__WIN32__) || (defined(linux) && defined(i386))
#define BTYPE "ieee_little_endian"
#endif
#
#if defined(_CRAY) && !defined (_CRAYMPP)
#define BTYPE "cray"
#endif /* ........... Cray Y-MP ........... */
#
#ifndef BTYPE
#define BTYPE "ieee_big_endian"
#endif /* default case in the absence of any other TYPE */
#endif /* ifndef TYPE */

void WritefieldF(FILE *fp, char *name, int step, double t,
		 int nfields, Element_List *U[]){
  register int n,i,j;
  time_t   tp;
  char     buf[128],state;
  Field    f;
  int      ntot;
  double   **u;
  int      nz     = U[0]->nz;
  int      nztot  = U[0]->nztot;
  int      procid = option("PROCID"), nprocs = option("NPROCS");

  if(nfields > MAXFIELDS)
    error_msg(too many fields -- must be less than MAXFIELDS);

  state = U[0]->fhead->state;
  /* check to see if all fields in same state */
  for(i = 1; i < nfields; ++i)
    if(state != U[i]->fhead->state)
      {error_msg(Writefield--Fields  in different space);}
      
  /* Get the current date and time from the system */
  
  tp = time((time_t*)NULL);
  strftime(buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  
  /* Set up the field structure */
  f.name      =  name;
  f.created   =  buf;
  f.state     =  state;
  f.dim       =  DIM;
  f.nel       =  U[0]->nel;
  f.lmax      =  LGmax;
  f.step      =  step;
  f.time      =  t;
  f.nz        =  nztot;
  f.lz        =  dparam("LZ");
  f.time_step =  dparam("DT");
  f.kinvis    =  (dparam("KINVIS-ORIG"))?
    dparam("KINVIS-ORIG") : dparam("KINVIS");   /* VARV */
  f.format    = (char *)((option("binary")) ? "binary-"BTYPE : "ascii");

  /* Generate the list of fields and assign data pointers */
  
  memset (f.type, '\0', MAXFIELDS);
  memset (f.data, '\0', MAXFIELDS * sizeof( double *));
  
  /* copyfield into matrix */
  if(state == 't'){
    ROOT {
      int cnt = count_facets(U[0]->fhead);
      f.size = ivector(0,cnt-1);
      for(i = 0,n=0; i < f.nel; ++i){
	for(j = 0; j < U[0]->flist[i]->Nedges; ++j,++n)
	  f.size[n] = U[0]->flist[i]->edge[j].l;
	
	for(j = 0; j < U[0]->flist[i]->Nfaces; ++j,++n)
	  f.size[n] = U[0]->flist[i]->face[j].l;
      }
    }
    ntot = U[0]->hjtot*nz;

    if(!procid)
      u = dmatrix(0,nfields-1,0,nprocs*ntot-1);    
    
    for(n = 0; n < nfields; ++n){
#ifdef PARALLEL
      if (procid) {
	MPI_Send (U[n]->base_hj,ntot,
		  MPI_DOUBLE, 0, W_MSG+procid, MPI_COMM_WORLD);
      }
      else {	
	dcopy(ntot, U[n]->base_hj, 1, u[n], 1);
	MPI_Status status;   
	for (i = 1; i < nprocs; ++i){
	  MPI_Recv (u[n]+i*ntot, ntot, MPI_DOUBLE, 
		    i, W_MSG+i, MPI_COMM_WORLD, &status);
	}
      }
#else
      dcopy(ntot, U[n]->base_hj, 1, u[n], 1);
#endif
    }
  }
  else{
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }

  ROOT {
    for(n = 0; n < nfields; ++n){
      f.type [n] = U[n]->fhead->type;
      f.data [n] = u[n];
    }
    writeFieldF(fp,&f,U[0]->fhead);
    free_dmatrix(u,0,0);
    free(f.size);
  }

  return;
}

#if 0
/*--------------------------------------------------------------------------*
 * This is a function that copies the field form the field structure 'f'    *
 * at the position given by pos to the Element 'U'. If in transformed space *
 * it allows the field to be copied at a different order. To do this it     *
 * either zeros the higher modes if U is higher than f else it ignores the  *
 * higher modes if f is higher than U                                       *
 *--------------------------------------------------------------------------*/
void copyfieldF(Field *f, int pos, Element_List *U){
  if(f->state == 'p'){ /* restart from fixed m physical field */
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
  else{ /* restart from transformed field at any order */
    register int n,i;
    int      cl,l;
    int      *size  = f->size;
    double   *data  = f->data[pos];
    
    register int j,k;
    int  nz = f->nz;
    Element *E;
    
    for(k = 0; k < nz; ++k){
      size  = f->size;
      for(j = 0; j < f->nel; ++j){
	E = U->flevels[k]->flist[j];
	
	/* copy vertices */
	for(i = 0; i < E->Nverts; ++i){
	  E->vert[i].hj[0] = *data;
	  ++data;
	}
	
	/* copy edges */
	for (i = 0; i < E->Nedges; ++i){
	  if(l = E->edge[i].l){
	    dzero(l,E->edge[i].hj,1);
	    cl = min(l,*size);
	    dcopy(cl,data,1,E->edge[i].hj,1);
	  }
	  data += *size;
	  ++size;
	}
	
	/* copy faces */
	for(n = 0; n < E->Nfaces; ++n){
	  if(E->identify()==Nek_Tri){
	    if(l = E->face[n].l){
	      dzero(l*(l+1)/2,*E->face[n].hj,1);
	      cl = min(l,*size);
	      for(i = 0; i < cl; ++i){
		dcopy(cl-i,data,1,E->face[n].hj[i],1);
		data += *size-i;
	      }
	      for(i = cl; i < *size; ++i)
		data += *size-i;
	    }
	    else
	      for(i = 0; i < *size; ++i)
		data += *size-i;
	    ++size;
	  }
	  else{
	    if(l = E->face[n].l){
	      dzero(l*l,*E->face[n].hj,1);
	      cl = min(l,*size);
	      for(i = 0; i < cl; ++i){
		dcopy(cl,data,1,E->face[n].hj[i],1);
		data += *size;
		}
	      for(i = cl; i < *size; ++i)
		data += *size;
	    }
	    else
	      for(i = 0; i < *size; ++i)
		data += *size;
	    ++size;
	  }
	}
	E->state = 't';
      }
    }
  }
}
#endif

/*--------------------------------------------------------------------------*
 * This is a function that copies the field form the field structure 'f'    *
 * at the position given by pos to the Element 'U'. If in transformed space *
 * it allows the field to be copied at a different order. To do this it     *
 * either zeros the higher modes if U is higher than f else it ignores the  *
 * higher modes if f is higher than U                                       *
 *--------------------------------------------------------------------------*/
void copyfieldF(Field *f, int pos, Element_List *U){
  if(f->state == 'p'){ /* restart from fixed m physical field */
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
  else{ /* restart from transformed field at any order */
    register int n,i;
    int      cl,l;
    int      *size  = f->size;
    double   *data  = f->data[pos];
    
    register int j,k;
    int  nz = min(f->nz,option("NZTOT"));
    Element *E;
    
    for(k = 0; k < nz; ++k){
      size  = f->size;
      for(j = 0; j < f->nel; ++j){
	E = U->flevels[k]->flist[j];
	
	/* copy vertices */
	for(i = 0; i < E->Nverts; ++i){
	  E->vert[i].hj[0] = *data;
	  ++data;
	}
	
	/* copy edges */
	for (i = 0; i < E->Nedges; ++i){
	  if(l = E->edge[i].l){
	    dzero(l,E->edge[i].hj,1);
	    cl = min(l,*size);
	    dcopy(cl,data,1,E->edge[i].hj,1);
	  }
	  data += *size;
	  ++size;
	}
	
	/* copy faces */
	for(n = 0; n < E->Nfaces; ++n){
	  if(E->identify()==Nek_Tri){
	    if(l = E->face[n].l){
	      dzero(l*(l+1)/2,*E->face[n].hj,1);
	      cl = min(l,*size);
	      for(i = 0; i < cl; ++i){
		dcopy(cl-i,data,1,E->face[n].hj[i],1);
		data += *size-i;
	      }
	      for(i = cl; i < *size; ++i)
		data += *size-i;
	    }
	    else
	      for(i = 0; i < *size; ++i)
		data += *size-i;
	    ++size;
	  }
	  else{
	    if(l = E->face[n].l){
	      dzero(l*l,*E->face[n].hj,1);
	      cl = min(l,*size);
	      for(i = 0; i < cl; ++i){
		dcopy(cl,data,1,E->face[n].hj[i],1);
		data += *size;
		}
	      for(i = cl; i < *size; ++i)
		data += *size;
	    }
	    else
	      for(i = 0; i < *size; ++i)
		data += *size;
	    ++size;
	  }
	}
	E->state = 't';
      }
    }
  }
}

static int gettypes (char *t, char *s)
{
  char    *p = strstr (s, "Fields");
  register int n = 0;
  if (!p){ fprintf(stderr,"invalid \"Fields\" string");exit(-1);}
  memset(t, '\0', MAXFIELDS);
  while (s != p && n < MAXFIELDS) {
    if (isalpha(*s)) t[n++] = *s;
    s++;
  }
  return n;
}


/* Check binary format compatibility */
#if defined(__sgi) || defined(cm5) || defined(_AIX) || defined (__hpux) || defined (__sparc) || defined (_CRAYMPP)
#define ieeeb
#define ieee
#endif

#if defined(i860) || defined (__alpha) || defined (__WIN32__)
#define ieeel
#endif

static int checkfmt (char *arch)
{
  char        **p;
  static char *fmtlist[] = {
#if defined(ieeeb)                  /* ... IEEE big endian machine ..... */
    "ieee_big_endian",
    "sgi", "iris4d", "SGI", "IRIX", /* kept for compatibility purposes   */
    "IRIX64",                       /* ........ Silicon Graphics ....... */
    "AIX",                          /* .......... IBM RS/6000 .......... */
    "cm5",                          /* ........ Connection Machine ..... */
#endif
#
#if defined(ieeel)                  /* ... IEEE little endian machine .. */
    "ieee_little_endian",
    "i860",                         /* ........... Intel i860 .......... */
#endif
#
#if defined(ieee)                   /* ...... Generic IEEE machine ..... */
    "ieee", "sim860",               /* kept for compatibility purposes   */
#endif                              /* same as IEEE big endian           */
#
#if defined(_CRAY) && !defined (_CRAYMPP) /* ...... Cray PVP ........... */
    "cray", "CRAY",                 /* kept for compatibility purposes   */
#endif                              /* no conflict with T3* as it        */
#                                   /* precedes this line                */
     0 };   /* a NULL string pointer to signal the end of the list */

  for (p = fmtlist; *p; p++)
    if (strncmp (arch, *p, strlen(*p)) == 0)
      return 0;

  return 1;
}

#undef ieee

void packBCf(Domain *Omega,int fields, double *data);
void unpackBCf(Domain *Omega,int fields, double *data);


void TransformBCs(Domain *Omega, Nek_Trans_Type ntt){
  int nprocs = option("NPROCS");
  int fields = 3;
#ifdef THERMO
  fields++;
#endif
  int nz     = option("NZ");
  int NZ     = nz*nprocs;
  int nel    = countbcs(Omega);
  int tot    = LGmax;/*Omega[0]->U->edge[0].l+2;*/
  int ntot   = tot*nel*fields;
  int nxy    = (ntot+nprocs-1)/nprocs;
  int ntotp  = nxy*nprocs;
  int ntotz  = nxy*NZ;

  dzero(ntotz, Pput, 1);
  dzero(ntotz, Pget, 1);

  packBCf (Omega,fields,Pput);
  Exchange(ntotz,Pput,Pget);
  unpackp ( nz, ntotp,Pget,Pput);

#ifdef OLDFFTS
  int i;
  int dir = (ntt==P_to_F) ? -1:1;
  for(i = 0; i < nxy; ++i)
    realft(NZ/2,Pput+i*NZ,dir);
#else
  int mskip = NZ/2;
  if (ntt == P_to_F) {
    rfftw(rplan, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
  } else
    rfftw(rplan_inv, nxy, (FFTW_COMPLEX *) Pput, 1, mskip, 0, 0, 0);
#endif

  packp   ( nz, ntotp,Pput,Pget);
  Exchange(ntotz,Pget,Pput);
  unpackBCf(Omega,fields,Pput);
  return;
}

void packBCf(Domain *Omega,int fields, double *data){
  int i, cnt;
  int nel    = countbcs(Omega);
  int nz     = option("NZ");
  int tot    = Omega->U->fhead->edge->l;
  int ntot   = (LGmax-2)*nel;
  int offset = fields*ntot*nz;
  Bndry *B;
  
  dzero(fields*LGmax*nel*nz, data, 1);
  for(i=0;i<nz;++i){
    for(cnt=0,B=Omega->Ubc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;
      dcopy(tot, B->bedge[0], 1, data+i+cnt*(LGmax-2)*nz,       nz);
    }

    for(cnt=0,B=Omega->Vbc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;	
      dcopy(tot, B->bedge[0], 1, data+i+cnt*(LGmax-2)*nz+ntot*nz  ,nz);
    }

    for(cnt=0,B=Omega->Wbc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;	
      dcopy(tot, B->bedge[0], 1, data+i+cnt*(LGmax-2)*nz+2*ntot*nz,nz);
    }
#ifdef THERMO
    for(cnt=0,B=Omega->Tbc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;	
      dcopy(tot, B->bedge[0], 1, data+i+cnt*(LGmax-2)*nz+3*ntot*nz,nz);
    }
#endif
    for(cnt=0, B=Omega->Ubc[i]; B; B=B->next,++cnt)
      dcopy(2, B->bvert, 1, data+i+offset+cnt*2*nz,         nz);
    for(cnt=0, B=Omega->Vbc[i]; B; B=B->next,++cnt)
      dcopy(2, B->bvert, 1, data+i+offset+cnt*2*nz+2*nel*nz,nz);
    for(cnt=0, B=Omega->Wbc[i]; B; B=B->next,++cnt)
      dcopy(2, B->bvert, 1, data+i+offset+cnt*2*nz+4*nel*nz,nz);
#ifdef THERMO
    for(cnt=0, B=Omega->Tbc[i]; B; B=B->next,++cnt)
      dcopy(2, B->bvert, 1, data+i+offset+cnt*2*nz+6*nel*nz,nz);
#endif
  }
  
}

void unpackBCf(Domain *Omega,int fields, double *data){
  int i, cnt;
  int nel    = countbcs(Omega);
  int nz     = option("NZ");
  int tot    = Omega->U->fhead->edge->l;
  int ntot   = (LGmax-2)*nel;
  int offset = fields*ntot*nz;
  Bndry *B;

  for(i=0;i<nz;++i){
    for(cnt=0,B=Omega->Ubc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;
      dcopy(tot,data+i+cnt*(LGmax-2)*nz          ,nz,B->bedge[0],1);
    }
    for(cnt=0,B=Omega->Vbc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;
      dcopy(tot,data+i+cnt*(LGmax-2)*nz+ntot*nz  ,nz,B->bedge[0],1);
    }
    for(cnt=0,B=Omega->Wbc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;
      dcopy(tot,data+i+cnt*(LGmax-2)*nz+2*ntot*nz,nz,B->bedge[0],1);
    }
#ifdef THERMO
    for(cnt=0,B=Omega->Tbc[i];B;B=B->next,++cnt){
      tot = B->elmt->edge[B->face].l;
      dcopy(tot,data+i+cnt*(LGmax-2)*nz+3*ntot*nz,nz,B->bedge[0],1);
    }
#endif
    
    for(cnt=0,B=Omega->Ubc[i];B;B=B->next,++cnt)
      dcopy(2,data+i+offset+cnt*2*nz,nz,B->bvert,1);
    for(cnt=0,B=Omega->Vbc[i];B;B=B->next,++cnt)
      dcopy(2,data+i+offset+cnt*2*nz+2*nel*nz,nz,B->bvert,1);
    for(cnt=0,B=Omega->Wbc[i];B;B=B->next,++cnt)
      dcopy(2,data+i+offset+cnt*2*nz+4*nel*nz,nz,B->bvert,1);
#ifdef THERMO
    for(cnt=0,B=Omega->Tbc[i];B;B=B->next,++cnt)
      dcopy(2,data+i+offset+cnt*2*nz+6*nel*nz,nz,B->bvert,1);
#endif
  }
}

int countbcs(Domain *Omega){
  int    i = 0;
  Bndry *B;
  for(B = Omega->Ubc[0]; B; B = B->next)
    ++i;
  return i;
}

#ifdef OLDFFTS
#ifdef _AIX
static double *factor = (double*)0;
static int    *bitrev =    (int*)0;
extern "C"
{
void    fftdf  (double *, int *, int *, int *, int *, int *, double *,
	    int *, int *, int *, int *);
}
#endif
#endif

void AtoA_FFT(Element_List *EL, double *base_to){
  if(EL->nz>1){
    double *base_from = EL->base_h;
    int     ntot      = EL->htot;
    int     nprocs    = option("NPROCS");
    int     nxy       = (ntot+nprocs-1)/nprocs;
    int     ntotp     = nxy*nprocs;
    int     nztot     = EL->nztot;
    int     ntotz     = nxy*nztot;
    
    if (option("dealias")) {
      dzero(3*ntotz/2, Pput, 1);
      dzero(3*ntotz/2, Pget, 1);
      nztot = 3*nztot/2;
    } else {
      dzero(ntotz, Pput, 1);
      dzero(ntotz, Pget, 1);
    }
    
    packf   ( EL->nz,   ntot,  base_from, Pput);
    Exchange( ntotz,    Pput,  Pget);
    FtPunpackp ( EL->nz,  ntotp,  Pget,      base_to);

#ifdef OLDFFTS
    int     dir       = 1;  // FFT -> physical
#ifndef _AIX    
    for(int i = 0; i < nxy; ++i)
      realft(nztot/2, base_to+i*nztot, 1);
#else
    double *Fu = base_to;
    int      n = nztot/2;
    int  nskip = 1;
    int  mskip = n;
    int      m = nxy;
    int  isign = dir;
    int   irev = 1;
    int   ierr = 1;
    int  ireal = 1;

    if (!factor) {                /* Check for initialization */
      factor = dvector(0, 6 * n);
      bitrev = ivector(0, 6 * n);
      fftdf_ (Fu, &n, &nskip, &m, &mskip, &isign, factor,
	      &irev, bitrev, &ierr, &ireal);
    }

    fftdf (Fu, &n, &nskip, &m, &mskip, &isign, factor,
	    &irev, bitrev, &ierr, &ireal);
#endif
#else
    int mskip = nztot/2;
    rfftw(rplan_inv32, nxy, (FFTW_COMPLEX *) base_to, 1, mskip, 0, 0, 0);
#endif
  }
}

void FFT_AtoA(Element_List *EL, double *base_from){
  if(EL->nz>1){
    double *base_to = EL->base_h;
    int     ntot    = EL->htot;
    int nprocs      = option("NPROCS");
    int    nxy      = (ntot+nprocs-1)/nprocs;
    int  ntotp      = nxy*nprocs;
    int   nztot     = EL->nztot;
    int  ntotz      = nxy*nztot;
    
    if (option("dealias")) nztot = 3*nztot/2;
#ifdef OLDFFTS
    int     dir     = -1;  // FFT -> fourier
#ifndef _AIX
    for(int i = 0; i < nxy; ++i)
      realft(nztot/2, base_from+i*nztot, dir);
#else
    double *Fu = base_from;
    int      n = nztot/2;
    int  nskip = 1;
    int  mskip = n;
    int      m = nxy;
    int  isign = dir;
    int   irev = 1;
    int   ierr = 1;
    int  ireal = 1;
    double scal   = 0.25 / n;

    fftdf  (Fu, &n, &nskip, &m, &mskip, &isign, factor,
	    &irev, bitrev, &ierr, &ireal);
    dscal (2*n*m, scal, Fu, 1);
#endif
#else
    int mskip = nztot/2;
    rfftw(rplan32, nxy, (FFTW_COMPLEX *) base_from, 1, mskip, 0, 0, 0);
#endif

    PtFpackp( EL->nz,   ntotp,  base_from, Pget);
    Exchange(  ntotz,    Pget,  Pput);
    unpackf ( EL->nz,    ntot,  Pput, base_to);
  }
}

void multifft(int nztot, int nxy, double *from, int dir){
#ifdef OLDFFTS
#ifndef _AIX
  for(int i = 0; i < nxy; ++i)
    realft(nztot/2, from+i*nztot, dir);
#else
  int      n = nztot/2;
  int  nskip = 1;
  int  mskip = n;
  int      m = nxy;
  int  isign = dir;
  int   irev = 1;
  int   ierr = 1;
  int  ireal = 1;
  double scal   = 0.25 / n;

  if (!factor) {                /* Check for initialization */
    factor = dvector(0, 6 * n);
    bitrev = ivector(0, 6 * n);
    fftdf_ (from, &n, &nskip, &m, &mskip, &isign, factor,
	    &irev, bitrev, &ierr, &ireal);
  }
  if(dir==1)  
    fftdf (from, &n, &nskip, &m, &mskip, &isign, factor,
	   &irev, bitrev, &ierr, &ireal);
  else{
    fftdf  (from, &n, &nskip, &m, &mskip, &isign, factor,
	    &irev, bitrev, &ierr, &ireal);
    dscal (2*n*m, scal, from, 1);
  }
#endif
#else
    int mskip = nztot/2;
    if (dir == -1) {
      rfftw(rplan, nxy, (FFTW_COMPLEX *) from, 1, mskip, 0, 0, 0);
    } else
      rfftw(rplan_inv, nxy, (FFTW_COMPLEX *) from, 1, mskip, 0, 0, 0);
#endif
}

#ifndef NEEDS_CLEANUP
double *Pput_new = (double*)0;
double *Pget_new = (double*)0;

double Global_Beta(int i){
  return  (i==1) ? 2.0*M_PI*(option("NZTOT")/2)/dparam("LZ"):
  2.0*M_PI*(i/2.0)/dparam("LZ");
}

void grad_z_new(int nztot, int  nxy, double *u, double *uz){
  double bet;
  int i;
  dzero(nxy, uz,   nztot);
  dzero(nxy, uz+1, nztot);
  
  for(i=2;i<nztot;i +=2){
    bet = Global_Beta(i);
    dsmul(nxy, bet, u+i,   nztot, uz+i+1, nztot);
    bet = -bet; 
    dsmul(nxy, bet, u+i+1, nztot, uz+i,   nztot);
  }
}

double **pstore = (double**)0;
static double *uz = (double*)0, *vz = (double*)0, *wz = (double*)0,
               *u = (double*)0, *v  = (double*)0, *w  = (double*)0,
	      *uF = (double*)0, *vF = (double*)0, *wF = (double*)0;

// assumes all elements have same number of d.o.f
void parallel_VdGradV(Domain *omega){

  Element_List *U  =  omega->U,  *V  =  omega->V,  *W  = omega->W,
              *Uf  =  omega->Uf, *Vf =  omega->Vf, *Wf = omega->Wf;

  Element     *eU, *eV, *eW;
  
  int nztot    = U->nztot;
  int nz       = U->nz;
  int procid   = option("PROCID");
  int nprocs   = option("NPROCS");  
  int nel_frac = (U->nel+nprocs-1)/nprocs;
  int Nm       = U->fhead->Nmodes;
  int nxy      = nel_frac*Nm;
  int ntotp    = nxy*nprocs;
  int ntotz    = nxy*nztot;
  int ntot     = Nm*U->nel;

  int Nq       = U->fhead->qa*U->fhead->qb;
  int nxy_q    = nel_frac*Nq;
  int ntotp_q  = nxy_q*nprocs;
  int ntotz_q  = nxy_q*nztot;
  int ntot_q   = Nq*U->nel;

  int  i,j,k;

  if(!pstore){
    pstore = dmatrix(0, 2, 0, U->hjtot-1);
    
    uz = dvector(0, nxy*nztot-1); dzero(nxy*nztot, uz, 1);
    vz = dvector(0, nxy*nztot-1); dzero(nxy*nztot, vz, 1);
    wz = dvector(0, nxy*nztot-1); dzero(nxy*nztot, wz, 1);
    				  
    u  = dvector(0, nxy*nztot-1); dzero(nxy*nztot, u, 1);
    v  = dvector(0, nxy*nztot-1); dzero(nxy*nztot, v, 1);
    w  = dvector(0, nxy*nztot-1); dzero(nxy*nztot, w, 1);
    
    uF = dvector(0, nxy_q*nztot-1); dzero(nxy_q*nztot, uF, 1);
    vF = dvector(0, nxy_q*nztot-1); dzero(nxy_q*nztot, vF, 1);
    wF = dvector(0, nxy_q*nztot-1); dzero(nxy_q*nztot, wF, 1);
  }    

  dcopy(U->hjtot, U->base_hj, 1, pstore[0], 1);
  dcopy(U->hjtot, V->base_hj, 1, pstore[1], 1);
  dcopy(U->hjtot, W->base_hj, 1, pstore[2], 1);


  double  *fu = dvector(0, Nq-1);
  double  *fv = dvector(0, Nq-1);
  double  *fw = dvector(0, Nq-1);

  double  *dx = dvector(0, Nq-1),   *dy = dvector(0, Nq-1);

  if(!Pput_new){
    Pput_new = dvector(0, nxy_q*nztot-1);
    Pget_new = dvector(0, nxy_q*nztot-1);
  }		       

  // packs u into ntot blocks of contiguous nz points
  packf   (      nz,   ntot,      U->base_hj, Pput_new);
  Exchange(   ntotz,   Pput_new,  Pget_new);
  unpackp (      nz,   ntotp,     Pget_new,   u);
  grad_z_new( nztot,   nxy,       u,          uz);

  multifft(nztot, nxy,  u, 1);
  multifft(nztot, nxy, uz, 1);

  packf   (      nz,   ntot,      V->base_hj, Pput_new);
  Exchange(   ntotz,   Pput_new,  Pget_new);
  unpackp (      nz,   ntotp,     Pget_new,   v);
  grad_z_new( nztot,   nxy,       v,          vz);

  multifft(nztot, nxy,  v, 1);
  multifft(nztot, nxy, vz, 1);
  
  packf   (      nz,   ntot,      W->base_hj, Pput_new);
  Exchange(   ntotz,   Pput_new,  Pget_new);
  unpackp (      nz,   ntotp,     Pget_new,   w);
  grad_z_new( nztot,   nxy,       w,          wz);

  multifft(nztot, nxy,  w, 1);
  multifft(nztot, nxy, wz, 1);
  
  int N = min(U->nel, (procid+1)*nel_frac);
  for(j=0;j<nztot;++j)
    for(i = procid*nel_frac,k=0; i < N; ++i,++k){
      eU = U->flist[i];      eV = V->flist[i];      eW = W->flist[i];

      dcopy(Nm, uz+j+k*Nm*nztot, nztot,  eU->vert[0].hj, 1);
      dcopy(Nm, vz+j+k*Nm*nztot, nztot,  eV->vert[0].hj, 1);
      dcopy(Nm, wz+j+k*Nm*nztot, nztot,  eW->vert[0].hj, 1);
      eU->Trans(eU, J_to_Q);
      eV->Trans(eV, J_to_Q);
      eW->Trans(eW, J_to_Q);
      dcopy(Nq, eU->h[0], 1, fu, 1);
      dcopy(Nq, eV->h[0], 1, fv, 1);
      dcopy(Nq, eW->h[0], 1, fw, 1);
      
      dcopy(Nm, u+j+k*Nm*nztot, nztot, eU->vert[0].hj, 1);
      dcopy(Nm, v+j+k*Nm*nztot, nztot, eV->vert[0].hj, 1);
      dcopy(Nm, w+j+k*Nm*nztot, nztot, eW->vert[0].hj, 1);
      eU->Trans(eU, J_to_Q);
      eV->Trans(eV, J_to_Q);
      eW->Trans(eW, J_to_Q);
      
      eU->Grad_d(dx, dy, 0, 'a');
      dvmul (Nq, eW->h[0], 1, fu, 1, fu, 1);
      dvvtvp(Nq, eU->h[0], 1, dx, 1, fu, 1, fu, 1);
      dvvtvp(Nq, eV->h[0], 1, dy, 1, fu, 1, fu, 1);
      dcopy (Nq, fu, 1, uF+j+k*Nq*nztot, nztot);

      eV->Grad_d(dx, dy, 0, 'a');
      dvmul (Nq, eW->h[0], 1, fv, 1, fv, 1);
      dvvtvp(Nq, eU->h[0], 1, dx, 1, fv, 1, fv, 1);
      dvvtvp(Nq, eV->h[0], 1, dy, 1, fv, 1, fv, 1);
      dcopy (Nq, fv, 1, vF+j+k*Nq*nztot, nztot);
      
      eW->Grad_d(dx, dy, 0, 'a');
      dvmul (Nq, eW->h[0], 1, fw, 1, fw, 1);
      dvvtvp(Nq, eU->h[0], 1, dx, 1, fw, 1, fw, 1);
      dvvtvp(Nq, eV->h[0], 1, dy, 1, fw, 1, fw, 1);
      dcopy (Nq, fw, 1, wF+j+k*Nq*nztot, nztot);
    }

  dscal(nxy_q*nztot, -1.0, uF, 1);
  
  multifft(nztot, nxy_q,  uF, -1);
  packp   (nz,  ntotp_q,       uF,  Pput_new);
  Exchange(ntotz_q,      Pput_new,  Pget_new);
  unpackf (nz,  ntot_q,  Pget_new,  Uf->base_h);

  dscal(nxy_q*nztot, -1.0, vF, 1);

  multifft(nztot, nxy_q,  vF, -1);
  packp   (nz,  ntotp_q,       vF,  Pput_new);
  Exchange(ntotz_q,      Pput_new,  Pget_new);
  unpackf (nz,  ntot_q,  Pget_new,  Vf->base_h);

  dscal(nxy_q*nztot, -1.0, wF, 1);

  multifft(nztot, nxy_q,  wF, -1);
  packp   (nz,  ntotp_q,       wF,  Pput_new);
  Exchange(ntotz_q,      Pput_new,  Pget_new);
  unpackf (nz,  ntot_q,  Pget_new,  Wf->base_h);

  Uf->Set_state('p');
  Vf->Set_state('p');
  Wf->Set_state('p');

  dcopy(U->hjtot, pstore[0], 1, U->base_hj, 1);   U->Trans(U, J_to_Q);
  dcopy(U->hjtot, pstore[1], 1, V->base_hj, 1);   V->Trans(V, J_to_Q);
  dcopy(U->hjtot, pstore[2], 1, W->base_hj, 1);   W->Trans(W, J_to_Q);


  free(dx); free(dy);
  free(fu); free(fv); free(fw);
}
#endif

static void pwriteHeaderF (FILE *fp, Field *f, Element *E)
{
  register int i;
  char     buf[128];
#ifdef _CRAY
  short *ssize;
  register int n;
#endif

  /* Write the header */
  fprintf (fp, fourier_hdr_fmt[0],  f->name);
  fprintf (fp, fourier_hdr_fmt[1],  f->created);
  fprintf (fp, fourier_hdr_fmt[2],  f->state);
  fprintf (fp, fourier_hdr_fmt[3],  f->nel, f->dim, f->lmax, f->nz);
  
  fprintf (fp, fourier_hdr_fmt[4],  f->step);
  fprintf (fp, fourier_hdr_fmt[5],  f->time);
  fprintf (fp, fourier_hdr_fmt[6],  f->time_step);
  
  fprintf (fp, fourier_hdr_fmt[7],  f->kinvis, f->lz);
  fprintf (fp, fourier_hdr_fmt[8],  f->type);
  fprintf (fp, fourier_hdr_fmt[9],  f->format);

  int cnt = count_facets(E);
  switch (tolower(*f->format)) {
  case 'a':
    for(i = 0; i < cnt; ++i)
      fprintf(fp,"%d ",f->size[i]);
    fputc('\n',fp);
    break;
    
  case 'b':
#ifdef _CRAY
    // int on Crays is 64b, short 32b - use short for portable fieldfiles
    ssize = (short *) calloc(cnt, sizeof(short));
    for (n = 0; n < cnt; n++)
      ssize[n] = (short) f->size[n];
    fwrite(ssize, sizeof(short),  cnt, fp);
    free (ssize);
#else
    fwrite(f->size, sizeof(int),  cnt, fp);
#endif
    break;
    
  default:
    sprintf  (buf, "unknown format -- %s", f->format);
    fprintf(stderr,buf);
    exit(-1);
    break;
  }

  fflush (fp);
  return;
}

static void pAwriteDataF (FILE *fp, Field *f, Element *E)
{
  int  nfields, ntot;
  register int i, n;

  /* Write the field files  */
  if(f->state == 't'){
    ntot = option("NZ")*data_len(f->nel,f->size,E);
  }
  else{
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
      
  nfields = (int) strlen (f->type);

  for (i = 0; i < ntot; i++) {
    for (n = 0; n < nfields; n++)
      fprintf (fp, "%#16.10g ", f->data[n][i]);
    fputc ('\n', fp);
  }
  
  return;
}

static void pBwriteDataF (FILE *fp, Field *f, int field, Element *E)
{
  int ntot;
  int procid = option("PROCID");

  /* Write the field files  */
  if(f->state == 't'){
    ntot = option("NZ")*data_len(f->nel,f->size,E);
  }
  else{
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
      
  if (fwrite (f->data[field], sizeof(double), ntot, fp) != ntot) {
    fprintf(stderr, "Proc %d: error writing field %c", procid, f->type[field]);
    exit(-1);
  }
  
  return;
}

void pWriteFieldF(FILE *fp, char *name, char *fname, int step, double t,
		  int nfields, Element_List *U[]){
  register int n,i,j;
  time_t   tp;
  char     buf[128],state;
  Field    f;
  int      nztot  = U[0]->nztot;
  
  if(nfields > MAXFIELDS)
    error_msg(too many fields -- must be less than MAXFIELDS);
  
  state = U[0]->fhead->state;
  /* check to see if all fields in same state */
  for(i = 1; i < nfields; ++i)
    if(state != U[i]->fhead->state)
      {error_msg(Writefield--Fields  in different space);}
      
  /* Get the current date and time from the system */
  
  tp = time((time_t*)NULL);
  strftime(buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  
  /* Set up the field structure */
  f.name      =  name;
  f.created   =  buf;
  f.state     =  state;
  f.dim       =  DIM;
  f.nel       =  U[0]->nel;
  f.lmax      =  LGmax;
  f.step      =  step;
  f.time      =  t;
  f.nz        =  nztot;
  f.lz        =  dparam("LZ");
  f.time_step =  dparam("DT");
  f.kinvis    =  (dparam("KINVIS-ORIG"))?
    dparam("KINVIS-ORIG") : dparam("KINVIS");   /* VARV */
  f.format    = (char *)((option("binary")) ? "binary-"BTYPE : "ascii");

  /* Generate the list of fields and assign data pointers */
  
  memset (f.type, '\0', MAXFIELDS);
  memset (f.data, '\0', MAXFIELDS * sizeof( double *));
  
  /* Set up sizes & data links */
  if(state == 't'){
    int cnt = count_facets(U[0]->fhead);
    f.size = ivector(0,cnt-1);
    for(i = 0,n=0; i < f.nel; ++i){
      for(j = 0; j < U[0]->flist[i]->Nedges; ++j,++n)
	f.size[n] = U[0]->flist[i]->edge[j].l;
      
      for(j = 0; j < U[0]->flist[i]->Nfaces; ++j,++n)
	f.size[n] = U[0]->flist[i]->face[j].l;
    }
    for(n = 0; n < nfields; ++n){
      f.type [n] = U[n]->fhead->type;
      f.data [n] = U[n]->base_hj;
    }
  } else {
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
  
  ROOT pwriteHeaderF(fp,&f,U[0]->fhead);
  
  FILE *Fp;
  char Fname[FILENAME_MAX];
  int procid = option("PROCID"), err;

  switch (tolower(*f.format)) {
  case 'a':
    sprintf (Fname, "%s.a%d", fname, procid);
    err = backup(Fname);
    if (err != 0 && step != 0) 
      fprintf(stderr,"Proc %d - WARNING: Analyser: failed to backup the chk file %s\n", procid, Fname);
    Fp = fopen(Fname,"w");
    pAwriteDataF (Fp, &f, U[0]->fhead);
    fclose (Fp);
    break;
    
  case 'b':
    for (i = 0; i < nfields; i++) {
      sprintf (Fname, "%s.%c%d", fname, f.type[i], procid);
      err = backup(Fname);
      if (err != 0 && step != 0) 
	fprintf(stderr,"Proc %d - WARNING: Analyser: failed to backup the chk file %s\n", procid, Fname);
      Fp = fopen(Fname,"w");
      pBwriteDataF (Fp, &f, i, U[i]->fhead);
      fclose (Fp);
    }
    break;
    
  default:
    sprintf  (buf, "unknown format -- %s", f.format);
    fprintf(stderr,buf);
    exit(-1);
    break;
  }

  free(f.size);

  return;
}

void pMWriteFieldF(FILE *fp, char *name, char *fname, int step, double t,
		  int nfields, Element_List *U, double **mat, char *types)
{
  register int n,i,j;
  time_t   tp;
  char     buf[128],state;
  Field    f;
  int      nztot  = U->nztot;
  
  if(nfields > MAXFIELDS)
    error_msg(too many fields -- must be less than MAXFIELDS);
  if(nfields != (int) strlen(types))
    error_msg(more fields than field types specified);
  
  state = U->fhead->state;
      
  /* Get the current date and time from the system */
  
  tp = time((time_t*)NULL);
  strftime(buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  
  /* Set up the field structure */
  f.name      =  name;
  f.created   =  buf;
  f.state     =  state;
  f.dim       =  DIM;
  f.nel       =  U->nel;
  f.lmax      =  LGmax;
  f.step      =  step;
  f.time      =  t;
  f.nz        =  nztot;
  f.lz        =  dparam("LZ");
  f.time_step =  dparam("DT");
  f.kinvis    =  (dparam("KINVIS-ORIG"))?
    dparam("KINVIS-ORIG") : dparam("KINVIS");   /* VARV */
  f.format    = (char *)((option("binary")) ? "binary-"BTYPE : "ascii");

  /* Generate the list of fields and assign data pointers */
  
  memset (f.type, '\0', MAXFIELDS);
  memset (f.data, '\0', MAXFIELDS * sizeof( double *));
  
  /* Set up sizes & data links */
  if(state == 't'){
    int cnt = count_facets(U->fhead);
    f.size = ivector(0,cnt-1);
    for(i = 0,n=0; i < f.nel; ++i){
      for(j = 0; j < U->flist[i]->Nedges; ++j,++n)
	f.size[n] = U->flist[i]->edge[j].l;
      
      for(j = 0; j < U->flist[i]->Nfaces; ++j,++n)
	f.size[n] = U->flist[i]->face[j].l;
    }
    for(n = 0; n < nfields; ++n){
      f.type [n] = types[n];
      f.data [n] = mat[n];
    }
  } else {
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }

  ROOT pwriteHeaderF(fp,&f,U->fhead);

  FILE *Fp;
  char Fname[FILENAME_MAX];
  int procid = option("PROCID"), err;

  switch (tolower(*f.format)) {
  case 'a':
    sprintf (Fname, "%s.a%d", fname, procid);
    err = backup(Fname);
    if (err != 0 && step != 0) 
      fprintf(stderr,"Proc %d - WARNING: Analyser: failed to backup the chk file %s\n", procid, Fname);
    Fp = fopen(Fname,"w");
    pAwriteDataF (Fp, &f, U->fhead);
    fclose (Fp);
    break;
    
  case 'b':
    for (i = 0; i < nfields; i++) {
      sprintf (Fname, "%s.%c%d", fname, f.type[i], procid);
      err = backup(Fname);
      if (err != 0 && step != 0) 
	fprintf(stderr,"Proc %d - WARNING: Analyser: failed to backup the chk file %s\n", procid, Fname);
      Fp = fopen(Fname,"w");
      pBwriteDataF (Fp, &f, i, U->fhead);
      fclose (Fp);
    }
    break;
    
  default:
    sprintf  (buf, "unknown format -- %s", f.format);
    fprintf(stderr,buf);
    exit(-1);
    break;
  }

  free(f.size);

  return;
}

void MWritefieldF(FILE *fp, char *name, int step, double t,
		 int nfields, Element_List *U, double **mat, char *types)
{
  register int n,i,j;
  time_t   tp;
  char     buf[128],state;
  Field    f;
  int      ntot;
  double   **u;
  int      nz     = U->nz;
  int      nztot  = U->nztot;
  int      procid = option("PROCID"), nprocs = option("NPROCS");

  if(nfields > MAXFIELDS)
    error_msg(too many fields -- must be less than MAXFIELDS);
  if(nfields != (int) strlen(types))
    error_msg(more fields than field types specified);

  state = U->fhead->state;

  /* Get the current date and time from the system */
  
  tp = time((time_t*)NULL);
  strftime(buf, 25, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  
  /* Set up the field structure */
  f.name      =  name;
  f.created   =  buf;
  f.state     =  state;
  f.dim       =  DIM;
  f.nel       =  U->nel;
  f.lmax      =  LGmax;
  f.step      =  step;
  f.time      =  t;
  f.nz        =  nztot;
  f.lz        =  dparam("LZ");
  f.time_step =  dparam("DT");
  f.kinvis    =  (dparam("KINVIS-ORIG"))?
    dparam("KINVIS-ORIG") : dparam("KINVIS");   /* VARV */
  f.format    = (char *)((option("binary")) ? "binary-"BTYPE : "ascii");

  /* Generate the list of fields and assign data pointers */
  
  memset (f.type, '\0', MAXFIELDS);
  memset (f.data, '\0', MAXFIELDS * sizeof( double *));
  
  /* copyfield into matrix */
  if(state == 't'){
    ROOT {
      int cnt = count_facets(U->fhead);
      f.size = ivector(0,cnt-1);
      for(i = 0,n=0; i < f.nel; ++i){
	for(j = 0; j < U->flist[i]->Nedges; ++j,++n)
	  f.size[n] = U->flist[i]->edge[j].l;
	
	for(j = 0; j < U->flist[i]->Nfaces; ++j,++n)
	  f.size[n] = U->flist[i]->face[j].l;
      }
    }
    ntot = U->hjtot*nz;

    if(!procid)
      u = dmatrix(0,nfields-1,0,nprocs*ntot-1);    
    
    for(n = 0; n < nfields; ++n){
#ifdef PARALLEL
      if (procid) {
	MPI_Send (mat[n], ntot, MPI_DOUBLE, 0, W_MSG+procid, MPI_COMM_WORLD);
      } else {	
	dcopy(ntot, mat[n], 1, u[n], 1);
	MPI_Status status;   
	for (i = 1; i < nprocs; ++i)
	  MPI_Recv (u[n]+i*ntot, ntot, MPI_DOUBLE, 
		    i, W_MSG+i, MPI_COMM_WORLD, &status);
      }
#else
      dcopy(ntot, mat[n], 1, u[n], 1);
#endif

    }
  } else {
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }

  ROOT {
    for(n = 0; n < nfields; ++n){
      f.type [n] = types[n];
      f.data [n] = u[n];
    }
    writeFieldF(fp,&f,U->fhead);
    free_dmatrix(u,0,0);
    free(f.size);
  }

  return;
}


