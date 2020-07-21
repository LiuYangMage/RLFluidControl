/*
 * Wannier Flow 
 *
 * This is benchmark problem for Prism.  This file sets up the solution
 * for Wannier flow, an exact solution to the Stokes equations for a cylinder
 * spinning above a wall.  The boundary conditions are complicated enough
 * to be hard-coded in this file.
 * 
 * ------------------------------------------------------------------------- */
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include "nektarF.h"
#include <polylib.h>

static char     dir;

#ifdef WANNIER
#if DIM == 2
static void wannier(int n, double *x, double *y, double *u, char dir);
#else
static void wannier(int n, double *x, double *y, double *z, double *u,
		    char dir);
#endif
/* 
 * Redefine vec_init() and set_vector() 
 * ---------------------------------------------------------------------- */
 
void vector_def(char *s, char *f)
{
  while (isspace(*f++));
  dir = *--f;
  return;
}

void vector_set(int n, ...)
{
#if DIM == 2
  double  *x, *y, *f;
#else
  double  *x, *y, *z, *f;
#endif
  va_list ap;

  va_start(ap, n);
  x = va_arg(ap, double *);
  y = va_arg(ap, double *);
#if DIM == 3
  z = va_arg(ap, double *);
#endif
  f = va_arg(ap, double *);
  va_end  (ap);

#if DIM == 2
  wannier(n, x, y, f, dir);
#else
  wannier(n, x, y, z, f, dir);
#endif
  return;
}

/* ---------------------------------------------------------------------- */

#if DIM == 2
static void wannier(int n, double *x, double *y, double *u, char dir)
#else
static void wannier(int n, double *x, double *y, double *z, double *u,
		    char dir)
#endif
{
  double K1, K2, Y1, Y2;

  static double U     = 1.0,
                V     = 0.0, /* 0.5,*/
                R     = 0.25,
                d     = 0.5,
                s     = 0.0;
  static double A, B, C, D, F, G;
  
  register int i;

  if (s == 0.0) {          /* initialization */
    s = sqrt(d*d - R*R);
    G = (d+s) / (d-s);
    A = -(U*d / log(G) + 0.5 * d*R*V/s);
    B = 2.0*(d+s)*U/log(G) + (d+s)*R*V/s;
    C = 2.0*(d-s)*U/log(G) + (d-s)*R*V/s;
    D = -U;
    F = U / log(G);
  }

  switch (dir) {
  case 'u':
    for (i = 0; i < n; i++) {
      Y1    = y[i] + d;
      Y2    = 2.0 * Y1;
      K1    = x[i]*x[i] + (s+Y1)*(s+Y1);
      K2    = x[i]*x[i] + (s-Y1)*(s-Y1);
      
      u[i]  = -2.0/K1*(A + F*Y1) * (s+Y1 + K1/K2*(s-Y1)) - F*log(K1/K2)
	      - B/K1*( s + Y2 - Y2 * (s+Y1)*(s+Y1)/K1)
	      - C/K2*( s - Y2 + Y2 * (s-Y1)*(s-Y1)/K2) - D;
    }
    break;

  case 'v':
    for (i = 0; i < n; i++) {
      Y1    = y[i] + d;
      Y2    = 2.0 * Y1;
      K1    = x[i]*x[i] + (s+Y1)*(s+Y1);
      K2    = x[i]*x[i] + (s-Y1)*(s-Y1);

      u[i]  = 2.0*x[i]/(K1*K2) * (A + F*Y1) * (K2 - K1) -
              x[i]*B*Y2*(s+Y1)/(K1*K1) - C*x[i]*Y2*(s-Y1)/(K2*K2);
    }
    break;
  case 'w':
    for (i = 0; i < n; i++) 
      u[i] = 0.0;
    break;
  case 'p':
    break;
  default:
    printf("unknown direction in wannier -- %c\n", dir);
    break;
  }
  
  return;
}
#endif


#ifdef BATCHELOR

#if DIM == 2
static void batchelor(int n, double *x, double *y, double *u, char dir);
#else
static void batchelor(int n, double *x, double *y, double *z, double *u,
                    char dir);
#endif
/* 
 * Redefine vec_init() and set_vector() 
 * ---------------------------------------------------------------------- */
 
void vector_def(char *s, char *f)
{
  while (isspace(*f++));
  dir = *--f;
  return;
}

void vector_set(int n, ...)
{
#if DIM == 2
  double  *x, *y, *f;
#else
  double  *x, *y, *z, *f;
#endif
  va_list ap;

  va_start(ap, n);
  x = va_arg(ap, double *);
  y = va_arg(ap, double *);
#if DIM == 3
  z = va_arg(ap, double *);
#endif
  f = va_arg(ap, double *);
  va_end  (ap);

#if DIM == 2
  batchelor(n, x, y, f, dir);
#else
  batchelor(n, x, y, z, f, dir);
#endif
  return;
}

/* ---------------------------------------------------------------------- */

#if DIM == 2
static void batchelor(int n, double *x, double *y, double *u, char dir)
#else
static void batchelor(int n, double *x, double *y, double *z, double *u,
                    char dir)
#endif
{
  double TOL = 1e-8;
  double L = dparam("L");
  double q = dparam("q");
  register int i;

  switch (dir) {
  case 'u':
    for (i = 0; i < n; i++) {
      if(((x[i]-1.0)*(x[i]-1.0)+y[i]*y[i]) < TOL){
        u[i] = 0.0;
      }
      else if(((x[i]+1.0)*(x[i]+1.0)+y[i]*y[i]) < TOL){
        u[i] = 0.0;
      }
      else{
        sxi = L/2/M_PI*sin(2*M_PI*x[i]/L) - 1;
        sxj = L/2/M_PI*sin(2*M_PI*x[i]/L) + 1;   
        sy  = L/2/M_PI*sin(2*M_PI*y[i]/L);
        u[i] = q*(1-exp(-(sxi*sxi + sy*sy)))*sy/(sxi*sxi + sy*sy) - 
            q*(1-exp(-(sxj*sxj + sy*sy)))*sy/(sxj*sxj + sy*sy);

      }
    }
    break;
  case 'v':
    for (i = 0; i < n; i++) {
      if(((x[i]-1.0)*(x[i]-1.0)+y[i]*y[i]) < TOL){
        u[i] = 0.0;
      }
      else if(((x[i]+1.0)*(x[i]+1.0)+y[i]*y[i]) < TOL){
        u[i] = 0.0;
      }
      else{
        sxi = L/2/M_PI*sin(2*M_PI*x[i]/L) - 1;
        sxj = L/2/M_PI*sin(2*M_PI*x[i]/L) + 1;   
        sy = L/2/M_PI*sin(2*M_PI*y[i]/L);
        u[i] = -q*(1-exp(-(sxi*sxi + sy*sy)))*sxi/(sxi*sxi + sy*sy) + 
            q*(1-exp(-(sxj*sxj + sy*sy)))*sxj/(sxj*sxj + sy*sy);
      }
    }
    break;
#if DIM == 3
  case 'w':
    for (i = 0; i < n; i++){
      sxi = L/2/M_PI*sin(2*M_PI*x[i]/L) - 1;
      sxj = L/2/M_PI*sin(2*M_PI*x[i]/L) + 1;   
      sy = L/2/M_PI*sin(2*M_PI*y[i]/L)*sin(2*M_PI*y[i]/L);
      u[i]  = exp(-(sxi*sxi + sy*sy))+ 
        exp(-(sxj*sxj + sy*sy));
    }
    break;
#endif
  case 'p':
    break;
  case '0':
    for (i = 0; i < n; i++) 
      u[i] = 0.0;
    break;
  default:
    printf("unknown direction in wannier -- %c\n", dir);
    break;
  }
  
  return;
}
#endif

#ifdef HEIMINEZ
/*
 * Heiminez Flow 
 *
 * ------------------------------------------------------------------------- */

#if DIM == 2
static void heiminez(int n, double *x, double *y, double *u, char dir);
#else
static void heiminez(int n, double *x, double *y, double *z, double *u,
		     char dir);
#endif
		     
/* 
 * Redefine vec_init() and set_vector() 
 * ---------------------------------------------------------------------- */
 
void vector_def(char *s, char *f)
{
  while (isspace(*f++));
  dir = *--f;
  return;
}

void vector_set(int n, ...)
{
#if DIM == 2
  double  *x, *y, *f;
#else
  double  *x, *y, *z, *f;
#endif
  va_list ap;

  va_start(ap, n);
  x = va_arg(ap, double *);
  y = va_arg(ap, double *);
#if DIM == 3
  z = va_arg(ap, double *);
#endif
  f = va_arg(ap, double *);
  va_end  (ap);

#if DIM == 2
  heiminez(n, x, y, f, dir);
#else
  heiminez(n, x, y, z, f, dir);
#endif
  return;
}

/* ---------------------------------------------------------------------- */

#if DIM == 2
static void heiminez(int n, double *x, double *y, double *u, char dir)
#else
static void heiminez(int n, double *x, double *y, double *z, double *u,
		    char dir)
#endif
{
  register int i,j;
  double *f, *fd, *eta, **imat;
  static double *zorig,*F,*Fd,len,scle;
  static int init = 1,npts;
  double Wvel = dparam("WVEL");
  
  if (init) {          /* initialization - read in solution */
    char buf[BUFSIZ];
    double ybeg,yend;
    FILE *fp = fopen("Heiminez","r");

    if(!fp) {
      fprintf(stderr,"Can not find file Heiminez in function heiminez\n");
      exit(1);
    }

    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%d %lf %lf",&npts,&ybeg,&yend);
    zorig = dvector(0,npts-1);
    F     = dvector(0,npts-1);
    Fd    = dvector(0,npts-1);

    /* Generate gauss-lobatto-chebychev points */
    zwglc(zorig,F,npts);

    for(i = 0; i < npts; ++i)
      fscanf(fp,"%*lf%lf%lf\n",F+i,Fd+i);

    len  = 0.1*(yend-ybeg);
    scle = 2.0*len/(yend-ybeg);

    init = 0;
  }

  if(dir != '0'){
    /* interpolate F and Fd from Y points to y points */
    eta = dvector(0,n-1);
    f   = dvector(0,n-1);
    fd  = dvector(0,n-1);
    
    /* calculate eta location from y */
    dsmul(n,1+scle,y,1,eta,1);
    dsadd(n,-len,eta,1,eta,1);
    dsadd(n, len,y,1,f,1);
    dvdiv(n,eta,1,f,1,eta,1);
    
    /* make interpolation matrix */
    imat = dmatrix(0,n-1,0,npts-1);
    iglcm(imat,zorig,eta,npts,n);
    
    for(i = 0; i < n; ++i){
      f [i] = ddot(npts,imat[i],1,F,1);
      fd[i] = ddot(npts,imat[i],1,Fd,1);
    }
  }

  switch (dir) {
  case 'u':
    for (i = 0; i < n; i++) 
      u[i]  = -x[i]*fd[i];
    break;
  case 'd': /* du/dn = du/dx */
    for(i = 0; i < n; i++) 
      u[i]  = -fd[i];
    break;
  case 'e': /* du/dn = -du/dx */
    for(i = 0; i < n; i++) 
      u[i]  = fd[i];
    break;
  case 'v':
    for (i = 0; i < n; i++) 
      u[i] = f[i];
    break;
  case 'w':
    for (i = 0; i < n; i++) 
      u[i] = Wvel;
    break;
  case 'p':
    break;
  case '0':
    for (i = 0; i < n; i++) 
      u[i] = 0.0;
    break;
  default:
    printf("unknown direction in heiminez -- %c\n", dir);
    break;
  }

  if(dir != '0'){
    free(f); free(fd); free(eta); free_dmatrix(imat,0,0);
  }

  return;
}
#endif

