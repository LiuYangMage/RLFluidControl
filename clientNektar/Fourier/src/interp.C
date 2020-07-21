#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nektarF.h"
#include <veclib.h>

/* Private functions */

static void Interp_pts(Intepts *I, Element_List **U, int nfields);
static void Fourier_interp(Element_List **U, Coord *X, double *ans, double *ufr, int nfields);

void interp(Domain *omega){
  int nfields;
  Element_List **E;

  nfields = DIM + 2;
  E = (Element_List**) malloc(nfields*sizeof(Element*));

  E[0] = omega->U;
  E[1] = omega->V;
  E[2] = omega->W;
  E[3] = omega->P;

  Interp_pts(omega->int_list, E, nfields);

  free(E); 

  return;
}

static void Interp_pts(Intepts *I, Element_List **U, int nfields){
  register int i,j,k;
  double  *ui, *ufr;
  Coord   X;
  int     nztot = option("NZTOT");
 
  ui  = dvector(0,nfields-1);
  ufr = dvector(0,nfields*nztot - 1);
  X.x = dvector(0,DIM+1);
  
  X.y = X.x+1;
  X.z = X.y+1;

  if (I->npts > 0){ 
    for(i = 0; i < I->npts; ++i){
      X.x[0] = I->X.x[i]; 
      X.y[0] = I->X.y[i];
      X.z[0] = I->X.z[i];

      Fourier_interp(U, &X, ui, ufr, nfields);
      
      ROOT{
	for(j = 0; j < nfields; ++j) {
	  I->ui[i][j] = ui[j];
	  for(k = 0; k < nztot; k++)
	    I->ufr[i][j*nztot+k] = ufr[j*nztot+k];
	}
      }
    }
  }

  free(ui);
  free(ufr);
  free(X.x); 
}

static void Fourier_interp(Element_List **U, Coord *X, double *ans, double *ufr, int nfields){
  int     nz = option("NZ"), nztot = option("NZTOT"), i, j;
  double* ui = dvector(0, nfields*nz), *tmp;
  Element_List **V = (Element_List**) malloc(U[0]->nz*nfields*sizeof(Element_List*));
  
  ROOT tmp = dvector(0, nfields*nztot-1);

  for(j=0;j<nfields;++j)
    for(i=0;i<U[0]->nz;++i)
      V[j+i*nfields] = U[j]->flevels[i];

  dzero(nfields, ans, 1);
  
  Interp_point(V, nfields*U[0]->nz, X, ui);

#ifdef PARALLEL
  MPI_Gather (ui, nfields*nz, MPI_DOUBLE, tmp, nfields*nz, MPI_DOUBLE, 0, 
			MPI_COMM_WORLD);
#else
  dcopy(nfields*nz,ui,1,tmp,1);
#endif

  ROOT{
    for(j = 0;j < nfields; ++j){
      register double angle;
      // This is supposedly zero
      // ans[j] = cos(Global_Beta(nztot)*X->z[0])*tmp[j+nfields];
      ufr[j*nztot+1] = 0.;
      ans[j] += tmp[j];
      ufr[j*nztot] = tmp[j];
      for(i = 2; i < nztot; i+=2) {
	angle = Beta(i)*X->z[0]; // on proc 0, Beta(i) is Global_Beta(i)
	ans[j] += 2.*(cos(angle)*tmp[j+i*nfields]
		      -sin(angle)*tmp[j+(i+1)*nfields]);
	ufr[j*nztot+i] = tmp[j+i*nfields];
	ufr[j*nztot+i+1] = tmp[j+(i+1)*nfields];
      }
    }
  
    free(tmp);
  } 
  
  free(ui);
  free(V);
}
