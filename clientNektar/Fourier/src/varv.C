/* ------------------------------------------------------------------------- *
 * 		Large Eddy Simulation Routines - Smagorinsky 	             *
 * The following functions calculate the Large Eddy Simulation additions     *
 * for the Smagorinsky modelling.					     * 
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
/* calculate the variable viscosity */
typedef  enum {
  M_Smag,
  M_Bardina
} ModelType;
static varv_variables *preprocessLES(Domain *omega, double *Diam);
void tijCalc(Domain *omega,varv_variables *vvexp, double *Diam);
extern double **u2nde;

varv_variables *CsCalc(Domain *omega){
  Element_List *U  =  omega->U,  *V  =  omega->V, *W  =  omega->W;
  int    nq =U->htot*U->nz;
  int model = (ModelType) iparam("MODEL");
  double  *Diam = dvector(0,U->htot);
  varv_variables *vvexp;
  
  if (option("dealias")) {
    nq = 3*U->htot*U->nz/2;
  }  
  
  /*  switch(model){
      case M_Smag:*/
  vvexp = preprocessLES(omega,Diam);
  tijCalc(omega,vvexp,Diam);
  /*    break;
	case M_Bardina::
	Bardina(omega,vvexp);
	break;
	}*/
  free(Diam);
  return vvexp;
}

void tijCalc(Domain *omega,varv_variables *vvexp, double *Diam)
{
  Element_List *U  =  omega->U,  *V  =  omega->V,  *W  = omega->W,
    *Ux =  omega->P,  *Uy = omega->Uf, *Uz = omega->Uf,
    *Vx =  omega->P,  *Vy = omega->Vf, *Vz = omega->Vf,
    *Wx =  omega->P,  *Wy = omega->Wf, *Wz = omega->Wf,
    *nul = (Element_List*)0, *Div = omega->P;
  FILE *fp = 0;
  char        trip = 'a';
  char      fname[FILENAME_MAX];
  static Nek_Trans_Type f_to_p = F_to_P,p_to_f = P_to_F;
  double *NLx  =  vvexp->NLx , *NLy  = vvexp->NLy,*NLz = vvexp->NLz,*dudz,
    *dvdz,*dwdz,*dwdx,*dwdy,Sijz,Sz,Miuz,D3,*dudx,*dudy,*dvdx,*dvdy,Sij,S,Miu,
    *D,cs,**tempo,*nu,*nux, *nuy,**tmp,MaxV=0.0,MinV = 0.0,AvV=0.0;	
  int nz= option("NZ"), nq=U->htot*U->nz;
  register int i,k;
  if (option("dealias")) {
    f_to_p = F_to_P32;
    p_to_f = P_to_F32;
    nq = 3*U->htot*U->nz/2;
  }
  tempo = dmatrix(0,16,0,nq); tmp  = dmatrix(0,8,0,nq);
  D = dvector(0,nq); nu  = dvector(0,nq);  nux = dvector(0,nq);
  nuy = dvector(0,nq); dudx= tmp[0]; dudy= tmp[1]; dvdx= tmp[2]; dvdy= tmp[3]; 
  dudz= tmp[4]; dvdz= tmp[5]; dwdz= tmp[6]; dwdy= tmp[7]; dwdx= tmp[8];

  if (option("dealias"))
    trip = 'A';
  
  U->Grad_h(U->base_h,dudx,dudy,(double *)NULL,trip);
  V->Grad_h(V->base_h,dvdx,dvdy,(double *)NULL,trip);
  W->Grad_h(W->base_h,dwdx,dwdy,(double *)NULL,trip);

  register int iq = 0;
  for (i = 0; i < nq; ++i)
    {
      if (iq == U->htot) iq = 0;

      Sij = 0.5 * (	2 * dudx[i] * dudx[i] + dudy[i] * dudy[i] + 
			Uz->base_h[i] * Uz->base_h[i] + dvdx[i] * dvdx[i] +
			2 * dvdy[i] * dvdy[i] + Vz->base_h[i] * 
			Vz->base_h[i] + dwdx[i] * dwdx[i] + dwdy[i] * dwdy[i] +
			2 * Wz->base_h[i]*Wz->base_h[i]+ 
			2 * Uz->base_h[i] * dwdx[i] + 2 * dudy[i] * dvdx[i] +
			2 * dwdy[i] * Vz->base_h[i] ) ; 
      S = sqrt(2 * Sij);  

      cs = dparam("SMAGS")*dparam("SMAGS");
      Miu = cs * Diam[iq]*Diam[iq] * S;  
      if (dparam("CAP") != 0.0)
	{
	  if (sqrt(Miu*Miu) > dparam("CAP")*dparam("KINVIS-ORIG")) 
	    Miu = dparam("CAP")*dparam("KINVIS-ORIG")*sqrt(Miu*Miu)/Miu;
	}
      Miu = Miu - dparam("KINVIS")+dparam("KINVIS-ORIG");
      //     Miu = Diam[iq];   REMOVE      
      nu[i] = Miu;

      if (option("STATAVG"))
	{
	  u2nde[1][i] = u2nde[1][i] + Diam[iq];
	  u2nde[2][i] = u2nde[2][i] + 2*nu[i]*0.5*(dudy[i]+dvdx[i]);
	}

      MinV = min(MinV,Miu);
      MaxV = max(MaxV,Miu);
      AvV  = AvV+Miu;
      iq++;
    }
  AvV=AvV/nq;
  if (dparam("RANGE") != 0.0)
    for (i=0;i<nq;i++)
      if (sqrt((nu[i]-AvV)*(nu[i]-AvV)) > dparam("RANGE")*AvV)
	nu[i] = dparam("RANGE")*AvV*sqrt((nu[i]-AvV)*(nu[i]-AvV))/(nu[i]-AvV);

  if (option("STATAVG"))
    dsvtvp(nq, 1/dparam("KINVIS-ORIG"), nu, 1, u2nde[0], 1, u2nde[0], 1);
    
  dscal(nq,2,dudx,1);
  dvmul(nq,nu,1,dudx,1,tempo[0],1);
  dvadd(nq,dudy,1,dvdx,1,tempo[1],1);
  dvmul(nq,nu,1,tempo[1],1,tempo[1],1);
  dvadd(nq,dvdx,1,dudy,1,tempo[2],1);
  dvmul(nq,nu,1,tempo[2],1,tempo[2],1);
  dscal(nq,2,dvdy,1);
  dvmul(nq,nu,1,dvdy,1,tempo[3],1);
  dvadd(nq,dwdx,1,Uz->base_h,1,dwdx,1);
  dvmul(nq,nu,1,dwdx,1,tempo[4],1);
  dvadd(nq,dwdy,1,Vz->base_h,1,dwdy,1);
  dvmul(nq,nu,1,dwdy,1,tempo[5],1);

  // Move the Uz data for safe keeping into tempo12
  dcopy(nq,Uz->base_h,1,tempo[12],1);

  // tempo13-15 to be kept so that d*dz may be determined
  dvadd(nq,Uz->base_h,1,dwdx,1,tempo[13],1);
  dvmul(nq,tempo[13],1,nu,1,tempo[13],1);
  dvadd(nq,Vz->base_h,1,dwdy,1,tempo[14],1);
  dvmul(nq,tempo[14],1,nu,1,tempo[14],1);
  dvmul(nq,nu,1,Wz->base_h,1,tempo[15],1);
  dscal(nq,2,tempo[15],1);

  dcopy(nq,tempo[0],1,Uz->base_h,1);
  Uz->Grad_h(Uz->base_h,tempo[6],(double *)NULL,(double *)NULL,trip);
  dcopy(nq,tempo[1],1,Uz->base_h,1);
  Uz->Grad_h(Uz->base_h,(double *)NULL,tempo[7],(double *)NULL,trip);
  dcopy(nq,tempo[2],1,Uz->base_h,1);
  Uz->Grad_h(Uz->base_h,tempo[8],(double *)NULL,(double *)NULL,trip);
  dcopy(nq,tempo[3],1,Uz->base_h,1);
  Uz->Grad_h(Uz->base_h,(double *)NULL,tempo[9],(double *)NULL,trip);
  dcopy(nq,tempo[4],1,Uz->base_h,1);
  Uz->Grad_h(Uz->base_h,tempo[10],(double *)NULL,(double *)NULL,trip);
  dcopy(nq,tempo[5],1,Uz->base_h,1);
  Uz->Grad_h(Uz->base_h,(double *)NULL,tempo[11],(double *)NULL,trip);

  dvadd(nq,tempo[6],1,tempo[7],1,vvexp->NLx,1);
  dvadd(nq,tempo[8],1,tempo[9],1,vvexp->NLy,1);
  dvadd(nq,tempo[10],1,tempo[11],1,vvexp->NLz,1);

  // Calculate the z- derivatives for vvexp->NL*
  // Transform z derivative of velocity to physical/quadrature space
  dcopy(nq,tempo[13],1,Uz->base_h,1);  Uz->Trans(Uz, p_to_f);
  // Take z derivatives in fourier/quadrature space
  Uz->Grad(nul,nul,Uz,'z'); Uz->Set_state('p');
  // Transform z derivative of velocity to physical/quadrature space
  Uz->Trans(Uz, f_to_p);  dcopy(nq,Uz->base_h,1,tempo[13],1);
  
  // Transform z derivative of velocity to physical/quadrature space
  dcopy(nq,tempo[14],1,Uz->base_h,1);  Uz->Trans(Uz, p_to_f);
  // Take z derivatives in fourier/quadrature space
  Uz->Grad(nul,nul,Uz,'z'); Uz->Set_state('p');
  // Transform z derivative of velocity to physical/quadrature space
  Uz->Trans(Uz, f_to_p);  dcopy(nq,Uz->base_h,1,tempo[14],1);
  
  // Transform z derivative of velocity to physical/quadrature space
  dcopy(nq,tempo[15],1,Uz->base_h,1);  Uz->Trans(Uz, p_to_f);
  // Take z derivatives in fourier/quadrature space
  Uz->Grad(nul,nul,Uz,'z'); Uz->Set_state('p');
  // Transform z derivative of velocity to physical/quadrature space
  Uz->Trans(Uz, f_to_p);  dcopy(nq,Uz->base_h,1,tempo[15],1);
  
  //Add the tempo13-15 to the vvexp->NL* terms
  dvadd(nq,tempo[13],1,vvexp->NLx,1,vvexp->NLx,1);
  dvadd(nq,tempo[14],1,vvexp->NLy,1,vvexp->NLy,1);
  dvadd(nq,tempo[15],1,vvexp->NLz,1,vvexp->NLz,1);
  
  dcopy(nq,tempo[12],1,Uz->base_h,1);
  ROOT printf("MinV=%e MaxV=%e AvV=%e \n",MinV,MaxV,AvV);
  
  free_dmatrix(tmp,0,0);free(D);    free_dmatrix(tempo,0,0); 
  free(nu);free(nux);free(nuy); 
  
  return;
}
static void Diameter(Domain *omega, double *Diam);

static varv_variables *preprocessLES(Domain *omega, double *Diam){
  varv_variables *vvexp;
  int nel=omega->U->nel;

  Diameter(omega,Diam);
    
  vvexp = (varv_variables*) malloc(sizeof(varv_variables));
  if (option("dealias"))
    {
      vvexp->NLx  = dvector(0,3*nel*QGmax*QGmax*option("NZ")/2);
      vvexp->NLy  = dvector(0,3*nel*QGmax*QGmax*option("NZ")/2);
      vvexp->NLz  = dvector(0,3*nel*QGmax*QGmax*option("NZ")/2);
    }
  else
    {
      vvexp->NLx  = dvector(0,nel*QGmax*QGmax*option("NZ"));
      vvexp->NLy  = dvector(0,nel*QGmax*QGmax*option("NZ"));
      vvexp->NLz  = dvector(0,nel*QGmax*QGmax*option("NZ"));
    }

  return vvexp;
}

static void Diameter(Domain *omega,double *Diam){
  Element_List *U  =  omega->U,  *V  =  omega->V, *W  =  omega->W;
  int    nq =U->htot,nprocs = option("NPROCS"),nztot= option("NZTOT");
  register int i,k,j;
  Coord X;
  double a,b,c,s,AC,AC1,Radius,D;

  X.x = dvector(0,nq); X.y = dvector(0,nq);

  /* 	Calculating sides a , b , c .		*/
  j=0;
  for (k=0; k<U->nel;k++) {
    Element *LU = U->flevels[0]->flist[k];
    LU->coord(&X);
    register double lx, ly;

    if(LU->identify() == Nek_Tri)
      {
	lx = LU->vert[1].x-LU->vert[0].x; ly = LU->vert[1].y-LU->vert[0].y;
	a = sqrt(lx*lx+ly*ly);
	
	lx = LU->vert[2].x-LU->vert[1].x; ly = LU->vert[2].y-LU->vert[1].y;
	b = sqrt(lx*lx+ly*ly);
	
	lx = LU->vert[0].x-LU->vert[2].x; ly = LU->vert[0].y-LU->vert[2].y;
	c = sqrt(lx*lx+ly*ly);
	
	s = 0.5 * ( a + b + c );
	
	AC  = asin(sqrt((LU->vert[1].y-LU->vert[0].y)*
			(LU->vert[1].y-LU->vert[0].y))/a);
	AC1 = asin(sqrt((LU->vert[2].y-LU->vert[0].y)*
			(LU->vert[2].y-LU->vert[0].y))/c);
	
	D = 0.5 * a * c * sin(AC-AC1) * M_PI * M_PI /(iparam("MODES")*   
						      iparam("MODES"));
      }

    if(LU->identify() == Nek_Quad)
      {
	lx = LU->vert[1].x-LU->vert[0].x; ly = LU->vert[1].y-LU->vert[0].y;
	a = sqrt(lx*lx+ly*ly);
	
	lx = LU->vert[2].x-LU->vert[1].x; ly = LU->vert[2].y-LU->vert[1].y;
	b = sqrt(lx*lx+ly*ly);
	
	lx = LU->vert[0].x-LU->vert[2].x; ly = LU->vert[0].y-LU->vert[2].y;
	c = sqrt(lx*lx+ly*ly);
	
	s = 0.5 * ( a + b + c );
	
	AC  = asin(sqrt((LU->vert[1].y-LU->vert[0].y)*
			(LU->vert[1].y-LU->vert[0].y))/a);
	AC1 = asin(sqrt((LU->vert[2].y-LU->vert[0].y)*
			(LU->vert[2].y-LU->vert[0].y))/c);

	D = 0.5 * a * c * sin(AC-AC1);
	
	lx = LU->vert[3].x-LU->vert[2].x; ly = LU->vert[3].y-LU->vert[2].y;
	a = sqrt(lx*lx+ly*ly);
	
	lx = LU->vert[0].x-LU->vert[3].x; ly = LU->vert[0].y-LU->vert[3].y;
	b = sqrt(lx*lx+ly*ly);
	
	lx = LU->vert[2].x-LU->vert[0].x; ly = LU->vert[3].y-LU->vert[0].y;
	c = sqrt(lx*lx+ly*ly);
	
	s = 0.5 * ( a + b + c );
	
	AC  = asin(sqrt((LU->vert[3].y-LU->vert[2].y)*
			(LU->vert[3].y-LU->vert[2].y))/a);
	AC1 = asin(sqrt((LU->vert[0].y-LU->vert[2].y)*
			(LU->vert[0].y-LU->vert[2].y))/c);


	D = (0.5 * a * c * sin(AC-AC1) + D) 
	  * M_PI * M_PI /(iparam("MODES")*iparam("MODES"));
      }

	D = pow(sqrt(D*D) * dparam("LZ") /(nztot),0.333);  
    for (i=0;i<LU->qa*LU->qb;i++)
      {
	Diam[j]=D;
	if (dparam("DNA") != 0.0)
	  {
	    Radius = sqrt(X.y[i]*X.y[i]+X.x[i]*X.x[i])-0.5;
	    if (Radius < 0.0) Radius = 0.0;
	    Radius = (2 / 3.141592654)*
	      atan(2* 0.41* Radius/(dparam("KINVIS-ORIG")*3.141592654))*
	      (1 -exp(- Radius/(dparam("KINVIS-ORIG")*dparam("DNA"))))*
	      (1 -  exp(- Radius /(dparam("KINVIS-ORIG")*dparam("DNA")))) ;
	    if (Radius > 1.0 ) Radius = 1.0;
	    Diam[j] = D * sqrt(Radius);
	  }	
	if (dparam("CNA") != 0.0)
	  {
	    if (sqrt(X.y[i]*X.y[i]) != 0.0 )
	      {
		Radius = 1-sqrt(X.y[i]*X.y[i]);
		if (Radius <= 0.0) Radius  = 0.0;
		Radius = 
		 (2 / 3.141592654)*atan(2* 0.41* Radius/(dparam("KINVIS-ORIG")*
		 3.141592654))*(1 -  exp(- Radius
                 /(dparam("KINVIS-ORIG")*dparam("CNA"))))*
                 (1 -  exp(- Radius /(dparam("KINVIS-ORIG")*dparam("CNA")))) ;
		 if (Radius > 1.0 ) Radius = 1.0;
		Diam[j] = D * sqrt(Radius);
	      }
	  }
	/* Diam[j]= -0.01*X.y[i]*X.y[i]; REMOVE*/
	j=j+1;
      }
  }
  free(X.x); free(X.y);
}

