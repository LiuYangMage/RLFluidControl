
#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include "nektarF.h"

#ifdef WALL_MODEL

double Newton(double c1, double c2, double y);

void SetWallBC(Domain *omega, Bndry **Ubc, int update)
{
  double  *u   =  omega->uk[1],     *v   =  omega->vk[1],     *w   = omega->wk[1],// velocity value at time step n
          *p   =  omega->pk[1];

  register int i,j,k;
  int     qt;
  double  *f;
  Coord   X,Y,XX;
  Bndry   *Bc;
  Element *eU, *eV, *eW, *eP;
	double  *za,*wa,*zb,*wb;

  int     htot = omega->U->nz*omega->U->htot;
  int     work_face;
  
  double dz = dparam("LZ")/option("NZTOT");

  double  **workspace = dmatrix(0,2,0,QGmax*QGmax-1);
  dzero(3*QGmax*QGmax,workspace[0],1);
  

  XX.x = dvector(0,QGmax*QGmax-1);
  XX.y = dvector(0,QGmax*QGmax-1);

  double  *ux = workspace[0], *uy = workspace[1], *uz = workspace[2];
  
  Element_List  *V[4];
  
  V[0] = omega->Ut;
  V[1] = omega->Vt;
  V[2] = omega->Wt;

  double niu = dparam("KINVIS");
  double karman_c = 0.41;
  double A = 17.;

  dcopy(htot, u, 1, V[0]->base_h, 1); //u
  dcopy(htot, v, 1, V[1]->base_h, 1); //v
  dcopy(htot, w, 1, V[2]->base_h, 1); //w

  Bc = Ubc[0];
  switch(Bc->type){
  case 'S':case 's': 
  }




  for(k = 0; k < V[0]->nz; ++k) 
    if(parid(k)!=1)
    {
     Bc = Ubc[k];

    int eid = Bc->elmt->id;
    int face = Bc->face;
  
  if(Bc->face == 0 || Bc->face ==2)
     qt = V[0]->flevels[k]->flist[eid]->qa;
  else 
     qt = V[0]->flevels[k]->flist[eid]->qb;
  
    X.x = dvector(0,qt-1);
    X.y = dvector(0,qt-1);

    Y.x = dvector(0,qt-1);
    Y.y = dvector(0,qt-1);
  
    f   = dvector(0,qt-1);
 
  switch(Bc->type){
  case 'S':case 's': 

//    if(update) //need to update slip velocity
     {
	     eU = V[update]->flevels[k]->flist[eid];
	     eU->coord(&XX);

       switch(Bc->face)
        {
          case 0:
          work_face = 2;
          break;
          case 1:
          work_face = 3;
          break;
          case 2:
          work_face = 0;
          break;
          case 3:
          work_face = 1;
          break;
           default:
          fprintf(stderr,"wrong element ! \n");
          exit(1);
          break;
        }
       
       eU->GetFaceCoord(work_face,&X); //opsite of the wall
       eU->GetFaceCoord(face,&Y); //wall
	     
       eU->GetFace(eU->h[0],face,ux);
	     eU->InterpToFace1(work_face,ux,uy); //uy has the aixal velocity
       
       getzw(eU->qa, &za, &wa, 'a');
       getzw(eU->qb, &zb, &wb, 'a');

        double um = 0., hm = 0.; //mean velocity and distance
       for(i=0; i<qt; ++i)
        {
         um += uy[i];
        
          double hh = dparam("RADIUS")-sqrt(X.x[i]*X.x[i]+X.y[i]*X.y[i]);
        if(hh<1e-10)
         {
           fprintf(stderr,"wrong value of distance to pipe wall...");
           exit(-1);
         }

        hm += hh;
        }
       um /= (double)qt;
       hm /= (double)qt;
     
       int m_hat = (int) omega->visc_max_indices[k*V[0]->nel+eid];
       
		   double u_hat = eU->h[0][m_hat];
       double h_hat = dparam("RADIUS")-sqrt(XX.x[m_hat]*XX.x[m_hat]+XX.y[m_hat]*XX.y[m_hat]);

       double ev = omega->visc_ave[eid*eU->qtot+m_hat] ; //entropy viscosity is constant within a elment!
       
       double k_hat = sqrt(ev/(u_hat*h_hat));

       double dxa = 0.5*(za[1]-za[0]);
       double dxb = 0.5*(zb[1]-zb[0]);

       double sca = (eU->vert[1].x-eU->vert[0].x)*(eU->vert[1].x-eU->vert[0].x)+
	            (eU->vert[1].y-eU->vert[0].y)*(eU->vert[1].y-eU->vert[0].y);
      
       double scb = (eU->vert[2].x-eU->vert[1].x)*(eU->vert[2].x-eU->vert[1].x)+
	             (eU->vert[2].y-eU->vert[1].y)*(eU->vert[2].y-eU->vert[1].y);
      
       double scc = (eU->vert[3].x-eU->vert[2].x)*(eU->vert[3].x-eU->vert[2].x)+
              (eU->vert[3].y-eU->vert[2].y)*(eU->vert[3].y-eU->vert[2].y);
      
        double scd = (eU->vert[3].x-eU->vert[0].x)*(eU->vert[3].x-eU->vert[0].x)+
	            (eU->vert[3].y-eU->vert[0].y)*(eU->vert[3].y-eU->vert[0].y);
      
       double Dx = min(sqrt(min(sca,scc)), sqrt(min(scb, scd)));
       double dx = min(dxa*sqrt(min(sca,scc)), dxb*sqrt(min(scb, scd)));
        
       double h_cr = 0.48*min(dz,dx);

       double denom = hm - h_cr;

       double factor = 0.;

       if(fabs(denom)>1e-8)
          factor = min((hm-h_hat)/denom,1.);
       
       double Ck = karman_c*(1.-factor)+k_hat*factor;
      
       double c1 = niu*um/hm;
       double c2 = Ck*um;
       double c3 = hm/A/niu;

//       fprintf(stdout,"c1 = %g  c2 = %g c3 = %g \n",c1,c2,c3);
//       double u_tau = Newton(c1,c2,c3);

       for(i=0; i<qt; ++i)
        {
         double tmp = dparam("RADIUS")-sqrt(X.x[i]*X.x[i]+X.y[i]*X.y[i]);
//         f[i] = uy[i]*tmp*u_tau*u_tau;  //slip velocity
//         f[i] = uy[i]*tmp*dparam("RE_TAU")*dparam("RE_TAU")*niu*niu;  //slip velocity
         f[i] = 0.;//um-tmp*dparam("RE_TAU")*dparam("RE_TAU")*niu*niu;  //slip velocity
  //       fprintf(stdout,"slip velocity =  %g \n",f[i]);
        }
    Bc->bvert[0] = f[0];
    eU->JtransEdge(Bc,Bc->face,0,f);

     }

    break;
    default: 
    break;
    }

    free(X.x);  free(X.y);  free(f);
    free(Y.x);  free(Y.y);
   }
    free(workspace);
    free(XX.x);free(XX.y);
}



double f_stress(double c1, double c2, double c3, double x);
double grad_stress(double c2, double c3, double x);

double Newton(double c1, double c2, double c3)
 {
  double abstol = 1.0e-10;
  int maxit = 200;
  double x0=dparam("KINVIS");
  double xnew;
 
  /* Newton's Method */
  for(int i=0;i<maxit;i++){
    xnew = x0 - f_stress(c1,c2,c3,x0)/grad_stress(c2,c3,x0); 
    if( fabs(f_stress(c1,c2,c3,xnew))<abstol)
      return(xnew);
    fprintf(stderr,"Newton root: %-25.6g Error: %-25.6g\n",xnew,fabs(f_stress(c1,c2,c3,xnew)));
    x0 = xnew;
  }

  /* Bisection Method */
  int counter = 0;
  double xleft  = xnew - 1.e-4;
  double xright = xnew + 1.e-4;
  while(fabs(f_stress(c1,c2,c3,xnew))>abstol){
    fprintf(stderr,"Bisection root: %-25.6g Error: %-25.6g\n",xnew,fabs(f_stress(c1,c2,c3,xnew)));
    if( f_stress(c1,c2,c3,xnew) * f_stress(c1,c2,c3,xleft) > 0.0e0){
      xleft = xnew;
      xnew = 0.5*(xleft+xright);
    }
    else{
      xright = xnew;
      xnew = 0.5*(xleft+xright);
    }
    counter++;
    if(counter>maxit)
      break;
  }

  fprintf(stderr, "Newton Iteration with Bisection did not converge: root: %-25.6g Error: %-25.6g \n",
	  xnew,fabs(f_stress(c1,c2,c3,xnew)));
  return(xnew);
}
 
 double f_stress(double c1, double c2, double c3, double x)
 {
    return c1+c2*x*pow(1.-exp(c3*x),2.)-x*x;
 }

 double grad_stress(double c2, double c3, double x)
 {
    return -2.*x+c2*(1.-exp(c3*x))*(1.-exp(c3*x)-2*c3*x*exp(c3*x));
 }

#endif
