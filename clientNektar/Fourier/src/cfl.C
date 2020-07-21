/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/
// tcew done
#include <stdio.h>
#include <time.h>
#include "nektarF.h"

static double *escales = (double*)0;

void   cfl_setup(Element_List *U, double dt){
  Element      *E;
  double *wa, *wb, *za, *zb;
  double dxa, dxb, sca, scb, scc, scd;
  
  escales = dvector(0, U->nel-1);
  
  dzero(U->nel, escales, 1);
  
  for(E=U->fhead;E;E=E->next){
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
      
      escales[E->id] = dt/min(dxa*sqrt(sca), dxb*sqrt(min(scb, scc)));
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
      
      escales[E->id] = dt/min(dxa*sqrt(min(sca,scc)), dxb*sqrt(min(scb, scd)));
    }
  }
}

double cfl_checker(Domain *omega, double dt){

  Element      *E, *F, *G;
  double cfl, umax;
  double *tmp = dvector(0, QGmax*QGmax-1);

  if(!escales)
    cfl_setup(omega->U,dt);
  
  cfl = 0.;

  for(E=omega->U->fhead,F=omega->V->fhead, G=omega->W->fhead;E;E=E->next,F=F->next, G=G->next){

    dvmul (E->qtot, E->h[0], 1, E->h[0], 1, tmp, 1);
    dvvtvp(E->qtot, F->h[0], 1, F->h[0], 1, tmp, 1, tmp, 1);
    dvvtvp(E->qtot, G->h[0], 1, G->h[0], 1, tmp, 1, tmp, 1);

    umax = sqrt(tmp[idmax(E->qtot, tmp, 1)])*escales[E->id];
    
    cfl = max(cfl,umax);
  }
  free(tmp);

  return cfl;
  
}

    

