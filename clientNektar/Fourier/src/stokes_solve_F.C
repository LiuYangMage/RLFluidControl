#include "nektarF.h"
#include "stokes_solve_F.h"
#include "BlkMat.h"
#include "SClevel.h"

#include <cassert>
#include <vector>

void StokesMatrix::GenMat(Element_List *U, Element_List *P, 
			  Bsystem *Ubsys, Bsystem *Pbsys, 
			  Metric *lambda, double beta ){
  
  LocMat    *helm;
  LocMat    *pse1, *pse2, *pse3, *pse4;      // MSB: Matrix for PSE terms
  LocMatDiv *div;
  Element *E;
  int eDIM = U->flist[0]->dim();
  int qt;
  
  if(!Ubsys->signchange){
    setup_signchange(U,Ubsys);
  }
  
  SClev = new SClevel;

  if(Ubsys->lambda->wave)                    // MSB: OSEEN Formulation
    SClev->Setup_BlkMat(U->nel,0);           // MSB: Full Matrix
  else                                       // MSB: STOKES Formulation
    SClev->Setup_BlkMat(U->nel,1);           // MSB: Symm Matrix       
    
  // set up matrices 
  for(E = U->fhead; E; E=E->next){ 
    helm = E->mat_mem ();
    div  = E->divmat_mem (P->flist[E->id]);

    // Generate Divergence Matrices 
    E->DivMat(div,P->flist[E->id]);
    
    // Generate Helmholtz/Laplacian matrix
    E->HelmMatC(helm,lambda+E->id);

    qt = E->qa*E->qb;    
#ifdef PSE
    double *DU;
    DU = dvector(0, qt - 1);

    if(Ubsys->lambda->wave){    
      pse1 = E->mat_mem ();        // MSB: Memory allocation of PSE terms
      pse2 = E->mat_mem ();
      pse3 = E->mat_mem ();
      pse4 = E->mat_mem ();
    }

    // MSB: Generate extra PSE terms matrices
    // MSB: Note that scaling by KINVIS is required because of the
    // MSB: way in which lambda->wave was scaled in Set_Oseen

    if(Ubsys->lambda->wave){   
      E->Grad_h((lambda+E->id)->wave[0], DU, NULL, NULL, 'x'); 
      dscal(qt, Ubsys->lambda[E->id].d, DU, 1); // MSB: Scale by KINVIS   
      E->PSE_Mat(E, pse1, DU);                  // MSB: u.dU/dx

      E->Grad_h((lambda+E->id)->wave[0], NULL, DU, NULL, 'y');  
      dscal(qt, Ubsys->lambda[E->id].d, DU, 1); // MSB: Scale by KINVIS   
      E->PSE_Mat(E, pse2, DU);                  // MSB: v.dU/dy

      E->Grad_h((lambda+E->id)->wave[1], DU, NULL, NULL, 'x');    
      dscal(qt, Ubsys->lambda[E->id].d, DU, 1); // MSB: Scale by KINVIS   
      E->PSE_Mat(E, pse3, DU);                  // MSB: u.dV/dx

      E->Grad_h((lambda+E->id)->wave[1], NULL, DU, NULL, 'y');    
      dscal(qt, Ubsys->lambda[E->id].d, DU, 1); // MSB: Scale by KINVIS   
      E->PSE_Mat(E, pse4, DU);                  // MSB: v.dV/dy
    }

    free(DU);
#endif

    // Multiply by viscosity 
    dscal(helm->asize*helm->asize,Ubsys->lambda[E->id].d,*helm->a,1);
    if(helm->csize){
      dscal(helm->asize*helm->csize,Ubsys->lambda[E->id].d,*helm->b,1);
      dscal(helm->csize*helm->csize,Ubsys->lambda[E->id].d,*helm->c,1);
      
      if(Ubsys->lambda->wave)                   // MSB: OSEEN Formulation
      	dscal(helm->asize*helm->csize,Ubsys->lambda[E->id].d,*helm->d,1); 
    }

    // MSB: PROJECT L0 MATRIX TERMS-----------------------
    Project(pse1,pse2,pse3,pse4,helm,div,E->id,Ubsys);
    // MSB: ----------------------------------------------  
    LocMatDiv *bet;
    //double beta   = dparam("BETA");             // MSB: BETA => bar(BETA)
    double sigma  = dparam("SIGMA");            // MSB: IMAGINARY cmpnt
    double *Beta  = dvector(0, qt - 1);
    double *Sigma = dvector(0, qt - 1);

#ifdef PSE_SLV // MSB: For PSE march in z-direction-------
    double idz = 1.0/dparam("DT");
    dfill(qt,beta,Beta,1);
    dfill(qt,sigma-idz,Sigma,1);   
#else // MSB: --------------------------------------------
    dfill(qt,beta,Beta,1);
    dfill(qt,sigma,Sigma,1);    
#endif

    bet  = E->divmat_mem (P->flist[E->id]);
    E->BET_Mat(P->flist[E->id], bet, Beta, Sigma);

    Project_Beta(helm, bet, E->id, Ubsys);

#ifdef COUPLED_TERMS
    LocMat    *amat, *bmat;

    double kinvis = Ubsys->lambda[E->id].d;
    double Omega  = dparam("OMEGA");            // MSB: Omega

    double *A     = dvector(0, qt - 1);
    double *B     = dvector(0, qt - 1);
    double *DWx   = dvector(0, qt - 1);
    double *DWy   = dvector(0, qt - 1);

    // MSB: Scaling of Wbar by KINVIS is due to scaling in Set_Oseen routine
    for(int i = 0; i < qt; i++){
      // A[i]      = -kinvis*sigma*sigma;
      // A[i]     -= sigma*kinvis*(lambda+E->id)->wave[2][i]; 
      A[i]      = -sigma*kinvis*(lambda+E->id)->wave[2][i];
      B[i]      = -2.0*(kinvis*beta*sigma);           
      B[i]     -= beta*kinvis*(lambda+E->id)->wave[2][i];  // MSB: -beta*Wbar 
      B[i]     += Omega;                                   // MSB: +Omega
      
#ifdef PSE_SLV
      A[i]     += idz*(kinvis*(lambda+E->id)->wave[2][i]+2.0*kinvis*sigma);
      B[i]     += idz*(2.0*kinvis*beta);
#endif
    }

    amat = E->mat_mem ();
    bmat = E->mat_mem ();
    
    E->PSE_Mat(E, amat, A);    
    E->PSE_Mat(E, bmat, B);    

    // MSB: PROJECT SOME OF THE L0 COUPLED TERMS----------
    // MSB: IF PSE_SLV THEN THIS INCLUDES L2/DZ-----------
    Project_CT(amat, bmat, helm, bet->rows-1, E->id, Ubsys);
    // MSB: ----------------------------------------------
    
    E->Grad_h((lambda+E->id)->wave[2], DWx, NULL, NULL, 'x'); 
    E->Grad_h((lambda+E->id)->wave[2], NULL, DWy, NULL, 'y'); 

    // MSB: Scale by KINVIS because of Set_Oseen 
    dscal(qt, kinvis, DWx, 1);
    dscal(qt, kinvis, DWy, 1);  
    E->PSE_Mat(E, amat, DWx);  
    E->PSE_Mat(E, bmat, DWy);    

    // MSB: PROJECT MORE OF THE L0 COUPLED TERMS----------
    Project_CT2(amat, bmat, bet, helm, E->id, Ubsys);
    // MSB: ----------------------------------------------

    dzero(qt, DWx, 1);     dzero(qt, A, 1); 
    dzero(qt, DWy, 1);     dzero(qt, B, 1);
    dzero(qt, Beta, 1);    dzero(qt, Sigma, 1);
    
#ifdef PSE_SLV
    LocMat  *cmat;
    cmat = E->mat_mem ();

    // MSB: Add a routine to find these terms.....
    // MSB: DWx => dU/dz, DWy => dV/dz, A => dW/dz
    // MSB: Add a routine to find these terms.....
    
    E->PSE_Mat(E, amat, DWx);  
    E->PSE_Mat(E, bmat, DWy);
    E->PSE_Mat(E, cmat, A);   

    // MSB: PROJECT THE L1 MATRIX-------------------------
    Project_L1(amat, bmat, cmat, helm, bet, E->id, Ubsys);
    // MSB: ----------------------------------------------

    E->mat_free   (cmat);
#endif    
    free(A), free(B), free(Beta), free(Sigma), free (DWx), free(DWy);

    E->mat_free   (amat);
    E->mat_free   (bmat);    
    E->divmat_free (bet);
#endif

    E->mat_free (helm);
    E->divmat_free(div);
  }
  
  Set_SignChange(U,Ubsys);
  SClev->Set_State(Filled);

 
  SClev->Condense();
  

  // need to set up assembly routine to interior matrix. 
  Ainv = new Gmat; 
     
  if(Pbsys->singular == 0){
    Set_locmap(Pbsys->bmap[0],Pbsys->nsolve,U);
    Ainv->set_lda(2*Pbsys->nsolve);             // MSB: Note change in lda
  }
  else{
    Set_locmap(Pbsys->bmap[0],Pbsys->nsolve-1,U);  
    Ainv->set_lda(2*(Pbsys->nsolve-1));         // MSB: Note change in lda
  }  
  
  Ainv->set_bwidth(SClev->bandwidth()); 

                                 
  if(Ubsys->lambda->wave)                       // MSB: OSEEN Formulation  
    Ainv->mat_form = General_Full;              // MSB: Full Matrix Solver
  else{                                         // MSB: STOKES Formulation
    if(3*Ainv->get_bwidth() > Ainv->get_lda())  // MSB: Symm Matrix Solver
      Ainv->mat_form = Symmetric;      
    else
      Ainv->mat_form = General_Banded; 
  }

  //assemble into global matrix system taking account of signchange if nec.
  Ainv->Assemble(SClev->_Am[0],SClev->_locmap,SClev->_signchange);
  
  // release memory associated with row matrix storage but set up offset !
  SClev->_Am->setup_offset();
  SClev->_Am->Reset_Rows();

  Ainv->Factor(); 

}
void StokesMatrix::Set_locmap(int *map,int nsolve, Element_List *U){

  SClev->_locmap = new int[SClev->_Asize];
  Element *E;
  
  int i;
  int *temp = new int[SClev->_Asize/2];
  vmath::vcopy(SClev->_Asize/2,map,1,temp,1);

  // reset values above nsolve to -1 for SC solver
  for(i=0; i < SClev->_Asize/2; ++i)
    if(temp[i] >= nsolve)
      temp[i] = -1;
  /*
  for(i=0; i < SClev->_Asize/2; ++i){
    SClev->_locmap[2*i]   = 2*temp[i];
    SClev->_locmap[2*i+1] = 2*temp[i] + 1;
  }
  */
  int cnt = 0, step = 0, delta = 0;
  int   j = 0,    k = 0;

  for(E=U->fhead; E; E=E->next, k++){             // Number Elements
    SClev->_locmap[step]   = 2*temp[cnt];         // P_real
    SClev->_locmap[step+1] = 2*temp[cnt++] + 1;   // P_imag
    for(i=0; i < 3; i++){                         // U, V, W => 3
      for(j=0; j < E->Nbmodes; j++){              // E->Nbmodes 
	SClev->_locmap[step+delta+j+2]            = 2*temp[cnt];
	SClev->_locmap[step+delta+E->Nbmodes+j+2] = 2*temp[cnt++] + 1;      
      }
      delta += 2*E->Nbmodes;
    }
    delta = 0;
    step += 2*(3*E->Nbmodes+1);                   // 3*(E->Nbmodes+1)
  }

  for(i=0; i < SClev->_Asize; ++i)
    if(SClev->_locmap[i] < -1)
      SClev->_locmap[i] = -1;
  
  /*

  int cnt = 0, step = 0, delta = 0;
  int   j = 0,    k = 0;

  for(k=0; k < 4; k++){                           // Number Elements
    SClev->_locmap[step]   = 2*temp[cnt];         // P_real
    SClev->_locmap[step+1] = 2*temp[cnt++] + 1;   // P_imag
    for(i=0; i < 3; i++){                         // U, V, W => 3
      for(j=0; j < 15; j++){                      // E->Nbmodes 
	SClev->_locmap[step+delta+j+2]    = 2*temp[cnt];
	SClev->_locmap[step+delta+15+j+2] = 2*temp[cnt++] + 1;      
      }
      delta += 2*15;
    }
    delta = 0;
    step += 2*(45+1);                             // 3*(E->Nbmodes+1)
  }
  */
   delete [] temp;
}

void StokesMatrix::Project(LocMat *pse1, LocMat *pse2, LocMat *pse3,
			   LocMat *pse4, LocMat *helm, LocMatDiv *div, 
			   int eid, Bsystem *Ubsys){
  int i,j,k,l,nr,nc,cnt;
  int psize = div->rows-1;   
             
  int Asize = 2*(3*helm->asize+1);
  int Csize = 2*(3*helm->csize+psize);    
  
  if((helm->asize != div->bsize)||(helm->csize != div->isize))
    NekError::error(fatal,"StokesMatrix::Project", 
		    "Matrices do not match in size");
  

  // FILL HELMHOLTZ SYSTEM-----------------------------------------------------

  //------------------------------------ 
  // fill Helm boundary boundary system   
  //------------------------------------ 

  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+k,
			 2+l,helm->a[k][l]);
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+2*helm->asize+k,
			 2+2*helm->asize+l,helm->a[k][l]);  
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
      			 2+4*helm->asize+l,helm->a[k][l]);  
    }

  // IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+helm->asize+k,
			 2+helm->asize+l,helm->a[k][l]);   
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+3*helm->asize+k,
			 2+3*helm->asize+l,helm->a[k][l]);  
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
      			 2+5*helm->asize+l,helm->a[k][l]);  
    }
  
  //------------------------------------ 
  // fill Helm boundary interior system   
  //------------------------------------  

  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+k,
			 2*psize+l,helm->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			 2*psize+2*helm->csize+l,helm->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
      			 2*psize+4*helm->csize+l,helm->b[k][l]);      
    }

  // IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			 2*psize+helm->csize+l,helm->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			 2*psize+3*helm->csize+l,helm->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
      			 2*psize+5*helm->csize+l,helm->b[k][l]);      
    }
  
  //------------------------------------ 
  // fill Helm interior interior system     
  //------------------------------------    
  
  // REAL
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+k,
			 2*psize+l,helm->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+2*helm->csize+k,
			 2*psize+2*helm->csize+l,helm->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+k,
			 2*psize+4*helm->csize+l,helm->c[k][l]);
    }  
  
  // IMAGINARY
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+helm->csize+k,
			 2*psize+helm->csize+l,helm->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+3*helm->csize+k,
			 2*psize+3*helm->csize+l,helm->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+k,
			 2*psize+5*helm->csize+l,helm->c[k][l]);
    }  
  
  // HELMHOLTZ SYSTEM FILLED---------------------------------------------------
  
  // FILL Dx,Dy SYSTEM---------------------------------------------------------

  //------------------------------------   
  // fill Dx,Dy boundary boundary system    
  //------------------------------------ 

  // REAL
  for(k = 0; k < div->bsize; ++k){
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,0,2+k,div->Dxb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,0,2+2*helm->asize+k,div->Dyb[0][k]);
    // transpose terms    
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+k,0,div->Dxb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+2*helm->asize+k,0,div->Dyb[0][k]);
  }
  
  // IMAGINARY
  for(k = 0; k < div->bsize; ++k){
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,1,2+helm->asize+k,div->Dxb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,1,2+3*helm->asize+k,div->Dyb[0][k]);
    // transpose terms    
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+helm->asize+k,1,div->Dxb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+3*helm->asize+k,1,div->Dyb[0][k]);
  } 
 
  //------------------------------------   
  // fill Dx,Dy boundary interior system 
  //------------------------------------   

  // REAL
  for(k = 0; k < div->isize; ++k){
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,0,2*psize+k,
		       div->Dxi[0][k]);
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,0,2*psize+2*helm->csize+k,
		       div->Dyi[0][k]);
  }
  // transpose terms
  for(k = 0; k < div->bsize; ++k)
    for(l = 0; l < psize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+k,l,
			 div->Dxb[l+1][k]); 
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,l,
			 div->Dyb[l+1][k]);  
    } 

  // IMAGINARY
  for(k = 0; k < div->isize; ++k){
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,1,2*psize+helm->csize+k,
		       div->Dxi[0][k]);
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,1,2*psize+3*helm->csize+k,
		       div->Dyi[0][k]);
  }
  // transpose terms
  for(k = 0; k < div->bsize; ++k)
    for(l = 0; l < psize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,l+psize,
			 div->Dxb[l+1][k]); 
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,l+psize,
			 div->Dyb[l+1][k]);  
    } 

  //------------------------------------   
  // fill Dx,Dy interior interior system 
  //------------------------------------     

  // REAL
  for(k = 0; k < psize; ++k)  
    for(l = 0; l < div->isize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,k,2*psize+l,
			 div->Dxi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,k,2*psize+2*helm->csize+l,
			 div->Dyi[k+1][l]);
      // transpose terms
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+l,k,
			 div->Dxi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+2*helm->csize+l,k,
			 div->Dyi[k+1][l]);
    }  

  // IMAGINARY

  for(k = 0; k < psize; ++k)  
    for(l = 0; l < div->isize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,psize+k,2*psize+helm->csize+l,
			 div->Dxi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,psize+k,2*psize+3*helm->csize+l,
			 div->Dyi[k+1][l]);
      // transpose terms
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+helm->csize+l,psize+k,
			 div->Dxi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+3*helm->csize+l,psize+k,
			 div->Dyi[k+1][l]);
    }  
  
  // Dx,Dy SYSTEM FILLED-------------------------------------------------------

  // ==========================================================================
  // MSB: OSEEN FORMULATION  
  // ==========================================================================
  if(Ubsys->lambda->wave){    // MSB: Fill _Dm because problem is not symmetric

    //------------------------------------ 
    // fill Helm interior boundary system     
    //------------------------------------  
  
    // REAL
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+k,
			   2*psize+l,helm->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			   2*psize+2*helm->csize+l,helm->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			   2*psize+4*helm->csize+l,helm->d[k][l]);      
      }
    
    // IMAGINARY
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			   2*psize+helm->csize+l,helm->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			   2*psize+3*helm->csize+l,helm->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			   2*psize+5*helm->csize+l,helm->d[k][l]);      
      }
 
    //------------------------------------   
    // fill Dx,Dy interior boundary system 
    //------------------------------------   

    // REAL
    for(k = 0; k < div->isize; ++k){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,0,2*psize+k,
			 div->Dxi[0][k]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,0,2*psize+2*helm->csize+k,
			 div->Dyi[0][k]);
    }
    // transpose terms
    for(k = 0; k < div->bsize; ++k)
      for(l = 0; l < psize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+k,l,
			   div->Dxb[l+1][k]); 
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,l,
			   div->Dyb[l+1][k]);  
      } 

    // IMAGINARY
    for(k = 0; k < div->isize; ++k){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,1,2*psize+helm->csize+k,
			 div->Dxi[0][k]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,1,2*psize+3*helm->csize+k,
			 div->Dyi[0][k]);
    }
    // transpose terms
    for(k = 0; k < div->bsize; ++k)
      for(l = 0; l < psize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,l+psize,
			   div->Dxb[l+1][k]); 
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,l+psize,
			   div->Dyb[l+1][k]); 
      } 
  }
  // ==========================================================================
  // MSB: OSEEN TERMS FILLED 
  // ==========================================================================

  // ==========================================================================
  // MSB: PSE FORMULATION  
  // pse1 => u.dU/dx
  // pse2 => v.dU/dy
  // pse3 => u.dV/dx
  // pse4 => v.dV/dy
  // ==========================================================================
#ifdef PSE
  // FILL u.grad(U) SYSTEM ----------------------------------------------------
  if(Ubsys->lambda->wave){
    
    //----------------------------------------- 
    // fill u.grad(U) boundary boundary system   
    //----------------------------------------- 
    
    // REAL
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->asize; ++l){
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+k,
			   2+l,pse1->a[k][l]);
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+2*helm->asize+k,
			   2+l,pse3->a[k][l]);
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+k,
			   2+2*helm->asize+l,pse2->a[k][l]);
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+2*helm->asize+k,
			   2+2*helm->asize+l,pse4->a[k][l]);  
      }

    // IMAGINARY
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->asize; ++l){
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+helm->asize+k,
			   2+helm->asize+l,pse1->a[k][l]);
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+3*helm->asize+k,
			   2+helm->asize+l,pse3->a[k][l]);
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+helm->asize+k,
			   2+3*helm->asize+l,pse2->a[k][l]);
	SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+3*helm->asize+k,
			   2+3*helm->asize+l,pse4->a[k][l]);  
      }    

    //----------------------------------------- 
    // fill u.grad(U) interior boundary system   
    //----------------------------------------- 

    // REAL
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+k,
			   2*psize+l,pse1->b[k][l]);
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			   2*psize+l,pse3->b[k][l]);	
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+k,
			   2*psize+2*helm->csize+l,pse2->b[k][l]);
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			   2*psize+2*helm->csize+l,pse4->b[k][l]);
      }

    // IMAGINARY
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			   2*psize+helm->csize+l,pse1->b[k][l]);
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			   2*psize+helm->csize+l,pse3->b[k][l]);	
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			   2*psize+3*helm->csize+l,pse2->b[k][l]);
	SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			   2*psize+3*helm->csize+l,pse4->b[k][l]);
      }

    //----------------------------------------- 
    // fill u.grad(U) interior interior system   
    //----------------------------------------- 

    // REAL
    for(k = 0; k < helm->csize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+k,
			   2*psize+l,pse1->c[k][l]);
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+2*helm->csize+k,
			   2*psize+l,pse3->c[k][l]);
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+k,
			   2*psize+2*helm->csize+l,pse2->c[k][l]);
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+2*helm->csize+k,
			   2*psize+2*helm->csize+l,pse4->c[k][l]);
      }

    // IMAGINARY
    for(k = 0; k < helm->csize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+helm->csize+k,
			   2*psize+helm->csize+l,pse1->c[k][l]);
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+3*helm->csize+k,
			   2*psize+helm->csize+l,pse3->c[k][l]);
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+helm->csize+k,
			   2*psize+3*helm->csize+l,pse2->c[k][l]);
	SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+3*helm->csize+k,
			   2*psize+3*helm->csize+l,pse4->c[k][l]);
      }

    //----------------------------------------- 
    // fill u.grad(U) boundary interior system   
    //----------------------------------------- 

    // REAL
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+k,
			   2*psize+l,pse1->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+k,
			   2*psize+2*helm->csize+l,pse3->d[k][l]);	
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			   2*psize+l,pse2->d[k][l]);	
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			   2*psize+2*helm->csize+l,pse4->d[k][l]);
      }

    // IMAGINARY
    for(k = 0; k < helm->asize; ++k)  
      for(l = 0; l < helm->csize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			   2*psize+helm->csize+l,pse1->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			   2*psize+3*helm->csize+l,pse3->d[k][l]);	
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			   2*psize+helm->csize+l,pse2->d[k][l]);	
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			   2*psize+3*helm->csize+l,pse4->d[k][l]);
      }
  }
  // u.grad(U) SYSTEM FILLED---------------------------------------------------
#endif
  // ========================================================================= 

  SClev->_Asize += Asize;
  SClev->_Csize += Csize;

}

// ===========================================================================
// PROJECT THE COUPLED TERMS INTO THE MATRIX SYSTEM
// =========================================================================== 
void StokesMatrix::Project_CT(LocMat *amat, LocMat *bmat, LocMat *helm,
			      int psize, int eid, Bsystem *Ubsys){

  
  int i,j,k,l,nr,nc,cnt;
  int Asize = 2*(3*helm->asize+1);
  int Csize = 2*(3*helm->csize+psize);      

  // FILL a,b,c,d TERMS-------------------------------------------------------

  //----------------------------------------- 
  // fill a,b,c,d boundary boundary system   
  //----------------------------------------- 
  
  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+k,
			 2+l,amat->a[k][l]);
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+2*helm->asize+k,
			 2+2*helm->asize+l,amat->a[k][l]);  
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
			 2+4*helm->asize+l,amat->a[k][l]);  

      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+helm->asize+k,
			 2+l,-bmat->a[k][l]);
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+3*helm->asize+k,
			 2+2*helm->asize+l,-bmat->a[k][l]);  
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
			 2+4*helm->asize+l,-bmat->a[k][l]);  

    }

  // IMAGINARY 
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+helm->asize+k,
			 2+helm->asize+l,amat->a[k][l]);   
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+3*helm->asize+k,
			 2+3*helm->asize+l,amat->a[k][l]);  
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
			 2+5*helm->asize+l,amat->a[k][l]);  

      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+k,
			 2+helm->asize+l,bmat->a[k][l]);   
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+2*helm->asize+k,
			 2+3*helm->asize+l,bmat->a[k][l]);  
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
			 2+5*helm->asize+l,bmat->a[k][l]);  
    }
  
  //----------------------------------------- 
  // fill a,b,c,d interior boundary system   
  //----------------------------------------- 
  
  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+k,
			 2*psize+l,amat->b[k][l]); 
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			 2*psize+2*helm->csize+l,amat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 2*psize+4*helm->csize+l,amat->b[k][l]);

      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			 2*psize+l,-bmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			 2*psize+2*helm->csize+l,-bmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 2*psize+4*helm->csize+l,-bmat->b[k][l]);
    }

  // IMAGINARY 
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			 2*psize+helm->csize+l,amat->b[k][l]);   
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			 2*psize+3*helm->csize+l,amat->b[k][l]); 
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 2*psize+5*helm->csize+l,amat->b[k][l]); 

      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+k,
			 2*psize+helm->csize+l,bmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			 2*psize+3*helm->csize+l,bmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 2*psize+5*helm->csize+l,bmat->b[k][l]);
    }

  //----------------------------------------- 
  // fill a,b,c,d interior interior system   
  //----------------------------------------- 
  
  // REAL
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+k,
			 2*psize+l,amat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+2*helm->csize+k,
			 2*psize+2*helm->csize+l,amat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+k,
			 2*psize+4*helm->csize+l,amat->c[k][l]);

      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+helm->csize+k,
			 2*psize+l,-bmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+3*helm->csize+k,
			 2*psize+2*helm->csize+l,-bmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+k,
			 2*psize+4*helm->csize+l,-bmat->c[k][l]);
    }  

  // IMAGINARY 
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+helm->csize+k,
			 2*psize+helm->csize+l,amat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+3*helm->csize+k,
			 2*psize+3*helm->csize+l,amat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+k,
			 2*psize+5*helm->csize+l,amat->c[k][l]);

      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+k,
			 2*psize+helm->csize+l,bmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+2*helm->csize+k,
			 2*psize+3*helm->csize+l,bmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+k,
			 2*psize+5*helm->csize+l,bmat->c[k][l]);
    }
  
  //----------------------------------------- 
  // fill a,b,c,d boundary interior system   
  //----------------------------------------- 
  
  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+k,
			 2*psize+l,amat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			 2*psize+2*helm->csize+l,amat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 2*psize+4*helm->csize+l,amat->d[k][l]);

      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,         // ***
			 2*psize+l,bmat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,       // ***
			 2*psize+2*helm->csize+l,bmat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,       // ***
			 2*psize+4*helm->csize+l,bmat->d[k][l]);
    }
  // *** INDICATES TRANSPOSE c.f. _Bm (CHECK IF THIS IS CORRECT)  
  // IMAGINARY 
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			 2*psize+helm->csize+l,amat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			 2*psize+3*helm->csize+l,amat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 2*psize+5*helm->csize+l,amat->d[k][l]);

      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+k,
			 2*psize+helm->csize+l,-bmat->d[k][l]);       // ***
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			 2*psize+3*helm->csize+l,-bmat->d[k][l]);     // ***
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 2*psize+5*helm->csize+l,-bmat->d[k][l]);     // ***
    }
  // *** INDICATES TRANSPOSE c.f. _Bm (CHECK IF THIS IS CORRECT)
  // a,b,c,d TERMS FILLED-----------------------------------------------------

}


// ===========================================================================
// PROJECT THE BETA and SIGMA  TERMS INTO THE MATRIX SYSTEM
// =========================================================================== 
void StokesMatrix::Project_Beta(LocMat *helm, LocMatDiv *bet, 
			      int eid, Bsystem *Ubsys){

  
  if((helm->asize != bet->bsize)||(helm->csize != bet->isize))
    NekError::error(fatal,"StokesMatrix::Project", 
		    "Matrices do not match in size");
  
  int i,j,k,l,nr,nc,cnt;
  int psize = bet->rows-1;   
  int Asize = 2*(3*helm->asize+1);
  int Csize = 2*(3*helm->csize+psize);      


  // FILL beta,sigma TERMS---------------------------------------------------

  //------------------------------------------   
  // fill beta,sigma boundary boundary system    
  //------------------------------------------ 
  // MSB: Note signs for consistency with div
  // REAL
  for(k = 0; k < bet->bsize; ++k){
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,0,
		       2+4*helm->asize+k,bet->Dyb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,1,
		       2+4*helm->asize+k,-bet->Dxb[0][k]);
    // transpose terms    
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
		       0,-bet->Dyb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
		       0,bet->Dxb[0][k]);
  }
  // MSB: Note signs for consistency with div  
  // IMAGINARY
  for(k = 0; k < bet->bsize; ++k){
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,0,
		       2+5*helm->asize+k,bet->Dxb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,1,
		       2+5*helm->asize+k,bet->Dyb[0][k]);
    // transpose terms    
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
		       1,-bet->Dxb[0][k]);
    SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
		       1,-bet->Dyb[0][k]);
  }

  //------------------------------------------   
  // fill beta,sigma boundary interior system    
  //------------------------------------------ 
  // MSB: Note signs for consistency with div
  // REAL
  for(k = 0; k < bet->isize; ++k){
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,0,
		       2*psize+4*helm->csize+k,bet->Dyi[0][k]);
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,1,
		       2*psize+4*helm->csize+k,-bet->Dxi[0][k]);    
  }
  // transpose terms
  for(k = 0; k < bet->bsize; ++k)
    for(l = 0; l < psize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 l,-bet->Dyb[l+1][k]); 
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 l,bet->Dxb[l+1][k]);  
    } 
  // MSB: Note signs for consistency with div
  // IMAGINARY
  for(k = 0; k < bet->isize; ++k){
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,0,
		       2*psize+5*helm->csize+k,bet->Dxi[0][k]);
    SClev->_Bm->GenBlk(eid,eid,Asize,Csize,1,
		       2*psize+5*helm->csize+k,bet->Dyi[0][k]);    
  }
  // transpose terms
  for(k = 0; k < bet->bsize; ++k)
    for(l = 0; l < psize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 psize+l,-bet->Dxb[l+1][k]); 
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 psize+l,-bet->Dyb[l+1][k]);  
    } 

  //------------------------------------------   
  // fill beta,sigma interior interior system    
  //------------------------------------------ 
  // MSB: Note signs for consistency with div
  // REAL
  for(k = 0; k < psize; ++k)  
    for(l = 0; l < bet->isize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,k,
			 2*psize+4*helm->csize+l,bet->Dyi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,psize+k,
			 2*psize+4*helm->csize+l,-bet->Dxi[k+1][l]);
      // transpose terms
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+l,
			 k,-bet->Dyi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+l,
			 k,bet->Dxi[k+1][l]);
    }  
  // MSB: Note signs for consistency with div
  // IMAGINARY
  for(k = 0; k < psize; ++k)  
    for(l = 0; l < bet->isize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,k,
			 2*psize+5*helm->csize+l,bet->Dxi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,psize+k,
			 2*psize+5*helm->csize+l,bet->Dyi[k+1][l]);
      // transpose terms
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+l,
			 psize+k,-bet->Dxi[k+1][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+l,
			 psize+k,-bet->Dyi[k+1][l]);
    }  

  //------------------------------------------   
  // fill beta,sigma interior boundary system    
  //------------------------------------------ 
  // MSB: Note signs for consistency with div
  // REAL
  for(k = 0; k < bet->isize; ++k){
    SClev->_Dm->GenBlk(eid,eid,Asize,Csize,0,
		       2*psize+4*helm->csize+k,-bet->Dyi[0][k]);
    SClev->_Dm->GenBlk(eid,eid,Asize,Csize,1,                         // ***
		       2*psize+4*helm->csize+k,-bet->Dxi[0][k]); 
  }
  // transpose terms
  for(k = 0; k < bet->bsize; ++k)
    for(l = 0; l < psize; ++l){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 l,bet->Dyb[l+1][k]); 
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,       // ***
			 l,bet->Dxb[l+1][k]); 
    } 
  // *** INDICATES TRANSPOSE c.f. _Bm (CHECK IF THIS IS CORRECT) 
  // MSB: Note signs for consistency with div   
  // IMAGINARY
  for(k = 0; k < bet->isize; ++k){
    SClev->_Dm->GenBlk(eid,eid,Asize,Csize,0,
		       2*psize+5*helm->csize+k,bet->Dxi[0][k]);       // ***
    SClev->_Dm->GenBlk(eid,eid,Asize,Csize,1,
		       2*psize+5*helm->csize+k,-bet->Dyi[0][k]);
  }
  // transpose terms
  for(k = 0; k < bet->bsize; ++k)
    for(l = 0; l < psize; ++l){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,       // ***
			 psize+l,-bet->Dxb[l+1][k]); 
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 psize+l,bet->Dyb[l+1][k]); 
    } 
  // *** INDICATES TRANSPOSE c.f. _Bm (CHECK IF THIS IS CORRECT)  
  // beta, sigma TERMS FILLED--------------------------------------------------
}




// ============================================================================
// PROJECT MORE OF THE COUPLED TERMS INTO THE MATRIX SYSTEM
// ============================================================================
void StokesMatrix::Project_CT2(LocMat *Wxmat, LocMat *Wymat, LocMatDiv *bet,
			       LocMat *helm, int eid, Bsystem *Ubsys){


  int i,j,k,l,nr,nc,cnt;
  int psize = bet->rows-1;   
  int Asize = 2*(3*helm->asize+1);
  int Csize = 2*(3*helm->csize+psize);      

  if((helm->asize != bet->bsize)||(helm->csize != bet->isize))
    NekError::error(fatal,"StokesMatrix::Project", 
		    "Matrices do not match in size");

  // FILL Wx,Wy SYSTEM---------------------------------------------------------

  //------------------------------------ 
  // fill Wx,Wy boundary boundary system   
  //------------------------------------ 

  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
			 2+l,Wxmat->a[k][l]);
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
			 2+2*helm->asize+l,Wymat->a[k][l]);
    }

  // IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
			 2+helm->asize+l,Wxmat->a[k][l]);   
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
			 2+3*helm->asize+l,Wymat->a[k][l]);   
    }
  
  //------------------------------------ 
  // fill Wx,Wy boundary interior system   
  //------------------------------------  

  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 2*psize+l,Wxmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 2*psize+2*helm->csize+l,Wymat->b[k][l]);
    }

  // IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 2*psize+helm->csize+l,Wxmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 2*psize+3*helm->csize+l,Wymat->b[k][l]);
    }

  //------------------------------------ 
  // fill Wx,Wy interior interior system     
  //------------------------------------    
  
  // REAL
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+k,
			 2*psize+l,Wxmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+k,
			 2*psize+2*helm->csize+l,Wymat->c[k][l]);
    }  
  
  // IMAGINARY
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+k,
			 2*psize+helm->csize+l,Wxmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+k,
			 2*psize+3*helm->csize+l,Wymat->c[k][l]);
    }  

  //------------------------------------ 
  // fill Wx,Wy interior boundary system     
  //------------------------------------  
  
  // REAL
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+k,
			 2*psize+4*helm->csize+l,Wxmat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			 2*psize+4*helm->csize+l,Wymat->d[k][l]);
    }
    
  // IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			 2*psize+5*helm->csize+l,Wxmat->d[k][l]);
      SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			 2*psize+5*helm->csize+l,Wymat->d[k][l]);
      }
  
  // Wx,Wy SYSTEM FILLED-------------------------------------------------------
}
// ============================================================================
// PROJECT  L1 INTO THE MATRIX SYSTEM
// ============================================================================
void StokesMatrix::Project_L1(LocMat *amat, LocMat *bmat, LocMat *cmat,
			      LocMat *helm, LocMatDiv *bet, 
			      int eid, Bsystem *Ubsys){

  int i,j,k,l,nr,nc,cnt;
  int psize = bet->rows-1;   
  int Asize = 2*(3*helm->asize+1);
  int Csize = 2*(3*helm->csize+psize);      

  if((helm->asize != bet->bsize)||(helm->csize != bet->isize))
    NekError::error(fatal,"StokesMatrix::Project", 
		    "Matrices do not match in size");
  
  // FILL Uz SYSTEM------------------------------------------------------------
  
  //------------------------------------ 
  // fill Uz boundary boundary system   
  //------------------------------------ 
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+k,
			 2+4*helm->asize+l,amat->a[k][l]);
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+helm->asize+k,
			 2+5*helm->asize+l,amat->a[k][l]);
    }
  
  //------------------------------------ 
  // fill Uz boundary interior system   
  //------------------------------------  

  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+k,
			 2*psize+4*helm->csize+l,amat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+helm->asize+k,
			 2*psize+5*helm->csize+l,amat->b[k][l]);
    }

  //------------------------------------ 
  // fill Uz interior interior system     
  //------------------------------------    
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+k,
			 2*psize+4*helm->csize+l,amat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+helm->csize+k,
			 2*psize+5*helm->csize+l,amat->c[k][l]);
    }  
  
  //------------------------------------ 
  // fill Uz interior boundary system     
  //------------------------------------  
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			   2*psize+l,amat->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			   2*psize+helm->csize+l,amat->d[k][l]);
    }

  // Uz SYSTEM FILLED----------------------------------------------------------

  // FILL Vz SYSTEM------------------------------------------------------------
  
  //------------------------------------ 
  // fill Vz boundary boundary system   
  //------------------------------------ 
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+2*helm->asize+k,
			 2+4*helm->asize+l,bmat->a[k][l]);
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+3*helm->asize+k,
			 2+5*helm->asize+l,bmat->a[k][l]);
    }
  
  //------------------------------------ 
  // fill Vz boundary interior system   
  //------------------------------------  

  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+2*helm->asize+k,
			 2*psize+4*helm->csize+l,bmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+3*helm->asize+k,
			 2*psize+5*helm->csize+l,bmat->b[k][l]);
    }

  //------------------------------------ 
  // fill Vz interior interior system     
  //------------------------------------    
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+2*helm->csize+k,
			 2*psize+4*helm->csize+l,bmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+3*helm->csize+k,
			 2*psize+5*helm->csize+l,bmat->c[k][l]);
    }  
  
  //------------------------------------ 
  // fill Vz interior boundary system     
  //------------------------------------  
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			   2*psize+2*helm->csize+l,bmat->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			   2*psize+3*helm->csize+l,bmat->d[k][l]);
    }

  // Vz SYSTEM FILLED----------------------------------------------------------

  // FILL Wz SYSTEM------------------------------------------------------------
  
  //------------------------------------ 
  // fill Wz boundary boundary system   
  //------------------------------------ 
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->asize; ++l){
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+4*helm->asize+k,
			 2+4*helm->asize+l,cmat->a[k][l]);
      SClev->_Am->GenBlk(eid,eid,Asize,Asize,2+5*helm->asize+k,
			 2+5*helm->asize+l,cmat->a[k][l]);
    }
  
  //------------------------------------ 
  // fill Wz boundary interior system   
  //------------------------------------  

  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			 2*psize+4*helm->csize+l,cmat->b[k][l]);
      SClev->_Bm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			 2*psize+5*helm->csize+l,cmat->b[k][l]);
    }

  //------------------------------------ 
  // fill Wz interior interior system     
  //------------------------------------    
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->csize; ++k)  
    for(l = 0; l < helm->csize; ++l){
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+4*helm->csize+k,
			 2*psize+4*helm->csize+l,cmat->c[k][l]);
      SClev->_Cm->GenBlk(eid,eid,Csize,Csize,2*psize+5*helm->csize+k,
			 2*psize+5*helm->csize+l,cmat->c[k][l]);
    }  
  
  //------------------------------------ 
  // fill Wz interior boundary system     
  //------------------------------------  
  
  // REAL/IMAGINARY
  for(k = 0; k < helm->asize; ++k)  
    for(l = 0; l < helm->csize; ++l){
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+4*helm->asize+k,
			   2*psize+4*helm->csize+l,cmat->d[k][l]);
	SClev->_Dm->GenBlk(eid,eid,Asize,Csize,2+5*helm->asize+k,
			   2*psize+5*helm->csize+l,cmat->d[k][l]);
    }

  // Wz SYSTEM FILLED----------------------------------------------------------
}
//=============================================================================
// Solve Au = f using statically condensed system 
void StokesMatrix::SolveSC(double *f){
  int    lda = Ainv->get_lda();
  double *u;
  double *rhs = new double [SClev->_Asize + SClev->_Csize];
  vmath::zero(SClev->_Asize + SClev->_Csize,rhs,1);
  
  // inner solve vector
  u  = new double[lda];
  vmath::zero(lda,u,1);

  // resort input vector into boundary dof followed by int dof
  SClev[0].Pack_BndInt(f,rhs);

  SClev->CondenseV(rhs);
  SClev->PackV(rhs,u);

  // Solve boundary system
  Ainv->Solve(u,1);

  SClev->UnPackV(u,rhs);
  SClev->SolveVd(rhs);
  
  // resort  vector into elemental form
  SClev->UnPack_BndInt(rhs,f);
    
  delete[] rhs;
  delete[] u;
}


void StokesMatrix::Set_SignChange(Element_List *U, Bsystem *Ubsys){
  int     Nbmodes;
  double  *sc,*sc1;
  Element *E;
  
  sc1 = SClev->_signchange = new double [SClev->_Asize];

  if(!Ubsys->signchange)
    setup_signchange(U,Ubsys);
  
  sc  = Ubsys->signchange;
  
  for(E=U->fhead;E;E = E->next){
    sc1[0] = 1; 
    ++sc1;

    sc1[0] = 1;
    ++sc1;    

    Nbmodes = E->Nbmodes;

    vmath::vcopy(Nbmodes,sc,1,sc1,1);
    sc1 += Nbmodes;
    vmath::vcopy(Nbmodes,sc,1,sc1,1);
    sc1 += Nbmodes;
    
    vmath::vcopy(Nbmodes,sc,1,sc1,1);
    sc1 += Nbmodes;
    vmath::vcopy(Nbmodes,sc,1,sc1,1);
    sc1 += Nbmodes;
    vmath::vcopy(Nbmodes,sc,1,sc1,1);
    sc1 += Nbmodes;
    vmath::vcopy(Nbmodes,sc,1,sc1,1);
    sc1 += Nbmodes;

    sc += Nbmodes;
  }
}


void StokesMatrix::Solve(Domain *Omega, int mode){
  int cnt;
  Element *E, *Ef, *EP, *EPf;
  Element *Ei,*Efi,*EPi,*EPfi;    
  Element_List *U = Omega->U;
  Element_List *V = Omega->V;
  Element_List *W = Omega->W;
  Element_List *P = Omega->P;
  double *rhs = new double [SClev->_Asize+SClev->_Csize];
  int nz = 2; // solve works accross two modes per call 
  int u_hjtot = 2*Omega->U->hjtot;
  int p_hjtot = 2*Omega->P->hjtot;

  Bsystem **Pbsys = Omega->Pressure_sys;
  int nsolve = Pbsys[0]->nsolve;
  double *ps  = new double [p_hjtot];
  
  if((mode%2)||(mode >= U->nz)){
    fprintf(stderr,"StokesMatrix::Solve: input parameter mode "
	    "must be a multiple of 2 and less than nz");
    exit(1);
  }

  //  currently not working with a  pressure initial condition prescribed. 
  dzero(p_hjtot,Omega->P->flevels[mode]->base_hj,1);

  //===========================================================================
  // store initial solution 
  //===========================================================================
  dcopy(u_hjtot, U->flevels[mode]->base_hj, 1, Omega->us+mode*u_hjtot, 1);
  dcopy(u_hjtot, V->flevels[mode]->base_hj, 1, Omega->vs+mode*u_hjtot, 1);
  dcopy(u_hjtot, W->flevels[mode]->base_hj, 1, Omega->ws+mode*u_hjtot, 1);  
  dcopy(p_hjtot, P->flevels[mode]->base_hj, 1, ps,  1);

  U->flevels[mode  ]->Trans(U->flevels[mode  ],J_to_Q);
  U->flevels[mode+1]->Trans(U->flevels[mode+1],J_to_Q);
  V->flevels[mode  ]->Trans(V->flevels[mode  ],J_to_Q);
  V->flevels[mode+1]->Trans(V->flevels[mode+1],J_to_Q);
  W->flevels[mode  ]->Trans(W->flevels[mode  ],J_to_Q);
  W->flevels[mode+1]->Trans(W->flevels[mode+1],J_to_Q);

#ifndef PCONTBASE
  Jtransbwd_Orth_Mode(P,P,mode);
#else
  P->flevels[mode  ]->Trans(P->flevels[mode  ],J_to_Q);
  P->flevels[mode+1]->Trans(P->flevels[mode+1],J_to_Q);
#endif

#ifdef PSE_SLV
  //===========================================================================
  // Weak PSE operator on U,V,W and Ui,Vi,Wi
  //===========================================================================
  PSE_Operator(U->flevels[mode  ],V->flevels[mode  ],W->flevels[mode  ],
  	       U->flevels[mode+1],V->flevels[mode+1],W->flevels[mode+1],
  	       P->flevels[mode  ],P->flevels[mode+1],Pbsys[mode]->lambda,
  	       Omega->Usys[mode]->lambda[0].d);

  //===========================================================================
#else
  //===========================================================================
  // Weak Stokes Operator on U,V,W assuming P is zero
  //===========================================================================
  Stokes3D_Operator(U->flevels[mode  ],V->flevels[mode  ],W->flevels[mode  ],
		    U->flevels[mode+1],V->flevels[mode+1],W->flevels[mode+1],
		    P->flevels[mode  ],P->flevels[mode+1],
		    Pbsys[0]->lambda,Omega->Usys[0]->lambda[0].d);

  //===========================================================================
#endif

  //===========================================================================
  // Subtract known solution from rhs
  //===========================================================================
  vmath::vsub(nz*Omega->U->hjtot,Omega->Uf->base_hj,1,Omega->U->base_hj,1,
	      Omega->Uf->base_hj,1);
  vmath::vsub(nz*Omega->V->hjtot,Omega->Vf->base_hj,1,Omega->V->base_hj,1,
	      Omega->Vf->base_hj,1);
  vmath::vsub(nz*Omega->W->hjtot,Omega->Wf->base_hj,1,Omega->W->base_hj,1,
	      Omega->Wf->base_hj,1);
  vmath::vsub(nz*Omega->P->hjtot,Omega->Pf->base_hj,1,Omega->P->base_hj,1,
	      Omega->Pf->base_hj,1);


  //===========================================================================
  // Assemble rhs vector so that solution is in Uf,Vf,Wf,Pf blocks for SolveSC
  //===========================================================================
  cnt = 0;
  for(Ef = Omega->Uf->flevels[mode]->fhead, 
	Efi = Omega->Uf->flevels[mode+1]->fhead,
	EPf = Omega->Pf->flevels[mode]->fhead, 
	EPfi = Omega->Pf->flevels[mode+1]->fhead;
      Ef; Ef = Ef->next, Efi = Efi->next, EPf = EPf->next, EPfi = EPfi->next){

    //=============================
    // Boundary modes
    //-----------------------------    
    // P real and imaginary
    //-----------------------------    
    rhs[cnt++] =  EPf->vert->hj[0];
    rhs[cnt++] = EPfi->vert->hj[0];
    //-----------------------------    
    // U real and imaginary
    //-----------------------------  
    vmath::vcopy(Ef->Nbmodes, Ef->vert->hj,1,rhs+cnt,1);
    cnt += Ef->Nbmodes;   
    vmath::vcopy(Ef->Nbmodes,Efi->vert->hj,1,rhs+cnt,1);
    cnt += Ef->Nbmodes; 
    //-----------------------------    
    // V real and imaginary
    //-----------------------------  
    vmath::vcopy(Ef->Nbmodes,Omega->Vf->flevels[0]->flist[Ef->id]->vert->hj,
		 1,rhs+cnt,1);
    cnt += Ef->Nbmodes;
    vmath::vcopy(Ef->Nbmodes,Omega->Vf->flevels[1]->flist[Ef->id]->vert->hj,
		 1,rhs+cnt,1);
    cnt += Ef->Nbmodes;
    //-----------------------------    
    // W real and imaginary
    //-----------------------------      
    vmath::vcopy(Ef->Nbmodes,Omega->Wf->flevels[0]->flist[Ef->id]->vert->hj,
		 1,rhs+cnt,1);
    cnt += Ef->Nbmodes;
    vmath::vcopy(Ef->Nbmodes,Omega->Wf->flevels[1]->flist[Ef->id]->vert->hj,
		 1,rhs+cnt,1);
    cnt += Ef->Nbmodes;

    //=============================
    // Interior modes
    //-----------------------------    
    // P real and imaginary
    //-----------------------------    
    vmath::vcopy(EPf->Nmodes-1, EPf->vert->hj+1,1,rhs+cnt,1);
    cnt += EPf->Nmodes-1; 
    vmath::vcopy(EPf->Nmodes-1,EPfi->vert->hj+1,1,rhs+cnt,1);
    cnt += EPf->Nmodes-1; 
    //-----------------------------    
    // U real and imaginary
    //-----------------------------  
    vmath::vcopy(Ef->Nmodes-Ef->Nbmodes, Ef->vert->hj+Ef->Nbmodes,1,
		 rhs+cnt,1);
    cnt += Ef->Nmodes-Ef->Nbmodes;
    vmath::vcopy(Ef->Nmodes-Ef->Nbmodes,Efi->vert->hj+Ef->Nbmodes,1,
		 rhs+cnt,1);
    cnt += Ef->Nmodes-Ef->Nbmodes;
    //-----------------------------    
    // V real and imaginary
    //-----------------------------  
    vmath::vcopy(Ef->Nmodes-Ef->Nbmodes,
		 Omega->Vf->flevels[0]->flist[Ef->id]->vert->hj+Ef->Nbmodes,
		 1,rhs+cnt,1);
    cnt += Ef->Nmodes-Ef->Nbmodes;
    vmath::vcopy(Ef->Nmodes-Ef->Nbmodes,
		 Omega->Vf->flevels[1]->flist[Ef->id]->vert->hj+Ef->Nbmodes,
		 1,rhs+cnt,1);
    cnt += Ef->Nmodes-Ef->Nbmodes;
    //-----------------------------    
    // W real and imaginary
    //-----------------------------      
    vmath::vcopy(Ef->Nmodes-Ef->Nbmodes,
		 Omega->Wf->flevels[0]->flist[Ef->id]->vert->hj+Ef->Nbmodes,
		 1,rhs+cnt,1);
    cnt += Ef->Nmodes-Ef->Nbmodes;
    vmath::vcopy(Ef->Nmodes-Ef->Nbmodes,
		 Omega->Wf->flevels[1]->flist[Ef->id]->vert->hj+Ef->Nbmodes,
		 1,rhs+cnt,1);
    cnt += Ef->Nmodes-Ef->Nbmodes;    
  }
  //===========================================================================
  //===========================================================================
  // Static condensation
  //===========================================================================
  SolveSC(rhs);

  //===========================================================================
  // Unpack solution into U,V,W,P and Ui,Vi,Wi,Pi
  //===========================================================================
  cnt = 0;
  for(E = U->flevels[mode]->fhead, Ei = U->flevels[mode+1]->fhead,
	EP = P->flevels[mode]->fhead, EPi = P->flevels[mode+1]->fhead;
      E; E = E->next, Ei = Ei->next, EP = EP->next, EPi = EPi->next){

    //=============================
    // Boundary modes
    //-----------------------------    
    // P real and imaginary
    //-----------------------------   
    EP->vert->hj[0]  = rhs[cnt++]; 
    EPi->vert->hj[0] = rhs[cnt++];     
    //-----------------------------    
    // U real and imaginary
    //-----------------------------     
    vmath::vcopy(E->Nbmodes,rhs+cnt,1, E->vert->hj,1);
    cnt += E->Nbmodes;
    vmath::vcopy(E->Nbmodes,rhs+cnt,1,Ei->vert->hj,1);
    cnt += E->Nbmodes;
    //-----------------------------    
    // V real and imaginary
    //-----------------------------     
    vmath::vcopy(E->Nbmodes,rhs+cnt,1,
		 Omega->V->flevels[0]->flist[E->id]->vert->hj,1);
    cnt += E->Nbmodes;    
    vmath::vcopy(E->Nbmodes,rhs+cnt,1,
		 Omega->V->flevels[1]->flist[E->id]->vert->hj,1);
    cnt += E->Nbmodes;        
    //-----------------------------    
    // W real and imaginary
    //-----------------------------  
    vmath::vcopy(E->Nbmodes,rhs+cnt,1,
		 Omega->W->flevels[0]->flist[E->id]->vert->hj,1);
    cnt += E->Nbmodes;    
    vmath::vcopy(E->Nbmodes,rhs+cnt,1,
		 Omega->W->flevels[1]->flist[E->id]->vert->hj,1);
    cnt += E->Nbmodes;  
    
    //=============================
    // Interior modes
    //-----------------------------    
    // P real and imaginary
    //-----------------------------   
    vmath::vcopy(EP->Nmodes-1,rhs+cnt,1, EP->vert->hj+1,1);
    cnt += EP->Nmodes-1; 
    vmath::vcopy(EP->Nmodes-1,rhs+cnt,1,EPi->vert->hj+1,1);
    cnt += EP->Nmodes-1;     
    //-----------------------------    
    // U real and imaginary
    //-----------------------------     
    vmath::vcopy(E->Nmodes-E->Nbmodes,rhs+cnt,1, E->vert->hj+E->Nbmodes,1);
    cnt += E->Nmodes-E->Nbmodes;
    vmath::vcopy(E->Nmodes-E->Nbmodes,rhs+cnt,1,Ei->vert->hj+E->Nbmodes,1);
    cnt += E->Nmodes-E->Nbmodes;
    //-----------------------------    
    // V real and imaginary
    //-----------------------------     
    vmath::vcopy(E->Nmodes-E->Nbmodes,rhs+cnt,1,
		 Omega->V->flevels[0]->flist[E->id]->vert->hj+E->Nbmodes,1);
    cnt += E->Nmodes-E->Nbmodes;    
    vmath::vcopy(E->Nmodes-E->Nbmodes,rhs+cnt,1,
		 Omega->V->flevels[1]->flist[E->id]->vert->hj+E->Nbmodes,1);
    cnt += E->Nmodes-E->Nbmodes;    
    //-----------------------------    
    // W real and imaginary
    //-----------------------------      
    vmath::vcopy(E->Nmodes-E->Nbmodes,rhs+cnt,1,
		 Omega->W->flevels[0]->flist[E->id]->vert->hj+E->Nbmodes,1);
    cnt += E->Nmodes-E->Nbmodes;    
    vmath::vcopy(E->Nmodes-E->Nbmodes,rhs+cnt,1,
		 Omega->W->flevels[1]->flist[E->id]->vert->hj+E->Nbmodes,1);
    cnt += E->Nmodes-E->Nbmodes;        

  }

  //===========================================================================
  // Add back us,vs,ws,ps to get final solution 
  //===========================================================================
  vmath::vadd(u_hjtot, Omega->us, 1, U->flevels[mode]->base_hj, 1, 
	      U->flevels[mode]->base_hj, 1);
  vmath::vadd(u_hjtot, Omega->vs, 1, V->flevels[mode]->base_hj, 1, 
	      V->flevels[mode]->base_hj, 1);
  vmath::vadd(u_hjtot, Omega->ws, 1, W->flevels[mode]->base_hj, 1, 
	      W->flevels[mode]->base_hj, 1);
  vmath::vadd(p_hjtot, ps, 1, P->flevels[mode]->base_hj, 1,
	      P->flevels[mode]->base_hj, 1);

  delete[] ps;
  delete[] rhs;
}

// ============================================================================
// WEAK PSE OPERATOR (REAL and IMAG cmpnt)
// ============================================================================
//
// Calculate the weak PSE operator on the first U,V,W
// assuming that P is zero and put back into U,V,W
//
// ============================================================================

void StokesMatrix::PSE_Operator(Element_List *U,   Element_List *V,
				Element_List *W,   Element_List *Ui,
				Element_List *Vi,  Element_List *Wi,
				Element_List *Pin, Element_List *Pini,
				Metric *lambda,    double kinvis){
  
  int qt;
  Element *E,*F,*G,*Ei,*Fi,*Gi,*P,*Pi;
  double *dx = new double [QGmax*QGmax];
  double *dy = new double [QGmax*QGmax];

  double beta   = dparam("BETA");             // MSB: BETA => bar(BETA)
  double sigma  = dparam("SIGMA");            // MSB: IMAGINARY cmpnt
  double Omega  = dparam("OMEGA");            // MSB: Omega

  // Declare temporary storage
  double *Ux, *Uy, *Vx, *Vy, *Wx, *Wy;
  double *Eh, *Fh, *Gh, *Ph; 
  double *u, *v, *w, *p, *store, *store1;
  double *A, *B, *Beta, *Sigma;
  
#ifdef PSE_SLV // MSB: For PSE march in z-direction-------
  double idz = 1.0/dparam("DT");
#endif // MSB: -------------------------------------------
  
  for(E=U->fhead,F=V->fhead,G=W->fhead,Ei=Ui->fhead,Fi=Vi->fhead,Gi=Wi->fhead;
      E;E=E->next,F=F->next,G=G->next,Ei=Ei->next,Fi=Fi->next,Gi=Gi->next){
    
    qt = E->qa*E->qb;
    
    P  =  Pin->flist[E->id];
    Pi = Pini->flist[E->id];
      
    //=========================================================================
    // Declare memory and store important initial data
    //=========================================================================

    // Store E->h[0], F->h[0], G->h[0] in temporary storage
    Eh = dvector(0, qt - 1);     
    Fh = dvector(0, qt - 1);     
    Gh = dvector(0, qt - 1);     

    Ph = dvector(0, qt - 1);
    dcopy(qt, P->h[0], 1, Ph, 1);

    dcopy(qt, E->h[0], 1, Eh, 1);
    dcopy(qt, F->h[0], 1, Fh, 1);
    dcopy(qt, G->h[0], 1, Gh, 1);

    A  = dvector(0, qt - 1);          Beta  = dvector(0, qt - 1);
    B  = dvector(0, qt - 1);          Sigma = dvector(0, qt - 1);
    
    store  = dvector(0, qt - 1);   
    store1 = dvector(0, qt - 1);   

    double *Phi = dvector(0, qt - 1); dcopy(qt, Pi->h[0], 1, Phi, 1);
    double *d   = dvector(0, qt - 1);

    // Scaling of Wbar by KINVIS is due to scaling in Set_Oseen routine
    for(int i = 0; i < qt; i++){
      // A[i]      = -kinvis*sigma*sigma;
      // A[i]     -= sigma*kinvis*(lambda+E->id)->wave[2][i]; 
      A[i]      = -sigma*kinvis*(lambda+E->id)->wave[2][i]; 
      B[i]      = -2.0*(kinvis*beta*sigma);           
      B[i]     -= beta*kinvis*(lambda+E->id)->wave[2][i];  
      B[i]     += Omega;
      // Beta[i]   = beta;
      // Sigma[i]  = sigma;

#ifdef PSE_SLV
      A[i]     += idz*(kinvis*(lambda+E->id)->wave[2][i]+2.0*kinvis*sigma);
      B[i]     += idz*(2.0*kinvis*beta);
#endif
    }

#ifdef PSE_SLV // MSB: For PSE march in z-direction-------
    dfill(qt,beta,Beta,1);
    dfill(qt,sigma-idz,Sigma,1);   
#else // MSB: --------------------------------------------
    dfill(qt,beta,Beta,1);
    dfill(qt,sigma,Sigma,1);    
#endif

    //=========================================================================
    // REAL COMPONENT OF PSE OPERATOR
    //=========================================================================
    // Divergence terms stored into P
    //=========================================================================
    E->Grad_d(dx,NULL,NULL,'x');
    F->Grad_d(NULL,dy,NULL,'y');
    vmath::vadd(qt,dx,1,dy,1,P->h[0],1);

#ifndef PCONTBASE
    P->Ofwd(*P->h,P->vert->hj,P->dgL);
#else
    P->Iprod(P);
#endif

    p = dvector(0, P->Nmodes - 1);
    dcopy(P->Nmodes, P->vert->hj, 1, p, 1); 

    //=========================================================================
    // Add on beta and sigma terms and store in P (nb: negate Beta, Sigma)
    //=========================================================================
    dcopy (qt, Ph, 1, P->h[0], 1);
    dneg  (qt, Beta,  1);      dneg(qt, Sigma, 1);
    dvmul (qt, Sigma, 1, G->h[0],  1, P->h[0], 1);  
    dvvtvp(qt, Beta,  1, Gi->h[0], 1, P->h[0], 1, P->h[0], 1);
    
#ifndef PCONTBASE
    P->Ofwd(*P->h,P->vert->hj,P->dgL);
#else
    P->Iprod(P);
#endif  

    dvadd (P->Nmodes, p, 1, P->vert->hj, 1, P->vert->hj, 1);  
    dneg(P->Nmodes, P->vert->hj, 1);       // MSB: Negate result as required!!!

    //=========================================================================
    // Helmholtz operator L acting on u,v,w (includes oseen terms)
    //=========================================================================
    E->HelmHoltz(lambda+E->id);
    blas::dscal(E->Nmodes,kinvis,E->vert->hj,1);
    F->HelmHoltz(lambda+E->id);
    blas::dscal(F->Nmodes,kinvis,F->vert->hj,1);
    G->HelmHoltz(lambda+E->id);
    blas::dscal(G->Nmodes,kinvis,G->vert->hj,1);      
      
    //=========================================================================
    // Helmholtz operator L stored in u,v,w respectively (includes oseen terms)
    //=========================================================================
    u = dvector(0, E->Nmodes - 1);
    v = dvector(0, F->Nmodes - 1);
    w = dvector(0, G->Nmodes - 1);

    dcopy(E->Nmodes, E->vert->hj, 1, u, 1);  
    dcopy(F->Nmodes, F->vert->hj, 1, v, 1); 
    dcopy(G->Nmodes, G->vert->hj, 1, w, 1);  

    // Recover E->h[0], F->h[0], G->h[0]
    dcopy(qt, Eh, 1, E->h[0], 1); 
    dcopy(qt, Fh, 1, F->h[0], 1); 
    dcopy(qt, Gh, 1, G->h[0], 1); 

    //=========================================================================
    // The terms associated with the Pressure field
    //=========================================================================
    E->GradT_h(Ph, d, NULL, NULL, 'x');
    dcopy(qt, d, 1, E->h[0], 1);
    E->Iprod(E);
    dcopy(qt, Eh, 1, E->h[0], 1);
    dvsub(E->Nmodes, u, 1, E->vert->hj, 1, u, 1);         // MSB: Check sign!!!

    F->GradT_h(Ph, NULL, d, NULL, 'y');
    dcopy(qt, d, 1, F->h[0], 1);
    F->Iprod(F);
    dcopy(qt, Fh, 1, F->h[0], 1);
    dvsub(F->Nmodes, v, 1, F->vert->hj, 1, v, 1);         // MSB: Check sign!!!

    dsmul(qt, -beta, Phi, 1, G->h[0], 1);                 // MSB: Check sign!!!
#ifdef PSE_SLV
#ifdef NODPDZ
    dsvtvp(qt, -sigma, Ph, 1, G->h[0], 1, G->h[0], 1);
#else
    dsvtvp(qt, idz-sigma, Ph, 1, G->h[0], 1, G->h[0], 1);
#endif
#else
    dsvtvp(qt, -sigma, Ph, 1, G->h[0], 1, G->h[0], 1);    // MSB: Check sign!!!
#endif
    G->Iprod(G);
    dcopy(qt, Gh, 1, G->h[0], 1);
    dvadd(G->Nmodes, w, 1, G->vert->hj, 1, w, 1);

    //=========================================================================
    // The u.dU/dx, v.dU/dy, u.dV/dx, v.dV/dy, u.dW/dx, v.dW/dy terms
    //=========================================================================
    Ux = dvector(0, qt - 1); 
    Uy = dvector(0, qt - 1); 
    Vx = dvector(0, qt - 1); 
    Vy = dvector(0, qt - 1); 
    Wx = dvector(0, qt - 1);
    Wy = dvector(0, qt - 1);

    E->Grad_h((lambda+E->id)->wave[0], Ux, NULL, NULL, 'x');
    E->Grad_h((lambda+E->id)->wave[0], NULL, Uy, NULL, 'y');
    E->Grad_h((lambda+E->id)->wave[1], Vx, NULL, NULL, 'x');
    E->Grad_h((lambda+E->id)->wave[1], NULL, Vy, NULL, 'y');
    E->Grad_h((lambda+E->id)->wave[2], Wx, NULL, NULL, 'x');
    E->Grad_h((lambda+E->id)->wave[2], NULL, Wy, NULL, 'y');

    // Scale by kinematic viscosity
    // Note that scaling by KINVIS is required because of the
    // way in which lambda->wave was scaled in Set_Oseen
    dscal(qt, kinvis, Ux, 1);
    dscal(qt, kinvis, Uy, 1);
    dscal(qt, kinvis, Vx, 1);
    dscal(qt, kinvis, Vy, 1);
    dscal(qt, kinvis, Wx, 1);
    dscal(qt, kinvis, Wy, 1);

    //------------------------------------------------------
    // u cmpnt
    //------------------------------------------------------
    // u.dU/dx + v.dU/dy   
    dvmul (qt, Ux, 1, E->h[0], 1, E->h[0], 1);
    dvvtvp(qt, Uy, 1, F->h[0], 1, E->h[0], 1, store, 1);

    // Recover E->h[0] from start
    dcopy(qt, Eh, 1, E->h[0], 1); 

    //------------------------------------------------------
    // v cmpnt
    //------------------------------------------------------
    // u.dV/dx + v.dV/dy
    dvmul (qt, Vx, 1, E->h[0], 1, E->h[0], 1);
    dvvtvp(qt, Vy, 1, F->h[0], 1, E->h[0], 1, store1, 1);

    // Recover E->h[0], F->h[0] from start
    dcopy(qt, Eh, 1, E->h[0], 1);     
    dcopy(qt, Fh, 1, F->h[0], 1);  

    //------------------------------------------------------
    // w cmpnt
    //------------------------------------------------------
    // u.dW/dx + v.dW/dy
    dvmul (qt, Wx, 1, E->h[0], 1, E->h[0], 1);
    dvvtvp(qt, Wy, 1, F->h[0], 1, E->h[0], 1, G->h[0], 1);  

    // Set E->h[0] to store and F->h[0] to store1
    dcopy(qt, store,  1, E->h[0], 1); 
    dcopy(qt, store1, 1, F->h[0], 1);    
    
    // Inner products
    E->Iprod(E);   
    F->Iprod(F);  
    G->Iprod(G);  

    //=========================================================================
    // The u.dU/dx, v.dU/dy, u.dV/dx, v.dV/dy, u.dW/dx, v.dW/dy terms stored
    //=========================================================================
    dvadd (E->Nmodes, u, 1, E->vert->hj, 1, u, 1);
    dvadd (F->Nmodes, v, 1, F->vert->hj, 1, v, 1);  
    dvadd (G->Nmodes, w, 1, G->vert->hj, 1, w, 1); 

    //=========================================================================
    // The a, b, c, d terms
    //=========================================================================
    
    // Recover E->h[0], F->h[0], G->h[0] from start    
    dcopy(qt, Eh, 1, E->h[0], 1); 
    dcopy(qt, Fh, 1, F->h[0], 1);    
    dcopy(qt, Gh, 1, G->h[0], 1);       

    //------------------------------------------------------
    // u cmpnt
    //------------------------------------------------------    
    dvmul (qt, A, 1, E->h[0],  1, E->h[0], 1);
    dvvtvp(qt, B, 1, Ei->h[0], 1, E->h[0], 1, E->h[0], 1);  

    //------------------------------------------------------
    // v cmpnt
    //------------------------------------------------------    
    dvmul (qt, A, 1, F->h[0],  1, F->h[0], 1);
    dvvtvp(qt, B, 1, Fi->h[0], 1, F->h[0], 1, F->h[0], 1);  

    //------------------------------------------------------
    // w cmpnt
    //------------------------------------------------------    
    dvmul (qt, A, 1, G->h[0],  1, G->h[0], 1);
    dvvtvp(qt, B, 1, Gi->h[0], 1, G->h[0], 1, G->h[0], 1);  

    // Inner products
    E->Iprod(E);   
    F->Iprod(F);  
    G->Iprod(G);      
    
    //=========================================================================
    // The a, b, c, d terms stored
    //=========================================================================
    // REAL COMPONENT OF OPERATOR STORED IN E, F, G
    //=========================================================================
    dvadd (E->Nmodes, u, 1, E->vert->hj, 1, E->vert->hj, 1);
    dvadd (F->Nmodes, v, 1, F->vert->hj, 1, F->vert->hj, 1);  
    dvadd (G->Nmodes, w, 1, G->vert->hj, 1, G->vert->hj, 1); 

    //=========================================================================
    // The L1 matrix terms
    //=========================================================================
    // Need to add Uz.u term to E->vert->hj (x-mom'm)
    // Need to add Vz.v term to F->vert->hj (y-mom'm)
    // Need to add Wz.w term to G->vert->hj (z-mom'm)

    //=========================================================================
    // IMAGINARY COMPONENT OF PSE OPERATOR
    //=========================================================================
    // Recover initial data
    //=========================================================================
    dcopy(qt, Eh, 1, E->h[0], 1);     
    dcopy(qt, Fh, 1, F->h[0], 1);      
    dcopy(qt, Gh, 1, G->h[0], 1);     
    dcopy(qt, Ph, 1, P->h[0], 1);  

    dcopy (qt, Ei->h[0], 1, Eh, 1);  
    dcopy (qt, Fi->h[0], 1, Fh, 1);  
    dcopy (qt, Gi->h[0], 1, Gh, 1);  
    dcopy (qt, Pi->h[0], 1, Ph, 1);    

    //=========================================================================
    // Divergence terms stored into Pi
    //=========================================================================
    Ei->Grad_d(dx,NULL,NULL,'x');
    Fi->Grad_d(NULL,dy,NULL,'y');
    vmath::vadd(qt,dx,1,dy,1,Pi->h[0],1);

#ifndef PCONTBASE
    Pi->Ofwd(*Pi->h,Pi->vert->hj,Pi->dgL);
#else
    Pi->Iprod(Pi);
#endif
    
    dcopy(Pi->Nmodes, Pi->vert->hj, 1, p, 1); 

    //=========================================================================
    // Add on beta and sigma terms and store in Pi (nb: negate Beta)
    //=========================================================================
    dcopy (qt, Ph, 1, Pi->h[0], 1);
    dneg  (qt, Beta,  1);
    dvmul (qt, Beta,  1, G->h[0],  1, Pi->h[0], 1);  
    dvvtvp(qt, Sigma, 1, Gi->h[0], 1, Pi->h[0], 1, Pi->h[0], 1);

#ifndef PCONTBASE
    Pi->Ofwd(*Pi->h,Pi->vert->hj,Pi->dgL);
#else
    Pi->Iprod(Pi);
#endif  

    dvadd (Pi->Nmodes, p, 1, Pi->vert->hj, 1, Pi->vert->hj, 1);  
    dneg(Pi->Nmodes, Pi->vert->hj, 1);     // MSB: Negate result as required!!!
    //=========================================================================
    // Helmholtz operator L acting on u,v,w (includes oseen terms)
    //=========================================================================
    Ei->HelmHoltz(lambda+E->id);
    blas::dscal(Ei->Nmodes,kinvis,Ei->vert->hj,1);
    Fi->HelmHoltz(lambda+E->id);
    blas::dscal(Fi->Nmodes,kinvis,Fi->vert->hj,1);
    Gi->HelmHoltz(lambda+E->id);
    blas::dscal(Gi->Nmodes,kinvis,Gi->vert->hj,1);      
      
    //=========================================================================
    // Helmholtz operator L stored in u,v,w respectively (includes oseen terms)
    //=========================================================================
    dcopy(Ei->Nmodes, Ei->vert->hj, 1, u, 1);  
    dcopy(Fi->Nmodes, Fi->vert->hj, 1, v, 1); 
    dcopy(Gi->Nmodes, Gi->vert->hj, 1, w, 1);  

    // Recover Ei->h[0], Fi->h[0], Gi->h[0]
    dcopy(qt, Eh, 1, Ei->h[0], 1); 
    dcopy(qt, Fh, 1, Fi->h[0], 1); 
    dcopy(qt, Gh, 1, Gi->h[0], 1); 

    //=========================================================================
    // The terms associated with the Pressure field
    //=========================================================================
    dcopy (qt, P->h[0], 1, Phi, 1);   // NB: Ph  stores imaginary (see earlier)
                                      // NB: Phi stores real
    Ei->GradT_h(Ph, d, NULL, NULL, 'x');
    dcopy(qt, d, 1, Ei->h[0], 1);
    Ei->Iprod(Ei);
    dcopy(qt, Eh, 1, Ei->h[0], 1);
    dvsub(Ei->Nmodes, u, 1, Ei->vert->hj, 1, u, 1);       // MSB: Check sign!!!

    Fi->GradT_h(Ph, NULL, d, NULL, 'y');
    dcopy(qt, d, 1, Fi->h[0], 1);
    Fi->Iprod(Fi);
    dcopy(qt, Fh, 1, Fi->h[0], 1);
    dvsub(Fi->Nmodes, v, 1, Fi->vert->hj, 1, v, 1);       // MSB: Check sign!!!

    dsmul(qt, beta, Phi, 1, Gi->h[0], 1);                 // MSB: Check sign!!!
#ifdef PSE_SLV
#ifdef NODPDZ
    dsvtvp(qt, -sigma, Ph, 1, Gi->h[0], 1, Gi->h[0], 1);
#else
    dsvtvp(qt, idz-sigma, Ph, 1, Gi->h[0], 1, Gi->h[0], 1);
#endif
#else
    dsvtvp(qt, -sigma, Ph, 1, Gi->h[0], 1, Gi->h[0], 1);  // MSB: Check sign!!!
#endif
    Gi->Iprod(Gi);
    dcopy(qt, Gh, 1, Gi->h[0], 1);
    dvadd(Gi->Nmodes, w, 1, Gi->vert->hj, 1, w, 1);

    //=========================================================================
    // The u.dU/dx, v.dU/dy, u.dV/dx, v.dV/dy, u.dW/dx, v.dW/dy terms
    //=========================================================================
    
    //------------------------------------------------------
    // u cmpnt
    //------------------------------------------------------
    // u.dU/dx + v.dU/dy   
    dvmul (qt, Ux, 1, Ei->h[0], 1, Ei->h[0], 1);
    dvvtvp(qt, Uy, 1, Fi->h[0], 1, Ei->h[0], 1, store, 1);

    // Recover Ei->h[0] from start
    dcopy(qt, Eh, 1, Ei->h[0], 1); 

    //------------------------------------------------------
    // v cmpnt
    //------------------------------------------------------
    // u.dV/dx + v.dV/dy
    dvmul (qt, Vx, 1, Ei->h[0], 1, Ei->h[0], 1);
    dvvtvp(qt, Vy, 1, Fi->h[0], 1, Ei->h[0], 1, store1, 1);

    // Recover Ei->h[0], Fi->h[0] from start
    dcopy(qt, Eh, 1, Ei->h[0], 1);     
    dcopy(qt, Fh, 1, Fi->h[0], 1);  

    //------------------------------------------------------
    // w cmpnt
    //------------------------------------------------------
    // u.dW/dx + v.dW/dy
    dvmul (qt, Wx, 1, Ei->h[0], 1, Ei->h[0], 1);
    dvvtvp(qt, Wy, 1, Fi->h[0], 1, Ei->h[0], 1, Gi->h[0], 1);  

    // Set Ei->h[0] to store and Fi->h[0] to store1
    dcopy(qt, store,  1, Ei->h[0], 1); 
    dcopy(qt, store1, 1, Fi->h[0], 1);    
    
    // Inner products
    Ei->Iprod(Ei);   
    Fi->Iprod(Fi);  
    Gi->Iprod(Gi);  

    //=========================================================================
    // The u.dU/dx, v.dU/dy, u.dV/dx, v.dV/dy, u.dW/dx, v.dW/dy terms stored
    //=========================================================================
    dvadd (Ei->Nmodes, u, 1, Ei->vert->hj, 1, u, 1);
    dvadd (Fi->Nmodes, v, 1, Fi->vert->hj, 1, v, 1);  
    dvadd (Gi->Nmodes, w, 1, Gi->vert->hj, 1, w, 1); 

    //=========================================================================
    // The a, b, c, d terms (nb: negate B)
    //=========================================================================
    
    // Recover Ei->h[0], Fi->h[0], Gi->h[0] from start    
    dcopy(qt, Eh, 1, Ei->h[0], 1); 
    dcopy(qt, Fh, 1, Fi->h[0], 1);    
    dcopy(qt, Gh, 1, Gi->h[0], 1);       

    //------------------------------------------------------
    // u cmpnt
    //------------------------------------------------------ 
    dneg  (qt, B, 1);
    dvmul (qt, A, 1, Ei->h[0], 1, Ei->h[0], 1);
    dvvtvp(qt, B, 1, E->h[0],  1, Ei->h[0], 1, Ei->h[0], 1);  

    //------------------------------------------------------
    // v cmpnt
    //------------------------------------------------------    
    dvmul (qt, A, 1, Fi->h[0], 1, Fi->h[0], 1);
    dvvtvp(qt, B, 1, F->h[0],  1, Fi->h[0], 1, Fi->h[0], 1);  

    //------------------------------------------------------
    // w cmpnt
    //------------------------------------------------------    
    dvmul (qt, A, 1, Gi->h[0], 1, Gi->h[0], 1);
    dvvtvp(qt, B, 1, G->h[0],  1, Gi->h[0], 1, Gi->h[0], 1);  

    // Inner products
    Ei->Iprod(Ei);   
    Fi->Iprod(Fi);  
    Gi->Iprod(Gi);      
    
    //=========================================================================
    // The a, b, c, d terms stored
    //=========================================================================
    // IMAGINARY COMPONENT OF OPERATOR STORED IN Ei, Fi, Gi
    //=========================================================================
    dvadd (Ei->Nmodes, u, 1, Ei->vert->hj, 1, Ei->vert->hj, 1);
    dvadd (Fi->Nmodes, v, 1, Fi->vert->hj, 1, Fi->vert->hj, 1);  
    dvadd (Gi->Nmodes, w, 1, Gi->vert->hj, 1, Gi->vert->hj, 1); 

    //=========================================================================
    // The L1 matrix terms
    //=========================================================================
    // Need to add Uz.u term to E->vert->hj (x-mom'm)
    // Need to add Vz.v term to F->vert->hj (y-mom'm)
    // Need to add Wz.w term to G->vert->hj (z-mom'm)

    //=========================================================================
    // Free all memory
    //=========================================================================
    free(Ux);    free(Uy);     free(Vx);   free(Vy);   free(Wx);   free(Wy);
    free(Eh);    free(Fh);     free(Gh);   free(Ph); 
    free(Beta);  free(Sigma);  free(A);    free(B);
    free(u);     free(v);      free(w);    free(p);
    free(store); free(store1);

    free(Phi);   free(d);

  }

  delete[] dx;
  delete[] dy;
}

void StokesMatrix::Stokes3D_Operator(Element_List *U,   Element_List *V,
				     Element_List *W,   Element_List *Ui,
				     Element_List *Vi,  Element_List *Wi,
				     Element_List *Pin, Element_List *Pini,
				     Metric *lambda,    double kinvis){  
  int qt;
  Element *E,*F,*G,*Ei,*Fi,*Gi,*P,*Pi;
  double *dx      = new double [QGmax*QGmax];
  double *dy      = new double [QGmax*QGmax];
  double *pstore  = new double [QGmax*QGmax];
  double *pistore = new double [QGmax*QGmax];
  double  beta    = dparam("BETA");  
  
  for(E=U->fhead,F=V->fhead,G=W->fhead,Ei=Ui->fhead,Fi=Vi->fhead,Gi=Wi->fhead;
      E;E=E->next,F=F->next,G=G->next,Ei=Ei->next,Fi=Fi->next,Gi=Gi->next){
    
    qt = E->qa*E->qb;
    
    P  =  Pin->flist[E->id];
    Pi = Pini->flist[E->id];
      
    dcopy     (qt,P->h[0],1,pstore,1);
    dcopy     (qt,Pi->h[0],1,pistore,1);
    //=========================================================================
    // Pressure evaluation
    //=========================================================================
    //=========================================================================
    // Divergence terms stored into P
    //=========================================================================

    E->Grad_d(dx,NULL,NULL,'x');
    F->Grad_d(NULL,dy,NULL,'y');
    vmath::vadd(qt,dx,1,dy,1,P->h[0],1);
    dsvtvp(qt, -beta, Gi->h[0], 1, P->h[0], 1, P->h[0], 1);

#ifndef PCONTBASE
    P->Ofwd(*P->h,P->vert->hj,P->lmax);
#else
    P->Iprod(P);
#endif
    dneg(P->Nmodes, P->vert->hj, 1);       // MSB: Negate result as required!!!

    //=========================================================================
    // Divergence terms stored into Pi
    //=========================================================================
    Ei->Grad_d(dx,NULL,NULL,'x');
    Fi->Grad_d(NULL,dy,NULL,'y');
    vmath::vadd(qt,dx,1,dy,1,Pi->h[0],1);
    dsvtvp (qt, beta, G->h[0], 1, Pi->h[0], 1, Pi->h[0], 1);

#ifndef PCONTBASE
    Pi->Ofwd(*Pi->h,Pi->vert->hj,Pi->lmax);
#else
    Pi->Iprod(Pi);
#endif
    dneg(Pi->Nmodes, Pi->vert->hj, 1);     // MSB: Negate result as required!!!


    //=========================================================================
    // REAL COMPONENT 
    //=========================================================================
    //=========================================================================
    // Helmholtz operator L acting on u,v,w (includes oseen terms)
    //=========================================================================
    double *store = dx;
    double *d     = dy;

    dcopy(qt,E->h[0],1,store,1);             // store solution 
    E->GradT_h(pstore,d,NULL,NULL,'x');      // calculate (phi,dp/dx)
    dcopy(qt,d,1,E->h[0],1);
    E->Iprod(E);
    dcopy(qt,store,1,E->h[0],1);             // copy back E->hj 
    dcopy(E->Nmodes,E->vert->hj,1,store,1);  // store E->vert->hj
    E->HelmHoltz(lambda+E->id);              // calculate Laplacian 
    blas::dscal(E->Nmodes,kinvis,E->vert->hj,1);

    dvsub(E->Nmodes,E->vert->hj,1,store,1,E->vert->hj,1);

    dcopy(qt,F->h[0],1,store,1);             // store solution 
    F->GradT_h(pstore,NULL,d,NULL,'y');      // calculate (phi,d/dy)
    dcopy(qt,d,1,F->h[0],1);
    F->Iprod(F);
    dcopy(qt,store,1,F->h[0],1);             // copy back E->hj 
    dcopy(F->Nmodes,F->vert->hj,1,store,1);  // store E->vert->hj
    F->HelmHoltz(lambda+E->id);
    blas::dscal(F->Nmodes,kinvis,F->vert->hj,1);

    dvsub(F->Nmodes,F->vert->hj,1,store,1,F->vert->hj,1);

    dcopy(qt,G->h[0],1,d,1);                 // store solution 
    dsmul(qt,-beta,pistore,1,G->h[0],1);     // calculate -beta (phi,pi)
    G->Iprod(G);
    dcopy(qt,d,1,G->h[0],1);                 // copy back E->hj 
    dcopy(G->Nmodes,G->vert->hj,1,store,1);  // store E->vert->hj

    G->HelmHoltz(lambda+E->id);
    blas::dscal(G->Nmodes,kinvis,G->vert->hj,1);      
    dvsub(G->Nmodes,G->vert->hj,1,store,1,G->vert->hj,1);
    //=========================================================================
    // IMAGINARY COMPONENT 
    //=========================================================================

    //=========================================================================
    // Helmholtz operator L acting on u,v,w (includes oseen terms)
    //=========================================================================
    dcopy(qt,Ei->h[0],1,store,1);              // store solution 
    Ei->GradT_h(pistore,d,NULL,NULL,'x');      // calculate (phi,dp/dx)
    dcopy(qt,d,1,Ei->h[0],1);
    Ei->Iprod(Ei);
    dcopy(qt,store,1,Ei->h[0],1);              // copy back E->hj 
    dcopy(Ei->Nmodes,Ei->vert->hj,1,store,1);  // store E->vert->hj
    Ei->HelmHoltz(lambda+E->id);
    blas::dscal(Ei->Nmodes,kinvis,Ei->vert->hj,1);
    dvsub(Ei->Nmodes,Ei->vert->hj,1,store,1,Ei->vert->hj,1);

    dcopy(qt,Fi->h[0],1,store,1);              // store solution 
    Fi->GradT_h(pistore,NULL,d,NULL,'y');      // calculate (phi,d/dy)
    dcopy(qt,d,1,Fi->h[0],1);
    Fi->Iprod(Fi);
    dcopy(qt,store,1,Fi->h[0],1);              // copy back E->hj 
    dcopy(Fi->Nmodes,Fi->vert->hj,1,store,1);  // store E->vert->hj
    Fi->HelmHoltz(lambda+E->id);
    blas::dscal(Fi->Nmodes,kinvis,Fi->vert->hj,1);
    dvsub(Fi->Nmodes,Fi->vert->hj,1,store,1,Fi->vert->hj,1);

    dcopy(qt,Gi->h[0],1,store,1);             // store solution 
    dsmul(qt,beta,pstore,1,Gi->h[0],1);       // calculate beta (phi,p)
    Gi->Iprod(Gi);
    dcopy(qt,store,1,Gi->h[0],1);              // copy back E->hj 
    dcopy(Gi->Nmodes,Gi->vert->hj,1,store,1); // store E->vert->hj
    Gi->HelmHoltz(lambda+E->id);
    blas::dscal(Gi->Nmodes,kinvis,Gi->vert->hj,1);      
    dvsub(Gi->Nmodes,Gi->vert->hj,1,store,1,Gi->vert->hj,1);
  }

  delete[] dx;
  delete[] dy;
  delete[] pstore;
  delete[] pistore;
}

// ============================================================================
// CALCULATE L2.q
// ============================================================================
void StokesMatrix::PSEForce(Domain *Omega){
  int k;

  Bsystem **Pbsys = Omega->Pressure_sys;

  dcopy(Omega->P->htot*Omega->P->nz,Omega->P->base_h, 1, Omega->Pf->base_h, 1);
  dcopy(Omega->U->htot*Omega->U->nz,Omega->U->base_h, 1, Omega->Uf->base_h, 1);
  dcopy(Omega->V->htot*Omega->V->nz,Omega->V->base_h, 1, Omega->Vf->base_h, 1);
  dcopy(Omega->W->htot*Omega->W->nz,Omega->W->base_h, 1, Omega->Wf->base_h, 1);

  //===========================================================================
  // CALL L2_Operator TO FORM L2.q OR PSE_Operator TO FORM IDENTITY MATRIX
  //===========================================================================

  for(k = 0; k < Omega->U->nz; k+=2){
    L2_Operator(Omega->Uf->flevels[k  ], Omega->Vf->flevels[k  ],
		Omega->Wf->flevels[k  ], Omega->Uf->flevels[k+1],
		Omega->Vf->flevels[k+1], Omega->Wf->flevels[k+1],
		Omega->Pf->flevels[k  ], Omega->Pf->flevels[k+1],
		Pbsys[k]->lambda, Omega->Usys[k]->lambda[0].d);
  }

} 

// ============================================================================
// WEAK L2 OPERATOR (REAL and IMAG cmpnt)
// ============================================================================
void StokesMatrix::L2_Operator(Domain *Omega){
  int k;

  for(k = 0; k < Omega->U->nz; k+=2){
    L2_Operator(Omega->Uf->flevels[k  ], Omega->Vf->flevels[k  ],
		Omega->Wf->flevels[k  ], Omega->Uf->flevels[k+1],
		Omega->Vf->flevels[k+1], Omega->Wf->flevels[k+1],
		Omega->Pf->flevels[k  ], Omega->Pf->flevels[k+1],
		Omega->Pressure_sys[k]->lambda,Omega->Usys[k]->lambda[0].d);
  }
}

void StokesMatrix::L2_Operator(Element_List *U,   Element_List *V,
			       Element_List *W,   Element_List *Ui,
			       Element_List *Vi,  Element_List *Wi,
			       Element_List *Pin, Element_List *Pini,
			       Metric *lambda,    double kinvis){

  int qt;
  Element *E,*F,*G,*Ei,*Fi,*Gi,*P,*Pi;
  
  double beta   = dparam("BETA");             // MSB: BETA => bar(BETA)
  double sigma  = dparam("SIGMA");            // MSB: IMAGINARY cmpnt
  double Omega  = dparam("OMEGA");            // MSB: Omega
  
  // Declare temporary storage
  double *Eh, *Fh, *Gh, *Ph, *B, *A; 
  double *IDZ;
  
  // MSB: For PSE march in z-direction--------------------
  double idz = 1.0/dparam("DT");  
  // MSB: ------------------------------------------------
   for(E=U->fhead,F=V->fhead,G=W->fhead,Ei=Ui->fhead,Fi=Vi->fhead,Gi=Wi->fhead;
      E;E=E->next,F=F->next,G=G->next,Ei=Ei->next,Fi=Fi->next,Gi=Gi->next){
    
    qt = E->qa*E->qb;

    P  =  Pin->flist[E->id];
    Pi = Pini->flist[E->id];

    //=========================================================================
    // Declare memory and store important initial data
    //=========================================================================

    // Store E->h[0], F->h[0], G->h[0] in temporary storage
    Eh = dvector(0, qt - 1);     
    Fh = dvector(0, qt - 1);     
    Gh = dvector(0, qt - 1);     

    Ph = dvector(0, qt - 1);
    dcopy(qt, P->h[0], 1, Ph, 1);

    dcopy(qt, E->h[0], 1, Eh, 1);
    dcopy(qt, F->h[0], 1, Fh, 1);
    dcopy(qt, G->h[0], 1, Gh, 1);

    A   = dvector(0, qt - 1);        
    B   = dvector(0, qt - 1);        
    IDZ = dvector(0, qt - 1);        
  
    // Scaling of Wbar by KINVIS is due to scaling in Set_Oseen routine
      
    for(int i = 0; i < qt; i++){
      A[i]     = idz*(2.0*kinvis*sigma);
      B[i]     = idz*(2.0*kinvis*beta);
      IDZ[i]   = idz;                    // MSB: Check the sign of this!!!
    }

    if(lambda[E->id].wave != NULL)
      for(int i = 0; i < qt; i++)
	A[i]     += idz*(kinvis*(lambda+E->id)->wave[2][i]);

    //=========================================================================
    // REAL COMPONENT OF L2 OPERATOR
    //=========================================================================
    // P-cmpnt
    //=========================================================================
    dcopy (qt, G->h[0], 1, P->h[0], 1);    
    blas::dscal(qt, -idz, P->h[0], 1);

#ifndef PCONTBASE
    P->Ofwd(*P->h,P->vert->hj,P->dgL);
#else
    P->Iprod(P);
#endif  

    dcopy(qt, Ph, 1, P->h[0], 1);  

    //=========================================================================
    // U-cmpnt
    //=========================================================================
    dvmul (qt, A, 1, E->h[0],  1, E->h[0], 1);    
    dvvtvp(qt, B, 1, Ei->h[0], 1, E->h[0], 1, E->h[0], 1);  
    
    E->Iprod(E);

    //=========================================================================
    // V-cmpnt
    //=========================================================================
    dvmul (qt, A, 1, F->h[0],  1, F->h[0], 1);    
    dvvtvp(qt, B, 1, Fi->h[0], 1, F->h[0], 1, F->h[0], 1);  
    
    F->Iprod(F);

    //=========================================================================
    // W-cmpnt
    //=========================================================================
    dvmul (qt, A,   1, G->h[0],  1, G->h[0], 1);    
    dvvtvp(qt, B,   1, Gi->h[0], 1, G->h[0], 1, G->h[0], 1);
#ifndef NODPDZ 
    dvvtvp(qt, IDZ, 1, P->h[0],  1, G->h[0], 1, G->h[0], 1); 
#endif
    G->Iprod(G);

    //=========================================================================
    // IMAGINARY COMPONENT OF PSE OPERATOR
    //=========================================================================
    // Recover initial data
    //=========================================================================
    dcopy(qt, Eh, 1, E->h[0], 1);     
    dcopy(qt, Fh, 1, F->h[0], 1);      
    dcopy(qt, Gh, 1, G->h[0], 1);     
    dcopy(qt, Ph, 1, P->h[0], 1);  

    dcopy (qt, Ei->h[0], 1, Eh, 1);  
    dcopy (qt, Fi->h[0], 1, Fh, 1);  
    dcopy (qt, Gi->h[0], 1, Gh, 1);  
    dcopy (qt, Pi->h[0], 1, Ph, 1); 
    //=========================================================================
    // P-cmpnt
    //=========================================================================
    dcopy (qt, Gi->h[0], 1, Pi->h[0], 1);    
    blas::dscal(qt, -idz, Pi->h[0], 1);

#ifndef PCONTBASE2
    Pi->Ofwd(*Pi->h,Pi->vert->hj,Pi->dgL);
#else
    Pi->Iprod(Pi);
#endif
  
    dcopy(qt, Ph, 1, Pi->h[0], 1);  

    //=========================================================================
    // U-cmpnt
    //=========================================================================
    dneg  (qt, B, 1);
    dvmul (qt, A, 1, Ei->h[0],  1, Ei->h[0], 1);    
    dvvtvp(qt, B, 1, E->h[0],   1, Ei->h[0], 1, Ei->h[0], 1);  
    
    Ei->Iprod(Ei);

    //=========================================================================
    // V-cmpnt
    //=========================================================================
    dvmul (qt, A, 1, Fi->h[0],  1, Fi->h[0], 1);    
    dvvtvp(qt, B, 1, F->h[0],   1, Fi->h[0], 1, Fi->h[0], 1);  
    
    Fi->Iprod(Fi);

    //=========================================================================
    // W-cmpnt
    //=========================================================================
    dvmul (qt, A,   1, Gi->h[0],  1, Gi->h[0], 1);    
    dvvtvp(qt, B,   1, G->h[0],   1, Gi->h[0], 1, Gi->h[0], 1); 
#ifndef NODPDZ  
    dvvtvp(qt, IDZ, 1, Pi->h[0],  1, Gi->h[0], 1, Gi->h[0], 1);
#endif    
    Gi->Iprod(Gi);

    //=========================================================================
    // Free all memory
    //=========================================================================
    free(Eh); free(Fh); free(Gh);  free(Ph); 
    free(B);  free(A);  free(IDZ);

  }
}



