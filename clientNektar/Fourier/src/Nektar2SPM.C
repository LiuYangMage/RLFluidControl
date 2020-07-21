
#include <mpi.h>
#include "nektarF.h"
#include "C_SPM.h"
#include <rfftw.h>

extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;

void compute_indicator_function(Domain *omega, C_SPM *SPMv);
void compute_hydrodynamic_force(Domain *Omega);
void compute_upf(Domain *omega);


void compute_indicator_function(Domain *omega, C_SPM *SPMv){

  Element_List  *EL = omega->U;
  Element *E;

  int nz = EL->nz;
  int htot = EL->htot*EL->nz;
   
  Coord X;
  X.x = dvector(0,QGmax*QGmax-1);
  X.y = dvector(0,QGmax*QGmax-1);

  int i,j,k,plane, nel, qtot,offset = 0;
  double lx, ly, lz, dt, r;
  double *ind_func;
  int part_ID;   

  nel = EL->nel;

  for (part_ID = 0; part_ID < SPMv[0].num_partls; ++part_ID){
     ind_func = SPMv[0].gaussian[part_ID];
 
     offset = 0; 
    for(k=0,plane=nz*option("PROCID"); k<nz; ++k,++plane)
      for (i = 0; i < nel; ++i){
         E=EL->flevels[k]->flist[i];
         qtot = E->qtot;

         E->coord(&X);
     
         SPMv[0].shape[part_ID].compute_indicator(qtot,plane, X.x,X.y,ind_func+offset); 
         offset += qtot;
      }
   }

   free(X.x);
   free(X.y);
}



void compute_hydrodynamic_force(Domain *Omega) {
//Uf,Vf,Wf are particle velocity in phsycial space
	double    dt    = dparam("DELT");
  Element   *el, *fl, *gl, *hl, *E;
  int       procid  = option("PROCID");
  int       nprocs  = option("NPROCS");
  int       nz = option("NZ"), nztot = option("NZTOT");
  int       htot = Omega->Uf->htot*nz;
  int       qtot = Omega->Uf->flevels[0]->flist[0]->qtot;
//  double *jacobian_weights= new double [Omega->CONCENTR->htot], *tmp;
  double    *za, *wa, *zb, *wb, *zc, *wc, weight, *r= new double [3];
  int       qa, qb, i, j, k, m,plane, counter, qab,qbc, qt;
  double    lz      = dparam("LZ");
  char      trip    = 'a';

  double **indicator = Omega->vSPM[0].gaussian;
  double ***x_o = Omega->vSPM[0].Xp;
  double ***y_o = Omega->vSPM[0].Yp;
  double ***FHx_o = Omega->vSPM[0].FHx;
  double ***FHy_o = Omega->vSPM[0].FHy;
  double ***NH_o  = Omega->vSPM[0].NH;
  double ***FHphix_o  = Omega->vSPM[0].FHphix;
  double ***FHphiy_o  = Omega->vSPM[0].FHphiy;
//add vy zwang
//  Element_List *Pp=Omega->Pf;
  Coord X;     
  X.x = dvector(0, QGmax*QGmax-1);
  X.y = dvector(0, QGmax*QGmax-1);
  
  int num_partls = Omega->vSPM[0].num_partls;

  double  *fhk[3], *fhpk[3],*nhk[3], *send_buff, *recv_buff;

  double    dz   = lz/(double) nztot;

  int send_len = 9*nz;
  send_buff = dvector(0,send_len-1);

  fhk[0]  = send_buff;
  fhk[1]  = send_buff+nz;
  fhk[2]  = send_buff+2*nz;
  fhpk[0] = send_buff+3*nz;
  fhpk[1] = send_buff+4*nz;
  fhpk[2] = send_buff+5*nz;
  nhk[0]  = send_buff+6*nz;
  nhk[1]  = send_buff+7*nz;
  nhk[2]  = send_buff+8*nz;
  
  int recv_len = send_len*nprocs;  
  recv_buff = dvector(0,recv_len-1);

  double *uback= new double [htot], *vback=new double [htot], *wback=new double [htot];
  double *dpdx= new double [htot], *dpdy=new double [htot], *dpdz=new double [htot];

  //trans back to physical space
  //old particle pressure
//  Omega->Pf->Trans(Omega->Pf, J_to_Q);
//  Omega->Pf->Set_state('p');
//  Omega->Pf->Trans(Omega->Pf, F_to_P);
//  Omega->Pf->Grad_d (dpdx,dpdy,dpdz,'a'); 
  Omega->Pf->Trans(Omega->Pf, J_to_Q);
  Omega->Pf->Grad(NULL,NULL,Omega->Pt,'z');
  Omega->Pf->Set_state('p');
  Omega->Pt->Trans(Omega->Pt, F_to_P);
  dcopy(htot, Omega->Pt->base_h, 1, dpdz, 1);
  Omega->Pf->Trans(Omega->Pf, F_to_P);
 
//  for(int i=0; i<200; ++i)
//    ROOT fprintf(stderr,"Pf = %lf \n", Omega->Pf->base_h[i]);

  Omega->Pf->Grad_h(Omega->Pf->base_h,dpdx,dpdy,NULL,trip);

  double *fh  = new double[3];
  double *fhp = new double[3];
  double *nh  = new double[3];


  int FLAG=0;
  double vsph=0.0, tmp =0.0;
//test
  double *wk=new double [3], test=1;
  
  for(int pid=0; pid < Omega->vSPM[0].num_partls; pid++){

     double  *fx = Omega->vSPM[0].mapx[pid].f;
     double  *fy = Omega->vSPM[0].mapy[pid].f;

     dcopy(htot, Omega->U->base_h, 1, uback, 1);
     dcopy(htot, Omega->V->base_h, 1, vback, 1);
     dcopy(htot, Omega->W->base_h, 1, wback, 1); //keep Omega->U safe , ow has error in calculating the second or more particle
     
     for(int ii=0;ii<htot;ii++) {
        Omega->U->base_h[ii]=indicator[pid][ii]*(Omega->Uf->base_h[ii]-Omega->U->base_h[ii]);
        Omega->V->base_h[ii]=indicator[pid][ii]*(Omega->Vf->base_h[ii]-Omega->V->base_h[ii]);
        Omega->W->base_h[ii]=indicator[pid][ii]*(Omega->Wf->base_h[ii]-Omega->W->base_h[ii]);
      }

   dzero(send_len,send_buff,1);
   dzero(recv_len,recv_buff,1);

   for(k=0,plane=nz*procid; k<nz; ++k,++plane){
     int offset = k*Omega->U->htot;

     double xp = x_o[plane][pid][0];
     double yp = y_o[plane][pid][0];


      fh[0]  =0.0;  fh[1]=0.0;   fh[2]=0.0;
      fhp[0] =0.0;  fhp[1]=0.0;  fhp[2]=0.0;
      nh[0]  =0.0;  nh[1]=0.0;   nh[2]=0.0;

     for (int eid=0; eid<Omega->U->nel; ++eid)  {  
      el = Omega->U->flevels[k]->flist[eid];
      fl = Omega->V->flevels[k]->flist[eid];
      gl = Omega->W->flevels[k]->flist[eid];
      hl = Omega->CONCENTR->flevels[k]->flist[eid];

      el->coord(&X);
      qa = el->qa;   qb = el->qb;
      el->GetZW(&za, &wa, &zb, &wb, &zc, &wc);

      for ( i=0; i<qa; i++ )
      for ( j=0; j<qb; j++ )
        {    
          m = i*qb + j;
          if( el->curvX )
            weight = wb[j]*wa[i]*el->geom->jac.p[m];
          else
            weight = wb[j]*wa[i]*el->geom->jac.d;

          fh[0] += el->h[0][m] * weight;
          fh[1] += fl->h[0][m] * weight;
          fh[2] += gl->h[0][m] * weight;

//          fhp[0] += dpdx[offset+eid*qtot+m] * hl->h[0][m] * weight;
//          fhp[1] += dpdy[offset+eid*qtot+m] * hl->h[0][m] * weight;
//          fhp[2] += dpdz[offset+eid*qtot+m] * hl->h[0][m] * weight;
          fhp[0] += dpdx[offset+eid*qtot+m] * indicator[pid][offset+eid*qtot+m]* weight;
          fhp[1] += dpdy[offset+eid*qtot+m] * indicator[pid][offset+eid*qtot+m]* weight;
          fhp[2] += dpdz[offset+eid*qtot+m] * indicator[pid][offset+eid*qtot+m]* weight;
     
          r[0]= X.x[m] - xp; 
          r[1]= X.y[m] - yp;
          r[2]= 0.;
          

          nh[0] += (r[1]*gl->h[0][m]-r[2]*fl->h[0][m])* weight;
          nh[1] += (r[2]*el->h[0][m]-r[0]*gl->h[0][m])* weight;
          nh[2] += (r[0]*fl->h[0][m]-r[1]*el->h[0][m])* weight;
        
//          if(k == 0){
//          vsph += indicator[pid][offset+eid*qtot+m]* weight;
          vsph += indicator[pid][offset+eid*qtot+m]* weight *dz;
//         }

      }
    }
     fh[0]=-fh[0]/dparam("DELT");
     fh[1]=-fh[1]/dparam("DELT");
     fh[2]=-fh[2]/dparam("DELT");

     nh[0]=-nh[0]/dparam("DELT");
     nh[1]=-nh[1]/dparam("DELT");
     nh[2]=-nh[2]/dparam("DELT");
     
 //no motion is allowed on z direction...
     fhk[0][k] = fh[0];
     fhk[1][k] = fh[1];
     fhk[2][k] = fh[2];

     nhk[0][k]  = nh[0];
     nhk[1][k]  = nh[1];
     nhk[2][k]  = nh[2];

     fhpk[0][k] = fh[0]-fhp[0]/getgamma();
     fhpk[1][k] = fh[1]-fhp[1]/getgamma();
     fhpk[2][k] = fh[2]-fhp[2]/getgamma();

//     fprintf(stdout,"k= %d  plane =%d  pid= %d  fhx= %g fhy= %g fhx_new= %g  fhy_new= %g \n",k,plane,pid,FHx[plane][pid][0],FHy[plane][pid][0],FHphix[plane][pid][0],FHphiy[plane][pid][0]);

   }
    

    MPI_Allgather (&send_buff[0], send_len, MPI_DOUBLE,
         &recv_buff[0], send_len, MPI_DOUBLE, MPI_COMM_WORLD);


    for(int id=0; id<nprocs;++id)
     for(k=0; k<nz;++k){

     plane = id*nz+k;

     FHx_o[plane][pid][1] = FHx_o[plane][pid][0];
     FHy_o[plane][pid][1] = FHy_o[plane][pid][0];
     NH_o[plane][pid][1] = NH_o[plane][pid][0];
     FHphix_o[plane][pid][1] = FHphix_o[plane][pid][0];
     FHphiy_o[plane][pid][1] = FHphiy_o[plane][pid][0];
      
     FHx_o[plane][pid][0] = recv_buff[id*send_len+k]; 
     FHy_o[plane][pid][0] = recv_buff[id*send_len+nz+k]; 
     NH_o[plane][pid][0]  = recv_buff[id*send_len+8*nz+k]; 
     FHphix_o[plane][pid][0] = recv_buff[id*send_len+3*nz+k]; 
     FHphiy_o[plane][pid][0] = recv_buff[id*send_len+4*nz+k]; 

//     fx[plane] = 1.5*FHphix_o[plane][pid][0]-0.5*FHphix_o[plane][pid][1];
//     fy[plane] = 1.5*FHphiy_o[plane][pid][0]-0.5*FHphiy_o[plane][pid][1];
   if(strcmp(Omega->vSPM[0].shape[pid].type,"wall") != 0) //for case 'wall', force from forces() is used
    {
     fx[plane] = FHphix_o[plane][pid][0];
     fy[plane] = FHphiy_o[plane][pid][0];
    }
      
//     ROOT fprintf(stdout,"Fx = %f Fy = %f Fxp = %f Fyp = %f , %d \n",FHx_o[plane][pid][0],FHy_o[plane][pid][0],FHphix_o[plane][pid][0],FHphiy_o[plane][pid][0],plane);
#ifdef MAP // we are in Physcial space
     if(pid >= 1){
       fprintf(stderr,"Errors: SPM and Mapping only work for the case that has one particle !!!!");
       exit(-1);
     }

//     Omega->mapx->f[plane] = 1.5*FHphix_o[plane][pid][0]-0.5*FHphix_o[plane][pid][1];
//     Omega->mapy->f[plane] = 1.5*FHphiy_o[plane][pid][0]-0.5*FHphiy_o[plane][pid][1];
   if(strcmp(Omega->vSPM[0].shape[pid].type,"wall") != 0) //for case 'wall', force from forces() is used
    {
     Omega->mapx->f[plane] = FHphix_o[plane][pid][0];
     Omega->mapy->f[plane] = FHphiy_o[plane][pid][0];
    }

#endif

     if( (strcmp(Omega->vSPM[0].shape[pid].type,"buoy") == 0)
         || (strcmp(Omega->vSPM[0].shape[pid].type,"cir") == 0) )
      {
        int offset = plane*9;
        fx[plane] /= Omega->vSPM[0].shape[pid].data[offset+6];
        fy[plane] /= Omega->vSPM[0].shape[pid].data[offset+6];
#ifdef MAP // we are in Physcial space
        Omega->mapx->f[plane] /= Omega->vSPM[0].shape[pid].data[offset+6];
        Omega->mapy->f[plane] /= Omega->vSPM[0].shape[pid].data[offset+6];
#endif
      }
    }

#ifdef MAP // we are in Physcial space
   if(strcmp(Omega->vSPM[0].shape[pid].type,"wall") != 0) //for case 'wall', force from forces() is used
     {
  #ifdef OLDFFTS
    realft (nztot/2, Omega->mapx->f, -1);
    realft (nztot/2, Omega->mapy->f, -1);
#else
    rfftw(rplan, 1, (FFTW_COMPLEX *) Omega->mapx->f, 1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) Omega->mapy->f, 1, 0, 0, 0, 0);
#endif
     
    int n = iparam ("NFILTER");
     if( plane>(nztot-n) ){  
       Omega->mapx->f[i] = 0.;
       Omega->mapy->f[i] = 0.;
     }
    }
#endif
    dcopy(htot, uback, 1, Omega->U->base_h, 1);
    dcopy(htot, vback, 1, Omega->V->base_h, 1);
    dcopy(htot, wback, 1, Omega->W->base_h, 1);
//    dcopy(Omega->Pf->htot,     Pfback, 1, Omega->Pf->base_h, 1);
   
    gsync ();
 } //end of num_partls


      
    if(dparam("VSPH")==0 )   
     { 
      MPI_Allreduce (&vsph, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dcopy(1,&tmp,1,&vsph,1);
      dparam_set("VSPH", vsph); 
      ROOT fprintf(stdout,"integrating the phi, get VSPH=%lf \n",vsph);
     }
    delete [] r;free(X.x);free(X.y);delete [] uback; delete [] vback; delete [] wback;
    free(fh); free(fhp); free(nh); free(wk);
    delete [] dpdx; delete dpdy; delete dpdz;
    free(send_buff);
    free(recv_buff);

 }
 



void compute_upf(Domain *omega){

  Element *el,*fl,*gl,*hl,*ml;
  int nz = omega->Uf->nz;
  int htot = omega->Uf->htot*nz;

  memset(omega->Uf->base_h,'\0',htot*sizeof(double));
  memset(omega->Vf->base_h,'\0',htot*sizeof(double));
  memset(omega->Wf->base_h,'\0',htot*sizeof(double));

  Coord X;
  X.x = dvector(0, QGmax*QGmax-1);
  X.y = dvector(0, QGmax*QGmax-1);


  double *r= new double [3];
  int ip;
 for (ip = 0; ip < omega->vSPM[0].num_partls; ++ip){

     dcopy(htot,omega->vSPM[0].gaussian[ip],1,omega->Pf->base_h,1);

//  if (strcmp(omega->vSPM[0].shape[ip].type,"cir") == 0){

 for(int k=0,plane=nz*option("PROCID"); k<nz; ++k,++plane){
//      int offset = k*omega->Uf->htot;

      double xp = omega->vSPM[0].Xp[plane][ip][0];
      double yp = omega->vSPM[0].Yp[plane][ip][0];
      double up = omega->vSPM[0].Up[plane][ip][0];
      double vp = omega->vSPM[0].Vp[plane][ip][0];

      double wp = 0.; // tranlational velocity on z direction is zero
      double ang0 = 0.; //angular velocity on x =and y directions are zero
      double ang1 = 0.;

      double ang2 = omega->vSPM[0].Ang[plane][ip][0];
      

    for(int i=0;i<omega->Uf->nel;++i){
       el = omega->Uf->flevels[k]->flist[i];
       fl = omega->Vf->flevels[k]->flist[i];
       gl = omega->Wf->flevels[k]->flist[i];
       hl = omega->CONCENTR->flevels[k]->flist[i];
       ml = omega->Pf->flevels[k]->flist[i];

       el->coord(&X);
  
      for(int m=0;m<el->qtot;++m)  {          
          
          r[0]= X.x[m] - xp; 
          r[1]= X.y[m] - yp;
          r[2]= 0.;
          
          if(hl->h[0][m]!=0) {

          el->h[0][m] += ml->h[0][m]/hl->h[0][m]*(up-(r[1]*ang2-r[2]*ang1));
          fl->h[0][m] += ml->h[0][m]/hl->h[0][m]*(vp-(r[2]*ang0-r[0]*ang2));
          gl->h[0][m] += ml->h[0][m]/hl->h[0][m]*(wp-(r[0]*ang1-r[1]*ang0));
          }
          else {
          el->h[0][m]+=0.0;
          fl->h[0][m]+=0.0;
          gl->h[0][m]+=0.0;
          }
        }
      }//end of Uf->nel
     }
//   } // end of case "cir"
//   else
//   {
//     fprintf(stderr,"compute_upf() Not implemented yet ! \n");
//     exit(-1);
//   }
 }

 omega->Uf->Set_state('p');
 omega->Vf->Set_state('p');
 omega->Wf->Set_state('p');

 free(X.x); 
 free(X.y);
 free(r);
}













