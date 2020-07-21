/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include "nektarF.h"
#ifdef MAP
#include "map.h"
#else
void Compare (Domain *omega, ACTION space);
#endif

#ifdef XMLRPC
#include <unistd.h>
#include "xmlrpc_nektar.h"
#endif

static int    Je;            /* Externals */
static double dt, Re;
#ifdef THERMO
static double Pr;
#endif

static void MakeF   (Domain *omega, ACTION act);
static void StartUp (Domain *, double *, int *);
void solve_w(Domain *omega, int step);

void switch_hj(Element_List *EL, double **hj);
void Bsystem_mem_free(Bsystem *Ubsys, Element_List *U);
void gsync();
void parallel_VdGradV(Domain *omega);
static void noise_filter (Element_List *EL);
void solve(Element_List *U, Element_List *Uf,
	   Bndry **Ubc,Bsystem **Ubsys,SolveType Stype,int step, double *scal);

#ifdef THERMO
void compute_temperature_viscosity(Domain *omega);
void setup_initial_temperature(Domain *omega);
#endif

#ifdef SPM
void compute_indicator_function(Domain *omega, C_SPM *SPMv);
void compute_hydrodynamic_force(Domain *Omega);
void compute_upf(Domain *omega);
void SPM_SetPBCs(Domain *omega);
#ifdef CONJUGATE_HEAT
void regenerate_temperature_system(Domain *omega) ;
#endif
static double **ps;
static double **grad_work;
#endif

int SPM_MAP_step = 0;

#ifdef WALL_MODEL
void SetWallBC(Domain *omega, Bndry **Ubc, int update);
#endif

#ifdef TIMING
#ifdef SYNC_TIMING
#define Timing(s) \
 gsync(); \
 ROOT fprintf(stdout,"%s Took %g cpu seconds\n",s,dclock()-st1); \
 st1 = dclock(); \
 ROOT fprintf(stdout,"%s Took %g wallclock seconds\n",s,MPI_Wtime()-st2); \
 st2 = MPI_Wtime()
#else
#define Timing(s) \
 ROOT fprintf(stdout,"%s Took %g cpu seconds\n",s,dclock()-st1); \
 st1 = dclock(); \
// ROOT fprintf(stdout,"%s Took %g wallclock seconds\n",s,MPI_Wtime()-st2); \
// st2 = MPI_Wtime()
#endif
#else
#define Timing(s) /* Nothing */
#endif

#if defined(VDEBUG) && defined(MAP)
#define Vdebug(myomega,s) \
  ROOT printf("%s\n",s); \
  Compare (myomega, Fourier)
#else
#define Vdebug(myomega,s) /* Nothing */
#endif

#ifdef TRY_ACC
static void Accel_term_implicit (Domain *omega, double dt);
#endif

#ifdef XMLRPC
  xmlrpc_env env;
  xmlrpc_value *resultP_init, 
               *resultP_start;
  xmlrpc_client * clientP;
  xmlrpc_value * action;
  int XMLRPC_BUFFSIZE = 100;

  const char * const serverUrl = "http://localhost:8160";
//  const char * const serverUrl = "http://10.2.91.20:8080";
  const char * const methodName = "request_stochastic_action";
  double current_u1 = 0.0, next_u1 = 0.0;
  double current_u2 = 0.0, next_u2 = 0.0;
  int action_counter = 0;
  void XMLRPC_update_boundary(int step);
#endif

void do_main(int argc, char *argv[]){
  register     int i;
  Domain      *Omega;             /* Solution Domain            */
  int          step, nsteps;      /* Discrete time (step)       */
  double       time = 0.0;        /* Continuous time            */
  ACTION       WaveProp;          /* Convective form            */
  double       begin_clock;
  double       *zerosv = (double *) NULL,
               *varx = zerosv, 
               *vary = zerosv;
#ifdef DEBUG
/*  mallopt(M_DEBUG,1); */
  init_debug();
#endif


  Omega     = PreProcess(argc,argv);
  Vdebug(Omega,"After PreProcess");

  nsteps    = iparam("NSTEPS");
  WaveProp  = (ACTION) iparam("EQTYPE");
  Je        = iparam("INTYPE");
  dt        = dparam("DELT");
  time      = dparam("STARTIME");
  step      = 0;
  Re        = 1./dparam("KINVIS");  
  int hisstep   =iparam("HISSTEP");
#ifdef THERMO
  Pr        = dparam("PRANDTL");   
#endif

  int cflstep = iparam("CFLSTEP");
  if (!cflstep) cflstep = option("hisstep");

#ifdef ZEROOUT
  ROOT {
    dzero(Omega->U->htot, Omega->U->base_h+Omega->U->htot, 1);
    dzero(Omega->V->htot, Omega->V->base_h+Omega->V->htot, 1);
    dzero(Omega->W->htot, Omega->W->base_h+Omega->W->htot, 1);
    dzero(Omega->P->htot, Omega->P->base_h+Omega->P->htot, 1);
  }
#endif

/*
#ifndef MAP
  Nek_Trans_Type f_to_p = F_to_P,
                 p_to_f = P_to_F;

  if (option("dealias")) {
      f_to_p = F_to_P32;
      p_to_f = P_to_F32;
  }

  Omega->U->Trans(Omega->U, f_to_p);
  Omega->V->Trans(Omega->V, f_to_p);
  Omega->W->Trans(Omega->W, f_to_p);
  
  for(i = 0; i < Omega->U->nz; ++i){
    dparam_set("z", zmesh(i));
    Omega->U->flevels[i]->Terror(Omega->soln[0]); 
    Omega->V->flevels[i]->Terror(Omega->soln[1]);
    Omega->W->flevels[i]->Terror(Omega->soln[2]);
  }

  Omega->U->Trans(Omega->U, p_to_f);
  Omega->V->Trans(Omega->V, p_to_f);
  Omega->W->Trans(Omega->W, p_to_f);
#else
#ifdef VDEBUG
  ROOT printf("After PreProcess\n");
#endif
  Compare (Omega, Fourier);
#endif
*/
#ifdef XMLRPC
    /* Initialize our error-handling environment. */

  sleep(300);

  int py_restart_episode = iparam("IPY_RESTART_EPISODE")? iparam("IPY_RESTART_EPISODE"):-1;

  ROOTONLY 
  {
   xmlrpc_env_init(&env);
    /* Create the global XML-RPC client object. */
   xmlrpc_client_init2(&env, XMLRPC_CLIENT_NO_FLAGS, NAME, VERSION, NULL, 0);
   dieIfFaultOccurred(&env);

   fprintf(stderr, "Making XMLRPC call to server url '%s' method '%s' "
           "to request the action \n", serverUrl, methodName);

   resultP_init  = xmlrpc_client_call(&env, serverUrl, "init",
                                 "(i)", py_restart_episode); //restart from Server/save/(-1->31)

   resultP_start = xmlrpc_client_call(&env, serverUrl, "start_episode",
                                 "(i)", -1);
  }
#endif

#ifdef SPM
    ROOT fprintf(stdout,"SPM: preprocessing - start\n");

    int  nztot = option("NZTOT");
    int num_partls;

    Omega->vSPM = new C_SPM[1];
    //set parameters
    if (dparam("INTERF_THICK") == 0)
       Omega->vSPM[0].buffer_width = 3.0/(iparam("MODES")-1);
    else
       Omega->vSPM[0].buffer_width = dparam("INTERF_THICK");

    num_partls = iparam("NUM_PARTLS");

    ROOT fprintf(stdout,"SPM: num_parts = %d\n",num_partls);

    int hjtot = Omega->U->hjtot*Omega->U->nz;
    int htot  = Omega->U->htot*Omega->U->nz;

    Omega->vSPM[0].SPM_init(num_partls,htot);
    Omega->vSPM[0].read_init_conditions(Omega->name);
    Omega->vSPM[0].read_in_particle_geometry(Omega->name);

#if defined(SPM) && defined(MAP)
    Omega->tens = dvector(0,nztot-1);
    Omega->stif = dvector(0,nztot-1);

    dcopy(nztot, Omega->vSPM[0].stif[0], 1, Omega->stif, 1);
    dcopy(nztot, Omega->vSPM[0].tens[0], 1, Omega->tens, 1);
#endif

#if defined(SPM) && !defined(MAP)
    Omega->vSPM[0].create_file_particle_out(Omega->name);
#endif

//    Omega->vSPM[0].ReadMap (Omega->name);
    
    ps = (double **) malloc(2*sizeof(double*));
    for(int dd=0; dd<2; ++dd)
     {
      ps[dd] = dvector(0, hjtot-1);
      dzero(hjtot, ps[dd], 1);
     }

    grad_work = (double **) malloc(7*sizeof(double*));
    for(int dd=0; dd<7; ++dd)
      grad_work[dd] = dvector(0, htot-1);
//generate SPM first
    compute_indicator_function(Omega, Omega->vSPM);
    dzero(htot, Omega->CONCENTR->base_h, 1); 
    for (int i = 0; i < Omega->vSPM[0].num_partls; ++i)
       dvadd(htot,Omega->vSPM[0].gaussian[i],1, Omega->CONCENTR->base_h, 1, Omega->CONCENTR->base_h, 1);
    //store the value in Q,P space to grad_work[6]
    dcopy(htot, Omega->CONCENTR->base_h, 1, grad_work[6], 1); 
    Omega->CONCENTR->Set_state('p');
//    //Trans for convective term
//    omega->CONCENTR->Trans(CONCENTR,Q_to_F);

//#ifdef CONJUGATE_HEAT
//    if(dparam("STARTIME") < 1e-8)
//    setup_initial_temperature(Omega);
//    regenerate_temperature_system(Omega);
//#endif

    ROOT fprintf(stdout,"SPM: preprocessing - done\n");
#endif

  StartUp (Omega,&time,&step);
  Vdebug(Omega,"After StartUp");
  /*------------------------------------------------------------------*
   *                    Main Integration Loop                         *
   *------------------------------------------------------------------*/
  begin_clock = MPI_Wtime();
#ifdef TIMING
  double st1, st2;
  st1 = dclock();
  st2 = MPI_Wtime();
#endif
  while (step < nsteps) {
    set_order(((step+1 < Je)? step+1: Je));
    dparam_set("t", time+dt); 
    MakeF (Omega, Prep);
    Timing("Prep.......");

    ROOT // approximate CFL number based on average flow
      if ((step % cflstep) == 0){
	double cfl = cfl_checker(Omega,dt);
	printf("Time %lf CFL = %lf\n", time, cfl);
	if (isnan(cfl)) exit(-1);
      }

    MakeF (Omega, WaveProp);
//    parallel_VdGradV(Omega);
    Timing("WavePropo..");
#ifdef XMLRPC
    XMLRPC_update_boundary(step);
#endif

  if(option("tvarying")){
    TransformBCs(Omega, F_to_P);
    Bndry *Bc;
    for(int k = 0; k < Omega->U->nz; ++k) 
//      if(parid(k)!=1)
   {
	for(Bc=Omega->Ubc[k];Bc;Bc=Bc->next)
	  Bc->elmt->update_bndry(Bc,1);
	for(Bc=Omega->Vbc[k];Bc;Bc=Bc->next)
	  Bc->elmt->update_bndry(Bc,1);
	for(Bc=Omega->Wbc[k];Bc;Bc=Bc->next)
	  Bc->elmt->update_bndry(Bc,1);
   }
    // Get the BCs in Fourier Space
    TransformBCs(Omega, P_to_F);
  }
    // Integration for first step of splitting scheme
    Integrate_SS(Je, dt, Omega->U, Omega->Uf, Omega->u, Omega->uf);
    Integrate_SS(Je, dt, Omega->V, Omega->Vf, Omega->v, Omega->vf);
    Integrate_SS(Je, dt, Omega->W, Omega->Wf, Omega->w, Omega->wf);
#ifdef THERMO
    Integrate_SS(Je, dt, Omega->T, Omega->Tf, Omega->Tt,Omega->tf);
#endif
    Timing("Integrate..");
    Vdebug(Omega,"After Integrate");
    
    // Solve pressure system for second step of splitting scheme
    MakeF (Omega, Pressure);
    Timing("Pressure...");
    Vdebug(Omega,"After MakeF Pressure");

#ifdef NOISE_FILTER
    for(i = 0; i < Omega->P->nz; ++i)
      noise_filter (Omega->P->flevels[i]);
    Timing("noise_filter....");
    Vdebug(Omega,"After noise_filter");
#endif

    solve (Omega->P,Omega->Uf,Omega->Pbc,Omega->Pressure_sys,Helm,step,zerosv);
    Timing("P_solve....");
    Vdebug(Omega,"After P_solve");
#ifdef SPM
    dcopy(hjtot, Omega->P->base_hj, 1, ps[0], 1);
#endif
    //    Omega->P->Trans     (Omega->P, Q_to_J); 

    // Solve velocity system for third step of splitting scheme
    MakeF (Omega, Viscous);
#ifdef ENTROPYVISCOSITY
   int ev_step = (dparam("EV_STEP") ? dparam("EV_STEP") : 2);
    if(step > ev_step)
     {
      compute_entropy_viscosity(Omega);
#ifdef THERMO
     compute_temperature_viscosity(Omega);
#endif
     }
#endif
    Timing("Viscous....");

#if defined(MAP) && defined(DIRNEW)
    varx = Omega->mapx->t;
    vary = Omega->mapy->t;
 #endif    
//zwang 20170105  
#ifdef WALL_MODEL
  if( option("wall_model") )
   {
    int update = (int) dparam("AXI");

    Bndry **Ubcs[3];
    Ubcs[0] = Omega->Ubc;
    Ubcs[1] = Omega->Vbc;
    Ubcs[2] = Omega->Wbc;

    SetWallBC(Omega,Ubcs[update],update);
   }
#endif
    solve (Omega->U,Omega->Uf,Omega->Ubc,Omega->Usys,Helm,step,varx);
    solve (Omega->V,Omega->Vf,Omega->Vbc,Omega->Vsys,Helm,step,vary);

#ifndef JAMES
    solve (Omega->W,Omega->Wf,Omega->Wbc,Omega->Wsys,Helm,step,zerosv);
#else
    solve_w(Omega,step);
#endif

#ifdef THERMO
    solve (Omega->T,Omega->Tf,Omega->Tbc,Omega->Tsys,Helm,step,zerosv);
#endif

    Timing("V_solve....");
    Vdebug(Omega,"After V_solve");


 // ROOT for (int k = 0; k <Omega->mapx->NZ; k++) fprintf (stdout, "%lf \n", Omega->mapx->z[k]);
#ifdef SPM 
    MakeF (Omega, SPM_Prep); //trans to physical space for u*
    Timing("SPM Prep  .......");

    if(iparam("MOVING_PARTLS")==1)
     {
       Omega->vSPM[0].UpdateMap(Omega->name);
       Timing("Updating Map .......");
     }
        //LG  make sure that compute_indicator_function is only called for (option("SPM")==1)  and for mobile particle
        //LG  consider changing it to update_"indicator_function"
  if( (step == 0) || (iparam("MOVING_PARTLS")==1) )
   {
//    compute_indicator_function(Omega, Omega->vSPM);
////    memcpy(Omega->CONCENTR->base_h,Omega->vSPM[0].gaussian[0],Omega->U->htot*Omega->U->nz*sizeof(double));
//    dzero(Omega->U->htot*Omega->U->nz, Omega->CONCENTR->base_h, 1); 
//    for (int i = 0; i < Omega->vSPM[0].num_partls; ++i)
//       dvadd(Omega->U->htot*Omega->U->nz,Omega->vSPM[0].gaussian[i],1, Omega->CONCENTR->base_h, 1, Omega->CONCENTR->base_h, 1);
//
//    Omega->CONCENTR->Set_state('p');
//    
//    Timing("Updating Concentration .......");
   //restore the (Q,P) value from grad_work[6]
    dcopy(htot, grad_work[6], 1, Omega->CONCENTR->base_h, 1);

    compute_upf(Omega);

    Timing("SPM upf .......");

    dcopy(htot, Omega->Uf->base_h, 1, grad_work[3], 1);
    dcopy(htot, Omega->Vf->base_h, 1, grad_work[4], 1);
    dcopy(htot, Omega->Wf->base_h, 1, grad_work[5], 1);

    //fix for the first step,ps[1] needs some pressure
    dcopy(hjtot, Omega->P->base_hj, 1, ps[1], 1);
   }
  else
   {
      dcopy(htot, grad_work[3], 1, Omega->Uf->base_h, 1);
      dcopy(htot, grad_work[4], 1, Omega->Vf->base_h, 1);
      dcopy(htot, grad_work[5], 1, Omega->Wf->base_h, 1);
   }

        //calculate up^n field,stored in Omega->Uf, Vf, Wf
      //integrate phi(Up^n-u*) to get F^H[i], hydro force for particle[i]
    /*Pf is safe to use*/
#ifndef MAP
  if ( (step % hisstep == 0) || ( (iparam("MOVING_PARTLS")) == 1 )) 
#endif
   if(strcmp(Omega->vSPM[0].shape[0].type,"wall") != 0) //for case 'wall', force from forces() is used
   {

    dcopy(hjtot, ps[1], 1, Omega->Pf->base_hj, 1);
    compute_hydrodynamic_force(Omega);
    Timing("SPM Compute Forces ........");
   }


#if 1
    MakeF (Omega, SPM_Pressure);
    Timing("SPM Pressure..");      
    SPM_SetPBCs(Omega);
    solve (Omega->P,Omega->Pf,Omega->Pbc,Omega->Pressure_sys,Helm,step,zerosv);
    dcopy(hjtot, Omega->P->base_hj, 1, ps[1], 1);
    Timing("SPM Solve...");

    MakeF (Omega, SPM_Post);
    Timing("SPM Postprocess ........");
#endif

#endif

    // "Back out" formula for new velocity values
    MakeF (Omega, Post);
    Timing("Post.......");

    Analyser(Omega, ++step, (time += dt));
    Timing("Analyser...");
    Vdebug(Omega,"After Analyzer");
  }
  printf("Wallclock time of solver (seconds):  %lf \n", MPI_Wtime()-begin_clock);

  #if defined(TEST) && !defined(VDEBUG)
  Compare (Omega, Fourier);
  #endif

  PostProcess(Omega, step, time);
  
#ifdef SPM
  ROOT for(int i=0; i<Omega->vSPM[0].num_partls; i++)
         fclose(Omega->vSPM[0].file_particleout[i]);
#endif

  return;
}
  
void solve(Element_List *U, Element_List *Uf, 
	   Bndry **Ubc, Bsystem **Ubsys,SolveType Stype, int step, 
	   double *scal)
{
  int k,i;
#ifdef DIRNEW
  static double *saveRHS = (double *) NULL;
  if (!saveRHS && U->fhead->type == 'u' && Ubc[0]->DirRHS) {
    saveRHS = dvector(0,Ubsys[0]->nsolve-1);
    ROOT dcopy (Ubsys[0]->nsolve, Ubc[0]->DirRHS, 1, saveRHS, 1);
    // distribute to all other processors
    MPI_Bcast (saveRHS, Ubsys[0]->nsolve, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif  
  
  if(step && step < Je  && U->fhead->type == 'u'){
     for(k = 0; k < U->nz; ++k){ 
       for(i=0;i<U->nel;++i)
	 Ubsys[k]->lambda[i].d  = Beta(k)*Beta(k)+Re*getgamma()/dt;
       ROOT fprintf(stdout,"Level: %d \n", k);
       
       if(!(parid(k) % 2)) {
	 Bsystem_mem_free(Ubsys[k], U->flevels[k]);
	 
	 ROOT fprintf(stdout,"Generating velocity system [."); fflush(stdout);
	 GenMat (U->flevels[k],Ubc[k],Ubsys[k],Ubsys[k]->lambda,Helm);
	 ROOT fprintf(stdout,"]\n");
       }
       else{
	 Ubsys[k]->Gmat = Ubsys[k-1]->Gmat;
	 Ubsys[k]->Pmat = Ubsys[k-1]->Pmat;
	 Ubsys[k]->rslv = Ubsys[k-1]->rslv;
       }
     }
  }

#ifdef THERMO
  Pr = dparam("PRANDTL");
//   double kr = dparam("RATIO_KAPPA") == 0? 1.:dparam("RATIO_KAPPA");
//   double kappa = dparam("DKAPPA") == 0? 1./Re:dparam("DKAPPA");

  if(step && step < Je  && U->fhead->type == 't'){
//  fprintf(stderr,"Pr = %g  Re = %g \n",Pr,Re);
     for(k = 0; k < U->nz; ++k){ 
       int nid = 0;
       for(i=0;i<U->nel;++i) {
       double Pr_tmp = Pr;
//       double Pr_tmp = Pr*kr;
//     #ifdef CONJUGATE_HEAT
//      int group_id =  U->flevels[k]->flist[i]->group_id; 
//      if(group_id > 1)
//       {
//        nid++;
//        Pr_tmp = Pr/kr;
////        ROOT fprintf(stdout,"Solid Kappa = %g \n",Pr_tmp);
//       }
//
//     #endif
    
//	 Ubsys[k]->lambda[i].d  = Beta(k)*Beta(k)+Pr_tmp*Re/kr*getgamma()/dt;
	 Ubsys[k]->lambda[i].d  = Beta(k)*Beta(k)+Re*Pr*getgamma()/dt;
 }
//       fprintf(stdout,"2 Number of solid elements = %d \n",nid);
       ROOT fprintf(stdout,"Level: %d \n", k);
       
       if(!(parid(k) % 2)) {
	 Bsystem_mem_free(Ubsys[k], U->flevels[k]);
	 
	 ROOT fprintf(stdout,"Generating temperature system [."); fflush(stdout);
	 GenMat (U->flevels[k],Ubc[k],Ubsys[k],Ubsys[k]->lambda,Helm);
	 ROOT fprintf(stdout,"]\n");
       }
       else{
	 Ubsys[k]->Gmat = Ubsys[k-1]->Gmat;
	 Ubsys[k]->Pmat = Ubsys[k-1]->Pmat;
	 Ubsys[k]->rslv = Ubsys[k-1]->rslv;
       }
     }
  }
#endif


  if(U->fhead->type != 'p' && step && step < Je)
    for(k = 0; k < U->nz; ++k){
      free(Ubc[k]->DirRHS); Ubc[k]->DirRHS = (double*) 0;
#ifdef DIRNEW
      if(!option("tvarying")) {
	DirBCs(U->flevels[k], Ubc[k], Ubsys[k], Helm);
	if ((k == 0) && (U->fhead->type == 'u')) {
	  ROOT dcopy (Ubsys[0]->nsolve, Ubc[0]->DirRHS, 1, saveRHS, 1);
	  // distribute to all other processors
	  MPI_Bcast (saveRHS, Ubsys[0]->nsolve, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
      }
#endif
    }

#ifdef DIRNEW
  if(!option("tvarying")) {
    if (U->fhead->type == 'u') {
      int kstart = 0;
      if (option("PROCID") == 0) {
	dsmul(Ubsys[0]->nsolve, 1.0-scal[0], saveRHS, 1, 
	      Ubc[0]->DirRHS, 1);
	kstart++;
      }
      for(k = kstart; k < U->nz; ++k)
	dsmul(Ubsys[0]->nsolve, -scal[parid(k)], saveRHS, 1, 
	      Ubc[k]->DirRHS, 1);
    } else if (U->fhead->type == 'v')
      for(k = 0; k < U->nz; ++k)
	dsmul(Ubsys[0]->nsolve, -scal[parid(k)], saveRHS, 1, 
	      Ubc[k]->DirRHS, 1);
  }
#endif
// zwang 05182020  time varying boundary condition
//#ifndef DIRNEW
//  if(option("tvarying")){
//    Bndry *Bc;
//    for(k = 0; k < U->nz; ++k) 
//      if(parid(k)!=1)
//	for(Bc=Ubc[k];Bc;Bc=Bc->next)
//	  Bc->elmt->update_bndry(Bc,0);
//  }
//#endif

  double soltol = (U->fhead->type == 'p' ? dparam("TOLFIX") : 0.);
  double rhs;
  for(k=0; k < U->nz; ++k){
    if(parid(k)!=1){
      SetBCs (U->flevels[k],Ubc[k],Ubsys[k]);
      Solve  (U->flevels[k],Uf->flevels[k],Ubc[k],Ubsys[k],Stype);
#ifdef NOISE_FILTER
      noise_filter (U->flevels[k]);
#endif
    }
    else {
      dzero(U->hjtot, U->base_hj+U->hjtot, 1);
      U->flevels[k]->Set_state('t');
    }
  }
}

static void noise_filter (Element_List *EL)
{
  Element *E;
  Nek_Facet_Type id;
  double vfilter = (EL->fhead->type == 'p' ? dparam("FILTERP") : dparam("FILTERV"));

  register double val, temp;
  register int i, j;

  if (EL->fhead->state == 't') { // Modal space
    for (E = EL->fhead; E; E=E->next) {
      id = E->identify();

      // find the maximum vertex value in the Element
      val = 0.;
      for (i = 0; i < E->Nverts; i++)
	val = max(val,fabs(E->vert[i].hj[0]));

      if (val > vfilter) 
	val *= vfilter;
      else
	val = vfilter;

      // filter out 
      for (i = 0; i < E->Nverts; i++) {
	temp = fabs(E->vert[i].hj[0]);
	if (temp < val)
	  E->vert[i].hj[0] = 0.;
      }

      for (i = 0; i < E->Nedges; i++) {
	// find the maximum edge value in the Element
	val = 0.;
	for (j = 0; j < E->edge[i].l; j++)
	  val = max(val,fabs(E->edge[i].hj[j]));

	if (val > vfilter) 
	  val *= vfilter;
	else
	  val = vfilter;

	// filter out 
 	for (j = 0; j < E->edge[i].l; j++) {
	  temp = fabs(E->edge[i].hj[j]);
	  if (temp < val)
	    E->edge[i].hj[j] = 0.;
	}
      }
	
      if (id == Nek_Tri) { // Triangle
	// find the maximum face value in the Element
	val = 0.;
	for (i = 0; i < E->face[0].l; i++)
	  for (j = 0; j < E->face[0].l-i-1; j++)
	    val = max(val,fabs(E->face[0].hj[i][j]));

	if (val > vfilter) 
	  val *= vfilter;
	else
	  val = vfilter;

	// filter out 
	for (i = 0; i < E->face[0].l; i++)
	  for (j = 0; j < E->face[0].l-i-1; j++) {
	    temp = fabs(E->face[0].hj[i][j]);
	    if (temp < val)
	      E->face[0].hj[i][j] = 0.;
	  } 
      } else { // Quad
	// find the maximum face value in the Element
	val = 0.;
	for (i = 0; i < E->face[0].l; i++)
	  for (j = 0; j < E->face[0].l; j++)
	    val = max(val,fabs(E->face[0].hj[i][j]));

	if (val > vfilter) 
	  val *= vfilter;
	else
	  val = vfilter;

	// filter out 
	for (i = 0; i < E->face[0].l; i++)
	  for (j = 0; j < E->face[0].l; j++) {
	    temp = fabs(E->face[0].hj[i][j]);
	    if (temp < val)
	      E->face[0].hj[i][j] = 0.;
	  }
      }   
    }
  } else { // Quadrature Space
    for (E = EL->fhead; E; E=E->next) {

      // find the maximum vertex value in the Element
      val = 0.;
      for (i = 0; i < E->Nverts; i++)
	val = max(val,fabs(E->vert[i].hj[0]));

      if (val > vfilter) 
	val *= vfilter;
      else
	val = vfilter;

      // filter out 
      for (i = 0; i < E->Nverts; i++) {
	temp = fabs(E->vert[i].hj[0]);
	if (temp < val)
	  E->vert[i].hj[0] = 0.;
      }

      // find the maximum quadrature value in the Element
      val = 0.;
      for (i = 0; i < E->qb; i++)
	for (j = 0; j < E->qa; j++)
	val = max(val,fabs(E->h[i][j]));

      if (val > vfilter) 
	val *= vfilter;
      else
	val = vfilter;
 
      // filter out 
      for (i = 0; i < E->qb; i++)
	for (j = 0; j < E->qa; j++) {
	  temp = fabs(E->h[i][j]);
	  if (temp < val)
	    E->h[i][j] = 0.;
	}
    }
  }
  return;
}

static void MakeF(Domain *omega, ACTION act){
  Element_List  *U    =  omega->U,  *V    =  omega->V,   *W    =  omega->W,
                *Uf   =  omega->Uf, *Vf   =  omega->Vf,  *Wf   =  omega->Wf,
                *P    =  omega->P;
 #ifdef THERMO
  Element_List *T     = omega->T,   *Tf = omega->Tf;
 #endif
 #ifdef SPM
  Element_List *CONCENTR     = omega->CONCENTR;
 #endif

  const    int  nel   =  U->nel;
  int           k;
  int           hjtot = U->nz*U->hjtot;
  int           htot = U->nz*U->htot;

  switch (act) {

  case Prep: /* put fields in physical space for waveprop */
#ifndef JAMES
    U->Trans(U,J_to_Q);
    V->Trans(V,J_to_Q);
    W->Trans(W,J_to_Q);
#ifdef THERMO
    T->Trans(T,J_to_Q);
#endif
//#ifdef SPM
//    CONCENTR->Trans(CONCENTR,J_to_Q);
//#endif
#else
    U->Set_state('p');
    V->Set_state('p');
    W->Set_state('p');
#ifdef THERMO
    T->Set_state('p');
#endif
//#ifdef SPM
//    CONCENTR->Set_state('p');
//#endif

#endif
    
//#if defined(MAP) && defined(DAVE)
//     UpdateMap (omega);
//     UpdateVbc (omega);
//#endif

    break;
    
  case Rotational:
    VxOmega (omega);
    goto AddForcing; 

  case Convective:
    VdgradV (omega); 
//     SkewSymmetric(omega);
    goto AddForcing; 

  case Stokes:
    StokesBC (omega);

  AddForcing: {
    int ntot;
    ROOT{
      double fx = dparam("FFX");
      double fy = dparam("FFY");
      double fz = dparam("FFZ");
 #ifdef THERMO
      double ft = dparam("FFT");
 #endif
      if(fx) dsadd(Uf->htot, fx, Uf->base_h, 1,Uf->base_h, 1);
      if(fy) dsadd(Uf->htot, fy, Vf->base_h, 1,Vf->base_h, 1);
      if(fz) dsadd(Uf->htot, fz, Wf->base_h, 1,Wf->base_h, 1);
 #ifdef THERMO
      if(ft) dsadd(Tf->htot, ft, Tf->base_h, 1,Tf->base_h, 1);

// #if defined(SPM) && defined(CONJUGATE_HEAT) 
//      dvadd(Tf->htot, omega->CONCENTR->base_h, 1, Tf->base_h, 1,Tf->base_h, 1); //dimensionless body heat generation
// #endif

 #endif

    if( dparam("ADJUSTFORCE") )
      {
       double ff = omega->adjust_force;
       if(dparam("PIPE"))
         dsadd(Wf->htot, ff, Wf->base_h, 1, Wf->base_h, 1);
       if(dparam("CHANNEL"))
         dsadd(Uf->htot, ff, Uf->base_h, 1, Uf->base_h, 1);
      }
      // The extra mapping terms are added in the advection phase
    }
    /* white-noise Forcing */
       
 if(omega->ForceFuncs){
	if(option("tvarying")) {
#if 1
	  int i;
     dvadd(htot, omega->ForceFuncs[0], 1, Uf->base_h, 1, Uf->base_h, 1);
	  dvadd(htot, omega->ForceFuncs[1], 1, Vf->base_h, 1, Vf->base_h, 1);
	  dvadd(htot, omega->ForceFuncs[2], 1, Wf->base_h, 1, Wf->base_h, 1);
	  omega->Pf->Set_field(omega->ForceStrings[0]);
	  dcopy(htot, omega->Pf->base_h, 1, omega->ForceFuncs[0], 1);
	  omega->Pf->Set_field(omega->ForceStrings[1]);
	  dcopy(htot, omega->Pf->base_h, 1, omega->ForceFuncs[1], 1);
	  omega->Pf->Set_field(omega->ForceStrings[2]);
	  dcopy(htot, omega->Pf->base_h, 1, omega->ForceFuncs[2], 1);
 
#else
	  omega->Pf->Set_field(omega->ForceStrings[0]);
	  daxpy(Uf->htot, 1.0, omega->Pf->base_h, 1, Uf->base_h, 1);
	  omega->Pf->Set_field(omega->ForceStrings[1]);
	  daxpy(Uf->htot, 1.0, omega->Pf->base_h, 1, Vf->base_h, 1);
	  omega->Pf->Set_field(omega->ForceStrings[2]);
	  daxpy(Uf->htot, 1.0, omega->Pf->base_h, 1, Wf->base_h, 1);
#endif
	}
	else {
	  dvadd(Uf->htot, omega->ForceFuncs[0], 1, Uf->base_h, 1, Uf->base_h, 1);
	  dvadd(Uf->htot, omega->ForceFuncs[1], 1, Vf->base_h, 1, Vf->base_h, 1);
	  dvadd(Uf->htot, omega->ForceFuncs[2], 1, Wf->base_h, 1, Wf->base_h, 1);
	}
      } 
    if (option("RAND")){
      Element *e;
      int knoise = iparam("KNOISE"),m;
      int nz = Uf->nz, procid = option("PROCID");
      double f;
      int i, j = max(1,procid*nz/2);

      for (m = 2*j-procid*nz; j <= knoise && m < nz; ++j, m = 2*j-procid*nz){
	if((f = 2.*dparam("XNOISE")) > 0.0){
	  for(k = 0; k < nel; ++k){
	    e = Uf->flevels[m]->flist[k];
	    ntot = e->qa*e->qb;
	    for(i = 0; i < ntot; ++i) e->h[0][i] += f*(drand() - .5);
	    e = Uf->flevels[m+1]->flist[k];
	    ntot = e->qa*e->qb;
	    for(i = 0; i < ntot; ++i) e->h[0][i] += f*(drand() - .5);
	  }
	}
	if((f = 2.*dparam("YNOISE")) > 0.0){
	  for(k = 0; k < nel; ++k){
	    e = Vf->flevels[m]->flist[k];
	    ntot = e->qa*e->qb;
	    for(i = 0; i < ntot; ++i) e->h[0][i] += f*(drand() - .5);
	     e = Vf->flevels[m+1]->flist[k];
	    ntot = e->qa*e->qb;
	     for(i = 0; i < ntot; ++i) e->h[0][i] += f*(drand() - .5);
	    }
	}
	if((f = 2.*dparam("ZNOISE")) > 0.0){
	  for(k = 0; k < nel; ++k){
	    e = Wf->flevels[m]->flist[k];
	    ntot = e->qa*e->qb;
	    for(i = 0; i < ntot; ++i) e->h[0][i] += f*(drand() - .5);
	    e = Wf->flevels[m+1]->flist[k];
	    ntot = e->qa*e->qb;
	    for(i = 0; i < ntot; ++i) e->h[0][i] += f*(drand() - .5);
	    }
	}
      }
      
       if(j=option("RAND"))
	option_set("RAND",j-1);
    }
    
    SetPBCs(omega);
    break;
  }

  case Pressure: { 
    double dtinv = 1./dt;
    double *nul  = (double*)0;
    double *ux   = dvector(0,3*htot-1);
    double *vy   = ux+htot;
    double *wz   = ux+2*htot;
    double *divv = omega->Uf->base_h;

    dzero(3*htot,ux,1);
    W->Set_state('p');
    U->Grad_d( ux, nul, nul, 'x');
    V->Grad_d(nul,  vy, nul, 'y');
    W->Grad_d(nul, nul,  wz, 'z');

    dvadd(htot, ux, 1, vy, 1, divv, 1);
    dvadd(htot, wz, 1, divv, 1, divv, 1);

    Uf->Set_state('p');
    Uf->Iprod(Uf);
    dscal(hjtot, dtinv, Uf->base_hj, 1);

    Uf->Set_state('t');
    free(ux);
// #if defined(MAP) && !defined(DAVE)
#ifdef MAP
//    if( SPM_MAP_step > iparam("ISPM_MAP_STEP") )
//     {
      UpdateMap (omega);
      UpdateVbc (omega);
    
#ifdef TRY_ACC
     if (option("acc_impl")) {
      // Add pseudo acceleration
       Accel_term_implicit (omega, dt);
       AddAccelPBCs(omega);
       }
#endif
//     }
//      SPM_MAP_step++;
#endif // defined(MAP) && !defined(DAVE)

    break;
  }
    
  case Viscous: {
    double dtinv = 1/dt;
 
    P->Trans     (P, J_to_Q); 

    P->Set_state ('p');
    P->Grad      (Uf,Vf,Wf,'a');
    P->Set_state ('t');

    daxpy(htot, -dtinv, U->base_h, 1, Uf->base_h, 1); 
    Uf->Iprod(Uf);    
    dscal(hjtot,    Re, Uf->base_hj, 1);
#ifndef SAVINGSD
    dcopy(hjtot, omega->us, 1, omega->U->base_hj, 1);
#endif

    daxpy(htot, -dtinv, V->base_h, 1, Vf->base_h, 1); 
    Vf->Iprod(Vf);
    dscal(hjtot,    Re, Vf->base_hj, 1);
#ifndef SAVINGSD
    dcopy(hjtot, omega->vs, 1, omega->V->base_hj, 1);
#endif

#ifdef JAMES
    ROOT{
      daxpy(W->htot, -dtinv, W->base_h, 1, Wf->base_h, 1); 
      Wf->flevels[0]->Iprod(Wf->flevels[0]);
      dscal(W->hjtot,    Re, Wf->base_hj, 1);    
#ifndef SAVINGSD
      dcopy(W->hjtot, omega->ws, 1, omega->W->base_hj, 1);
#endif
    }
#else
    daxpy(htot, -dtinv, W->base_h, 1, Wf->base_h, 1); 

    Wf->Iprod(Wf);
    dscal(hjtot,    Re, Wf->base_hj, 1);
#ifndef SAVINGSD
    dcopy(hjtot, omega->ws, 1, omega->W->base_hj, 1);
#endif
#endif

#ifdef THERMO   
    Pr = dparam("PRANDTL");
//    double kr = dparam("RATIO_KAPPA") == 0? 1.:dparam("RATIO_KAPPA");
//    double kappa = dparam("DKAPPA") == 0? 1./Re:dparam("DKAPPA");
//    double Pr = dparam("D") == 0? 1./Re:dparam("DKAPPA");

    dzero(htot, Tf->base_h, 1);
    daxpy(htot, -dtinv, T->base_h, 1, Tf->base_h, 1); 

    Tf->Set_state('p');
    Tf->Iprod(Tf);    
 
//    double kr = dparam("RATIO_KAPPA") == 0? 1.:dparam("RATIO_KAPPA");
//    double Pr_tmp = Pr*kr;
//#ifdef CONJUGATE_HEAT
//    Element *el;
//    
//    for(int k=0; k<Tf->nz; ++k)
//    {
//      int nid=0;
//     for(int eid=0; eid<nel; ++eid)
//      {
//       el = Tf->flevels[k]->flist[eid];
//       int group_id =  el->group_id; 
//
//
//        if(group_id > 1)
//        {
//         Pr_tmp = Pr/kr;
//         nid++;
//        }
//
//         dscal(el->Nmodes, Re*Pr_tmp, Tf->base_hj+k*Tf->hjtot+eid*el->Nmodes, 1);
//      }
//       fprintf(stdout,"4 Number of solid elements = %d \n",nid);
//
//    }
//#else
//    dscal(hjtot,  Re*Pr/kr, Tf->base_hj, 1);
    dscal(hjtot,  Re*Pr, Tf->base_hj, 1);
//    dscal(hjtot,  Re*Pr_tmp, Tf->base_hj, 1);
//#endif

	  Tf->Set_state('t');
#endif
    
    break;
  }
  case Post: {
#ifndef SAVINGSD
    dcopy(hjtot, U->base_hj, 1, omega->us, 1);
    dcopy(hjtot, V->base_hj, 1, omega->vs, 1);
    
    //fixed for james
#ifdef JAMES
    ROOT
      dcopy(W->hjtot, W->base_hj, 1, omega->ws, 1);
#else
    dcopy(hjtot, W->base_hj, 1, omega->ws, 1);
#ifdef THERMO
    dcopy(hjtot, T->base_hj, 1, omega->ts, 1);
#endif

#endif
#endif

    break;
  }
#ifdef SPM  
  case SPM_Prep: /* put fields in physical space for waveprop */
//#ifndef JAMES
    U->Set_state('p');
    V->Set_state('p');
    W->Set_state('p');
//    CONCENTR->Set_state('p');

    U->Trans(U,J_to_Q);
    U->Trans(U, F_to_P); 
    V->Trans(V,J_to_Q);
    V->Trans(V, F_to_P); 
    W->Trans(W,J_to_Q);
    W->Trans(W, F_to_P); // here we are in phsyical/quadrature space 
//    CONCENTR->Trans(CONCENTR, F_to_P); // here we are in phsyical/quadrature space 
//#ifdef THERMO
//    T->Trans(T,J_to_Q);
//    T->Trans(T, F_to_P); // here we are in phsyical/quadrature space
//#endif
//#else
//#ifdef THERMO
//    T->Set_state('p');
//#endif

//#endif
    
    break;
    case SPM_Pressure: { 
    double    dtinv = 1./dt;
    double    *dudx   = grad_work[0];
    double    *dvdy   = grad_work[1];
    double    *dwdz   = grad_work[2];

    double *nul= (double*)0;
    double *divU = omega->Pf->base_h;
    double *phi= omega->CONCENTR->base_h;
//    MakeF (omega, Prep, gaussian, Up, Vp, Wp); //recover U*

    for(int i=0;i<htot;i++) {
    Uf->base_h[i]=phi[i]*(Uf->base_h[i]-U->base_h[i]);
    Vf->base_h[i]=phi[i]*(Vf->base_h[i]-V->base_h[i]);
    Wf->base_h[i]=phi[i]*(Wf->base_h[i]-W->base_h[i]);
    }    
    
    Uf->Set_state('p');
    Vf->Set_state('p');
    Wf->Set_state('p');

    Uf->Trans(Uf, P_to_F); //go to Fourier space 
    Vf->Trans(Vf, P_to_F); //go to Fourier space 
    Wf->Trans(Wf, P_to_F); //go to Fourier space 
    
    Uf->Grad_d( dudx,  nul,  nul, 'x');
    Vf->Grad_d(  nul, dvdy,  nul, 'y');
    dvadd(htot, dudx, 1, dvdy, 1, divU, 1);

 //   Wf->Set_state('p');
    Wf->Grad_d(  nul,  nul,  dwdz, 'z');
    dvadd(htot, dwdz, 1, divU, 1, divU, 1);
    
    dscal(htot,getgamma(),divU,1);
    omega->Pf->Set_state('p');
    omega->Pf->Iprod(omega->Pf);
    omega->Pf->Set_state('t');
    dscal(hjtot, dtinv, omega->Pf->base_hj, 1); 

//    free(dudx);

    break;
  }
case SPM_Post: {
#if 1
    static int SPM_step = 0;
//    double *dpdx = omega->uf[0],
//           *dpdy = omega->vf[0], 
//           *dpdz = omega->wf[0]; 
    double    *dpdx   = grad_work[0];
    double    *dpdy   = grad_work[1];
    double    *dpdz   = grad_work[2];
    
    P->Trans(P, J_to_Q);
    P->Set_state('p');
    P->Grad_d (dpdx,dpdy,dpdz,'a'); 
    P->Set_state ('t');
    dscal(htot, dt/getgamma(),dpdx,1);
    dscal(htot, dt/getgamma(),dpdy,1);
    dscal(htot, dt/getgamma(),dpdz,1);
    dneg(htot, dpdx, 1); 
    dneg(htot, dpdy, 1);
    dneg(htot, dpdz, 1);
    
//    dcopy(htot, P->base_h, 1, omega->Pf->base_h, 1); //solid pressure save to Pf
#endif

    U->Set_state('p');
    V->Set_state('p');
    W->Set_state('p');

    U->Trans(U, P_to_F); 
    V->Trans(V, P_to_F); 
    W->Trans(W, P_to_F); 

//    dvadd(htot, omega->uf[0], 1, U->base_h, 1, U->base_h, 1 );  //U=U^*+ -dt*grad(p_p)/gama0
//    dvadd(htot, omega->vf[0], 1, V->base_h, 1, V->base_h, 1 );  //
//    dvadd(htot, omega->wf[0], 1, W->base_h, 1, W->base_h, 1 );  //
    dvadd(htot, dpdx, 1, U->base_h, 1, U->base_h, 1 );  //U=U^*+ -dt*grad(p_p)/gama0
    dvadd(htot, dpdy, 1, V->base_h, 1, V->base_h, 1 );  //
    dvadd(htot, dpdz, 1, W->base_h, 1, W->base_h, 1 );  //
    
    dvadd(htot, Uf->base_h, 1, U->base_h, 1, U->base_h, 1 );  //U= +phi*(up-u*)
    dvadd(htot, Vf->base_h, 1, V->base_h, 1, V->base_h, 1 );
    dvadd(htot, Wf->base_h, 1, W->base_h, 1, W->base_h, 1 );

    U->Trans(U,Q_to_J); U->Set_state('t');
    V->Trans(V,Q_to_J); V->Set_state('t');
    W->Trans(W,Q_to_J); W->Set_state('t');

    
    dvadd(hjtot, P->base_hj, 1, ps[0], 1, P->base_hj, 1); //total pressure
    //P->base_h has total pressure in (Q,F) space
    P->Trans(P, J_to_Q);

//    free(dpdx);

    break;
  }
#endif
  default:
    error_msg(MakeF--unknown type of action);
    break;
  }

  return;
}


/* Do initial time step assuming for case where startup field is in physical *
 * space or copy V field to Vs if it is a restart                           */

static void StartUp(Domain *omega, double *time_1, int *step_1){
  int       hjtot    = omega->U->nz*omega->U->hjtot;
  //  int       htot     = omega->U->nz*omega->U->htot;
  int       step     = *step_1;
  double    time     = *time_1;
  ACTION    WaveProp = (ACTION) iparam("EQTYPE");
  double    *zerosv = (double *) NULL,
            *varx = zerosv, 
            *vary = zerosv;
#ifdef TIMING
  double st1, st2;
  st1 = dclock();
  st2 = MPI_Wtime();
#endif

  if(omega->U->fhead->state == 'p'){
    set_order(step+1);

    MakeF (omega, WaveProp);

    Integrate_SS(Je, dt, omega->U, omega->Uf, omega->u, omega->uf);
    Integrate_SS(Je, dt, omega->V, omega->Vf, omega->v, omega->vf);
    Integrate_SS(Je, dt, omega->W, omega->Wf, omega->w, omega->wf);
#ifdef THERMO
    Integrate_SS(Je, dt, omega->T, omega->Tf, omega->Tt, omega->tf);
#endif
    Vdebug(omega,"After Integrate");

    MakeF (omega, Pressure);
    /* zero initial guess in transformed space */
    dzero (hjtot, omega->P->base_hj, 1);
    solve (omega->P,omega->Uf,omega->Pbc,omega->Pressure_sys,Helm,step,zerosv);
    Vdebug(omega,"After P_solve");

#ifdef NOISE_FILTER
    for(int i = 0; i < omega->P->nz; ++i)
      noise_filter (omega->P->flevels[i]);
    Timing("noise_filter....");
    Vdebug(omega,"After noise_filter");
#endif

    MakeF(omega, Viscous);
    /* zero initial guess in transformed space */
#if defined(MAP) && defined(DIRNEW)
     varx = omega->mapx->t;
    vary = omega->mapy->t;
#endif    

    dzero(hjtot, omega->U->base_hj, 1);
    solve(omega->U,omega->Uf,omega->Ubc,omega->Usys,Helm,step,varx);

    dzero(hjtot, omega->V->base_hj, 1);
    solve(omega->V,omega->Vf,omega->Vbc,omega->Vsys,Helm,step,vary);

    dzero(hjtot, omega->W->base_hj, 1);
#ifndef JAMES
    solve(omega->W,omega->Wf,omega->Wbc,omega->Wsys,Helm,step,zerosv);

#ifdef THERMO
    dzero(hjtot, omega->T->base_hj, 1);
    solve(omega->T,omega->Tf,omega->Tbc,omega->Tsys,Helm,step,zerosv);
#endif

#else
    solve_w(omega, step);
#endif
    Vdebug(omega,"After V_solve");

    MakeF(omega, Post);

#ifdef SPM 
    MakeF (omega, SPM_Prep); //trans to physical space for u*
    Vdebug(omega,"SPM Prep  .......");

    compute_indicator_function(omega, omega->vSPM);
    dzero(omega->U->htot*omega->U->nz, omega->CONCENTR->base_h, 1); 
    for (int i = 0; i < omega->vSPM[0].num_partls; ++i)
       dvadd(omega->U->htot*omega->U->nz,omega->vSPM[0].gaussian[i],1, omega->CONCENTR->base_h, 1, omega->CONCENTR->base_h, 1);

     omega->CONCENTR->Set_state('p');
//    for(int i=0; i<omega->U->htot*omega->U->nz; ++i)
//     if(omega->CONCENTR->base_h[i] > 0.5)
//        fprintf(stdout," error indicator = %g  i = %d \n",omega->CONCENTR->base_h[i],i);
        //calculate up^n field,stored in Omega->Uf, Vf, Wf
    compute_upf(omega);
    Vdebug(omega,"SPM upf .......");

    MakeF (omega, SPM_Pressure);
    Vdebug(omega,"SPM Pressure..");      
    SPM_SetPBCs(omega);
    solve (omega->P,omega->Pf,omega->Pbc,omega->Pressure_sys,Helm,step,zerosv);
    Vdebug(omega,"SPM Solve...");

    MakeF (omega, SPM_Post);
    Vdebug(omega,"SPM Postprocess ........");
#endif

    Analyser(omega, ++step, (time += dt));
    Vdebug(omega,"After Analyser");

    *step_1 = step;
    *time_1 = time;
  }
  else{ /* Put V into Vs */ 
#ifndef SAVINGSD
    dcopy(hjtot, omega->U->base_hj, 1, omega->us, 1); 
    dcopy(hjtot, omega->V->base_hj, 1, omega->vs, 1); 
    dcopy(hjtot, omega->W->base_hj, 1, omega->ws, 1);
#ifdef THERMO
    dcopy(hjtot, omega->T->base_hj, 1, omega->ts, 1);

#endif
#endif
///
//debug SPM
#if 0
    Integrate_SS(Je, dt, omega->U, omega->Uf, omega->u, omega->uf);
    Integrate_SS(Je, dt, omega->V, omega->Vf, omega->v, omega->vf);
    Integrate_SS(Je, dt, omega->W, omega->Wf, omega->w, omega->wf);
#ifdef THERMO
    Integrate_SS(Je, dt, omega->T, omega->Tf, omega->Tt, omega->tf);
#endif
    Vdebug(omega,"After Integrate");

    MakeF (omega, Pressure);
    /* zero initial guess in transformed space */
    dzero (hjtot, omega->P->base_hj, 1);
    solve (omega->P,omega->Uf,omega->Pbc,omega->Pressure_sys,Helm,step,zerosv);
    Vdebug(omega,"After P_solve");

#ifdef NOISE_FILTER
     for(int i = 0; i < omega->P->nz; ++i)
      noise_filter (omega->P->flevels[i]);
    Timing("noise_filter....");
    Vdebug(omega,"After noise_filter");
#endif

    MakeF(omega, Viscous);
    /* zero initial guess in transformed space */
#if defined(MAP) && defined(DIRNEW)
     varx = omega->mapx->t;
    vary = omega->mapy->t;
#endif    

    dzero(hjtot, omega->U->base_hj, 1);
    solve(omega->U,omega->Uf,omega->Ubc,omega->Usys,Helm,step,varx);

    dzero(hjtot, omega->V->base_hj, 1);
    solve(omega->V,omega->Vf,omega->Vbc,omega->Vsys,Helm,step,vary);

    dzero(hjtot, omega->W->base_hj, 1);
#ifndef JAMES
    solve(omega->W,omega->Wf,omega->Wbc,omega->Wsys,Helm,step,zerosv);

#ifdef THERMO
    dzero(hjtot, omega->T->base_hj, 1);
     solve(omega->T,omega->Tf,omega->Tbc,omega->Tsys,Helm,step,zerosv);
#endif

#else
    solve_w(omega, step);
#endif
    Vdebug(omega,"After V_solve");

    MakeF(omega, Post);
//////
//
#endif
    // ce107
    Analyser(omega, step, time);
    Vdebug(omega,"After Analyser");

#ifdef JAMES
    omega->U->Trans(omega->U, J_to_Q);
    omega->V->Trans(omega->V, J_to_Q);
     omega->W->Trans(omega->W, J_to_Q);
#endif
  }
}

void switch_hj(Element_List *EL, double **hj){
  double *d = EL->base_hj;
  EL->Mem_shift(EL->base_h, *hj);
  *hj = d;
}


void solve_w(Domain *omega, int step){
  Element_List *W    = omega->W;
  Element_List *Wf   = omega->Wf;
  Bndry       **Wbc  = omega->Wbc;
  Bsystem    **Wbsys = omega->Wsys;
  int              i = 0, k = 0;
  double         bet;

  omega->U->Trans(omega->U, J_to_Q);
  omega->V->Trans(omega->V, J_to_Q);
  omega->U->Set_state('p');
  omega->V->Set_state('p');

  ROOT{
    if(step && step < Je){
      free(Wbc[0]->DirRHS); Wbc[0]->DirRHS = (double*) 0;
      DirBCs(W->flevels[0], Wbc[0], Wbsys[0], Helm);
    }
    
    SetBCs (W->flevels[0],Wbc[0],Wbsys[0]);
    Solve  (W->flevels[0],Wf->flevels[0],Wbc[0],Wbsys[0],Helm);
    W->flevels[0]->Trans(W->flevels[0], J_to_Q);
    dzero  (W->htot, W->base_h+W->htot, 1);
    W->Set_state('p');
    if(W->nz == 2) return;
    i=2;
  }
  
  double *tmp = dvector(0, W->htot-1);
  
  for(k=i;k<W->nz;k+=2){
    bet = 1.0/Beta(k);
    omega->U->flevels[k]->Grad_d(W->flevels[k+1]->base_h, NULL, NULL, 'x');
    dscal (W->htot, bet,  W->flevels[k+1]->base_h, 1);
    omega->V->flevels[k]->Grad_d(NULL, tmp, NULL, 'y');
    daxpy (W->htot, bet, tmp, 1, W->flevels[k+1]->base_h, 1);
    W->flevels[k+1]->Set_state('p');
    
    omega->U->flevels[k+1]->Grad_d(W->flevels[k]->base_h, NULL, NULL, 'x');
    dscal (W->htot, -bet,  W->flevels[k]->base_h, 1);
    omega->V->flevels[k+1]->Grad_d(NULL, tmp, NULL, 'y');
    daxpy (W->htot, -bet, tmp, 1, W->flevels[k]->base_h, 1);
    W->flevels[k]->Set_state('p');
  }

  free(tmp);
}    

#ifdef TRY_ACC
static void Accel_term_implicit(Domain *omega, double dt)
{
  Element_List *U = omega->U,
               *V = omega->V;
  Element      *fU;
  double f;
  register int i, j, k, plane;
  double nu = dparam("KINVIS");
  int NZ    = U->nz;
  const int nel = omega->U->nel;

  // Adding "acceleration" terms
  for (k = 0, plane = NZ*option("PROCID"); k < NZ; k++, plane++) {
    if ((f = (omega->mapx->tzz[plane]*nu-omega->mapx->tt[plane])*dt) != 0.0) {
      for (i = 0; i < nel; i++) {
	fU = U->flevels[k]->flist[i]; 
	for (j = 0; j < fU->Nverts; j++)
	  fU->vert->hj[j] += f;
      }
    }
    if ((f = (omega->mapy->tzz[plane]*nu-omega->mapy->tt[plane])*dt) != 0.0) {
      for (i = 0; i < nel; i++) {
	fU = V->flevels[k]->flist[i];
	for (j = 0; j < fU->Nverts; j++)
	  fU->vert->hj[j] += f;
      }
    }
  }
  return;
}
#endif

#ifdef XMLRPC
void XMLRPC_update_boundary(int step)
 {
  int hisstep   =iparam("HISSTEP");
  char arr[XMLRPC_BUFFSIZE];
  char *eptr;
  char *token;
//  xmlrpc_value * action;
  const char *action;
  char *action_new;
  xmlrpc_value * resultP; 
  double *action_values   = dvector(0,2-1);
  action_values[0] = 0.0;
  action_values[1] = 0.0;
  int train_freq = iparam("NUM_TRAIN_FREQ");
  xmlrpc_value *resultP_train; 
  xmlrpc_value *resultP_restart; 

//  double current_u1 = dparam("D_U1");
//  double current_u2 = dparam("D_U2");

  double Cd = dparam("D_DRAGCOEFF");
  double Cf = dparam("D_LIFTCOEFF");
  if (step % hisstep == 0){
   ROOTONLY
   {
     sprintf(arr, "%lf_%lf", Cf,Cd);
     resultP = xmlrpc_client_call(&env, serverUrl, methodName,
                                 "(s)", arr);

     dieIfFaultOccurred(&env);
    /* Get our action and print it out. */
     xmlrpc_read_string(&env, resultP, &action);
     dieIfFaultOccurred(&env);
    
//     size_t action_size =  sizeof(action);

//     memcpy(action_new, action, action_size);
     action_new = strdup(action);
     token = strtok(action_new, "_");
     int ii=0;
     while(token != NULL) {
      action_values[ii] = strtod(token,&eptr);
//      fprintf(stderr,"Returned action value %d: %g \n",ii, action_value);
      token = strtok(NULL,"_");
      ii++;
      }

     } //end of root
    MPI_Bcast(action_values,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier  (MPI_COMM_WORLD);                  /* sync before work */

    current_u1 = next_u1;
    current_u2 = next_u2;

    next_u1 = action_values[0];
    next_u2 = action_values[1];

    action_counter ++;
   
    ROOTONLY
    if(action_counter%train_freq == 0)
    {
     resultP_train = xmlrpc_client_call(&env, serverUrl, "train",
                                 "(i)", train_freq*10);

     resultP_train = xmlrpc_client_call(&env, serverUrl, "save",
                                 "(i)", train_freq);

     resultP_restart = xmlrpc_client_call(&env, serverUrl, "start_episode",
                                 "(i)", -1);
    }
   }

  double intermediate_u1 = current_u1 +double(step% hisstep)/hisstep*(next_u1-current_u1);
  double intermediate_u2 = current_u2 +double(step% hisstep)/hisstep*(next_u2-current_u2);

  double scaling = dparam("D_SCALING_XMLRPC");
  dparam_set("D_U1", intermediate_u1*scaling);
  dparam_set("D_U2", intermediate_u2*scaling);

  ROOTONLY fprintf(stderr,"rotating speed: %lf  %lf \n",intermediate_u1*scaling,intermediate_u2*scaling);

  free(action_values);

 }
#endif

