/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source:
 * $Revision:
 * $Date:
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include "nektarF.h"
#ifndef OLDFFTS
#include <rfftw.h>
extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;
#endif
#ifdef MAP
#include "map.h"
#endif

void compute_local_vorticity_flux(Domain *omega, int step, double time);
#ifdef THERMO
void compute_local_nusselt(Domain *omega, double time);
void output_temperature(Domain *omega, double time);
#endif

#ifdef SPM
void compute_hydrodynamic_force(Domain *Omega, double **indicator, 
                                double ***x_o, double ***y_o, double ***FHx, double ***FHy, double ***NH, double ***FHphix, double ***FHphiy);
#endif

double       **u2nd = (double **) NULL;
double       **u2nde= (double **) NULL;
static double **avg = (double **) NULL;

extern int SPM_MAP_step ;
/*
 * Run-time Analyzer ... called every time step
 */
static int    History       (Domain *omega, double time);
static void   hisHeader     (Domain *omega);
static void   intHeader     (Domain *omega);
static void   SaveInt       (Domain *omega, double time);
//static void   addfields     (Element_List *V[]);
       void   addfields     (Element_List *V[]);
static void   MEnergy       (Domain *omega, FILE *fpenergy, double time);
static int check_number=0;

static int init = 1, verbose, iostep, hisstep, mstep, nsteps, timeavg;
static double last_utime, last_wtime, clock_skip;
static FILE *fpenergy;

double zmesh(int);

void Analyser (Domain *omega, int step, double time){
  FILE     *fp = 0;
  int      nfields = DIM+2;
  int      i,j;
  Element_List  **V;
  char      fname[FILENAME_MAX];
  double    step_length;

  dparam_set("t", time);
  
  int nf = nfields;
//#ifdef MAP
//  nf++;
//#endif

#ifdef THERMO
  nf++;
  nfields++;
#endif

#ifdef SPM
  nf++;
  nfields++;
#endif

#ifdef OUTPUT_VISCOSITY
   nfields=nfields+2;
   nf +=2; 
#endif

  V = (Element_List**) malloc((nf)*sizeof(Element*));

  if(init){
    verbose   = option("verbose");
    iostep    = iparam("IOSTEP");
    hisstep   = option("hisstep");
    mstep     = iparam("MSTEP");
    if (mstep == 0) mstep = hisstep;
    nsteps    = iparam("NSTEPS");
    timeavg   = option("timeavg");
    init = 0;
    ROOT {
      if (omega->his_list) hisHeader (omega);
      if (omega->int_list) intHeader (omega);
      sprintf (fname, "%s.mea", omega->name);
      fpenergy = fopen(fname, "w");
      fprintf (fpenergy,"#Time Ek[k=0,...,k=%d], Total\n",option("NZTOT")/2-1);
      fflush (fpenergy);
    }
    last_utime = dclock();
    clock_skip = pow(2.0,sizeof(clock_t)*8.0)*1E-6;
#ifdef PARALLEL
    last_wtime = MPI_Wtime();
#endif
  } else {
    /* .......... General Output ......... */
    step_length = dclock()-last_utime;
    /* a negative increment implies we are using clock() which rolls over
       to negative values when an the capacity of the variable it returns 
       is reached. Here we assume that the rollover doesn't happen twice
       in the same timestep */
    if (step_length < 0.0) step_length += clock_skip;
#ifdef PARALLEL
    ROOT
      printf("Time step = %d, Time = %g User-Time = %g Wallclock-Time = %g\n",
	     step, time, step_length, MPI_Wtime()-last_wtime);
      last_wtime = MPI_Wtime();
#else
#ifndef FLOK
      printf("Time step = %d, Time = %g User-Time = %g\n",
	     step, time, step_length);
#endif
#endif
    
    last_utime = dclock();
  }

  if (step == 0){                        /* Everything else is for step > 0 */
    if (omega->his_list)                 /* Do initial history point        */
      History (omega, time);
#ifndef FLOK
    MEnergy (omega, fpenergy, time);
#endif

#ifndef SPM
    if (omega->fce_file){
      forces(omega,time);
      ROOT fflush(omega->fce_file);
    }
#else
   if(strcmp(omega->vSPM[0].shape[0].type,"wall") == 0) //for 'wall'
    if (omega->fce_file){
      forces(omega,time);
      ROOT fflush(omega->fce_file);
    }
#endif
    return;
  }

#ifndef SPM
  if (omega->fce_file){
    forces(omega,time);
    ROOT fflush(omega->fce_file);
  }  
#else
   if(strcmp(omega->vSPM[0].shape[0].type,"wall") == 0) //for 'wall'
    if (omega->fce_file){
      forces(omega,time);
      ROOT fflush(omega->fce_file);
    }
#endif
 
  if (step % hisstep == 0){
    if (omega->his_list){
      History (omega, time);
      ROOT fflush(omega->his_file);
    }
    if (omega->int_list){
      interp (omega);
      ROOT {
	SaveInt (omega, time);
	fflush(omega->int_file);
      }
    }

 if(iparam("IOUTPUT_WALL_VORTICITY") && step > 1)
    compute_local_vorticity_flux(omega,step,time);

#if !defined(MAP) && !defined(SPM)
//#ifndef MAP
    if (omega->fce_file){
      forces(omega,time);
      ROOT fflush(omega->fce_file);
    }   
#endif
 
#if !defined(MAP) && defined(SPM)
  ROOT omega->vSPM[0].hist_data_out(time);
#endif

//#ifdef THERMO
//  if(dparam("DNUSSELT"))
//    compute_local_nusselt(omega,time);
//
//  if(dparam("DTEMPOUTPUT"))
//  output_temperature(omega,time);
//#endif
    /* flush stdout at step as well */
    fflush(stdout);
  }
  
#ifndef FLOK
  if (step % mstep == 0) MEnergy (omega, fpenergy, time);
#endif
#if defined(MAP) && !defined(TEST)
  /* .......... Cylinder variables ..... */
//  if( SPM_MAP_step > iparam("ISPM_MAP_STEP") )
  run_obj_vars(omega);
#endif

//#ifdef SPM
//  omega->vSPM[0].writedog2();
//#endif

  /* ..........  Field Files   ......... */
  V[0]   = omega->U;  V[1]   = omega->V;  V[2]   = omega->W;
  V[3]   = omega->P;

#if defined(THERMO) && defined(SPM)
  V[4]   = omega->T;
  omega->CONCENTR->Trans(omega->CONCENTR,P_to_F);
  omega->CONCENTR->Trans(omega->CONCENTR,Q_to_J);
  omega->CONCENTR->Set_state('t');

  V[5]   = omega->CONCENTR;

#ifdef OUTPUT_VISCOSITY
  V[6]   = omega->visc1;
  V[7]   = omega->visc2;
#endif

#endif

#if !defined(THERMO) && defined(SPM)

  omega->CONCENTR->Trans(omega->CONCENTR,P_to_F);
  omega->CONCENTR->Trans(omega->CONCENTR,Q_to_J);
  omega->CONCENTR->Set_state('t');

  V[4]   = omega->CONCENTR;

#ifdef OUTPUT_VISCOSITY
  V[5]   = omega->visc1;
  V[6]   = omega->visc2;
#endif

#endif

#if defined(THERMO) && !defined(SPM)
  V[4]   = omega->T;
#endif

#if !defined(THERMO) && !defined(SPM) 

#ifdef OUTPUT_VISCOSITY
  V[4]   = omega->visc1;
  V[5]   = omega->visc2;
#endif

#endif

 // if( (timeavg)&&(!dparam("PIPE")) ) addfields (V);

  if (step % iostep == 0 && step < nsteps) {         

#if 0
    if(verbose)
#ifndef MAP
    { double *tmp = dvector(0,V[0]->htot*V[0]->nz-1);
      for(j=0;j<nfields-1;++j){
//      for(j=0;j<nfields;++j){

	V[j]->Trans(V[j], J_to_Q);
 
	dcopy(V[j]->htot*V[j]->nz, V[j]->base_h, 1, tmp, 1);
	if(V[j]->fhead->type == 'p')
	  V[j]->Trans(V[j], J_to_Q);
	V[j]->Trans(V[j], F_to_P);

	for(i = 0; i < V[0]->nz; ++i){
	  dparam_set("z", zmesh(i));
	  V[j]->flevels[i]->Terror(omega->soln[j]); 
	}
	dcopy(V[j]->htot*V[j]->nz, tmp, 1, V[j]->base_h, 1);

	V[j]->Set_state('t');
      }
      free(tmp);
    }
#else
    Compare (omega, Fourier);
#endif

#endif



    if (option ("checkpt")) {
      if(option("SLICES")&&check_number<1000){
	sprintf (fname, "%s_%d.chk",  omega->name,check_number);
	++check_number;
      }
      else {
	sprintf (fname, "%s.chk", omega->name);
	ROOT {
	  int err = backup(fname);
	  if (err != 0 && step != 0) 
	    fprintf(stderr,"WARNING: Analyser: failed to backup the chk file\n");
	}
      }
	
      ROOT fp = fopen(fname,"w");

      for(j=0;j<nfields;++j)
	V[j]->Set_state('t');

#ifdef MAP
      ROOT WriteMap (omega);
//      nfields++; // to output Wlast as well
//
//      Nek_Trans_Type f_to_p = F_to_P,
//	             p_to_f = P_to_F;
//      if (option("dealias")) {
//	f_to_p = F_to_P32;
//	p_to_f = P_to_F32;
//      }
//
//      V[nfields-1]->Trans(V[nfields-1], p_to_f); // Get us back to Fourier Space
//      V[nfields-1]->Trans(V[nfields-1], Q_to_J); // Get us back to modal space
//

#endif 

#ifdef SPM
     ROOT omega->vSPM[0].WriteInit(omega->name);
#endif

      if (option("parts"))
	pWriteFieldF (fp, omega->name, fname, step, time, nfields, V);
      else
	WritefieldF (fp, omega->name, step, time, nfields, V);
      ROOT fclose (fp);
#if 0
#ifdef MAP // Return Wlast to quadrature/physical space
      V[nfields]->Set_state('p'); // no need for J_to_Q
      V[nfields]->Trans(V[nfields], f_to_p);
#endif
#endif      
      // save the modal storage of P
      dcopy (V[0]->hjtot*V[0]->nz, omega->P->base_hj, 1, omega->U->base_h, 1);

      if (timeavg) averagefields(omega, step, time); //U,V,W,P...
      if(option("STATAVG"))
	{
	  average_u2_avg(omega, step, time); //u*u v*v....
	  if (option("VARV"))
	    average_u2e_avg(omega, step, time);
	}
      
      // restore the modal storage of P
      dcopy (V[0]->hjtot*V[0]->nz, omega->U->base_h, 1, omega->P->base_hj, 1);
    }
  }
  
#ifdef SPM
  omega->CONCENTR->Trans(omega->CONCENTR,J_to_Q);
  omega->CONCENTR->Set_state('p');
  omega->CONCENTR->Trans(omega->CONCENTR,F_to_P);
#endif
  
  free (V);
  return;		
}
    

/* ------------------------------------------------------------------------- *
 * History() -- Process history points                                       *
 *                                                                           *
 * This function processes the history point list and transforms data points *
 * if necessary to physical space.                                           *
 * ------------------------------------------------------------------------- */
  static int gatherPts     (HisPoint *hp, Element_List **V, double **vbuf);

  static int History (Domain *omega, double time)
{
  FILE     *fp = omega->his_file;
  HisPoint *hp = omega->his_list;
   Element_List  **V;
  double   *vbuf[HP_MAX];
  register  int n, cnt;

  V = (Element_List**) malloc((DIM+2)*sizeof(Element_List*));
  if (!hp) return 0;

  V[0] = omega->U; V[1] = omega->V; V[2] = omega->W; V[3] = omega->P;

  gatherPts (hp, V, vbuf);
  cnt     = 0;
  
  ROOT do 
    { fprintf (fp, "%#13.6g ", time);
      for (n = 0; n < strlen(hp->flags); n++)
	fprintf (fp, "%#13.11g ", vbuf[cnt][n]);
      fprintf (fp, ":%d\n", cnt+1);
      if (iparam("ISO") !=0) fprintf (fp, "%#13.6g ", time);     /* VARV */
      free (vbuf[cnt++]);    } 
  while (hp = hp->next);
  
  free(V);
  return cnt;
}

/* Collect history points from the Elements */
static int vFe[][2] = {{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
static double *modecenter(Element *E, Edge *e);

static int gatherPts (HisPoint *hp, Element_List **V, double **vbuf)
{
  register  int i, j, n, pos;
  int       nz = option("NZ"), nztot = option("NZTOT");
  double    *tmp, *tmps;

  tmps = dvector(0,nz-1);
  ROOT tmp = dvector(0,nztot-1);

  for (i = 0; hp ; ++i, hp = hp->next) {
    vbuf[i] = dvector (0, strlen (hp->flags)-1);
    // hard wired for (DIM+1) V and P
    for  (n = pos = 0; n < DIM+2; ++n) {
      if (strchr (hp->flags, V[n]->fhead->type)) {
	switch (hp->mode) {
	case TransVert:
	  if(V[n]->fhead->state == 't'){
	    for(j = 0; j < nz; ++j)
	      tmps[j] = V[n]->flevels[j]->flist[hp->id]->vert[hp->i].hj[0];
	     
#ifdef PARALLEL
	    MPI_Gather (tmps, nz, MPI_DOUBLE, tmp, nz, MPI_DOUBLE, 0, 
			MPI_COMM_WORLD); 
#else
	    dcopy(nz,tmps,1,tmp,1);
#endif


	    ROOT{
#ifdef OLDFFTS
	      realft(nztot/2,tmp,1);
#else
	      rfftw(rplan_inv, 1, (FFTW_COMPLEX *) tmp, 1, 0, 0, 0, 0);
#endif
	      vbuf[i][pos++] = tmp[hp->k];
	    }
	  }
	  break;
	case TransEdge:
	  if(V[n]->fhead->state == 't'){
	    int    k;
	    Edge   *e = V[n]->flist[hp->id]->edge+hp->i;
	    double *mode = modecenter(V[n]->flist[hp->id],e);
	    
	    for(k = 0; k < nz; ++k){
	      tmps[k] =
	0.5*(V[n]->flevels[k]->flist[hp->id]->vert[vFe[hp->i][0]].hj[0] + 
             V[n]->flevels[k]->flist[hp->id]->vert[vFe[hp->i][1]].hj[0]);

	      e = V[n]->flevels[k]->flist[hp->id]->edge+hp->i;
	      for(j = 0; j < e->l; ++j) tmps[k] += e->hj[j]*mode[j];
	    }
#ifdef PARALLEL
	    MPI_Gather (tmps, nz, MPI_DOUBLE, tmp, nz, MPI_DOUBLE, 0, 
			MPI_COMM_WORLD); 
#else
	    dcopy(nz,tmps,1,tmp,1);
#endif

	    ROOT{
#ifdef OLDFFTS
	      realft(nztot/2,tmp,1);
#else
	      rfftw(rplan_inv, 1, (FFTW_COMPLEX *) tmp, 1, 0, 0, 0, 0);
#endif
	      vbuf[i][pos++] = tmp[hp->k];
	    }
	   }
	  break;
	case Physical:
	  break;
	default:
	  error_msg (History -- unknown history point mode);
	  break;
	}
      }
    }
  }
  free(tmps);
  ROOT
    free(tmp);

  return i;
}

/* find the center of modes  and store */
static double *modecenter(Element *E, Edge *edg){
  static double *mode;
  static int Lmode;

  if(!(mode&&(edg->l <= Lmode))){
    int i,qa = E->qa;
    Mode *e = E->getbasis()->edge[0];
    
    if(!mode) free(mode);
    mode = dvector(0,edg->l);

    if(qa%2 == 0){ /* if even spacing interpolate to center point */
      double **im;
      getim(qa,qa+1,&im,a2a);
      
      for(i = 0; i < edg->l; ++i)
	mode[i] = ddot(qa,im[qa/2],1,e[i].a,1);
    }      
    else           /* else use center value which is always at center point */
      for(i = 0; i < edg->l; ++i)
	mode[i] = e[i].a[qa/2];
  }

  return mode;
}

/* Write the header for the history point file */

static void hisHeader (Domain *omega)
 {
  FILE      *fp = omega->his_file;
  HisPoint  *hp = omega->his_list;
  Element_List   *U  = omega->U;
  int        n  = 1;

  if (!fp) return;

  fputs ("# Nektar history point file\n"
	 "# \n"
	 "# History points:\n", fp);

  do 
    {
#if DIM == 3
      if(hp->mode == TransVert){
	fprintf (fp, "#  %3d: x = %#6.2lf, y = %#6.2lf, z = %#6.2lf, ", n++, 
		 U->flist[hp->id]->vert[hp->i].x,
		 U->flist[hp->id]->vert[hp->i].y,
		 U->[hp->id]->vert[hp->i].z);
      }
      else if(hp->mode == TransEdge){ /* assumes edge is straight */
	fprintf (fp, "#  %3d: x = %#6.2lf, y = %#6.2lf, z = %#6.2lf, ", n++,
		 0.5*(U->flist[hp->id]->vert[vFe[hp->i][0]].x + 
		      U->flist[hp->id]->vert[vFe[hp->i][1]].x),
		 0.5*(U->flist[hp->id]->vert[vFe[hp->i][0]].y + 
		      U->flist[hp->id]->vert[vFe[hp->i][1]].y),
		 0.5*(U->flist[hp->id]->vert[vFe[hp->i][0]].z + 
		      U->flist[hp->id]->vert[vFe[hp->i][1]].z));
      }	
#else
      if(hp->mode == TransVert){
	fprintf (fp, "#  %3d: x = %#6.2lf, y = %#6.2lf, z = %#6.2lf, ", n++, 
		 U->flist[hp->id]->vert[hp->i].x,
		 U->flist[hp->id]->vert[hp->i].y,
		 hp->k*dparam("LZ")/option("NZTOT"));
      }
      else if(hp->mode == TransEdge){
	fprintf (fp, "#  %3d: x = %#6.2lf, y = %#6.2lf, z = %#6.2lf, ", n++,
		 0.5*(U->flist[hp->id]->vert[vFe[hp->i][0]].x + 
		      U->flist[hp->id]->vert[vFe[hp->i][1]].x),
		 0.5*(U->flist[hp->id]->vert[vFe[hp->i][0]].y + 
		      U->flist[hp->id]->vert[vFe[hp->i][1]].y),
		   hp->k*dparam("LZ")/option("NZTOT"));
      }	
#endif
      fprintf (fp, "fields = %4s, [%2d %3d %4d]\n",
	       hp->flags,  hp->id+1, hp->i+1, hp->k+1);  
     }
  while
    (hp = hp->next);
  fputs ("#\n", fp);
  return;
}

static void intHeader (Domain *omega)
{
  FILE           *fp = omega->int_file;
  intepts        *I  = omega->int_list;
  int             i;
 
  if (!fp) return;
 
  fputs ("# Nektar interpolation point file\n"
         "# \n"
         "# interpolation points:\n", fp);
 
  for (i = 0; i < I->npts; i++)
    fprintf (fp, "#  %3d: x = %#6.2lf, y = %#6.2lf, z = %#6.2lf, fields = uvwp Fu, Fv, Fw, Fp\n", i+1, 
             I->X.x[i], I->X.y[i], I->X.z[i]);
  fputs ("#\n", fp);
 
  return;
}
 
/* Save interpolated points' date */
 
static void SaveInt (Domain *omega, double time)
{
  FILE    *fp = omega->int_file;
  Intepts *I = omega->int_list;
  int      i, j, n, nfields, nztot = option("NZTOT");
  
  nfields = DIM + 2;
 
  for (i = 0; i < I->npts; i++){ 
    fprintf (fp, "%lf ", time);
    for (n = 0; n < nfields; n++)
      fprintf (fp, "%#13.6g ", I->ui[i][n]);
    for (n = 0; n < nfields; n++)
      for (j = 0; j < nztot; j++)
	fprintf (fp, "%#13.6g ", I->ufr[i][n*nztot+j]);
    fprintf (fp, ":%d\n", i+1); 
  }
  
  return;
 }

//static void addfields (Element_List *V[]){
void addfields (Element_List *V[]){
  register int i;
  int nf    = DIM+2;
#ifdef THERMO
  nf++;
#endif
  int njtot = V[0]->nz*V[0]->hjtot;

  int          qa,qb,m;
  Element      *E,*F,*G,*H;

  qa = V[0]->flevels[0]->fhead->qa,
  qb = V[0]->flevels[0]->fhead->qb;

  //convert to cylindric coordinate first!
  int nz = option("NZ"), nztot = option("NZTOT");
  if (option("dealias")) nz=3*nz/2;

//zwang
  if(dparam("PIPE"))
   {
  	Coord X;
	  X.x  = dvector(0, qa*qb-1);
	  X.y  = dvector(0, qa*qb-1);
    //convert to cylindric coordinate first!
  //  int nz = option("NZ"), nztot = option("NZTOT");
  //  if (option("dealias")) nz=3*nz/2;

    for (int k = 0; k < nz; k++)
    {
     for(int eid =0; eid<V[0]->nel; ++eid)
      {

	     E   = V[0]->flevels[0]->flist[eid]; //u
//	     F   = V[1]->flevels[0]->flist[eid]; //v
//	     G   = V[2]->flevels[0]->flist[eid]; //w
//	     H   = V[3]->flevels[0]->flist[eid]; //p

       
       qa = E->qa;   qb = E->qb;
	     E->coord(&X);
       
       int offset = V[0]->htot*k+eid*qa*qb;

       for ( int i=0; i<qa; i++ )
        for ( int j=0; j<qb; j++ )
         {
            m = i*qb + j;
            double phi = atan2(X.x[m],X.y[m]);
            
        //    double ur   = cos(phi)*E->h[0][m]+sin(phi)*F->h[0][m];
        //    double uphi = -sin(phi)*E->h[0][m]+cos(phi)*F->h[0][m];
        //    double uz   = G->h[0][m];
        //    double up   = H->h[0][m];
            double ur   = cos(phi)*V[0]->base_h[offset+m]+sin(phi)*V[1]->base_h[offset+m];
            double uphi = -sin(phi)*V[0]->base_h[offset+m]+cos(phi)*V[1]->base_h[offset+m];
            double uz   = V[2]->base_h[offset+m];
            double up   = V[3]->base_h[offset+m];
        #ifdef THERMO
            double ut   = V[4]->base_h[offset+m];
        #endif
            
        //    avg[0][k*V[0]->htot+eid*E->qtot+m] += ur;
        //    avg[1][k*V[0]->htot+eid*E->qtot+m] += uphi;
        //    avg[2][k*V[0]->htot+eid*E->qtot+m] += uz;
        //    avg[3][k*V[0]->htot+eid*E->qtot+m] += up;
        //    u2nd[0][k*V[0]->htot+eid*E->qtot+m] += ur*ur;
        //    u2nd[1][k*V[0]->htot+eid*E->qtot+m] += uphi*uphi;
        //    u2nd[2][k*V[0]->htot+eid*E->qtot+m] += uz*uz;
        //    u2nd[3][k*V[0]->htot+eid*E->qtot+m] += uz*ur;
        //    u2nd[4][k*V[0]->htot+eid*E->qtot+m] += up*up;

            avg[0][offset+m] += uz;
            avg[1][offset+m] += ur;
            avg[2][offset+m] += uphi;
            avg[3][offset+m] += up;
         #ifdef THERMO
            avg[4][offset+m] += ut;
         #endif

            u2nd[0][offset+m] += uz*uz;
            u2nd[1][offset+m] += ur*ur;
            u2nd[2][offset+m] += uphi*uphi;
            u2nd[3][offset+m] += uz*ur;
            u2nd[4][offset+m] += up*up;
         #ifdef THERMO
            u2nd[5][offset+m] += ut*ut;
            u2nd[6][offset+m] += ut*uz;
         #endif

         }


      }
 
    }
     free(X.x);
     free(X.y);
  }
  else if(dparam("CYLINDER"))
  {

    for (int k = 0; k < nz; k++)
    {
     for(int eid =0; eid<V[0]->nel; ++eid)
      {

	     E   = V[0]->flevels[0]->flist[eid]; //u
       
       qa = E->qa;   qb = E->qb;
       
       int offset = V[0]->htot*k+eid*qa*qb;

       for ( int i=0; i<qa; i++ )
        for ( int j=0; j<qb; j++ )
         {
            m = i*qb + j;
            
            avg[0][offset+m] += V[0]->base_h[offset+m]; //u
            avg[1][offset+m] += V[1]->base_h[offset+m]; //v
            avg[2][offset+m] += V[2]->base_h[offset+m]; //w
            avg[3][offset+m] += V[3]->base_h[offset+m]; //p
         #ifdef THERMO
            avg[4][offset+m] += V[4]->base_h[offset+m]; //T
         #endif

            u2nd[0][offset+m] += V[0]->base_h[offset+m]*V[0]->base_h[offset+m];
            u2nd[1][offset+m] += V[1]->base_h[offset+m]*V[1]->base_h[offset+m];
            u2nd[2][offset+m] += V[2]->base_h[offset+m]*V[2]->base_h[offset+m];
            u2nd[3][offset+m] += V[0]->base_h[offset+m]*V[1]->base_h[offset+m];
            u2nd[4][offset+m] += V[3]->base_h[offset+m]*V[3]->base_h[offset+m];
         #ifdef THERMO
            u2nd[5][offset+m] += V[4]->base_h[offset+m]*V[4]->base_h[offset+m];
            u2nd[6][offset+m] += V[0]->base_h[offset+m]*V[4]->base_h[offset+m];
         #endif

         }

      }
 
    }

  }
 else
  {
   for(i = 0; i < nf; ++i)
    dvadd(njtot, V[i]->base_hj, 1, avg[i], 1, avg[i], 1);
  }

  return;
}

#if defined(NEW_OUTPUT)
void averagefields (Domain *omega, int nsteps, double time){
  register int i;
  double fac;
  int nf = DIM+2;
  Element_List *U = omega->U;
  
  int njtot = U->nz*U->hjtot;
  fac = 1.0/(double)nsteps;
  
  for(i = 0; i < nf; ++i)
    dscal(njtot, fac, avg[i], 1);

  FILE *fp;
  char fname[BUFSIZ];
  
  sprintf(fname, "%s.1st", omega->name);
  ROOT {
    int err = backup(fname);
    if (err) 
      fprintf(stderr,"WARNING: Analyser: failed to backup the chk file %s\n",
	      fname);
    fp = fopen(fname, "w");
  }
  if (option("parts"))
    pMWriteFieldF (fp, omega->name, fname, nsteps, time, nf, U, avg, "uvpwT");
  else 
    MWritefieldF(fp, omega->name, nsteps, time, nf, U, avg, "uvpwT");
  ROOT fclose (fp);

  for(i = 0; i < nf; ++i)
    dscal(njtot, ((double) nsteps), avg[i], 1);
  
  return;
}
#else
void averagefields (Domain *omega, int nsteps, double time){
  register int i,j;
  double       *ut0  =  omega->ut[0],    *vt0   =  omega->vt[0],    *wt0 = omega->wt[0], // used as temps arrary for gradients
               *ut1  =  omega->ut[1],    *vt1   =  omega->vt[1],    *wt1 = omega->wt[1]; // used for calculating laplacians
  double fac;
  int nf = DIM+2;
 #ifdef THERMO
   nf++;
 #endif
  Element_List **V;
   
  V = (Element_List**) calloc(nf, sizeof(Element_List*));
  
  V[0] = omega->Uf;  V[1] = omega->Vf; 
  V[2] = omega->Wf;  V[3] = omega->P;
 #ifdef THERMO
  V[4] = omega->Tf;
 #endif

  int njtot = V[0]->nz*V[0]->hjtot;
  int ntot = V[0]->nz*V[0]->htot;

  if (option("dealias")) ntot=3*ntot/2;
 
#ifdef THERMO
  dcopy(ntot, omega->Tf->base_h, 1, ut0, 1);
#endif

  fac = 1.0/(double)nsteps;
 
   if(dparam("PIPE") || dparam("CYLINDER"))
   {
    for(i = 0; i < nf; ++i)
      dsmul(ntot, fac, avg[i], 1, V[i]->base_h, 1);//<u>,<v>,<w>,<p>


 #ifndef THERMO
    for(i = 0; i < nf+1; ++i)
      dsmul(ntot, fac, u2nd[i], 1, u2nde[i], 1); //<uu>,<vv>,<ww>,<pp>,<uw>

    for(j = 0; j < ntot; ++j)
     {
      u2nde[0][j] =u2nde[0][j]-V[0]->base_h[j]*V[0]->base_h[j]; // <u*u> - <u>*<u> 
      u2nde[1][j] =u2nde[1][j]-V[1]->base_h[j]*V[1]->base_h[j]; // <v*v> - <v>*<v>
      u2nde[2][j] =u2nde[2][j]-V[2]->base_h[j]*V[2]->base_h[j]; // <w*w> - <w>*<w>
      u2nde[3][j] =u2nde[3][j]-V[3]->base_h[j]*V[3]->base_h[j]; // <p*p> - <p>*<p>
      u2nde[4][j] =u2nde[4][j]-V[2]->base_h[j]*V[0]->base_h[j]; // <u*w> - <u>*<w>
     }
 #else
    for(i = 0; i < nf+2; ++i)
      dsmul(ntot, fac, u2nd[i], 1, u2nde[i], 1); //<uu>,<vv>,<ww>,<pp>,<uw> <tt> <ut>

    for(j = 0; j < ntot; ++j)
     {
      u2nde[0][j] =u2nde[0][j]-V[0]->base_h[j]*V[0]->base_h[j]; // <u*u> - <u>*<u> 
      u2nde[1][j] =u2nde[1][j]-V[1]->base_h[j]*V[1]->base_h[j]; // <v*v> - <v>*<v>
      u2nde[2][j] =u2nde[2][j]-V[2]->base_h[j]*V[2]->base_h[j]; // <w*w> - <w>*<w>
      u2nde[3][j] =u2nde[3][j]-V[3]->base_h[j]*V[3]->base_h[j]; // <p*p> - <p>*<p>
      u2nde[4][j] =u2nde[4][j]-V[2]->base_h[j]*V[0]->base_h[j]; // <u*w> - <u>*<w>
      u2nde[5][j] =u2nde[5][j]-V[4]->base_h[j]*V[4]->base_h[j]; // <t*t> - <p>*<p>
      u2nde[6][j] =u2nde[6][j]-V[0]->base_h[j]*V[4]->base_h[j]; // <u*t> - <u>*<w>
     }

 #endif
     
    for(i = 0; i < nf; ++i)
     {
      V[i]->Set_state('p');
      V[i]->Trans(V[i], P_to_F);
      V[i]->Trans(V[i], Q_to_J);
     }

   }
   else if(dparam("CYLINDER"))
   {
#ifndef THERMO
    for(i = 0; i < nf; ++i)
#else
    for(i = 0; i < nf+1; ++i)
#endif
      dsmul(ntot, fac, avg[i], 1, V[i]->base_h, 1);//<u>,<v>,<w>,<p>


    for(i = 0; i < nf+1; ++i)
      dsmul(ntot, fac, u2nd[i], 1, u2nde[i], 1); //<uu>,<vv>,<ww>,<pp>,<uw>

    for(j = 0; j < ntot; ++j)
     {
      u2nde[0][j] =u2nde[0][j]-V[0]->base_h[j]*V[0]->base_h[j]; // <u*u> - <u>*<u> 
      u2nde[1][j] =u2nde[1][j]-V[1]->base_h[j]*V[1]->base_h[j]; // <v*v> - <v>*<v>
      u2nde[2][j] =u2nde[2][j]-V[2]->base_h[j]*V[2]->base_h[j]; // <w*w> - <w>*<w>
      u2nde[3][j] =u2nde[3][j]-V[3]->base_h[j]*V[3]->base_h[j]; // <p*p> - <p>*<p>
      u2nde[4][j] =u2nde[4][j]-V[2]->base_h[j]*V[0]->base_h[j]; // <u*w> - <u>*<w>
     }
     
    for(i = 0; i < nf; ++i)
     {
      V[i]->Set_state('p');
      V[i]->Trans(V[i], P_to_F);
      V[i]->Trans(V[i], Q_to_J);
     }

   }
   else //(!dparam("PIPE"))
   {
    for(i = 0; i < nf; ++i)
      dsmul(njtot, fac, avg[i], 1, V[i]->base_hj, 1);
   }

  if(!dparam("PIPE"))
  for(j=0;j<nf;++j)
	V[j]->Set_state('t');

  FILE *fp;
  char fname[BUFSIZ];
  
  sprintf(fname, "%s.1st", omega->name);
  ROOT {
     int err = backup(fname);
     if (err) 
      fprintf(stderr,"WARNING: Analyser: failed to backup the chk file %s\n",
	      fname);
    fp = fopen(fname, "w");
  }
  if (option("parts"))
    pWriteFieldF (fp, omega->name, fname, nsteps, time, nf, V);
  else 
    WritefieldF(fp, omega->name, nsteps, time, nf, V);
  ROOT fclose (fp);

#ifdef THERMO
  dcopy(ntot, ut0, 1, omega->Tf->base_h, 1);
#endif
  free(V);
  return;
}
#endif

static void MEnergy(Domain *omega, FILE *fpenergy, double time)
{
  int nz = option("NZ"), nztot = option("NZTOT");
  register int j, k;
  static double *Em, *GEm;
  double E;
  Element_List *V[3];
  V[0]   = omega->U;  V[1]   = omega->V;  V[2]   = omega->W;
#ifdef MAP
  MapField(omega, 0);
#endif
  
   if (!Em) {
    Em = dvector(0, nz/2-1);
    ROOT GEm = dvector(0, nztot/2-1);
  }

  for (k = 0; k < nz; k += 2) {
    E = 0.;
    for (j = 0; j < 3; j++) {
      E += L2norm(V[j]->flevels[k]);
      if (parid(k+1) != 1) 
	E += L2norm(V[j]->flevels[k+1]);
    }
    Em[k/2] = E;
  }
  Em[0] *= 0.5; // all other modes have a corresponding -N mode 

#ifdef MAP
  for (j = 0; j < 3; j++)
    V[j]->Set_state('t');
 #endif
 
#ifdef PARALLEL
  MPI_Gather (Em, nz/2, MPI_DOUBLE, GEm, nz/2, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
#else
  dcopy(nz/2,Em,1,GEm,1);
#endif

  ROOT {
    fprintf (fpenergy, "%#13.6g", time);
    for(k = 0; k < nztot/2; k++)
      fprintf (fpenergy, "%#13.6g", GEm[k]);
    fprintf (fpenergy, "%#13.6g", dsum(nztot/2, GEm, 1));
    fprintf (fpenergy, "\n");
    
    fflush  (fpenergy);
  }
} 

#ifdef MAP
double VolInt (Element_List *U, double shift)
#else
double VolInt (Element_List *U)
#endif
{
  Element *E = U->fhead;
  double volint = 0., area = 0.;

  for (E = U->fhead; E; E = E->next) {
     const int qa = E->qa, qb = E->qb, qt = E->qtot;
    double *H    = E->h[0];
    double *b    = dvector(0,qt-1),
           *z,*wa,*wb;
    register int i,j;
    
    getzw(qa,&z,&wa,'a');
    switch (E->identify()) {
    case Nek_Tri:
      getzw(qb,&z,&wb,'b');
      break;
    case Nek_Quad:
      getzw(qb,&z,&wb,'a');
      break;
    default:
      fprintf(stderr,"VolInt() not intented for this type of Element!\n");
      break;
    }
    
 #ifndef NEEDS_CHECKING
     for(i = 0; i < qb; ++i)
      dcopy(qa,wa,1,b+i*qa,1);
    for(j = 0; j < qa; ++j)
      dvmul(qb,wb,1,b+j,qa,b+j,qa);
#else
    for (i = 0; i < qb; ++i)
      for(j = 0; j < qa; ++j)
	b[i*qa+j] = wa[j]*wb[i];
#endif
    
    if(E->curvX)
      dvmul(qt, E->geom->jac.p, 1, b, 1, b, 1);
    else
      dscal(qt, E->geom->jac.d, b, 1);
    
#ifdef MAP // this shouldn't change
    area   += dsum(qt, b, 1);
#endif
    volint += ddot(qt, b, 1, H, 1);
    free(b);
  }
  
#ifdef MAP
  return volint + area*shift;
#else
  return volint;
#endif
}
 
double L2norm(Element_List *V){
  int   qt,trip=0;
  double l2=0.0,areat = 0.0,l2a,area;
  double *s,*utmp;
  Element *E;

  utmp = dvector(0,QGmax*QGmax-1);
    
  for(E=V->fhead;E;E = E->next){
    qt = E->qtot;
    s  = E->h[0];
    
    if(E->state == 't'){ E->Trans(E,J_to_Q); trip = 1;}
    dcopy(qt,s,1,utmp,1);
    
    E->Norm_l2m(&l2a,&area);
    
    l2    += l2a; 
     areat += area;
 
    if(trip){ E->state = 't'; trip = 0;}
    else     dcopy(qt,utmp,1,s,1);
  }
  free(utmp);

  return l2/areat;
}

double GL2(Element_List *Div)
{
  int nz = option("NZ"), start = 0;
  register int k;
  double l2 = 0., Gl2 = 1., dz = dparam("LZ")/option("NZTOT");
  
  ROOT { // handle mode 0 & 1 separately
    l2 = L2norm(Div->flevels[0]);
    start = 2;
  }
  for (k = start; k < nz; k++)
    l2 += 2*L2norm(Div->flevels[k]);
  
#ifdef PARALLEL
  MPI_Reduce (&l2, &Gl2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  Gl2 = l2;
#endif
    
   return sqrt(dz*Gl2);
}

void init_avg(Domain *omega)
{
//  fprintf(stderr, "Warning: there is a problem with this function.. \n");
  int nprocs = option("NPROCS");
  int ntot = omega->U->htot;
  int jtot = omega->U->hjtot*omega->U->nz;
  int nf   = (DIM+2); //nf=4
#ifdef THERMO
   nf++; //nf=5;
#endif
  int ntotz, nwork;

  switch (iparam("EQTYPE")) {
  case Rotational:
    ntotz = (ntot+nprocs-1)/nprocs;
    nwork = ntotz*omega->U->nztot;
    break;
  case Convective:
    nwork = ntot*omega->U->nz;
     break;
   default:
    error_msg (2nd moments not implemented for this EQTYPE.);
    break;
  }
  if (option("dealias")) nwork=3*nwork/2;
#if 0
#ifdef MAP
   int npts, i;
  if (npts = omega->mstat->n) {
    if (option("STATAVG")) {
      u2nd = (double **) calloc(npts, sizeof (double *));
      for (i = 0; i < npts; i++) {
	u2nd[i] = dmatrix(0, nf-1, 0, nwork-1);
	dzero(nwork*nf, u2nd[i][0], 1);
      }
    }
    if (option("timeavg")) {
      avg = (double **) calloc(npts, sizeof (double *));
      for (i = 0; i < npts; i++) {
	avg[i] = dmatrix(0, nf-1, 0, jtot-1);
	dzero(jtot*nf, avg[i][0], 1);
     }
    }
  } else
#endif
#endif
  {
  if (option("STATAVG")) {
  if( (dparam("PIPE")) || (dparam("CYLINDER")))
   {
 #ifndef THERMO
    u2nd = dmatrix(0, nf, 0, nwork-1); //nf=4
    dzero(nwork*(nf+1), u2nd[0], 1);
    u2nde= dmatrix(0, nf, 0, nwork-1);
    dzero(nwork*(nf+1), u2nde[0], 1);   
 #else
    u2nd = dmatrix(0, nf+1, 0, nwork-1); //nf=5
    dzero(nwork*(nf+2), u2nd[0], 1);
    u2nde= dmatrix(0, nf+1, 0, nwork-1);
    dzero(nwork*(nf+2), u2nde[0], 1);   
 #endif

   }
  else
   {
    u2nd = dmatrix(0, nf-1, 0, nwork-1);
    dzero(nwork*nf, u2nd[0], 1);
   }
  }
  if (option("timeavg")) {

  if(dparam("PIPE") || dparam("CYLINDER"))
  {
    avg = dmatrix(0, nf-1, 0, nwork-1);
    dzero(nwork*nf, avg[0], 1);
  }
  else
  {
    avg = dmatrix(0, nf-1, 0, jtot-1);
    dzero(jtot*nf, avg[0], 1);
  }


   }
  }
  return;
}

void FFT_AtoA(Element_List *EL, double *base_from);
void AtoA_FFT(Element_List *EL, double *base_to);

void average_u2_avg(Domain *omega, int step, double time)
{
//  fprintf(stderr, "Warning: there is a problem with this function.. \n");
 double       *ut0  =  omega->ut[0],    *vt0   =  omega->vt[0],    *wt0 = omega->wt[0], // used as temps arrary for gradients
              *ut1  =  omega->ut[1],    *vt1   =  omega->vt[1],    *wt1 = omega->wt[1]; // used for calculating laplacians
 double       *us0  =  omega->uss[0],   *vs0   =  omega->vss[0],   *ws0 = omega->wss[0]; // backup arrary

  int i, nf     = DIM+3;

//  if(dparam("PIPE"))
 #ifdef THERMO
    nf = nf+2;
 #endif

  int nsteps    = step;
  Element_List **V;
   Nek_Trans_Type AVGp_to_f = P_to_F;
  
  V = (Element_List**) calloc(nf, sizeof(Element_List*));
  
  V[0] = omega->Uf;  V[1] = omega->Vf; 
  V[2] = omega->Wf;  V[3] = omega->P;
  V[4] = omega->U;

 #ifdef THERMO
  V[5] = omega->V;
  V[6] = omega->Tf;
 #endif

//  if(dparam("PIPE"))
//  {
//    V[4] = omega->Uf;
//  }

  double fac  = 1./((double)nsteps);
  int    ntot = omega->U->hjtot*omega->U->nz,
         nq   = omega->U->htot*omega->U->nz;


  if (option("dealias")) {
    nq   = 3*nq/2;
    AVGp_to_f = P_to_F32;
  }

  //backup
  dcopy(ntot, omega->U->base_hj, 1, us0, 1);
  dcopy(nq,   omega->U->base_h, 1, ut0, 1);
 #ifdef THERMO
  dcopy(ntot, omega->V->base_hj, 1, vs0, 1);
  dcopy(nq, omega->V->base_h, 1, vt0, 1);
  dcopy(ntot, omega->Tf->base_hj, 1, ws0, 1);
  dcopy(nq, omega->Tf->base_h, 1, wt0, 1);
 #endif

#if 0
#ifdef MAP
  int npts, nstats;
  if (npts = omega->mstat->n)
    for (int k = 0; k < npts; k++) {
      if (nstats = omega->mstat->nstat[k]) {
	fac = 1.0/(double) nstats;
	
	switch (iparam("EQTYPE")) {
	case Rotational:
	  for (i = 0; i < nf; i++)
	    FFT_AtoA(V[i], u2nd[k][i]);    
	  break;
	case Convective:
	   for (i = 0; i < nf; i++) {
	    dcopy (nq, u2nd[k][i], 1, V[i]->base_h, 1);
	    V[i]->Set_state('p');
	    V[i]->Trans(V[i], AVGp_to_f);
	  }
	  break;
	default:
	  error_msg (2nd moments not implemented for this EQTYPE.);
	  break;
	}

	for (i = 0; i < nf; i++) {
	  V[i]->Trans(V[i], Q_to_J);
	  dscal(ntot, fac, V[i]->base_hj, 1);
	}
	
	 FILE *fp;
	char fname[BUFSIZ];
	
	sprintf(fname, "%s.p%d.2nd", omega->name, k);
	ROOT {
	  int err = backup(fname);
	  if (err) 
	    fprintf(stderr,
		    "WARNING: Analyser: failed to backup the chk file %s\n", 
		    fname);
	  fp = fopen(fname, "w");
	}
	if (option("parts"))
	  pWriteFieldF (fp, omega->name, fname, nstats, time, nf, V);
	else 
	  WritefieldF(fp, omega->name, nstats, time, nf, V);
	ROOT fclose (fp);
      }
    }
  else
#endif
 #endif
    {
      fac = 1.0/(double) nsteps;

  switch (iparam("EQTYPE")) {
  case Rotational:
    fprintf(stderr, "Warning: there is a problem with Rotational form \n");
    for (i = 0; i < nf; i++)
	FFT_AtoA(V[i], u2nd[i]);  
    break;
  case Convective:
 //zwang 20170101
/*
    for(int j=0; j<nq; ++j)
      {
        u2nd[0][j] =u2nd[0][j]-avg[0][j]*avg[0][j]; // <u*u> - <u>*<u> 
        u2nd[1][j] =u2nd[1][j]-avg[1][j]*avg[1][j]; // <v*v> - <v>*<v>
        u2nd[2][j] =u2nd[2][j]-avg[2][j]*avg[2][j]; // <w*w> - <w>*<w>
        u2nd[3][j] =u2nd[3][j]-avg[3][j]*avg[3][j]; // <p*p> - <p>*<p>
        u2nd[4][j] =u2nd[4][j]-avg[2][j]*avg[1][j]; // <p*p> - <p>*<p>
         
      }
*/
    for (i = 0; i < nf; i++) {

//     if(!dparam("PIPE"))
      dcopy (nq, u2nd[i], 1, V[i]->base_h, 1);
//     else
//      dcopy (nq, u2nde[i], 1, V[i]->base_h, 1);

      V[i]->Set_state('p');
      V[i]->Trans(V[i], AVGp_to_f);
     }
    break;
  default:
    error_msg (2nd moments not implemented for this EQTYPE.);
    break;
  }

  for (i = 0; i < nf; i++) {
    V[i]->Trans(V[i], Q_to_J);
//zwang 20170101
//   if(!dparam("PIPE"))
    dscal(ntot, fac, V[i]->base_hj, 1);
  }

  FILE *fp;
  char fname[BUFSIZ];
 
  sprintf(fname, "%s.2nd", omega->name);
  ROOT {
    int err = backup(fname);
    if (err) 
      fprintf(stderr,"WARNING: Analyser: failed to backup the chk file %s\n", 
	      fname);
    fp = fopen(fname, "w");
  }
  if (option("parts"))
    pWriteFieldF (fp, omega->name, fname, nsteps, time, nf, V);
  else 
    WritefieldF(fp, omega->name, nsteps, time, nf, V);
  ROOT fclose (fp);
 }

  //backup
  dcopy(ntot, us0, 1, omega->U->base_hj, 1);
  dcopy(nq, ut0, 1, omega->U->base_h, 1);
 #ifdef THERMO
  dcopy(ntot, vs0, 1, omega->V->base_hj, 1);
  dcopy(nq, vt0, 1, omega->V->base_h, 1);
  dcopy(ntot, ws0, 1, omega->Tf->base_hj, 1);
  dcopy(nq, wt0, 1, omega->Tf->base_h, 1);
 #endif

  free(V);
  return;
}

void average_u2e_avg(Domain *omega, int step, double time)
{
  int i, nf     = DIM+2;
  int nsteps    = step;
  Element_List **V;
  Nek_Trans_Type AVGp_to_f = P_to_F;
  
  V = (Element_List**) calloc(nf, sizeof(Element_List*));
  
  V[0] = omega->Uf;  V[1] = omega->Vf; 
  V[2] = omega->Wf;  V[3] = omega->P;

  double fac  = 1./((double)nsteps);
  int    ntot = omega->U->hjtot*omega->U->nz,
         nq   = omega->U->htot*omega->U->nz;

  if (option("dealias")) {
    nq   = 3*nq/2;
    AVGp_to_f = P_to_F32;
  }

  switch (iparam("EQTYPE")) {
  case Rotational:
    printf("2nd moment e not implemented for Rotational form.\n");  
    break;
  case Convective:
    for (i = 0; i < nf-1; i++) {
      dcopy (nq, u2nde[i], 1, V[i]->base_h, 1);
       V[i]->Set_state('p');
      V[i]->Trans(V[i], AVGp_to_f);
    }
    break;
  default:
    error_msg (2nd moments not implemented for this EQTYPE.);
    break;
  }

  for (i = 0; i < nf-1; i++) {
    V[i]->Trans(V[i], Q_to_J);
    dscal(ntot, fac, V[i]->base_hj, 1);
  }

  FILE *fp;
  char fname[BUFSIZ];

  sprintf(fname, "%s.sfs", omega->name);
  ROOT {
    int err = backup(fname);
    if (err) 
      fprintf(stderr,"WARNING: Analyser: failed to backup the chk file %s\n", 
	      fname);
    fp = fopen(fname, "w");
  }
  if (option("parts"))
    pWriteFieldF (fp, omega->name, fname, nsteps, time, nf-1, V);
  else 
    WritefieldF(fp, omega->name, nsteps, time, nf-1, V);
  ROOT fclose (fp);

  free(V);
  return;
}
