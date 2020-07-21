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
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "nektarF.h"
#ifdef MAP
#include "map.h"
#endif

/* local functions */
void   BlasiusInterp(Element_List *U, FILE *rea_file, Metric *lambda, 
		     double z);

static void Stokes_Pbsys(Bsystem **Pbsys, Bsystem *Ubsys, 
			 Element_List *U, Bndry *Ubc);

#ifdef THERMO
Bndry **ReadBCs_Temp (FILE *fp, Element_List *EL);
#endif

#ifdef PSE_SLV//modified nadir 15.06.06
static void Set_Oseen(FILE *rea_file,Bsystem *PB, Bsystem *B, 
		      Element_List *U, Domain *omega);
#endif
static Bsystem *gen_bsystem(Element_List *E, Bndry *Ebc);
static void     Summary      (void);
int bandwidth(Element *E, Bsystem *Bsys);

struct ssn_tag {
  char name[FILENAME_MAX];    /* name of session without postfix   */
  char fld [FILENAME_MAX];    /* Output (field) file               */
  char his [FILENAME_MAX];    /* History point file                */
  char rea [FILENAME_MAX];    /* Input file (from pre-processor)   */
  char itp [FILENAME_MAX];    /* Interp. point file                */
  char fce [FILENAME_MAX];    /* Forces file                       */
} session;

static struct {              /* Default options for Nekscal         */
  char *name;
  int   val;
  char *descrip;
} nektar_ops[] = {
  "order"     ,0, "-n #   ... run with # of modes ",
  "binary"    ,1, "-a     ... write output in ascii format",
  "verbose"   ,0, "-v     ... plot divergence error",
  "variable"  ,0, "-var   ... run with variable order (uses [].ord file)",
  "checkpt"   ,0, "-chk   ... dump checkpoint fields at \'iostep\' steps",
#ifndef SAVINGSD
  "iterative" ,0, "-i     ... iterative solve",
  "mixediter" ,0, "-m     ... mixed iterative, press. iterative, vel. direct",
#endif
#ifdef FLOK
  "recursive" ,0, "-r #  ... recursive static condesation with # recursions\n \t -r0 uses a reverse cuthill Mckee option (defaults is -r0)",
#else
  "recursive" ,10, "-r #  ... recursive static condesation with # recursions\n \t -r0 uses a reverse cuthill Mckee option (defaults is -r10)",
#endif
  "timeavg"   ,0, "-T     ... evaluate time average field",
  "tvarying"  ,0, "-V     ... time dependent boundary conditions",
  "dealias"   ,0, "-deal  ... use 3/2 rule dealiasing in the z-direction",
  "STATAVG"   ,0, "-A     ... calculate statistical averages at runtime",
  "SLICES"    ,0, "-S     ... output several time slices for a movie",
  "RAND"      ,0, "-R #   ... run for first # steps with random noise",
  "NZTOT"     ,2, "-z #   ... number of z planes  (power of 2)",
  "THETA"     ,0, "-t #   ... theta C-N time integration parameter",
  "spanval"   ,0, "-span  ... calculate & output forces along the span",
  "parts"     ,0, "-part  ... partial output on each processor",
  "VARV"       ,0, "-LES   ... Large Eddy Simulation",
  "wall_model" ,0, "-wmod  ... Large Eddy Simulation with near-wall model",
  "scalar"     ,0, "",
#ifdef MAP
  "oldupdate"  ,0, "-ou    ... using old update_free() ",
  "acc_impl"   ,0, "-acim  ... implicit treatment of the acceleration terms",
#endif
   0, 0, 0
};

void parse_args(int argc, char *argv[])
{
  int     i;
  char    c;

  /* initialize the parser and install the defaults */
  
  manager_init();

  for (i = 0; nektar_ops[i].name; i++) 
    option_set(nektar_ops[i].name, nektar_ops[i].val);
  
  if(argc == 1) goto showUsage;

  option_set("PROCID", mynode());
  option_set("NPROCS", numnodes());

  while (--argc && (*++argv)[0] == '-') {
    if (strcmp (*argv+1, "chk") == 0) {
      option_set("checkpt",1);
      continue;
    }
    else if (strcmp (*argv+1, "var") == 0) {
      option_set("variable",1);
      continue;
    }
    else if (strcmp (*argv+1, "deal") == 0) {
      option_set("dealias",1);
      continue;
    }
    else if (strcmp (*argv+1, "part") == 0) {
      option_set("parts",1);
      continue;
    }
    else if (strcmp (*argv+1, "LES") == 0) {
      option_set("VARV",1);
      continue;
    }
    else if (strcmp (*argv+1, "span") == 0) {
      option_set("spanval",1);
      continue;
    }
    else if (strcmp (*argv+1, "wmod") == 0) {
      option_set("wall_model",1);
      continue;
    }
#ifdef MAP
    else if (strcmp (*argv+1,"ou") == 0){
	   option_set("oldupdate",1);
 	    continue;
    }
    else if (strcmp (*argv+1, "acim") == 0) {
      option_set("acc_impl",1);
      continue;
    }
#endif
    while (c = *++argv[0])
      switch (c) {
      case 'A':
	option_set("STATAVG",1);
	break;
      case 'a':
	option_set("binary",0);
	break;
      case 'b':
	fprintf(stderr,"Option b is now automatic\n");
	break;
      case 'c':
	{ int n; 
	  if (*++argv[0]) 
	    n = atoi (*argv);
	  else {
	    n = atoi (*++argv);
	    argc--;
	  }
	  option_set("TORDER",n);
	  (*argv)[1] = '\0'; 
	break;
	}
#ifndef SAVINGSD
      case 'i':
	option_set("iterative",1);
	break;
      case 'm':
	option_set("mixediter",1);
	break;
#endif
      case 'n':
	{ int n; 
	  if (*++argv[0]) 
	    n = atoi (*argv);
	  else {
	    n = atoi (*++argv);
	    argc--;
	  }
	  option_set("NORDER.REQ",n);
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'r':
	{ int n; 
	  if (*++argv[0]) 
	    n = atoi (*argv);
	  else {
	    n = atoi (*++argv);
	    argc--;
	  }
	  option_set("recursive",n);
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'R':
	{ int n; 
	  if (*++argv[0]) 
	    n = atoi (*argv);
	  else {
	    n = atoi (*++argv);
	    argc--;
	  }
	  option_set("RAND",n);
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'z':
	{ int nz; 
	  if (*++argv[0]) 
	    nz = atoi (*argv);
	  else {
	    nz = atoi (*++argv);
	    argc--;
	  }
	  option_set("NZTOT.REQ",nz);

	  /* check to see that nz is power of 2 */
#ifdef OLDFFTS
	  while (nz > 2){
	    if(nz%2 != 0) error_msg(-z  must be power of 2);
	    nz /= 2;
	  }
#else
#if !(defined(FLOK) || defined(LNS))
	  if (nz%4 != 0) 
	    ROOT fprintf(stderr, "For dealiasing -z must be a multiple of 4\n");
	  if(nz%(2*option("NPROCS")) != 0) 
	    error_msg(-z must be an multiple of 2 & the number of processors);
#endif
#endif
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'S':
	{
	  option_set("SLICES",1);
	  break;
	}
      case 't':
	{ double n; 
	  if (*++argv[0]) 
	    n = atof (*argv);
	  else {
	    n = atof (*++argv);
	    argc--;
	  }
	  dparam_set("THETA",n);
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'T':
	option_set("timeavg",1);
	break;
      case 'v':
	option_set("verbose",2);
	break;
      case 'V':
	// ce107 note: see note in drive.C
	// ROOT fprintf(stderr, "-V flag not currently supported, see code!\n");
	option_set("tvarying",1);
	break;
      default:
	goto showUsage;
      }
  }
  return;

 showUsage:
  fputs("usage: nektar [options] file[.rea]\n\n"
	"options:\n", stderr);
  for (i = 0; nektar_ops[i].name; i++)
    fprintf(stderr, "%s\n", nektar_ops[i].descrip);
  exit(-1);
  
}
  
/* ------------------------------------------------------------------------ *
 * PreProcess() - Create a new domain                                       *
 *                                                                          *
 * This function creates a fully independent problem domain.                *
 *                                                                          *
 * ------------------------------------------------------------------------ */
void CorrectRobinFluxes(Bndry *B);
void ReadSetLink_T(FILE *fp, Element *U);

Domain *PreProcess(int argc, char **argv)
{
  FILE     *rea_file;
  Element_List   *U,*V,*W,*Uf,*Vf,*Wf,*P;
#ifdef THERMO
  Element_List   *T, *Tf;
#endif
#ifdef SPM
  Element_List   *Pf, *CONCENTR;
#endif
  Element_List   *Ut,*Vt,*Wt,*Pt;
#ifdef OUTPUT_VISCOSITY
  Element_List   *Visc1, *Visc2;
#endif
#ifdef MAP
  Element_List *Wlast;
#endif
  Bndry         **Ubc,**Vbc, **Wbc,**Pbc;
  Bsystem       **Ubsys, **Vbsys, **Wbsys, **Pbsys;
#ifdef THERMO
  Bndry         **Tbc;
  Bsystem       **Tbsys;
#endif 
  Element        *E;
  Domain         *omega;
  int            k,Je,nz,i;
  int            hjtot,htot;
  double         Re,dt;
#ifdef THERMO
  double         Pr;
#endif
  double         **u,**v,**w;
  double         **uf,**vf,**wf;
#ifdef THERMO
  double         **Tt;
  double         **tf;
  double         *ts;
#endif
#ifndef SAVINGSD
  double        *us,*vs,*ws;
#endif
  SLVTYPE        SolveType; 
  ACTION         Eqtype;
#ifdef THERMO
  Bndry          *TBc;	
  Gmap           *gmap_Temp;	
#endif
  /* Create the new domain and open the input files */

  parse_args(argc, argv);

  omega = (Domain*) calloc(1,sizeof(Domain));

  sprintf(session.name,"%s"     , strtok(argv[argc-1],"."));
  sprintf(session.rea, "%s.rea" , argv[argc-1]);
  sprintf(session.fld, "%s.fld" , argv[argc-1]);
  sprintf(session.his, "%s.his" , argv[argc-1]);
  sprintf(session.itp, "%s.itp" , argv[argc-1]);
  sprintf(session.fce, "%s.fce" , argv[argc-1]);

#ifdef DEBUG
#if PARALLEL
  char *buf = (char*) calloc(BUFSIZ,sizeof(char));
  int procid = mynode();
  sprintf(buf, "%s.dbx.%d", argv[argc-1], procid);
  debug_out = fopen(buf, "w");
#else
  debug_out = stdout;
#endif
#endif

  omega->name     = argv[argc-1];
  rea_file        = fopen(session.rea,"r");
  if (rea_file == (FILE*) NULL)
    error_msg(PreProcess(): Could not open input file(s))
  else 
    ROOT printf ("Reading from %s\n\n", session.rea);

  ROOT {
    omega->fld_file = fopen(session.fld,"w");
    omega->his_file = fopen(session.his,"w");
    omega->int_file = fopen(session.itp,"w");
    omega->fce_file = fopen(session.fce,"w"); 
#ifdef PARALLEL
     MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  else {
#ifdef PARALLEL
     MPI_Barrier(MPI_COMM_WORLD);
#endif
    omega->fce_file = fopen(session.fce,"r");
  }

  /* Read the input parameters */
  ReadParams  (rea_file);
  ReadPscals  (rea_file);
  ReadLogics  (rea_file);

  SolveType  = (SLVTYPE)iparam("SLVTYPE"); 
  
  if(nz = iparam("NZTOT"))     option_set("NZTOT",nz); //set nz from param list
  if(nz = option("NZTOT.REQ")) option_set("NZTOT",nz); //set nz from comm line

  init_rfft_struct(); // initialize RFFT structures

#ifdef FLOK3D
  //  if((option("NZTOT") == 2) && (dparam("BETA") == 0)){
  //  fprintf(stdout,"Resetting NZ = 1 since Beta = 0\n");
  //  option_set("NZTOT",1);
  //}
#endif    

  nz = (int)(option("NZTOT")/option("NPROCS"));
  option_set("NZ",nz);

  if(option("tvarying"))
    dparam_set("t",0);

  dparam_set("KINVIS-ORIG",dparam("KINVIS"));

  /* reset implicit viscosity to average VARV */
  if (option("VARV")) 
    {
      if (dparam("AVKINVIS") == 0)
	dparam_set("AVKINVIS",dparam("KINVIS"));
      dparam_set("KINVIS",dparam("AVKINVIS"));
      if (dparam("SMAGS")==0.0)
	{
	  printf("Smagorinsky model not set, setting to 0.1\n");
	  dparam_set("SMAGS",0.1);
	}
    }

//  option_set("recursive", 0);

  Je = iparam("INTYPE");
  Eqtype  = (ACTION) iparam("EQTYPE");
  
  int *size = ivector(0, MAXFACETS-1); //modified nadir 15.06.06
  switch(SolveType){
  case Splitting:
           
    /* Build the mesh */
    U  = ReadMesh(rea_file, session.name);
    
    hjtot = U->nz*U->hjtot;
    htot = U->nz*U->htot;
  
    for(k = 0; k < U->nz; ++k)
      ReadSetLink      (rea_file,U->flevels[k]);
    
    Summary     ();

#ifdef MAP
    InitMap (omega);
#endif
    
    /* set up velocity system using U */
#ifdef MAP
     Ubc = ReadBCs (rea_file,U, omega->mapx->d, omega->mapy->d);
#else
     Ubc = ReadBCs (rea_file,U);      
#endif
  
    /* set up velocity matrix structure  */
    Ubsys = (Bsystem**) malloc(U->nz*sizeof(Bsystem*));
    for(k = 0; k < U->nz; ++k)
      Ubsys[k] = gen_bsystem(U->flevels[k],Ubc[k]); 

    P  = U->gen_aux_field ('p');    
  break;
  case StokesSlv:

    {
#if 0 
      int *size = ivector(0, MAXFACETS-1);
       

#ifdef INVPOW
      iparam_set("NSTEPS",1);
#endif
      option_set("Oseen", 1);
      option_set("Stokes",1);
      option_set("FAMOFF",1);
      
      /* Build the mesh */
      /* set velocity to L+1 */
      iparam_set("MODES",iparam("MODES")+1);  
      U  = ReadMesh(rea_file, session.name);
      
      hjtot = U->nz*U->hjtot;
      htot = U->nz*U->htot;
      
      for(k = 0; k < U->nz; ++k)
	ReadSetLink      (rea_file,U->flevels[k]);
      
      Summary();
      init_ortho_basis();
      
      P = U->gen_aux_field('p');    
      
      option_set("variable",1);
      
      /* set pressure field to L-1 */
      if (P->base_h) {
	free(P->base_h);
	P->base_h = 0;
      }
      if (P->base_hj) {
	free(P->base_hj);
	P->base_hj = 0;
      }
      for(k = 0; k < U->nz; ++k){
	for(E=P->flevels[k]->fhead;E;E=E->next){
	  for(i=0;i<E->Nedges;++i)
	    size[i] = max(0,E->edge[i].l-2);    // MSB: Set dimension of pressure
	  size[i] = max(0,E->face[0].l-2);      // MSB: space to get stability
	  
	  E->Mem_J(size, 'n');
	  E->dgL = E->lmax;
	}
      }
      P->Cat_mem();
      
      
      free(size);
      
      /* set up velocity system using U */
      Ubc = ReadBCs (rea_file,U);      
      
      /* set up velocity matrix structure  */
      Ubsys = (Bsystem**) malloc(U->nz*sizeof(Bsystem*));
      for(k = 0; k < U->nz; ++k){
	Ubsys[k] = gen_bsystem(U->flevels[k],Ubc[k]); 
	Ubsys[k]->lambda = (Metric *) calloc(U->nel, sizeof(Metric)); 
      }
      
      omega->Pf = P->gen_aux_field('p');
#endif
    }
    break;
  default:
    fprintf(stderr,"SolveType not recognised in prepost");
    exit(1);
    break;
  }

  /* set up pressure system */
  Pbc = (Bndry**) malloc(P->nz*sizeof(Bndry*));
  for(k = 0; k < P->nz; ++k)
    Pbc[k]   = BuildPBCs     (P->flevels[k],Ubc[k]);

  Pbsys = (Bsystem**) malloc(P->nz*sizeof(Bsystem*));
  for(k = 0; k < P->nz; ++k)
    Pbsys[k] = gen_bsystem(P->flevels[k],Pbc[k]); 

  // Set up U
  
  u  = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
    u[k] = dvector(0,htot-1);    dzero(htot, u[k], 1);
  }

  Uf    = U->gen_aux_field('u');
  uf    = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
    uf[k] = dvector(0,htot-1);    dzero(htot, uf[k], 1);
  }
#ifndef SAVINGSD
  us = dvector(0, hjtot-1);      dzero(hjtot, us, 1);
#endif  
  // Set up V
  V   = U->gen_aux_field ('v');
  v  = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
    v[k] = dvector(0,htot-1);    dzero(htot, v[k], 1);
  }

  if(option("REFLECT"))
    for(E=V->fhead;E;E=E->next)
      for(i=0;i<E->Nverts;++i)
	E->vert[i].solve = 1;


#ifdef MAP
  Vbc = ReadBCs       (rea_file,V, omega->mapx->d, omega->mapy->d);
#else
  Vbc = ReadBCs       (rea_file,V);  
#endif

  /* set up velocity matrix structure  */
  if(option("REFLECT")){
    Vbsys = (Bsystem**) malloc(U->nz*sizeof(Bsystem*));
    for(k = 0; k < U->nz; ++k)
      Vbsys[k] = gen_bsystem(V->flevels[k],Vbc[k]); 
  }
  else
    Vbsys = Ubsys;

  Vf  = V->gen_aux_field('v');
  vf = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
    vf[k] = dvector(0,htot-1);   dzero(htot, vf[k], 1);
  }
#ifndef SAVINGSD
  vs = dvector(0, hjtot-1);      dzero(hjtot, vs, 1);
#endif
  // Set up W
  W   = U->gen_aux_field ('w');
#ifdef MAP
  Wlast = U->gen_aux_field ('W');
#endif
  w  = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
     w[k] = dvector(0,htot-1);    dzero(htot, w[k], 1);
  }
  
#ifdef MAP
  Wbc = ReadBCs       (rea_file,W, omega->mapx->d, omega->mapy->d);
#else
  Wbc   = ReadBCs(rea_file,W);   
#endif
  Wbsys = Ubsys; 
  
  Wf  = W->gen_aux_field('w');
  wf = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
    wf[k] = dvector(0,htot-1);    dzero(htot, wf[k], 1);
  }

#ifdef THERMO
  T = U->gen_aux_field('t');

  Tt  = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
    Tt[k] = dvector(0,htot-1);    dzero(htot, Tt[k], 1);
  }

  Tf    = U->gen_aux_field('t');
  tf    = (double**) malloc(Je*sizeof(double*));
  for(k = 0; k < Je; ++k){
    tf[k] = dvector(0,htot-1);    dzero(htot, tf[k], 1);
  }

  ts = dvector(0, T->hjtot*T->nz-1);      dzero(T->hjtot*T->nz, ts, 1);
  
  Tbc = ReadBCs_Temp(rea_file, T);      

   /* set up temperature matrix structure  */
   Tbsys  = (Bsystem**) malloc(T->nz*sizeof(Bsystem*));
   for(k = 0; k < T->nz; ++k)
     Tbsys[k] = gen_bsystem(T->flevels[k],Tbc[k]); 
#endif

#ifdef SPM
  Pf        = U->gen_aux_field ('p');  /*set Pf: temperal field*/
  CONCENTR  = U->gen_aux_field ('c');  /*set up phi: SPM particle concentration field*/
#endif

#ifdef OUTPUT_VISCOSITY
  Visc1     = U->gen_aux_field ('f');  /*set Pf: temperal field*/
  Visc2     = U->gen_aux_field ('g');  /*set Pf: temperal field*/
#endif

#ifndef SAVINGSD
  ws = dvector(0, hjtot-1);  dzero(hjtot, ws, 1);
#endif
  omega->soln = NULL;
  ReadSoln(rea_file, omega);

  dt = dparam("DELT");
  Re = 1.0 / dparam("KINVIS");
#ifdef THERMO
  Pr = dparam("PRANDTL");
#endif
  if (iparam("MN"))
    option_set("FAMOFF", 1);
  if (iparam("KFILZ"))
    ROOT printf("Filter applied in the Fourier direction, Delta=%i \n",
	   option("NZTOT")-2*iparam("KFILZ"));
  if (iparam("MN"))
    ROOT printf("********Spectral Viscosity activated MN = %i eps = %lf \n",
	   iparam("MN"),dparam("EPSILON"));
  if (iparam("MNF"))
    ROOT printf("********Spectral Viscosity activated MNF= %i epsF= %lf \n",
	   int(iparam("MNF")),dparam("EPSILONF"));

  for(k = 0; k < U->nz; ++k){
    Ubsys[k]->lambda = (Metric*) calloc(U->nel, sizeof(Metric));

 #ifdef THERMO
    Tbsys[k]->lambda = (Metric*) calloc(T->nel, sizeof(Metric));
 #endif

    if(option("REFLECT"))
      Vbsys[k]->lambda = (Metric*) calloc(U->nel, sizeof(Metric));

    if (iparam("MN"))
      if (dparam("EPSILON") == 0.0)
	if (!(parid(k) % 2))
	  SetEpsilon(rea_file, U->flevels[k], Ubc[k], Ubsys[k]->lambda); // SV
    
    for(i=0;i<U->nel;++i)  {
    if (iparam("MN"))
	      Ubsys[k]->lambda[i].MN      = iparam("MN");
    if (iparam("MNF"))
     {
	     if ((int)(parid(k)/2) <= iparam("MNF"))
	      Ubsys[k]->lambda[i].d = Beta(k)*Beta(k)+ Re*getgamma(1)/dt;
	     else
	      Ubsys[k]->lambda[i].d = dparam("EPSILONF")*
	        exp(-1.0*((int)(parid(k)/2)-(int)(option("NZTOT")/2))*
		      ((int)(parid(k)/2)-(int)(option("NZTOT")/2))/
		      (((int)(parid(k)/2)-iparam("MNF"))*
		      ((int)(parid(k)/2)-iparam("MNF"))))*
	        Beta(k)*Beta(k) + Beta(k)*Beta(k) + Re*getgamma(1)/dt;
     }
   else //No SVV
    {
 	   Ubsys[k]->lambda[i].d = Beta(k)*Beta(k) + Re*getgamma(1)/dt;
    }
      
      Vbsys[k]->lambda[i].d = Ubsys[k]->lambda[i].d;

    #ifdef THERMO
//     double Pr_tmp = Pr;
//     double kr = dparam("RATIO_KAPPA") == 0? 1.:dparam("RATIO_KAPPA");
//     double kappa = dparam("DKAPPA") == 0? 1./Re:dparam("DKAPPA");
//     Pr_tmp = Pr/kr;

//     #ifdef CONJUGATE_HEAT
//
//      int group_id =  T->flevels[k]->flist[i]->group_id; 
//      if(group_id > 1)
//     #endif

 	   Tbsys[k]->lambda[i].d = Beta(k)*Beta(k) + Re*Pr*getgamma(1)/dt;
// 	   Tbsys[k]->lambda[i].d = Beta(k)*Beta(k) + 1./kappa*getgamma(1)/dt;
    #endif
    }

    Pbsys[k]->lambda = (Metric*) calloc(U->nel, sizeof(Metric));
    for(i=0;i<U->nel;++i)
      Pbsys[k]->lambda[i].d = Beta(k)*Beta(k);

    ROOT fprintf(stdout,"Level: %d \n", k);

    if(!(parid(k) % 2)) {
      ROOT fprintf(stdout,"Generating pressure system [."); 
      ROOT fflush(stdout);
      GenMat (P->flevels[k],Pbc[k],Pbsys[k],Pbsys[k]->lambda,Helm);
      ROOT fprintf(stdout,"]\n");
      
      ROOT fprintf(stdout,"Generating velocity system [."); 
      ROOT fflush(stdout);
      GenMat (U->flevels[k],Ubc[k],Ubsys[k],Ubsys[k]->lambda,Helm);
      ROOT fprintf(stdout,"]\n");
 #ifdef THERMO
      ROOT fprintf(stdout,"Generating temperature system [."); 
      ROOT fflush(stdout);
      GenMat (T->flevels[k],Tbc[k],Tbsys[k],Tbsys[k]->lambda,Helm);
      ROOT fprintf(stdout,"]\n");
 #endif

      if(option("REFLECT")){
	ROOT fprintf(stdout,"Generating velocity system (y)[.");
	ROOT fflush(stdout);
	GenMat (V->flevels[k],Vbc[k],Vbsys[k],Vbsys[k]->lambda,Helm);
	ROOT fprintf(stdout,"]\n"); 
      }
      ROOT fflush(stdout);
    }
    else{
      Ubsys[k]->Gmat = Ubsys[k-1]->Gmat; 
      Ubsys[k]->Pmat = Ubsys[k-1]->Pmat;
      Ubsys[k]->rslv = Ubsys[k-1]->rslv;
#ifdef THERMO
      Tbsys[k]->Gmat = Tbsys[k-1]->Gmat; 
      Tbsys[k]->Pmat = Tbsys[k-1]->Pmat;
      Tbsys[k]->rslv = Tbsys[k-1]->rslv;
#endif

      Vbsys[k]->Gmat = Vbsys[k-1]->Gmat; // necc. for symm bc's
      Vbsys[k]->Pmat = Vbsys[k-1]->Pmat;
      Vbsys[k]->rslv = Vbsys[k-1]->rslv;

      Pbsys[k]->Gmat = Pbsys[k-1]->Gmat;
      Pbsys[k]->Pmat = Pbsys[k-1]->Pmat;
      Pbsys[k]->rslv = Pbsys[k-1]->rslv;
    }
  }

  omega->U            = U;    omega->V   = V;   omega->W   = W; 
  omega->Uf           = Uf;   omega->Vf  = Vf;  omega->Wf  = Wf;
  omega->Ubc          = Ubc;  omega->Vbc = Vbc; omega->Wbc = Wbc;
  omega->u            = u;    omega->v   = v;   omega->w   = w;
#ifndef SAVINGSD
  omega->us           = us;   omega->vs  = vs;  omega->ws  = ws;
#endif
  omega->uf           = uf;   omega->vf  = vf;  omega->wf  = wf;	
  omega->Usys         = Ubsys;  
  omega->Vsys         = Vbsys;  
  omega->Wsys         = Wbsys;  

#ifdef MAP
  omega->Wlast        = Wlast;
#endif

  omega->P            = P;
  omega->Pbc          = Pbc;
  omega->Pressure_sys = Pbsys;

#ifdef THERMO
  omega->T            = T; 
  omega->Tf           = Tf;
  omega->Tbc          = Tbc;
  omega->Tt           = Tt;
  omega->tf           = tf;
  omega->ts           = ts;
  omega->Tsys         = Tbsys;
   
#endif

#ifdef SPM
  omega->Pf           = Pf;
  omega->CONCENTR     = CONCENTR;
#endif

#ifdef OUTPUT_VISCOSITY
  omega->visc1        = Visc1;
  omega->visc2        = Visc2;
#endif
  if(SolveType == StokesSlv){
/*
    Bsystem  bsys[3];  // MSB: Modified from 2 to 3 
    
    int l;
    double kinvis = dparam("KINVIS");
    double sigma = dparam("SIGMA");
    
    omega->StkSlv = new StokesMatrix * [nz];

    //set up the Fourier modal Stokes solves
    for(k = 0; k < U->nz; k+=2){
      Vbsys[k]->lambda = Ubsys[k]->lambda;
      Wbsys[k]->lambda = Ubsys[k]->lambda;     // MSB: Wbsys required for 3D
    
      memcpy(bsys  ,Ubsys[k],sizeof(Bsystem));
      memcpy(bsys+1,Vbsys[k],sizeof(Bsystem));
      memcpy(bsys+2,Wbsys[k],sizeof(Bsystem)); // MSB: Wbsys required for 3D
    
      Stokes_Pbsys(&Pbsys[k],bsys,U,Ubc[k]);
    
      omega->Pbc[k] = (Bndry *)NULL;
      Pbsys[k]->lambda = (Metric *) calloc(P->nel, sizeof(Metric)); 
  
#ifdef PSE_SLV
      Set_Oseen(rea_file,Pbsys[k],Ubsys[k],U,omega);  
#endif
      
      for(l = 0; l < U->nel; ++l)
	Ubsys[k]->lambda[l].d = kinvis; 
      
#ifdef INVPOW
      double theta = dparam("THETA");
      for(l = 0; l < U->nel; ++l)
	Pbsys[k]->lambda[l].d = Beta(k)*Beta(k)-sigma*sigma;  
      
#else
      if((Eqtype == StokesS)||(Eqtype == Oseen)){
	for(l = 0; l < U->nel; ++l) Pbsys[k]->lambda[l].d = 0.0;      
      }
      else{
	double theta = dparam("THETA");
	
	for(l = 0; l < U->nel; ++l)
	  Pbsys[k]->lambda[l].d = (1.0/dt/kinvis/(1.0-theta))+Beta(k)*Beta(k)-
	    sigma*sigma;
      }
#endif 
    
      fprintf(stdout,"Generating stokes system [(%d)",k); fflush(stdout);
      
      omega->StkSlv[k] = new StokesMatrix();
      omega->StkSlv[k]->GenMat(U,P,Ubsys[k],Pbsys[k],Pbsys[k]->lambda,Beta(k));
      
      fprintf(stdout,"]\n"); fflush(stdout);
      rewind(rea_file);
    }
 */
  }

#ifndef FLOK
  setup_transfer_space(omega);
  TransformBCs(omega, P_to_F);
#endif

  ReadICs     (rea_file, omega);


#ifdef MAP
  if (fabs(omega->mapx->time-dparam("STARTIME")) > 1e-6) {
    ROOT fprintf(stderr, "readMap and Restart gave different start time!\n");
    if (dparam("STARTIME") != 0.0) {
      // use the restart time from the field restart file
      omega->mapx->time = dparam("STARTIME");
      omega->mapy->time = omega->mapx->time;
    } else
      dparam_set("STARTIME", omega->mapx->time);
  }
#endif

  ReadDF      (rea_file, DIM+1);
  ReadDFunc   (rea_file, omega);
//#ifdef MAP
#if defined(MAP) && !defined(SPM)
  int forcx = (int) dparam("FORCX");
  int forcy = (int) dparam("FORCY");
  if(!option("oldupdate") && ((forcx<0) || (forcy<0)) )
   {
     ReadTension  (omega);  
     ReadStiffne  (omega);
   }
#endif

#ifdef FLOW_CONTROL 
  ReadHilbertCoeffs (omega);
#endif

  ReadHisData (rea_file, omega);
#ifdef FLOK
  ReadIntData (rea_file, omega);
#endif

  if(option("STATAVG") || option("timeavg"))
    init_avg(omega);
  
  if(option("GivenICs")){
#ifdef MAP
    ResetICs  (omega); // reset the ics using the mapping

    if (option("dealias")) { // Need to get back nz planes
      U->Trans(U, P_to_F32);
      V->Trans(V, P_to_F32);
      W->Trans(W, P_to_F32);
 #ifdef THERMO
      T->Trans(T, P_to_F32);
 #endif
      
      Wlast->Trans(Wlast, P_to_F);
      Wlast->Trans(Wlast, F_to_P32);
    } else
#endif
     {
	    U->Trans(U, P_to_F);
	    V->Trans(V, P_to_F);
 	    W->Trans(W, P_to_F);
 #ifdef THERMO
      T->Trans(T, P_to_F);
 #endif
     }
    P->Trans(P, P_to_F);
  }

#ifdef ZEROOUT
  ROOT {
    dzero(U->htot, U->base_h+U->htot, 1);
    dzero(V->htot, V->base_h+V->htot, 1);
    dzero(W->htot, W->base_h+W->htot, 1);
 #ifdef THERMO
    dzero(T->htot, T->base_h+T->htot, 1);
 #endif
  }
#endif

//#ifdef ENTROPYVISCOSITY
  omega->velocity_range = dmatrix(0, 3-1, 0, 2-1);

	Ut = U->gen_aux_field('r');
  Vt = U->gen_aux_field('s');
  Wt = U->gen_aux_field('t');
  Pt = U->gen_aux_field('q');
  
  omega->Ut = Ut;   omega->Vt  = Vt;  omega->Wt  = Wt; omega->Pt = Pt;
//#endif
//
   int worksize = Ut->htot*Ut->nz;
   if (option("dealias")) {
     worksize = 3*htot/2;
   }

#ifdef ENTROPYVISCOSITY
   omega->visc_tmp = dvector(0,worksize-1);
   omega->visc_ave = dvector(0,U->htot-1);
   omega->visc_bak = dvector(0,U->htot-1);
   omega->visc_max_indices = dvector(0,worksize-1);
#endif
//  if(Je<2)
//   {
//    fprintf(stdout,"We are using second-order time-splitting scheme ... \n");
//    exit(1);
//   }

   int JJ = 2;
	 omega->uk  = dmatrix(0,JJ-1,0,worksize-1);
	 omega->vk  = dmatrix(0,JJ-1,0,worksize-1);
	 omega->wk  = dmatrix(0,JJ-1,0,worksize-1);
	 omega->pk  = dmatrix(0,JJ-1,0,worksize-1);
 #ifdef THERMO
	 omega->tk  = dmatrix(0,JJ-1,0,worksize-1);
 #endif

	 omega->ut  = dmatrix(0,JJ-1,0,worksize-1);
	 omega->vt  = dmatrix(0,JJ-1,0,worksize-1);
	 omega->wt  = dmatrix(0,JJ-1,0,worksize-1);
	    
   omega->up  = dmatrix(0,JJ-1,0,worksize-1);
	 omega->vp  = dmatrix(0,JJ-1,0,worksize-1);
	 omega->wp  = dmatrix(0,JJ-1,0,worksize-1);
	    
   omega->uss  = dmatrix(0,JJ-1,0,hjtot-1);
	 omega->vss  = dmatrix(0,JJ-1,0,hjtot-1);
	 omega->wss  = dmatrix(0,JJ-1,0,hjtot-1);
	    
//#if defined(OUTPUT_VISCOSITY) || defined(WALE)
   omega->st  = dmatrix(0,9-1,0,worksize-1);
//#endif
    omega->afx  = dmatrix(0,JJ-1,0,worksize-1); // here we use this array store the accerlating term relating to coordinate transformation
    omega->afy  = dmatrix(0,JJ-1,0,worksize-1); // here we use this array store the accerlating term relating to coordinate transformation
    omega->afz  = dmatrix(0,JJ-1,0,worksize-1); // here we use this array store the accerlating term relating to coordinate transformation
#ifdef THERMO
    omega->aft  = dmatrix(0,JJ-1,0,worksize-1); // here we use this array store the accerlating term relating to coordinate transformation
#endif

//   omega->Q = dmatrix(0,11,0,htot-1);
   if(dparam("STATISTICS"))
    {
     if(dparam("CHANNEL"))
        statistics_setup(omega);
      
     if(dparam("CYLINDER")||dparam("PIPE")||dparam("SPHERE"))
        point_output_prepare(omega);

    }


  fclose(rea_file);
  return omega;
}

void PostProcess(Domain *omega, int step, double time){
  Element_List **V;
  FILE  *fld_file = omega->fld_file;
  int nfields;
  nfields = DIM+1;

  ++nfields;  // add W
#ifdef MAP
  nfields++; // add Wlast
#endif  
#ifdef THERMO
  nfields++;
#endif

  V = (Element_List**) malloc(nfields*sizeof(Element_List*));
  
  V[0] = omega->U;
  V[1] = omega->V;
  V[2] = omega->W;
  V[3] = omega->P;
#ifdef MAP
  V[4] = omega->Wlast;
#endif
#ifdef THERMO
  V[5] = omega->T;
#endif
  
#ifdef MAP
  ROOT WriteMap (omega);
  Nek_Trans_Type p_to_f = P_to_F;
  if (option("dealias"))
    p_to_f = P_to_F32;

  V[4]->Trans(V[4], p_to_f); // Get us back to Fourier Space
  V[4]->Trans(V[4], Q_to_J); // Get us back to modal space
#endif
  if (option("parts"))
    pWriteFieldF (fld_file, omega->name, session.fld, step, time, nfields, V);
  else 
    WritefieldF(fld_file, omega->name, step, time, nfields, V);

  ROOT {
    fclose(fld_file);
    fclose(omega->his_file);
    fclose(omega->int_file);
    fclose(omega->fce_file);
  }

  if(option("timeavg"))
    averagefields(omega, step, time);
  if(option("STATAVG"))
    {
      average_u2_avg(omega, step, time);
#ifdef SMAGORINSKY
      if (option("VARV"))
	average_u2e_avg(omega, step, time);
#endif
    }
  return;
}

/* ----------------------------------------------------------------------- *
 * Summary() - Print a summary of the input data                           *
 *                                                                         *
 * This function collects the echo of the input data which used to be      *
 * scattered throughout this file.                                         *
 *                                                                         *
 * ----------------------------------------------------------------------- */

static void Summary (void){
  ROOT{  
    printf ("Input File          : %s\n",session.name);
    printf ("Reynolds number     : %g\n",1./dparam("KINVIS-ORIG"));
#ifdef THERMO
    printf ("Prandtl number      : %g\n",dparam("PRANDTL"));
    printf ("Gravity direction   : %d\n",iparam("DIRGRAVITY"));
#endif
    printf ("NZ                  : %d\n",option("NZ"));
    printf ("NZTOT               : %d\n",option("NZTOT"));
#ifdef FLOK
    printf ("Beta                : %lf\n",dparam("BETA"));
#else
    printf ("LZ                  : %lf\n",dparam("LZ"));
#endif
    printf ("Time step           : %g\n",dparam("DELT"));
    printf ("Integration order   : %d\n",iparam("INTYPE"));
    if(!option("variable"))
      printf("Number of modes     : %d\n",iparam("MODES"));
    else
      printf("Number of modes     : variable\n");
    printf ("Number of elements  : %d\n",iparam("ELEMENTS"));
    printf ("Number of Families  : %d\n",iparam("FAMILIES"));    
    {
      fputs ("Equation type       : ", stdout);
#ifdef FLOK
      puts("Linearised Eigenvalue Solver");
#elif defined(LNS)
      puts("Linearised Solver");
      iparam_set("EQTYPE",3);
#else
      switch (iparam("EQTYPE")) {
      case Rotational:
	puts("Navier-Stokes (rotational)");
	break;
      case Convective:
	puts("Navier-Stokes (convective)");
	break;
      case Stokes:
	puts("Stokes flow");
	break;
      default:
	puts("undefined");
	break;
      }
#endif
    }
    if(option("RAND"))
      printf("Random forcing      : %d steps for k < %d, amp = (%g %g %g)\n",
	     option("RAND"), iparam("KNOISE"), dparam("XNOISE"), 
	     dparam("YNOISE"), dparam("ZNOISE"));
    
    printf("Integration time    : %g, or %d steps\n", 
	   dparam("FINTIME"), iparam("NSTEPS"));
    printf("I/O time for saves  : %g, or %d steps",
	   dparam("IOTIME"),  iparam("IOSTEP"));
    
    if  (option("checkpt"))
      fputs (" [checkpoint]\n", stdout);
    else { 
      putchar ('\n');
      if(iparam("NSTEPS"))
	if (iparam("NSTEPS") / iparam("IOSTEP") > 10)
	  fputs ("Summary: " "You have more than 10 dumps..."
		 "did you want to use -chk?\n", stderr);
    }
#ifdef MAP
    int forcx = (int) dparam("FORCX");
    int forcy = (int) dparam("FORCY");
    double PI = 4.*atan(1.);
    double m_a = PI/4.;
   if( (forcx<0) || (forcy<0) )
    {
     printf ("Ma                  : %lf\n",m_a);
    double ee = dparam("ZMASS");
    printf ("Mc/Ma               : %lf\n",ee);
    double stiff = dparam("WNB")*dparam("WNB")*ee;
    printf ("Stiffness           : %lf\n",stiff);
    double tens = dparam("WNC")*dparam("WNC")*ee;
    printf ("Tension             : %lf\n",tens);
    double u_bulk = dparam("U_BULK");
    printf ("Mean velocity       : %lf\n",u_bulk);

    double Lz = dparam("LZ");
    printf ("u                   : %lf\n",sqrt(m_a/stiff)*Lz*u_bulk);
    printf ("Ip                  : %lf\n",m_a*u_bulk*u_bulk/(tens*(1+ee)/(ee-0.125)));

    for(int n=1; n<5; ++n)
     {
      double  freq = (double)n*(double)n*PI*sqrt(stiff/(m_a+m_a*ee))/2/Lz/Lz
         *sqrt(1.+tens*Lz*Lz/stiff/PI/PI/(double)n/(double)n);

      printf("Natrual frequency %d mode : %lf \n", n, freq);
     }
   }
#endif

#ifdef ENTROPYVISCOSITY
    double  CA    = dparam("STABILIZATION_ALPHA");
    double  CB    = dparam("STABILIZATION_BETA");
    printf ("Entropy Viscosity Method:  \n");
    printf ("alpha: %lf  beta: %lf \n", CA,CB);

    if(dparam("VAN_DRIEST"))
     {
       printf ("Van Driest damping function:  \n");
       printf ("Damping length: %lf \n",dparam("DAMPING_LENGTH"));
     }
#endif

  }
  return;
}



/* This is a function to sort out the global numbering scheme for      *
 * vertices, edges and faces. It also sets up the list of cumulative   *
 * indices for edges and vertices.                                     */

static void set_bmap       (Element *, Bsystem *);
static int  suminterior    (Element *E);
static void BasicNumScheme (Element_List *E, Bndry *Ebc, Bsystem *Bsys);

static Bsystem *gen_bsystem(Element_List *UL, Bndry *Ebc){
  Bsystem  *Bsys;
  Element  *E=UL->fhead;

  Bsys       = (Bsystem *)calloc(1,sizeof(Bsystem));

  if(option("iterative")) Bsys->smeth = iterative;

  /* basis numbering scheme for vertices, edges and faces */
  BasicNumScheme(UL,Ebc,Bsys);

  if(option("recursive"))
    Mlevel_SC_decom(UL,Bsys);
  else  if(Bsys->smeth == direct){
    bandwidthopt(E,Bsys);
    ROOT
      fprintf(stdout,"rcm bandwidth (%c)   : %d [%d (%d)] \n",E->type,
	      bandwidth(E,Bsys),Bsys->nsolve,suminterior(E));
    fflush(stdout);
  }
  
  /* set up boundary maps for element edges */
  set_bmap(E,Bsys);
  Bsys->families = iparam("FAMILIES");
  return Bsys;
}

static void  setGid(Element_List *E);

static void BasicNumScheme(Element_List *UL, Bndry *Ebc, Bsystem *Bsys){
  register int i;
  int       nvg,nvs,neg,nes,nfg,nfs,scnt,ncnt;
  int      *gsolve,*gbmap,l,l1;
  Bndry    *Ubc; 
  Vert     *v;
  Edge     *e;
#if DIM == 3
  Face     *f;
  extern   int ednum[][3];
#endif
  Element *E, *U=UL->fhead;

  nfg = 0; /* compiler trick */

  setGid (UL);  /* setup a global numbering scheme i.e. without boundaries*/

  /* This part of the routine re-orders the vertices, edges and then
     the faces so that the knowns are listed first. For the edges and
     the faces it also initialises a cummalative list stored in Bsys. */
  
  /*--------------------*/
  /* Vertex re-ordering */
  /*--------------------*/

  /* find maximum number of global vertices; */
  for(E=U,nvg=0; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      nvg = (E->vert[i].gid > nvg)? E->vert[i].gid:nvg;
  ++nvg;
  

  gsolve = ivector(0,nvg-1);
  gbmap  = ivector(0,nvg-1);

  ifill(nvg, 1, gsolve,1);
  /* Assemble vertex solve mask to sort out multiplicity */
  for(E=U; E; E = E->next){
    for(i = 0; i < E->Nverts; ++i)
      gsolve[E->vert[i].gid] &= E->vert[i].solve;
  }
  
  /* copy back mask */
  for(E=U; E; E = E->next){
    v = E->vert;
    for(i = 0; i < E->Nverts; ++i)
     v[i].solve = gsolve[E->vert[i].gid];
  }
  
  scnt = 0;  ncnt = 0;
  for(i = 0; i < nvg; ++i)
     gbmap[i] = gsolve[i]? scnt++:ncnt++;
  nvs = scnt;

  /* place unknowns at the end of the pile */
  for(i = 0; i < nvg; i++) 
    gbmap[i] += gsolve[i]?  0:nvs;

  /* replace vertices numbering into vertex structures */
  for(E=U; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      E->vert[i].gid = gbmap[E->vert[i].gid];

  free(gsolve); free(gbmap);

  /*------------------*/
  /* Edge re-ordering */
  /*------------------*/

  /* find maximum number of global edges; */
  for(E = U,neg=0; E; E = E->next)
    for(i = 0; i < E->Nedges; ++i)
      neg = (E->edge[i].gid > neg)? E->edge[i].gid:neg;
  ++neg;

  gsolve = ivector(0,neg-1);
  gbmap  = ivector(0,neg-1);

  /* form edge gsolve and gbmap */
  ifill(neg, 1, gsolve,1);
  for(Ubc = Ebc; Ubc; Ubc = Ubc->next)
    if((Ubc->type == 'V')||(Ubc->type == 'W')||(Ubc->type == 'o')){
      l = Ubc->face;
#if DIM == 2
      gsolve[Ubc->elmt->edge[l].gid] = 0;
#else
      for(i = 0; i < 3; ++i)
	gsolve[Ubc->elmt->edge[ednum[l][i]].gid] = 0;
#endif
    }
  
  scnt = 0;  ncnt = 0;
  for(i = 0; i < neg; ++i)
    gbmap[i] = gsolve[i]? (scnt++):(ncnt++);
  nes   = scnt;

  /* place unknowns at the end of the pile */
  for(i = 0; i < neg; ++i)
    gbmap[i] += gsolve[i]?  0:nes;

  /* replace sort gid's */
  for(E = U; E; E = E->next){
    e = E->edge;
    for(i = 0; i < E->Nedges; ++i)
      e[i].gid = gbmap[E->edge[i].gid];
  }
  
  /* set up cumulative edge list */
  Bsys->edge = ivector(0,neg);
  for(E = U; E; E = E->next)
    for(i = 0; i < E->Nedges; ++i)
      Bsys->edge[E->edge[i].gid] = E->edge[i].l;
  
  l = Bsys->edge[0];
  Bsys->edge[0] = 0;
  for(i = 1; i < neg; ++i){
    l1 = Bsys->edge[i];
    Bsys->edge[i] =  Bsys->edge[i-1]+l;
    l = l1;
  }
  Bsys->edge[i] =  Bsys->edge[i-1]+l;

  free(gsolve); free(gbmap);

#if DIM == 3

  /*------------------*/
  /* Face re-ordering */
  /*------------------*/

  /* find maximum number of global faces; */
  for(k = 0,nfg =0; k < nel; ++k)
    for(i = 0; i < Nfaces; ++i)
      nfg = (E[k].face[i].gid > nfg)? E[k].face[i].gid:nfg;
  ++nfg;
  
  gsolve = ivector(0,nfg-1);
  gbmap  = ivector(0,nfg-1);

  /* form faces part of gsolve, gbmap */
  ifill(nfg, 1, gsolve,1);
  for(Ubc = Ebc; Ubc; Ubc = Ubc->next)
    if((Ubc->type == 'V')||(Ubc->type == 'W')||(Ubc->type == 'o'))
      gsolve[Ubc->elmt->face[Ubc->face].gid] = 0;

  scnt = 0; ncnt = 0;
  for(i = 0; i < nfg; ++i)
    gbmap[i] = gsolve[i]? (scnt++):(ncnt++);
  nfs = scnt;

  /* place unknowns at the end of the pile */
  for(i = 0; i < nfg; ++i)
    gbmap[i] += gsolve[i]?  0:nfs;

  /* replace sorted gid's */
  for(k = 0; k < nel; ++k){
    f=E[k].face;
    for(i = 0; i < Nface; ++i)
      f[i].gid = gbmap[f[i].gid];
  }

  /* set up cumulative face list */
  Bsys->face = ivector(0,nfg);
  for(k = 0;k < nel; ++k){
    f = E[k].face;
    for(i = 0; i < Nface; ++i)
      Bsys->face[f[i].gid] = f[i].l;
  }

  l = Bsys->face[0];
  Bsys->face[0] = 0;
  for(i = 1; i < nfg; ++i){
    l1 = Bsys->face[i];
    Bsys->face[i] =  Bsys->face[i-1] + l*(l+1)/2;
    l = l1;
  }
  Bsys->face[i] = Bsys->face[i-1] + l*(l+1)/2;

  free(gsolve); free(gbmap);
#endif

  Bsys->nv_solve  = nvs;
  Bsys->ne_solve  = nes;
  Bsys->nel       = countelements(U);
  Bsys->nsolve    = nvs + Bsys->edge[nes];
  Bsys->nglobal   = nvg + Bsys->edge[neg];
#if DIM == 3
  Bsys->nsolve   += Bsys->face[nfs];
  Bsys->nglobal  += Bsys->face[nfg];
  Bsys->nf_solve  = nfs;
#endif

  /* make numbering scheme consequative the unknowns are listed by
     vertices, edges and then faces */

  /* place known degree of freedom at end of list */
  for(E = U; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      if(E->vert[i].gid >= nvs)
	E->vert[i].gid += Bsys->nsolve-nvs;
  
#if DIM == 3  
  for(i = 0; i < nfs; ++i)
    Bsys->face[i] += nvs + Bsys->edge[nes];
  
  ncnt = Bsys->nsolve + (nvg-nvs) +
    (Bsys->edge[neg]-Bsys->edge[nes]) - Bsys->face[nfs];
  for(i = nfs; i <= nfg; ++i)
    Bsys->face[i] += ncnt;
#endif
  
  for(i = 0; i < nes; ++i)
    Bsys->edge[i] += nvs;
  
  ncnt = Bsys->nsolve + nvg-nvs - Bsys->edge[nes];
  for(i = nes; i <= neg; ++i)
    Bsys->edge[i] += ncnt;

}

typedef struct vertnum {
  int    id; 
  struct vertnum *base;
  struct vertnum *link;
} Vertnum;


#if DIM == 2
  int tri_Vnum  [][2] = {{0,1},{1,2},{2,0}};
  int tri_Vnum1 [][2] = {{1,0},{2,1},{0,2}};

  int quad_Vnum  [][2] = {{0,1},{1,2},{2,3},{3,0}};
  int quad_Vnum1 [][2] = {{1,0},{2,1},{3,2},{0,3}};
#else
  int Vnum  [][3] = {{1,0,2},{0,1,3},{1,2,3},{2,0,3}};
  int Vnum1 [][3] = {{0,1,2},{1,0,3},{2,1,3},{0,2,3}};
#endif

static void setGid(Element_List *UL){
  register int i,j;
  int      nvg, face, eid, edgeid, faceid, edgmax;
  const int    nel = UL->nel;
  Edge     *e,*ed;
  Face     *f;
  Vertnum  *V,*vb,*v;
  Element *E, *U=UL->fhead;

  /* set vector of consequative numbers */
  /* set up vertex list */
  edgmax = 4;

  V = (Vertnum *) calloc(edgmax*nel,sizeof(Vertnum));

  for(E = U; E; E = E->next)
    for(i = 0; i < E->Nedges; ++i){
      for(j = 0; j < DIM; ++j){
	if(E->Nedges == 3)
	  v   = V + E->id*edgmax +  tri_Vnum[i][j];
	if(E->Nedges == 4)
	  v   = V + E->id*edgmax + quad_Vnum[i][j];

	if(E->edge[i].base){
	  if(E->edge[i].link){
	    eid  = E->edge[i].link->eid;
	    face = E->edge[i].link->id; 
	  }
	  else{
	    eid  = E->edge[i].base->eid;
	    face = E->edge[i].base->id;
	  }

	  if(UL->flist[eid]->Nedges == 3)
	    vb  = V[eid*edgmax +  tri_Vnum1[face][j]].base;
	  if(UL->flist[eid]->Nedges == 4)
	    vb  = V[eid*edgmax + quad_Vnum1[face][j]].base;
	  
	  if(eid < E->id){  /* connect to lower element */
	    if(!v->base) v->base = v;
	    
	    /* search through all points and assign to same base */
	    for(;vb->link;vb = vb->link);
	    if(vb->base != v->base) vb->link = v->base;
	    for(v = v->base;v; v = v->link) v->base = vb->base;
	  }
	  else if(!v->base) v->base = v;
	}
	else if(!v->base) v->base = v;
      }
    }

  /* number vertices consequatively */
  for(E = U, nvg = 1; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      if(!V[E->id*edgmax+i].base->id) 
	V[E->id*edgmax+i].id = V[E->id*edgmax+i].base->id = nvg++;
      else                       
	V[E->id*edgmax+i].id = V[E->id*edgmax+i].base->id;
  nvg--;

  for(E = U; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      E->vert[i].gid = V[E->id*edgmax+i].id-1;

  
  /* set gid's to -1 */
  for(E = U; E; E = E->next){
    for(i = 0; i < E->Nedges; ++i) E->edge[i].gid = -1;
    for(i = 0; i < E->Nfaces; ++i) E->face[i].gid = -1;
  }

  /* at present just number edge and faces consequatively */
  faceid = 0; edgeid = 0;
  for(E = U; E; E = E->next){
    e = E->edge;
    for(i = 0; i < E->Nedges; ++i){
      if(e[i].gid==-1)
	if(e[i].base){
	  for(ed = e[i].base; ed; ed = ed->link)
	    ed->gid = edgeid;
	  ++edgeid;
	}
	else
	  e[i].gid = edgeid++;
    }
    f = E->face;
    for(i = 0; i < E->Nfaces; ++i)
      if(f[i].gid==-1)
	if(f[i].link){
	  f[i].gid       = faceid;
	  f[i].link->gid = faceid++;
	}
	else
	  f[i].gid       = faceid++;
  }

  free(V);
  return;
}


static int suminterior(Element *E){
  int sum=0;
  for(;E;E = E->next){
#if DIM == 2
    if(E->Nedges == 3)
      sum += E->face->l*(E->face->l+1)/2;
    else
      sum += E->face->l*E->face->l;
#else
    sum += E->l*(E->l+1)*(E->l+2)/6;
#endif
  }
  return sum;
}

static void set_bmap(Element *U, Bsystem *B){
  register int i,j,k,n;
  int   l;
  const int nel = B->nel;
  int  **bmap;
  Element *E;

  /* declare memory */
  bmap = (int **) calloc(nel,sizeof(int *));
  for(E=U,l=0;E;E=E->next) l += E->Nbmodes;

  bmap[0] = ivector(0,l-1);
  for(i = 0, E=U; i < nel-1; ++i, E=E->next)
    bmap[i+1] = bmap[i] + E->Nbmodes;

  /* fill with bmaps */
  for(E=U;E;E=E->next){
    for(j = 0; j < E->Nverts; ++j)
      bmap[E->id][j] = E->vert[j].gid;
    
    for(j = 0,n = E->Nverts; j < E->Nedges; ++j,n+=l){
      l = E->edge[j].l;
      for(k = 0; k < l; ++k)
	bmap[E->id][n+k] = B->edge[E->edge[j].gid] + k;
    }
  }

  B->bmap = bmap;
}

/* ---------------------------------------------------------------------- *
 * backup() -- Create a backup copy of a file                             *
 * ---------------------------------------------------------------------- */

int backup (char *path1)
{
  int  stat;
  char path2[FILENAME_MAX];

  sprintf (path2, "%s.bak", path1);
  unlink  (path2);                    /* unlink path2 regardless    */
  if (!(stat = link(path1, path2)))   /* try to link path1 -> path2 */
    unlink (path1);                   /* unlink path1 only if the   */
  return stat;                        /* link was sucessful         */
}


static void Stokes_Pbsys(Bsystem **Sys, Bsystem *Ubsys, Element_List *U,
			 Bndry *Ubc){
  
  register int i,j,d;
  int   l,cnt;
  const int nel = Ubsys->nel;
  int   **bmap;
  Element *E;
  Bsystem *Pbsys;

  Sys[0] = Pbsys = (Bsystem *)calloc(1,sizeof(Bsystem));
  Pbsys->Gmat    = (MatSys  *)malloc(  sizeof(MatSys) );

  Pbsys->nsolve  = Ubsys->nsolve  + Ubsys[1].nsolve  + Ubsys[2].nsolve  + nel;
  Pbsys->nglobal = Ubsys->nglobal + Ubsys[1].nglobal + Ubsys[2].nglobal + nel;

  Pbsys->nel     = Ubsys->nel;
  Pbsys->smeth   = Ubsys->smeth;
  Pbsys->Precon  = Ubsys->Precon;

  /* declare memory */
  bmap = (int **) malloc(nel*sizeof(int *));
  for(E=U->fhead,l=0;E;E=E->next) l += E->Nbmodes;
  l = (DIM+1)*l + nel;
  bmap[0] = ivector(0,l-1);

  for(i = 0, E=U->fhead; i < nel-1; ++i, E=E->next)
    bmap[i+1] = bmap[i] + (DIM+1)*E->Nbmodes + 1;

  int offset = 1; // offset bmap so that pressure is listed first

  // MSB: Modified for 3D Stokes Solve <= DIM
  
  for(E=U->fhead;E;E=E->next){
    int nslv_space = 0;
    int nkwn_space = 0;
    for(j = 0; j <= DIM; ++j){
      for(i = 0; i < E->Nbmodes; ++i){
	if(Ubsys[j].bmap[E->id][i] < Ubsys[j].nsolve)
	  bmap[E->id][i +j*E->Nbmodes+offset] =
	    nslv_space + Ubsys[j].bmap[E->id][i];
	else
	  bmap[E->id][i+j*E->Nbmodes+offset] = 
	    (Ubsys[j].bmap[E->id][i]-Ubsys[j].nsolve) +
	    nkwn_space + Pbsys->nsolve;
      }
      nslv_space += Ubsys[j].nsolve;
      nkwn_space += Ubsys[j].nglobal - Ubsys[j].nsolve;
    }

    /* add pressure dof */ 
    bmap[E->id][0] = Ubsys[0].nsolve + Ubsys[1].nsolve
                   + Ubsys[2].nsolve + E->id;
  }
  
  /* Check to see if an outflow is defined - if so turn off
     singular pressure mode */
  Pbsys->singular = 1;
  for(;Ubc;Ubc = Ubc->next)
    if((Ubc->type == 'O')||(Ubc->type == 'F')){
      Pbsys->singular = 0;
      break;
    }
#ifdef UNSTEADY_PSE
  Pbsys->singular = 0;
#endif  
  if(Pbsys->singular){
    
    // put iesing element at end of global matrix system 
    for(i = U->nel-1; i > iparam("IESING"); --i)
      bmap[i][0] = bmap[i-1][0];
    bmap[iparam("IESING")][0] = Ubsys[0].nsolve + Ubsys[1].nsolve 
                              + Ubsys[2].nsolve + nel - 1; 
  }
  
  Pbsys->bmap = bmap;

  /* reset recursive numbering system and set up in Pbsys */
  if(Ubsys->rslv){
    register int k,n;
    Rsolver *Vrslv = Ubsys->rslv;
    Rsolver *Prslv;
    int nrecur     = Vrslv->nrecur;
    Recur *Vrdat   = Vrslv->rdata;
    Recur *Prdat;
    int len,cstart,top,alen,maxid,eid,maplen,*nelmt,*elmtid,*nelmt1,*elmtid1;
    int *mapping,*nptchold,*n2optch, *swap, *swap1;

    elmtid    = ivector(0,nel-1);
    nelmt     = ivector(0,nel);
    elmtid1   = ivector(0,nel-1);
    nelmt1    = ivector(0,Vrdat[0].npatch);

    n2optch   = ivector(0,nel-1);
    nptchold  = ivector(0,Vrdat[0].npatch);

    mapping = ivector(0,Pbsys->nsolve-1);

    /* set up elmt numbering for original patch i.e. individual elements */
    top = 0; eid = 1;
    iramp(nel+1,&top,&eid,nelmt,1);
    iramp(nel  ,&top,&eid,elmtid,1);

    Pbsys->rslv   = Prslv = (Rsolver *)calloc(1,sizeof(Rsolver));
    Prslv->nrecur = nrecur; 

    Prdat = Prslv->rdata = (Recur *)calloc(nrecur,sizeof(Recur));
    
    /* declare rdata in Pbsys for new system */
    for(i = 0; i < nrecur; ++i){
      ifill(Pbsys->nsolve,-1,mapping,1);

      Prdat[i].npatch = Vrdat[i].npatch;
      Prdat[i].patchlen_a = ivector(0,Prdat[i].npatch-1);
      Prdat[i].patchlen_c = ivector(0,Prdat[i].npatch-1);
      Prdat[i].map = (int **)malloc(Vrdat[i].npatch*sizeof(int *));
      Prdat[i].pmap = Vrdat[i].pmap;

      if(i){
	top = Prdat[i-1].cstart;
	maplen = Prdat[i-1].npatch;
      }
      else{
	top = Pbsys->nsolve;
	maplen = nel;
      }

      /* determine new to old mapping of patches  */
      nptchold[0] = 0;
      for(j = 0; j < Vrdat[i].npatch; ++j){
	for(k = 0,cnt=0; k < maplen; ++k)
	  if(Vrdat[i].pmap[k] == j)
	    n2optch[nptchold[j]+ cnt++] = k;
	nptchold[j+1] = nptchold[j] + cnt;
      }

      nelmt1[0] = 0;
      for(j = 0; j < Vrdat[i].npatch; ++j){
	nelmt1[j+1] = nelmt1[j];
	for(k=nptchold[j]; k < nptchold[j+1]; ++k){
	  len = nelmt[n2optch[k]+1]-nelmt[n2optch[k]];
	  icopy(len,elmtid+nelmt[n2optch[k]],1,elmtid1+nelmt1[j+1],1);
	  nelmt1[j+1] +=  len;
	}
      }
      
      swap    = nelmt;
      swap1   = elmtid;
      nelmt   = nelmt1;
      elmtid  = elmtid1;
      nelmt1  = swap;
      elmtid1 = swap1;
      
      /* reshuffle bmaps so that pressure dof are lumped with patch blocks */
      cstart = DIM*Vrdat[i].cstart;
      
      /*!!!!!!!!!!!!!!!!

      Note that bmap should probably be Prdat->bmap in here. This needs
      checking throught propertly. 

      !!!!!!!!!!*/

      cnt = Vrdat[i].npatch;
      for(j = 0; j < Vrdat[i].npatch; ++j){

	for(k = nelmt[j],maxid=0; k < nelmt[j+1]; ++k){
	  eid  = elmtid[k];
	  alen = DIM*U->flist[eid]->Nbmodes+1;
	
	  for(n = 0; n < alen-1; ++n)
	    if((bmap[eid][n] >= cstart)&&(bmap[eid][n] < top)){
	      mapping[bmap[eid][n]] = bmap[eid][n] + cnt; 
	      bmap[eid][n] += cnt;
	      maxid = max(bmap[eid][n],maxid);
	    }
	}
	/* Redefine pressure dof putting first element at
	   beginning of cstart so it is carried over to boundary
	   solve. Otherwise define the dof to be in the interior
	   block */
	
	eid  = elmtid1[nelmt1[n2optch[nptchold[j]]]];
	alen = DIM*U->flist[eid]->Nbmodes+1;
	mapping[bmap[eid][alen-1]] = cstart+j;
	bmap[eid][alen-1] = cstart+j;
	
	for(k = nptchold[j]+1,n=1; k < nptchold[j+1]; ++k,++n){
	  eid  = elmtid1[nelmt1[n2optch[k]]];
	  alen = DIM*U->flist[eid]->Nbmodes+1;
	  mapping[bmap[eid][alen-1]] = maxid+n;
	  bmap[eid][alen-1] = maxid+n;
	}
                    
	len  = nptchold[j+1]-nptchold[j]; 
	cnt += len-1;

	Prdat[i].patchlen_a[j] = DIM*Vrdat[i].patchlen_a[j] + 1;
	Prslv->max_asize       = max(Prslv->max_asize,Prdat[i].patchlen_a[j]);
	Prdat[i].map[j]        = ivector(0,Prdat[i].patchlen_a[j]-1);
	Prdat[i].patchlen_c[j] = DIM*Vrdat[i].patchlen_c[j] + len-1;
      }
      
      Prdat[i].cstart = cstart + Prdat[i].npatch;
      
      /* fill out mapping for each patch to new local A system */
      cstart = Prdat[i].cstart;

      /* go through and redefine previous Prdat[i].map which will have
	 been altered due to new recursion */
      for(j = 0; j < i; ++j){
	for(k = 0; k < Prdat[j].npatch; ++k)
	  for(n = 0; n < Prdat[j].patchlen_a[k]; ++n)
	    if(mapping[Prdat[j].map[k][n]]+1)
	      Prdat[j].map[k][n] = mapping[Prdat[j].map[k][n]];
      }

      for(j = 0; j < Prdat[i].npatch; ++j){
	ifill(cstart,-1,mapping,1);
	for(k = nelmt[j],maxid=0; k < nelmt[j+1]; ++k){
	  eid  = elmtid[k];
	  alen = DIM*U->flist[eid]->Nbmodes+1;
	
	  for(n = 0; n < alen; ++n)
	    if(bmap[eid][n] < cstart)
	      mapping[bmap[eid][n]] = bmap[eid][n];
	}

	for(k = 0,cnt=0; k < cstart; ++k)
	  if(mapping[k]+1)
	    Prdat[i].map[j][cnt++] = mapping[k];

	if(cnt != Prdat[i].patchlen_a[j])
	  error_msg(error in reset_pbmap);
      }
    }
    /* find the patch id that bmap[IEsing] belongs to */
    k = iparam("IESING");
    for(i = 0; i < nrecur; ++i)
      k = Prdat[i].pmap[k];
      /* work out  new patch id */

    if(k >= Prdat[nrecur-1].npatch) 
      error_msg(error sorting singular pressure in Pbsys_bmap);
					      
    /* set pressure singular mode */
    Pbsys->singular = Prdat[nrecur-1].cstart-k;

    free(mapping); 
    free(nelmt);    free(elmtid);
    free(nelmt1);   free(elmtid1);
    free(nptchold); free(n2optch);
  }
}


#ifdef PSE_SLV //Nadir modified nadir 15.06.06
  /* setup Oseen solver */
static void Set_Oseen(FILE *rea_file,Bsystem *PB, Bsystem *B, 
		      Element_List *U, Domain *omega){
  
  register int i,j;
  int      eDIM = U->fhead->dim(),qtot;
  eDIM++;
  
//#ifdef BLASIUS    // MSB:  1 for current, 2 for d/dz, 3 for initial f(\eta)
  eDIM = 3*eDIM+1;  // MSB: +1 for storage
//#endif
  
  double   **wave,*s;

  wave = (double **)malloc(eDIM*B->nel*sizeof(double*));
  s = dvector(0,eDIM*U->htot);
  dzero(eDIM*U->htot,s,1);

  for(i = 0; i < B->nel; ++i){
    qtot = U->flist[i]->qtot;
    for(j = 0; j < eDIM; ++j, s+= qtot){
      wave[i*eDIM   + j] = s;
    }
    B ->lambda[i].wave = wave + i*eDIM;
    PB->lambda[i].wave = wave + i*eDIM;
  }
  
  // MSB: Note for Blasius eDIM in ReadWave has changed also
  ReadWave(rea_file,wave,U,omega);

  /* divide through by the viscosity to account for way that 
     matrix is formed */

  dscal(eDIM*U->htot,1.0/dparam("KINVIS"),wave[0],1);

#ifdef BTVPG  // MSB: Copy initial data into whole structure
  double kinvis =  dparam("KINVIS");

  for(i = 0; i < B->nel; ++i){
    qtot = U->flist[i]->qtot;
    
    dcopy(qtot, B->lambda[i].wave[0], 1, B->lambda[i].wave[6], 1);
    dcopy(qtot, B->lambda[i].wave[1], 1, B->lambda[i].wave[7], 1);
    dcopy(qtot, B->lambda[i].wave[2], 1, B->lambda[i].wave[8], 1);
 
    dcopy(qtot, PB->lambda[i].wave[0], 1, PB->lambda[i].wave[6], 1);
    dcopy(qtot, PB->lambda[i].wave[1], 1, PB->lambda[i].wave[7], 1);
    dcopy(qtot, PB->lambda[i].wave[2], 1, PB->lambda[i].wave[8], 1);
  }

  Element *NE;
  Coord X;

  X.x      = dvector(0, QGmax*QGmax-1); 
  X.y      = dvector(0, QGmax*QGmax-1); 

  double *zvec, *ivec, *zzpyy2, *zzpyy, *zzmyy, R, offset, z;

  R        = dparam("RADIUS");  // MSB: Radius of circle for pressure gradient
  offset   = dparam("OFFSET");  // MSB: Offset of circle for pressure gradient
  z        = dparam("Z0");
 
  zvec    = dvector(0, QGmax*QGmax-1); 
  ivec    = dvector(0, QGmax*QGmax-1); 
  zzpyy2  = dvector(0, QGmax*QGmax-1);
  zzpyy   = dvector(0, QGmax*QGmax-1);
  zzmyy   = dvector(0, QGmax*QGmax-1);

  dfill(QGmax*QGmax, z + offset, zvec, 1);
  dfill(QGmax*QGmax, 1.0, ivec, 1);  

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;
    NE->coord(&X);
   
    // MSB: Form y^2 coordinate vector
    dvmul(qtot, X.y, 1, X.y, 1, X.x, 1);
    
    // MSB: Store (z^2 - y^2)
    dvvtvm(qtot, zvec, 1, zvec, 1, X.x, 1, zzmyy, 1);
    // MSB: Store (z^2 + y^2)
    dvvtvp(qtot, zvec, 1, zvec, 1, X.x, 1, zzpyy, 1);    
    // MSB: Store (z^2 + y^2)^2
    dvmul(qtot, zzpyy, 1, zzpyy, 1, zzpyy2, 1);
   
    // MSB: Store V
    dvmul(qtot, zvec, 1, X.y, 1, B->lambda[NE->id].wave[9], 1);
    dvdiv(qtot, B->lambda[NE->id].wave[9], 1, zzpyy2, 1,
	  B->lambda[NE->id].wave[1], 1);
    dscal(qtot, -2.0*R*R/kinvis, B->lambda[NE->id].wave[1], 1); // NB: KINVIS

    dvadd(qtot, B->lambda[NE->id].wave[1], 1, B->lambda[NE->id].wave[7], 1,
	  B->lambda[NE->id].wave[1], 1);	  
    
    // MSB: Store W
    dvdiv(qtot, zzmyy, 1, zzpyy2, 1, B->lambda[NE->id].wave[2], 1);
    dscal(qtot, -R*R/kinvis, B->lambda[NE->id].wave[2], 1); // NB: KINVIS 

    dvadd(qtot, B->lambda[NE->id].wave[2], 1, B->lambda[NE->id].wave[8], 1,
	  B->lambda[NE->id].wave[2], 1);
    
    // MSB: Store dV/dz   
    dvmul(qtot, zvec, 1, zvec, 1, B->lambda[NE->id].wave[9], 1);
    dvdiv(qtot, B->lambda[NE->id].wave[9], 1, zzpyy, 1, 
	  B->lambda[NE->id].wave[9], 1);
    dscal(qtot, 4.0, B->lambda[NE->id].wave[9], 1); 
    dvsub(qtot, ivec, 1, B->lambda[NE->id].wave[9], 1,
	  B->lambda[NE->id].wave[9], 1);
    dvmul(qtot, X.y, 1, B->lambda[NE->id].wave[9], 1,
	  B->lambda[NE->id].wave[9], 1);
    dvdiv(qtot, B->lambda[NE->id].wave[9], 1, zzpyy2, 1,
	  B->lambda[NE->id].wave[4], 1);
    dscal(qtot, -2.0*R*R/kinvis, B->lambda[NE->id].wave[4], 1); // NB: KINVIS
    
    // MSB: Store dW/dz   
    dvdiv(qtot, zzmyy, 1, zzpyy, 1, B->lambda[NE->id].wave[9], 1);
    dscal(qtot, 2.0, B->lambda[NE->id].wave[9], 1); 
    dvsub(qtot, ivec, 1, B->lambda[NE->id].wave[9], 1,
	  B->lambda[NE->id].wave[9], 1);
    dvmul(qtot, zvec, 1, B->lambda[NE->id].wave[9], 1,
	  B->lambda[NE->id].wave[9], 1);
    dvdiv(qtot, B->lambda[NE->id].wave[9], 1, zzpyy2, 1,
	  B->lambda[NE->id].wave[5], 1);
    dscal(qtot, -2.0*R*R/kinvis, B->lambda[NE->id].wave[5], 1); // NB: KINVIS
 
    dcopy(qtot, B->lambda[NE->id].wave[1], 1, PB->lambda[NE->id].wave[1], 1); 
    dcopy(qtot, B->lambda[NE->id].wave[2], 1, PB->lambda[NE->id].wave[2], 1); 
    dcopy(qtot, B->lambda[NE->id].wave[4], 1, PB->lambda[NE->id].wave[4], 1); 
    dcopy(qtot, B->lambda[NE->id].wave[5], 1, PB->lambda[NE->id].wave[5], 1); 

    dzero(qtot,  B->lambda[NE->id].wave[3], 1);
    dzero(qtot, PB->lambda[NE->id].wave[3], 1);    
  }

  free(X.x); free(X.y);
#endif
  
  // MSB: Note on BLASIUS -----------------------------------------------------
  //
  // lambda[nel].wave[0] <=> U (mean cmpnt)
  // lambda[nel].wave[1] <=> V (mean cmpnt)
  // lambda[nel].wave[1] <=> W (mean cmpnt)
  //
  // lambda[nel].wave[0] <=> dU/dz (mean cmpnt)
  // lambda[nel].wave[1] <=> dV/dz (mean cmpnt)
  // lambda[nel].wave[1] <=> dW/dz (mean cmpnt)
  //
  // lambda[nel].wave[0] <=> U(eta) (initial Blasius profile)
  // lambda[nel].wave[1] <=> V(eta) (initial Blasius profile)
  // lambda[nel].wave[1] <=> W(eta) (initial Blasius profile)
  // --------------------------------------------------------------------------
  
#ifdef BLASIUS // MSB: Copy initial data into whole structure
  for(i = 0; i < B->nel; ++i){
    qtot = U->flist[i]->qtot;
    
    dcopy(qtot, B->lambda[i].wave[0], 1, B->lambda[i].wave[3], 1);
    dcopy(qtot, B->lambda[i].wave[1], 1, B->lambda[i].wave[4], 1);
    dcopy(qtot, B->lambda[i].wave[2], 1, B->lambda[i].wave[5], 1);
 
    dcopy(qtot, B->lambda[i].wave[0], 1, B->lambda[i].wave[6], 1);
    dcopy(qtot, B->lambda[i].wave[1], 1, B->lambda[i].wave[7], 1);
    dcopy(qtot, B->lambda[i].wave[2], 1, B->lambda[i].wave[8], 1);
 
    dcopy(qtot, PB->lambda[i].wave[0], 1, PB->lambda[i].wave[3], 1);
    dcopy(qtot, PB->lambda[i].wave[1], 1, PB->lambda[i].wave[4], 1);
    dcopy(qtot, PB->lambda[i].wave[2], 1, PB->lambda[i].wave[5], 1);
 
    dcopy(qtot, PB->lambda[i].wave[0], 1, PB->lambda[i].wave[6], 1);
    dcopy(qtot, PB->lambda[i].wave[1], 1, PB->lambda[i].wave[7], 1);
    dcopy(qtot, PB->lambda[i].wave[2], 1, PB->lambda[i].wave[8], 1);

    // MSB:************************************************
    /*
    dscal(qtot, dparam("KINVIS"), B->lambda[i].wave[6], 1);
    dscal(qtot, dparam("KINVIS"), B->lambda[i].wave[7], 1);    
    dscal(qtot, dparam("KINVIS"), B->lambda[i].wave[8], 1);
    */
    // MSB:************************************************
  }

  // MSB: i.e. wave[6] contains f(\eta) and wave[7] contains f'(\eta)
  // MSB: and wave[8] contains f''(\eta) - these should NOT be modified

  // MSB: Now make sure that wave[1] contains V and wave[2] W
  // MSB: and wave[4] contains dV/dz and wave[5] dW/dz
  // MSB: and wave[9] acts as storage for manipulation
  // MSB: this is assuming that W_\infty = 1
  
  double z = dparam("Z0");
  double kinvis =  dparam("KINVIS");
  double scal1  =  0.5*sqrt(kinvis/z);
  double scal2  = -0.25*sqrt(kinvis)*pow(z, -1.5); 
 
  fprintf(stderr, " \n scal1 = %f \n\n", scal1);
  fprintf(stderr, " \n scal2 = %f \n\n", scal2);

  Element *NE;
  Coord X, detadz; 

  // MSB: For these purposes X.y stores \eta co-ordinates
  // MSB: and detadz.y stores d\eta/dz

  X.x      = dvector(0, QGmax*QGmax-1); 
  X.y      = dvector(0, QGmax*QGmax-1); 

  detadz.x = dvector(0, QGmax*QGmax-1); 
  detadz.y = dvector(0, QGmax*QGmax-1); 

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;
    
    NE->coord(&X);
    NE->coord(&detadz);
    dscal(qtot, -0.5/z, detadz.y, 1);

    // MSB: Store W = W_\infty*f'(\eta)
    dcopy(qtot, B->lambda[NE->id].wave[7], 1,  B->lambda[NE->id].wave[2], 1);
    dcopy(qtot,PB->lambda[NE->id].wave[7], 1, PB->lambda[NE->id].wave[2], 1); 
 
    // MSB: Store V (see notebook for equation)
    dvvtvm(qtot, X.y, 1, B->lambda[NE->id].wave[7], 1,
	   B->lambda[NE->id].wave[6], 1, B->lambda[NE->id].wave[1], 1);   
    dscal(qtot, scal1, B->lambda[NE->id].wave[1], 1);   
    
    dvvtvm(qtot, X.y, 1, PB->lambda[NE->id].wave[7], 1,
	   PB->lambda[NE->id].wave[6], 1, PB->lambda[NE->id].wave[1], 1);   
    dscal(qtot, scal1, PB->lambda[NE->id].wave[1], 1); 

    // MSB: Zero U component
    dzero(qtot,  B->lambda[NE->id].wave[0], 1);
    dzero(qtot, PB->lambda[NE->id].wave[0], 1); 

    // MSB: Zero dU/dz
    dzero(qtot,  B->lambda[NE->id].wave[3], 1);
    dzero(qtot, PB->lambda[NE->id].wave[3], 1); 

    // MSB: Store dV/dz
    dvvtvm(qtot, X.y, 1, B->lambda[NE->id].wave[7], 1,
	   B->lambda[NE->id].wave[6], 1, B->lambda[NE->id].wave[4], 1);
    dscal(qtot, scal2, B->lambda[NE->id].wave[4], 1);   
    dvmul(qtot, X.y, 1, B->lambda[NE->id].wave[8], 1,
	  B->lambda[NE->id].wave[9], 1);    
    dvmul(qtot, detadz.y, 1, B->lambda[NE->id].wave[9], 1,
	  B->lambda[NE->id].wave[9], 1);    
    dsvtvp(qtot, scal1, B->lambda[NE->id].wave[9], 1,
	   B->lambda[NE->id].wave[4], 1, B->lambda[NE->id].wave[4], 1);

    dvvtvm(qtot, X.y, 1, PB->lambda[NE->id].wave[7], 1,
	   PB->lambda[NE->id].wave[6], 1, PB->lambda[NE->id].wave[4], 1);
    dscal(qtot, scal2, PB->lambda[NE->id].wave[4], 1); 
    dvmul(qtot, X.y, 1, PB->lambda[NE->id].wave[8], 1,
	  PB->lambda[NE->id].wave[9], 1);    
    dvmul(qtot, detadz.y, 1, PB->lambda[NE->id].wave[9], 1,
	  PB->lambda[NE->id].wave[9], 1);    
    dsvtvp(qtot, scal1, PB->lambda[NE->id].wave[9], 1,
	   PB->lambda[NE->id].wave[4], 1, PB->lambda[NE->id].wave[4], 1);

    // MSB: Store dW/dz
    dvmul(qtot, detadz.y, 1, B->lambda[NE->id].wave[8], 1,
	  B->lambda[NE->id].wave[5], 1);  
    dvmul(qtot, detadz.y, 1, PB->lambda[NE->id].wave[8], 1,
	  PB->lambda[NE->id].wave[5], 1);       


    // MSB:********************************************************
    /*
    dscal(qtot, 1./dparam("KINVIS"), B->lambda[NE->id].wave[0], 1);
    dscal(qtot, 1./dparam("KINVIS"), B->lambda[NE->id].wave[1], 1);    
    dscal(qtot, 1./dparam("KINVIS"), B->lambda[NE->id].wave[2], 1);

    dscal(qtot, 1./dparam("KINVIS"), B->lambda[NE->id].wave[3], 1);
    dscal(qtot, 1./dparam("KINVIS"), B->lambda[NE->id].wave[4], 1);    
    dscal(qtot, 1./dparam("KINVIS"), B->lambda[NE->id].wave[5], 1);
    */
    // MSB:********************************************************


    // MSB: IMPORTANT========================================================
    // Divide through by the viscosity to account for 
    // the way that the matrices are formed 
    /*
    dscal(qtot,1.0/dparam("KINVIS"),B->lambda[NE->id].wave[0],1);
    dscal(qtot,1.0/dparam("KINVIS"),B->lambda[NE->id].wave[1],1);
    dscal(qtot,1.0/dparam("KINVIS"),B->lambda[NE->id].wave[2],1);
    dscal(qtot,1.0/dparam("KINVIS"),B->lambda[NE->id].wave[3],1);
    dscal(qtot,1.0/dparam("KINVIS"),B->lambda[NE->id].wave[4],1);
    dscal(qtot,1.0/dparam("KINVIS"),B->lambda[NE->id].wave[5],1);

    dscal(qtot,1.0/dparam("KINVIS"),PB->lambda[NE->id].wave[0],1);
    dscal(qtot,1.0/dparam("KINVIS"),PB->lambda[NE->id].wave[1],1);
    dscal(qtot,1.0/dparam("KINVIS"),PB->lambda[NE->id].wave[2],1);
    dscal(qtot,1.0/dparam("KINVIS"),PB->lambda[NE->id].wave[3],1);
    dscal(qtot,1.0/dparam("KINVIS"),PB->lambda[NE->id].wave[4],1);
    dscal(qtot,1.0/dparam("KINVIS"),PB->lambda[NE->id].wave[5],1);
    */
    // MSB: IMPORTANT========================================================
    
  }
  
  free(X.x); free(X.y); free(detadz.x); free(detadz.y);

  // MSB:======================================================================
  // INTERPOLATION DONE HERE
  // MSB:======================================================================

  Element_List *Uw;

  // MSB:------------------------------------------------------------
  // Interpolate V
  // MSB:------------------------------------------------------------

  Uw  = U->gen_aux_field('u');
  
  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot, B->lambda[NE->id].wave[1], 1, Uw->flist[NE->id]->h[0], 1);
  }
  
  // MSB:-------------------------------------------
  BlasiusInterp(Uw, rea_file, B->lambda, z);
  // MSB:-------------------------------------------

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot,  B->lambda[NE->id].wave[9], 1,  B->lambda[NE->id].wave[1], 1);
  }

  // MSB:------------------------------------------------------------
  // Interpolate W
  // MSB:------------------------------------------------------------

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot, B->lambda[NE->id].wave[2], 1, Uw->flist[NE->id]->h[0], 1);
  }
  
  // MSB:-------------------------------------------
  BlasiusInterp(Uw, rea_file, B->lambda, z);
  // MSB:-------------------------------------------

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot,  B->lambda[NE->id].wave[9], 1,  B->lambda[NE->id].wave[2], 1);
  }

  // MSB:------------------------------------------------------------
  // Interpolate dV/dz
  // MSB:------------------------------------------------------------

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot, B->lambda[NE->id].wave[4], 1, Uw->flist[NE->id]->h[0], 1);
  }
  
  // MSB:-------------------------------------------
  BlasiusInterp(Uw, rea_file, B->lambda, z);
  // MSB:-------------------------------------------

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot, B->lambda[NE->id].wave[9], 1, B->lambda[NE->id].wave[4], 1);
  }

  // MSB:------------------------------------------------------------
  // Interpolate dW/dz
  // MSB:------------------------------------------------------------

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot, B->lambda[NE->id].wave[5], 1, Uw->flist[NE->id]->h[0], 1);
  }
  
  // MSB:-------------------------------------------
  BlasiusInterp(Uw, rea_file, B->lambda, z);
  // MSB:-------------------------------------------

  for(NE = U->fhead; NE; NE = NE->next){
    qtot = NE->qa*NE->qb;   
    dcopy(qtot, B->lambda[NE->id].wave[9], 1, B->lambda[NE->id].wave[5], 1);
  }

  fprintf(stderr, "------------------------------------------------- \n"  );
  fprintf(stderr, "Blasius profile interpolated \n"   );
  fprintf(stderr, "------------------------------------------------- \n"  );
#else
#ifndef BTVPG
  for(i = 0; i < B->nel; ++i){
    qtot = U->flist[i]->qtot;

    dzero(qtot, B->lambda[i].wave[3], 1);
    dzero(qtot, B->lambda[i].wave[4], 1);
    dzero(qtot, B->lambda[i].wave[5], 1);
    dzero(qtot, B->lambda[i].wave[6], 1);
    dzero(qtot, B->lambda[i].wave[7], 1);
    dzero(qtot, B->lambda[i].wave[8], 1);

    dzero(qtot, PB->lambda[i].wave[3], 1);
    dzero(qtot, PB->lambda[i].wave[4], 1);
    dzero(qtot, PB->lambda[i].wave[5], 1);
    dzero(qtot, PB->lambda[i].wave[6], 1);
    dzero(qtot, PB->lambda[i].wave[7], 1);
    dzero(qtot, PB->lambda[i].wave[8], 1);

  }
#endif
#endif

}
#endif


/* ------------------------------------------------------------------------ */
/* Blasius solution interpolation                                           */
/* ------------------------------------------------------------------------ */

void BlasiusInterp(Element_List *U, FILE *rea_file, Metric *lambda, double z){

  /* Scale Blasius solution and store in lambda->wave */

  Element *E, *NE; 
  Coord Xe, *A;
  double *hr = dvector(0, QGmax-1);
  double *hs = dvector(0, QGmax-1);
  int npts, *eids;

  Xe.x = dvector(0, QGmax*QGmax-1);
  Xe.y = dvector(0, QGmax*QGmax-1); 

  double scal = 1./sqrt(dparam("KINVIS")*z); 
  // MSB: CHECK SCALING: SHOULD IT BE 1./???

  fprintf(stderr, " \n scal  = %f \n\n", scal);
  
  for(NE = U->fhead; NE; NE = NE->next){
    NE->coord(&Xe);  
    fprintf(stderr,"Interpolating element %d:",NE->id+1);
    npts = NE->qa*NE->qb;
    
    for(int n = 0; n < npts; n++){
      Xe.x[n] = 1.0; // MSB: Constant value to remove interpolation issues
      Xe.y[n] = scal*Xe.y[n];
    }
    
    Find_local_coords(U,&Xe,npts,&eids,&A);
    
    for(int k = 0; k < npts; ++k){
      if(eids[k] != -1){
	E = U->flist[eids[k]];	
	get_point_shape_2d(E, A->x[k], A->y[k], hr, hs);
	/* Copy data into wave */
	(lambda+NE->id)->wave[9][k] =
	  eval_field_at_pt_2d(E->qa,E->qb, U->flist[E->id]->h[0], hr, hs);
      }
    }
  }
  free(Xe.x); free(Xe.y); 
}

/* ------------------------------------------------------------------------ */
