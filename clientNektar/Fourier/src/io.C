/*
 * Input functions based on the NEKTON spectral element data model
 */
/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $RCSfile: io.C,v $
 * $Revision: 1.6 $
 * $Author: ssherw $
 * $Date: 2006/06/07 18:11:14 $
 * $State: Exp $
 * ------------------------------------------------------------------------- */

// set up for 2d only


#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <veclib.h>
#include "nektarF.h"
#include "Quad.h"
#include "Tri.h"
#include <rfftw.h>
extern rfftw_plan rplan, rplan_inv;
                                              /* .....  Section names .....  */
#define  SEC_PARAMS   sections[0]             /* Parameters                  */
#define  SEC_PSCALS   sections[1]             /* Passive scalars             */
#define  SEC_LOGICS   sections[2]             /* Logical switches            */
#define  SEC_MESH     sections[3]             /* Mesh definition             */
#define  SEC_CURVE    sections[4]             /* Curved sides                */
#define  SEC_BCS      sections[5]             /* Boundary conditions         */
#define  SEC_ICS      sections[6]             /* Initial conditions          */
#define  SEC_DF       sections[7]             /* Drive force                 */
#define  SEC_HIST     sections[8]             /* History points              */
#define  SEC_T_BCS    sections[9]             /* Thermo boundary conditions  */
#define  SEC_INT      sections[10]            /* Interp. points              */
#define  SEC_WAVE     sections[11]            /* Wave Prop in Oseen eqn      */
#define  SEC_DFUNCS   sections[12]            /* Drive force                 */
#define  SEC_TENS     sections[13]            /* Tension                     */
#define  SEC_STIF     sections[14]            /* Bending stiffness           */
#define  SEC_LAMB     sections[15]            /* Heat Conductivities         */
#define  SEC_RHO_CP   sections[16]            /* Density times Heat Capacity */

static int   checklist   (char *name, char *list[]);
static char *findSection (char *name, char *buf, FILE *fp);
static int   bedgcmp     (const void *b1, const void *b2);

static char *sections[] = { 
  "Paramaters",  "Passive Scalars",     "LOGICAL",     "MESH", "Curved sides",
  "BOUNDARY",    "INITIAL CONDITIONS",  "DRIVE FORCE",  "HISTORY", "THERMAL",
  "INTERPOLATION", "WAVE","DRIVE FUNCT","TENSION",  "STIFFNESS", "LAMBDAS","RHO_CP"};

/* ------------------------------------------------------------------------ *
 * Spectral Element Data Model                                              *
 *                                                                          *
 * This file contains functions to read a spectral element input file in    *
 * the NEKTON format.  The file format is fairly specific, expecting the    *
 * following sections in order:                                             *
 *                                                                          *
 *     PARAMETERS -> Passive scalars -> Logical switches -> MESH      ->    *
 *     BC's       -> INITIAL COND.'s -> DRIVE FORCE      -> Variables ->    *
 *     HISTORY POINTS                                                       *
 *                                                                          *
 * The uppercase sections are the ones which are actually used.  In the     *
 * parameters section, the following should be defined if you are solving   *
 * Stokes or Navier-Stokes:                                                 *
 *                                                                          *
 *     KINVIS              Viscosity, or 1/Re                               *
 *     NSTEP | FINTIME     Number of steps or simulation time.              *
 *     DELT                Time step                                        *
 *                                                                          *
 * Return value: none                                                       *
 * ------------------------------------------------------------------------ */

void ReadParams (FILE *fp)
{
  int  n;
  char buf[BUFSIZ], value[25], name[25], c;
  double val, dt;

  static char *dspecial[] = { "KINVIS", "LAMBDA", "IOTIME", "KC", "LZ",
                               "FRNDX","FRNDY","FRNDZ" , "PECLET", 0 };
  static char *ispecial[] = { "TEQN", "EQTYPE", "SLVTYPE",
			      "CFLSTEP", "HISSTEP","STASTEP","DIRGRAVITY", 0 };
  rewind (fp);
  for (n = 0; n < 4; n++) fgets(buf, BUFSIZ, fp);

  if (sscanf(buf, "%d", &n) != 1) 
    {fputs("ReadParams: can't read # of parameters", stderr);exit(-1);}
  
  while (n--) {
    fgets (buf, BUFSIZ, fp);
    if(sscanf(buf, "%25s%25s", value, name) == 2) {
      if (checklist (name, dspecial))
	 dparam_set (name, scalar(value));
      else if (checklist (name, ispecial))
	iparam_set (name, (int) scalar(value));
      else if (isupper(c = *name) && 'I' <= c && c <= 'N')
	iparam_set(name, (int) scalar(value));
      else
	dparam_set(name, scalar(value));
    }
  }
  

  /* The following section is only for unsteady problems */

  if (!option("scalar")) {
    if ((dt=dparam("DELT")) == 0.) {
      if ((dt=dparam("DT")) != 0.) 
	dparam_set("DELT", dt);
      else
	{fputs("ReadParams: no time step specified", stderr);exit(-1);}    
    } else
      dparam_set("DT", dt);

    /* FINTIME overrides NSTEPS */

    if ((val=dparam("FINTIME")) > 0.)
      iparam_set("NSTEPS", n = (int) ((val+.5*dt)/dt));
    else
      dparam_set("FINTIME", val = (n=iparam("NSTEPS"))*dt);
    
    /* IOTIME overrides IOSTEP */

    n   = clamp(iparam("IOSTEP"), 0 , n);
    val = clamp(dparam("IOTIME"), 0., val);
  
    if (val == 0. && n == 0) {
      n   = iparam("NSTEPS");
      val = dparam("FINTIME");
    } else if (val > 0.)
      n   = (int) ((val + .5 * dt) / dt);
    else
      val = n * dt;

    iparam_set("IOSTEP", n);
    dparam_set("IOTIME", val);
  }

  /* check to see if there is a command line input for norder */
  /* set modes to be equal to norder if it is not already set */

  if(n = option("NORDER.REQ")) iparam_set("MODES",n);
  if(!iparam("MODES")        ) iparam_set("MODES",iparam("NORDER"));

  return;
}

static int checklist (char *name, char *list[])
{
  do
    if (strcmp(name,*list) == 0)
      return 1;
  while 
    (*++list);

  return 0;
}

/* The next two functions are optional */

void ReadPscals (FILE *fp)
{
  int  n;
  char buf[BUFSIZ];
  
  fgets(buf, BUFSIZ, fp);
  if (sscanf(buf,"%d", &n) != 1)
    {fputs("ReadPscals: no passive scalar data", stderr);exit(-1);}
  
  while (n--) fgets(buf, BUFSIZ, fp);
  return;
}

void ReadLogics (FILE *fp)
{
  int  n;
  char buf[BUFSIZ];

  fgets(buf, BUFSIZ, fp);
  if (sscanf(buf, "%d", &n) != 1 || !strstr(buf, SEC_LOGICS)) 
    {fputs("ReadLogics: no logical switches found\n", stderr); exit(-1);}

  while (n--) fgets(buf, BUFSIZ, fp);
  return;
}

static void ReadCurve(FILE *fp, Element_List *new_E);

Element_List *ReadMesh (FILE *fp, char *session_name)
{
  int     nel, L, qa, qb,  k;
  char    buf[BUFSIZ];
  Element **new_E;
  register int i;
  double  xscal = dparam("XSCALE");
  double  yscal = dparam("YSCALE");
  double  xmove = dparam("XMOVE");
  double  ymove = dparam("YMOVE");

	char dummy1[BUFSIZ];
	char dummy2[BUFSIZ];
	char dummy3[BUFSIZ];
  int  group_id = 1; //default value is 1
  int  elmt_id;

  /* set up modes and quadrature points */
  
  if(!( L = iparam("MODES")))
    {fputs("ReadMesh: Number of modes not specified\n",stderr);exit(-1);}

  /* note quadrature order reset for variable order runs */
  if(qa = iparam("LQUAD"));
  else qa = L + 1;

  if(qb = iparam("MQUAD"));
  else qb = L;

  if (!findSection (SEC_MESH, buf, fp))
    {fputs("ReadMesh: Section not found\n", stderr);exit(-1);}

  if (sscanf(fgets(buf,BUFSIZ,fp), "%d%d", &nel,&k) != 2)
   {fputs("ReadMesh: unable to get the number of elements\n",stderr);exit(1);}

  if (k != DIM)
    {fputs("ReadMesh: Read file is wrong dimension\n", stderr);exit(1);}


  iparam_set("ELEMENTS", nel);
  
  /* Set up a new element vector */
  QGmax = qa;
  LGmax = L;

  new_E = (Element**) malloc(nel*sizeof(Element*));

  Coord X;
  X.x = dvector(0,3);
  X.y = dvector(0,3);

  int nvs;
  
  /* Read in mesh information */
  for(k = 0; k < nel; k++) {
    
    fgets(buf, BUFSIZ, fp);  /* element header */
    nvs = (strstr(buf,"Qua") || strstr(buf,"qua")) ? 4 : 3;

    sscanf(buf, "%s %d %s %s %d", dummy1,&elmt_id, dummy2, dummy3, &group_id);

    /* -X- Coordinates */
    for (i = 0; i < nvs; i++)
      fscanf (fp, "%lf", X.x+i);

    fgets(buf, BUFSIZ, fp);  /* finish the line */

    /* -Y- Coordinates */
    for (i = 0; i < nvs; i++)
      fscanf (fp, "%lf", X.y+i);

    fgets(buf, BUFSIZ, fp);  /* finish the line */

    if(xscal)  dscal(nvs,xscal,X.x,1);
    if(yscal)  dscal(nvs,yscal,X.y,1);
    if(xmove)  dsadd(nvs,xmove,X.x,1,X.x,1);
    if(ymove)  dsadd(nvs,ymove,X.y,1,X.y,1);

    if(nvs == 4)
      new_E[k] = new Quad(k,'u', L,qa,qa,0, &X);
    else
      new_E[k] = new  Tri(k,'u', L,qa,qb,0, &X);

      new_E[k]->group_id = group_id;
  }

  for(k = 0; k < nel-1; k++)     new_E[k]->next = new_E[k+1];
  new_E[k]->next = (Element*) NULL;
  
  Element_List* E_List =  new Element_List(new_E, nel);			  
  
  if(option("variable"))
    ReadOrderFile (session_name,E_List);
  
  Tri_work();
  Quad_work();

  E_List->nz = 1;
  E_List->nztot = 1;
  E_List->Cat_mem();

  /* Read the curved side information */
  ReadCurve(fp,E_List);
  
  for(k = 0; k < nel; ++k){
    new_E[k]->set_curved_elmt(E_List);
    new_E[k]->set_geofac();
  }

  // Free the  memory of the first level
  free(E_List->base_h);
  free(E_List->base_hj);
  
  int nz = option("NZ");
  Fourier_List *F_List = new Fourier_List();
  
  F_List->flevels = (Element_List**) malloc(nz*sizeof(Element_List*));
  F_List->flevels[0] = E_List;
  
  if(nz > 1)
    for(k = 1; k < nz; ++k)
      F_List->flevels[k] = E_List->gen_aux_field('u');

  F_List->fhead = E_List->fhead;
  F_List->nel   = E_List->nel;
  F_List->flist = E_List->flist;
  
  F_List->nz    = nz;   // set this afterwards to stop infinite recursion 
  F_List->nztot = option("NZTOT");
  F_List->Cat_mem();
  
  return (Element_List*) F_List;
}


static void ReadCurve(FILE *fp, Element_List *new_E){
  register int i,j;
  int      k;
  int    eid,face;
  char   type, *p, buf[BUFSIZ];
  Curve  *curve,*ctmp;
  
  i = 2; while (i--) fgets (buf, BUFSIZ, fp);
  if(strstr(buf,"type")){
    int    ntags;
    char   *tagid;

    sscanf(buf, "%d", &ntags);

    if(ntags){
      ctmp  = (Curve *) calloc(ntags,sizeof(Curve));
      tagid = (char *)  calloc(ntags+1,sizeof(char));
      
      /* read in curved info */
      for(i = 0; i < ntags; ++i){
	fgets (buf, BUFSIZ, fp);

	if(strstr(buf,"Str")||strstr(buf,"str")){  /* straight sided */
	  fgets (buf, BUFSIZ, fp);
	  sscanf(buf,"%c",tagid+i);
	  ctmp[i].type = T_Straight;      
	}
        else if(strstr(buf,"Cir")||strstr(buf,"cir")){ /* circle */
	  fgets  (buf, BUFSIZ, fp);
	  sscanf (buf, "%*lf %*lf %lf %c",&(ctmp[i].info.arc.radius),tagid+i);
	  ctmp[i].type = T_Arc;
	}	  
	else if(strstr(buf,"Fil")||strstr(buf,"fil")){ /* spline fit - file */
	  char file[BUFSIZ];
	  fgets  (buf, BUFSIZ, fp);
	  sscanf (buf, "%s",file);
	  ctmp[i].info.file.name = strdup(file);
	  sscanf (buf, "%*s %c",tagid+i);
	  ctmp[i].type = T_File;
	}	
	else if(strstr(buf,"Sin")||strstr(buf,"sin")){
	  fgets  (buf, BUFSIZ, fp);
	  sscanf (buf, "%lf %lf %lf %lf %c",&(ctmp[i].info.sin.xo),
		  &(ctmp[i].info.sin.yo), &(ctmp[i].info.sin.amp),
		  &(ctmp[i].info.sin.wavelength),tagid+i);
	  ctmp[i].type = T_Sin;
	}
      }
      
      fgets (buf, BUFSIZ, fp);
      sscanf(buf, "%d", &k);
      
      if (k) {
	if(!ntags) error_msg(No curved type information specified);
	
	curve =  (Curve *) calloc(k,sizeof(Curve));
	
	for(i = 0; i < k; ++i){
	  if (fscanf(fp, "%d %d %c", &face, &eid, &type) != 3)
	    error_msg(ReadMesh -- unable to read curved information);
	  
	  new_E->flist[--eid]->curve = curve + i;
	  
	  for(j = 0; j < ntags; ++j)
	    if(tagid[j] == type){
	      memcpy(curve+i,ctmp+j,sizeof(Curve)); 
	      break;
	    }
	  if(j == ntags){ 
	    fprintf(stderr,"can not find tag type `%c\' called in element"
		    " %d face %d\n" ,type,eid+1,face);
	    exit(-1);
	  }
	  
	  curve[i].face    = --face;
	}
      }
      if(ntags){
	free(ctmp); free(tagid);
      }
    }
  }
  else{
  
    sscanf(buf, "%d", &k);
    
    if (k) {
      curve =  (Curve *) calloc(k,sizeof(Curve));
      
      for(i = 0; i < k; ++i){
	fgets(buf, BUFSIZ, fp);
	if (sscanf(buf, "%d%d", &face, &eid) != 2)
	  error_msg(ReadMesh -- unable to read curved information);

	new_E->flist[--eid]->curve = curve + i;	

	curve[i].next = new_E->flist[--eid]->curve;
	new_E->flist[eid]->curve = curve + i;
	curve[i].face    = --face;
	
	p    = buf + strlen(buf); while (isspace(*--p));
	type = toupper(*p);
	
	switch (type) {
	case 'S':                                   /* Straight side */
	  curve[i].type = T_Straight;
	  break;
	case 'C':                                   /* Circular arc */
	  sscanf (buf, "%*d%*d%lf",&(curve[i].info.arc.radius));
	  curve[i].type = T_Arc;
	  break;
	default: 
	  sprintf (buf, "unknown curve type -- %c", type);
	  fputs(buf,stderr);
	  exit(-1);
	  break;
	}
      }
    }
  }
}

static char *findSection (char *name, char *buf, FILE *fp)
{
  char *p;

  while (p = fgets (buf, BUFSIZ, fp))
    if (strstr (p, name))
      break;

  return p;
}


/* ----------------------------------------------------------------------- *
 * ReadICs() - Read Initial Conditions                                     *
 *                                                                         *
 * This function reads in the initial conditions for the velocity fields   *
 * from the .rea file.  The form of the velocity initial conditions is     *
 * determined by one of the following keywords in the line following the   *
 * beginning of the INITIAL CONDITIONS section:                            *
 *                                                                         *
 *     Default            Initial velocity field is at rest                *
 *     Restart            File name to read IC's from is given             *
 *     Given              V(x,y,z) (t = 0) specified                       *
 * ----------------------------------------------------------------------- */
static void Restart (char *name, int nfields, Element_List *U[]);
double zmesh(int);

void ReadICs (FILE *fp, Domain *omega)
{
  register int i;
  int      n, type, nfields;
  char     buf[BUFSIZ];
  Element_List *U [MAXFIELDS];
  const    int nel=omega->U->nel;
  static  char *ictypes[] = { "PrepOnly", "Restart", "Default", "Given", 
			      "Hardwired", "Minion"};
  
  /* initialize the Element pointers */
  
  U[0] = omega->U; U[1] = omega->V; U[2] = omega->W; 
//  U[3] = omega->Uf;
  U[3] = omega->P;
  nfields = DIM+2;
#ifdef THERMO
  U[4] = omega->T;
  nfields++;
#endif

//#if !defined(THERMO) && defined(MAP)
//  U[4] = omega->Wlast;
//  nfields++;
//#endif

//#if defined(THERMO) && defined(MAP)
// U[5] = omega->Wlast;
//  nfields++;
//#endif

  omega->Uf->fhead->type = 'p';

//  nfields = DIM+1;
  
  /* zero pressure fields */
  dzero(omega->P->nz*omega->P->hjtot,  omega->P->base_hj, 1);
  dzero(omega->P->nz*omega->P->htot,  omega->P->base_h, 1);
  
  if (!findSection (SEC_ICS, buf, fp))
    {fprintf(stderr,"ReadICs: no initial conditions specified "
	     "(last line read below) %s", buf); exit(-1);}
  if (sscanf(buf, "%d", &n) != 1)
    {fprintf(stderr,"ReadICs: # of initial conditions not specified "
	     "(last line read below) %s", buf);exit(-1);}
  
  if (option("prep"))              /* Preprocessor */
    type = 0;
  else if (n == 0)                 /* Default      */
    type = 2;
  else {                           /* Find it...   */
    n--; fgets (buf, BUFSIZ, fp);
    for (i = 1; i < 6; i++) 
      if (strstr(buf, ictypes[type = i])) break;    
  }
  
  switch (type) {
  case 0: 
    break;                      /* Pre-processor only */
  case 1: {                           /* Restart from a file */
#ifdef FLOK // skip past ISLICE lines 
    for (i = 0; i < iparam("ISLICE"); ++i,n--)
      fgets   (buf, BUFSIZ, fp);
#endif

    fgets   (buf, BUFSIZ, fp); n--;
    // restore previous pressure as well
    Restart (buf, nfields, U);


    break;
  }
  case 2: {                           /* Default Initial Conditions */
    for (i = 0; i < nfields; i++){
      U[i]->zerofield();
      U[i]->Set_state('t');
    }
    ROOT printf("Initial solution   : Default\n");
    break;
  }
  case 3: {                          /* Function(x,y,z) Initial conditions */
    char    *s;
    
    ROOT printf("Initial solution   :");
    //U[3] = omega->P;
    
    for (i = 0; i < nfields && n--; i++) {
      if(!(s = strchr(fgets(buf,BUFSIZ,fp),'=')))
	{fprintf(stderr,"ReadICs: the following function is invalid %s",s);
	 exit(-1);}
      while (isspace(*++s));
#ifdef MAP
      Set_mapped_field(U[i], s, omega->mapx, omega->mapy);
#else
      U[i]->Set_field(s);
#endif
      ROOT 
	if(!i) printf (" %c = %s", U[i]->fhead->type, s);
	else   printf ("                     %c = %s", U[i]->fhead->type, s);
    }
#ifdef FLOK
    //    if(omega->U->nz == 2){
    //      printf("Zeroing w component initial conditions\n");
    //      U[2]->zerofield();      
    //    }
#endif

    option_set("GivenICs",1);
    break;
  }
  case 5: {

    // Not set up for Fourier
    int j;
    double x, y, u, v,
      eps   = dparam("SHEAR_EPS"),
      delta = dparam("SHEAR_DELTA");
    Element* E;
    Coord X; 
    
    X.x = dvector(0,QGmax*QGmax-1);
    X.y = dvector(0,QGmax*QGmax-1);
    
    for(i = 0; i < nel; ++i){
      E=U[0]->flist[i];
      E->coord(&X);
      for(j = 0; j < E->qa*E->qb; ++j){
	x = X.x[j];
	y = X.y[j];

	if (y < 0.5)  u = tanh (eps*(y-0.25));
	else          u = tanh (eps*(0.75-y));
	v = delta * cos (2*M_PI*x);

	(*U[0]->flist[i]->h)[j] = u;
	(*U[1]->flist[i]->h)[j] = v;
      }
    }
    set_state(U[0]->fhead,'p');
    set_state(U[1]->fhead,'p');
    ROOT printf("Initial solution   : Minion\n");
    break;
  }

  default:
    fputs("ReadICs: invalid initial conditions", stderr);
    exit(-1);
    break;
  }
  
  omega->Uf->fhead->type = 'u';

  while (n--&&n>0) fgets(buf, BUFSIZ, fp); /* read to the end of the section */
  
  return;
}

int data_len(int nel,int *size, Element *E){
  int ntot = 0;
  int l,i;
  
  for(i=0;i<nel;E=E->next, ++i){
    ntot += E->data_len(size);
    l = E->Nedges + E->Nfaces;
    if(E->dim() == 3) ++l;
    size += l;
  }
  
  return ntot;
}
#if 0
static void Restart (char *name, int nfields, Element_List *U[]){
  Field    f;
  FILE    *fp;
  char    *p, fname[FILENAME_MAX];
  register int i, pos, nel;
  int      nz = U[0]->nz,k;
  int      procid = option("PROCID");
  
  if (sscanf(name, "%s", fname) != 1)
    {fputs("ReadICs: couldn't read the restart file name", stderr);exit(-1);}

  if ((fp = fopen(fname,"r")) == (FILE *) NULL)
    {fprintf(stderr,"ReadICs: restart file %s not found",fname); exit(-1);}

  nel   = U[0]->nel;
  memset(&f, '\0', sizeof (Field));
  //  int dump  = 0;
  //  while (readFieldF (fp, &f, U[0])) dump++;
  //  if (!dump)
  //    {fputs("Restart: no dumps read from restart file", stderr);exit(-1);}
  if (!(readFieldF (fp, &f, U[0])))
    {fputs("Restart: no dumps read from restart file", stderr);exit(-1);}
  if (f.dim != DIM)
    {fputs("Restart: file if wrong dimension", stderr); exit(-1);}

  if (f.nel != nel)
    {fputs("Restart: file is the wrong size", stderr); exit(-1);}
  
  ROOT printf("Initial solution   : %s at t = %g", fname, f.time);
  putc ('\n', stdout);

  if(procid*nz < f.nz){
    f.nz = min(f.nz - procid*nz,nz);
#ifdef FULL_FIELD
    // ce107 - use strlen so as not to have problems with 2d restarts    
    int fnfields = (int) strlen(f.type);
    int ntot = data_len(f.nel,f.size,U[0]->fhead);
    // ce107 - avoid the copy to oneself.
    if (procid)
      for(i = 0; i < fnfields; ++i) /* move relevant data to beginning */
	dcopy(f.nz*ntot,f.data[i]+ntot*procid*nz,1,f.data[i],1);
#endif
  }
  else
    f.nz = 0;

  for (i = 0; i < nfields; i++) {
    if (!(p = strchr(f.type, U[i]->fhead->type))) {
      ROOT fprintf (stderr, "Restart: variable %c is not in file; field set to zero", 
		    U[i]->fhead->type);
      for(k = 0; k < nz; ++k)
	U[i]->flevels[k]->zerofield();
      U[i]->Set_state('t');
    }
    else{
      pos = (int) (p - f.type);
      /* copy all planes into field and zero any others*/
      for(k = 0; k < nz; ++k)
	U[i]->flevels[k]->zerofield();
      U[i]->Set_state('t');

      copyfieldF(&f,pos,U[i]); 
    }
  }
  
  dparam_set("STARTIME", f.time);
  dparam_set("FINTIME", dparam("FINTIME") + f.time);

  freeField (&f);
  fclose    (fp);
  return;
}
#endif
#if 1
static void Restart (char *name, int nfields, Element_List *U[]){
  Field    f;
  FILE    *fp;
  char    *p, fname[FILENAME_MAX];
  register int i, pos, nel;
  int      nz = U[0]->nz,k;
  int      procid = option("PROCID");
  
  if (sscanf(name, "%s", fname) != 1)
    {fputs("ReadICs: couldn't read the restart file name", stderr);exit(-1);}

  if ((fp = fopen(fname,"r")) == (FILE *) NULL)
    {fprintf(stderr,"ReadICs: restart file %s not found",fname); exit(-1);}

  if(1||option("oldhybrid")) // at present need this to read in old files
    set_nfacet_list(U[0]);

  nel   = U[0]->nel;
  memset(&f, '\0', sizeof (Field));

//  if (!(readField (fp, &f)))
  if (!(readFieldF (fp, &f, U[0])))
    {fputs("Restart: no dumps read from restart file", stderr);exit(-1);}
  if (f.dim != DIM)
    {fputs("Restart: file if wrong dimension", stderr); exit(-1);}

  if (f.nel != nel)
    {fputs("Restart: file is the wrong size", stderr); exit(-1);}
  
  ROOT printf("Initial solution   : %s at t = %g", fname, f.time);
  putc ('\n', stdout);

#if 1 
  if(procid*nz < f.nz){
    f.nz = min(f.nz - procid*nz,nz);
#ifdef FULL_FIELD
    // ce107 - use strlen so as not to have problems with 2d restarts    
    int fnfields = (int) strlen(f.type);
    int ntot = data_len(f.nel,f.size,U[0]->fhead);
    // ce107 - avoid the copy to oneself.
    if (procid)
      for(i = 0; i < fnfields; ++i) /* move relevant data to beginning */
	dcopy(f.nz*ntot,f.data[i]+ntot*procid*nz,1,f.data[i],1);
#endif
  }
  else
    f.nz = 0;
#else

  if(procid*nz < f.nz){
    int fnfields = (int) strlen(f.type);
    int ntot = data_len(f.nel,f.size,f.nfacet,f.dim);

    f.nz = min(f.nz - procid*nz,nz);

    if (procid)
      for(i = 0; i < fnfields; ++i) /* move relevant data to beginning */
	dcopy(f.nz*ntot,f.data[i]+ntot*procid*nz,1,f.data[i],1);
  }
  else
    f.nz = 0;


#endif

  for (i = 0; i < nfields; i++) {
    if (!(p = strchr(f.type, U[i]->fhead->type))) {
      ROOT fprintf (stderr, "Restart: variable %c is not in file; "
		    "field set to zero\n", U[i]->fhead->type);
      for(k = 0; k < nz; ++k)
	U[i]->flevels[k]->zerofield();
      U[i]->Set_state('t');
    }
    else{
      pos = (int) (p - f.type);
      /* copy all planes into field and zero any others*/
      for(k = 0; k < nz; ++k)
	U[i]->flevels[k]->zerofield();
      U[i]->Set_state('t');

      copyfieldF(&f,pos,U[i]); 
    }
  }
  
  dparam_set("STARTIME", f.time);
  dparam_set("FINTIME", dparam("FINTIME") + f.time);
  freeField (&f);
  fclose    (fp);
  return;
}
#endif

/* ----------------------------------------------------------------------- *
 * ReadDF() - Read Drive Force data                                        *
 *                                                                         *
 * This is really one function designed to handle two different types of   *
 * drive force specification.  If the "steady" flag is active, the drive   *
 * force can be an arbitary function.  Otherwise, it must be constant but  *
 * can have any direction.                                                 *
 *                                                                         *
 * If the flowrate option is active, this section is skipped entirely.     *
 *                                                                         * 
 * Example:                                                                *
 *                                                                         *
 *   ***** DRIVE FORCE DATA ***** PRESSURE GRAD, FLOW, Q                   *
 * 4                   Lines of Drive force data follow                    *
 *      FFX = 0.                                                           *
 *      FFY = 0.                                                           *
 *      FFZ = 2. * KINVIS                                                  *
 * C                                                                       *
 *                                                                         *
 * The names given to these lines does not matter, but their order is      *
 * taken as (x,y,z).  The '=' MUST be present, and everything to the right *
 * of it determines the forcing functions.                                 *
 * ----------------------------------------------------------------------- */

void ReadDF(FILE *fp, int nforces, ...)
{
  register int i;
  int  nlines;
  char buf[BUFSIZ], *p;
#ifndef THERMO
  static char *ftypes[] = { "FFX", "FFY", "FFZ" };
#else
  static char *ftypes[] = { "FFX", "FFY", "FFZ", "FFT" };
#endif
  Element *E;


  if (!findSection (SEC_DF, buf, fp))
    {fputs("ReadDF: section not found", stderr); exit(-1);}
  if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &nlines) != 1)
    {fputs("ReadDF: can't read the number of lines", stderr); exit(-1);}

  /* Don't process this section if running in pre-processor mode */

  if (option("prep")) {
    while (nlines--) fgets(buf, BUFSIZ, fp);
    return;
  }
  
  /* Check to see if the problem is steady scalar or vector */

  if (option("scalar")) {
    va_list ap;
    Element_List *U [MAXFIELDS];
    
    va_start (ap, nforces);
    for (i = 0; i < nforces; i++) U[i] = va_arg(ap, Element_List *);
    va_end (ap);

    fputs  ("Forcing function    : ", stdout);
    switch (nlines ? iparam("EQTYPE") : Laplace) {
    case Poisson:  case Helmholtz:  default:
      for (i = 0; i < nforces && nlines--; i++) {
	if (p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
	  while (isspace(*++p));
	  if (i) fputs ("                     ", stdout); 
	  printf ("%c   = %s", U[i]->fhead->type, p);
	  U[i]->Set_field(p);
	}
	else
	  fprintf (stderr,"ReadDF -- invalid function (below) %s", buf);
      }
      break;
    case Laplace:
      puts ("none");
      for  (i = 0; i < nforces; i++){
	for(E=U[i]->fhead;E;E=E->next)
	  dzero(E->qa * E->qb , E->h[0], 1);
      }
      break;
    }
  } else {    /* ..... Force specification for unsteady problems ..... */
    if (dparam("FLOWATE") != 0.)       /* check for flowrate */
      option_set ("flowrate", 1);    
    else if (nlines) {                             /* read the list... */
      for (i = 0; (i < nforces) && nlines--; i++) {
	for(p = fgets(buf, BUFSIZ, fp); isspace(*p); p++);
	if (toupper(*p) == 'C') {
	  /* do nothing */
	} else {
	  if (p = strchr(buf,'='))
	    dparam_set(ftypes[i], scalar(++p));
	  else
	    fprintf(stderr,"ReadDF -- invalid function (below) %s", buf);
	}
      }
    }
    ROOT
      if (option("flowrate"))
	printf("Flowrate           : %g\n", dparam("FLOWRATE"));
      else {
	printf("Drive Force        : FFx = %g\n", dparam("FFX"));
	printf("                     FFy = %g\n", dparam("FFY"));
	printf("                     FFz = %g\n", dparam("FFZ"));
#ifdef THERMO
	printf("                     FFt = %g\n", dparam("FFT"));
#endif
      }
  }
  
  while (nlines--&&nlines>0) fgets(buf, BUFSIZ, fp);    /* finish the section */
  return;
}

void ReadDFunc(FILE *fp, Domain *Omega){

  int  nlines, i;
  char buf[BUFSIZ], *p;
  double alpha = dparam("FCESCAL");
  
  if (!findSection (SEC_DFUNCS, buf, fp))
    {
      ROOTONLY fputs("ReadDFunc: section not found\n", stderr);
      Omega->ForceFuncs = (double**)0;
      rewind(fp);
      return;
    }
  if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &nlines) != 1)
    {fputs("ReadDFunc: can't read the number of lines\n", stderr); exit(0);}
  
  int htot=Omega->U->htot*Omega->U->nz;
  int nfields = 3;
#ifdef THERMO
  nfields++;
#endif
  Omega->ForceFuncs = 
    dmatrix(0, nfields-1, 0, htot-1);
          
  Omega->ForceStrings=(char **)malloc(sizeof(char*)*(3));
  for(i=0;i<nfields;i++)
     Omega->ForceStrings[i] =(char *) malloc(sizeof(char)*(1025));

 if (p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
    while (isspace(*++p));
    strcpy(Omega->ForceStrings[0], p);
    Omega->Uf->Set_field(p);
    dcopy(htot, Omega->Uf->base_h, 1, Omega->ForceFuncs[0],1);
  }

  if (p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
    while (isspace(*++p));
    strcpy(Omega->ForceStrings[1], p);
    Omega->Uf->Set_field(p);
    dcopy(htot, Omega->Uf->base_h, 1, Omega->ForceFuncs[1],1);
  }

 if (p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
    while (isspace(*++p));
    strcpy(Omega->ForceStrings[2], p);
    Omega->Uf->Set_field(p);
    dcopy(htot, Omega->Uf->base_h, 1, Omega->ForceFuncs[2],1);
 }
#ifdef THERMO
 if (p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
    while (isspace(*++p));
    strcpy(Omega->ForceStrings[3], p);
    Omega->Tf->Set_field(p);
    dcopy(htot, Omega->Tf->base_h, 1, Omega->ForceFuncs[3],1);
 }
#endif

 return;

}


/* ------------------------------------------------------------------------ *
 * ReadBCs() - Read Boundary Conditions                                     *
 *                                                                          *
 * This function reads boundary information for a spectral element field    *
 * and sets up the appropriate boundary structures.  The following boundary *
 * conditions are currently recognized:                                     *
 *                                                                          *
 *   Flag      Type         Field #1       Field #2       Field #3          *
 *  -----    --------      ----------     ----------     ----------         *
 *    V      Velocity      U-velocity     V-Velocity     W-Velocity         *
 *    W      Wall          = 0            = 0            = 0                *
 *    v      v             U(x,y,z)       V(x,y,z)       W(x,y,z)           *
 *    F      Flux          U' = f1        V' = f2        W' = f3            *
 *    f      flux          U'(x,y,z)      V'(x,y,z)      W'(x,y,z)          *
  *    O      Outflow       U' = 0         V' = 0         W' = 0             *
 *    E      Element       w/ ID          w/ EDGE                           *
 *    P      Periodic      w/ ID          w/ EDGE                           *
 *                                                                          *
 * NOTES:                                                                   *
 *                                                                          *
 *   - On the first call, this function records the beginning of the BC     *
 *     section.  Subsequent calls will re-read the BC's and assign data     *
 *     from the next available parameter.  There should be a total of DIM   *
 *     parameters available (DIM <= 3).                                     *
 *                                                                          *
 *   - For boundaries with a given function of (x,y,z),  the first line     *
  *     following is taken to be U(x,y,z), the second is V(x,y,z), etc.      *
 *     Here is an example (element 1, edge 1):                              *
 *                                                                          *
 *            1  1  v   0.0     0.0     0.0                                 *
 *                 ux = (1 - y^2)*sin(t)                                    *
 *                 uy = (1 - y^2)*cos(t)                                    *
 *                 uz = 0.                                                  *
 *            1  2  E   2.0     1.0     0.0                                 *
 *                                                                          *
 *     The "=" MUST be present, and the function is defined by the string   *
 *     to the right of it.                                                  *
 *                                                                          *
 *   - NO TIME DEPENDENT BOUNDARY CONDITIONS.                               *
 *                                                                          *
 *   - This routine will read any type of boundary conditions, FLUID or     *
 *     otherwise, as long as the file is positioned correctly.              *
 *                                                                          *
 * ------------------------------------------------------------------------ */

#define on  1
#define off 0

#ifdef THERMO
static fpos_t bc_start_temp, bc_stop_temp;
static void   prescan_temp (FILE *fp),
              posscan_temp (FILE *fp);

/* this routine reads all the global mesh boundary conditions to
 determine the type of boundaries and generates a list of dirichlet
 boundary conditions . No extra memory is assign to set the boundary
 conditions */
Bndry **ReadBCs_Temp (FILE *fp, Element_List *EL);

Bndry **ReadBCs_Temp (FILE *fp, Element_List *EL)
{
  register int   n,i,j;
  char       bc, buf[BUFSIZ];
  int        iel, iside, nbcs=0;

  Element    *E;
  Bndry      **bndry_list, **new_bndry_list, 
             *new_bndry  = (Bndry *) NULL;
  double     f[3];
  static int fnum_temp = 0;         /* this tracks the field numbers read */
  char       *temp_ptr;
  double     hc, T_infty, sigma;

  if (!fnum_temp) prescan_temp(fp);

  bndry_list = (Bndry**) calloc(EL->nz,sizeof(Bndry*));
  new_bndry_list = (Bndry**) calloc(EL->nz,sizeof(Bndry*));

  //clear all vertex mode
  for(E = EL->fhead; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      E->vert[i].solve = 1;

  /* .. Begin Reading Boundary Conditions .. */
  for(i = 0; i < EL->nz; ++i){
    rewind(fp);
    fsetpos(fp, &bc_start_temp);
    dparam_set("z", zmesh(i));

    for(E=EL->flevels[i]->fhead;E;E=E->next){ 
      for (n = 0; n < E->Nedges; ++n){
	    char *p = fgets (buf, BUFSIZ, fp);
	/* vertex solve mask is on by default */
	
	while (isspace(*p)) p++; bc = *p++;
	sscanf(p, "%d%d%lf%lf%lf", &iel, &iside, f, f+1, f+2); 
	
	if (iel - E->id != 1 || iside - n != 1) {
	  sprintf(buf, "Mismatched element/side -- got %d,%d, expected %d,%d",
		  iel, iside, E->id+1, n+1);
	}
	
	// Boundaries set by value
	
	switch (bc) {
	case 'W':
	  f[fnum_temp] = 0.;
	case 'V':
	  E->set_solve(n,off);
	  new_bndry = E->gen_bndry(bc, n, f[fnum_temp]);
	  break;

	case 'Z': // this is only set up for reflection in x-axis
	  f[fnum_temp] = 0.;
	  if(E->type == 'v'){
	    E->set_solve(n,off);
	    new_bndry = E->gen_bndry('W', n, 0.0);
	  }
	  else
	     new_bndry = E->gen_bndry('F', n, 0.0 );
	  
	  
	  new_bndry->usrtype = bc;
	  option_set("REFLECT", 1);
	  break;
	case 'I': case 'O':
	   E->set_solve(n,on);

	  f[fnum_temp] = 0.;
	case 'F': 
	  new_bndry = E->gen_bndry(bc, n, f[fnum_temp]);
	  break;
  case 'R': // Robin Boundary condition
	  fgets(buf, BUFSIZ, fp); // line of "h = ..." convective heat transfer coefficient
	  temp_ptr = strchr(buf, '=');
	  if(!temp_ptr) {
	    fprintf(stderr, "ERROR: ReadBCs_Temperature: Robin BC, hc line format error, "
		    "element %d, face %d\n", iel, iside);
	    exit(-1);
	  }
	  temp_ptr++;
	  hc = scalar(temp_ptr); // get h value

	  fgets(buf, BUFSIZ, fp); // line of "T_infty = ..." Wll temperature
	  temp_ptr = strchr(buf, '=');
	  if(!temp_ptr) {
	    fprintf(stderr, "ERROR: ReadBCs_Temperature: Robin BC, Tz line format error, "
		    "element %d, face %d\n", iel, iside);
	    exit(-1);
	  }
	  temp_ptr++;
	  T_infty = scalar(temp_ptr);


	  fgets(buf, BUFSIZ, fp); // line of "sigma = ..." Robin coefficient
	  temp_ptr = strchr(buf, '=');
	  if(!temp_ptr) {
	    fprintf(stderr, "ERROR: ReadBCs_Temperature: Robin BC, Tz line format error, "
		    "element %d, face %d\n", iel, iside);
	    exit(-1);
	  }
	  temp_ptr++;
	  sigma = scalar(temp_ptr);

	  f[fnum_temp] = 0;
	  new_bndry = E->gen_bndry(bc, n, f[fnum_temp]);
	  new_bndry->aux_data[0] = hc;
	  new_bndry->aux_data[1] = T_infty; 
    new_bndry->sigma       = sigma;

	  E->Surface_geofac(new_bndry);

    break;
	  
	  /* Functional boundaries */
	  
	case 't': case 'v':
	  E->set_solve(n,off);
	case 'f': case 'r': {
	  fpos_t pos;
	
	  for (j = 0; j <= fnum_temp; j++) fgets(buf, BUFSIZ, fp);
	  new_bndry = E->gen_bndry(bc, n, buf);
	  
	  do {                         
	    fgetpos (fp, &pos);         /* Search through the following    */
	    fgets   (buf, BUFSIZ, fp);  /* lines for the first one without */
	  } while (strchr (buf, '='));  /* an '=' (function specification) */
	  
	  fsetpos (fp, &pos);
	  break;
	}
	  
	  /* Element and Periodic boundaries */
	  
	case 'E': case 'P':
	  new_bndry = (Bndry *) NULL;
	  break;
	default:
	  fputs("ReadBCs: read a strange boundary condition", stderr);
	  exit(-1);
	  break;
	}
	
	if (new_bndry) {
	  nbcs++;
	  new_bndry->next = bndry_list[i];
	  bndry_list[i]   = new_bndry;
	}
      }
 #if DIM == 2
      if(E->Nedges!=4)
	fgets (buf, BUFSIZ, fp); /* skip extra side */
#endif
    }
  }
  
  if (fnum_temp++ == 0) posscan_temp(fp); fsetpos(fp, &bc_stop_temp);
  
  for(i = 0; i < EL->nz; ++i)
    new_bndry_list[i] = bsort(bndry_list[i],nbcs/EL->nz);
      
  return new_bndry_list;      /* Return the sorted list          */
}

static void prescan_temp(FILE *fp)
{
  char buf[BUFSIZ];

  rewind(fp);
  findSection(SEC_T_BCS,buf,fp);
  if (!strstr(buf, SEC_T_BCS))
    error_msg(ReadBCs_Temp: do not see any thermal boundary conditions);

//  do
//    fgets(buf, BUFSIZ, fp);
//  while
//    (strstr(buf, "NO"));
  
  fgetpos (fp, &bc_start_temp);
  return;
}

static void posscan_temp (FILE *fp)
{
  char buf[BUFSIZ];
  fgetpos (fp, &bc_stop_temp);
  while (strstr(fgets(buf,BUFSIZ,fp), SEC_T_BCS))
    fgetpos (fp, &bc_stop_temp);
  return;
}

#endif

static fpos_t bc_start, bc_stop;
static void   prescan (FILE *fp),
              posscan (FILE *fp);

#ifdef MAP
Bndry **ReadBCs (FILE *fp, Element_List *EL, double *x, double *y)
#else
Bndry **ReadBCs (FILE *fp, Element_List *EL)
#endif
{
  register int   n,i,j;
  char       bc, buf[BUFSIZ];
  int        iel, iside, nbcs=0;

   Element    *E;
  Bndry      **bndry_list, **new_bndry_list, 
             *new_bndry  = (Bndry *) NULL;
  double     f[3];
  static int fnum = 0;         /* this tracks the field numbers read */

  if (!fnum) prescan(fp);

  bndry_list = (Bndry**) calloc(EL->nz,sizeof(Bndry*));
  new_bndry_list = (Bndry**) calloc(EL->nz,sizeof(Bndry*));
#ifdef MAP
  // invFFT the maps to physical space
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) x, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) y, 1, 0, 0, 0, 0);

  int nztot=option("NZTOT");
  double *fD;
  rewind(fp);
  if (findSection("INFLOW PROFILE", buf, fp)){
    int mpnts=nztot;
    fD=new double[mpnts];
    dzero(mpnts,fD,1);
    for(i=0; i<nztot; i++)
      if(sscanf(fgets(buf, BUFSIZ, fp), "%lf", fD+i) != 1){
	{ROOT fputs("ReadBCs: can't read fft value", stderr);
	exit(-1);}
      }
  }
#endif

  /* .. Begin Reading Boundary Conditions .. */
  for(i = 0; i < EL->nz; ++i){
    rewind(fp);
    fsetpos(fp, &bc_start);
    dparam_set("z", zmesh(i));
#ifdef MAP
    int iloc = parid(i);
     dparam_set("ORIGINX", x[iloc]);
    dparam_set("ORIGINY", y[iloc]);
#endif
    for(E=EL->flevels[i]->fhead;E;E=E->next){ 
      for (n = 0; n < E->Nedges; ++n){
	char *p = fgets (buf, BUFSIZ, fp);
	
	/* vertex solve mask is on by default */
	
	while (isspace(*p)) p++; bc = *p++;
	sscanf(p, "%d%d%lf%lf%lf", &iel, &iside, f, f+1, f+2); 
	
	if (iel - E->id != 1 || iside - n != 1) {
	  sprintf(buf, "Mismatched element/side -- got %d,%d, expected %d,%d",
		  iel, iside, E->id+1, n+1);
	  fprintf(stderr,"ReadBCs : %s\n", buf); exit(-1);
	}
	
	// Boundaries set by value
	
	switch (bc) {
#ifdef MAP
	case 'D':
	  E->set_solve(n,off);
	  if(!fnum)
	  	new_bndry = E->gen_bndry('V', n, fD[iloc]);
	  else
	  	new_bndry = E->gen_bndry('V', n, (double)0);
	  break;
#endif
	case 'W':
	  f[fnum] = 0.;
	case 'V':
	  E->set_solve(n,off);
	  new_bndry = E->gen_bndry(bc, n, f[fnum]);
 //zwang 05042020 
	  new_bndry->usrtype = 'W';
	  break;
	case 'S': //slip velocity boundary condition
	  E->set_solve(n,off);
	  new_bndry = E->gen_bndry(bc, n, f[fnum]);
	  break;

	case 'Z': // this is only set up for reflection in x-axis
	  f[fnum] = 0.;
	  if(E->type == 'v'){
	    E->set_solve(n,off);
	    new_bndry = E->gen_bndry('W', n, 0.0);
	  }
	  else
	     new_bndry = E->gen_bndry('F', n, 0.0 );
	  
	  
	  new_bndry->usrtype = bc;
	  option_set("REFLECT", 1);
	  break;
	case 'I': case 'O': 
	  f[fnum] = 0.;
	case 'F': case 'R':
	  new_bndry = E->gen_bndry(bc, n, f[fnum]);
	  break;
	  
	  /* Functional boundaries */
	  
  case 'v':
	  E->set_solve(n,off);
	case 'f': case 'r': {
	  fpos_t pos;
	
	  for (j = 0; j <= fnum; j++) fgets(buf, BUFSIZ, fp);
	  new_bndry = E->gen_bndry(bc, n, buf);
	  
	  do {                         
	    fgetpos (fp, &pos);         /* Search through the following    */
	    fgets   (buf, BUFSIZ, fp);  /* lines for the first one without */
	  } while (strchr (buf, '='));  /* an '=' (function specification) */
	  
	  fsetpos (fp, &pos);
	  break;
	}
	  
	case 't': 
  {
//zwang 05042020
    bc = 'v';
	  E->set_solve(n,off);
//    fprintf(stderr,"I get here \n");
	  fpos_t pos;
	
	  for (j = 0; j <= fnum; j++) fgets(buf, BUFSIZ, fp);
	  new_bndry = E->gen_bndry(bc, n, buf);
	  
	  do {                         
	    fgetpos (fp, &pos);         /* Search through the following    */
	    fgets   (buf, BUFSIZ, fp);  /* lines for the first one without */
	  } while (strchr (buf, '='));  /* an '=' (function specification) */
	  
	  fsetpos (fp, &pos);
 //zwang 05042020 
	  new_bndry->usrtype = 'T';
	  break;
	}

	case 'x': 
  {
//zwang 05042020
    bc = 'v';
	  E->set_solve(n,off);
//    fprintf(stderr,"I get here \n");
	  fpos_t pos;
	
	  for (j = 0; j <= fnum; j++) fgets(buf, BUFSIZ, fp);
	  new_bndry = E->gen_bndry(bc, n, buf);
	  
	  do {                         
	    fgetpos (fp, &pos);         /* Search through the following    */
	    fgets   (buf, BUFSIZ, fp);  /* lines for the first one without */
	  } while (strchr (buf, '='));  /* an '=' (function specification) */
	  
	  fsetpos (fp, &pos);
 //zwang 05042020 
	  new_bndry->usrtype = 'X';
	  break;
	}
	  /* Element and Periodic boundaries */
	  
	case 'E': case 'P':
	  new_bndry = (Bndry *) NULL;
	  break;
	default:
	  fputs("ReadBCs: read a strange boundary condition", stderr);
	  exit(-1);
	  break;
	}
	
	if (new_bndry) {
	  nbcs++;
	  new_bndry->next = bndry_list[i];
	  bndry_list[i]   = new_bndry;
	}
      }
 #if DIM == 2
      if(E->Nedges!=4)
	fgets (buf, BUFSIZ, fp); /* skip extra side */
#endif
    }
  }
#ifdef MAP
  // FFT the maps to fourier space
  rfftw(rplan, 1, (FFTW_COMPLEX *) x, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) y, 1, 0, 0, 0, 0);
#endif
  
  if (fnum++ == 0) posscan(fp); fsetpos(fp, &bc_stop);
  
  for(i = 0; i < EL->nz; ++i)
    new_bndry_list[i] = bsort(bndry_list[i],nbcs/EL->nz);
      
  return new_bndry_list;      /* Return the sorted list          */
}

/* ------------------------------------------------------------------------ *
 * bsort() - Sort Boundary Conditions                                       *
 *                                                                          *
 * The following function sorts boundary conditions so they are set accord- *
 * ing to precedence and not mesh order.  The following is the "hierarchy"  *
 * of boundary conditions :                                                 *
 *                                                                          *
 *    W          : Walls (zero velocity) are the most binding               *
 *    V          : Velocity boundary conditions (constant)                  *
 *    v          : Velocity function boundary conditions                    *
 *    all others                                                            *
 *                                                                          *
 * This function returns the sordid boundary list (ha ha).                  *
 *                                                                          *
 * ------------------------------------------------------------------------ */
Bndry *bsort(Bndry *bndry_list, int nbcs){
  Bndry   **tmp  = (Bndry**) malloc( nbcs * sizeof(Bndry*) );
  Bndry   *bedg  = bndry_list;
  register int i;

  if(!nbcs) return bndry_list;

  /* Copy the linked list into a regular array and sort it */

  for(i = 0 ; i < nbcs; i++, bedg = bedg->next) tmp[i] = bedg;
  qsort(tmp, nbcs, sizeof(Bndry*), bedgcmp);

  /* Create a new version of the ordered list */
  
  bedg  = (Bndry*) calloc (nbcs, sizeof(Bndry));
  for(i = 0; i < nbcs; i++) {
    memcpy (bedg + i, tmp[i], sizeof(Bndry));
    bedg[i].id         = i;
    bedg[i].next = bedg + i + 1;
  }
  bedg[nbcs-1].next   = (Bndry*) NULL;
  
  /* Free up the space occupied by the old list */

  bndry_list = bedg;
#if 0
  for(i = 0; i < nbcs; i++) free (tmp[i]);
  free (tmp);
#endif
  return bndry_list;
}

/*
 * Boundaries are sorted by their ASCII character values except for
 * types 'D' and 'N', which are sorted by their element ID numbers.
 */

static int bedgcmp(const void *b1, const void *b2)
{
  Bndry  *bedg_1 = *((Bndry**) b1),
         *bedg_2 = *((Bndry**) b2);
  char  btype_1 = bedg_1 -> type,
        btype_2 = bedg_2 -> type;

  /* Convert {D,N} -> {G,H} to get the ASCII precedence right */

  if (btype_1 == 'D') btype_1 = 'G';
  if (btype_2 == 'D') btype_2 = 'G';
  if (btype_1 == 'N') btype_1 = 'H';
  if (btype_2 == 'N') btype_2 = 'H';
  
  if      (btype_1 > btype_2)              /* Check ASCII code */
    return  1;
  else if (btype_1 < btype_2)
    return -1;
  else                                     /* Check element ID */
    {
      int id_1 = bedg_1 -> elmt -> id,
          id_2 = bedg_2 -> elmt -> id;

      if       (id_1 > id_2)
	return  1;
      else if  (id_1 < id_2)
	return -1;
      else                                   /* Check edge ID    */
	{
#if DIM == 3
	  id_1 = bedg_1 -> face;
	  id_2 = bedg_2 -> face;
#else
	  id_1 = bedg_1 -> face;
	  id_2 = bedg_2 -> face;
#endif
	  if      (id_1 > id_2)
	    return  1;
	  else if (id_1 < id_2)
	    return -1;
	}
    }

  /* If we get this far there is an error in the boundary conditions */
  
  fputs("bedgcmp(): " "Duplicate boundary\n",stderr); exit(-1);
  return 0;
}

static void prescan(FILE *fp)
{
  char buf[BUFSIZ];

  rewind(fp);
  findSection(SEC_BCS,buf,fp);
  if (!strstr(buf, SEC_BCS))
    error_msg(ReadBCs: do not see any boundary conditions);

  do
    fgets(buf, BUFSIZ, fp);
  while
    (strstr(buf, "NO"));
  
  fgetpos (fp, &bc_start);
  return;
}

static void posscan (FILE *fp)
{
  char buf[BUFSIZ];
  fgetpos (fp, &bc_stop);
  while (strstr(fgets(buf,BUFSIZ,fp), SEC_BCS))
    fgetpos (fp, &bc_stop);
  return;
}

/* set up inter element links using the coundary connectivity information */

static void set_link (Element *U, int eid1, int id1, int eid2, int id2);

void ReadSetLink(FILE *fp, Element_List *UL)
{
  register int   n;
  char       bc, buf[BUFSIZ];
  int        iel, iside;
  Element    *E;
  double     f[3];
  Element    *U = UL->fhead;

  rewind(fp);
  findSection(SEC_BCS,buf,fp);
  if (!strstr(buf, SEC_BCS))
    error_msg(ReadSetLink: do not see any boundary conditions);
  fgets(buf,BUFSIZ,fp);
  
  /* ....... Search through Boundary Conditions ...... */

  char *p = fgets(buf,BUFSIZ, fp);

  for(E=U;E;E=E->next){
    for (n = 0; n < E->Nedges; ++n){

      while (isspace(*p)) p++; bc = *p++;
      sscanf(p,"%d%d%lf%lf%lf", &iel, &iside, f, f+1, f+2); 
 
      if (iel - E->id != 1 || iside - n != 1) {
	sprintf(buf, "Mismatched element/side -- got %d,%d, expected %d,%d",
		iel, iside, E->id+1, n+1);
	fprintf(stderr,"ReadSetLink : %s\n", buf); exit(-1);
      }

      if ((bc == 'E')||(bc == 'P')) /* Element and Periodic boundaries */
	set_link(U,iel-1,iside-1,(int)f[0]-1,(int)f[1]-1);
//      if ((bc == 'v')||(bc == 'f')||(bc == 'r'))
//      zwang 05042020
      if ((bc == 'v')||(bc == 'f')||(bc == 'r')||(bc == 't') || (bc == 'x'))
	while(strchr(p=fgets(buf,BUFSIZ,fp),'='));
      else
	p=fgets(buf,BUFSIZ,fp);
    }
#if DIM == 2
    if(E->Nverts==3){
      p=fgets(buf,BUFSIZ,fp);
    }
#endif
  }
  return;
}

int  con[4][4] = {{1,1,0,0},{1,1,0,0},{0,0,1,1},{0,0,1,1}};


static void set_link(Element *U, int eid1, int id1, int eid2, int id2){
  /* set up connectivity inverses if required */
  /* like faces meet */
  Element *E=U,*F=U;
  while(E->id !=eid1)
    E=E->next;
  while(F->id !=eid2)
    F=F->next;
  
  if(eid1>eid2)
    E->edge[id1].con = con[id1][id2];
  else
    F->edge[id2].con = con[id1][id2];
  
  /* set links */
  if(eid1 < eid2){
    E->edge[id1].base = E->edge + id1;
    E->edge[id1].link = F->edge + id2;  
    F->edge[id2].base = E->edge + id1;
  }
}

void  ReadOrderFile(char *name, Element_List *E){
  register int i,k;
  int *size,n;
  const int nel = E->nel;
  char buf[BUFSIZ],*p;
  FILE *fp;
  
  sprintf(buf,"%s.ord",name);
  if(!(fp = fopen(buf,"r"))){
    fprintf(stderr,"ERROR: file \"%s\" does not exist\n",buf); 
    exit(-1);
  }

#if DIM == 2
  size = ivector(0,5);
#else
  size = ivector(0,10);
#endif
  
  /* skip comments */
  while(p=fgets(buf,BUFSIZ,fp))
    if(strstr(p,"Modes")){
      for(k=0;k<nel;++k){
	fscanf(fp,"%*d");
#if DIM == 2
	for(i = 0; i < E->flist[k]->Nedges; ++i)
	  fscanf(fp,"%d",size+i);
	fscanf(fp,"%d",size+i);
#else
	for(i = 0; i < 11; ++i)
	  fscanf(fp,"%d",size+i);
#endif
	E->flist[k]->Mem_J(size,'n');  
	E->flist[k]->Mem_Q();
      }
      break;
    }
    else if (strstr(p,"Function")){
#if DIM == 2
      double x[5],y[5],f[5];
#else
      double x[11],y[11],z[11],f[11];
      extern int ednum[][3];
#endif
      p = fgets(buf,BUFSIZ,fp);
      
      if((p = (strchr(p,'='))) == (char *) NULL)
	error_msg(CheckOrderFIle: Illegal function definition);
#if DIM == 2
      vector_def("x y",++p);
#else
      vector_def("x y z",++p);
#endif
      int max_f;
      Element *F;

      for(k=0;k<nel;++k){
	F=E->flist[k];
#if DIM == 2
	for(i = 0; i < F->Nverts; ++i){
	  x[i] = 0.5*(F->vert[i].x + F->vert[(i+1)%F->Nverts].x);
	  y[i] = 0.5*(F->vert[i].y + F->vert[(i+1)%F->Nverts].y);
	}
	vector_set(F->Nverts,x,y,f);
	/* set fourth side to be consistent with highest edge */
	max_f=0;
	for(i = 0; i < F->Nedges; ++i)
	  max_f = (f[i]>max_f) ? (int)f[i]:max_f;
	f[F->Nverts] = max_f-1;
	
	if(n=option("aposterr")) 
	  for(i = 0; i < F->Nedges+1; ++i) size[i] = (int)max(f[i],0.0) + n;
	else
	  for(i = 0; i < F->Nedges+1; ++i) size[i] = (int)max(f[i],0.0);
#else
	for(i = 0; i < 3; ++i){
	  x[i] = 0.5*(E[k].vert[i].x + E[k].vert[(i+1)%3].x);
	  y[i] = 0.5*(E[k].vert[i].y + E[k].vert[(i+1)%3].y);
	  z[i] = 0.5*(E[k].vert[i].z + E[k].vert[(i+1)%3].z);
	}
	for(i = 0; i < 3; ++i){
	  x[i+3] = 0.5*(E[k].vert[i].x + E[k].vert[3].x);
	  y[i+3] = 0.5*(E[k].vert[i].y + E[k].vert[3].y);
	  z[i+3] = 0.5*(E[k].vert[i].z + E[k].vert[3].z);
	}

	vector_set(6,x,y,z,f);
	/* set faces to be maximun of surrounding edges */
	/* and edges to be maximun of faces             */
	f[10] = 0;
	for(i = 0; i < 4; ++i){
	  f[i+6] = max(max(f[ednum[i][0]],f[ednum[i][1]]),f[ednum[i][2]])-1;
	  f[10]  = max(f[10],f[i+6]);
	}
	--f[10];

	for(i = 0; i < 11; ++i) size[i] = (int)max(f[i],0.0);
#endif
	F->Mem_J(size,'n');  
	F->Mem_Q();
      }
      break;
    }
  option_set("variable",1);

  free(size);  
}


/* ----------------------------------------------------------------------- *
 * ReadHisData() - Read History and Integral data                          *
 *                                                                         *
 * This function reads in and sets up structures for the processing of     *
 * history points.  This is comprised of a linked list of structures which *
 * are processed by the Analyzer().                                        *
 *                                                                         *
 * ----------------------------------------------------------------------- */

static HisPoint *appendHisPoint (HisPoint *list, HisPoint *hp);

void ReadHisData (FILE *fp, Domain *omega)
{
  char      buf[BUFSIZ], *p;
  int       cnt, npts, nel,nz = option("NZTOT");
  HisPoint *hp, *his_list = (HisPoint*)NULL;
  
  if (!findSection(SEC_HIST, buf, fp))
    {ROOT fputs("ReadHisData: Section not found", stderr); exit(-1);}
  if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &npts) != 1)
    {ROOT fputs("ReadHisData: can't read the number of points", stderr);
     exit(-1);}

  nel    = iparam("ELEMENTS");

  while (npts--) {
    hp = (HisPoint*) calloc(1, sizeof(HisPoint));
    for (p = fgets(buf, BUFSIZ, fp); !isdigit(*p); p++)
      if (isalpha(*p=tolower(*p)))
	switch(*p){
	case 'h': /* treat as vertex point */
	  hp->mode = TransVert;
	  break;
	case 'e': /* treat as mid point of edge */
	  hp->mode = TransEdge;
	  break;
	default:
	  strncat(hp->flags, p, 1);
	}
    
    hp->i = (int) strtol (p, &p, 10);
    hp->j = (int) strtol (p, &p, 10);
    hp->k = (int) strtol (p, &p, 10);

    if (hp->k >= nz)
      error_msg(Bad history points: nz < k);

    if(hp->k) --hp->k; /* allow points with k = 0 to remain be first points */

    hp->id = atoi (p);
    
    if ((--hp->i) > omega->U->flist[hp->id-1]->Nedges || (hp->id--) > nel ) 
      { free(hp); continue; }
    
    his_list = appendHisPoint (his_list, hp);
  }


  /* This history step can be overridden by setting HISSTEP */

    option_set ("hisstep", iparam ("HISSTEP") > 0 ?
		iparam ("HISSTEP") : max(1, iparam("NSTEPS") / 1000));

    omega->his_list = his_list;

    ROOT if (omega->his_list) {
      /* Echo the history points */

      puts("\nHistory points:");
      for (hp = his_list, cnt = 0; hp; hp = hp->next) {
	if(hp->mode == TransVert)
	  printf ("%2d: vertice = %2d, elmt = %3d, k = %3d, fields = %s\n", 
		  ++cnt, hp->i + 1, hp->id + 1, hp->k + 1, hp->flags);
	else if(hp->mode == TransEdge)
	  printf ("%2d: edge    = %2d, elmt = %3d, k = %3d, fields = %s\n", 
		  ++cnt, hp->i + 1, hp->id + 1, hp->k + 1, hp->flags);
      }
    }
  
  putchar ('\n');

  return;
}

static HisPoint *appendHisPoint (HisPoint *list, HisPoint *hp)
{
  HisPoint *h = list;

  if (h) {
    while (h->next) h = h->next;
    h->next = hp;
  } else
    list = hp;

  return list;
}

void ReadLambda (FILE *fp, Domain *omega)
{
  char    buf[BUFSIZ];
  double *lambdas;
  int     i,nlams;
 
  if (!findSection(SEC_LAMB, buf, fp)) {
    fputs("ReadLambda: Section not found", stderr); 
  }
  
  fgets (buf, BUFSIZ, fp);
  sscanf(buf, "%d", &nlams);

  lambdas = dvector(0,nlams-1);

  for (i = 0; i < nlams; i++)
      fscanf(fp,"%lf",lambdas+i);

  omega->lambdas = lambdas;

  return;
}

void ReadRhoCp (FILE *fp, Domain *omega)
{
  char    buf[BUFSIZ];
  double *RhoCp;
  int     i,nrc;
 
  if (!findSection(SEC_RHO_CP, buf, fp)) {
    fputs("ReadRhoCp: Section not found", stderr); 
  }
  
  fgets (buf, BUFSIZ, fp);
  sscanf(buf, "%d", &nrc);

  RhoCp = dvector(0,nrc-1);

  for (i = 0; i < nrc; i++)
      fscanf(fp,"%lf",RhoCp+i);

  omega->RhoCp = RhoCp;

  return;
}

//REMI 
#if 0
void ReadTension (FILE *fp, Domain *omega)
{
  char     buf[BUFSIZ];
  int      i, nztot = option("NZTOT");
  double *Tension = (double *) calloc (nztot, sizeof(double));
 
//  if (option("dealias")) nztot = 3*nztot/2;
 ROOTONLY fprintf(stderr,"nztot=%d\n",nztot);
  if (!findSection(SEC_TENS, buf, fp)) {
    fputs("ReadTension: Section not found", stderr); 
  }
  Tension = dvector(0,nztot-1);
  for (i = 0; i < nztot; i++)
      fscanf(fp,"%lf",Tension+i);

  omega->tens = Tension;
 // for (i = 0; i < nztot; i++)
//  ROOTONLY fprintf(fp,"%lf",Tension+i);
  return;
}

void ReadStiffne (FILE *fp, Domain *omega)
{
  char     buf[BUFSIZ];
  int      i, nztot = option("NZTOT");
  double *Stiffness = (double *) calloc (nztot, sizeof(double));

//  if (option("dealias")) nztot = 3*nztot/2;

  if (!findSection(SEC_STIF, buf, fp)) {
    fputs("ReadStiffne: Section not found", stderr); 
  }
  Stiffness = dvector(0,nztot-1);
  for (i = 0; i < nztot; i++)
      fscanf(fp,"%lf",Stiffness+i);

  omega->stif = Stiffness;
 // for (i = 0; i < nztot; i++)
  // ROOTONLY fprintf(fp,"%lf",Stiffness+i);
  return;
}
#else

void ReadTension (Domain *omega)
{
  char     buf[BUFSIZ];
  int      i, nztot = option("NZTOT");
  int      Nplanes;
  double *Tension = (double *) calloc (nztot, sizeof(double));
 
  char Fname[FILENAME_MAX];
  sprintf(Fname, "%s_TensStiff.geo", omega->name); // stiffness and tension file
  FILE *fp;
  fp = fopen(Fname,"r");
  if (fp == NULL){
    fprintf(stderr,"ReadTension: file %s is not found\n",Fname);
    exit(-1);
  }
//  if (option("dealias")) nztot = 3*nztot/2;
// ROOTONLY fprintf(stderr,"nztot=%d\n",nztot);
  if (!findSection(SEC_TENS, buf, fp)) {
    fputs("ReadTension: Section not found", stderr); 
  }
  Tension = dvector(0,nztot-1);

  fscanf(fp,"%d",&Nplanes);

  if(Nplanes != nztot)
   {
    fprintf(stderr,"error in ReadTension: nztot != Nplanes \n");
    exit(-1);
   }

  for (i = 0; i < nztot; i++)
      fscanf(fp,"%lf",Tension+i);

  omega->tens = Tension;
 // for (i = 0; i < nztot; i++)
//  ROOTONLY fprintf(fp,"%lf",Tension+i);
  return;
}

//void ReadStiffne (FILE *fp, Domain *omega)
void ReadStiffne (Domain *omega)
{
  char     buf[BUFSIZ];
  int      i, nztot = option("NZTOT");
  int      Nplanes;
  double *Stiffness = (double *) calloc (nztot, sizeof(double));

//  if (option("dealias")) nztot = 3*nztot/2;
  char Fname[FILENAME_MAX];
  sprintf(Fname, "%s_TensStiff.geo", omega->name); // stiffness and tension file
  FILE *fp;
  fp = fopen(Fname,"r");
  if (fp == NULL){
    fprintf(stderr,"ReadStiffne: file %s is not found\n",Fname);
    exit(-1);
  }

  if (!findSection(SEC_STIF, buf, fp)) {
    fputs("ReadStiffne: Section not found", stderr); 
  }
  Stiffness = dvector(0,nztot-1);

  fscanf(fp,"%d",&Nplanes);

  if(Nplanes != nztot)
   {
    fprintf(stderr,"error in ReadStiffne: nztot != Nplanes \n");
    exit(-1);
   }
  for (i = 0; i < nztot; i++)
      fscanf(fp,"%lf",Stiffness+i);

  omega->stif = Stiffness;
 // for (i = 0; i < nztot; i++)
  // ROOTONLY fprintf(fp,"%lf",Stiffness+i);
  fclose(fp); 
  return;
}

#endif

#ifdef FLOW_CONTROL
//Hilbert Coefficients for P control
void ReadHilbertCoeffs (Domain *omega)
{
  char     buf[BUFSIZ];
  int      i, ncoeff = iparam("NHF_COEFF");
  double   *HF_Coeff = (double *) calloc (ncoeff, sizeof(double));
 
  char Fname[FILENAME_MAX];
  sprintf(Fname, "%s_hilbert.coeff", omega->name); // hilbert coeffs file
  FILE *fp;
  fp = fopen(Fname,"r");
  if (fp == NULL){
    fprintf(stderr,"ReadHilbertCoeffs: file %s is not found\n",Fname);
    exit(-1);
  }

  for (i = 0; i < ncoeff; i++)
      fscanf(fp,"%lf",HF_Coeff+i);

  omega->hilbert_coeff = HF_Coeff;

  return;
}
#endif

void ReadIntData (FILE *fp, Domain *omega)
{
  char     buf[BUFSIZ];
  int      i, npts, nfields, nztot = option("NZTOT");
  Intepts *I;
  
  nfields = DIM + 2;
  I = (Intepts *)malloc(sizeof(Intepts));
 
  if (!findSection(SEC_INT, buf, fp)) {
    fputs("ReadIntData: Section not found", stderr); 
    npts = 0;
  } else
    if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &npts) != 1) {
      fputs("ReadIntData: can't read the number of points", stderr);
      exit(-1);
    }
 
  I->npts = npts;
  if (npts > 0){ 
    I->X.x  = dvector(0,I->npts-1);
    I->X.y  = dvector(0,I->npts-1);
    I->X.z  = dvector(0,I->npts-1);
    I->ui   = dmatrix(0,I->npts-1,0,nfields-1);
    I->ufr  = dmatrix(0,I->npts-1,0,nfields*nztot-1);
    for (i = 0; i < I->npts; i++)
      fscanf(fp,"%lf%lf%lf",I->X.x+i,I->X.y+i,I->X.z+i); 
    omega->int_list = I;
  } else
    omega->int_list = (Intepts *) NULL;
 
  return;
}

void ReadSoln(FILE* fp, Domain* omega){
  char buf[BUFSIZ], *s;
  int  i;
  int nfields = DIM+2;

  rewind(fp);
  omega->soln = (char**) malloc(nfields*sizeof(char*));

  if(findSection("Exact", buf, fp)){
    for(i=0;i<nfields;++i){
      if(s = strchr(fgets(buf,BUFSIZ,fp),'=')){
	while (isspace(*++s));
	omega->soln[i] = (char*) malloc(BUFSIZ*sizeof(char));
	strcpy(omega->soln[i],s);
      }
      else {
	omega->soln[i] = (char*) malloc(BUFSIZ*sizeof(char));
	sprintf(omega->soln[i],"0.0");
      }
    }
    ROOT{
      fprintf(stderr,"Exact U : %s",omega->soln[0]);
      fprintf(stderr,"Exact V : %s",omega->soln[1]);
      fprintf(stderr,"Exact W : %s",omega->soln[2]);
      fprintf(stderr,"Exact P : %s",omega->soln[3]);
    }
  }
  else
    for(i=0;i<nfields;++i){
      omega->soln[i] = (char*) malloc(BUFSIZ*sizeof(char));
      sprintf(omega->soln[i],"0.0");
    }

  rewind(fp);
}

// This is for Spectral Viscosity  SV
void SetEpsilon(FILE* fp, Element_List *U, Bndry *Ubc, Metric *lambda){
  int  i,k,skip;
  char buf[BUFSIZ], *s;
  Element *E;

  option_set("FAMOFF", 1);
  for (k=0;k<U->nel;k++)
    lambda[k].epsilon = dvector(0, QGmax*QGmax);
  printf("Sorry..... this does not work!!!!!! \n");
  rewind(fp);
  if(findSection("Epsilon", buf, fp)){
    if(s = strchr(fgets(buf,BUFSIZ,fp),'=')){
      while (isspace(*++s));
      fprintf(stdout, "Epsilon set by function: %s \n", s);
      
      U->Set_field(s);
      dcopy(U->htot, U->base_h, 1, lambda->epsilon, 1);
      //dsmul(U->htot,dparam("EPSILON"),lambda->epsilon, 1,lambda->epsilon, 1);
      for(skip=0, k=0;k<U->nel;skip+=E->qa*E->qb,++k){
	E = U->flist[k];
	lambda[k].epsilon = lambda->epsilon + skip;
      }
    }
  }
  else
    {
      for (k=0;k<U->nel;k++)
	{
	  E = U->flist[k];
	  for (i=0;i<E->qa*E->qb;i++)
	    lambda[k].epsilon[i] = dparam("EPSILON"); 
	}
    }

  rewind(fp);
}

void ReadWave(FILE *fp, double **wave, Element_List *U, Domain *omega){
  int  i;
  int eDIM = U->fhead->dim();

  eDIM++;

//#ifdef BLASIUS    // MSB: 1 for current, 2 for d/dz, 3 for initial f(\eta)
  eDIM = 3*eDIM+1;  // MSB: +1 for storage
//#endif
  char buf[BUFSIZ], *p;
  
  if (!findSection (SEC_WAVE, buf, fp)){
    rewind(fp);
    if(!findSection (SEC_WAVE, buf, fp)){
      fputs("ReadWave: section not found\n", stderr);
      return;
    }
  }

  int      n, type, nfields;
  //  char     buf[BUFSIZ];
  Element_List *Uw [MAXFIELDS];
  const    int nel=omega->U->nel;
  static  char *ictypes[] = { "PrepOnly", "Restart", "Default", "Given", 
			      "Hardwired", "Minion"};
  
  /* initialize the Element pointers */
  
  Uw[0] = omega->U; Uw[1] = omega->V; Uw[2] = omega->W; Uw[3] = omega->P;
  nfields = DIM+1;

  /* zero pressure fields */

  dzero(omega->P->nz*omega->P->hjtot,  omega->P->base_hj, 1);
  dzero(omega->P->nz*omega->P->htot,  omega->P->base_h, 1);

  if (!findSection (SEC_ICS, buf, fp))
    {fprintf(stderr,"ReadWave: no initial conditions specified "
	     "(last line read below) %s", buf); exit(-1);}
  if (sscanf(buf, "%d", &n) != 1)
    {fprintf(stderr,"ReadWave: # of initial conditions not specified "
	     "(last line read below) %s", buf);exit(-1);}

  if (option("prep"))              /* Preprocessor */
    type = 0;
  else if (n == 0)                 /* Default      */
    type = 2;
  else {                           /* Find it...   */
    n--; fgets (buf, BUFSIZ, fp);
    for (i = 1; i < 6; i++) 
      if (strstr(buf, ictypes[type = i])) break;    
  }
  
  switch (type) {
  case 0: 
    break;                            /* Pre-processor only  */
  case 1: {                           /* Restart from a file */
    //#ifdef FLOK // skip past ISLICE lines 
    //    for (i = 0; i < iparam("ISLICE"); ++i,n--)
    //      fgets   (buf, BUFSIZ, fp);
    //#endif
    
    fgets   (buf, BUFSIZ, fp); n--;
    // restore previous pressure as well
    ROOT printf("Base Flow   : ");
    Restart (buf, DIM+2, Uw);
    for (i = 0; i < nfields; i++) Uw[i]->Trans(Uw[i],J_to_Q);
    break;
  }
  case 2: {                           /* Default Initial Conditions */
    for (i = 0; i < nfields; i++){
      Uw[i]->zerofield();
      Uw[i]->Set_state('t');
    }
    ROOT printf("Base Flow   : Default\n");
    break;
  }
  case 3: {                          /* Function(x,y,z) Initial conditions */
    char    *s;
    
    ROOT printf("Base Flow   : Given");
    
    for (i = 0; i < nfields && n--; i++) {
      if(!(s = strchr(fgets(buf,BUFSIZ,fp),'=')))
	{fprintf(stderr,"ReadWave: the following function is invalid %s",s);
	 exit(-1);}
      while (isspace(*++s));

      Uw[i]->Set_field(s);
      ROOT 
	if(!i) printf ("  %c = %s", Uw[i]->fhead->type, s);
	else   printf ("                     %c = %s", Uw[i]->fhead->type, s);
    }
    
    option_set("GivenICs",1);
    break;
  }
    return;
  }

  for(int j = 0; j < nfields; j++){
    for(i = 0; i < U->nel; ++i)
      dcopy(Uw[j]->flist[i]->qtot,Uw[j]->flist[i]->h[0],1,wave[i*eDIM+j],1);
  }

  
  /*
  
  if(p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
    while (isspace(*++p));
    U->Set_field(p);
    for(i = 0; i < U->nel; ++i)
      dcopy(U->flist[i]->qtot,U->flist[i]->h[0],1,wave[i*eDIM],1);
  }
  
  if(p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
    while (isspace(*++p));
    U->Set_field(p);
    for(i = 0; i < U->nel; ++i)
      dcopy(U->flist[i]->qtot,U->flist[i]->h[0],1,wave[i*eDIM+1],1);
  }

  if(p = strchr (fgets(buf, BUFSIZ, fp), '=')) {
    while (isspace(*++p));
    U->Set_field(p);
    for(i = 0; i < U->nel; ++i)
      dcopy(U->flist[i]->qtot,U->flist[i]->h[0],1,wave[i*eDIM+2],1);
  }
  */

}

#ifdef MAP
void ReadMStatPoints(FILE* fp, Domain* omega){
  char buf[BUFSIZ];
  int  i, npts;
  
  MStatStr *mstat = (MStatStr *) calloc(1,sizeof(MStatStr));

  rewind(fp);

  if (findSection("Moving Statistics", buf, fp)) {
    if (sscanf(fgets(buf, BUFSIZ, fp), "%d", &npts) != 1)
    {fputs("ReadMStatPoints: can't read the number of points", stderr);
     exit(-1);}
 
    mstat->n = npts;
    if (npts > 0) {
      mstat->x     = dvector(0,npts-1);
      mstat->y     = dvector(0,npts-1);
      mstat->sgnx  = ivector(0,npts-1);
      mstat->sgny  = ivector(0,npts-1);
      mstat->nstat = (int *) calloc(npts,sizeof(int));
      for (i = 0; i < npts; i++)
	fscanf(fp,"%lf%lf%d%d",mstat->x+i, mstat->y+i, 
	       mstat->sgnx+i, mstat->sgny+i); 
    }
  }

  omega->mstat = mstat;
  rewind(fp);
  return;
}
#endif

