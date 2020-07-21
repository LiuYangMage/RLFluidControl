/**************************************************************************/
//                                                                        //
//   Author:    S.Sherwin                                                 //
//   Design:    T.Warburton && S.Sherwin                                  //
//   Date  :    12/4/96                                                   //
//                                                                        //
//   Copyright notice:  This code shall not be replicated or used without //
//                      the permission of the author.                     //
//                                                                        //
/**************************************************************************/

#include <stdio.h>
#include <veclib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "hotel.h"

#ifndef BTYPE
#if defined(i860) || defined (__alpha) || defined (__WIN32__) || (defined(linux) && defined(i386))
#define BTYPE "ieee_little_endian"
#endif
#
#if defined(_CRAY) && !defined (_CRAYMPP)
#define BTYPE "cray"
#endif /* ........... Cray Y-MP ........... */
#
#ifndef BTYPE
#define BTYPE "ieee_big_endian"
#endif /* default case in the absence of any other TYPE */
#endif /* ifndef TYPE */

#define DESCRIP 25
static char *hdr_fmt[] = { 
  "%-25s "            "Session\n",
  "%-25s "            "Created\n",
  "%-5c Hybrid              " "State 'p' = physical, 't' transformed\n",
  "%-7d %-7d %-7d   " "Number of Elements; Dim of run; Lmax\n",
  "%-25d "            "Step\n",
  "%-25.6g "          "Time\n",
  "%-25.6g "          "Time step\n",
  "%-25.6g "          "Kinvis;\n",
  "%-25s "            "Fields Written\n",
  "%-25s "            "Format\n"
  };

static char *hdr_fmt_comp[] = { 
  "%-25s "            "Session\n",
  "%-25s "            "Created\n",
  "%-5c HyCompress          " "State 'p' = physical, 't' transformed\n",
  "%-7d %-7d %-7d   " "Number of Elements; Dim of run; Lmax\n",
  "%-25d "            "Step\n",
  "%-25.6g "          "Time\n",
  "%-25.6g "          "Time step\n",
  "%-25.6g "          "Kinvis;\n",
  "%-25s "            "Fields Written\n",
  "%-25s "            "Format\n"
  };


#if 0
static char *fourier_hdr_fmt[] = { 
  "%-25s "            "Session\n",
  "%-25s "            "Created\n",
  "%-5c Fourier Hybrid      " "State 'p' = physical, 't' transformed\n",
  "%-5d %-5d %-5d %-5d   "    "Number of Elements; Dim of run; Lmax; NZ\n",
  "%-25d "            "Step\n",
  "%-25.6g "          "Time\n",
  "%-25.6g "          "Time step\n",
  "%-11.6g %-11.6g   " "Kinvis; LZ\n",
  "%-25s "            "Fields Written\n",
  "%-25s "            "Format\n"
  };
#endif

#define BINARY   strings[0]	   /* 0. Binary format file            */ 
#define ASCII    strings[1]	   /* 1. ASCII format file             */ 
#define KINVIS   strings[2]	   /* 2. Kinematic viscosity           */ 
#define DELT     strings[3]	   /* 3. Time step                     */ 
#define CHECKPT  strings[4]        /* 4. Check-pointing option         */
#define BINFMT   strings[5]        /* 5. Binary format string          */

static char *strings [] = { 
  "binary", "ascii", "KINVIS", "DELT", "checkpt",
  "binary-"BTYPE
};

static int data_len(int nel,int *size, Element *E);
int gettypes    (char *t, char *s);
static int checkfmt    (char *format);

int  count_facets(Element *U){
  Element *E;
  int cnt = 0;  // interior mode length

  for(E=U;E; E= E->next)
    cnt += E->Nedges+E->Nfaces;

  if(U->dim() == 3)
    for(E=U;E; E= E->next)
      ++cnt;

  return cnt;
}

int  count_facets(Element_List *U, int *emap, int nel){
  Element *E;
  int i, cnt = 0;  // interior mode length
  
  for(i=0;i<nel;++i){
    E = U->flist[emap[i]];
    
    cnt += E->Nedges+E->Nfaces;

    if(E->dim() == 3)
      ++cnt;
  }

  return cnt;
}



void Writefield(FILE *fp, char *name, int step, double t,
		 int nfields, Element_List *U[]){

  register int n,i,j;
  time_t   tp;
  char     buf[128],state;
  Field    f;
  int      ntot;
  double   **u;

  //  fprintf(stderr,"max= %d; nfi= %d\n",MAXFIELDS,nfields);

  if(nfields > MAXFIELDS)
    error_msg(too many fields -- must be less than MAXFIELDS);

  state = U[0]->fhead->state;
  /* check to see if all fields in same state */
  for(i = 1; i < nfields; ++i)
    if(state != U[i]->fhead->state)
      {error_msg(Writefield--Fields  in different space);}
      
  /* Get the current date and time from the system */
  
  tp = time((time_t*)NULL);
  strftime(buf, DESCRIP, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  
  /* Set up the field structure */
  f.name      =  name;
  f.created   =  buf;
  f.state     =  state;
  f.dim       =  U[0]->fhead->dim();
  f.nel       =  U[0]->nel;
  f.lmax      =  LGmax;
  f.step      =  step;
  f.time      =  t;
  f.time_step =  dparam("DT");
  f.kinvis    =  (dparam("KINVIS-ORIG"))?
    dparam("KINVIS-ORIG") : dparam("KINVIS");   /* VARV */ 
  f.format    = (option(BINARY)) ? BINFMT : ASCII;

  /* Generate the list of fields and assign data pointers */
  
  memset (f.type, '\0', MAXFIELDS);
  memset (f.data, '\0', MAXFIELDS * sizeof( double *));
  
  int cnt = count_facets(U[0]->fhead);

  /* copyfield into matrix */
  if(state == 't'){
    f.size = ivector(0,cnt-1);

    for(i = 0,n=0; i < f.nel; ++i){
      for(j = 0; j < U[0]->flist[i]->Nedges; ++j,++n)
	f.size[n] = U[0]->flist[i]->edge[j].l;
      
      for(j = 0; j < U[0]->flist[i]->Nfaces; ++j,++n)
	f.size[n] = U[0]->flist[i]->face[j].l;
      
      if(U[0]->fhead->dim() == 3){
	f.size[n] = U[0]->flist[i]->interior_l;
	++n;
      }
    }
    
    ntot = U[0]->hjtot;
    u = dmatrix(0,nfields-1,0,ntot-1);    
    dzero(nfields*ntot, u[0], 1);
    for(n = 0; n < nfields; ++n){
      f.type [n] = U[n]->fhead->type;
      f.data [n] = u[n];
      dcopy(U[0]->hjtot, U[n]->fhead->vert[0].hj, 1, u[n], 1);
    }
  }
  else{
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }

  writeField(fp,&f, U[0]->fhead);
  
  free(f.size);
  free_dmatrix(u,0,0);

  return;
}

/* Transfer a field-type array from s -> t.            *
 *                                                     *
 * NOTE: we assume that t is at least MAXFIELDS long!! */

int gettypes (char *t, char *s)
{
  char    *p = strstr (s, "Fields");
  int n = 0;
  fprintf(stderr,"max= %d \n",MAXFIELDS);

  if (!p){ fprintf(stderr,"invalid \"Fields\" string");exit(-1);}
  memset(t, '\0', MAXFIELDS);
  while (s != p && n < MAXFIELDS) {
    if (isalpha(*s)) t[n++] = *s;
    s++;
  }
  return n;
}


/* Check binary format compatibility */
#if defined(__sgi) || defined(cm5) || defined(_AIX) || defined (__hpux) || defined (__sparc) || defined (_CRAYMPP)
#define ieeeb
#define ieee
#endif

#if defined(i860) || defined (__alpha) || defined (__WIN32__) || (defined(linux) && defined(i386))
#define ieeel
#endif

static int checkfmt (char *arch)
{
  char        **p;
  static char *fmtlist[] = {
#if defined(ieeeb)                  /* ... IEEE big endian machine ..... */
    "ieee_big_endian",
    "sgi", "iris4d", "SGI", "IRIX", /* kept for compatibility purposes   */
    "IRIX64",                       /* ........ Silicon Graphics ....... */
    "AIX",                          /* .......... IBM RS/6000 .......... */
    "cm5",                          /* ........ Connection Machine ..... */
#endif
#
#if defined(ieeel)                  /* ... IEEE little endian machine .. */
    "ieee_little_endian",
    "i860",                         /* ........... Intel i860 .......... */
#endif
#
#if defined(ieee)                   /* ...... Generic IEEE machine ..... */
    "ieee", "sim860",               /* kept for compatibility purposes   */
#endif                              /* same as IEEE big endian           */
#
#if defined(_CRAY) && !defined (_CRAYMPP) /* ...... Cray PVP ........... */
    "cray", "CRAY",                 /* kept for compatibility purposes   */
#endif                              /* no conflict with T3* as it        */
#                                   /* precedes this line                */
     0 };   /* a NULL string pointer to signal the end of the list */

  for (p = fmtlist; *p; p++)
    if (strncmp (arch, *p, strlen(*p)) == 0)
      return 0;

  return 1;
}

#undef ieee

/* ---------------------------------------------------------------------- *
 * writeField() -- Write a field file                                     *
 *                                                                        *
 * This function writes a field file from a Field structure.  The format  *
 * is a simple array, stored in 64-bit format.                            *
 *                                                                        *
 * Return value: number of variables written, -1 for error.               *
 * ---------------------------------------------------------------------- */
  
int writeField (FILE *fp, Field *f, Element *E)
{
  int  nfields, ntot;
  char buf[BUFSIZ];
  register int i, n;
#ifdef _CRAY
  short *ssize, *spllinfo;
#endif
  char **fmt;

  if(option("COMPRESSWRITE"))
    fmt = hdr_fmt_comp;
  else
    fmt = hdr_fmt;
  
  /* Write the header */
  fprintf (fp, fmt[0],  f->name);
  fprintf (fp, fmt[1],  f->created);
  fprintf (fp, fmt[2],  f->state);
#ifdef FOURIER 
  fprintf (fp, fmt[3],  f->nel, f->dim, f->lmax, f->nz);
#else
  fprintf (fp, fmt[3],  f->nel, f->dim, f->lmax);
#endif
  fprintf (fp, fmt[4],  f->step);
  fprintf (fp, fmt[5],  f->time);
  fprintf (fp, fmt[6],  f->time_step);
#ifdef FOURIER
  fprintf (fp, fmt[7],  f->kinvis, f->lz);
#else
  fprintf (fp, fmt[7],  f->kinvis);
#endif
  fprintf (fp, fmt[8],  f->type);
  fprintf (fp, fmt[9],  f->format);

  /* Write the field files  */
  if(f->state == 't'){
#ifdef FOURIER
    ntot = f->nz*data_len(f->nel,f->size,E);
#else
    ntot = data_len(f->nel,f->size,E);
#endif
  }
  else{
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
      
  nfields = (int) strlen (f->type);

  int cnt = count_facets(E);
  switch (tolower(*f->format)) {
  case 'a':
    for(i = 0; i < cnt; ++i)
      fprintf(fp,"%d ",f->size[i]);
    fputc('\n',fp);
    for (i = 0; i < ntot; i++) {
      for (n = 0; n < nfields; n++)
	fprintf (fp, "%#16.10g ", f->data[n][i]);
      fputc ('\n', fp);
    }
    
#ifdef PARALLEL
    for(i = 0; i < pllinfo.nloop; ++i)
      fprintf(fp,"%d ",pllinfo.eloop[i]);
    fputc ('\n', fp);
#endif

    break;
    
  case 'b':
#ifdef _CRAY
    // int on Crays is 64b, short 32b - use short for portable fieldfiles
    ssize = (short *) calloc(cnt, sizeof(short));
    for (n = 0; n < cnt; n++)
      ssize[n] = (short) f->size[n];
    fwrite(ssize, sizeof(short),  cnt, fp);
    free (ssize);
#else
    fwrite(f->size, sizeof(int),  cnt, fp);
#endif
    
    for (n = 0; n < nfields; n++)
      if (fwrite (f->data[n], sizeof(double), ntot, fp) != ntot) {
	fprintf  (stderr, "error writing field %c", f->type[n]);
	exit(-1);
      }
#ifdef PARALLEL
#ifdef _CRAY
    // int on Crays is 64b, short 32b - use short for portable fieldfiles
    spllinfo = (short *) calloc(pllinfo.nloop, sizeof(short));
    for (n = 0; n < pllinfo.nloop; n++)
      spllinfo[n] = (short) pllinfo.eloop[n];
    fwrite(spllinfo, sizeof(short), pllinfo.nloop, fp);
    free (spllinfo);
#else
    fwrite(pllinfo.eloop, sizeof(int),  pllinfo.nloop, fp);
#endif
#endif
    break;
    
  default:
    sprintf  (buf, "unknown format -- %s", f->format);
    fprintf(stderr,buf);
    exit(-1);
    break;
  }

  fflush (fp);
  return nfields;
}

/* ---------------------------------------------------------------------- *
 * readField() -- Load a field file                                       *
 *                                                                        *
 * This function loads a field file into a Field structure.  The format   *
 * is a simple array, stored in 64-bit format.                            *
 *                                                                        *
 * Return value: number of variables loaded, 0 for EOF, -1 for error.     *
 * ---------------------------------------------------------------------- */

#define READLINE(fmt,arg) \
  if (fscanf(fp,fmt,arg)==EOF) return -1; fgets(buf,BUFSIZ,fp)

int readField (FILE *fp, Field *f, Element *E)
{
  register int i, n, ntot, nfields;
  char     buf[BUFSIZ];
#ifdef _CRAY
  short *ssize; // semap
#endif
  if (fscanf (fp, "%s", buf) == EOF) return 0;            /* session name */
  f->name    = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               /* creation date */
  f->created = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);


  READLINE ("%c"  , &f->state);                  /* simulation parameters */

  if(strstr(buf, "HyCompress"))
    option_set("COMPRESSREAD",1);

  
  fscanf(fp,"%d%d%d",&f->nel,&f->dim,&f->lmax);  fgets(buf,BUFSIZ,fp);

  READLINE ("%d"  , &f->step);                  
  READLINE ("%lf" , &f->time);
  READLINE ("%lf" , &f->time_step);
  READLINE ("%lf" , &f->kinvis);

  nfields = gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = (char *)strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("Warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);

  /* allocate memory and load the data */
  
#ifdef PARALLEL
  f->emap = ivector(0,f->nel-1);
#endif
  
  int cnt = count_facets(E);
  switch (tolower(*f->format)) {
  case 'a':
    if(f->state == 't'){
      f->size = ivector(0,cnt-1);
      for(i = 0; i < cnt ; ++i) fscanf(fp,"%d",f->size+i);
    }
    else {
      fputs("Error: reading field file in physical space is "
	    "not implemented \n",stderr);
      exit(-1);
    }
    
    ntot = data_len(f->nel,f->size, E);
    
    for(i = 0; i < nfields; ++i)
      if(!(f->data [i]))
	f->data[i] = dvector (0, ntot-1);

    for (i = 0; i < ntot; i++) {
      for (n = 0; n < nfields; n++)
	if (fscanf (fp, "%lf", f->data[n] + i) != 1) {
	  fprintf(stderr,"Error: reading field %c, point %d\n",f->type[n],i+1);
	  exit(-1);
	}
    }
    fgets (buf, BUFSIZ, fp);
#if 0
#ifdef PARALLEL
    for(i = 0; i < f->nel; ++i)
	if(fscanf(fp,"%d",f->emap + i) != 1){
	  fprintf(stderr,"Error: reading emap, point %d\n",i+1);
	  exit(-1);
	}
#endif
#endif
    break;
    
  case 'b':
    if(f->state == 't'){
      f->size = ivector(0,cnt-1);
#ifdef _CRAY
      // int on Crays is 64b, short 32b - use short for portable fieldfiles
      ssize = (short *) calloc(cnt, sizeof(short));
      fread(ssize, sizeof(short), (size_t) cnt, fp);
      for (n = 0; n < cnt; n++)
	f->size[n] = (int) ssize[n];
      free (ssize);
#else
      fread(f->size, sizeof(int), (size_t) cnt,fp);
#endif
    }
    else {
      fputs("Error: reading field file in physical space is "
	    "not implemented \n",stderr);
      exit(-1);
    }

    ntot = data_len(f->nel,f->size,E);
    
    for(i = 0; i < nfields; ++i)
      if(!(f->data [i]))
	f->data[i] = dvector (0, ntot-1);
    
    for (n = 0; n < nfields; n++) {
      if (fread (f->data[n], sizeof(double), (size_t) ntot, fp) != ntot) {
	fprintf (stderr, "error reading field %c", f->type[n]);
	exit(-1);
      }
    }
#if 0
#ifdef PARALLEL
#ifdef _CRAY
    // int on Crays is 64b, short 32b - use short for portable fieldfiles
    semap = (short *) calloc(f->nel, sizeof(short));
    if(fread(semap, sizeof(short), (size_t) f->nel, fp) != f->nel)
      fprintf (stderr, "error reading emap");
    for (n = 0; n < f->nel; n++)
      f->emap[n] = (int) semap[n];
    free (semap);
#else
    if(fread(f->emap, sizeof(int), (size_t) f->nel,fp) != f->nel)
      fprintf (stderr, "error reading emap");
#endif
#endif
#endif
    break;

  default:
    fprintf (stderr, "unknown field format -- %s\n", f->format);
    exit(-1);
    break;
  }

  return nfields;
}

/* ---------------------------------------------------------------------- *
 * freeField() -- Free the memory associated with a Field                 *
 * ---------------------------------------------------------------------- */

void freeField (Field *f)
{
  int      n = (int) strlen (f->type);
  double **d = f->data;

  if (f->name)    free (f->name);
  if (f->created) free (f->created);
  if (f->format)  free (f->format);
  
  while (*d && n--) free (*d++);
  
  memset (f, '\0', sizeof(Field));
  return;
}

/*--------------------------------------------------------------------------*
 * This is a function that copies the field form the field structure 'f'    *
 * at the position given by pos to the Element 'U'. If in transformed space *
 * it allows the field to be copied at a different order. To do this it     *
 * either zeros the higher modes if U is higher than f else it ignores the  *
 * higher modes if f is higher than U                                       *
 *--------------------------------------------------------------------------*/
void copyfield(Field *f, int pos, Element *U){

  if(f->state == 'p'){ /* restart from fixed m physical field */
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }
  else{ /* restart from transformed field at any order */
    int      l;
    int      *size  = f->size;
    double   *data  = f->data[pos];
    
    while(U){
      U->Copy_field(data, size);
      
      l = U->Nedges + U->Nfaces;
      if(U->dim() == 3)
	++l;
      
      data += U->data_len(size);
      size += l;
      
      U->state = 't';
      U=U->next;
      
    }
  }
}

static int data_len(int nel,int *size, Element *E){
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


int readHeader(FILE* fp, Field *f){
  register int  nfields;
  char     buf[BUFSIZ];
  if (fscanf (fp, "%s", buf) == EOF) return 0;            /* session name */
  f->name    = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               /* creation date */
  f->created = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  
  READLINE ("%c"  , &f->state);                  /* simulation parameters */
#ifdef FOURIER
  /* check to see if it is standard 2d file */
  if(strstr(buf,"Fourier")){
    trip = 1;
    fscanf(fp,"%d%d%d%d",&f->nel,&f->dim,&f->lmax,&f->nz);fgets(buf,BUFSIZ,fp);
  }
  else{
    fscanf(fp,"%d%d%d",&f->nel,&f->dim,&f->lmax);fgets(buf,BUFSIZ,fp);
    f->nz = 1;
  }
#else
  fscanf(fp,"%d%d%d",&f->nel,&f->dim,&f->lmax);  fgets(buf,BUFSIZ,fp);
  f->nz = 1;
#endif
  
  READLINE ("%d"  , &f->step);                  
  READLINE ("%lf" , &f->time);
  READLINE ("%lf" , &f->time_step);
#ifdef FOURIER
  if(trip){
    fscanf(fp,"%lf%lf " ,&f->kinvis,&f->lz);
    fgets(buf,BUFSIZ,fp);
  }
  else{
    READLINE ("%lf" , &f->kinvis);
    f->lz = 1.0;
  }
#else
  READLINE ("%lf" , &f->kinvis);
#endif
  
  nfields = gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = (char *)strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("Warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);
  rewind(fp);
  return nfields;
}
static void writeHeader(FILE *fp, Field *f, int Shuffled, Element_List *EL);
static int  writeData(FILE *fp, Field *f, Element_List *EL);


static char *p_hdr_fmt[] = { 
  "%-25s ",
  "%-25s ",
  "%-5c ",
  "%-5d %-5d %-5d ",
  "%-25d ", 
  "%-25.6g ", 
  "%-25.6g ", 
  "%-11.6g ",
  "%-25s ",
  "%-25s "
};


static char *p_hdr_fmt_exp[] = { 
  "Session\n",
  "Created\n",
  "State 'p' = physical, 't' transformed\n",
  "Number of Elements; Dim of run; Lmax; nz\n",
  "Step\n",
  "Time\n",
  "Time step\n",
  "Kinvis; LZ\n",
  "Fields Written\n",
  "Format\n"
};



void Writefld(FILE **fp, char *name, int step, double t,
		 int nfields, Element_List *UL[]){
  register int n,i,j;
  time_t   tp;
  char     buf[128],state;
  Field    f;
  int      ntot;
  double   **u;

  if(nfields > MAXFIELDS)
    error_msg(too many fields -- must be less than MAXFIELDS);

  state = UL[0]->fhead->state;
  /* check to see if all fields in same state */
  for(i = 1; i < nfields; ++i)
    if(state != UL[i]->fhead->state)
      {error_msg(Writefield--Fields  in different space);}
      
  /* Get the current date and time from the system */
  
  tp = time((time_t*)NULL);
  strftime(buf, DESCRIP, "%a %b %d %H:%M:%S %Y", localtime(&tp));
  
  /* Set up the field structure */
  f.name      =  name;
  f.created   =  buf;
  f.state     =  state;
  f.dim       =  UL[0]->fhead->dim();
  f.nel       =  UL[0]->nel;
  f.lmax      =  LGmax;
  f.step      =  step;
  f.time      =  t;
  f.time_step =  dparam("DT");
  f.kinvis    =  (dparam("KINVIS-ORIG"))?
    dparam("KINVIS-ORIG") : dparam("KINVIS");   /* VARV */
  f.format    = (option(BINARY)) ? BINFMT : ASCII;
  f.lz        =  dparam("LZ");
  f.nz        =  option("NZTOT");

  /* Generate the list of fields and assign data pointers */
  
  memset (f.type, '\0', MAXFIELDS);
  memset (f.data, '\0', MAXFIELDS * sizeof( double *));

  int cnt = count_facets(UL[0]->fhead);
  
  /* copyfield into matrix */
  if(state == 't'){
    f.size = ivector(0,cnt-1);

    for(i = 0,n=0; i < f.nel; ++i){
      for(j = 0; j < UL[0]->flist[i]->Nedges; ++j,++n)
	f.size[n] = UL[0]->flist[i]->edge[j].l;
      
      for(j = 0; j < UL[0]->flist[i]->Nfaces; ++j,++n)
	f.size[n] = UL[0]->flist[i]->face[j].l;
      
      if(UL[0]->fhead->dim() == 3){
	f.size[n] = UL[0]->flist[i]->interior_l;
	++n;
      }
    }
    
    ntot = UL[0]->hjtot;
    u = dmatrix(0,nfields-1,0,ntot-1);    
    dzero(nfields*ntot, u[0], 1);
    for(n = 0; n < nfields; ++n){
      f.type [n] = UL[n]->fhead->type;
      f.data [n] = u[n];
      dcopy(UL[0]->hjtot, UL[n]->fhead->vert[0].hj, 1, u[n], 1);
    }
  }
  else{
    fprintf(stderr,"Not implemented\n");
    exit(-1);
  }

  /* if parallel the header and data are written out seperately */
#ifdef PARALLEL
  DO_PARALLEL
    { 
      if(fp[0]) writeHeader(fp[0],&f,1,UL[0]);
      else fprintf(stderr,"No header file pointer\n");
      if(fp[1]) writeData  (fp[1],&f,UL[0]);
      else fprintf(stderr,"No data file pointer\n");
    }
    else
#endif
      { /* else write it all to one file */
	if(fp[0]){
	  writeHeader(fp[0],&f,0,UL[0]);
	  writeData  (fp[0],&f,UL[0]);
	}
	else
	  fprintf(stderr,"No file pointer \n");
      }
  
  free(f.size);
  free_dmatrix(u,0,0);

  return;
}


/* Write the header and expansion order information */
static void writeHeader(FILE *fp, Field *f, int Shuffled, Element_List *EL){
#ifdef PARALLEL
  char buf[BUFSIZ];
  register int i;
#ifdef _CRAY
  short *ssize, *spllinfo;
  register int n;
#endif
  
  if(f->state == 'p'){
    fprintf(stderr,"Can't write field in physical space\n");
    exit(-1);
  }

  sprintf (buf,"%s%s",p_hdr_fmt[0],p_hdr_fmt_exp[0]);
  fprintf (fp, buf,  f->name);
  sprintf (buf,"%s%s",p_hdr_fmt[1],p_hdr_fmt_exp[1]);
  fprintf (fp, buf,  f->created);

  
  if(Shuffled)
      sprintf (buf,"%s    Shuffled Hybrid       %s",p_hdr_fmt[2],p_hdr_fmt_exp[2]);
  else
    if(f->nz > 1) /* fourier case */
      sprintf (buf,"%s    Fourier         %s",p_hdr_fmt[2],p_hdr_fmt_exp[2]);
    else
      sprintf (buf,"%s                    %s",p_hdr_fmt[2],p_hdr_fmt_exp[2]);
  
  fprintf (fp, buf,  f->state);
  
  if(f->nz > 1){ /* fourier case */
    sprintf (buf,"%s%s%s",p_hdr_fmt[3],"%-5d   ",p_hdr_fmt_exp[3]);
    fprintf (fp, buf,  f->nel, f->dim, f->lmax, f->nz);
  }
  else{
    sprintf (buf,"%s        %s",p_hdr_fmt[3],p_hdr_fmt_exp[3]);
    fprintf (fp, buf,  f->nel, f->dim, f->lmax);
  }
 
  sprintf (buf,"%s%s",p_hdr_fmt[4],p_hdr_fmt_exp[4]);
  fprintf (fp, buf, f->step);
  sprintf (buf,"%s%s",p_hdr_fmt[5],p_hdr_fmt_exp[5]);
  fprintf (fp, buf, f->time);
  sprintf (buf,"%s%s",p_hdr_fmt[6],p_hdr_fmt_exp[6]);
  fprintf (fp, buf, f->time_step);
  if(f->nz > 1){
    sprintf (buf,"%s%s%s",p_hdr_fmt[7],"%-11.6g   ",p_hdr_fmt_exp[7]);
    fprintf (fp, buf, f->kinvis, f->lz);
  }
  else{
    sprintf (buf,"%s              %s",p_hdr_fmt[7],p_hdr_fmt_exp[7]);
    fprintf (fp, buf, f->kinvis);
  }
  sprintf (buf,"%s%s",p_hdr_fmt[8],p_hdr_fmt_exp[8]);
  fprintf (fp, buf, f->type);
  sprintf (buf,"%s%s",p_hdr_fmt[9],p_hdr_fmt_exp[9]);
  fprintf (fp, buf, f->format);

  int cnt = count_facets(EL->fhead);

  switch (tolower(*f->format)){
  case 'a':
    if(Shuffled){
      for(i = 0; i < pllinfo.nloop; ++i)
	fprintf(fp,"%d ",pllinfo.eloop[i]);
      fputc('\n',fp);
    }

    for(i = 0; i < cnt; ++i)
      fprintf(fp,"%d ",f->size[i]);
    fputc('\n',fp);
    break;
  case 'b':
#ifdef _CRAY
    if(Shuffled) {
      // int on Crays is 64b, short 32b - use short for portable fieldfiles
      spllinfo = (short *) calloc(pllinfo.nloop, sizeof(short));
      for (n = 0; n < pllinfo.nloop; n++)
	spllinfo[n] = (short) pllinfo.eloop[n];
      fwrite(spllinfo, sizeof(short), pllinfo.nloop, fp);
      free (spllinfo);
    }
    ssize = (short *) calloc(cnt, sizeof(short));
    for (n = 0; n < cnt; n++)
      ssize[n] = (short) f->size[n];
    fwrite(ssize, sizeof(short),  cnt, fp);
    free (ssize);
#else
    if(Shuffled)
      fwrite(pllinfo.eloop, sizeof(int),  pllinfo.nloop, fp);
    fwrite(f->size, sizeof(int),  cnt, fp);
#endif
    break;
  default:
    break;
  }
  
  fflush(fp);
#endif
}

/* Write the header and expansion order information */
static int writeData(FILE *fp, Field *f, Element_List *EL){
#ifdef PARALLEL
  int      ntot,nfields;
  char     buf[BUFSIZ];
  register int i, n;

  if(f->state == 'p'){
    fprintf(stderr,"Can't write field in physical space\n");
    exit(-1);
  }

  /* Write the field files  */
  if(f->nz > 1) /* fourier case */
    ntot = f->nz*data_len(f->nel,f->size, EL->fhead);
  else
    ntot = data_len(f->nel,f->size, EL->fhead);
      
  nfields = (int) strlen (f->type);
  
  switch (tolower(*f->format)) {
  case 'a':
    for (i = 0; i < ntot; i++) {
      for (n = 0; n < nfields; n++)
	fprintf (fp, "%#16.10g ", f->data[n][i]);
      fputc ('\n', fp);
    }
    break;
  case 'b':
    for (n = 0; n < nfields; n++)
      if (fwrite (f->data[n], sizeof(double), ntot, fp) != ntot) {
	fprintf  (stderr, "error writing field %c\n", f->type[n]);
	exit(-1);
      }
    break;
  default:
    sprintf  (buf, "unknown format -- %s", f->format);
    fprintf(stderr,buf);
    exit(-1);
    break;
  }

  fflush (fp);

  return nfields;
#else
  return 0;
#endif
}

// Expects the full mesh
/* read in a parallel field using the local to global mapping in pllinfo */
int readFieldP (FILE *fp, Field *f, Element_List *Mesh){
#ifdef PARALLEL
  register int i, j, n, ntot, nfields;
  int      Fourier = 0, skip = 0, *gsize, gnel;
  int      dskip, cnt,el;
  char     buf[BUFSIZ];
#ifdef _CRAY
  short *sgsize;
#endif

  if (fscanf (fp, "%s", buf) == EOF) return 0;            /* session name */
  f->name    = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               /* creation date */
  f->created = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  
  READLINE ("%c"  , &f->state);                  /* simulation parameters */
  /* check to see if it is standard 2d file */
  if(strstr(buf,"Fourier")){
    Fourier = 1;
    fscanf(fp,"%d%d%d%d",&gnel,&f->dim,&f->lmax,&f->nz);fgets(buf,BUFSIZ,fp);
  }
  else{
    fscanf(fp,"%d%d%d",&gnel,&f->dim,&f->lmax);fgets(buf,BUFSIZ,fp);
    f->nz = 1;
  }
  f->nel = pllinfo.nloop;

  READLINE ("%d"  , &f->step);                  
  READLINE ("%lf" , &f->time);
  READLINE ("%lf" , &f->time_step);
  if(Fourier){
    fscanf(fp,"%lf%lf " ,&f->kinvis,&f->lz);
    fgets(buf,BUFSIZ,fp);
  }
  else{
    READLINE ("%lf" , &f->kinvis);
    f->lz = 1.0;
  }

  nfields = gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = (char *)strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("Warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);

  /* allocate memory and load the data */


  int *cfs = pllinfo.cumfacets;
  switch (tolower(*f->format)) {
  case 'a':
    if(f->state == 't'){
      double tmp;

      gsize   = ivector(0,cfs[gnel]-1);
      for(i = 0; i < cfs[gnel]; ++i) fscanf(fp,"%d",gsize+i);
      
      /* copy in relevant data from gsize to f->size */
      f->size = ivector(0,cfs[gnel]-1);
      
      for(i = 0,skip = 0; i < pllinfo.nloop; ++i){
	el = pllinfo.eloop[i];
	icopy(pllinfo.efacets[el], gsize + cfs[el],1,f->size + skip,1);
	skip += pllinfo.efacets[el];
      }
      
      // calculate local total given global size vector
      ntot = 0;
      for(i=0; i<pllinfo.nloop;++i){
	el = pllinfo.eloop[i];
	ntot += f->nz*data_len(1,gsize+cfs[el], Mesh->flist[el]);
      }

      for(i = 0; i < nfields; ++i)
	f->data[i] = dvector (0, ntot-1);
      
      dskip = 0; cnt = 0;
      for(i = 0; i < gnel; ++i){
	ntot = data_len(1,gsize + cfs[i], Mesh->flist[i]);
	if((cnt < pllinfo.nloop)&&(i == pllinfo.eloop[cnt])){
	  for(j = 0; j < ntot; ++j, ++dskip)
	    for(n = 0; n < nfields; ++n)	
	      if (fscanf (fp, "%lf", f->data[n] + dskip) != 1) {
		fprintf(stderr,"Error: reading field %c\n",f->type[n]);
		exit(-1);
	      }
	    ++cnt;
	}
	else{
	  for(j = 0; j < ntot; ++j)
	    for(n = 0; n < nfields; ++n)
	      if (fscanf (fp, "%lf", &tmp) != 1) {
		fprintf(stderr,"Error: reading field %c\n",f->type[n]);
		exit(-1);
	      }
	}
      }
      free(gsize);
    }
    else {
      fputs("Error: reading field file in physical space is "
	    "not implemented \n",stderr);
      exit(-1);
    }
    break;
  case 'b':
    if(f->state == 't'){
      double *tmp;

      gsize   = ivector(0,cfs[gnel]-1);
      tmp = dvector(0,f->lmax*f->lmax*f->lmax-1);


#ifdef _CRAY
      // int on Crays is 64b, short 32b - use short for portable fieldfiles
      sgsize = (short *) calloc(cfs[gnel], sizeof(short));
      if(fread(sgsize, sizeof(short), (size_t) cfs[gnel], fp) != cfs[gnel]){
	fprintf (stderr, "error reading size of field ");
	exit(-1);
      }
      const int gcfs = cfs[gnel];
      for (n = 0; n < gcfs; n++)
	gsize[n] = (int) sgsize[n];
      free (sgsize);
#else
      if(fread(gsize, sizeof(int),(size_t)cfs[gnel],fp) != cfs[gnel]){
	fprintf (stderr, "error reading size of field ");
	exit(-1);
      }
#endif
      
      /* copy in relevant data from gsize to f->size */
      f->size = ivector(0,cfs[gnel]-1);
      for(i = 0, skip = 0; i < pllinfo.nloop; ++i){
	el = pllinfo.eloop[i];
	icopy(pllinfo.efacets[el], gsize + cfs[el], 1, f->size + skip,1);
	skip += pllinfo.efacets[el];
      }

      // calculate local total given global size vector
      ntot = 0;
      for(i=0; i<pllinfo.nloop;++i) {
	el = pllinfo.eloop[i];
	ntot += f->nz*data_len(1,gsize+cfs[el],Mesh->flist[el]);
      }

      for(i = 0; i < nfields; ++i)
	if(!(f->data [i]))
	  f->data[i] = dvector (0, ntot-1);
    
      for(n = 0; n < nfields; ++n){
	dskip = cnt = 0;
	for(i = 0; i < gnel; ++i){
	  ntot = data_len(1,gsize + cfs[i],Mesh->flist[i]);
	  if((cnt < pllinfo.nloop)&&(i == pllinfo.eloop[cnt])){
	    if(fread(f->data[n]+dskip,sizeof(double),(size_t)ntot,fp) != ntot){
	      fprintf (stderr, "error reading field %c", f->type[n]);
	      exit(-1);
	    }
	    dskip += ntot;
	    ++cnt;
	  }
	  else{
	    if (fread (tmp, sizeof(double), (size_t) ntot, fp) != ntot){
	      fprintf (stderr, "error reading field %c", f->type[n]);
	      exit(-1);
	    }
	  }
	}
      }
      free(gsize); free(tmp);
    }
    else {
      fputs("Error: reading field file in physical space is "
	    "not implemented \n",stderr);
      exit(-1);
    }
  break;
  default:
    fprintf (stderr, "unknown field format -- %s\n", f->format);
    exit(-1);
    break;
  }
  
  return nfields;
#else
  return 0;
#endif
}


int readHeaderP(FILE* fp, Field *f, Element_List *Mesh){
  register int  nfields;
  char     buf[BUFSIZ];
#ifdef _CRAY
  short *ssize, *semap;
#endif
  int      shuffled = 0,i, cnt;

  if (fscanf (fp, "%s", buf) == EOF) return 0;            /* session name */
  
  f->name    = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               /* creation date */
  f->created = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  
  READLINE ("%c"  , &f->state);                  /* simulation parameters */

  if(strstr(buf,"Shuffled")) shuffled = 1; /*check for non-sequential setup*/
  fscanf(fp,"%d%d%d",&f->nel,&f->dim,&f->lmax);fgets(buf,BUFSIZ,fp);
  f->nz = 1;
  
  READLINE ("%d"  , &f->step);                  
  READLINE ("%lf" , &f->time);
  READLINE ("%lf" , &f->time_step);
  READLINE ("%lf" , &f->kinvis);

  nfields = gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = (char *)strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("Warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);
  
  f->emap = ivector(0,f->nel-1);
  switch (tolower(*f->format)) {
  case 'a':
    if(option("ALLTET")){

      cnt = f->nel*4;
      f->size = ivector(0,cnt-1);
      for(i = 0; i < cnt; ++i) fscanf(fp,"%d",f->size+i);

      if(shuffled) for(i = 0; i < f->nel; ++i) fscanf(fp,"%d",f->emap + i);
      else for(i = 0; i < f->nel; ++i) f->emap[i] = i;
      
    }
    else{
      if(shuffled) for(i = 0; i < f->nel; ++i) fscanf(fp,"%d",f->emap + i);
      else for(i = 0; i < f->nel; ++i) f->emap[i] = i;
      
      cnt = count_facets(Mesh, f->emap, f->nel);
      f->size = ivector(0,cnt-1);
      for(i = 0; i < cnt; ++i) fscanf(fp,"%d",f->size+i);
    }
    break;
  case 'b':
#ifdef _CRAY
    if(option("ALLTET"))
      fprintf(stderr, "Shuffled file not implemented on cray routine\n");

    register int n;
    if(shuffled) {
      // int on Crays is 64b, short 32b - use short for portable fieldfiles
      semap = (short *) calloc(f->nel, sizeof(short));
      fread(semap,sizeof(short),(size_t)f->nel,fp);
      for (n = 0; n < f->nel; n++)
	f->emap[n] = (int) semap[n];
      free (semap);
    } else
      for(i = 0; i < f->nel; ++i) f->emap[i] = i;

    cnt = count_facets(Mesh, f->emap, f->nel);
    f->size = ivector(0,cnt-1);

    ssize = (short *) calloc(cnt, sizeof(short));
    if(fread(ssize, sizeof(short),(size_t)cnt,fp) != cnt){
      fprintf (stderr, "error reading size of field \n");
      exit(-1);
    }
    for (n = 0; n < cnt; n++)
      f->size[n] = (int) ssize[n];
    free (ssize);    
#else
    if(option("ALLTET")){
      cnt = f->nel*4; 
      f->size = ivector(0,cnt-1);
      
      if(fread(f->size, sizeof(int),(size_t)cnt,fp) != cnt){
	fprintf (stderr, "error reading size of field \n");
	exit(-1);
      } 
      if(shuffled)
	fread(f->emap,sizeof(int),(size_t)f->nel,fp);
      else
	for(i = 0; i < f->nel; ++i) f->emap[i] = i;
    }
    else{
      if(shuffled)
	fread(f->emap,sizeof(int),(size_t)f->nel,fp);
      else
	for(i = 0; i < f->nel; ++i) f->emap[i] = i;
      
      cnt = count_facets(Mesh, f->emap, f->nel);
      f->size = ivector(0,cnt-1);
      
      if(fread(f->size, sizeof(int),(size_t)cnt,fp) != cnt){
	fprintf (stderr, "error reading size of field \n");
	exit(-1);
      }
    }
#endif
    break;
  default:
    fprintf (stderr, "unknown field format -- %s\n", f->format);
    exit(-1);
    break;
  }
  rewind(fp);

  return nfields;
}

// Assumes global Mesh
// This is for the utilites to assemble global fld

Field *readFieldFiles(int nfiles, char *name, Element_List *Mesh){
  int i,j,k,n;
  int nfields, ntot;

  Field **Fld = (Field**) calloc(nfiles, sizeof(Field*));

  char hdr_name[BUFSIZ];
  char dat_name[BUFSIZ];
  FILE *hdr_fp, *dat_fp;
  Element *E;

  // setup global cumulative lists for facet count
  int *Nfcts = ivector(0, Mesh->nel-1);
  for(i=0;i<Mesh->nel;++i){
    E = Mesh->flist[i];
      
    Nfcts[i] = E->Nedges + E->Nfaces;
    if(E->dim() == 3)
      ++Nfcts[i];
  }
  
  int *cumNfcts = ivector(0, Mesh->nel);
  cumNfcts[0] = 0;
  for(i=1;i<Mesh->nel+1;++i)
    cumNfcts[i] = cumNfcts[i-1] + Nfcts[i-1];
  
  // read header files
  Fld = (Field**) calloc(nfiles, sizeof(Field*));

  for(i=0;i<nfiles;++i){

    sprintf(hdr_name, "%s.hdr.%d", name, i);
    
    hdr_fp = fopen(hdr_name, "r");
    if(!hdr_fp){
      fprintf(stderr, "WARNING readFieldFiles: error opening hdr %d\n",i);
      //      exit(-1);
      Fld[i] = NULL;
    }
    else{
      Fld[i]  = (Field*) calloc(1, sizeof(Field));
      nfields = readHeaderP(hdr_fp, Fld[i], Mesh);
    }
  }

  // find out data length for each global element
  int *datalen = ivector(0, Mesh->nel-1);
  int skip,elmt,gskip;
  //  izero(Mesh->nel, datalen, 1);
  for(k=0;k<Mesh->nel;++k)
    datalen[k] = Mesh->flist[k]->Nverts;
  
  for(i=0;i<nfiles;++i){
    skip = 0;
    
    if(Fld[i])
      for(k=0;k<Fld[i]->nel;++k){
	datalen[Fld[i]->emap[k]] = data_len(1, Fld[i]->size + skip, 
					    Mesh->flist[Fld[i]->emap[k]]);
	skip += Nfcts[Fld[i]->emap[k]];
      }
  }
  
  // make cumalative length for global elements
  int *cumdatalen = ivector(0, Mesh->nel);
  cumdatalen[0] = 0;
  for(k=1;k<Mesh->nel+1;++k)
    cumdatalen[k] = cumdatalen[k-1] + datalen[k-1];
  
  // set up global field storage
  Field *GFld = (Field*) calloc(1, sizeof(Field));
  memcpy(GFld, *Fld, sizeof(Field));

  for(i=0;i<nfields;++i){
    GFld->data[i] = dvector(0, cumdatalen[Mesh->nel]-1);
    dzero(cumdatalen[Mesh->nel], GFld->data[i], 1);
  }
  
  GFld->size  = ivector(0, cumNfcts[Mesh->nel]-1);
  GFld->nel   = Mesh->nel;

  skip = 0;

  izero(cumNfcts[Mesh->nel], GFld->size, 1);

    
  for(i=0;i<nfiles;++i){
    if(Fld[i]){
    sprintf(dat_name, "%s.dat.%d", name, i);
    
    dat_fp = fopen(dat_name, "r");
    if(!dat_fp){
      fprintf(stderr, "WARNING readFieldFiles: error opening dat %d\n",i);
      exit(-1);
    }
    else{
      
      skip = 0;
      for(k=0;k<Fld[i]->nel;++k){
	elmt = Fld[i]->emap[k];
	icopy(Nfcts[elmt], Fld[i]->size+skip, 1,  
	      GFld->size + cumNfcts[elmt], 1);
	skip += Nfcts[elmt];
      }
      switch (tolower(*Fld[i]->format)) {
      case 'a':
	for(k=0;k<Fld[i]->nel;++k){
	  elmt   = Fld[i]->emap[k];
	  ntot   = datalen[elmt];
	  gskip  = cumdatalen[elmt];
	  for (j = 0; j < ntot; j++) {
	    for (n = 0; n < nfields; n++)
	      if (fscanf (dat_fp, "%lf", GFld->data[n] + gskip + j) != 1) {
		fprintf(stderr,"Error: reading field %c, point %d\n",
			Fld[i]->type[n],j+1);
		exit(-1);
	      }
	  }
	}
	break;
      case 'b':
	for (n = 0; n < nfields; n++) {
	  for(k=0;k<Fld[i]->nel;++k){
	    elmt   = Fld[i]->emap[k];
	    ntot   = datalen[elmt];
	    gskip  = cumdatalen[elmt];
	    if(fread(GFld->data[n]+gskip, sizeof(double), (size_t) ntot, dat_fp)!=ntot){
	      fprintf (stderr, "error reading field %c", GFld->type[n]);
	      exit(-1);
	    }
	  }
	}
	break;
      }
    }
  }
  }
  return GFld;
}

/* ---------------------------------------------------------------------- *
 * readField() -- Load a field file                                       *
 *                                                                        *
 * This function loads a field file into a Field structure.  The format   *
 * is a simple array, stored in 64-bit format.                            *
 *                                                                        *
 * Return value: number of variables loaded, 2 for EOF, -1 for error.     *
 * ---------------------------------------------------------------------- */

#define READLINE(fmt,arg) \
  if (fscanf(fp,fmt,arg)==EOF) return -1; fgets(buf,BUFSIZ,fp)

int readHeader_old (FILE *fp, Field *f){
  register int i;
  int      Fourier = 0,lskip;
  int      shuffled = 0;
  char     buf[BUFSIZ];

  if (fscanf (fp, "%s", buf) == EOF) return 2;           

  /* check to see if name is defined and if so free the field */
  /* this can happen if reading in multiple dumps             */
  if(f->name) freeField(f);

  /* session name */
  f->name    = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               
  /* creation date */
  f->created = (char *)strdup (buf); fgets (buf, BUFSIZ, fp);
  
  READLINE ("%c"  , &f->state);                  /* simulation parameters */

  if(f->state == 'p'){
    fputs("Error: reading field file in physical space is "
	  "not implemented \n",stderr);
    exit(-1);
  }

  /* check to see if it is standard 2d file */
  if(strstr(buf,"Fourier")){
    Fourier = 1;
    fscanf(fp,"%d%d%d%d",&f->nel,&f->dim,&f->lmax,&f->nz);fgets(buf,BUFSIZ,fp);
  }
  else{
    if(strstr(buf,"Shuffled")) shuffled = 1; /*check for non-sequential setup*/
    fscanf(fp,"%d%d%d",&f->nel,&f->dim,&f->lmax);fgets(buf,BUFSIZ,fp);
    f->nz = 1;
  }
  
  lskip = (f->dim == 2) ?  4:11;

  READLINE ("%d"  , &f->step);                  
  READLINE ("%lf" , &f->time);
  READLINE ("%lf" , &f->time_step);

  if(Fourier){
    fscanf(fp,"%lf%lf " ,&f->kinvis,&f->lz);
    fgets(buf,BUFSIZ,fp);
  }
  else{
    READLINE ("%lf" , &f->kinvis);
    f->lz = 1.0;
  }

  gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = (char *)strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("Warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);

  if(shuffled) f->emap = ivector(0,f->nel-1);
    
  switch (tolower(*f->format)) {
  case 'a':
    f->size = ivector(0,lskip*f->nel-1);
    for(i = 0; i < lskip*f->nel; ++i) fscanf(fp,"%d",f->size+i);
    if(shuffled) for(i = 0; i < f->nel; ++i) fscanf(fp,"%d",f->emap + i);
    break;
  case 'b':
    f->size = ivector(0,lskip*f->nel-1);
    if(fread(f->size, sizeof(int),(size_t)lskip*f->nel,fp)!= lskip*f->nel){
      fprintf (stderr, "error reading size of field \n");
      exit(-1);
    }
    if((shuffled)&&(fread(f->emap,sizeof(int),(size_t)f->nel,fp) != f->nel)){
      fprintf (stderr, "error reading emap \n ");
      exit(-1);
    }
    break;
  default:
    fprintf (stderr, "unknown field format -- %s\n", f->format);
    exit(-1);
    break;
  }
  return shuffled;
}

int readField_old (FILE *fp, Field *f, Element_List *Mesh){
  register int i, n, ntot, nfields;
  int      dlen,lskip;
  int      shuffled = 0;
  char     buf[BUFSIZ];


  lskip = (Mesh->fhead->dim() == 2) ? 4:11;
  
  shuffled = readHeader_old(fp,f);
  if(shuffled == 2) return 0; /* if shuffled == 2 the at EOF */

  nfields  = strlen(f->type);

  if(!shuffled){  /* standard input routine format */
#ifdef PARALLEL
    DO_PARALLEL{ /* read in relevant elements */
      register int j;
      int *gsize,dskip,cnt;

      switch (tolower(*f->format)) {
      case 'a':
	if(f->state == 't'){
	  double tmp;
	  
	  gsize  = f->size;
	  
	  /* copy in relevant data from gsize to f->size */
	  f->size = ivector(0,lskip*pllinfo.nloop-1);
	  for(i = 0; i < pllinfo.nloop; ++i)
	    icopy(lskip,gsize + pllinfo.eloop[i]*lskip,1,f->size + i*lskip,1);
	  
	  ntot = data_len(pllinfo.nloop,f->size, Mesh->fhead);
	  
	  for(i = 0; i < nfields; ++i)
	    f->data[i] = dvector (0, ntot-1);
	  
	  dskip = 0; cnt = 0;
	  for(i = 0; i < f->nel; ++i){
	    ntot = data_len(1,gsize + i*lskip, Mesh->fhead);
	    if((cnt < pllinfo.nloop)&&(i == pllinfo.eloop[cnt])){
	      for(j = 0; j < ntot; ++j, ++dskip)
		for(n = 0; n < nfields; ++n)	
		  if (fscanf (fp, "%lf", f->data[n] + dskip) != 1) {
		    fprintf(stderr,"Error: reading field %c\n",f->type[n]);
		    exit(-1);
		  }
	      ++cnt;
	    }
	    else
	      for(j = 0; j < ntot; ++j)
		for(n = 0; n < nfields; ++n)
		  if (fscanf (fp, "%lf", &tmp) != 1) {
		    fprintf(stderr,"Error: reading field %c\n",f->type[n]);
		    exit(-1);
		  }
	  }
	  free(gsize);
	}
	else {
	  fputs("Error: reading field file in physical space is "
		"not implemented \n",stderr);
	  exit(-1);
	}
	break;
      case 'b':
	if(f->state == 't'){
	  double *tmp;

	  if(Mesh->fhead->dim() == 2)
	    tmp = dvector(0,f->lmax*(f->lmax+1)/2-1);
	  else
	    tmp = dvector(0,f->lmax*(f->lmax+1)*(f->lmax+2)/6-1);

	  
	  gsize = f->size;
	  
	  /* copy in relevant data from gsize to f->size */
	  f->size = ivector(0,lskip*pllinfo.nloop-1);
	  for(i = 0; i < pllinfo.nloop; ++i)
	    icopy(lskip,gsize + pllinfo.eloop[i]*lskip,1,f->size + i*lskip,1);
	  
	  ntot = data_len(pllinfo.nloop,f->size, Mesh->fhead);
	  
	  for(i = 0; i < nfields; ++i)
	    if(!(f->data [i]))
	      f->data[i] = dvector (0, ntot-1);
	  
	  for(n = 0; n < nfields; ++n){
	    dskip = cnt = 0;
	    for(i = 0; i < f->nel; ++i){
	      ntot = data_len(1,gsize + i*lskip, Mesh->fhead);
	      if((cnt < pllinfo.nloop)&&(i == pllinfo.eloop[cnt])){
		
		if(fread(f->data[n]+dskip,sizeof(double),
			 (size_t)ntot,fp) != ntot){
		  fprintf (stderr, "error reading field %c", f->type[n]);
		  exit(-1);
		}
		dskip += ntot;
		++cnt;
	      }
	      else{
		if (fread (tmp, sizeof(double), (size_t) ntot, fp) != ntot){
		  fprintf (stderr, "error reading field %c", f->type[n]);
		  exit(-1);
		}
	      }
	    }
	  }
	  free(gsize); free(tmp);
	}
	f->nel = pllinfo.nloop;
      }
    }
    else
#endif
      {
	/* allocate memory and load the data */
      dlen = data_len(f->nel,f->size, Mesh->fhead);
      ntot = f->nz*dlen;
      
      for(i = 0; i < nfields; ++i)
	if(!(f->data [i]))
	  f->data[i] = dvector (0, ntot-1);
      
      
      switch (tolower(*f->format)) {
      case 'a':
	for (i = 0; i < ntot; i++) {
	  for (n = 0; n < nfields; ++n)
	  if (fscanf (fp, "%lf", f->data[n] + i) != 1) {
	    fprintf(stderr,"Error: reading field %c, point %d\n",
		    f->type[n],i+1);
	    exit(-1);
	  }
	}
	
	if(f->nz > 1){ /* zero out second fourier plane */
	  for(n = 0; n < nfields; ++n)
	    dzero(dlen,f->data[n] + dlen,1);
	}
	fgets (buf, BUFSIZ, fp);
	break;
      case 'b':
	for (n = 0; n < nfields; n++) {
	  if (fread (f->data[n], sizeof(double), (size_t) ntot, fp) != ntot) {
	    fprintf (stderr, "error reading field %c", f->type[n]);
	  }
	}
	
	if(f->nz > 1){ /* zero out second fourier plane */
	  for(n = 0; n < nfields; ++n)
	    dzero(dlen,f->data[n] + dlen,1);
	}
	break;
      default:
	fprintf (stderr, "unknown field format -- %s\n", f->format);
	exit(-1);
	break;
      }
    }
  }
  else{ /* element storage is non-sequential and stored on elemental basis */
    
    /* NOTE: not set up for fourier with parallel split in each plane */
#ifdef PARALLEL
    DO_PARALLEL{ /* just read in relevant elements */
      int j;
      int *gskip,cnt,eid,skip;
      int *peid,*dskip,*len,*shuf;
    
      gskip   = f->size;   
      f->size = ivector(0,lskip*pllinfo.nloop-1);
   
      /* determine the length of every element as stored */
      /* and generate map of unshuffled indices          */
      len   = ivector(0,f->nel-1);
      shuf  = ivector(0,f->nel-1);
      
      for(i = 0; i < f->nel; ++i){
	len  [i] = data_len(1,gskip + lskip*i, Mesh->fhead);
	shuf [f->emap[i]] = i;
      }

      dskip = ivector(0,pllinfo.nloop-1);
      peid  = ivector(0,pllinfo.nloop-1);

      /* assemble the storage order of the elements we want 
         and the location that data should be stored       */
      for(i = 0,ntot = 0; i < pllinfo.nloop; ++i){
	peid[i]  = shuf[pllinfo.eloop[i]];
	dskip[i] = ntot;
	ntot += len[peid[i]];
      }
      
      /* reform f->size */
      for(i = 0; i < pllinfo.nloop; ++i)
	icopy(lskip,gskip + peid[i]*lskip,1,f->size + i*lskip,1);
      
      /* declare data memory */
      for(i = 0; i < nfields; ++i)
	f->data[i] = dvector (0, ntot-1);
      
      /* at this stage we have a local index peid[i] of the elements 
	 we want to extract in terms of the order they are stored 
	 and a position vector as to where they need to be stored.
	 To make the extraction easier we now need to re-order
	 these vectors so they are sequential in terms of the 
	 ordering of the stored elements                         */
      
      for(i = 0; i < pllinfo.nloop; ++i){
	/* find lowest index between i and nloop */
	for(cnt = j = i; j < pllinfo.nloop; ++j)
	  cnt = (peid[j] < peid[cnt])? j:cnt;
	eid  = peid[cnt];
	skip = dskip[cnt];
	
	/* put this id at position i and push all others up the list */
	for(j = cnt; j > i; --j){
	  peid [j] = peid [j-1];
	  dskip[j] = dskip[j-1];
	}
	peid[i]  = eid;
	dskip[i] = skip;
      }

      switch(tolower(*f->format)){
      case 'a':
	{
	  double tmp;
	  cnt = 0;
	  for(i = 0; i < f->nel; ++i) 
	    if((cnt < pllinfo.nloop)&&(i == peid[cnt])){
	      for(j = 0; j < len[i]; ++j) 
		for(n = 0; n < nfields; ++n)
		  if (fscanf (fp, "%lf", f->data[n] + dskip[cnt]+j) != 1) {
		    fprintf(stderr,"Error: reading field %c, element %d\n",
			    f->type[n],i+1);
		    exit(-1);
		  }
	      ++cnt;
	    }
	    else{
	      for(j = 0; j < len[i]; ++j)
		for(n = 0; n < nfields; ++n)
		  if (fscanf (fp, "%lf", &tmp) != 1) {
		    fprintf(stderr,"Error: reading field %c, element %d\n",
			    f->type[n],i+1);
		    exit(-1);
		  }
	    }
	}
	break;
      case 'b':
	{
	  double *tmp;
	  if(Mesh->fhead->dim() == 2)
	    tmp = dvector(0,f->lmax*(f->lmax+1)/2-1);
	  else
	    tmp = dvector(0,f->lmax*(f->lmax+1)*(f->lmax+2)/6-1);
	  
	  
	  cnt = 0;
	  for(i = 0; i < f->nel; ++i) 
	    if((cnt < pllinfo.nloop)&&(i == peid[cnt])){
	      for(n = 0; n < nfields; ++n)
		if(fread (f->data[n]+dskip[cnt],sizeof(double),
			  (size_t)len[i],fp) != len[i]){
		  fprintf(stderr, "error reading field %c, element %d", 
			  f->type[n],i+1);
		  exit(-1);
		}
	      ++cnt;
	    }
	    else{
	      for(n = 0; n < nfields; ++n)
		if(fread (tmp,sizeof(double),(size_t)len[i],fp) != len[i]){
		  fprintf(stderr, "error reading field %c, element %d", 
			  f->type[n],i+1);
		  exit(-1);
		}
	    }
	}
	break;
      default:
	fprintf (stderr, "unknown field format -- %s\n", f->format);
	exit(-1);
	break;
      }
      
      f->nel = pllinfo.nloop;
      free(gskip); free(dskip); free(peid); free(len); free(shuf);
    }
    else 
#endif
      { /* read in all elements unsorting as we go! */
      int j;
      int *dskip,*len;
    
      /* determine the length and staring position of every piece of data */
      dskip = ivector(0,f->nel-1);
      len   = ivector(0,f->nel-1);
      
      /* find the length of the i th element and put in unshuffled order */
      for(i = 0; i < f->nel; ++i)
	len  [f->emap[i]] = data_len(1,f->size + lskip*i, Mesh->fhead);

      /* find starting location for every data point */
      for(i = 0; i < f->nel; ++i){
	dskip[i] = 0;
	for(j = 0; j < f->emap[i]; ++j)
	  dskip[i] += len[j];
      }
      ntot = data_len(f->nel,f->size,Mesh->fhead);
      
      for(i = 0; i < nfields; ++i)
	if(!(f->data [i]))
	  f->data[i] = dvector (0, ntot-1);
      
      switch(tolower(*f->format)){
      case 'a':
	for (i = 0; i < f->nel; ++i) 
	  for (j = 0; j < len[f->emap[i]]; ++j) 
	    for (n = 0; n < nfields; ++n)
	      if (fscanf (fp, "%lf", f->data[n] + dskip[i]+j) != 1) {
		fprintf(stderr,"Error: reading field %c, element %d\n",
			f->type[n],i+1);
		exit(-1);
	      }
	break;
      case 'b':
	for (i = 0; i < f->nel; ++i) 
	  for (n = 0; n < nfields; ++n)
	    if(fread (f->data[n]+dskip[i], sizeof(double), 
		      (size_t)len[f->emap[i]],fp) != len[f->emap[i]]) {
	      fprintf(stderr, "Error reading field %c, element %d\n", 
		      f->type[n],i+1);
	      exit(-1);
	    }
	break;
      default:
	fprintf (stderr, "unknown field format -- %s\n", f->format);
	exit(-1);
	break;
      }

      free(len); free(dskip);
    }
  }
  
  return nfields;
}

/* 

Function name: Element::Copy_field

Function Purpose:
 Copy field coefficients from field file format modal storage into element modal storage.

Argument 1: double *data
Purpose: 
 Contains the modal storage read from file.

Argument 2: int *size
Purpose: 
 Contains the degree vector for the modal data read from file.

Function Notes:

*/

void Tri::Copy_field(double *data, int *size){
  int i,j,l,cl;
  
  for(i = 0; i < Nverts; ++i){
    vert[i].hj[0] = *data;
    ++data;
  }
   
  /* copy edges */
  for (i = 0; i < Nedges; ++i){
    if(l = edge[i].l){
      dzero(l,edge[i].hj,1);
      cl = min(l,*size);
      dcopy(cl,data,1,edge[i].hj,1);
    }
    data += *size;
    ++size;
  }
  
  /* copy faces */
  for(j = 0; j < Nfaces; ++j){
    if(l = face[j].l){
      dzero(l*(l+1)/2,*face[j].hj,1);
      cl = min(l,*size);
      for(i = 0; i < cl; ++i){
	dcopy(cl-i,data,1,face[j].hj[i],1);
	data += *size-i;
      }
      for(i = cl; i < *size; ++i)
	data += *size-i;
    }
    else
      for(i = 0; i < *size; ++i)
	data += *size-i;
    ++size;
  }
}




void Quad::Copy_field(double *data, int *size){
  int i,j,l,cl;
  
  for(i = 0; i < Nverts; ++i){
    vert[i].hj[0] = *data;
    ++data;
  }
   
  /* copy edges */
  for (i = 0; i < Nedges; ++i){
    if(l = edge[i].l){
      dzero(l,edge[i].hj,1);
      cl = min(l,*size);
      dcopy(cl,data,1,edge[i].hj,1);
    }
    data += *size;
    ++size;
  }
  
  /* copy faces */
  for(j = 0; j < Nfaces; ++j){
    if(l = face[j].l){
      dzero(l*l,*face[j].hj,1);
      cl = min(l,*size);
      for(i = 0; i < cl; ++i){
	dcopy(cl,data,1,face[j].hj[i],1);
	data += *size;
      }
      for(i = cl; i < *size; ++i)
	data += *size;
    }
    else
      for(i = 0; i < *size; ++i)
	data += *size;
    ++size;
  }
}




void Tet::Copy_field(double *data, int *size){
  int i,j,l,cl;
  
  dzero(Nmodes, vert[0].hj, 1);
  
  for(i = 0; i < Nverts; ++i){
    vert[i].hj[0] = *data;
    ++data;
  }
   
  /* copy edges */
  for (i = 0; i < Nedges; ++i){
    if(l = edge[i].l){
      cl = min(l,*size);
      dcopy(cl,data,1,edge[i].hj,1);
    }
    data += *size;
    ++size;
  }
  
  /* copy faces */
  for(j = 0; j < Nfaces; ++j){
    if(l = face[j].l){
      cl = min(l,*size);
      for(i = 0; i < cl; ++i){
	dcopy(cl-i,data,1,face[j].hj[i],1);
	data += *size-i;
      }
      for(i = cl; i < *size; ++i)
	data += *size-i;
    }
    else
      for(i = 0; i < *size; ++i)
	data += *size-i;
    ++size;
  }

  if(l = interior_l){ 
    cl = min(l,*size);
    for(i = 0; i < cl; ++i){
      for(j = 0; j < cl-i; ++j){
	dcopy(cl-i-j,data,1,hj_3d[i][j],1);
	data += *size-i-j;
      }
      for(j = cl-i; j < *size-i; ++j)
	data += *size-i-j;
    }
	
    for(i = cl; i < *size; ++i)
      for(j = 0; j < *size-i; ++j)
	data += *size-i-j;
  }
}




void Pyr::Copy_field(double *data, int *size){
  int i,j,l,cl;
  // clear modes
  dzero(Nmodes, vert[0].hj, 1);
   
  for(i = 0; i < NPyr_verts; ++i){
    vert[i].hj[0] = *data;
    ++data;
  }
   
  /* copy edges */
  for (i = 0; i < NPyr_edges; ++i){
    if(l = edge[i].l){
      cl = min(l,*size);
      dcopy(cl,data,1,edge[i].hj,1);
    }
    data += *size;
    ++size;
  }
  
  /* copy faces */
  for(j = 0; j < NPyr_faces; ++j){
    if(Nfverts(j) == 3){
      if(l = face[j].l){
	cl = min(l,*size);
	for(i = 0; i < cl; ++i){
	  dcopy(cl-i,data,1,face[j].hj[i],1);
	  data += *size-i;
	}
	for(i = cl; i < *size; ++i)
	  data += *size-i;
      }
      else
	for(i = 0; i < *size; ++i)
	  data += *size-i;
    }
    else{
      if(l = face[j].l){
	cl = min(l,*size);
	for(i = 0; i < cl; ++i){
	  dcopy(cl,data,1,face[j].hj[i],1);
	  data += *size;
	}
	for(i = cl; i < *size; ++i)
	  data += *size;
      }
      else
	for(i = 0; i < *size; ++i)
	  data += *size;
    }
    ++size;
  }

  if(l = interior_l){
    cl = min(l,*size);
    for(i = 0; i < cl; ++i){
      for(j = 0; j < cl-i; ++j){
	dcopy(cl-i-j,data,1,hj_3d[i][j],1);
	data += *size-i-j;
      }
      for(j = cl-i; j < *size-i; ++j)
	data += *size-i-j;
    }
	
    for(i = cl; i < *size; ++i)
      for(j = 0; j < *size-i; ++j)
	data += *size-i-j;
  }
}




void Prism::Copy_field(double *data, int *size){
  int i,j,l,cl;

  dzero(Nmodes, vert[0].hj, 1);
  
  for(i = 0; i < NPrism_verts; ++i){
    vert[i].hj[0] = *data;
    ++data;
  }
   
  /* copy edges */
  for (i = 0; i < NPrism_edges; ++i){
    if(l = edge[i].l){
      cl = min(l,*size);
      dcopy(cl,data,1,edge[i].hj,1);
    }
    data += *size;
    ++size;
  }
  
  /* copy faces */
  for(j = 0; j < NPrism_faces; ++j){
    if(Nfverts(j) == 3){
      if(l = face[j].l){
	cl = min(l,*size);
	for(i = 0; i < cl; ++i){
	  dcopy(cl-i,data,1,face[j].hj[i],1);
	  data += *size-i;
	}
	for(i = cl; i < *size; ++i)
	  data += *size-i;
      }
      else
	for(i = 0; i < *size; ++i)
	  data += *size-i;
    }
    else{
      if(l = face[j].l){
	cl = min(l,*size);
	for(i = 0; i < cl; ++i){
	  dcopy(cl,data,1,face[j].hj[i],1);
	  data += *size;
	}
	for(i = cl; i < *size; ++i)
	  data += *size;
      }
      else
	for(i = 0; i < *size; ++i)
	  data += *size;
    }
    ++size;
  }

  if((l = interior_l)>1){
    cl = min(l,*size);
    for(i = 0; i < cl-1; ++i){
      for(j = 0; j < cl; ++j){
	dcopy(cl-i-1,data,1,hj_3d[i][j],1);
	data += *size-i-1;
      }
      for(j = cl; j < *size; ++j)
	data += *size-i-1;
    }
	
    for(i = cl; i < *size-1; ++i)
      for(j = 0; j < *size; ++j)
	data += *size-i-1;
  }
}




void Hex::Copy_field(double *data, int *size){
  int i,j,l,cl;

  dzero(Nmodes, vert[0].hj, 1);
  
  for(i = 0; i < NHex_verts; ++i){
    vert[i].hj[0] = *data;
    ++data;
  }
   
  /* copy edges */
  for (i = 0; i < NHex_edges; ++i){
    if(l = edge[i].l){
      cl = min(l,*size);
      dcopy(cl,data,1,edge[i].hj,1);
    }
    data += *size;
    ++size;
  }
  
  /* copy faces */
  for(j = 0; j < NHex_faces; ++j){
    if(l = face[j].l){
      cl = min(l,*size);
      for(i = 0; i < cl; ++i){
	dcopy(cl,data,1,face[j].hj[i],1);
	data += *size;
      }
      for(i = cl; i < *size; ++i)
	data += *size;
    }
    else
      for(i = 0; i < *size; ++i)
	data += *size;

    ++size;
  }

  if((l = interior_l)>1){
    cl = min(l,*size);
    for(i = 0; i < cl; ++i){
      for(j = 0; j < cl; ++j){
	dcopy(cl,data,1,hj_3d[i][j],1);
	data += *size;
      }
      for(j = cl; j < *size; ++j)
	data += *size;
    }
	
    for(i = cl; i < *size; ++i)
      for(j = 0; j < *size; ++j)
	data += *size;
  }
}

void Element::Copy_field(double *, int *){ERR;}
