/*
 * Mapping stuff
 *
 * RCS Information
 * ---------------
 * $Author: ce107 $
 * $Date: 1997/03/07 18:21:23 $
 * $Source: /users/ce107/Prism/flex/RCS/flex.c,v $
 * $Revision: 1.6 $
 * -------------------------------------------------------------- */

#include <mpi.h>
#include <time.h>
#include <string.h>
#include "nektarF.h"
#include "map.h"
#include <rfftw.h>
extern rfftw_plan rplan, rplan_inv;
extern rfftw_plan rplan32, rplan_inv32;
extern "C"
{
  void sinft (int, double *, int);
  void cosft (int, double *, int);
}

static int forcx, forcy;

void writeheader (FILE *fp, char *name, int lmax, double lz, int nz, 
			 int nel, int step, double t, char *typelist);

void writedog (FILE *fp, Mapping *mapx, Mapping *mapy);
void writedog2 (FILE *fp, Mapping *mapx, Mapping *mapy, double *basep);
void writecab (FILE *fp, Mapping *mapx, Mapping *mapy);
void setMap (Mapping *mapx, Mapping *mapy); /* update the map using mapping() */
			 
static double *oldd_1, *oldt_1, *oldf_1;
static double *Mat_plus, *Mat_inv, *RHS_f;
static int  map_step = 0;
void generate_structure_matrix( double wn, double zeta, double *Tension, double *Stiffness );
void Derivatives_PinnedPinned(double *oldd, double *oldzz, double *oldzzzz, int NZ, int NZM, double dz);

void newmark (double *x, double *v, double *a,
	      double  f, double dt, double  m, double wn, double zeta);
#ifdef FLOW_CONTROL
double *yd, *td;
double current_amp, next_amp;
#endif

#ifdef MAP
/* ------------------------------------------------------------------------ *
 *
 * INITIALIZATION
 * ==============
 *
 * InitMap()
 * allocate_map()
 *
 * ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ *
 * InitMap() -- allocate map and set map or read in map file
 * ------------------------------------------------------------------------ */

void InitMap (Domain *omega)
{
  int  NZ = option("NZTOT");

  omega->mapx = allocate_map (NZ);
  omega->mapy = allocate_map (NZ);

  #ifdef FLOW_CONTROL
   int ncoeff = iparam("NHF_COEFF");
   yd = dvector(0, ncoeff-1);
   td = dvector(0, ncoeff-1);
   current_amp = 0.0;
   next_amp = 0.0;
  #endif

  forcx = (int) dparam("FORCX");
  forcy = (int) dparam("FORCY");
 
  if (forcx < -7 || forcx > 7) {
    fprintf(stderr, "Acceptable values for FORCX are -7 to 7! Using 0!\n");
    forcx = 0; dparam_set("FORCX",0.0);
  }
  if (forcy < -7 || forcy > 7) {
    fprintf(stderr, "Acceptable values for FORCY are -7 to 7! Using 0!\n");
    forcy = 0; dparam_set("FORCY",0.0);
  }
 
#ifdef TEST
  setMap  (omega->mapx, omega->mapy);     /* use mapping() */
#else
  if (!ReadMap (omega)) {                  /* read in .map.rst file */
    update_forc (omega->mapy, dparam("AMPY"), dparam("FREQY"), dparam("PHIZY"), dparam("PHITY"), dparam("BETAY"), forcy);
  #ifndef FLOW_CONTROL
    update_forc (omega->mapx, dparam("AMPX"), dparam("FREQX"), dparam("PHIZX"), dparam("PHITX"), dparam("BETAX"), forcx);
  #else
    update_forc (omega->mapx, dparam("AMPX"), dparam("FREQX"), dparam("PHIZX"), dparam("PHITX"), dparam("BETAX"), forcx,omega->hilbert_coeff,
                 omega->mapy->d, omega->mapy->t);
  #endif
	
  }
#endif
  
  gradz_map (omega->mapx);                /* compute z-derivatives */
  gradz_map (omega->mapy);
  
  dparam_set ("ORIGINX", omega->mapx->d[0]);
  dparam_set ("ORIGINY", omega->mapy->d[0]);

  #ifdef FLOW_CONTROL
   current_amp = next_amp = omega->mapx->d[0];
  #endif

  return;
}

/* ------------------------------------------------------------------------ *
 * allocate_map() - memory for map
 * ------------------------------------------------------------------------ */

Mapping *allocate_map (int NZ)
{
  Mapping *newmap    = (Mapping *) calloc (1, sizeof(Mapping));
  
  newmap->NZ   = NZ;                                     /* number of z-plns */
  newmap->time = 0.0;                                    /* time             */

  if (option("dealias")) NZ = 3*NZ/2; // make space for dealiasing

  newmap->d    = (double *) calloc (NZ, sizeof(double)); /* displacement     */
  newmap->z    = (double *) calloc (NZ, sizeof(double)); /* z   - derivative */
  newmap->zz   = (double *) calloc (NZ, sizeof(double)); /* zz  - derivative */
  newmap->t    = (double *) calloc (NZ, sizeof(double)); /* t   - derivative */
  newmap->tt   = (double *) calloc (NZ, sizeof(double)); /* tt  - derivative */
  newmap->tz   = (double *) calloc (NZ, sizeof(double)); /* tz  - derivative */
  newmap->tzz  = (double *) calloc (NZ, sizeof(double)); /* tzz - derivative */
  newmap->f    = (double *) calloc (NZ, sizeof(double)); /* forcing          */
  
  return newmap;
}

/* ------------------------------------------------------------------------ *
 * 
 * UPDATING 
 * ========
 * 
 * UpdateMap()
 * writedog()
 * writedog2()
 * writecab()
 * setMap() (for TEST)
 * update_forc()
 * update_free()
 * newmark()
 * MapField()
 * 
 * ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ *
 * UpdateMap() -- update map using (1) mapping (2) forced (3) free 
 * ------------------------------------------------------------------------ */

static FILE  *Fp, *Fp2;
static double *basep;
static int step = 0, nz, nztot, mstep;

void UpdateMap (Domain *omega)
{
  Mapping  *mapx = omega->mapx,
       *mapy = omega->mapy;
  static double ampx, freqx, phizx, phitx, ampy, freqy, phizy, phity,
                betax, betay, wn, wnc, wnb, mass_, zeta_, dt, ramp, zramp;
  static double wnx,wny; //for the case of different IL and CF natural frequencies 
  double        mass, zeta;

  if (step == 0) {                               /* first time through */
    char fname[FILENAME_MAX];
    ampx = dparam("AMPX"); freqx = dparam("FREQX"); phizx = dparam("PHIZX"); phitx = dparam("PHITX"); betax = dparam("BETAX");
    ampy = dparam("AMPY"); freqy = dparam("FREQY"); phizy = dparam("PHIZY"); phity = dparam("PHITY"); betay = dparam("BETAY");

    wn   = dparam("WN");   wnc   = dparam("WNC");   wnb   = dparam("WNB");   mass_ = dparam("ZMASS");
    zeta_= dparam("ZETA"); dt    = dparam("DELT");  ramp  = dparam("RAMP");  zramp = dparam("ZRAMP");

    wnx  = dparam("WNX");  wny  = dparam("WNY");

    nz = option("NZ"); nztot = option("NZTOT"); mstep = max(iparam("MSTEP"),1);
    basep = dvector (0, nztot-1);
    ROOT {
      sprintf (fname, "%s.dog", omega->name);
      if ((Fp = fopen (fname,"w")) == (FILE *) NULL)
	ErrorHandler ("UpdateMap", "failed to open the force file", ERROR);
      writeheader (Fp, omega->name, LGmax, dparam("LZ"), nztot,
		   omega->U->nel, step, omega->mapx->time, "txyDLpk");
      sprintf (fname, "%s.cab", omega->name);
      if ((Fp2 = fopen (fname,"w")) == (FILE *) NULL)
	ErrorHandler ("UpdateMap", "failed to open the cable energy file", ERROR);
      writeheader (Fp2, omega->name, LGmax, dparam("LZ"), nztot,
		   omega->U->nel, step, omega->mapx->time, "t Ex Ey");
    }
    

     oldd_1 = (double *) calloc (nztot+1, sizeof(double));
     oldt_1 = (double *) calloc (nztot+1, sizeof(double));
     oldf_1 = (double *) calloc (nztot+1, sizeof(double));
    
   if( (forcy == -3) || (forcx == -3) )
    {
     generate_structure_matrix(wn, zeta_, omega->tens, omega->stif);
    }
   
    if( (forcx == -1) && (fabs(wnx)<1e-9) ) wnx = wn;
    if( (forcy == -1) && (fabs(wny)<1e-9) ) wny = wn;
  }
  
  step++;
  mapx->time += dt;
  mapy->time += dt;

#ifdef TMAP /* Do this only if we have time-dependent mapping */
#ifdef TEST
  setMap  (mapx, mapy);                           /* use mapping()        */
#else

  mass = mass_;                                   /* variable mass/zeta at startup */
  zeta = zeta_;
  if (step*dt < ramp) {
    mass = mass_ / (0.5 * (1.0 - cos(M_PI*step*dt/ramp)) + 1e-6);
    ROOT printf ("ramp: mass = %g\n", mass);
  }
  if (step*dt < zramp) {
    double expon = exp(-7.0*step*dt/zramp);
    zeta = 1.0 * expon + zeta_ * (1.0 - expon);
    ROOT printf ("ramp: zeta = %g\n", zeta);
  }

  if (forcy > 0)
    update_forc (mapy, ampy, freqy, phizy, phity, betay, forcy);
  else {
//    filter      (mapy, 'y');
  if(option("oldupdate"))
   update_free (mapy, dt, mass, wny, wnc, wnb, zeta, forcy);
  else
    update_free (mapy, dt, mass, wny, wnc, wnb, zeta, forcy, omega->tens, omega->stif);
  }

  /* update motion */
  if (forcx > 0) 
 #ifndef FLOW_CONTROL
    update_forc (mapx, ampx, freqx, phizx, phitx, betax, forcx); 
 #else
    update_forc (omega->mapx, dparam("AMPX"), dparam("FREQX"), dparam("PHIZX"), dparam("PHITX"), dparam("BETAX"), forcx,omega->hilbert_coeff,
                 omega->mapy->d, omega->mapy->t);
 #endif
  else {
//    filter      (mapx, 'x');
  if(option("oldupdate"))
    update_free (mapx, dt, mass, wnx, wnc, wnb, zeta, forcx);
  else
    update_free (mapx, dt, mass, wnx, wnc, wnb, zeta, forcx, omega->tens, omega->stif);	
  }
#endif /* ifdef TEST */

 if(!dparam("FIRST_ORDER"))
  {
    gradz_map (mapx);                  /* compute z-derivatives */
    gradz_map (mapy);
  }
 
  dparam_set ("ORIGINX", mapx->d[0]);
  dparam_set ("ORIGINY", mapy->d[0]);

#endif /* ifdef TMAP */
  return;
}

/* ------------------------------------------------------------------------ *
 * run_obj_vars() -- Runtime cylinder/cable/beam motion variables
 * ------------------------------------------------------------------------ */

void run_obj_vars(Domain *omega)
{
  if (step % mstep == 0) {  /* gather base pressure: */
    gather_baseP (omega, basep);

    ROOT 
      {
//	 writedog  (Fp,  omega->mapx, omega->mapy);
	writedog2 (Fp,  omega->mapx, omega->mapy, basep);
	writecab  (Fp2, omega->mapx, omega->mapy);  

      }
  }
  return;
}

/* ------------------------------------------------------------------------ *
 * writedog() - write displacement and forces to session.dog
 * ------------------------------------------------------------------------ */

void writedog (FILE *fp, Mapping *mapx, Mapping *mapy)
{
  int i, NZ = mapx->NZ;
  
#ifdef OLDFFTS
  realft (NZ/2, mapx->d,  1);              /* transform to Physical space  */
  realft (NZ/2, mapx->f,  1);
  realft (NZ/2, mapy->d,  1);
  realft (NZ/2, mapy->f,  1);
#else
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->f, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->f, 1, 0, 0, 0, 0);
#endif
  for (i = 0; i < NZ; i++)
    fprintf (fp, "%#10.4f %#10.4f %#10.4f %#10.4f %#10.4f %#2d\n",
	     mapx->time, mapx->d[i], mapy->d[i], mapx->f[i], mapy->f[i], i);
  fflush (fp);
#ifdef OLDFFTS
  realft (NZ/2, mapx->d, -1);              /* transform to Fourier space  */
  realft (NZ/2, mapx->f, -1);
  realft (NZ/2, mapy->d, -1);
  realft (NZ/2, mapy->f, -1);
#else
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->f, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->f, 1, 0, 0, 0, 0);
#endif

  return;
}
/* ------------------------------------------------------------------------ *
 * writedog2() - write displacement and forces and base pressure to session.dog
 * ------------------------------------------------------------------------ */

void writedog2 (FILE *fp, Mapping *mapx, Mapping *mapy, double *basep)
{
  int i, NZ = mapx->NZ;
  
#ifdef OLDFFTS
  realft (NZ/2, mapx->d,  1);              /* transform to Physical space  */
  realft (NZ/2, mapx->f,  1);
  realft (NZ/2, mapy->d,  1);
  realft (NZ/2, mapy->f,  1);
#else
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->tt, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->tt, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->f, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->f, 1, 0, 0, 0, 0);
#endif
  for (i = 0; i < NZ; i++)
//    fprintf (fp, "%#8.3f %#8.4f %#8.4f %#8.4f %#8.4f %#8.4f %#2d\n",
	fprintf (fp, "%#10.5f  %#10.6f %#10.6f  %#10.6f %#10.6f  %#10.6f %#10.6f  %#10.6f %#10.6f  %#10.6f %#10.6f %#10.6f  %#2d\n",
	     mapx->time, mapx->d[i], mapy->d[i],  mapx->z[i], mapy->z[i], mapx->t[i], mapy->t[i], 
       mapx->tt[i], mapy->tt[i], mapx->f[i], mapy->f[i],  basep[i], i);
//	fprintf (fp, "%#10.5f  %#10.6f %#10.6f  %#10.6f %#10.6f  %#10.6f %#10.6f  %#10.6f  %#2d\n",
//	     mapx->time, mapx->d[i], mapy->d[i],  mapx->z[i], mapy->z[i],
//       mapx->f[i], mapy->f[i],  basep[i], i);
  fflush (fp);
#ifdef OLDFFTS
  realft (NZ/2, mapx->d, -1);              /* transform to Fourier space  */
  realft (NZ/2, mapx->f, -1);
  realft (NZ/2, mapy->d, -1);
  realft (NZ/2, mapy->f, -1);
#else
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->tt, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->tt, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->f, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->f, 1, 0, 0, 0, 0);
#endif

  return;
}
void writecab (FILE *fp, Mapping *mapx, Mapping *mapy)
{
  int i, NZ = mapx->NZ;
  double wn   = dparam("WN"),
         wnc  = dparam("WNC"),
         wnb  = dparam("WNB"),
	 mass = dparam("ZMASS"),
         Ex   = 0.0,
         Ey   = 0.0;

  /* for free cables, may not want to include mean x displacement */

  // This should be do-able in Fourier Space.
#ifdef OLDFFTS
  realft (NZ/2, mapx->d,  1);              /* transform to Physical space  */
  realft (NZ/2, mapx->t,  1);
  realft (NZ/2, mapx->z,  1);
  realft (NZ/2, mapy->d,  1);
  realft (NZ/2, mapy->t,  1);
  realft (NZ/2, mapy->z,  1);
#else
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(rplan_inv, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
#endif

  for (i = 0; i < NZ; i++) {
#ifndef ENERGYX
    Ex += mapx->t[i]*mapx->t[i] + (wnc*wnc + wnb*wnb) * mapx->z[i]*mapx->z[i] + wn*wn * mapx->d[i]*mapx->d[i];
#else
    Ex += mapx->t[i]*mapx->t[i] + (wnc*wnc + wnb*wnb) * mapx->z[i]*mapx->z[i];
#endif /* ifndef ENERGYX */
#ifndef ENERGYY
    Ey += mapy->t[i]*mapy->t[i] + (wnc*wnc + wnb*wnb) * mapy->z[i]*mapy->z[i] + wn*wn * mapy->d[i]*mapy->d[i];
#else
    Ey += mapy->t[i]*mapy->t[i] + (wnc*wnc + wnb*wnb) * mapy->z[i]*mapy->z[i];
#endif /* ifndef ENERGYY */
  }

  Ex *= 0.5*mass/NZ;
  Ey *= 0.5*mass/NZ;

  fprintf (fp, "%#8.3f %#8.4f %#8.4f\n", mapx->time, Ex, Ey);
  fflush  (fp);

#ifdef OLDFFTS
  realft (NZ/2, mapx->d, -1);              /* transform to Fourier space  */
  realft (NZ/2, mapx->t, -1);
  realft (NZ/2, mapx->z, -1);
  realft (NZ/2, mapy->d, -1);
  realft (NZ/2, mapy->t, -1);
  realft (NZ/2, mapy->z, -1);
#else
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->t, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->t, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->z, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->z, 1, 0, 0, 0, 0);
#endif
  return;
}


/* ------------------------------------------------------------------------ *
 * writecab() - write x and y cable energy to session.cab
 * ------------------------------------------------------------------------ */

#ifdef TEST
/* ------------------------------------------------------------------------ *
 * setMap() -- update the map using mapping()
 * ------------------------------------------------------------------------ */

void setMap (Mapping *mapx, Mapping *mapy)
{
  int     NZ = mapx->NZ;
  register int k;
  double  *z = Zmesh (NZ);
  for (k = 0; k < NZ; k++)
    mapping (z[k], mapx->time,
	     mapx->d  + k, mapy->d  + k,
	     mapx->t  + k, mapy->t  + k,
	     mapx->tt + k, mapy->tt + k);
  free(z);

#ifdef OLDFFTS
  realft (NZ/2, mapx->d,  -1);              /* transform to Fourier space  */
  realft (NZ/2, mapx->t,  -1);
  realft (NZ/2, mapx->tt, -1);
  realft (NZ/2, mapy->d,  -1);
  realft (NZ/2, mapy->t,  -1);
  realft (NZ/2, mapy->tt, -1);
#else
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->d,  1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->d,  1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->t,  1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->t,  1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapx->tt, 1, 0, 0, 0, 0);
  rfftw(rplan, 1, (FFTW_COMPLEX *) mapy->tt, 1, 0, 0, 0, 0);
#endif
  return;
}
#endif

/* ------------------------------------------------------------------------ *
 * gradz_map() -- compute Mz, Mzz, Mtz, Mtzz
 * ------------------------------------------------------------------------ */

void gradz_map (Mapping *map)
{
  const int    NZ   = map->NZ;
  double   beta = 2.*M_PI/dparam("LZ"),
           fac, save;
  register int m;
  
#ifndef DOESNT_WORK
  dcopy  (NZ,   map->d,   1, map->z,   1);
  dcopy  (NZ,   map->z,   1, map->zz,  1);

  dcopy  (NZ,   map->t,   1, map->tz,  1);
  dcopy  (NZ,   map->tz,  1, map->tzz, 1);

  for (m = 0; m < NZ; m += 2) {
    fac = beta * (m / 2);

    save          =            map->z[m];
    map->z[m]     = -fac     * map->z[m+1];
    map->z[m+1]   =  fac     * save;              /* map->z   */

    map->zz[m]    = -fac*fac * map->zz[m];
    map->zz[m+1]  = -fac*fac * map->zz[m+1];      /* map->zz  */

    save          =            map->tz[m];
    map->tz[m]    = -fac     * map->tz[m+1];
    map->tz[m+1]  =  fac     * save;              /* map->tz  */

    map->tzz[m]   = -fac*fac * map->tzz[m];
    map->tzz[m+1] = -fac*fac * map->tzz[m+1];     /* map->tzz */
  }
  
//  int n = iparam ("NFILTER");
//   
//   if(n>0)  
//   for (int i = NZ-n; i < NZ; i++) {
//      map->z[i]   = 0.;
//      map->zz[i]  = 0.;
//      map->tz[i]  = 0.;
//      map->tzz[i] = 0.;
//     }
#else 
/* ---------------------------------------------------------------------*
 * Unless we can come up with a DCT which takes a input of n datapoints
 * this is not a viable option - both the ESSL as well as the Numerical 
 * recipes (v.2) DCT require n+1 datapoints - the v.1 Numerical recipes
 * routine doesn't BUT gived wrong results - so does the FFTPACK/BIHAR
 * routine that is also very slow for n-1 prime wtc. 
 * ---------------------------------------------------------------------*/
  dcopy (NZ,   map->d,   1, map->z,   1);
  sinft (NZ,   map->z,  -1);                     /* to fourier */
  dcopy (NZ,   map->z,   1, map->zz,  1);

  dcopy (NZ,   map->t,   1, map->tz,  1);
  sinft (NZ,   map->tz, -1);                     /* to fourier */
  dcopy (NZ,   map->tz,  1, map->tzz, 1);

  beta *= 0.5;                 /* halve beta for this expansion */
  for (m = 0; m < NZ; m++) {
    fac = beta * m;

    map->z[m]     =  fac     * map->z[m];         /* map->z   */

    map->zz[m]    = -fac*fac * map->zz[m];        /* map->zz  */

    map->tz[m]    =  fac     * map->tz[m];        /* map->tz  */

    map->tzz[m]   = -fac*fac * map->tzz[m];       /* map->tzz */
  }
 
  
  cosft (NZ, map->z,   1);                        /* to phys  */
  sinft (NZ, map->zz,  1);
  cosft (NZ, map->tz,  1);
  sinft (NZ, map->tzz, 1);
#endif /* ifdef FIXED */
  return;
}

#endif //end of map 

/* ------------------------------------------------------------------------ *
 * update_forc() -- update the map using amp freq phi
 * ------------------------------------------------------------------------ */

void update_forc (Mapping *map, double amp, double freq, double phiz, double phit,
		  double BETA, int forc)
{
  const int NZ = map->NZ;
  const double beta = (BETA > 0.0 ? BETA : 2.*M_PI/dparam("LZ")),
               timec = dparam("TIMEC");
  double *z = Zmesh (NZ);

  register int k;
  double t, CoS, SiN, denom, fac, dfac, ddfac, startime;
  register double CoSZ, loop_phase;

  double  init_position = 0.;

  if (timec != 0.0) {
    if (iparam("IGRSTART")) /* restarting needs .rea file info */
      startime = dparam("STARTING");
    else
      startime = dparam("STARTIME");
    t = map->time - startime;
    denom = (t*t + timec);
    fac = t*t/denom;
    dfac = 2.*t*timec/(denom*denom);
    ddfac = -2.*timec*(3.*t*t - timec)/(denom*denom*denom);
  } else
    t = map->time;

  if (forc < 0) forc = 0; /* Handle restart properly for these cases */
  switch (forc) {
  case 0: /* Only gets called if structure is uninitialized, sets it to zero */
    dzero(NZ,map->d,1);
    dzero(NZ,map->t,1);
    dzero(NZ,map->tt,1);
    break;    
  case 1: /* standing wave */
    if (timec != 0.0) {
      CoS = amp * cos(freq*t+phit);
      SiN = amp * sin(freq*t+phit);
      for (k = 0; k < NZ; k++) {
	CoSZ = cos(beta*z[k]+phiz);
	map->d[k]  = CoSZ * (fac * CoS);
	map->t[k]  = CoSZ * (fac * (-freq * SiN) + dfac * CoS);
	map->tt[k] = CoSZ * (fac * (-freq*freq * CoS) + ddfac * CoS
			     -2. * dfac * freq * SiN);
      }
    } else {
      CoS = amp * cos(freq*t+phit);
      SiN = amp * sin(freq*t+phit);
      for (k = 0; k < NZ; k++) {
	CoSZ = cos(beta*z[k]+phiz);
	map->d[k]  =              CoSZ * CoS;
	map->t[k]  = -freq      * CoSZ * SiN;
	map->tt[k] = -freq*freq * CoSZ * CoS;
      }
    }
#ifdef MAP
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->d,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->t,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->tt, 1, 0, 0, 0, 0);
#endif
    break;
  case 2: /* traveling wave */
    if (timec != 0.0) {
      for (k = 0; k < NZ; k++) {
	loop_phase = beta*z[k]+phiz + freq*t+phit;
	CoS = amp * cos(loop_phase);
	SiN = amp * sin(loop_phase);
	map->d[k]  = fac * CoS;
	map->t[k]  = fac * (-freq * SiN) + dfac * CoS;
	map->tt[k] = fac * (-freq*freq * CoS) + ddfac * CoS
	  -2. * dfac * freq * SiN;
      }
    } else {
      for (k = 0; k < NZ; k++) {
	loop_phase = beta*z[k]+phiz + freq*t+phit;
	CoS = amp * cos(loop_phase);
	SiN = amp * sin(loop_phase);
	map->d[k]  =              CoS;
	map->t[k]  = -freq      * SiN;
	map->tt[k] = -freq*freq * CoS;
      }
    }
#ifdef MAP
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->d,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->t,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->tt, 1, 0, 0, 0, 0);
#endif
    break;
  case 3: /* rigid motion */
   
#ifdef SPM
    init_position = map->d[0];
#endif
    if (timec != 0.0) {
      CoS = amp * cos(freq*t+phit);
      SiN = amp * sin(freq*t+phit);
      denom = (t*t + timec);
      fac = t*t/denom;
      dfac = 2.*t*timec/(denom*denom);
      ddfac = -2.*timec*(3.*t*t - timec)/(denom*denom*denom);
      map->d[0]  = fac * CoS;
      map->t[0]  = fac * (-freq * SiN) + dfac * CoS;
      map->tt[0] = fac * (-freq*freq * CoS) + ddfac * CoS
	-2. * dfac * freq * SiN;
    } else {
      CoS = amp * cos(freq*t+phit);
      SiN = amp * sin(freq*t+phit);
      map->d[0]  =  CoS;
      map->t[0]  = -freq      * SiN;
      map->tt[0] = -freq*freq * CoS;
    }

#ifdef MAP
    dzero(NZ-1,map->d+1,1);
    dzero(NZ-1,map->t+1,1);
    dzero(NZ-1,map->tt+1,1);
#endif

#if !defined(MAP) && defined(SPM)
    for(k=1;k<NZ; ++k) //copy all map stuffs from zero plane
     {
      map->d[k]  = map->d[0];
      map->t[k]  = map->t[0];
      map->tt[k] = map->tt[0];
     }
#endif

    break;
  case 4: /* wavy pipe but not move */
    if (timec != 0.0) {
      CoS = amp * cos(freq*t+phit);
      SiN = amp * sin(freq*t+phit);
      for (k = 0; k < NZ; k++) {
	CoSZ = cos(beta*z[k]+phiz);
	map->d[k]  = CoSZ * (fac * CoS);
	map->t[k]  = CoSZ * (fac * (-freq * SiN) + dfac * CoS);
	map->tt[k] = CoSZ * (fac * (-freq*freq * CoS) + ddfac * CoS
			     -2. * dfac * freq * SiN);
      }
    } else {
      for (k = 0; k < NZ; k++) {
	CoSZ = amp*cos(beta*z[k]+phiz);
	map->d[k]  =              CoSZ;
	map->t[k]  = 0.;
	map->tt[k] = 0.;
      }
    }
 #ifdef MAP
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->d,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->t,  1, 0, 0, 0, 0);
    rfftw(rplan, 1, (FFTW_COMPLEX *) map->tt, 1, 0, 0, 0, 0);
  #endif
    break;
  default:
    fprintf(stderr,"Invalid case: Coding error!\n");
    exit (-1);
    break;
  }

  free(z);

  return;
}

#ifdef FLOW_CONTROL
void update_forc (Mapping *map, double amp, double freq, double phiz, double phit,
		  double BETA, int forc,  double *hilbert_coeff, double *dy, double *ty)
{
  int mstep  = iparam("MSTEP");
  int ncoeff = iparam("NHF_COEFF");
  double dt = dparam("DELT");

  double Kp = dparam("DKP");
  double phi_d = dparam("DPHI_D");
  double wny  = dparam("WNY");
  //double delta_phi = dparam("DELTA_PHI");//2.0*wny*dt*mstep; 
  double delta_phi = 2.0*wny*dt*mstep; 
  double max_disp  = dparam("DISP_CONTROL")? dparam("DISP_CONTROL"):0.3 ;

  const int NZ = map->NZ;
  const double beta = (BETA > 0.0 ? BETA : 2.*M_PI/dparam("LZ")),
               timec = dparam("TIMEC");
  double *z = Zmesh (NZ);

  register int k;
  double t, CoS, SiN, denom, fac, dfac, ddfac, startime;
  register double CoSZ, loop_phase;

  double  init_position = 0.;

  if (timec != 0.0) {
    if (iparam("IGRSTART")) /* restarting needs .rea file info */
      startime = dparam("STARTING");
    else
      startime = dparam("STARTIME");
    t = map->time - startime;
    denom = (t*t + timec);
    fac = t*t/denom;
    dfac = 2.*t*timec/(denom*denom);
    ddfac = -2.*timec*(3.*t*t - timec)/(denom*denom*denom);
  } else
    t = map->time;

  if (forc < 0) forc = 0; /* Handle restart properly for these cases */
  switch (forc) {
  case 0: /* Only gets called if structure is uninitialized, sets it to zero */
    dzero(NZ,map->d,1);
    dzero(NZ,map->t,1);
    dzero(NZ,map->tt,1);
    break;    
  case 1: /* rigid motion but under P control*/
    {

      double old_d  = map->d[0];
      double old_t  = map->t[0];
      //update action
      if (step % mstep == 0)
       {
        for(int i=0; i<ncoeff-1; ++i)
         {
          yd[i] = yd[i+1];
          td[i] = td[i+1];
         }
          yd[ncoeff-1] = dy[0];
          td[ncoeff-1] = ty[0];
       
       if(step >= mstep*ncoeff)
        {
         current_amp = next_amp;
        //next_amp = HilbertNextAmp (yd, td, hilbert_coeff);
         double Hzt = 0.0;
        for(int i=0; i<ncoeff; ++i)
          Hzt += hilbert_coeff[i]*yd[i]*td[i];
          Hzt = -Hzt;

         double Hzt_sign =1.0;
         if(Hzt<0.0) Hzt_sign = -1.0;
         double ZHzt = Hzt_sign*sqrt(2.0*fabs(Hzt)/wny);


         double zt = dy[0]*ty[0];
         double zt_sign =1.0;
         if(zt<0.0) zt_sign = -1.0;
         double Zt = zt_sign*sqrt(2.0*fabs(zt)/wny);
          
          next_amp = Kp*(cos(phi_d-0.5*M_PI+delta_phi)*Zt+sin(phi_d-0.5*M_PI+delta_phi)*ZHzt);

        }
        
       }

      double amp_disp = next_amp-current_amp;
      double sign_disp = 1.0;
      if(amp_disp<0.0) sign_disp=-1.0;
      if(fabs(amp_disp)>max_disp)
      {
        amp_disp = sign_disp*max_disp;
        next_amp = current_amp+amp_disp;
      }

      double intermediate_amp = current_amp +double(step% mstep)/mstep*(next_amp-current_amp);

     ROOTONLY  fprintf(stderr,"current = %lf  next =  %lf  intermediate =  %lf \n",
      current_amp, next_amp, intermediate_amp );

      map->d[0]  = intermediate_amp;
      map->t[0]  = (map->d[0]-old_d)/2.0/dt; //firdst order
      map->tt[0] = (map->t[0]-old_t)/2.0/dt; //first order

#ifdef MAP
    dzero(NZ-1,map->d+1,1);
    dzero(NZ-1,map->t+1,1);
    dzero(NZ-1,map->tt+1,1);
#endif

#if !defined(MAP) && defined(SPM)
    for(k=1;k<NZ; ++k) //copy all map stuffs from zero plane
     {
      map->d[k]  = map->d[0];
      map->t[k]  = map->t[0];
      map->tt[k] = map->tt[0];
     }
#endif
 
    break;
  }
  default:
    fprintf(stderr,"Invalid case: Coding error!\n");
    exit (-1);
    break;
  }

  free(z);

  return;
}
#endif

/* ------------------------------------------------------------------------ *
 * update_free() -- update the map using dt, mass, wn and zeta
 * ------------------------------------------------------------------------ */

void update_free (Mapping *map, double dt, double mass, double wn, double wnc, double wnb, double zeta, int forc)
{
  
  int      NZ   = map->NZ;
  int NZM = iparam("MNL")  ;
 // ROOTONLY fprintf(stderr,"NZM=%d", NZM);
  double   beta = 2.*M_PI/dparam("LZ"),
           lambda, zeta_m, frq, keep;
  double   dz = dparam("LZ")/option("NZTOT"),
	   AM = dparam("AMOR"),k10,k20,k30,k40,
	   savd,savt,savdz0,savdzN,keepmapd, keepmapt;

  double addmass = dparam("DADDMASS");
  double zetan   = dparam("DZETAN");
  double frhs = 0.;
  double massa = mass+addmass;

#ifdef AFMASS
  double af_residual = 0.0;
  double af_mass = 0.0;
  double af_alpha = dparam("DAF_ALPHA");
#endif
  double   KSI = 500;
  double *oldd = (double *) calloc (NZ+1, sizeof(double));
  double *oldt = (double *) calloc (NZ+1, sizeof(double));
  double *oldtt = (double *) calloc (NZ+1, sizeof(double));
  double *olddz = (double *) calloc (NZ+1, sizeof(double));
  double *olddzz = (double *) calloc (NZ+1, sizeof(double));
  double *olddzzz = (double *) calloc (NZ+1, sizeof(double));
  double *olddzzzz = (double *) calloc (NZ+1, sizeof(double));
  double *subd = (double *) calloc (NZ, sizeof(double));
  double *subt = (double *) calloc (NZ, sizeof(double));
  double *subtt = (double *) calloc (NZ, sizeof(double));
  double *k1d = (double *) calloc (NZ, sizeof(double));
  double *k2d = (double *) calloc (NZ, sizeof(double));
  double *k3d = (double *) calloc (NZ, sizeof(double));
  double *k4d = (double *) calloc (NZ, sizeof(double));
  double *k1t = (double *) calloc (NZ, sizeof(double));
  double *k2t = (double *) calloc (NZ, sizeof(double));
  double *k3t = (double *) calloc (NZ, sizeof(double));
  double *k4t = (double *) calloc (NZ, sizeof(double)); 

//input for external force
  double   amp = dparam("AMM_Force");
  double   freq = dparam("freq_Force");
  double   phit = dparam("phit_Force"); 
  double   phiz = dparam("phiz_Force"); 

  register int m;
  
  const double   timec = dparam("TIMEC");
  double *z = Zmesh (NZ);
  register int k;
  double t, CoS, SiN, denom, fac, dfac, ddfac, startime;
  register double CoSZ, loop_phase;

  if (timec != 0.0) {
    if (iparam("IGRSTART")) /* restarting needs .rea file info */
      startime = dparam("STARTING");
    else
      startime = dparam("STARTIME");
    t = map->time - startime;
    denom = (t*t + timec);
    fac = t*t/denom;
    dfac = 2.*t*timec/(denom*denom);
    ddfac = -2.*timec*(3.*t*t - timec)/(denom*denom*denom);
  } else
    t = map->time;
	

  switch (forc) {
  case 0: /* periodic in all derivatives, only constrained by virtual sping */
    if (option("nomean")) { /* Allow no motion in the mean */
      /* We explicitly force all of the variables below to be zero */
      map->d[0] = map->t[0] = map->tt[0] = 0.;
      map->d[1] = map->t[1] = map->tt[1] = 0.;
    } else {
      /* We explicitly force all of the variables below to be zero */
      lambda = wn;
      zeta_m = zeta;
      newmark (map->d, map->t, map->tt, map->f[0], dt, mass, lambda, zeta_m);
      map->d[1] = map->t[1] = map->tt[1] = 0.;
    }
    for (m = 2; m < NZ; m++) {
      frq    = beta * (m/2);
      lambda = frq*frq * (wnc*wnc + frq*frq  * wnb*wnb);
      lambda = sqrt (wn*wn + lambda);
      zeta_m = (lambda*zeta*wn != 0.0 ? zeta * wn/lambda : 0.0);
      newmark (map->d+m, map->t+m, map->tt+m,
	       map->f[m], dt, mass, lambda, zeta_m);
    }
    break;
  case -1: /* Allowed to freely move rigidly */
    lambda = wn;
#ifdef AFMASS
    af_residual = map->f[0]/mass-wn*wn*map->d[0]-2.0*zeta*wn*map->t[0]-map->tt[0];
    af_mass = af_alpha*fabs(af_residual)/fabs(map->d[0])*dt*dt;
    ROOTONLY fprintf(stderr,"artificial mass: %g  %g  %g  \n",af_alpha,af_residual,af_mass);
#endif
//added mass zwang
    lambda = lambda*sqrt(mass/massa);
    frhs = map->f[0]+addmass*map->tt[0]+(2.*zetan*wn)*map->t[0];
#ifdef AFMASS
    frhs -= af_mass*map->tt[0];
#endif
//    zeta_m = zeta;
    zeta_m = zeta+zetan;
//    newmark (map->d, map->t, map->tt, map->f[0], dt, mass, lambda, zeta_m);
    newmark (map->d, map->t, map->tt, frhs, dt, massa, lambda, zeta_m);

    dzero(NZ-1,map->d+1,1);
    dzero(NZ-1,map->t+1,1);
    dzero(NZ-1,map->tt+1,1);
    break;
  case -2: /* Allowed to move but with endpoints fixed to zero */

    fprintf(stderr, "there is some problems in this option....please use forcx =-3 and without -ou \n");
    exit(-1);

    realft (NZ/2, map->d,  1);              /* transform to Physical space  */
    realft (NZ/2, map->t,  1);
    realft (NZ/2, map->tt, 1);
    map->d[0] = map->t[0] = map->tt[0] = 0.0;    /* make sure this is true  */
//    sinft (NZ, map->d,  -1);                 /* transform to Fourier space  */
//    sinft (NZ, map->t,  -1);
//    sinft (NZ, map->tt, -1);
    realft (NZ/2, map->f, 1);        /* transform forcing to Physical space */
//	for (m = NZ/2+1; m < NZ; m++) {
//	map->f[m] =0;	
//		}
    keep = map->f[0];         /* Remember the value of the endpoint forcing */
    map->f[0]=0.0;                         /* zero out the endpoint forcing */
//    sinft (NZ, map->f,  -1);     /* sine transform forcing to Fourier space */
    beta *= 0.5;                           /* halve beta for this expansion */
  
    for (m = 1; m < NZ; m++) {
      frq = beta * m;
      lambda = frq*frq * (wnc*wnc + frq*frq  * wnb*wnb);
      lambda = sqrt (wn*wn + lambda);
//added mass zwang
      lambda = lambda*sqrt(mass/massa);

      zeta_m = (lambda*zeta*wn != 0.0 ? (zeta+zetan) * wn/lambda : 0.0);

      frhs = map->f[m]+addmass*map->tt[m]+(2.*zetan*wn)*map->t[m];
//      zeta_m = (lambda*zeta*wn != 0.0 ? zeta * wn/lambda : 0.0);
//      newmark (map->d+m, map->t+m, map->tt+m,
//	       map->f[m], dt, mass, lambda, zeta_m);
      newmark (map->d+m, map->t+m, map->tt+m,
	       frhs, dt, massa, lambda, zeta_m);
    }
#if SHOW_MAPD
    ROOT {
      for (m = 0; m < NZ; m++) 
	printf (" %G", map->d[m]);
      putchar ('\n');
    }
#endif
    sinft (NZ, map->d,  1);            /* sine transform to Physical space  */
#if SHOW_MAPD
    ROOT {
      for (m = 0; m < NZ; m++) 
	printf (" %G", map->d[m]);
      putchar ('\n');
    }
#endif
//    sinft (NZ, map->t,  1);            /* sine transform to Physical space  */
//    sinft (NZ, map->tt, 1);
//    
//    sinft (NZ, map->f,  1);     /* sine transform forcing to Physical space */

    realft (NZ/2, map->d, -1);                /* transform to Fourier space */
    realft (NZ/2, map->t, -1);
    realft (NZ/2, map->tt,-1);

    map->f[0] = keep;           /* Recall the value of the endpoint forcing */
    realft (NZ/2, map->f, -1);        /* transform forcing to Fourier space */
    break;
  case -3: /* periodic in all derivatives, with damping at the end */
 
    realft (NZ/2, map->d,  1);              /* transform to Physical space  */
    realft (NZ/2, map->t,  1);
    realft (NZ/2, map->tt, 1);
    realft (NZ/2, map->f, 1);       


    for (m = 0; m < NZ; m += 1) {
     subd[m]=map->d[m];
     subt[m]=map->t[m];
     subtt[m]=map->tt[m];
     oldd[m]=map->d[m];
     oldt[m]=map->t[m];
     oldtt[m]=map->tt[m];}

     oldd[NZ]=oldd[0];
     oldt[NZ]=oldt[0];
     oldtt[NZ]=oldtt[0];

  fprintf(stdout,"KSI= %g\n",KSI);

// Premiere etape //

    for (m = 1; m < NZ; m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz); 
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

    for (m = 1; m < NZ; m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}

      k10=wnc*wnc/AM*(olddz[0]-olddz[NZ])-KSI/AM*oldd[0];
      oldd[0]=subd[0]+k10*dt/2;
for (m = 1; m < NZ; m++) {
      k1d[m]=subt[m];
      k1t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m];
      oldd[m]=subd[m]+k1d[m]*dt/2;
      oldt[m]=subt[m]+k1t[m]*dt/2;}

    oldd[NZ]=oldd[0];

// Deuxieme etape //

    for (m = 1; m < NZ; m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz);
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

    for (m = 1; m < NZ; m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}

      k20=wnc*wnc/AM*(olddz[0]-olddz[NZ])-KSI/AM*oldd[0];
      oldd[0]=subd[0]+k20*dt/2;
for (m = 1; m < NZ; m++) {
      k2d[m]=subt[m]+k1t[m]*dt/2;
      k2t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m];
      oldd[m]=subd[m]+k2d[m]*dt/2;
      oldt[m]=subt[m]+k2t[m]*dt/2;}

    oldd[NZ]=oldd[0];

// Troisieme etape //

    for (m = 1; m < NZ; m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz);
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

    for (m = 1; m < NZ; m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}

      k30=wnc*wnc/AM*(olddz[0]-olddz[NZ])-KSI/AM*oldd[0];
      oldd[0]=subd[0]+k30*dt;
for (m = 1; m < NZ; m++) {
      k3d[m]=subt[m]+k2t[m]*dt/2;
      k3t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m];
      oldd[m]=subd[m]+k3d[m]*dt;
      oldt[m]=subt[m]+k3t[m]*dt;}

    oldd[NZ]=oldd[0];

// Quatrieme etape //

    for (m = 1; m < NZ; m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz);
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

    for (m = 1; m < NZ; m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}

      k40=wnc*wnc/AM*(olddz[0]-olddz[NZ])-KSI/AM*oldd[0];
      oldd[0]=subd[0]+(k10+2*k20+2*k30+k40)*dt/6;
for (m = 1; m < NZ; m++) {
      k4d[m]=subt[m]+k3t[m]*dt;
      k4t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m];
      oldd[m]=subd[m]+(k1d[m]+2*k2d[m]+2*k3d[m]+k4d[m])*dt/6;
      oldt[m]=subt[m]+(k1t[m]+2*k2t[m]+2*k3t[m]+k4t[m])*dt/6;}

    oldd[NZ]=oldd[0];

// Derniere etape pour l'acceleration  //

for (m = 1; m < NZ; m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz);
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

    for (m = 1; m < NZ; m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}

    for (m = 1; m < NZ; m++) {
     oldtt[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m];}
     
    oldt[0]=wnc*wnc/AM*(olddz[0]-olddz[NZ])-KSI/AM*oldd[0];
    oldtt[0]=wnc*wnc/(2*dz*AM)*(4*oldt[1]+4*oldt[NZ-1]-oldt[2]-oldt[NZ-2]-6*oldt[0])-KSI/AM*oldt[0];  

for (m = 0; m < NZ; m += 1) {
     map->d[m]=oldd[m];
     map->t[m]=oldt[m];
     map->tt[m]=oldtt[m];}

      realft (NZ/2, map->d, -1);                /* transform to Fourier space */
      realft (NZ/2, map->t, -1);
      realft (NZ/2, map->tt,-1);
      realft (NZ/2, map->f,-1);
      map->d[1] = map->t[1] = map->tt[1] = 0.;

    break;
 case -4: /* Allowed to move with prescibed map->f */
    realft (NZ/2, map->d,  1);              /* transform to Physical space  */
    realft (NZ/2, map->t,  1);
    realft (NZ/2, map->tt, 1);
    map->d[0] = map->t[0] = map->tt[0] = 0.0;    /* make sure this is true  */
    sinft (NZ, map->d,  -1);                 /* transform to Fourier space  */
    sinft (NZ, map->t,  -1);
    sinft (NZ, map->tt, -1);
	
   if (timec != 0.0) {
      CoS = amp * cos(freq*t+phit);
      for (k = 0; k < NZ; k++) {
	CoSZ = cos(beta*z[k]+phiz);
	map->f[k]  = CoSZ * (fac * CoS);
      }
    } else {
      CoS = amp * cos(freq*t+phit);
      for (k = 0; k < NZ; k++) {
	CoSZ = cos(beta*z[k]+phiz);
	map->f[k]  =              CoSZ * CoS;
      }
    }

    realft (NZ/2, map->f, 1);        /* transform forcing to Physical space */
    keep = map->f[0];         /* Remember the value of the endpoint forcing */
    map->f[0]=0.0;                         /* zero out the endpoint forcing */
    sinft (NZ, map->f,  -1);     /* sine transform forcing to Fourier space */
  
    for (m = 1; m < NZ; m++) {
      lambda = 0;
      zeta_m = 0;
      newmark (map->d+m, map->t+m, map->tt+m,
	       map->f[m], dt, mass, lambda, zeta_m);
    }
	
#if SHOW_MAPD
    ROOT {
      for (m = 0; m < NZ; m++) 
	printf (" %G", map->d[m]);
      putchar ('\n');
    }
#endif
    sinft (NZ, map->d,  1);            /* sine transform to Physical space  */
#if SHOW_MAPD
    ROOT {
      for (m = 0; m < NZ; m++) 
	printf (" %G", map->d[m]);
      putchar ('\n');
    }
#endif
    sinft (NZ, map->t,  1);            /* sine transform to Physical space  */
    sinft (NZ, map->tt, 1);
    
    sinft (NZ, map->f,  1);     /* sine transform forcing to Physical space */

    realft (NZ/2, map->d, -1);                /* transform to Fourier space */
    realft (NZ/2, map->t, -1);
    realft (NZ/2, map->tt,-1);

    map->f[0] = keep;           /* Recall the value of the endpoint forcing */
    realft (NZ/2, map->f, -1);        /* transform forcing to Fourier space */
    break;
	//////////////////////////////////////////////
	case -5:
	 realft (NZ/2, map->d,  1);           
    realft (NZ/2, map->t,  1);
    realft (NZ/2, map->tt, 1);
    realft (NZ/2, map->f, 1);    
	
	map->d[0] = map->t[0] = map->tt[0] = 0.0; 
     for (m = 0; m < (NZ); m += 1) 
	 {
     subd[m]=map->d[m];
     subt[m]=map->t[m];
     subtt[m]=map->tt[m];
     oldd[m]=map->d[m];
     oldt[m]=map->t[m];
     oldtt[m]=map->tt[m];
	 }

	 oldd[NZ]=oldd[0];
     oldt[NZ]=oldt[0];
     oldtt[NZ]=oldtt[0];

 //  ROOTONLY fprintf(stderr,"go to first step\n");

// Premiere etape //

     for (m = 1; m < (NZ); m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz); 
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

     for (m = 1; m < (NZ); m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}
		olddzz[0] = olddzz[NZ] = 0.0; //both endpoints of beam are pinned. 
     for (m = 01; m < (NZ); m += 1) {
     olddzzz[m]= (olddzz[m+1]-olddzz[m-1])/(2*dz);}
   olddzzz[0] = (-3*olddzz[0]+4*olddzz[1]-olddzz[2])/(2*dz); 
   olddzzz[NZ]= (3*olddzz[NZ]-4*olddzz[NZ-1]+olddzz[NZ-2])/(2*dz);
     for (m = 1; m < (NZ); m += 1) {
     olddzzzz[m]= (olddzzz[m+1]-olddzzz[m-1])/(2*dz);}
     for (m = 0; m < (NZM); m += 1) {
      k1d[m]=subt[m];
      oldd[m]=subd[m]+k1d[m]*dt/2;
	  k1t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
      oldt[m]=subt[m]+k1t[m]*dt/2;
	  }

	 for (m = NZM; m < (NZ); m += 1) {
      k1d[m]=subt[m];
      oldd[m]=subd[m]+k1d[m]*dt/2;
      k1t[m]=(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
	  oldt[m]=subt[m]+k1t[m]*dt/2;
	  }
	  
 //  ROOTONLY fprintf(stderr,"go to second step\n");
// Deuxieme etape //
    for (m = 1; m < (NZ); m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz); 
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

     for (m = 1; m < (NZ); m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}
		olddzz[0] = olddzz[NZ] = 0.0; //both endpoints of beam are pinned. 
     for (m = 01; m < (NZ); m += 1) {
     olddzzz[m]= (olddzz[m+1]-olddzz[m-1])/(2*dz);}
   olddzzz[0] = (-3*olddzz[0]+4*olddzz[1]-olddzz[2])/(2*dz); 
   olddzzz[NZ]= (3*olddzz[NZ]-4*olddzz[NZ-1]+olddzz[NZ-2])/(2*dz);
     for (m = 1; m < (NZ); m += 1) {
     olddzzzz[m]= (olddzzz[m+1]-olddzzz[m-1])/(2*dz);}
	 
	for (m = 1; m < (NZM); m++) {
      k2d[m]=subt[m]+k1t[m]*dt/2;
	k2t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
      oldd[m]=subd[m]+k2d[m]*dt/2;
      oldt[m]=subt[m]+k2t[m]*dt/2;}
	for (m = NZM; m < (NZ); m++) {
      k2d[m]=subt[m]+k1t[m]*dt/2;
	k2t[m]=(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
      oldd[m]=subd[m]+k2d[m]*dt/2;
      oldt[m]=subt[m]+k2t[m]*dt/2;}
  
 //  ROOTONLY fprintf(stderr,"go to third step\n");
// Troisieme etape //
    for (m = 1; m < (NZ); m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz); 
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

     for (m = 1; m < (NZ); m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}
		olddzz[0] = olddzz[NZ] = 0.0; //both endpoints of beam are pinned. 
     for (m = 1; m < (NZ); m += 1) {
     olddzzz[m]= (olddzz[m+1]-olddzz[m-1])/(2*dz);}
   olddzzz[0] = (-3*olddzz[0]+4*olddzz[1]-olddzz[2])/(2*dz); 
   olddzzz[NZ]= (3*olddzz[NZ]-4*olddzz[NZ-1]+olddzz[NZ-2])/(2*dz);
     for (m = 1; m < (NZ); m += 1) {
     olddzzzz[m]= (olddzzz[m+1]-olddzzz[m-1])/(2*dz);}
	for (m = 1; m < (NZM); m++) {
      k3d[m]=subt[m]+k2t[m]*dt/2;
      k3t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
        oldd[m]=subd[m]+k3d[m]*dt;
      oldt[m]=subt[m]+k3t[m]*dt;}
	for (m = NZM; m < (NZ); m++) {
      k3d[m]=subt[m]+k2t[m]*dt/2;
      k3t[m]=(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
        oldd[m]=subd[m]+k3d[m]*dt;
      oldt[m]=subt[m]+k3t[m]*dt;}

//   ROOTONLY fprintf(stderr,"go to fourth step\n");
// Quatrieme etape //
    for (m = 1; m < (NZ); m += 1) {
     olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz); 
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);

     for (m = 1; m < (NZ); m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}
		olddzz[0] = olddzz[NZ] = 0.0; //both endpoints of beam are pinned. 
     for (m = 01; m < (NZ); m += 1) {
     olddzzz[m]= (olddzz[m+1]-olddzz[m-1])/(2*dz);}
   olddzzz[0] = (-3*olddzz[0]+4*olddzz[1]-olddzz[2])/(2*dz); 
   olddzzz[NZ]= (3*olddzz[NZ]-4*olddzz[NZ-1]+olddzz[NZ-2])/(2*dz);
     for (m = 1; m < (NZ); m += 1) {
     olddzzzz[m]= (olddzzz[m+1]-olddzzz[m-1])/(2*dz);}
	for (m = 1; m < (NZM); m++) {
      k4d[m]=subt[m]+k3t[m]*dt;
      k4t[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
      oldd[m]=subd[m]+(k1d[m]+2*k2d[m]+2*k3d[m]+k4d[m])*dt/6;
      oldt[m]=subt[m]+(k1t[m]+2*k2t[m]+2*k3t[m]+k4t[m])*dt/6;}
	for (m = NZM; m < (NZ); m++) {
      k4d[m]=subt[m]+k3t[m]*dt;
      k4t[m]=(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
      oldd[m]=subd[m]+(k1d[m]+2*k2d[m]+2*k3d[m]+k4d[m])*dt/6;
      oldt[m]=subt[m]+(k1t[m]+2*k2t[m]+2*k3t[m]+k4t[m])*dt/6;}


//   ROOTONLY fprintf(stderr,"go to final step\n");
// Derniere etape pour l'acceleration  //

 for (m = 1; m < (NZ); m += 1) {
    olddz[m]=	(oldd[m+1]-oldd[m-1])/(2*dz);}

    olddz[0] = (-3*oldd[0]+4*oldd[1]-oldd[2])/(2*dz);
    olddz[NZ]= (3*oldd[NZ]-4*oldd[NZ-1]+oldd[NZ-2])/(2*dz);
     for (m = 1; m < (NZ); m += 1) {
     olddzz[m]=	(olddz[m+1]-olddz[m-1])/(2*dz);}
		olddzz[0] = olddzz[NZ] = 0.0; //both endpoints of beam are pinned. 	 
        for (m = 1; m < (NZ); m += 1) {
     olddzzz[m]= (olddzz[m+1]-olddzz[m-1])/(2*dz);}
   olddzzz[0] = (-3*olddzz[0]+4*olddzz[1]-olddzz[2])/(2*dz); 
   olddzzz[NZ]= (3*olddzz[NZ]-4*olddzz[NZ-1]+olddzz[NZ-2])/(2*dz);
        for (m = 1; m < (NZ); m += 1) {
     olddzzzz[m]= (olddzzz[m+1]-olddzzz[m-1])/(2*dz);}
    for (m = 1; m < (NZM); m++) {
     oldtt[m]=map->f[m]/mass+(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
    }  
    for (m = NZM; m < (NZ); m++) {
     oldtt[m]=(wnc*wnc)*olddzz[m]-(wnb*wnb)*olddzzzz[m];
    } 
    oldt[0]=0;//wnc*wnc/AM*(olddz[0]-olddz[NZ])-KSI/AM*oldd[0];
    oldtt[0]=0;//wnc*wnc/(2*dz*AM)*(4*oldt[1]+4*oldt[NZ-1]-oldt[2]-oldt[NZ-2]-6*oldt[0])-KSI/AM*oldt[0];  
 //  ROOTONLY fprintf(stderr,"the calculation is finished\n");
    
for (m = 0; m < NZ; m += 1) {
     map->d[m]=oldd[m];
     map->t[m]=oldt[m];
     map->tt[m]=oldtt[m];}
	 
      realft (NZ/2, map->d, -1);               
      realft (NZ/2, map->t, -1);
      realft (NZ/2, map->tt,-1);
      realft (NZ/2, map->f,-1);
	  
      map->d[1] = map->t[1] = map->tt[1] = 0.;


    break;
  case -6: //rigid vibration
  {
     double oldd, oldt, oldtt;
     double subd, subt, subtt, k1d, k2d, k3d, k4d, k1t, k2t, k3t, k4t;
 	
    lambda = wn*wn;

	   {
      subd=map->d[0];
      subt=map->t[0];
      subtt=map->tt[0];
      oldd=map->d[0];
      oldt=map->t[0];
      oldtt=map->tt[0];
	   }

/**RK4 first stage***/
   {
     k1d=subt;
  	 k1t=map->f[0]/massa-lambda*sqrt(mass/massa)*oldd
              -2.*zeta*wn*mass/massa*oldt+addmass/massa*oldtt;

     oldd=subd+k1d*dt/2;
     oldt=subt+k1t*dt/2;
	  }

//
//
//
/**RK4 second stage***/
   {
        k2d=subt+k1t*dt/2;
  	    k2t=map->f[0]/massa-lambda*sqrt(mass/massa)*oldd
                -2.*zeta*wn*mass/massa*oldt+addmass/massa*oldtt;

        oldd=subd+k2d*dt;
        oldt=subt+k2t*dt;
   }
  
//
//
//
/**RK4 third stage***/

   {
      k3d=subt+k2t*dt/2;
  	  k3t=map->f[0]/massa-lambda*sqrt(mass/massa)*oldd
                -2.*zeta*wn*mass/massa*oldt+addmass/massa*oldtt;

      oldd=subd+k3d*dt;
      oldt=subt+k3t*dt;
   }

//
//
//
/**RK4 fourth stage***/
    {
      k4d=subt+k3t*dt;
  	  k4t=map->f[0]/massa-lambda*sqrt(mass/massa)*oldd
                -2.*zeta*wn*mass/massa*oldt+addmass/massa*oldtt;

      oldd=subd+(k1d+2*k2d+2*k3d+k4d)*dt/6;
      oldt=subt+(k1t+2*k2t+2*k3t+k4t)*dt/6;
    }


/**RK4 fifth stage***/


  	  oldtt=map->f[0]/massa-lambda*sqrt(mass/massa)*oldd
                -2.*zeta*wn*mass/massa*oldt+addmass/massa*oldtt;


     map->d[0]=oldd;
     map->t[0]=oldt;
     map->tt[0]=oldtt;

//rigid cylinder only mean plane
     dzero(NZ-1,map->d+1,1);
     dzero(NZ-1,map->t+1,1);
     dzero(NZ-1,map->tt+1,1);

  break;
   }
  case -7: //rigid vibration backward Euler
  {
     double oldd, oldt, oldtt;
     double newd, newt, newtt;
     double kd, fd;
 	
     lambda = wn*wn;

     oldd=map->d[0];
     oldt=map->t[0];
     oldtt=map->tt[0];

     kd = lambda+(1.0/dt+2.0*zeta*wn)/dt;
     fd = map->f[0]/mass+oldt/dt+oldd/dt*(1.0/dt+2.0*zeta*wn);

     newd = fd/kd;
     newt = (newd-oldd)/dt;
  	 newtt=map->f[0]/mass-lambda*oldd-2.*zeta*wn*oldt;

     map->d[0]=newd;
     map->t[0]=newt;
     map->tt[0]=newtt;

//rigid cylinder only mean plane
     dzero(NZ-1,map->d+1,1);
     dzero(NZ-1,map->t+1,1);
     dzero(NZ-1,map->tt+1,1);

  break;
   }
  default:
    fprintf(stderr,"Invalid case: Coding error!\n");
    exit (-1);
    break;
  }

free(oldd);
free(oldt);
free(oldtt);
free(olddz);
free(olddzz);
free(olddzzz);
free(olddzzzz);
free(subd);
free(subt);
free(subtt);
free(k1d);
free(k2d);
free(k3d);
free(k4d);
free(k1t);
free(k2t);
free(k3t);
free(k4t);


  return;
}

void update_free (Mapping *map, double dt, double mass, double wn, double wnc, double wnb, double zeta, int forc, double *Tension, double *Stiffness)
{
  int      NZ   = map->NZ;
#ifdef DISC_FORCING
  int 	   NZM = iparam("MNL")  ;
#else
  int 	   NZM = iparam("NZTOT")  ;
#endif
  double addEI = dparam("addEI");
  double freqqq = dparam("freqqq");
  double Amplitude = dparam("Amplitude"); 
  int nztot= option("NZTOT");
  double   beta = 2.*M_PI/dparam("LZ"), fac,
           lambda, zeta_m, frq, keep, save, 
	   keepmapd, keepmapt,
	   Amnk, Bmnk,
//	   k1d, k2d, k3d, k4d, 
//	   k1t, k2t, k3t, k4t,
	   newd, newt, newa;
  double   dz = dparam("LZ")/option("NZTOT");

  double addmass = dparam("DADDMASS");
  double zetan   = dparam("DZETAN");
  double frhs = 0.;
  double massa = mass+addmass;

  int ipiv[2*(NZ-1)];
  int info=0; 

  register int m, n, k,j,i;
  const double   timec = dparam("TIMEC");
  double *z = Zmesh (NZ);
  double t, CoS, SiN, denom,  dfac, ddfac, startime;
  register double CoSZ, loop_phase;
//#if !defined(MAP) && defined(SPM)
  double fmean = 0.;
//#endif

  #ifdef DISC_FORCING
   int diff = NZ-NZM;
   if( (NZM<3) || (diff<3) )
    {
      fprintf(stderr, "NZM = %d, it  has to be in a range of (3, NZ-3) ... \n",NZM);
      exit(-1);
    }
  #endif

  if (timec != 0.0) {
    if (iparam("IGRSTART")) /* restarting needs .rea file info */
      startime = dparam("STARTING");
    else
      startime = dparam("STARTIME");
    t = map->time - startime;
    denom = (t*t + timec);
    fac = t*t/denom;
    dfac = 2.*t*timec/(denom*denom);
    ddfac = -2.*timec*(3.*t*t - timec)/(denom*denom*denom);
  } else
    t = map->time;
  switch (forc) {
  case 0: /* periodic in all derivatives, only constrained by virtual sping */
    if (option("nomean")) { /* Allow no motion in the mean */
      /* We explicitly force all of the variables below to be zero */
      map->d[0] = map->t[0] = map->tt[0] = 0.;
      map->d[1] = map->t[1] = map->tt[1] = 0.;
    } else {
      /* We explicitly force all of the variables below to be zero */
      lambda = wn;
      zeta_m = zeta;
      newmark (map->d, map->t, map->tt, map->f[0], dt, mass, lambda, zeta_m);
      map->d[1] = map->t[1] = map->tt[1] = 0.;
    }
    for (m = 2; m < NZ; m++) {
      frq    = beta * (m/2);
      lambda = frq*frq * (wnc*wnc + frq*frq  * wnb*wnb);
      lambda = sqrt (wn*wn + lambda);
      zeta_m = (lambda*zeta*wn != 0.0 ? zeta * wn/lambda : 0.0);
      newmark (map->d+m, map->t+m, map->tt+m,
	       map->f[m], dt, mass, lambda, zeta_m);
    }
    break;
  case -1: /* Allowed to freely move rigidly */
    lambda = wn;
    lambda = lambda*sqrt(mass/massa);
    fmean = map->f[0]; 

#if !defined(MAP) && defined(SPM)
    fmean = 0.;
    for(m=0; m<NZ; ++m)
     fmean += map->f[m];
    fmean /= (double)NZ; //mean force
#endif

#if defined(MAP) && defined(SPM)
    fmean = map->f[0]; 
#endif

    frhs = fmean+addmass*map->tt[0]+(2.*zetan*wn)*map->t[0];
//    zeta_m = zeta;
    zeta_m = zeta+zetan;

    newmark (map->d, map->t, map->tt, frhs, dt, mass, lambda, zeta_m);

#if !defined(MAP) && defined(SPM)
    for(m=1;m<NZ; ++m) //copy all map stuffs from zero plane
     {
      map->d[m]  = map->d[0];
      map->t[m]  = map->t[0];
      map->tt[m] = map->tt[0];
     }
#endif

#ifdef MAP
    dzero(NZ-1,map->d+1,1);
    dzero(NZ-1,map->t+1,1);
    dzero(NZ-1,map->tt+1,1);
#endif

    break;
 
  case -3: /* Allowed to move but with endpoints fixed to zero */
   {
     
    const int RZ = 2*(NZ-1); 

    double *second_d = (double *) calloc (NZ, sizeof(double));
    double *fourth_d = (double *) calloc (NZ, sizeof(double));
    double *mid_d = (double *) calloc (RZ, sizeof(double));

//  #ifdef OLDFFTS
    realft (NZ/2, map->d,  1);           
    realft (NZ/2, map->t,  1);
    realft (NZ/2, map->tt, 1);
    realft (NZ/2, map->f, 1);    
//  #else  
//    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->d, 1, 0, 0, 0, 0);
//    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->t, 1, 0, 0, 0, 0);
//    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->tt, 1, 0, 0, 0, 0);
//    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->f, 1, 0, 0, 0, 0);
//  #endif
    
    map->d[0] = map->d[NZ]=0.;    /* make sure this is true  */
    map->t[0] = map->t[NZ]=0.;    /* make sure this is true  */
    map->tt[0] = map->tt[NZ]=0.;    /* make sure this is true  */

    for (int i = 0; i < RZ*RZ; i += 1)  Mat_inv[i] = Mat_plus[i];//restore the matrix

     for (int i = 0; i < (NZ-1); i += 1) //for i th row for displacement
//        RHS_f[i] = 2.*map->d[i+1]-0.5*oldd_1[i+1];
        RHS_f[i] = map->d[i+1];

     for (int i = NZ-1, j=0; i < RZ; i += 1, j += 1) //for i th row for displacement
//        RHS_f[i] = dt*map->f[j+1]/mass+2.*map->t[j+1]-0.5*oldt_1[j+1];
        RHS_f[i] = dt*map->f[j+1]/mass+map->t[j+1];

//note following lapack function will change the matrix 
	     dgesv(RZ, 1, &(Mat_inv[0]), RZ, &(ipiv[0]), &(RHS_f[0]), RZ, info);	
       
//       for (int i = 0; i < (NZ-1); i += 1) //for i th row for displacement
//        {
//         oldd_1[i+1]=map->d[i+1];
//         oldt_1[i+1]=map->t[i+1];
//        }

       for (int i = 0; i < (NZ-1); i += 1) //for i th row for displacement
         map->d[i+1] = RHS_f[i];
       
       for (int i = NZ-1, j=0; i < RZ; i += 1, j += 1) //
         map->t[j+1] = RHS_f[i];

    map->d[0] = map->d[NZ]=0.;    /* make sure this is true  */
    map->t[0] = map->t[NZ]=0.;    /* make sure this is true  */
    map->tt[0] = map->tt[NZ]=0.;    /* make sure this is true  */
#if 1
 //d^2 y/ dz^2
   for (int m = 1; m < (NZ); m += 1)   second_d[m]=	(map->d[m+1]-2.*map->d[m]+map->d[m-1])/(dz*dz);
		second_d[0] = 0.; //pinned
    second_d[NZ] = 0.0; //both pinned and spring-supported ends have zero second derivative. 

   for (int m = 2; m < (NZ-1); m += 1) 
     fourth_d[m]= (map->d[m+2]-4.*map->d[m+1]+6*map->d[m]-4.*map->d[m-1]+map->d[m-2])/(dz*dz*dz*dz);

//    fourth_d[1] = (map->d[3]-4.*map->d[2]+5.*map->d[1]-4.*map->d[0])/(dz*dz*dz*dz);
    fourth_d[1] = (map->d[3]-4.*map->d[2]+6.*map->d[1]-4.*map->d[0]+map->d[NZ-1])/(dz*dz*dz*dz);
    fourth_d[0] = fourth_d[1]; //assume 4th order gradient doesn't change we don't have a formula for gradient at this point

//    fourth_d[NZ-1] = (-4.*map->d[NZ]+5.*map->d[NZ-1]-4.*map->d[NZ-2]+map->d[NZ-3])/(dz*dz*dz*dz); //pinned
    fourth_d[NZ-1] = (map->d[1]-4.*map->d[NZ]+6.*map->d[NZ-1]-4.*map->d[NZ-2]+map->d[NZ-3])/(dz*dz*dz*dz); //pinned
    fourth_d[NZ]   = -fourth_d[NZ-1];

//upate d^2/dt^2
    for (int m = 1; m < (NZ); m++) map->tt[m]=map->f[m]/mass+Tension[m]*second_d[m]-Stiffness[m]*fourth_d[m]-2.*zeta*wn*map->t[m];
#else
 	  dgemv('N',RZ,RZ,1,Mat_plus,RZ,RHS_f,1,0,mid_d,1)	;	
    
    for (int i = 0; i < NZ-1; i += 1) //
		 map->tt[i+1] =  -mid_d[i]+map->f[i+1]/mass - 2.*zeta*wn*map->t[i+1];	

#endif  
     free(second_d);
     free(fourth_d);
//     free(Mat_plus);
//     free(RHS_f);
     free(mid_d);
//#ifdef OLDFFTS
     realft (NZ/2, map->d, -1);               
     realft (NZ/2, map->t, -1);
     realft (NZ/2, map->tt,-1);
     realft (NZ/2, map->f,-1);
//#else
//     rfftw(rplan, 1, (FFTW_COMPLEX *) map->d, 1, 0, 0, 0, 0);
//     rfftw(rplan, 1, (FFTW_COMPLEX *) map->t, 1, 0, 0, 0, 0);
//     rfftw(rplan, 1, (FFTW_COMPLEX *) map->tt, 1, 0, 0, 0, 0);
//     rfftw(rplan, 1, (FFTW_COMPLEX *) map->f, 1, 0, 0, 0, 0);
//#endif
//     map->d[1] = map->t[1] = map->tt[1] = 0.; //make it be periodic ? //zwang 02272018

//    n = iparam ("NFILTER");
//
//  for (m = 0; m < NZ; m ++)
//     if( m>(NZ-n) ){  
//       map->f[m] = 0.;
//       map->t[m] = 0.;
//       map->tt[m] = 0.;
//     }
     
//     map_step ++;
    break;
   }
  case -4:
  { 
     double *oldd, *oldt, *oldtt, *olddz, *olddzz, *olddzzz, *olddzzzz;
     double *subd, *subt, *subtt, *k1d, *k2d, *k3d, *k4d, *k1t, *k2t, *k3t, *k4t;

     oldd = (double *) calloc (NZ+1, sizeof(double));
     oldt = (double *) calloc (NZ+1, sizeof(double));
     oldtt = (double *) calloc (NZ+1, sizeof(double));
     olddz = (double *) calloc (NZ+1, sizeof(double));
     olddzz = (double *) calloc (NZ+1, sizeof(double));
     olddzzz = (double *) calloc (NZ+1, sizeof(double));
     olddzzzz = (double *) calloc (NZ+1, sizeof(double));
     subd = (double *) calloc (NZ, sizeof(double));
     subt = (double *) calloc (NZ, sizeof(double));
     subtt = (double *) calloc (NZ, sizeof(double));
     k1d = (double *) calloc (NZ, sizeof(double));
     k2d = (double *) calloc (NZ, sizeof(double));
     k3d = (double *) calloc (NZ, sizeof(double));
     k4d = (double *) calloc (NZ, sizeof(double));
     k1t = (double *) calloc (NZ, sizeof(double));
     k2t = (double *) calloc (NZ, sizeof(double));
     k3t = (double *) calloc (NZ, sizeof(double));
     k4t = (double *) calloc (NZ, sizeof(double)); 

  #ifdef OLDFFTS
    realft (NZ/2, map->d,  1);           
    realft (NZ/2, map->t,  1);
    realft (NZ/2, map->tt, 1);
    realft (NZ/2, map->f, 1);    
  #else  
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->d, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->t, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->tt, 1, 0, 0, 0, 0);
    rfftw(rplan_inv, 1, (FFTW_COMPLEX *) map->f, 1, 0, 0, 0, 0);
  #endif
 	
   	map->d[0] = map->t[0] = map->tt[0] = 0.0; 
     for (m = 0; m < (NZ); m += 1) 
	 {
     subd[m]=map->d[m];
     subt[m]=map->t[m];
     subtt[m]=map->tt[m];
     oldd[m]=map->d[m];
     oldt[m]=map->t[m];
     oldtt[m]=map->tt[m];
	 }
	 
   oldd[NZ]=oldd[0];
   oldt[NZ]=oldt[0];
   oldtt[NZ]=oldtt[0];

#ifdef DISC_FORCING
   oldd[NZM]= 0.;
   oldt[NZM]= 0.;
   oldtt[NZM]= 0.;
#endif
/**RK4 first stage***/

   Derivatives_PinnedPinned(oldd, olddzz, olddzzzz, NZ, NZM, dz);


#ifdef DISC_FORCING
   for (int m = 1; m < (NZM); m += 1)
#else
   for (int m = 1; m < (NZ); m += 1) 
#endif
   {
     k1d[m]=subt[m];
  	 k1t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
              -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

     oldd[m]=subd[m]+k1d[m]*dt/2;
     oldt[m]=subt[m]+k1t[m]*dt/2;
	  }

#ifdef DISC_FORCING
   for (int m = NZM+1; m < (NZ); m += 1) {
      k1d[m]=subt[m];
//  	  k1t[m]=Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
//              -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];
  	 k1t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
              -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

      oldd[m]=subd[m]+k1d[m]*dt/2;
      oldt[m]=subt[m]+k1t[m]*dt/2;
	  }
#endif
//
//
//
/**RK4 second stage***/
   oldd[NZ]=oldd[0]=0.;
   oldt[NZ]=oldt[0]=0.;

#ifdef DISC_FORCING
   oldd[NZM]=oldt[NZM]=0.;
#endif
	  
   Derivatives_PinnedPinned(oldd, olddzz, olddzzzz, NZ, NZM, dz);

#ifdef DISC_FORCING
   for (int m = 1; m < (NZM); m += 1) 
#else
	 for (int m = 1; m < (NZ); m++)
#endif
   {
        k2d[m]=subt[m]+k1t[m]*dt/2;
  	    k2t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

        oldd[m]=subd[m]+k2d[m]*dt;
        oldt[m]=subt[m]+k2t[m]*dt;
   }
  
#ifdef DISC_FORCING
   for (int m = NZM+1; m < (NZ); m += 1) {
        k2d[m]=subt[m]+k1t[m]*dt/2;
//  	    k2t[m]=Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
//                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];
  	    k2t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

        oldd[m]=subd[m]+k2d[m]*dt;
        oldt[m]=subt[m]+k2t[m]*dt;
     }
#endif
//
//
//
/**RK4 third stage***/
   oldd[NZ]=oldd[0]=0.;
   oldt[NZ]=oldt[0]=0.;
#ifdef DISC_FORCING
   oldd[NZM]=oldt[NZM]=0.;
#endif

   Derivatives_PinnedPinned(oldd, olddzz, olddzzzz, NZ, NZM, dz);
#ifdef DISC_FORCING
   for (int m = 1; m < (NZM); m += 1) 
#else
	 for (int m = 1; m < (NZ); m++) 
#endif
   {
      k3d[m]=subt[m]+k2t[m]*dt/2;
  	  k3t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

      oldd[m]=subd[m]+k3d[m]*dt;
      oldt[m]=subt[m]+k3t[m]*dt;
   }

#ifdef DISC_FORCING
   for (int m = NZM+1; m < (NZ); m += 1) {
      k3d[m]=subt[m]+k2t[m]*dt/2;
//  	  k3t[m]=Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
//                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];
  	  k3t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

      oldd[m]=subd[m]+k3d[m]*dt;
      oldt[m]=subt[m]+k3t[m]*dt;
    }
#endif
//
//
//
/**RK4 fourth stage***/
   oldd[NZ]=oldd[0]=0.;
   oldt[NZ]=oldt[0]=0.;

#ifdef DISC_FORCING
   oldd[NZM]=oldt[NZM]=0.;
#endif
   Derivatives_PinnedPinned(oldd, olddzz, olddzzzz, NZ, NZM, dz);

#ifdef DISC_FORCING
   for (int m = 1; m < (NZM); m += 1) 
#else
	  for (int m = 1; m < (NZ); m++) 
#endif
    {
      k4d[m]=subt[m]+k3t[m]*dt;
  	  k4t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

      oldd[m]=subd[m]+(k1d[m]+2*k2d[m]+2*k3d[m]+k4d[m])*dt/6;
      oldt[m]=subt[m]+(k1t[m]+2*k2t[m]+2*k3t[m]+k4t[m])*dt/6;
    }

#ifdef DISC_FORCING
   for (int m = NZM+1; m < (NZ); m += 1) {
      k4d[m]=subt[m]+k3t[m]*dt;
//  	  k4t[m]=Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
//                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];
  	  k4t[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

      oldd[m]=subd[m]+(k1d[m]+2*k2d[m]+2*k3d[m]+k4d[m])*dt/6;
      oldt[m]=subt[m]+(k1t[m]+2*k2t[m]+2*k3t[m]+k4t[m])*dt/6;
      }
#endif


/**RK4 fifth stage***/
   oldd[NZ]=oldd[0]=0.;
   oldt[NZ]=oldt[0]=0.;

#ifdef DISC_FORCING
   oldd[NZM]=oldt[NZM]=0.;
#endif
   Derivatives_PinnedPinned(oldd, olddzz, olddzzzz, NZ, NZM, dz);

#ifdef DISC_FORCING
   for (int m = 1; m < (NZM); m += 1)
#else
   for (int m = 1; m < (NZ); m++) 
#endif
  	  oldtt[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];

#ifdef DISC_FORCING
   for (int m = NZM+1; m < (NZ); m += 1) 
//  	  oldtt[m]=Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
//                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];
  	  oldtt[m]=map->f[m]/massa+Tension[m]*mass/massa*olddzz[m]-Stiffness[m]*mass/massa*olddzzzz[m]
                -2.*zeta*wn*mass/massa*oldt[m]+addmass/massa*oldtt[m];
#endif
   
    oldt[0]=0;
    oldtt[0]=0;
#ifdef DISC_FORCING
   oldt[NZM]=oldtt[NZM]=0.;
#endif


for (m = 0; m < NZ; m += 1) {
     map->d[m]=oldd[m];
      map->t[m]=oldt[m];
     map->tt[m]=oldtt[m];}
	 
  #ifdef OLDFFTS
     realft (NZ/2, map->d, -1);               
     realft (NZ/2, map->t, -1);
     realft (NZ/2, map->tt,-1);
     realft (NZ/2, map->f,-1);
  #else
     rfftw(rplan, 1, (FFTW_COMPLEX *) map->d, 1, 0, 0, 0, 0);
     rfftw(rplan, 1, (FFTW_COMPLEX *) map->t, 1, 0, 0, 0, 0);
     rfftw(rplan, 1, (FFTW_COMPLEX *) map->tt, 1, 0, 0, 0, 0);
     rfftw(rplan, 1, (FFTW_COMPLEX *) map->f, 1, 0, 0, 0, 0);
  #endif
/*	  
    n = iparam ("NFILTER");

  for (m = 0; m < NZ; m ++)
     if( m>(NZ-n) ){  
       map->d[m] = 0.;
       map->t[m] = 0.;
       map->tt[m] = 0.;
     }
*/
//     map->d[1] = map->t[1] = map->tt[1] = 0.;

     free(oldd);
     free(oldt);
     free(oldtt);
     free(olddz);
     free(olddzz);
     free(olddzzz);
     free(olddzzzz);
     free(subd);
     free(subt);
     free(subtt);
     free(k1d);
     free(k2d);
     free(k3d);
     free(k4d);
     free(k1t);
     free(k2t);
     free(k3t);
     free(k4t);

  break;
   }
  default:
    fprintf(stderr,"Invalid case: Coding error!\n");
    exit (-1);
    break;
  }
  return;
}


void Derivatives_PinnedPinned(double *oldd, double *olddzz, double *olddzzzz, int NZ, int NZM, double dz)
 {


 //d^2 y/ dz^2
  #ifdef DISC_FORCING
   for (int m = 1; m < (NZM); m += 1)   olddzz[m]=	(oldd[m+1]-2.*oldd[m]+oldd[m-1])/(dz*dz);

   for (int m = NZM+1; m < (NZ); m += 1)   olddzz[m]=	(oldd[m+1]-2.*oldd[m]+oldd[m-1])/(dz*dz);
   
   olddzz[NZM] = 0.0; //pinned end have zero second derivative. 
  #else
   for (int m = 1; m < (NZ); m += 1)    olddzz[m]=	(oldd[m+1]-2.*oldd[m]+oldd[m-1])/(dz*dz);
  #endif

		olddzz[0] = 0.; //pinned
    olddzz[NZ] = 0.0; //both pinned and spring-supported ends have zero second derivative. 


  #ifdef DISC_FORCING
   for (int m = 2; m < (NZM-1); m += 1) 
     olddzzzz[m]= (oldd[m+2]-4.*oldd[m+1]+6*oldd[m]-4.*oldd[m-1]+oldd[m-2])/(dz*dz*dz*dz);

   for (int m = NZM+2; m < (NZ-1); m += 1) 
     olddzzzz[m]= (oldd[m+2]-4.*oldd[m+1]+6*oldd[m]-4.*oldd[m-1]+oldd[m-2])/(dz*dz*dz*dz);

    olddzzzz[NZM-1] = (-4.*oldd[NZM]+5.*oldd[NZM-1]-4.*oldd[NZM-2]+oldd[NZM-3])/(dz*dz*dz*dz); 
    olddzzzz[NZM] = 0.;//we don't solve this point 
    olddzzzz[NZM+1] = (-4.*oldd[NZM]+5.*oldd[NZM+1]-4.*oldd[NZM+2]+oldd[NZM+3])/(dz*dz*dz*dz); 
  #else
   for (int m = 2; m < (NZ-1); m += 1) 
     olddzzzz[m]= (oldd[m+2]-4.*oldd[m+1]+6*oldd[m]-4.*oldd[m-1]+oldd[m-2])/(dz*dz*dz*dz);
  #endif

    olddzzzz[0] = 0.;  //don't know the fourth derivative for this point but it's OK since we don't solve this point
    olddzzzz[1] = (oldd[3]-4.*oldd[2]+5.*oldd[1]-4.*oldd[0])/(dz*dz*dz*dz);

    olddzzzz[NZ-1] = (-4.*oldd[NZ]+5.*oldd[NZ-1]-4.*oldd[NZ-2]+oldd[NZ-3])/(dz*dz*dz*dz); //pinned
    olddzzzz[NZ]   = 0.;//pinned, we don't solve this point

 }

void generate_structure_matrix( double wn, double zeta, double *Tension, double *Stiffness )
 {
      const int      NZ = option("NZTOT");
      const int      RZ = 2*(NZ-1); 

      const double   dz = dparam("LZ")/option("NZTOT");
      const double   dt = dparam("DELT");

      Mat_plus = (double *) calloc  (RZ*RZ, sizeof(double));
      Mat_inv = (double *) calloc  (RZ*RZ, sizeof(double));
      RHS_f = (double *) calloc  (RZ, sizeof(double));

    for (int i = 0; i < RZ*RZ; i += 1) {Mat_plus[i] = 0.; Mat_inv[i] = 0.;}//clear;

//first step gennerate the matrix... 
    for (int i = 0; i < (NZ-1); i += 1) //for i th row
      {
         int j= i;
//         Mat_plus[i+j*RZ] = 1.5;  //top-left identity matrix
         Mat_plus[i+j*RZ] = 1.;  //top-left identity matrix
      }

    for (int i = 0; i < (NZ-1); i += 1) //for i th row
      {
         int j = NZ-1+i; 
         Mat_plus[i+j*RZ] = -dt;  //top-right identity matrix
      }
      
    for (int i = NZ-1; i < RZ; i += 1) //for i th row
     {
        int j=i;
//         Mat_plus[i+j*RZ] = 1.5+2.*zeta*wn*dt;  //bottom-right identity matrix
         Mat_plus[i+j*RZ] = 1.+2.*zeta*wn*dt;  //bottom-right identity matrix
     }
    
#if 0   
    //low-left part
    for (int i = NZ+1; i < RZ-2; i += 1) //for i th row
      {
        int j = i-(NZ-1); //major diagonal
        //main diagonal term
         Mat_plus[i+j*RZ] = dt/dz/dz*(15./12.*(Tension[j+2]+Tension[j])+3.*(Stiffness[j+2]+Stiffness[j])/dz/dz);  // j^th row of the matrix corresponds to j+1 plane ! 
        //left two terms
         Mat_plus[i+(j-1)*RZ] = dt/dz/dz*(-8./12.*(Tension[j+2]+Tension[j])-2.*(Stiffness[j+2]+Stiffness[j])/dz/dz);  //
         Mat_plus[i+(j-2)*RZ] = dt/dz/dz*(1./24.*(Tension[j+2]+Tension[j])+0.5*(Stiffness[j+2]+Stiffness[j])/dz/dz);  //
        //right two terms
         Mat_plus[i+(j+1)*RZ] = dt/dz/dz*(-8./12.*(Tension[j+2]+Tension[j])-2.*(Stiffness[j+2]+Stiffness[j])/dz/dz);  //
         Mat_plus[i+(j+2)*RZ] = dt/dz/dz*(1./24.*(Tension[j+2]+Tension[j])+0.5*(Stiffness[j+2]+Stiffness[j])/dz/dz);  //

      }
 
    //first row j=0
//        Mat_plus[(NZ-1)+0*RZ] = dt/dz/dz*(29./24.*(Tension[2]+Tension[0])+2.5*(Stiffness[2]+Stiffness[0])/dz/dz);  //
        Mat_plus[(NZ-1)+0*RZ] = dt/dz/dz*(15./12.*(Tension[2]+Tension[0])+3.0*(Stiffness[2]+Stiffness[0])/dz/dz);  // periodic  1 == NZ-1
        Mat_plus[(NZ-1)+1*RZ] = dt/dz/dz*(-8./12.*(Tension[2]+Tension[0])-2.*(Stiffness[2]+Stiffness[0])/dz/dz);  //
        Mat_plus[(NZ-1)+2*RZ] = dt/dz/dz*(1./24.*(Tension[2]+Tension[0])+0.5*(Stiffness[2]+Stiffness[0])/dz/dz);  //
        Mat_plus[(NZ-1)+(NZ-2)*RZ] = dt/dz/dz*(1./24.*(Tension[2]+Tension[0])+0.5*(Stiffness[2]+Stiffness[0])/dz/dz);  // periodic NZ-1 == 1
    //second row j=1
        Mat_plus[(NZ)+0*RZ] = dt/dz/dz*(-8./12.*(Tension[1]+Tension[3])-2.*(Stiffness[1]+Stiffness[3])/dz/dz);  //
        Mat_plus[(NZ)+1*RZ] = dt/dz/dz*(15./12*(Tension[1]+Tension[3])+3.*(Stiffness[1]+Stiffness[3])/dz/dz);  //
        Mat_plus[(NZ)+2*RZ] = dt/dz/dz*(-8./12.*(Tension[1]+Tension[3])-2.*(Stiffness[1]+Stiffness[3])/dz/dz);  //
        Mat_plus[(NZ)+3*RZ] = dt/dz/dz*(1./24.*(Tension[1]+Tension[3])+0.5*(Stiffness[1]+Stiffness[3])/dz/dz);  //

   //row NZ-2  second last
        Mat_plus[(RZ-2)+(NZ-5)*RZ] = dt/dz/dz*(1./24.*(Tension[NZ-1]+Tension[NZ-3])+0.5*(Stiffness[NZ-1]+Stiffness[NZ-3])/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-4)*RZ] = dt/dz/dz*(-8./12.*(Tension[NZ-1]+Tension[NZ-3])-2.*(Stiffness[NZ-3]+Stiffness[NZ-1])/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-3)*RZ] = dt/dz/dz*(15./12.*(Tension[NZ-3]+Tension[NZ-1])+3.*(Stiffness[NZ-3]+Stiffness[NZ-1])/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-2)*RZ] = dt/dz/dz*(-8./12.*(Tension[NZ-3]+Tension[NZ-1])-2.*(Stiffness[NZ-3]+Stiffness[NZ-1])/dz/dz);  //
  //last row
        Mat_plus[(RZ-1)+0*RZ] = dt/dz/dz*(1./24.*(Tension[NZ]+Tension[NZ-2])+0.5*(Stiffness[NZ]+Stiffness[NZ-2])/dz/dz);  // periodic NZ-1 == 1
        Mat_plus[(RZ-1)+(NZ-4)*RZ] = dt/dz/dz*(1./24.*(Tension[NZ]+Tension[NZ-2])+0.5*(Stiffness[NZ]+Stiffness[NZ-2])/dz/dz);  //
        Mat_plus[(RZ-1)+(NZ-3)*RZ] = dt/dz/dz*(-8./12.*(Tension[NZ]+Tension[NZ-2])-2.*(Stiffness[NZ]+Stiffness[NZ-2])/dz/dz);  //
//        Mat_plus[(RZ-1)+(NZ-2)*RZ] = dt/dz/dz*(29./24.*(Tension[NZ]+Tension[NZ-2])+2.5*(Stiffness[NZ]+Stiffness[NZ-2])/dz/dz);   //
        Mat_plus[(RZ-1)+(NZ-2)*RZ] = dt/dz/dz*(15./12.*(Tension[NZ]+Tension[NZ-2])+3.0*(Stiffness[NZ]+Stiffness[NZ-2])/dz/dz);   //
#endif
    //low-left part
    for (int i = NZ+1; i < RZ-2; i += 1) //for i th row
     {
        int j = i-(NZ-1); //major diagonal
        //main diagonal term
         Mat_plus[i+j*RZ] = dt/dz/dz*(2.*Tension[j+1]+6.*Stiffness[j+1]/dz/dz);  // j^th row of the matrix corresponds to j+1 plane ! 
        //left two terms
         Mat_plus[i+(j-1)*RZ] = dt/dz/dz*(-Tension[j+1]-4.*Stiffness[j+1]/dz/dz);  //
         Mat_plus[i+(j-2)*RZ] = dt/dz/dz*(Stiffness[j+1]/dz/dz);  //
        //right two terms
         Mat_plus[i+(j+1)*RZ] = dt/dz/dz*(-Tension[j+1]-4.*Stiffness[j+1]/dz/dz);  //
         Mat_plus[i+(j+2)*RZ] = dt/dz/dz*(Stiffness[j+1]/dz/dz);  //
      }
    //first row j=0
        Mat_plus[(NZ-1)+0*RZ] = dt/dz/dz*(2.*Tension[1]+6.*Stiffness[1]/dz/dz);  //
        Mat_plus[(NZ-1)+1*RZ] = dt/dz/dz*(-Tension[1]-4.*Stiffness[1]/dz/dz);  //
        Mat_plus[(NZ-1)+2*RZ] = dt/dz/dz*(Stiffness[1]/dz/dz);  //
        Mat_plus[(NZ-1)+(NZ-2)*RZ] = dt/dz/dz*(1.*(Stiffness[1])/dz/dz);  // periodic NZ-1 == 1
    //second row j=1
        Mat_plus[(NZ)+0*RZ] = dt/dz/dz*(-Tension[2]-4.*Stiffness[2]/dz/dz);  //
        Mat_plus[(NZ)+1*RZ] = dt/dz/dz*(2.*Tension[2]+6.*Stiffness[2]/dz/dz);  //
        Mat_plus[(NZ)+2*RZ] = dt/dz/dz*(-Tension[2]-4.*Stiffness[2]/dz/dz);  //
        Mat_plus[(NZ)+3*RZ] = dt/dz/dz*(Stiffness[2]/dz/dz);  //

   //row NZ-2  second last
        Mat_plus[(RZ-2)+(NZ-5)*RZ] = dt/dz/dz*(Stiffness[NZ-2]/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-4)*RZ] = dt/dz/dz*(-Tension[NZ-2]-4.*Stiffness[NZ-2]/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-3)*RZ] = dt/dz/dz*(2*Tension[NZ-2]+6.*Stiffness[NZ-2]/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-2)*RZ] = dt/dz/dz*(-Tension[NZ-2]-4.*Stiffness[NZ-2]/dz/dz);  //
  //last row
        Mat_plus[(RZ-1)+0*RZ] = dt/dz/dz*(1.*(Stiffness[NZ-1])/dz/dz);  // periodic NZ-1 == 1
        Mat_plus[(RZ-1)+(NZ-4)*RZ] = dt/dz/dz*(Stiffness[NZ-1]/dz/dz);  //
        Mat_plus[(RZ-1)+(NZ-3)*RZ] = dt/dz/dz*(-Tension[NZ-1]-4.*Stiffness[NZ]/dz/dz);  //
        Mat_plus[(RZ-1)+(NZ-2)*RZ] = dt/dz/dz*(2.*Tension[NZ-1]+6.*Stiffness[NZ]/dz/dz);   //
 #if 0
     //low-left part
    for (int i = NZ+1; i < RZ-2; i += 1) //for i th row
      {
        int j = i-(NZ-1); //major diagonal
        //main diagonal term
         Mat_plus[i+j*RZ] = dt/dz/dz*(15./12.*(Tension[j+1]+Tension[j+1])+3.*(Stiffness[j+1]+Stiffness[j+1])/dz/dz);  // j^th row of the matrix corresponds to j+1 plane ! 
        //left two terms
         Mat_plus[i+(j-1)*RZ] = dt/dz/dz*(-8./12.*(Tension[j+1]+Tension[j+1])-2.*(Stiffness[j+1]+Stiffness[j+1])/dz/dz);  //
         Mat_plus[i+(j-2)*RZ] = dt/dz/dz*(1./24.*(Tension[j+1]+Tension[j+1])+0.5*(Stiffness[j+1]+Stiffness[j+1])/dz/dz);  //
        //right two terms
         Mat_plus[i+(j+1)*RZ] = dt/dz/dz*(-8./12.*(Tension[j+1]+Tension[j+1])-2.*(Stiffness[j+1]+Stiffness[j+1])/dz/dz);  //
         Mat_plus[i+(j+2)*RZ] = dt/dz/dz*(1./24.*(Tension[j+1]+Tension[j+1])+0.5*(Stiffness[j+1]+Stiffness[j+1])/dz/dz);  //

      }
 
    //first row j=0
        Mat_plus[(NZ-1)+0*RZ] = dt/dz/dz*(29./24.*(Tension[1]+Tension[1])+2.5*(Stiffness[1]+Stiffness[1])/dz/dz);  // pinned
//        Mat_plus[(NZ-1)+0*RZ] = dt/dz/dz*(31./24.*(Tension[1]+Tension[1])+3.*(Stiffness[1]+Stiffness[1])/dz/dz);  // fixed
        Mat_plus[(NZ-1)+1*RZ] = dt/dz/dz*(-8./12.*(Tension[1]+Tension[1])-2.*(Stiffness[1]+Stiffness[1])/dz/dz);  //
        Mat_plus[(NZ-1)+2*RZ] = dt/dz/dz*(1./24.*(Tension[2]+Tension[0])+0.5*(Stiffness[1]+Stiffness[1])/dz/dz);  //
    //second row j=1
        Mat_plus[(NZ)+0*RZ] = dt/dz/dz*(-8./12.*(Tension[2]+Tension[2])-2.*(Stiffness[2]+Stiffness[2])/dz/dz);  //
        Mat_plus[(NZ)+1*RZ] = dt/dz/dz*(15./12*(Tension[2]+Tension[2])+3.*(Stiffness[2]+Stiffness[2])/dz/dz);  //
        Mat_plus[(NZ)+2*RZ] = dt/dz/dz*(-8./12.*(Tension[2]+Tension[2])-2.*(Stiffness[2]+Stiffness[2])/dz/dz);  //
        Mat_plus[(NZ)+3*RZ] = dt/dz/dz*(1./24.*(Tension[2]+Tension[2])+0.5*(Stiffness[2]+Stiffness[2])/dz/dz);  //

   //row NZ-2  second last
        Mat_plus[(RZ-2)+(NZ-5)*RZ] = dt/dz/dz*(1./24.*(Tension[NZ-2]+Tension[NZ-2])+0.5*(Stiffness[NZ-2]+Stiffness[NZ-2])/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-4)*RZ] = dt/dz/dz*(-8./12.*(Tension[NZ-2]+Tension[NZ-2])-2.*(Stiffness[NZ-2]+Stiffness[NZ-2])/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-3)*RZ] = dt/dz/dz*(15./12.*(Tension[NZ-2]+Tension[NZ-2])+3.*(Stiffness[NZ-2]+Stiffness[NZ-2])/dz/dz);  //
        Mat_plus[(RZ-2)+(NZ-2)*RZ] = dt/dz/dz*(-8./12.*(Tension[NZ-2]+Tension[NZ-2])-2.*(Stiffness[NZ-2]+Stiffness[NZ-2])/dz/dz);  //
  //last row
        Mat_plus[(RZ-1)+(NZ-4)*RZ] = dt/dz/dz*(1./24.*(Tension[NZ-1]+Tension[NZ-1])+0.5*(Stiffness[NZ-1]+Stiffness[NZ-1])/dz/dz);  //
        Mat_plus[(RZ-1)+(NZ-3)*RZ] = dt/dz/dz*(-8./12.*(Tension[NZ-1]+Tension[NZ-1])-2.*(Stiffness[NZ-1]+Stiffness[NZ-1])/dz/dz);  //
        Mat_plus[(RZ-1)+(NZ-2)*RZ] = dt/dz/dz*(29./24.*(Tension[NZ-1]+Tension[NZ-1])+2.5*(Stiffness[NZ-1]+Stiffness[NZ-1])/dz/dz);   //pinned
//        Mat_plus[(RZ-1)+(NZ-2)*RZ] = dt/dz/dz*(31./24.*(Tension[NZ]+Tension[NZ-1])+3.*(Stiffness[NZ-1]+Stiffness[NZ-1])/dz/dz);   //fixed
   #endif
 }
/* ------------------------------------------------------------------------ *
 * newmark() - Newmark integration of:
 *             Xtt + 2*zeta*wn * Xt + wn*wn * X = f/m
 * ------------------------------------------------------------------------ */

void newmark (double *x, double *v, double *a,
		     double  f, double dt, double  m, double wn, double zeta)
{
  double x_, v_, a_,
         d, k, c,
         a11, a12, a13, a21, a22, a23, a31, a32, a33, det;

  d   =  dt/2;
  k   =  wn*wn;
  c   =  2*zeta*wn;
  a11 =  1 + d*c;
  a12 = -d*d*c;
  a13 =  d*d;
  a21 = -d*k;
  a22 =  1 + d*d*k;
  a23 =  d;
  a31 = -k;
  a32 = -c;
  a33 =  1;
  det =  1 + d*d*k + d*c;
  
  x_ = (*x) + dt*(*v) + (dt*dt/4)*(*a);
  v_ =           (*v) +    (dt/2)*(*a);
  a_ = f/m;
  
  *x = (a11 * x_ + a12 * v_ + a13 * a_) / det;
  *v = (a21 * x_ + a22 * v_ + a23 * a_) / det;
  *a = (a31 * x_ + a32 * v_ + a33 * a_) / det;

  return;
}


#ifdef MAP
/* ------------------------------------------------------------------------ *
 *
 * IO
 * ==
 *
 * ReadMap()
 * readMap()
 * WriteMap()
 * writeMap()
 * 
 * ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ *
 * ReadMap() - read the map file
 * ------------------------------------------------------------------------ */

int ReadMap (Domain *omega)
{
  char  fname[FILENAME_MAX];
  FILE *fp;
  
  sprintf (fname, "%s.map.rst", omega->name);
  if ((fp = fopen (fname,"r")) == (FILE *) NULL) {
    ErrorHandler ("ReadMap", "failed to open the map file", WARNING);
    return 0;
  } else {
    readMap (fp, omega->mapx);
    readMap (fp, omega->mapy);
    fclose  (fp);
  }
   // check for consistent times 

  if ((omega->mapx->time != omega->mapy->time)) {
      ROOT fprintf(stderr, "Different start times for mapx and mapy!\n");
    // make them agree
    omega->mapy->time = omega->mapx->time;
 	
  }




  return 1;
}


/* ------------------------------------------------------------------------ *
 * readMap() - read a single map
 * ------------------------------------------------------------------------ */

void readMap (FILE *fp, Mapping *map)
{
  char   buf[BUFSIZ];
  int    k, NZ, nZ;
  double time, dummy;

  /* read up to -- and then past -- next set of comment lines */
  
  while ((*buf = getc(fp)) != '#') fgets (buf, BUFSIZ, fp); ungetc (*buf, fp); 
  while ((*buf = getc(fp)) == '#') fgets (buf, BUFSIZ, fp); ungetc (*buf, fp); 
  
  fscanf (fp, "%lf", &time);
  map->time = time;

  fscanf (fp, "%d",  &NZ);
  nZ = min(NZ,map->NZ); // To avoid overwrites
  
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->d   + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy); // skip over rest of line
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->z   + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->zz  + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->t   + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->tt  + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->tz  + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->tzz + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);
  for (k = 0; k < nZ; k++) fscanf (fp, "%lf", map->f   + k);
  for (; k < NZ; k++) fscanf (fp, "%lf", &dummy);

  if (option("nomean")) {
    /* We explicitly force all of the variables below to be zero */
    map->d[0] = map->t[0] = map->tt[0] = 0.;
  }

  return;
}

/* ------------------------------------------------------------------------ *
 * WriteMap() - write the map file (for restarts etc.)
 * ------------------------------------------------------------------------ */

void WriteMap (Domain *omega)
{
  char  fname[FILENAME_MAX];
  FILE *fp;
  static int mapstep = 0, map_number = 0;
  
  if (option("SLICES")) {
    sprintf (fname, "%s_%d.map",  omega->name, map_number);
    ++map_number;
    if ((fp = fopen (fname,"w")) == (FILE *) NULL)
      ErrorHandler ("WriteMap", "failed to open the map file", ERROR);
  } else {
    sprintf (fname, "%s.map", omega->name);
    int err = backup (fname);
    if (err != 0 && mapstep != 0) 
      ErrorHandler ("WriteMap", "failed to backup the map file", WARNING);
    if ((fp = fopen (fname,"w")) == (FILE *) NULL)
      ErrorHandler ("WriteMap", "failed to open the map file", ERROR);
  }
  /* Write the header just to be safe */
  
  mapstep = min(iparam("NSTEPS"),mapstep+iparam("IOSTEP"));
  writeheader (fp,
	       omega->name,
	       LGmax, 
	       dparam("LZ"),
	       option("NZTOT"),
	       omega->U->nel,
	       mapstep,                   /* ignore number of steps */
	       omega->mapx->time,
	       "maps");
  
  /* Write additional information */
  
  fprintf (fp, "%-25.6g Forcx\n", dparam("FORCX"));
  fprintf (fp, "%-25.6g Ampx \n", dparam("AMPX"));
  fprintf (fp, "%-25.6g Freqx\n", dparam("FREQX"));
  fprintf (fp, "%-25.6g Phizx\n", dparam("PHIZX"));
  fprintf (fp, "%-25.6g Phitx\n", dparam("PHITX"));
  fprintf (fp, "%-25.6g Forcy\n", dparam("FORCY"));
  fprintf (fp, "%-25.6g Ampy \n", dparam("AMPY"));
  fprintf (fp, "%-25.6g Freqy\n", dparam("FREQY"));
  fprintf (fp, "%-25.6g Phizy\n", dparam("PHIZY"));
  fprintf (fp, "%-25.6g Phity\n", dparam("PHITY"));
 
  fprintf (fp, "%-25.6g Natural Frequency  \n", dparam("WN"));
  fprintf (fp, "%-25.6g Mass               \n", dparam("ZMASS"));
  fprintf (fp, "%-25.6g Damping Coefficient\n", dparam("ZETA"));
  fprintf (fp, "%-25.6g Phase Speed cable  \n", dparam("WNC"));
  fprintf (fp, "%-25.6g Phase Speed beam   \n", dparam("WNB"));

  
  /* Write out the maps */
  
  fprintf (fp, "#\n# mapx \n#\n"); writeMap (fp, omega->mapx );
  fprintf (fp, "#\n# mapy \n#\n"); writeMap (fp, omega->mapy );
 
  fclose (fp);
  
  return;
}

/* ------------------------------------------------------------------------ *
 * writeMap() - write a single map
 * ------------------------------------------------------------------------ */

void writeMap (FILE *fp, Mapping *map)
{
  int k, NZ = map->NZ;
  
  fprintf (fp, "%-25.6g\n", map->time);
  fprintf (fp, "%-25d  \n", map->NZ);

  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->d[k]);   fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->z[k]);   fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->zz[k]);  fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->t[k]);   fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->tt[k]);  fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->tz[k]);  fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->tzz[k]); fprintf (fp, "\n");
  for (k = 0; k < NZ; k++) fprintf (fp, "%#20.18E ", map->f[k]);   fprintf (fp, "\n");

  return;
}

#ifdef TO_DELETE
/* ------------------------------------------------------------------------ *
 * readmap() - for sem2tec
 * ------------------------------------------------------------------------ */

void readmap (char *name, int NZ, Mapping *mapx, Mapping *mapy)
{
  char  fname[FILENAME_MAX];
  FILE *fp;
  
  sprintf (fname, "%s.map", name);
  if ((fp = fopen (fname,"r")) == (FILE *) NULL)
    ErrorHandler ("readmap", "failed to open the map file", ERROR);
  
  readMap (fp, mapx);
  readMap (fp, mapy);
  return;
}
#endif

#endif //end of map
/* ------------------------------------------------------------------------ *
 * filter() - n-point linear filter (2-d only)
 * ------------------------------------------------------------------------ */

void filter (Mapping *map, char type)
{
  int            n = iparam ("NFILTER");
  static double *x, *y;

  register int       i = n;

  if (n == 0) return;

  /* out[i] = 1/n * (in[i] + in[i-1] + in[i-2] + ... + in[i-n+1]) */

  switch (type) {
  case 'x':
    if (!x) {
      printf ("filter: using filter of %d points\n", n);
      x = dvector (0, n-1);
      dfill (n, map->f[0], x, 1);
    }
    while (--i)
      x[i] = x[i-1];
    x[0] = map->f[0];
    map->f[0] = dsum(n, x, 1) / n;
    break;
  case 'y':
    if (!y) {
      y = dvector (0, n-1);
      dfill (n, map->f[0], y, 1);
    }
    while (--i)
      y[i] = y[i-1];
    y[0] = map->f[0];
    map->f[0] = dsum(n, y, 1) / n;
    break;
	
  default:
    puts ("no such case");
    break;
  }

  return;
}

/*
 * Generate a Z-mesh 
 */

double *Zmesh (int nz)
{
  double *z = dvector(0, nz),
       zmin = 0.0,
         lz = dparam("LZ"),
         dz = lz/nz;

  dramp (nz, &zmin, &dz, z, 1);
  z[nz] = zmin + lz;

  return z;
}

static char *hdr_fmt[] = { 
  "%-25s "            "Session\n",
  "%-25s "            "Created\n",
  "%-5c Hybrid              " "State 'p' = physical, 't' transformed\n",
  "%-5d %-5d %-5d         "   "Number of Elements; Dim of run; Lmax\n",
  "%-25d "            "Step\n",
  "%-25.6g "          "Time\n",
  "%-25.6g "          "Time step\n",
  "%-25.6g "          "Kinvis;\n",
  "%-25s "            "Fields Written\n",
  "%-25s "            "Format\n"
  };

static char *fourier_hdr_fmt[] = { 
  "%-25s "            "Session\n",
  "%-25s "            "Created\n",
  "%-5c Fourier Hybrid        " "State 'p' = physical, 't' transformed\n",
  "%-5d %-5d %-5d %-5d   "      "Number of Elements; Dim of run; Lmax; NZ\n",
  "%-25d "            "Step\n",
  "%-25.6g "          "Time\n",
  "%-25.6g "          "Time step\n",
  "%-11.6g %-11.6g   " "Kinvis; LZ\n",
  "%-25s "            "Fields Written\n",
  "%-25s "            "Format\n"
  };

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

void writeheader (FILE *fp, char *name, int lmax, double lz, int nz, 
			 int nel, int step, double t, char *typelist)
{
#define DESCRIP   25
#define BINARY   strings[0]        /* 0. Binary format file            */ 
#define ASCII    strings[1]        /* 1. ASCII format file             */ 
#define KINVIS   strings[2]        /* 2. Kinematic viscosity           */ 
#define DELT     strings[3]        /* 3. Time step                     */ 
#define CHECKPT  strings[4]        /* 4. Check-pointing option         */
#define BINFMT   strings[5]        /* 5. Binary format string          */

static char *strings [] = { 
  "binary", "ascii", "KINVIS", "DELT", "checkpt",
#ifdef BTYPE
  "binary-"BTYPE
#else
  "binary"
#endif
};
  time_t   tp;
  char     buf[128];

  /* set up parameters */

  tp = time((time_t*) NULL);
  strftime (buf, DESCRIP, "%a %b %d %H:%M:%S %Y", localtime(&tp));

  /* write out the file header */
  char **Hdr_fmt = fourier_hdr_fmt;

  fprintf (fp, Hdr_fmt[0], name);
  fprintf (fp, Hdr_fmt[1], buf);
  fprintf (fp, Hdr_fmt[2], 't');
  fprintf (fp, Hdr_fmt[3], nel, 2, lmax, nz);
  fprintf (fp, Hdr_fmt[4], step);
  fprintf (fp, Hdr_fmt[5], t);
  fprintf (fp, Hdr_fmt[6], dparam(DELT));
  fprintf (fp, Hdr_fmt[7], dparam(KINVIS), lz);
  fprintf (fp, Hdr_fmt[8], typelist);
  fprintf (fp, Hdr_fmt[9], ASCII);

  fflush (fp);

  return;
}

void ErrorHandler(char *section, char *message, FLAG type)
{
  char buf [BUFSIZ];

  switch (type) {
  case WARNING:
    sprintf (buf, "\n"
             "* * * * *  WARNING  * * * * *\n"
             "Section : %s\n"
             "Reason  : %s\n", section, message);
    sprintf (buf + strlen(buf), "node    : %d\n", option("PROCID"));
    fputs   (buf, stderr);
    break;

  case ERROR:
    sprintf (buf, "\n"
             "* * * * *  ERROR  * * * * *\n"
             "Section : %s\n"
             "Reason  : %s\n", section, message);
    sprintf (buf + strlen(buf), "node    : %d\n", option("PROCID"));
    fputs   (buf, stderr);
    exit    (-1);
    break;
  }
  return;
}


void four1(int nn, double *data, int isign);

/* ------------------------------------------------------------------------ *
 * orig_realft() - (Numerical Recipes') FFT of a single real function       *
 *                                                                          *
 * The real component of the first mode and the real component of the last  *
 * mode are packed into the first two elements of data.  This convention is *
 * expected regardless of the direction of the transform.                   *
 *                                                                          *
 * ------------------------------------------------------------------------ */

void orig_realft(double *data, int n, int isign)
{
  int i,i1,i2,i3,i4,n2p3;
  double c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=M_PI/(double) n;
  if (isign == 1) {                        /* FORWARD TRANSFORM */
    c2 = -0.5;
    four1(n,data,1);
  } else {                                 /* INVERSE TRANFORM */
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  n2p3=2*n+3;
  for (i=2;i<=n/2;i++) {           /* case i=1 done separately below */
    i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
    /* The two separate transforms are generated out of the data */
    h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    /* squeeze the first and last data elements together */
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    /* squeeze the first and last data elements together */
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(n,data,-1);
  }
}

void four1(int nn, double *data, int isign)
{
  int    n,mmax,m,j,istep,i;
  double wtemp, wr, wpr, wpi, wi, theta;  
  double tempr,tempi;

  n = nn << 1;    
  j = 1;

  /* 
   * BIT REVERSEAL SECTION
   */

  for(i = 1;i < n; i += 2) { 
    if (j > i) {                       /* Exchange the two complex numbers */
      tempr     = data[j];
      tempi     = data[j+1];
      data[j]   = data[i];
      data[j+1] = data[i+1];
      data[i]   = tempr;
      data[i+1] = tempi;
    }
    m = n >> 1;
    while(m >= 2 && j > m) {
      j  -= m;
      m >>= 1;
    }
    j += m;
  }

  /* 
   * DANIELSON-LANCZOS SECTION            
   */

  mmax = 2;
  while(n > mmax) {                 /* Outer loop executed log_2(nn) times */

    istep = mmax << 1;
    theta = M_PI * 2. / ( isign * mmax );
    wtemp = sin(0.5 * theta); 
    wpr   = -2. * wtemp * wtemp;
    wpi   = sin( theta );
    wr    = 1.;
    wi    = 0.;


    /*         FFT INNER NESTED LOOPS         *
     * -------------------------------------- */

    for(m = 1;m < mmax;m += 2) { 
        for(i = m;i <= n;i += istep) {
            j          = i + mmax;           /* Danielson-Lanczos formula */
            tempr      = wr * data[j]   - wi * data[j+1];
            tempi      = wr * data[j+1] + wi * data[j];
            data[j]    = data[i]   - tempr;
            data[j+1]  = data[i+1] - tempi;
            data[i]   += tempr;
            data[i+1] += tempi;
	}                                        /* trigonometric recurrence */
        wr = (wtemp=wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
      }

    /* -------------------------------------- */

    mmax = istep;
  }
  return;
}


/* ------------------------------------------------------------------------ *
 * sinft() - Discrete Fourier Sine Transform                                *
 *                                                                          *
 * Replaces y[] by its discrete Fourier sine transform.                     *
 *                                                                          *
 * ------------------------------------------------------------------------ */

void sinft(int n, double *y, int isign)
{
  int j,m=n/2,n2=n+2;
  double sum,y1,y2;
  double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;
  
  y     -=  1;      /* realign the data vector */
  isign *= -1;      /* switch the sign         */
  
  theta=M_PI/(double) n;
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  y[1]=0.0;
  for (j=2;j<=m+1;j++) {
    wr=(wtemp=wr)*wpr-wi*wpi+wr;    /* calculate the sine for the aux array */
    wi=wi*wpr+wtemp*wpi+wi;      /* the cosine is needed for the recurrence */
    y1=wi*(y[j]+y[n2-j]);                        /* construct the aux array */
    y2=0.5*(y[j]-y[n2-j]);
    y[j]=y1+y2;                              /* terms j and n-j are related */
    y[n2-j]=y1-y2;
  }
  orig_realft(y,m,1);                            /* transform the aux array */
  y[1]*=0.5;                 /* initialize the sum used for odd terms below */
  sum=y[2]=0.0;
  for (j=1;j<=n-1;j+=2) {
    sum += y[j];
    y[j]=y[j+1];                          /* even terms determined directly */
    y[j+1]=sum;                 /* odd terms determined by this running sum */
  }
  if (isign == 1) dscal(n, 2.0/n, y+1, 1);                       /* scaling */
}

/* ------------------------------------------------------------------------ *
 * cosft() - Discrete Fourier Cosine Transform                              *
 *                                                                          *
 * Replaces y[] by its discrete Fourier cosine transform.                   *
 *                                                                          *
 * ------------------------------------------------------------------------ */

void cosft(int n, double *y, int isign)
{
  int j,m,n2;
  double enf0,even,odd,sum,sume,sumo,y1,y2;
  double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;
  
  y     -=  1;      /* realign the data vector */
  isign *= -1;      /* switch the sign         */
  
  theta=M_PI/(double) n;
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  sum=y[1];
  m=n >> 1;
  n2=n+2;
  for (j=2;j<=m;j++) {
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
    y1=0.5*(y[j]+y[n2-j]);
    y2=(y[j]-y[n2-j]);
    y[j]=y1-wi*y2;
    y[n2-j]=y1+wi*y2;
    sum += wr*y2;
  }
  orig_realft(y,m,1);
  y[2]=sum;
  for (j=4;j<=n;j+=2) {
    sum += y[j];
    y[j]=sum;
  }
  if (isign == -1) {
    even=y[1];
    odd=y[2];
    for (j=3;j<=n-1;j+=2) {
      even += y[j];
      odd += y[j+1];
    }
    enf0=2.0*(even-odd);
    sumo=y[1]-enf0;
    sume=(2.0*odd/n)-sumo;
    y[1]=0.5*enf0;
    y[2] -= sume;
    for (j=3;j<=n-1;j+=2) {
      y[j] -= sumo;
      y[j+1] -= sume;
    }
  } else dscal(n, 2.0/n, y+1, 1);                                /* scaling */
}


