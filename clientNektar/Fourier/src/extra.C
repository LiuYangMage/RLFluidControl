#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <polylib.h>
#include "veclib.h"
#include "hotel.h"
#include "nekstruct.h"
#include "Tri.h"
#include "Quad.h"
#include "nektarF.h"
#ifdef MAP
#include "map.h"
#endif

/* pack from groups of nztot into groups of nz */
void FtPpackp (int nz, int np, double *data, double *output){
  register int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  for(i = 0;i < nxy; ++i)
    for(j = 0; j < nprocs; ++j) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++) 
	output[i*nz+j*block+SLn] = data[i*nztot+j*nz+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nztot+j*nz, *SL_B = output+i*nz+j*block; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nztot+j*nz,output+i*nz+j*block);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nztot+j*nz,1,output+i*nz+j*block,1);
#endif // ifdef SL
    }
  return;
}

void PtFpackp (int nz, int np, double *data, double *output){
  register int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  if (option("dealias")) nztot = 3*nztot/2;
  for(i = 0;i < nxy; ++i)
    for(j = 0; j < nprocs; ++j) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++)
	output[i*nz+j*block+SLn] = data[i*nztot+j*nz+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nztot+j*nz, *SL_B = output+i*nz+j*block; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nztot+j*nz,output+i*nz+j*block);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nztot+j*nz,1,output+i*nz+j*block,1);
#endif // ifdef SL
    }
  return;
}

/* unpack vector of local groups of nz and put into groups of nztot */
void FtPunpackp (int nz, int np, double *data, double *output){
  int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  if (option("dealias")) {
    nztot = 3*nztot/2;
    dzero (nxy*nztot, output, 1);
  }
  for(j = 0; j < nprocs; ++j)
    for(i = 0;i < nxy; ++i) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++)
	output[i*nztot+j*nz+SLn] = data[i*nz+j*block+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nz+j*block, *SL_B = output+i*nztot+j*nz; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nz+j*block,output+i*nztot+j*nz);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nz+j*block,1,output+i*nztot+j*nz,1);
#endif // ifdef SL
    }
  return;
}

void PtFunpackp (int nz, int np, double *data, double *output){
  register int i,j;
  int nprocs = option("NPROCS");
  int nztot  = nprocs*nz;
  int nxy    = np/nprocs;
  int block  = nxy*nz;

  for(j = 0; j < nprocs; ++j)
    for(i = 0;i < nxy; ++i) {
#ifdef SL
#ifndef TRY_OUT
      for (register int SLn = 0; SLn < nz; SLn++)
	output[i*nztot+j*nz+SLn] = data[i*nz+j*block+SLn];
#else // ifndef TRY_OUT
#ifdef INLINE_SL
      { register int SLn; double *SL_A = data+i*nz+j*block, *SL_B = output+i*nztot+j*nz; for (SLn = 0; SLn < nz; SLn++) SL_B[SLn] = SL_A[SLn]; }
#else // ifdef INLINE_SL
      SLcopy(nz, data+i*nz+j*block,output+i*nztot+j*nz);
#endif // ifdef INLINE_SL
#endif // ifndef TRY_OUT
#else // ifdef SL
      dcopy(nz, data+i*nz+j*block,1,output+i*nztot+j*nz,1);
#endif // ifdef SL
    }
  return;
}
