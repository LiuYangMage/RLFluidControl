//**********************************************************************//
// DMath class
// common routines
//**********************************************************************//

#ifdef __FFTW__
#include "fftw.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

#ifndef __DMATH_H__
#define __DMATH_H__

class DMath
{
 public:
  DMath() {};
  ~DMath() {};

 public:
  static int checkEndian()
  {
    // return value: -1 ---> error in function call
    //                0 ---> little endian
    //                1 ---> big endian
    //                2 ---> unknown endian type
    int num = 0x0102;
    int len_int = sizeof(int);
    if(len_int<2) return -1;
    char *ptr = (char*)&num;
    
    if(ptr[0]==0x02 && ptr[1]==0x01) return 0; // little endian
    else if(ptr[len_int-1]==0x02 && ptr[len_int-2]==0x01) return 1; // big endian
    else return 2; // unknown endian type
  }

 public:
  static double* newD(int dim)
  {
    if(dim<=0) return NULL;

    double *temp = NULL;
    temp = new double[dim];
    return temp;
  };
  static double** newD(int dimx, int dimy)
  {
    if(dimx<=0 || dimy<=0) return NULL;

    double **temp = NULL;
    temp = new double*[dimx];
    *temp = new double[dimx*dimy];
    for(int i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;

    return temp;
  };
  static double*** newD(int dimx, int dimy, int dimz)
  {
    if(dimx<=0 || dimy<=0 || dimz<=0) return NULL;

    double ***temp = NULL;
    temp = new double**[dimx];
    *temp = new double*[dimx*dimy];
    **temp = new double[dimx*dimy*dimz];

    int i, j;
    for(j=0;j<dimx*dimy;j++)
      (*temp)[j] = **temp + j*dimz;
    for(i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;

    return temp;
  };
  static void del(double *a)
  {
    if(!a) return;
    delete []a;
  };
  static void del(double** a)
  {
    if(!a) return;
    delete []*a;
    delete []a;
  };
  static inline void del(double ***a)
  {
    if(!a) return;
    delete []**a;
    delete []*a;
    delete []a;
  };

 public:
  static int* newI(int dim)
  {
    if(dim<=0) return NULL;

    int *temp = NULL;
    temp = new int[dim];
    return temp;
  };
  static int** newI(int dimx, int dimy)
  {
    if(dimx<=0 || dimy<=0) return NULL;

    int **temp = NULL;
    temp = new int*[dimx];
    *temp = new int[dimx*dimy];
    for(int i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;

    return temp;
  };
  static int*** newI(int dimx, int dimy, int dimz)
  {
    if(dimx<=0 || dimy<=0 || dimz<=0) return NULL;

    int ***temp = NULL;
    temp = new int**[dimx];
    *temp = new int*[dimx*dimy];
    **temp = new int[dimx*dimy*dimz];

    int i,j;
    for(j=0;j<dimx*dimy;j++)
      (*temp)[j] = **temp + j*dimz;
    for(i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;

    return temp;
  };
  static void del(int* a)
  {
    if(!a) return;
    delete []a;
  };
  static void del(int **a)
  {
    if(!a) return;
    delete []*a;
    delete []a;
  };
  static void del(int ***a)
  {
    if(!a) return;
    delete []**a;
    delete []*a;
    delete []a;
  };

 public:
  static char* newCH(int dim)
  {
    if(dim<=0) return NULL;

    char *temp = NULL;
    temp = new char[dim];
    return temp;
  };
  static char** newCH(int dimx, int dimy)
  {
    if(dimx<=0 || dimy<=0) return NULL;

    char **temp = NULL;
    temp = new char*[dimx];
    *temp = new char[dimx*dimy];
    for(int i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;

    return temp;
  };
  static char*** newCH(int dimx, int dimy, int dimz)
  {
    if(dimx<=0 || dimy<=0 || dimz<=0) return NULL;

    char ***temp = NULL;
    temp = new char**[dimx];
    *temp = new char*[dimx*dimy];
    **temp = new char[dimx*dimy*dimz];

    int i,j;
    for(j=0;j<dimx*dimy;j++)
      (*temp)[j] = **temp + j*dimz;
    for(i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;

    return temp;
  };
  static void del(char* a)
  {
    if(!a) return;
    delete []a;
  };
  static void del(char **a)
  {
    if(!a) return;
    delete []*a;
    delete []a;
  };
  static void del(char ***a)
  {
    if(!a) return;
    delete []**a;
    delete []*a;
    delete []a;
  };  

 public:
#ifdef __FFTW__
  static FFTW_COMPLEX* newC(int dim)
  {
    if(dim<=0) return NULL;

    FFTW_COMPLEX *temp = NULL;
    temp = new FFTW_COMPLEX[dim];
    return temp;
  };
  static FFTW_COMPLEX** newC(int dimx, int dimy)
  {
    if(dimx<=0 || dimy<=0) return NULL;

    FFTW_COMPLEX **temp = NULL;
    temp = new FFTW_COMPLEX*[dimx];
    *temp = new FFTW_COMPLEX[dimx*dimy];

    for(int i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;
    return temp;
  };
  static FFTW_COMPLEX*** newC(int dimx, int dimy, int dimz)
  {
    if(dimx<=0 || dimy<=0 || dimz<=0) return NULL;

    FFTW_COMPLEX ***temp = NULL;
    temp = new FFTW_COMPLEX**[dimx];
    *temp = new FFTW_COMPLEX*[dimx*dimy];
    **temp = new FFTW_COMPLEX[dimx*dimy*dimz];

    int i,j;
    for(j=0;j<dimx*dimy;j++)
      (*temp)[j] = **temp + j*dimz;
    for(i=0;i<dimx;i++)
      temp[i] = *temp + i*dimy;

    return temp;
  };
  static void del(FFTW_COMPLEX *a)
  {
    if(!a) return;
    delete []a;
  };
  static inline void del(FFTW_COMPLEX **a)
  {
    if(!a) return;
    delete []*a;
    delete []a;
  };
  static void del(FFTW_COMPLEX ***a)
  {
    if(!a) return;
    delete []**a;
    delete []*a;
    delete []a;
  };
#endif // __FFTW__

};

#endif // __DMATH_H__
