
CXX  = mpic++
CC   = mpicc
FC   = mpif77

MPICXX  = mpic++
MPICC   = mpicc

MPILIB	= -L/gpfs/runtime/opt/openmpi/1.8/lib
MPIINC	= -I/gpfs/runtime/opt/openmpi/1.8/include

#LLAPACK = -L/opt/intel/mkl/10.0.3.020/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -lmkl_core
#LLAPACK = -L/usr/local/lib -lscalapack 
#LBLAS   = -lm 
LBLAS   = -L/users/zwang/Applications/lapack-3.6.0 -lrefblas 

#LLAPACK = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack 
#LLAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lirc #-lguide
#LLAPACK = -L/gpfs/runtime/opt/intel/2013.1.106/mkl/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
LLAPACK = -L/users/zwang/Applications/lapack-3.6.0 -llapack

OPTXX	= -O3 
MOPTXX	= -O3 
OPT	= $(OPTXX)
MOPT	= $(MOPTXX)
OPTF	= -O3 
MOPTF	= -O3 
DEBUG	= -g 
DEBUGF	= -g

LD	= $(CXX)
DEBUGLD	= $(DEBUG)
MOPTLD	= $(MOPTXX)
OPTLD	= $(OPTXX)

#FFTLIB	= -L/users/zwang/Applications/fftw_2.1.5_install/lib/ -lrfftw -lfftw -lm
FFTLIB	= -lfftw -lm

LIBS	= -lvec -lfftw $(FFTLIB) $(MPILIB) $(LLAPACK) $(LBLAS) -lgfortran #-lsvml -limf -lirc


