COMPATH	= /usr/bin
CXX	= $(COMPATH)/g++
CC	= $(COMPATH)/gcc 
FC	= $(COMPATH)/gfortran

#CXX	= mpic++
#CC	= mpicc 
#FC	= mpif90


#MPILIB	= -L/gpfs/runtime/opt/openmpi/1.8/lib
#MPIINC	= -I/gpfs/runtime/opt/openmpi/1.8/include

ifdef PARALLEL
CXX  = mpic++
CC   = mpicc
FC   = mpif90
endif

OPTXX	= -O2 $(MCPU) -ffast-math 
MOPTXX	= -O3 $(MCPU) -ffast-math -funroll-loops -fstrict-aliasing 
OPT	= -O2 $(MCPU) -ffast-math
MOPT	= -O3 $(MCPU) -ffast-math -funroll-loops -fstrict-aliasing 
OPTF	= -O2 $(MCPU) -ffast-math
MOPTF	= -O3 $(MCPU) -ffast-math -funroll-loops
DEBUG	= -g
DEBUGF	= -g 
LD	= $(CXX)
DEBUGLD	= $(DEBUG)
MOPTLD	= $(MOPTXX)
OPTLD	= $(OPTXX)
ifdef FFTPRG
#FFTLIB	= -L/users/zwang/Applications/fftw_2.1.5_install/lib/ -lrfftw -lfftw -lm
FFTLIB	= -lfftw -lm
endif
#LLAPACK = -L/usr/lib64 -llapack 
#LBLAS   = -L/usr/lib64 -lblas
#LLAPACK = /usr/lib64/liblapack.so
#LBLAS   = /usr/lib64/libblas.so
LLAPACK = -L/users/zwang/Applications/lapack-3.6.0 -llapack
LBLAS   = -L/users/zwang/Applications/lapack-3.6.0 -lrefblas 
ifdef UMFPACKSLV
UMFPACKLIB = -lumfpack -lamd
endif
LG2C  = -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/ -lg2c
#LPG   	= -L/usr/local/lib -lpgftnrtl -lpgc   
#LIBS	= -lvec $(FFTLIB) $(UMFPACKLIB) $(LLAPACK) $(LBLAS) -lgfortran #$(LG2C) #-lg2c
LIBS = -lvec $(FFTLIB) -L$(MKLROOT)/lib/intel64 -lmkl_core -lmkl_sequential -lmkl_intel_lp64 -lsvml -lifcore -lgfortran
