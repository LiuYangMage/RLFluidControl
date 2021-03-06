message = Extra build flags:\n\n
#COMPATH        = /opt/xt-pe/1.3.10/bin/snos64/
CXX     = CC
CC      = cc
FC      = ftn
ifdef MPIPRG
# CFM default rule
message += CFM default uses LAM,\n to use MPICH define MPICH.\n\n
ifdef MPICH
MPICXX  = CC
MPICC   = cc
MPILIB  =
else
MPICXX  = CC
MPICC   = cc
MPILIB  =
endif
endif
LLAPACK = -llapack #-lsci
LBLAS   = -lblas
LPDIO    = #-L/usr/local/packages/pdio/nektar/lib -lpdio
INCIOBUF = #-I${IOBUF_INC}
LIOBUF   = #-liobuf


F90LIB = -lmpichf90 -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl



message += You can specify MCPU=-mcpu=pentiumpro,ev5 etc.\n or MCPU=-march=pentium. See the gcc info pages for more information.
OPTXX   = -O2 #$(MCPU)  -fastsse   -malign-double -ffast-math
MOPTXX  = -g -O3  -fastsse -Munroll -tp barcelona-64 #-tp k8-64 $(INCIOBUF)  #$(MCPU)  -fastsse   -malign-double -ffast-math -funroll-loops -fstrict-aliasing
OPT     = $(OPTXX)
MOPT    = $(MOPTXX)
OPTF    = -O2 #$(MCPU) -malign-double -ffast-math # -fno-globals
MOPTF   =  -g -O3 -fastsse -Munroll -tp barcelona-64 #-tp k8-64 #$(MCPU) -malign-double -ffast-math -funroll-loops # -fno-globals
DEBUG   = -g -Mprof=func
DEBUGF  = -g -Mprof=func
LD      = $(CXX)
DEBUGLD = $(DEBUG)
MOPTLD  = $(MOPTXX)
OPTLD   = $(OPTXX)
ifdef FFTPRG
FFTLIB  = -lfftw
endif
#LPDIO    = #-L/usr/local/packages/pdio/nektar/lib -lpdio
#INCIOBUF = -I/${IOBUF_INC}
#LIOBUF  = -liobuf

LIBS    = $(INCIOBUF)  $(FFTLIB) $(MPILIB) $(LLAPACK) $(LBLAS) $(F90LIB) -lvec  $(LPDIO) $(LIOBUF)
#LIBS   = $(FFTLIB) $(MPILIB) $(LLAPACK) $(LBLAS) -lvec
XLIBDIRS = -L/usr/X11/lib
XINCDIRS = -I/usr/include/X11

