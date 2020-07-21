ARCH	:= $(shell uname)
PROG	= nektar
HYBRID	= ../../
VPATH   = ../src
HLIBDIR	= $(HYBRID)/Hlib/$(ARCH)
HDRDIR	= $(HYBRID)/include
TARFILE	= $(PROG).tar

ALL     = mpinektar drive io prepost bwoptim convective comm\
	  forces pressure rotational stokes wannier analyser\
	  mlevel varv convarv Fourier_List interp cfl viscosity postprocess nse_util # stokes_solve_F

defines = DIM=2  _ALGOBASE_H TIMING FORCES DYNAMIC_ALPHA #SVV_KERNEL #NS_RESIDUAL  #CONJUGATE_HEAT  #BETA0  #OLDFFTS
CPPFLAGS = -I$(HDRDIR)
FFTPRG	= yes
#MPIPRG  = yes

ifeq ($(ARCH),SunOS)			# SunOS 5.*
include $(HYBRID)/Flags/SunOS.inc
#OPT += -Qoption iropt -O4+vector
endif

ifeq ($(ARCH),ap3000)			# ap3000
include $(HYBRID)/Flags/ap3000.inc
endif

ifeq ($(ARCH),IRIX)			# Silicon Graphics IRIX
include $(HYBRID)/Flags/IRIX.inc
endif

ifeq ($(ARCH),IRIX64)			# Silicon Graphics IRIX64
include $(HYBRID)/Flags/IRIX64.inc
endif

ifeq ($(ARCH),AIX)			# IBM AIX
include $(HYBRID)/Flags/AIX.inc
endif

ifeq ($(ARCH),cray-pvp)			# Cray PVP
include $(HYBRID)/Flags/CRAYPVP.inc
endif

ifeq ($(ARCH),cray-t3d)			# Cray T3D
include $(HYBRID)/Flags/T3D.inc
endif

ifeq ($(ARCH),unicosmk)			# Cray T3E
include $(HYBRID)/Flags/T3E.inc
endif

ifeq ($(ARCH),alpha)			# DEC Alpha AXP
include $(HYBRID)/Flags/DECAXP.inc
endif

ifeq ($(ARCH),HP-UX)			# HP-UX
include $(HYBRID)/Flags/HP-UX.inc
endif

ifeq ($(ARCH),GCC)			# Generic GCC platform
include $(HYBRID)/Flags/GCC.inc
endif

ifeq ($(ARCH),Linux)			# Linux using  GCC 
include $(HYBRID)/Flags/Linux.inc
endif

ifeq ($(ARCH),Darwin)			# Max OSX
include ../../Flags/Darwin.inc
endif

ifdef PARALLEL
override CPPFLAGS :=  -DPARALLEL  $(CPPFLAGS) 
endif

ifeq ($(MAP), 1)
override defines := $(defines) MAP
override ALL	 := flex map $(ALL)
override PROG	 := flex
endif

ifeq ($(TEST), 1)
override defines := $(defines) TEST MAP
override ALL	 := flex map rea $(ALL)
endif

ifdef    defines
DEF     = $(foreach def, $(defines) $(ADDONS), -D$(def))
override CPPFLAGS :=  $(DEF) $(CPPFLAGS)
endif

ifeq ($(THERMO), 1)
override CPPFLAGS := -DTHERMO $(CPPFLAGS)
override ALL	 := temperature_viscosity $(ALL)
override PROG	 := nektar_thermo
endif

ifeq ($(NATURAL_CONVECTION), 1)
override CPPFLAGS := -DNATURAL_CONVECTION $(CPPFLAGS)
override ALL	 := temperature_etc $(ALL)
endif

ifeq ($(SPM), 1)
override CPPFLAGS := -DSPM $(CPPFLAGS)
override ALL	 := C_SPM C_SHAPE Nektar2SPM $(ALL)
override PROG	 := nektar_SPM
endif

ifeq ($(DISC_FORCING), 1)
override CPPFLAGS :=  -DDISC_FORCING $(CPPFLAGS) 
endif

ifeq ($(WALE), 1)
override CPPFLAGS := -DWALE $(CPPFLAGS)
endif

ifeq ($(EVM), 1)
override CPPFLAGS := -DENTROPYVISCOSITY $(CPPFLAGS)
endif

ifeq ($(OUTPUT_VISC), 1)
override CPPFLAGS :=  -DOUTPUT_VISCOSITY $(CPPFLAGS) 
endif

ifeq ($(EVM_NOT_APPLIED), 1)
override CPPFLAGS := -DEVM_NOT_APPLIED $(CPPFLAGS)
endif

ifeq ($(XMLRPC), 1)
override CPPFLAGS := -DXMLRPC  -I/gpfs_home/zwang197/xmlrpc-c-1.51.06/build/include  $(CPPFLAGS)
override ALL	 := xmlrpc_c $(ALL)
endif

export

help:
	@echo ''
	@echo ' Nektar Makefiles:'
	@echo ''
	@echo ' Run make with one of the following options:'	
	@echo '	a) dbx  ...Fourier dbx version'
	@echo '	b) opt  ...Fourier optimised version'
	@echo '	c) mopt ...Fourier maximally (?) optimised version'
	@echo ''
	@echo '	need PARALLEL=1 in command line for parallel version'
	@echo ''

dbx: 
	$(MAKE) LTYPE=g -f MakeNek
opt:
	$(MAKE) LTYPE=opt -f MakeNek
mopt:
	$(MAKE) LTYPE=mopt -f MakeNek

tar: Makefile
	@echo 'creating a tape archive in' $(TARFILE)
	@tar cvf $(TARFILE) Makefile MakeNek src/*.C
clean:
	$(MAKE) -f MakeNek clean

