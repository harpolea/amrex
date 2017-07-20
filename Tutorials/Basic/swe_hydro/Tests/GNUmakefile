AMREX_HOME ?= ../../../../../amrex

ADR_DIR  ?= $(AMREX_HOME)/Tutorials/Basic/swe_hydro

TOP := $(ADR_DIR)

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
#DEBUG      = FALSE

USE_PARTICLES = TRUE

#DIM        = 2
DIM       = 3

COMP	   = g++

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_IPM      = FALSE
PROFILE      = TRUE
TRACE_PROFILE = TRUE
COMM_PROFILE = TRUE
USE_DUMPI    = FALSE
PRECISION    = DOUBLE
EBASE        = test


CEXE_sources	+= test.cpp
f90EXE_sources	+= test_f.f90 #$(AMREX_HOME)/Tutorials/Basic/swe_hydro/Source/Src_3d/compute_flux_3d.f90

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include $(AMREX_HOME)/Src/Base/Make.package

Bdirs 	:= Source/Src_nd Source/Src_$(DIM)d
Bpack	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs   += $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)


Pdirs 	:= Base AmrCore Boundary Particle
Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
