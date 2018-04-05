#==============================================================================
# Macros file for DMI cray XT5, intel compiler
#==============================================================================

INCLDIR := -I. -I/usr/include
SLIBS   :=

#--- Compiler/preprocessor ---
FC         := ftn
CC         := cc
CPP        := /usr/bin/cpp
CPPFLAGS   := -P -traditional  # ALLOW fortran double backslash "\\" 

#--- Flags ---
CFLAGS     := -c -O2
FFLAGS     := -O2 -ffree-line-length-none -finit-real=nan -g -fcheck=all -Wall -fbacktrace

# Preprocessor flags
CPPDEFS    := -DLINUX $(ICE_CPPDEFS)
CPPDEFS    := -DLINUX -DDMI -DDMI_alloc -DDMI_Sanity -DDMI_TEST $(ICE_CPPDEFS)
CPPDEFS    := -DLINUX -DDMI -DDMI_alloc -DDMI_Sanity -DDMI_TEST $(ICE_CPPDEFS)
CPPDEFS    := -DLINUX -DDMI_alloc -DDMI_nml $(ICE_CPPDEFS)
CPPDEFS    := -DLINUX -DDMI_alloc -DDMI_nml -DDMI_Sanity -DDMI_Sanity_MHRI $(ICE_CPPDEFS)

# Linker flags
LD         := $(FC)
LDFLAGS    := $(FFLAGS) -v

# Additional flags
FIXEDFLAGS := -132
FREEFLAGS  := -FR

#--- TEST flags ---
#ifeq ($(ICE_BLDDEBUG), true)
#  FFLAGS     += -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created
#else
#  FFLAGS     += -O2
#endif


ifeq ($(DITTO), yes)
   CPPDEFS :=  $(CPPDEFS) -DREPRODUCIBLE
endif


#--- OpenMP flags ---
ifeq ($(ICE_THREADED), true) 
   LDFLAGS += -fopenmp
   CFLAGS += -fopenmp
   FFLAGS += -fopenmp
endif

#--- NetCDF ---
ifeq ($(IO_TYPE), netcdf)
   CPPDEFS :=  $(CPPDEFS) -Dncdf
endif

ifeq ($(IO_TYPE), netcdf_bin)
   CPPDEFS :=  $(CPPDEFS) -Dncdf
endif

### if using parallel I/O, load all 3 libraries.  PIO must be first!
ifeq ($(ICE_IOTYPE), pio)
   PIO_PATH:=/usr/projects/climate/SHARED_CLIMATE/software/conejo/pio/1.7.2/intel-13.0.1/openmpi-1.6.3/netcdf-3.6.3-parallel-netcdf-1.3.1/include
   INCLDIR += -I$(PIO_PATH)
   SLIBS   := $(SLIBS) -L$(PIO_PATH) -lpio
endif

