#-------------------------------------------------------------------------------
# Common Makefile for CICE
#
# Command-line variables
#   MACFILE=<file> ~ the macros definition file to use/include
#   EXEC=<name>    ~ name given to executable, default is a.out
#   VPATH=<vpath>  ~ VPATH            , default is . (cwd only)
#   SRCS=<files>   ~ list of src files, default is all .c .F .F90 files in VPATH
#   VPFILE=<file>  ~ file with list of dirs, used to create VPATH
#   SRCFILE=<file> ~ file with list of src files, used to create SRCS
#
#   <macro defns>  ~ any macro definitions found in this file or the included 
#                    MACFILE will be over-riden by cmd-line macro definitions
#   MODEL=<model>  ~ a standard macro definition, often found in the included 
#                    MACFILE, used to trigger special compilation flags
#
# Usage examples:
#   % gmake -j 8 VPFILE=Filepath EXEC=${ICE_RUNDIR}/cice \
#      -f  ${ICE_CASEDIR}/Makefile MACFILE=${ICE_CASEDIR}/Macros.conrad_intel \
#      DEPFILE=${ICE_CASEDIR}/makdep.c cice
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# parse cmd-line and establish values for EXEC, VPATH, SRCS, OBJS, etc
#-------------------------------------------------------------------------------

EXEC     := a.out
MACFILE  := NONE 
DEPFILE  := NONE 
MODEL    := NONE 
VPFILE   := NONE
VPATH    := .
SRCFILE  := NONE
SRCS     := NONE

# dependency generator
DEPGEN   := ./makdep
OBJS_DEPGEN := $(DEPFILE)

ifneq ($(ESMFMKFILE),)
  -include $(ESMFMKFILE)
  INCLDIR += $(ESMF_F90COMPILEPATHS)
  SLIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

ifneq ($(VPATH),.)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(VPFILE),NONE)
    # explicit list of VPATH dirs is provided
    VPATH := $(wildcard . $(shell cat $(VPFILE) ) )
  endif
endif

ifneq ($(SRCS),NONE)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(SRCFILE),NONE)
    # explicit list of src files is provided
    SRCS := $(shell cat $(SRCFILE) )
  else
    # list of src files is all .F90 .F .c files in VPATH
    SRCS := $(wildcard $(addsuffix /*.F90 , $(VPATH)) \
		       $(addsuffix /*.[cF], $(VPATH)) )
  endif
endif

OBJS  := $(addsuffix .o, $(sort $(basename $(notdir $(SRCS)))))
DEPS  := $(addsuffix .d, $(sort $(basename $(notdir $(SRCS)))))
INCS  := $(patsubst %,-I%, $(VPATH) )
MODDIR:= -I.
RM    := rm
AR    := ar

.SUFFIXES:
.SUFFIXES: .F90 .F .c .o

.PHONY: all cice libcice targets target db_files db_flags clean realclean helloworld calchk
all: $(EXEC)

cice: $(EXEC)

#-------------------------------------------------------------------------------
# include the file that provides macro definitions required by build rules
#-------------------------------------------------------------------------------

-include $(MACFILE)

#-------------------------------------------------------------------------------
# echo supported targets, file names, paths, compile flags, etc. used during build
#-------------------------------------------------------------------------------

targets:
	@echo " "
	@echo "Supported Makefile Targets are: cice, libcice, makdep, depends, clean, realclean"
	@echo "                   Diagnostics: targets, db_files, db_flags"
	@echo "                   Unit Tests : helloworld, calchk"
target: targets

db_files:
	@echo " "
	@echo "* EXEC    := $(EXEC)"
	@echo "* MACFILE := $(MACFILE)"
	@echo "* VPFILE  := $(VPFILE)"
	@echo "* VPATH   := $(VPATH)"
	@echo "* SRCFILE := $(SRCFILE)"
	@echo "* INCS    := $(INCS)"
	@echo "* SRCS    := $(SRCS)"
	@echo "* OBJS    := $(OBJS)"
	@echo "* DEPS    := $(DEPS)"
	@echo "* ULIBS   := $(ULIBS)"
	@echo "* SLIBS   := $(SLIBS)"
	@echo "* INCLDIR := $(INCLDIR)"
	@echo "* DEPFILE := $(DEPFILE)"
	@echo "* OBJS_DEPGEN := $(OBJS_DEPGEN)"
db_flags:
	@echo " "
	@echo "* $(DEPGEN) := $(SCC) $(CFLAGS_HOST)"
	@echo "* %.d : %.c   := $(DEPGEN) $(INCS)"
	@echo "* %.d : %.F   := $(DEPGEN) $(INCS)"
	@echo "* %.d : %.F90 := $(DEPGEN) $(INCS)"
	@echo "* %.d : %.H   := $(DEPGEN) $(INCS)"
	@echo "* cpp     := $(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCLDIR)"
	@echo "* .c.o    := $(CC) $(CFLAGS) $(CPPDEFS) $(INCLDIR)"
	@echo "* .F.o    := $(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(CPPDEFS) $(INCLDIR)"
	@echo "* .F90.o  := $(FC) -c $(FFLAGS) $(FREEFLAGS) $(CPPDEFS) $(MODDIR) $(INCLDIR)"
	@echo "* libcice := $(AR) -r $(EXEC) "
	@echo "* $(notdir $(EXEC)) := $(LD) $(LDFLAGS) $(ULIBS) $(SLIBS)"

#-------------------------------------------------------------------------------
# build rule for makdep: MACFILE, cmd-line, or env vars must provide
# the needed macros
#-------------------------------------------------------------------------------

$(DEPGEN): $(OBJS_DEPGEN)
	@ echo "Building makdep"
	$(SCC) -o $@ $(CFLAGS_HOST) $<

#-------------------------------------------------------------------------------
# unit tests
#-------------------------------------------------------------------------------

# this builds all dependent source code automatically even though only a subset might actually be used
# this is no different than the cice target and in fact the binary is called cice
# it exists just to create separation as needed for unit tests
calchk: $(EXEC)

# this builds just a subset of source code specified explicitly and requires a separate target
HWOBJS := helloworld.o
helloworld: $(HWOBJS)
	$(LD) -o $(EXEC) $(LDFLAGS) $(HWOBJS) $(ULIBS) $(SLIBS)

#-------------------------------------------------------------------------------
# build rules: MACFILE, cmd-line, or env vars must provide the needed macros
#-------------------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(LD) -o $(EXEC) $(LDFLAGS) $(OBJS) $(ULIBS) $(SLIBS)

libcice: $(OBJS)
	@ echo "$(AR) -r $(EXEC) $(OBJS)"
	$(AR) -r $(EXEC) $(OBJS)

.c.o:
	$(CC) $(CFLAGS) $(CPPDEFS) $(INCLDIR) $<

.F.o:
	$(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(CPPDEFS) $(INCLDIR) $<

.F90.o:
	$(FC) -c $(FFLAGS) $(FREEFLAGS) $(CPPDEFS) $(MODDIR) $(INCLDIR) $<

clean:
	$(RM) -f *.o *.d *.mod $(EXEC)

realclean: clean
	$(RM) -f $(DEPGEN)

#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------
# ASSUMPTIONS:
# o the dependency generator, $(DEPGEN), can be built, 
#   its cmd line syntax is compatible with the build rules below.   Eg, for 
#   each .o file, there is a corresponding .d (dependency) file, and both
#   will be dependent on the same src file, eg.    foo.o foo.d : foo.F90
#   Also, the dependancy genorator's capabilities, limitations, and assumptions
#   are understood & accepted.
#-------------------------------------------------------------------------------

depends: $(DEPS)

%.d : %.c $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $< | head -3  > $@
%.d : %.F $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.F90 $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.H $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@

# the if-tests prevent DEPS files from being created when they're not needed
ifneq ($(MAKECMDGOALS), db_files)
ifneq ($(MAKECMDGOALS), db_flags)
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), realclean)
ifneq ($(MAKECMDGOALS), targets)
ifneq ($(MAKECMDGOALS), target)
ifneq ($(MAKECMDGOALS), makdep)
ifneq ($(MAKECMDGOALS), depends)
    -include $(DEPS)
endif
endif
endif
endif
endif
endif
endif
endif
