#################################
### Choose whether to use LHAPDF
### LHAPDF = on,off
LHAPDF = on
LHADIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/lhapdf/6.1.6-giojec/lib
#################################

#################################
FCC = gfortran
CPP=g++
CC=gcc
#################################

ifeq (${FCC},g77)
LIBS=-lg2c -lm
CCLINK=$(CPP)
FFLAGS = -ffixed-line-length-none -w
endif

ifeq ($(FCC),gfortran)
LIBS=-lstdc++ -lm
CCLINK=$(FCC)
FFLAGS = -ffixed-line-length-none -w -O1
endif

ifeq ($(FCC),ifort)
LIBS=-lstdc++ 
CCLINK=$(FCC)
LDFLAGS=-nofor-main 
FFLAGS= -extend_source -fpp
endif

BASEDIR := $(shell pwd)

CUBADIR = $(BASEDIR)/Cuba-3.0
LIBS += -L$(CUBADIR) -lcuba

ifeq ($(LHAPDF),on)
LIBS += -L$(LHADIR) -lLHAPDF
LDFLAGS += -Wl,-rpath,$(LHADIR)
endif

PDFSRC = $(BASEDIR)/pdfSrc
SRCDIR = $(BASEDIR)/src
SRCDIR_Z = $(SRCDIR)/z
SRCDIR_W = $(SRCDIR)/w
EXECDIR = $(BASEDIR)/bin
PDFBUILD = $(BASEDIR)/pdfBuild
BUILDDIR = $(BASEDIR)/build
BUILDDIR_Z = $(BUILDDIR)/z
BUILDDIR_W = $(BUILDDIR)/w

INCLUDEDIR = $(BASEDIR)/include

# add include flag
FFLAGS += -I${INCLUDEDIR}

MAINSRC = $(BASEDIR)/mainSrc

MAINBUILD = $(BASEDIR)/mainBuild

EXEC_Z = $(EXECDIR)/fewzz
CONDOR_EXEC_Z = $(EXECDIR)/condor_fewzz
EXEC_W = $(EXECDIR)/fewzw
CONDOR_EXEC_W = $(EXECDIR)/condor_fewzw
targets = $(EXEC_Z) $(EXEC_W) $(CONDOR_EXEC_W) $(CONDOR_EXEC_Z)

.PHONY: force clean distclean cuba fewzz fewzw condor_fewzw condor_fewzz

all: cuba $(targets)

cuba: $(CUBADIR)/libcuba.a
fewzz: cuba $(EXEC_Z)
condor_fewzz: cuba $(CONDOR_EXEC_Z)
fewzw: cuba $(EXEC_W)
condor_fewzw: cuba $(CONDOR_EXEC_W)

$(BUILDDIR_Z)/%.o : $(SRCDIR_Z)/%.F
	@mkdir -p $(BUILDDIR_Z)
	$(FCC) $(FFLAGS)  -c $< -o $@
	
$(BUILDDIR_W)/%.o : $(SRCDIR_W)/%.F
	@mkdir -p $(BUILDDIR_W)
	$(FCC) $(FFLAGS)  -c $< -o $@

$(BUILDDIR)/%.o : $(SRCDIR)/%.F
	@mkdir -p $(BUILDDIR)
	$(FCC) $(FFLAGS)  -c $< -o $@

$(PDFBUILD)/%.o : $(PDFSRC)/%.F
	@mkdir -p $(PDFBUILD)
	$(FCC) $(FFLAGS)  -c $< -o $@

$(MAINBUILD)/%.o : $(MAINSRC)/%.F
	@mkdir -p $(MAINBUILD)
	$(FCC) $(FFLAGS) -c $< -o $@

#$(MAINBUILD)/%.o : $(MAINSRC)/%.cpp
#	@mkdir -p $(MAINBUILD)
#	$(CPP) -c $< -o $@

$(MAINBUILD)/fewzwrap.o : $(MAINSRC)/fewzwrap.cc
	@mkdir -p $(MAINBUILD)
	$(CPP) -c $(MAINSRC)/fewzwrap.cc -o $(MAINBUILD)/fewzwrap.o

$(MAINBUILD)/condor_fewzwrap.o : $(MAINSRC)/fewzwrap.cc
	@mkdir -p $(MAINBUILD)
	$(CPP) -c $(MAINSRC)/fewzwrap.cc -o $(MAINBUILD)/condor_fewzwrap.o

objects := $(patsubst $(SRCDIR)/%.F,$(BUILDDIR)/%.o,$(wildcard $(SRCDIR)/*.F))

ifneq ($(LHAPDF),on)
objects += $(patsubst $(PDFSRC)/%.F,$(PDFBUILD)/%.o,$(wildcard $(PDFSRC)/*.F))
endif

#execobj := $(patsubst $(MAINSRC)/%.cc,$(MAINBUILD)/%.o,$(wildcard $(MAINSRC)/*.cc))
#execobj += $(patsubst $(MAINSRC)/%.cpp,$(MAINBUILD)/%.o,$(wildcard $(MAINSRC)/*.cpp))
#execobj += $(patsubst $(MAINSRC)/%.F,$(MAINBUILD)/%.o,$(wildcard $(MAINSRC)/*.F))
execobj += $(MAINBUILD)/cuba.o
execobj += $(MAINBUILD)/FEWZ.o
ifeq ($(LHAPDF),on)
execobj += $(MAINBUILD)/InitPDFsLHA.o $(MAINBUILD)/ObtainPDFsLHA.o $(MAINBUILD)/fitpdf.o
else
execobj += $(MAINBUILD)/InitPDFs.o $(MAINBUILD)/ObtainPDFs.o $(MAINBUILD)/fitpdf.o
endif
execobj += ${objects}
condor_execobj := ${execobj}
condor_execobj += $(MAINBUILD)/condor_fewzwrap.o
execobj += $(MAINBUILD)/fewzwrap.o
execobj_z := ${execobj}
execobj_w := ${execobj}
z_obj := $(patsubst $(SRCDIR_Z)/%.F,$(BUILDDIR_Z)/%.o,$(wildcard $(SRCDIR_Z)/*.F))
z_obj += $(MAINBUILD)/FEWZz.o
w_obj := $(patsubst $(SRCDIR_W)/%.F,$(BUILDDIR_W)/%.o,$(wildcard $(SRCDIR_W)/*.F))
w_obj += $(MAINBUILD)/FEWZw.o
condor_execobj_z := ${condor_execobj}
condor_execobj_w := ${condor_execobj}
execobj_z += $(z_obj)
execobj_w += $(w_obj)
condor_execobj_z += $(z_obj)
condor_execobj_w += $(w_obj)

# EXECUTABLE

$(CUBADIR)/libcuba.a:
	cd $(CUBADIR) && ./configure && make lib

$(EXEC_Z): $(execobj_z) 
	@mkdir -p $(EXECDIR)
	$(CCLINK) $(LDFLAGS) -o $@ $(execobj_z) $(LIBS)

$(CONDOR_EXEC_Z): $(condor_execobj_z)
	@mkdir -p $(EXECDIR)
	condor_compile $(CCLINK) $(LDFLAGS) -o $@ $(condor_execobj_z) $(LIBS)

$(EXEC_W): $(execobj_w) 
	@mkdir -p $(EXECDIR)
	$(CCLINK) $(LDFLAGS) -o $@ $(execobj_w) $(LIBS)

$(CONDOR_EXEC_W): $(condor_execobj_w)
	@mkdir -p $(EXECDIR)
	condor_compile $(CCLINK) $(LDFLAGS) -o $@ $(condor_execobj_w) $(LIBS)

wipe:
	rm -rf $(MAINBUILD)

clean:
	rm -rf ${MAINBUILD} ${BUILDDIR} $(BUILDDIR_Z) $(BUILDDIR_W) $(PDFBUILD) && cd ${CUBADIR} && make clean

distclean:
	rm -rf ${MAINBUILD} ${BUILDDIR} $(BUILDDIR_Z) $(BUILDDIR_W) $(PDFBUILD) ${EXEC_Z} ${CONDOR_EXEC_Z} ${EXEC_W} ${CONDOR_EXEC_W} && cd ${CUBADIR} && make clean && rm libcuba.a

