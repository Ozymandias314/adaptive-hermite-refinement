# Makefile for Nuno
#
# type "make" to build an executable file 
#
########################################################################

#############################
# WITH_HAC is used always.
#############################
#WITH_HAC=1 # comment this line for a HAC-free build on helios and nino

#############################
# Defaults
# (These get used when 
# VIRIATO_SYSTEM is not set)
############################

ifeq "$(findstring nino,$(shell hostname))" "nino"
ifdef WITH_HAC
#ADIOS_DIR=/usr/local
ADIOS_DIR=/opt/adios-1.8.0
ADIOS_INC = $(shell $(ADIOS_DIR)/bin/adios_config -c -f)
ADIOS_FLIB = $(shell $(ADIOS_DIR)/bin/adios_config -l -f)
HAC_OBJ = hlst_adios_checkpoint.o
HAC_INC=-I.
HAC_FLAGS=-D''HLST_HAC_USE_REAL=1 -D''HLST_HAC_USE_CMSK=0 -D''HLST_HAC_USE_ADIM=4 $(HAC_INC)
VIRIATIO_FLAGS= -D WITH_HAC=1 $(ADIOS_INC) $(HAC_FLAGS)
VIRIATIO_LIBS = $(ADIOS_FLIB)
endif
VIRIATIO_FLAGS=
VIRIATIO_LIBS =
FC_ =  mpif90
FLIBS_ = -llapack /opt/fftw-2.1.5/lib/*.a $(VIRIATIO_LIBS)
#F90FLAGS_ = -O0 -ggdb -pipe -fdefault-real-8 -fdefault-double-8  -Dgasca3d $(VIRIATIO_FLAGS) -Wall -ffree-line-length-none
F90FLAGS_ = -O0 -ggdb -pg -pipe -fdefault-real-8 -fdefault-double-8  -Dgasca3d $(VIRIATIO_FLAGS) -Wall
#F90FLAGS_ = -O3 -pg -pipe -fdefault-real-8 -fdefault-double-8  -Dgasca3d $(VIRIATIO_FLAGS) -Wall
# -fdefault-integer-8 breaks MPI calls
else
FC_ =  mpif90
FLIBS_ = -I/opt/apps/lahey/L8.10b/fftw/2.1.5/include -L/opt/apps/lahey/L8.10b/fftw/2.1.5/lib -ldfftw -ldrfftw
F90FLAGS_ = -O3  --dbl --ml cdecl
endif


###############
# Bass
##############
FC_bass =  mpif90
FLIBS_bass = -I/opt/apps/lahey/L8.10b/fftw/2.1.5/include -L/opt/apps/lahey/L8.10b/fftw/2.1.5/lib -ldfftw -ldrfftw
F90FLAGS_bass = -O3  --dbl --ml cdecl

###############
#Linux
##############

FC_linux =  mpif90
FLIBS_linux =   -lfftw -lrfftw
F90FLAGS_linux = -O3  -fdefault-real-8

###############
# helios
##############
#
#HAC_LIB = 
HAC_INC=
FC_ =  mpif90
#FC_ =  $(MPIFC)
ifdef WITH_HAC
ADIOS_INC = $(shell $(ADIOS_DIR)/bin/adios_config -c -f)
ADIOS_FLIB = $(shell $(ADIOS_DIR)/bin/adios_config -l -f)
HAC_FLAGS=-D''HLST_HAC_USE_REAL=1 -D''HLST_HAC_USE_CMSK=0 -D''HLST_HAC_USE_ADIM=4 $(HAC_INC)
VIRIATIO_FLAGS= -D WITH_HAC $(ADIOS_INC) $(HAC_FLAGS)
VIRIATIO_LIBS = $(ADIOS_FLIB)
HAC_OBJ = hlst_adios_checkpoint.o
else
VIRIATIO_FLAGS=
VIRIATIO_LIBS =
HAC_OBJ =
endif

#FC_helios = mpif90
FC_helios = $(MPIFC)
#F90FLAGS_helios = -O2 -r8 -i8 -g -C
#F90FLAGS_helios = -O3 -r8 -i8 -xAVX -ipo  
#ttr F90FLAGS_helios = -O3 -r8 -i8 -xhost  
F90FLAGS_helios = -O3 -r8 -i8 $(VIRIATIO_FLAGS)
###F90FLAGS_helios = -O3 -r8
FLIBS_helios = -lfftw -lrfftw -mkl $(VIRIATIO_LIBS)
#
ifdef GASCA2D
FPPFLAGS_helios = -Dgasca2d
endif
#
ifdef GASCA3D
FPPFLAGS_helios = -Dgasca3d
endif
#
ifdef GASCA4D
FPPFLAGS_helios = -Dgasca4d -Dgasca3d
endif
#
FPPFLAGS_helios += -Dmpi1
#
ifdef PERFLIB
#need to include the following lines in the .bash_profile for PERF to work
#    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/csc/home0/dannert/helios_inst/lib (OLD!)
#    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/csc/home3/ttr/perf/lib
#    export LD_LIBRARY_PATH
FPPFLAGS_helios += -Dperf
PERFLIB_HOME_helios = /csc/home3/ttr/perf
FLIBS_helios  += -L$(PERFLIB_HOME_helios)/lib -looperf -lpfm -lstdc++
#
ifdef FORTDIAGS
F90FLAGS_helios:=$(F90FLAGS_helios) -check all -warn all,nodec,interfaces -gen_interfaces -traceback
endif
#
endif

###############
# hydra (rzg)
##############

FC_hydra = mpiifort
F90FLAGS_hydra = -O3 -r8 -i8
FLIBS_hydra = -L$(FFTW_HOME)/lib -ldfftw -ldrfftw \
$(MKL_HOME)/lib/intel64/libmkl_intel_lp64.a  -Wl,--start-group \
        $(MKL_HOME)/lib/intel64/libmkl_sequential.a \
        $(MKL_HOME)/lib/intel64/libmkl_core.a \
        -Wl,--end-group -lpthread
#
ifdef GASCA2D
FPPFLAGS_hydra = -Dgasca2d
endif
#
ifdef GASCA3D
FPPFLAGS_hydra = -Dgasca3d
endif
#
ifdef GASCA4D
FPPFLAGS_hydra = -Dgasca4d -Dgasca3d
endif
#
FPPFLAGS_hydra += -Dmpi1
#
ifdef PERFLIB
#need to issue this at command line for PERF to work
#   module load perflib
#   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$PERFLIB_HOME/lib
FPPFLAGS_hydra += -Dperf
FLIBS_hydra  += -L$(PERFLIB_HOME)/lib -looperf
#
ifdef FORTDIAGS
F90FLAGS_hydra:=$(F90FLAGS_hydra) -check all -warn all,nodec,interfaces -gen_interfaces -traceback
endif
#
endif

###############
# hpcff
##############

FC_hpcff = mpif90
F90FLAGS_hpcff = -O3 -ipo -axSSE4.2 -r8 -i8 #  -default64  
FLIBS_hpcff = -L /usr/local/fftw/2.1.5/lib -ldfftw -ldrfftw

###############                                                                                                                                                                                        
# engaging                                                                                                                                                                                             
#############                                                                                 

H5_FCOMPILEFLAGS = /home/software/psfc/pkg/hdf5-1.10.0-patch1/intel-17/include
H5_FLINKFLAGS = /home/software/psfc/pkg/hdf5-1.10.0-patch1/intel-17/lib -lhdf5_fortran -lhdf5 -lz -lm

FC_engaging = mpiifort
F90FLAGS_engaging = -O3 -i8 -r8 -mkl=cluster -I$(FFTWINCLUDE) -I$(H5_FCOMPILEFLAGS)
FLIBS_engaging = -L$(FFTWLIB) -ldfftw -ldrfftw -L$(H5_FLINKFLAGS)  # -lreflapack -lscalapack -lrefblas -ltmg


ifdef GASCA2D
FPPFLAGS_engaging = -Dgasca2d
endif
#
ifdef GASCA3D
FPPFLAGS_engaging = -Dgasca3d
endif
#
ifdef GASCA4D
FPPFLAGS_engaging = -Dgasca4d -Dgasca3d
endif
#
###############
# Hopper
##############

FC_hopper = ftn
F90FLAGS_hopper =-target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I/opt/fftw/2.1.5.3/include
FLIBS_hopper = -L/opt/fftw/2.1.5.3/lib -ldfftw -ldrfftw

###############
# Edison
##############

FC_edison = ftn
F90FLAGS_edison =-target=linux -r8 -fast -no-ipo -I. -Iutils -I/opt/fftw/2.1.5.4/include
FLIBS_edison = -L/opt/fftw/2.1.5.4/lib -ldfftw -ldrfftw

###############
# Kraken
##############
FC_kraken =  ftn
F90FLAGS_kraken = -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include
FLIBS_kraken = -L/opt/fftw/2.1.5.3/lib -ldfftw -ldrfftw

###############
# Stampede
##############
HAC_INC=
#FC_ =  mpif90
ifdef WITH_HAC
ADIOS_INC = $(shell $(ADIOS_DIR)/bin/adios_config -c -f)
ADIOS_FLIB = $(shell $(ADIOS_DIR)/bin/adios_config -l -f)
HAC_FLAGS=-D''HLST_HAC_USE_REAL=1 -D''HLST_HAC_USE_CMSK=0 -D''HLST_HAC_USE_ADIM=4 $(HAC_INC)
VIRIATIO_FLAGS= -D WITH_HAC $(ADIOS_INC) $(HAC_FLAGS)
VIRIATIO_LIBS = $(ADIOS_FLIB)
HAC_OBJ = hlst_adios_checkpoint.o
else
VIRIATIO_FLAGS=
VIRIATIO_LIBS =
HAC_OBJ =
endif

FC_stampede =  mpif90
FLIBS_stampede = -L/opt/apps/intel13/impi_4_1/fftw2/2.1.5/lib -ldrfftw -ldfftw -mkl $(VIRIATIO_LIBS)
F90FLAGS_stampede = -O3 -r8 -i8 -xhost -I/opt/apps/intel13/impi_4_1/fftw2/2.1.5/include  $(VIRIATIO_FLAGS)
#
ifdef GASCA2D
FPPFLAGS_stampede = -Dgasca2d
endif
#
ifdef GASCA3D
FPPFLAGS_stampede = -Dgasca3d
endif
#
ifdef GASCA4D
FPPFLAGS_stampede = -Dgasca4d -Dgasca3d
endif
#
FPPFLAGS_stampede += -Dmpi1
#
ifdef PERFLIB
#need to include the following lines in the .bash_profile for PERF to work
#    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/csc/home0/dannert/helios_inst/lib (OLD!)
#    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/csc/home3/ttr/perf/lib
#    export LD_LIBRARY_PATH
#FPPFLAGS_stampede += -Dperf
PERFLIB_HOME_stampede = /scratch/00295/tg455668/perf/
FLIBS_stampede  += -L$(PERFLIB_HOME_stampede)/lib -looperf -lpfm -lstdc++
#
ifdef FORTDIAGS
F90FLAGS_stampede:=$(F90FLAGS_stampede) -check all -warn all,nodec,interfaces -gen_interfaces -traceback
endif
#
endif

###############
# Mac amd64_darwin_120
##############

FC_mac =  mpif90-openmpi-gcc5
#FLIBS_mac =   -L/opt/local/lib -lfftw3 -lrfftw3
FLIBS_mac =   -L/opt/local/lib -lfftw3
F90FLAGS_mac = -O3  -fdefault-real-8
#
ifdef GASCA2D
FPPFLAGS_mac = -Dgasca2d
endif
#
ifdef GASCA3D
FPPFLAGS_mac = -Dgasca3d
endif
#
ifdef GASCA4D
###FPPFLAGS_mac += -Dgasca4d
FPPFLAGS_mac = -Dgasca4d -Dgasca3d
endif
#
FPPFLAGS_mac += -Dmpi1


##############
# Set actual variables
###########

FC = ${FC_${VIRIATO_SYSTEM}}
FLIBS = ${FLIBS_${VIRIATO_SYSTEM}}
F90FLAGS = ${F90FLAGS_${VIRIATO_SYSTEM}}
FPPFLAGS = ${FPPFLAGS_${VIRIATO_SYSTEM}}

#TTR
ifdef GPROF
F90FLAGS:=$(F90FLAGS) -g -p
endif

#TTR
ifdef MAP
F90FLAGS:=$(F90FLAGS) -g
endif

#FC =  h5pfc
#FLIBS = -L $$HOME/FFTW/lib -lfftw -lrfftw

#F90FLAGS = -r8 -i8 #  -default64  
#F90FLAGS = -O2 -r8 -i8 -C -g
#F90FLAGS = -g -C -ffortran-bounds-check -default64 -zerouv

#F90FLAGS = -i8 -r8 -O2
#F90FLAGS = -g -qfree=f90 -qsuffix=f=f90  -qrealsize=8 -qinitauto -C -qmoddir=/tmp/nfl -I/tmp/nfl 

#ifeq ($(debug),on) 
#   F90FLAGS += -g -qsuffix=f=f90 -qstrict
#else
#   F90FLAGS += -O3 -qsuffix=f=f90 -qstrict
#endif

OBJS = REGK.o fft_work_fftw.o constants.o transforms.o grid.o \
        mpi_mod.o mp_mpi_r8.o redistribute_mpi.o diag.o forcing.o \
        functions.o initialize.o brackets.o fluxes.o aux.o stepping.o $(HAC_OBJ)

OBJS2 = test_par.o fft_work_fftw.o constants.o transforms.o grid.o \
        mp_mpi_r8.o mpi_mod.o redistribute_mpi.o
OBJS3 =contours.o fft_work_fftw.o constants.o transforms.o grid.o \
        mp_mpi_r8.o mpi_mod.o redistribute_mpi.o diag.o

OBJS_PP = postproc.o fft_work_fftw.o constants.o transforms.o grid.o \
        mp_mpi_r8.o mpi_mod.o redistribute_mpi.o diag.o

FFLAGS = $(F90FLAGS)

.SUFFIXES: .F90 .f90 .f

.F90.o:
	$(FC) $(F90FLAGS) $(FPPFLAGS) -c $<

.f90.o: 
	$(FC) $(F90FLAGS) -c $<

.f.o: 
	$(FC) $(FFLAGS) -c $<

%.o : %.mod

all:    REGK h2v

$(HAC_OBJ): $(patsubst %o, %F90 ,$(HAC_OBJ))
	$(FC) $(ADIOS_INC) $(HAC_FLAGS) -c $<

REGK:   $(OBJS)
	$(FC) $(FFLAGS) -o viriato $(OBJS) $(FLIBS)

test2:  test2.o fft_work_fftw.o
	$(FC) $(FFLAGS) -o test2 test2.o fft_work_fftw.o $(FLIBS)

test3:  test3.o fft_work_fftw.o constants.mod
	$(FC) $(FFLAGS) -o test3 test3.o fft_work_fftw.o constants.mod  $(FLIBS)

test:   test.o fft_work_fftw.o
	$(FC) $(FFLAGS) -o test test.o fft_work_fftw.o $(FLIBS)

testpar: $(OBJS2)
	$(FC) $(FFLAGS) -o testpar $(OBJS2) $(FLIBS)
clean:
	rm -f *.o *.mod *~ regk.e* regk.o* exit debugregi tomates contours.e* contours.o*

contours: $(OBJS3)
	$(FC) $(FFLAGS) -o contours $(OBJS3) $(FLIBS) 

postproc: $(OBJS_PP)
	$(FC) $(FFLAGS) -o postproc $(OBJS_PP) $(FLIBS)

h2v: h2v.o $(filter-out REGK.o,$(OBJS))
	$(FC) $(FFLAGS) -o $@ $+ $(FLIBS) 
test_make:
	@echo $(USE_FFT)
	@echo $(HOME)
# dependencies:

REGK.o: constants.o mp_mpi_r8.o transforms.o grid.o diag.o forcing.o \
        functions.o initialize.o brackets.o stepping.o fluxes.o aux.o $(HAC_OBJ)
test3.o: constants.o fft_work_fftw.o
grid.o: constants.mod mp_mpi_r8.o
transforms.o: constants.mod fft_work_fftw.o grid.o redistribute_mpi.o
diag.o: mp_mpi_r8.o constants.mod grid.o transforms.o $(HAC_OBJ)
constants.o: mp_mpi_r8.o
test_par.o: constants.mod mp_mpi_r8.o transforms.o grid.o
contours.o: constants.mod mp_mpi_r8.o transforms.o grid.o diag.o
forcing.o: constants.mod grid.o
postproc.o: constants.mod mp_mpi_r8.o transforms.o grid.o diag.o
mp_mpi_r8.o: mpi_mod.o
functions.o: constants.mod grid.o
initialize.o: constants.mod mp_mpi_r8.o transforms.o grid.o forcing.o
brackets.o: constants.mod mp_mpi_r8.o grid.o functions.o transforms.o
fluxes.o: constants.mod mp_mpi_r8.o
aux.o: constants.mod mp_mpi_r8.o grid.o transforms.o
stepping.o: aux.o constants.mod mp_mpi_r8.o grid.o functions.o fluxes.o


