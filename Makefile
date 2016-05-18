# /usr/sbin/smake
#
FC	= mpifort 
CC	= mpic++
LD	= $(CC)
F77	= $(FC)
RM	= /bin/rm -f
#
FCODE = fsrc
CCODE = csrc
FS = build/fsrc
CS = build/csrc
FMOD = $(FS)/modules
#

# Flag for FORTRAN compiler
FIXED   =  -extend_source 
FFLAGS	=  -O3 -fdefault-real-8 -fcheck=all -fbounds-check -fbacktrace -finit-real=zero -J $(FMOD) -I $(FMOD)
#FFLAGS	=  -fdefault-real-8 -fcheck=all -fbounds-check -fbacktrace -finit-real=zero -J $(FMOD) -I $(FMOD) -Wall

# -Wunused
#FFLAGS  = $(FIXED) -align -O2  -ftz 

#FFLAGS = -O0 -g -traceback -check bounds  -check bounds -fpe0 -i_dynamic -r8

# Flag for PROFILING: -g -pg
# after compilation it's necessary to run the code and use the txt file as:
# gprof <executable> gmon.out

# Flag for large files (ENZO): -bmaxdata:0x80000000
# Flag for large files (webSP4): -bmaxdata:0x70000000


# Flag for C++ compiler #-Wall
CFLAGS	=  -O0 -fmessage-length=0
#CFLAGS	=  -O0 -fmessage-length=0

# Used libraries
LIBS	= -lmpi_mpifh -lgfortran -lm -lmpi -lquadmath #-lmpi_usempi #-L/cineca/lib/mass -lmass -lmpi_mpifh

# Flag for linking
LDFLAGS	= $(FFLAGS) $(CFLAGS) 


# List of Fortran code files 
FOBJS   =  \
 scala3.o mysending.o period.o mysettings.o subgrid.o velpar.o geometricRoutines.o\
 \
 myarrays_velo3.o myarrays_metri3.o wind_module.o \
 orlansky_module.o turbo_module.o wallmodel_module.o trilinear.o\
 ibm_module.o tridag_module.o particle_module.o  multigrid_module.o  \
 contour_module.o ricerca_module.o output_module.o buffer_bodyforce_module.o inflow_module.o jord_module.o \
 filtro_module.o filter_module.o flu_module.o nesting_module.o \
 \
 ada_rho.o adams.o cellep.o check_divergence.o voronoi.o\
 contra.o contrin.o correggi_ib.o courant.o diver.o \
 fmassa.o gradie.o indy.o inverse_para2.o \
 metrica.o mix_para.o mixrho_para.o \
 restart.o rhs1_rho.o turbo_lagrdin.o \
 vel_up.o turbo_statico.o \
 \
 stratiParticle.o

# List of C++ code files 

COBJS	= myvector.o utils.o elmt.o IO.o DEM.o wrapper.o

OBJS	= $(FOBJS:%.o=$(FS)/%.o) $(COBJS:%.o=$(CS)/%.o)

EXEC	= ./strati_p.exe

# all Target 
all: $(EXEC)

$(EXEC):	$(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(OBJS) $(LIBS)

# clean Target
clean:
	$(RM) $(EXEC) $(FS)/*.o $(CS)/*.o $(FMOD)/*.mod fort.* 

.SUFFIXES: 
	.o .f03 .cpp .h

# Code in Fortran
$(FS)/%.o: $(FCODE)/%.f03
	$(FC) -c $(FFLAGS) $< -o $@
	
$(FS)/%.mod: $(FCODE)/%.f03 $(FS)/%.o
	@true
	
$(FS)/%.o: $(FCODE)/%.f03 $(FS)/%.mod
	$(FC) -c $(FFLAGS) $< -o $@
	
# Code in c++
$(CS)/%.o: $(CCODE)/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(CS)/%.o: $(CCODE)/%.h
	$(CC) -c $(CFLAGS) $< -o $@
	


