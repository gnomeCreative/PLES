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
FFLAGS	=  -O3 -fdefault-real-8 -fcheck=all -fbounds-check -fbacktrace -finit-real=zero -J $(FMOD) -I $(FMOD) -Wall
#FFLAGS	=  -fdefault-real-8 -fcheck=all -fbounds-check -fbacktrace -finit-real=zero -J $(FMOD) -I $(FMOD)

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
LIBS	= -lmpi_mpifh -lgfortran -lm -lmpi -lmpi_usempi #-L/cineca/lib/mass -lmass -lmpi_mpifh

# Flag for linking
LDFLAGS	= $(FFLAGS) $(CFLAGS) 


# List of Fortran code files 
FOBJS   =  \
 scala3.o subgrid.o period.o convex.o print.o tipologia.o parti.o parete.o orl.o velpar.o mysettings.o mysending.o \
 \
 mysettings_boundary.o myarrays_metri3.o myarrays_density.o myarrays_cor.o myarrays_WB.o myarrays_LC.o myarrays_moisture.o \
 myarrays_ibm.o myarrays_velo3.o myarrays_wallmodel.o myarrays2.o myarrays_nesting.o myarrays_buffer_bodyforce.o\
 turbo2_data.o turbo3bis.o multigrid.o ricercaGeomRoutines.o rotation.o trilinear.o ricerca.o \
 \
 contour_module.o flucn_module.o jord_module.o output_module.o ada_rho.o adams.o angolo.o aree_parziali.o \
 filter_module.o buffer_bodyforce.o carico_immb.o cellep.o check_divergence.o coed1.o coed2.o coed3_tr.o\
 coef1_par.o coef2_par.o condi1.o  condi2.o condi3.o \
 contra.o contra_infout.o contrin.o contrin_lat.o contrin_pot.o correggi_ib.o courant.o \
 diver.o drift.o disturbance_face1.o disturbance_face2.o disturbance_face5.o disturbance_face6.o \
 facce.o filtro.o flu_turbo.o  flucnesp.o flucrho.o flucrhoesp.o flud1.o flud2.o flud3.o \
 flux1.o flux2.o flux3.o fmassa.o gauss_random.o genero_random.o gprima.o gradie.o gseconda.o gterza.o\
 indy.o inflow.o iniz.o iniz_metrica.o inverse_para2.o interp3.o interpola_pareti_pom.o \
 leggi_pareti_pom.o langmuir2.o leggivento.o metrica.o mix_para.o mixrho_para.o \
 nesting.o orlansky_generale.o passo_ibm.o passo_solide.o \
 periodic.o prepare_nesting.o prolong.o \
 readflux.o readtau.o redistribuzione.o resid_sndrcv.o resid_sndrcv_cart.o restart.o \
 restrict_sndrcv.o rhostress.o rhs1_rho.o sett.o set_transpose_implicit.o tridag.o \
 tridiag_trasp_para.o tridiag_trasp_para_rho.o triper.o turbo_lagrangian.o update.o \
 vel_up.o vorticitag.o wall_function_bodyfitted.o turbo_dinamico.o turbo_statico.o wernerwengle.o\
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
	
# Code in c++
$(CS)/%.o: $(CCODE)/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@



