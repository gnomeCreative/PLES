# /usr/sbin/smake
#
FC	= mpif90 
CC	= mpic++
LD	= $(CC)
F77	= $(FC)
RM	= /bin/rm -f
#
FCODE=../fsrc
FS = ./fsrc
CS = ./csrc
FMOD = $(FS)/modules
#

# Flag for FORTRAN compiler
FIXED   =  -extend_source -fixed
FFLAGS	=  -fdefault-real-8 -fcheck=all -fbacktrace -finit-real=zero -J $(FMOD) -I $(FMOD)
# -Wunused
#FFLAGS  = $(FIXED) -align -O2  -ftz 
#FFLAGS = -O0 -g -traceback -check bounds  -check bounds -fpe0 -i_dynamic -r8

# Flag for PROFILING: -g -pg
# after compilation it's necessary to run the code and use the txt file as:
# gprof <executable> gmon.out

# Flag for large files (ENZO): -bmaxdata:0x80000000
# Flag for large files (webSP4): -bmaxdata:0x70000000


# Flag for C++ compiler
CFLAGS	=  -O0 -Wall -fmessage-length=0

# Used libraries
LIBS	= -lmpi_mpifh -lgfortran -lm -lmpi -lmpi_usempi #-L/cineca/lib/mass -lmass

# Flag for linking
LDFLAGS	= $(FFLAGS) $(CFLAGS) 


# List of Fortran code files 
FOBJS   =  \
 $(FS)/scala3.o $(FS)/subgrid.o $(FS)/period.o $(FS)/convex.o \
 $(FS)/print.o $(FS)/tipologia.o $(FS)/parti.o $(FS)/parete.o $(FS)/orl.o $(FS)/velpar.o\
 $(FS)/mysettings.o $(FS)/mysending.o \
 \
 $(FS)/mysettings_boundary.o\
 $(FS)/myarrays_metri3.o $(FS)/myarrays_density.o\
 $(FS)/myarrays_cor.o $(FS)/myarrays_WB.o $(FS)/myarrays_LC.o $(FS)/myarrays_moisture.o \
 $(FS)/myarrays_ibm.o $(FS)/myarrays_velo3.o $(FS)/myarrays_wallmodel.o\
 $(FS)/myarrays2.o $(FS)/myarrays_nesting.o $(FS)/myarrays_buffer_bodyforce.o\
 $(FS)/turbo2_data.o $(FS)/turbo3bis.o \
 $(FS)/multigrid.o \
 \
 $(FS)/contourp_se_module.o \
 $(FS)/ts_module.o $(FS)/flucn_module.o $(FS)/jord_module.o $(FS)/output_module.o \
 $(FS)/ada_rho.o $(FS)/adams.o $(FS)/angolo.o $(FS)/aree_parziali.o \
 $(FS)/buffer_bodyforce.o $(FS)/buffer1g.o $(FS)/buffer1d_par.o $(FS)/buffer2d_par.o \
 $(FS)/buffer2.o $(FS)/buffer2g.o $(FS)/buffer2gg.o \
 $(FS)/buffer1old_par.o $(FS)/buffer2old_par.o $(FS)/buffer1old_par_nscal.o $(FS)/buffer2old_par_nscal.o \
 $(FS)/buffvect1d.o $(FS)/buffvect2d.o \
 $(FS)/carico_immb.o $(FS)/cellep.o $(FS)/celle_vicine_immb.o $(FS)/check_divergence.o \
 $(FS)/coed1.o $(FS)/coed2.o $(FS)/coed3_tr.o $(FS)/compute_ustar.o \
 $(FS)/coef1_par.o $(FS)/coef2_par.o $(FS)/condi1.o  $(FS)/condi2.o $(FS)/condi3.o  \
 $(FS)/contour.o $(FS)/contour_se.o $(FS)/contour_se_nesting.o $(FS)/contourp.o  \
 $(FS)/contra.o $(FS)/contra_infout.o $(FS)/contrin.o $(FS)/contrin_lat.o $(FS)/contrin_pot.o $(FS)/correggi_ib.o $(FS)/courant.o \
 $(FS)/diver.o $(FS)/drift.o \
 $(FS)/disturbance_face1.o $(FS)/disturbance_face2.o $(FS)/disturbance_face5.o $(FS)/disturbance_face6.o \
 $(FS)/facce.o $(FS)/filterb_csieta.o $(FS)/filterb_zita.o  $(FS)/filtro.o \
 $(FS)/filter01.o $(FS)/filter02.o $(FS)/filter03.o $(FS)/filter01np.o \
 $(FS)/filter02np.o $(FS)/filter04g.o $(FS)/filter05g.o \
 $(FS)/filter04b.o $(FS)/filter06b.o $(FS)/filter06.o \
 $(FS)/flu_turbo.o  $(FS)/flucnesp.o $(FS)/flucrho.o $(FS)/flucrhoesp.o \
 $(FS)/flud1.o $(FS)/flud2.o $(FS)/flud3.o \
 $(FS)/flux1.o $(FS)/flux2.o $(FS)/flux3.o $(FS)/fmassa.o \
 $(FS)/gauss_random.o $(FS)/genero_random.o $(FS)/gprima.o $(FS)/gradie.o $(FS)/gseconda.o $(FS)/gterza.o\
 $(FS)/incrementi.o $(FS)/indy.o $(FS)/inflow.o $(FS)/iniz.o $(FS)/iniz_metrica.o \
 $(FS)/inverse_para2.o $(FS)/interp3.o $(FS)/interpola_pareti_pom.o \
 $(FS)/leggi_pareti_pom.o $(FS)/langmuir2.o $(FS)/leggivento.o \
 $(FS)/metrica.o $(FS)/mix_para.o $(FS)/mixrho_para.o \
 $(FS)/nesting.o $(FS)/orlansky_generale.o $(FS)/output_ibm.o \
 $(FS)/passo_ibm.o $(FS)/passo_solide.o \
 $(FS)/periodic.o $(FS)/potenziale_ibm.o\
 $(FS)/prepare_nesting.o $(FS)/prolong.o \
 $(FS)/readflux.o $(FS)/readtau.o \
 $(FS)/read_simulation_setting.o $(FS)/redistribuzione.o $(FS)/resid_sndrcv.o $(FS)/resid_sndrcv_cart.o $(FS)/restart.o \
 $(FS)/restrict_sndrcv.o $(FS)/rhostress.o $(FS)/rhs1_rho.o \
 $(FS)/sett.o $(FS)/set_transpose_implicit.o $(FS)/tridag.o \
 $(FS)/tridiag_trasp_para.o $(FS)/tridiag_trasp_para_rho.o \
 $(FS)/triper.o $(FS)/turbo_lagrangian.o $(FS)/update.o \
 $(FS)/vel_up.o $(FS)/vorticitag.o \
 $(FS)/wall_function_bodyfitted.o $(FS)/turbo_dinamico.o \
 $(FS)/wernerwengle.o\
 \
 $(FS)/stratiParticle.o

# List of C++ code files 
COBJS	= $(CS)/wrapper.o

OBJS	= $(FOBJS) $(COBJS) 

EXEC	= ../strati_p.exe

# all Target 
all: $(EXEC)

$(EXEC):	$(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(OBJS) $(LIBS)

# clean Target
clean:
	$(RM) $(EXEC) $(OBJS) *.lst *.mod ./$(FMOD)/*.mod

.SUFFIXES: 
	.o .F .c .f .f03 .cpp

# Code in Fortran
$(FS)/%.o: ../fsrc/%.f03
	$(FC) -c $(FFLAGS) -o "$@" "$<"

# Code in c++
$(CS)/%.o: ../csrc/%.cpp
	$(CC) -c $(CFLAGS) -o "$@" "$<"


