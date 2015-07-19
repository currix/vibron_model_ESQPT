#
# Makefile for endohedral H2 project
#
# by Currix TM
#
BINPATH = ../bin/
#
U4_vibron_HAM_SO4_SRC = nrtype.f90  defparam_vibron.f90 \
vibron_model_eigensystem_SO4.f90
#
U4_vibron_HAM_SO4_OBJ = $(U4_vibron_HAM_SO4_SRC:.f90=.o)
#
U4_vibron_HAM_SO4_IPR_SRC = nrtype.f90  defparam_vibron.f90 \
vibron_model_eigensystem_SO4_ipr.f90
#
U4_vibron_HAM_SO4_IPR_OBJ = $(U4_vibron_HAM_SO4_IPR_SRC:.f90=.o)
#
U4_vibron_HAM_U3_SRC = nrtype.f90  defparam_vibron.f90	\
vibron_model_eigensystem_U3.f90
#
U4_vibron_HAM_U3_OBJ = $(U4_vibron_HAM_U3_SRC:.f90=.o)
#
U4_vibron_HAM_U3_IPR_SRC = nrtype.f90  defparam_vibron.f90	\
vibron_model_eigensystem_U3_ipr.f90
#
U4_vibron_HAM_U3_IPR_OBJ = $(U4_vibron_HAM_U3_IPR_SRC:.f90=.o)
#####################################################################################################
##
FC = gfortran
FOPT	= -c 
FOPT   += -Wall
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
FOPT   += $(MODLIB)
#FOPT   += -O0
#
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
## Cluster
LAPACK77 = -llapack
BLAS     = -lblas
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
######################################################################################################
all: build_Vibron_HAM_SO4 build_Vibron_HAM_U3  build_Vibron_HAM_SO4_ipr build_Vibron_HAM_U3_ipr
.PHONY : all
#
# Standard compilation rules
#
%.o: %.f90
	$(FC) $(FOPT) -c $<
#
%.o: %.F90
	$(FC) $(FOPT) -c $<
#
%.o: %.f
	$(FC) $(FOPT) $<
#
#######################
build_Vibron_HAM_SO4: $(U4_vibron_HAM_SO4_OBJ) 
	$(FC) -o $(BINPATH)/vibron_u4_so4_energies  $(U4_vibron_HAM_SO4_OBJ) $(LIBS) 
#######################
build_Vibron_HAM_U3: $(U4_vibron_HAM_U3_OBJ) 
	$(FC) -o $(BINPATH)/vibron_u4_u3_energies  $(U4_vibron_HAM_U3_OBJ) $(LIBS) 
#######################
build_Vibron_HAM_SO4_ipr: $(U4_vibron_HAM_SO4_IPR_OBJ) 
	$(FC) -o $(BINPATH)/vibron_u4_so4_ipr  $(U4_vibron_HAM_SO4_IPR_OBJ) $(LIBS) 
#######################
build_Vibron_HAM_U3_ipr: $(U4_vibron_HAM_U3_IPR_OBJ) 
	$(FC) -o $(BINPATH)/vibron_u4_u3_ipr  $(U4_vibron_HAM_U3_IPR_OBJ) $(LIBS) 
#######################
clean:
	@rm -f  *.o *.mod




