#
# Makefile for ESQPT Dynamics Calculations in U(4)
#
# by Currix TM
#
BINPATH = ../bin/
#
U4_vibron_HAM_SO4_SRC = nrtype.f90  defparam_vibron.f90 \
vibron_model_eigensystem_SO4.f90
#
U4_vibron_HAM_SO4_IPR_SRC = nrtype.f90  defparam_vibron.f90 \
vibron_model_eigensystem_SO4_ipr.f90
#
U4_vibron_HAM_U3_SRC = nrtype.f90  defparam_vibron.f90	\
vibron_model_eigensystem_U3.f90
#
U4_vibron_HAM_U3_IPR_SRC = nrtype.f90  defparam_vibron.f90	\
vibron_model_eigensystem_U3_ipr.f90
#####################################################################################################
##
FC = gfortran
FOPT   = -Wall
#FOPT   += -O0
#
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
#
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
#
LAPACK77 = -llapack
BLAS     = -lblas
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
######################################################################################################
all: vibron_u4_so4_energies vibron_u4_u3_energies vibron_u4_so4_ipr vibron_u4_u3_ipr
#
.PHONY : all
#
#######################
vibron_u4_so4_energies: $(U4_vibron_HAM_SO4_SRC) Makefile 
	$(FC) -o $(BINPATH)/$@  $(MODLIB) $(U4_vibron_HAM_SO4_SRC) $(LIBS) 
#######################
vibron_u4_u3_energies: $(U4_vibron_HAM_U3_SRC)  Makefile 
	$(FC) -o $(BINPATH)/$@  $(MODLIB) $(U4_vibron_HAM_U3_SRC) $(LIBS) 
#######################
vibron_u4_so4_ipr: $(U4_vibron_HAM_SO4_IPR_SRC) Makefile 
	$(FC) -o $(BINPATH)/$@ $(MODLIB) $(U4_vibron_HAM_SO4_IPR_SRC) $(LIBS) 
#######################
vibron_u4_u3_ipr: $(U4_vibron_HAM_U3_IPR_SRC) Makefile 
	$(FC) -o $(BINPATH)/$@ $(MODLIB) $(U4_vibron_HAM_U3_IPR_SRC) $(LIBS) 
#######################
clean:
	@rm -f  *.o *.mod vibron_u4_u3_ipr vibron_u4_so4_ipr vibron_u4_u3_energies vibron_u4_so4_energies




