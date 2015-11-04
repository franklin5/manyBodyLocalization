
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = 
EXAMPLESC        = main.cpp polaron.cpp 
EXAMPLESF        = 
MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
OBJ=$(EXAMPLESC:.cpp=.o)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

GSL_DIR=/work/02946/ld7/bin/gsl-2.0/installation
GSL_INCLUDE_DIR=-I$(GSL_DIR)/include 
GSL_LINK_DIR=-L$(GSL_DIR)/lib -lgsl -lgslcblas -limf -lm 

COMPILE = $(GSL_INCLUDE_DIR) $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
LINK = $(GSL_LINK_DIR) $(LINK_FLAGS)
%.o : %.cpp
	 -${CLINKER} $(COMPILE) -c $<

EXE = mbl
all: ${OBJ}  chkopts
	-${CLINKER} -o ${EXE} ${OBJ}  $(LINK) ${PETSC_KSP_LIB} ${SLEPC_EPS_LIB}
#	${RM} ${OBJ}

touch:
	touch *.cpp *.h




