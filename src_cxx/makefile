
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

GSL_DIR=/home/ld7/bin/gsl-1.16/installation
GSL_INCLUDE_DIR=-I$(GSL_DIR)/include
GSL_LINK_DIR=-L$(GSL_DIR)/lib -lgsl -lgslcblas -limf -lm 

COMPILE = $(GSL_INCLUDE_DIR) $(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
LINK = $(GSL_LINK_DIR) $(LINK_FLAGS)
%.o : %.cpp
	 -${CLINKER} $(COMPILE) -c $<

EXE = mbl
all: ${OBJ}  chkopts
	-${CLINKER} -o ${EXE} ${OBJ}  $(LINK) ${PETSC_KSP_LIB} 
#	${RM} ${OBJ}

touch:
	touch *.cpp *.h




