
static char help[] = "Localiazation and delocalization transition of a polaron on a 1D disordered lattice.\n\n";

/*T
   Concepts:
   Processors: n
   mpirun -n 6 mbl
   or
   ./mbl
T*/


#include "polaron.h"
#undef __FUNCT__
#define __FUNCT__ "main"
#define root 0
int main(int argc,char **args){
  PetscErrorCode ierr;
  SlepcInitialize(&argc,&args,(char*)0,help);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     "======================================================================\n"
		     "The purpose of this program is to study the Localiazation and \n"
		     "  delocalization transition of a polaron on a 1D disordered lattice.\n"
		     " Initiated and researched the matlab code by Shangshun Zhang,  \n"
		     " proposed and implemented the parallel C++ code by Lin Dong  \n"
		     " at Rice University. Run at "  __TIME__  ", on "  __DATE__  "\n"
		     "Petsc is initialized and program starts from " __FILE__  "\n"
		     "======================================================================\n");CHKERRQ(ierr);
  cHamiltonianMatrix Hpolaron;
  Hpolaron.input();
  Hpolaron.fock();
  Hpolaron.timeEvolutaion();
  Hpolaron.destruction();
  ierr = SlepcFinalize();
  return 0;
}
