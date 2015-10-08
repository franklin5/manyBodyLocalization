
static char help[] = "Localiazation and delocalization transition of a polaron on a 1D disordered lattice.\n\n";

/*T
   Concepts: KSP^basic parallel example; Sparse matrix construction;
   Processors: n
   mpirun -n 6 SteadyState -ksp_monitor_short   -pc_type jacobi   -ksp_type gmres -ksp_gmres_restart 200
   or
   ./SteadyState -ksp_monitor_short   -pc_type jacobi   -ksp_type gmres -ksp_gmres_restart 200
T*/

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

*/
#include <petscksp.h>
#include "polaron.h"
#undef __FUNCT__
#define __FUNCT__ "main"
#define root 0
int main(int argc,char **args){
  PetscErrorCode ierr;
  PetscInitialize(&argc,&args,(char*)0,help);
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
  ierr = PetscFinalize();
  return 0;
}
