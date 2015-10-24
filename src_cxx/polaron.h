#ifndef POLARON_H
#define POLARON_H
/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners

*/
#include <slepceps.h>
//#include <petscksp.h>
//#include <petscsys.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <cmath> 
//#include <ctgmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_bessel.h>
#include <sys/time.h>
using namespace std;
const PetscInt __MAXNOZEROS__ = 100; // TODO: This is the max number in a row --> theoretically largest recursion relation index given by the Hamiltonian.

class cHamiltonianMatrix{
private:
  Vec            X1,X2,X3;
  Mat            Hpolaron;
  Vec*			 WFt;
  EPS			 eps;
  EPSType		 type;
  gsl_matrix     *basis1, *basis2;
  gsl_vector	 *randV;
  PetscInt       ROW,COLUMN,rstart,rend,nlocal,col[__MAXNOZEROS__],nev,its,maxit,nconv;
  PetscScalar    value[__MAXNOZEROS__],HpolaronMax,HpolaronMin;
  double		 a_scaling, b_scaling;
  int			 set_gsl_under_flow_ratio;
  gsl_vector *	 rr;
protected:
  PetscErrorCode ierr;
  PetscInt       N,N2,L,position;
  PetscInt		 judge, boundary; // judge=0;	% Clean background!---This background serves as a heat bath!
  	  	  	  	  	  	  	  	  // judge=1;    % Disordered background!---'energy' bath due to interaction
  	  	  	  	  	  	  	  	  // boundary=1; % @Open Boundary
  	  	  	  	  	  	  	  	  // boundary=0; % @Periodic Boundary
  PetscInt		 dim, dim2,DIM;   // Derived parameters: too large L or N will give non-number: NaN or Inf.
  double	     W,U,tmax, dt,tol,error,re,im;
  PetscMPIInt    rank, size;
  int 			 _jdim1, _jdim2,Nt;
public:
  cHamiltonianMatrix(){}
  ~cHamiltonianMatrix(){}
  PetscErrorCode destruction();
  PetscErrorCode input();
  PetscErrorCode fock();
  void construct_basis(gsl_matrix *, const int);
  PetscErrorCode timeEvolution();
  void randomPotential(gsl_vector*);
  PetscErrorCode hamiltonianConstruction();
  PetscErrorCode hamiltonianRescaling();
  PetscErrorCode assemblance();
  void block(int, int &, int &);
  PetscReal compute_diag();
  unsigned long int random_seed();
  unsigned long int random_seed_devrandom();
  void compute_hopping(int &,const int);
  void spin_flag_condition(int &, int &, const int);
  void compute_col_index(const int, const int, const int );
  int myindex(gsl_vector *, const int);
  PetscErrorCode initial_state();
  PetscErrorCode KernalPolynomialMethod();
  PetscErrorCode measurement();
  PetscErrorCode WaveFunctionUpdate(const int);
};
#endif
