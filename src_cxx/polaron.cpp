#include <petscksp.h>
#include "polaron.h"
#undef __FUNCT__
#define __FUNCT__ "polaron"

PetscErrorCode cHamiltonianMatrix::input(){
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
//	cout << rank << '\t' << size << endl;
	for (int ig = 0; ig < size; ++ig) {
	    if (ig ==rank){
	    	char dummyname[100];
	    	double dummyvalue;
	    	int intdummyvalue;
	        FILE *input;
	        input = fopen("input.txt","r");
	        assert(input != NULL);
	        if (ig == 0)  {
	        	cout << "Starting to read in parameters from file input.txt" << endl;
			    cout << " --- Time Dependent Evolution of Disordered Impurity fermion ---" << endl;
			    cout << "Parameters of System:" << endl;
			    cout << "---------------------------------------------------------- " << endl;
	        }
	        fscanf(input,"%s %d", dummyname, &intdummyvalue);
	        L = intdummyvalue;    if (ig == 0) cout << "SystemSize " << dummyname << "=" << L << endl;
	        fscanf(input,"%s %d", dummyname, &intdummyvalue);
	        N = intdummyvalue;    if (ig == 0) cout << "Majority fermion number " << dummyname << "=" << N << endl;
	        // N = L/2; // --> L has to be even number here!!! TODO: odd number case?
//	        cout << "Number of Particle N=" << N << endl;
	        fscanf(input,"%s %d", dummyname, &intdummyvalue);
	        N2 = intdummyvalue;    if (ig == 0) cout << dummyname << "=" << N2 << endl;
	        fscanf(input,"%s %lf", dummyname, &dummyvalue);
	        dim = PetscInt ( gsl_sf_choose(L,N));
	        dim2= PetscInt ( gsl_sf_choose(L,N2));
	        DIM = dim*dim2;
	        if (ig == 0) cout << "Dimension of Hilbert Space: " << dim*dim2 << endl;
	        W = dummyvalue;    if (ig == 0) cout << "Disorder " << dummyname << "=" << W << endl;
	        fscanf(input,"%s %lf", dummyname, &dummyvalue);
	        U = dummyvalue;    if (ig == 0) cout << "Interaction" << dummyname << "=" << U << endl;
	        fscanf(input,"%s %lf", dummyname, &dummyvalue);
	        tmax = dummyvalue;    if (ig == 0) cout << "Total time tmax=" << dummyname << "=" << tmax << endl;
	        fscanf(input,"%s %d", dummyname, &intdummyvalue);
	        Nt = intdummyvalue;    if (ig == 0) cout << dummyname << "=" << Nt << endl;
	    	dt=log(tmax)/Nt; if (ig == 0) cout << "time steps in log time scale = " << dt << endl;
	        fscanf(input,"%s %d", dummyname, &intdummyvalue);
	        judge = intdummyvalue;    if (ig == 0) cout << dummyname << "=" << judge << endl;
	        fscanf(input,"%s %d", dummyname, &intdummyvalue);
	        boundary = intdummyvalue;    if (ig == 0) cout << dummyname << "=" << boundary << endl;
	        fscanf(input,"%s %d", dummyname, &intdummyvalue);
	        position = intdummyvalue;    if (ig == 0) cout << dummyname << "=" << position << endl;
	        fclose(input);
	        if (ig == 0) {
	        	if (judge==0) {
					cout << "Clean background!---This background serves as a heat bath!  " << endl;
				} else {
					cout << "Disordered background!---'energy' bath due to interaction" << endl;
				}
				if (boundary==0) {
					cout << "#Periodic Boundary Conditions" << endl;
				} else {
					cout << "#Open Boundary Conditions" << endl;
				}
				cout << "---------------------------------------------------------- " << endl;
	        }
	    }
	}
    return ierr;
}

PetscErrorCode cHamiltonianMatrix::fock(){
	int spin_flag;
	spin_flag = 1;
	basis1 = gsl_matrix_alloc(N,dim);
	construct_basis(basis1, spin_flag);
	spin_flag = 2;
	basis2 = gsl_matrix_alloc(N2,dim2);
	construct_basis(basis2, spin_flag);
		if (rank == 0){ // only root CPU prints the result to screen.
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < dim; ++j) {
					cout << gsl_matrix_get(basis1,i,j) << '\t';
				}
				cout << "EOL" << endl;
			}
			for (int i = 0; i < N2; ++i) {
				for (int j = 0; j < dim2; ++j) {
					cout << gsl_matrix_get(basis2,i,j) << '\t';
				}
				cout << "EOL" << endl;
			}
		}
	return ierr;
}

PetscErrorCode cHamiltonianMatrix::timeEvolutaion(){
	randV = gsl_vector_alloc(L);
	randomPotential(randV);
	for (int ig = 0; ig < size; ++ig) { // I have to make sure the generated random number potential is the same across all CPUs.
		    if (ig ==rank){
		    	cout << "rank " << rank << " has rand numbers:" << '\t';
				cout.precision(16);
		    	for (int i = 0; i < L; i++) {
					cout << gsl_vector_get(randV,i) << '\t';
				}
				cout << endl;
		    }
	}

	ierr = hamiltonianConstruction();CHKERRQ(ierr);

	ierr = hamiltonianRescaling();CHKERRQ(ierr);

	// initial state preperation
	ierr = initial_state();CHKERRQ(ierr);
	//	%==============================================================
		//	    % -------- time-dependent evolution --------------
		//	    %-------------------------------------------------
	// Kernal Polynomial Method

	return ierr;
}

PetscErrorCode cHamiltonianMatrix::hamiltonianRescaling(){

//	ierr = MatCreateVecs(Hpolaron,NULL,&xr);CHKERRQ(ierr);
//	ierr = MatCreateVecs(Hpolaron,NULL,&xi);CHKERRQ(ierr);
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	ierr = EPSSetOperators(eps,Hpolaron,NULL);CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
	PetscBool her;
	ierr = EPSIsHermitian(eps, &her); CHKERRQ(ierr); if (rank==0) {cout << endl;cout << "is hermitian? " << her << endl;}
	PetscBool pos;
	ierr = EPSIsPositive(eps, &pos); CHKERRQ(ierr); if (rank==0) {cout << "is positive? " << pos << endl;cout << endl;}
//	ierr = EPSSetType(eps,EPSPOWER);CHKERRQ(ierr);
	ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr);
//	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
	ierr = EPSSolve(eps);CHKERRQ(ierr);
//	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
//	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
//	ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
//	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);
	// ------- maximum and minimum of energy spectrum ----------
	ierr = EPSGetEigenpair(eps,0,&HpolaronMax,NULL,NULL,NULL);CHKERRQ(ierr);
	ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
	ierr = EPSSolve(eps);CHKERRQ(ierr);
	ierr = EPSGetEigenpair(eps,0,&HpolaronMin,NULL,NULL,NULL);CHKERRQ(ierr);
	double epsilon_cut_off = 0.01;
	a_scaling = PetscRealPart(HpolaronMax-HpolaronMin)/(2.0-epsilon_cut_off);
	b_scaling = PetscRealPart(HpolaronMax+HpolaronMin)/2.0;
//	if (rank==0) {
//		cout << "HpolaronMax is " << HpolaronMax << endl;
//		cout << "HpolaronMin is " << HpolaronMin << endl;
//		cout << "a_rescaling is " << a_rescaling << endl;
//		cout << "b_rescaling is " << b_rescaling << endl;
//	}
	/*
	% --------- rescaled Hamiltonian -------------
	 */
	ierr = MatShift(Hpolaron,-b_scaling);CHKERRQ(ierr);
	ierr = MatScale(Hpolaron,1.0/a_scaling);
//	ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_DENSE  );CHKERRQ(ierr);
//	//      ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_MATLAB  );CHKERRQ(ierr);
//	ierr = MatView(Hpolaron,	PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr);

//	ierr = EPSSetOperators(eps,Hpolaron,NULL);CHKERRQ(ierr);
//	ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
//	ierr = EPSSolve(eps);CHKERRQ(ierr);
//	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
//	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
//	ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
//	ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//	ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	if (rank==0)    cout << "Characteristic Coherance Time: t_c=" << 1.0/(2.0*a_scaling/dim/L) << endl;
	/*
	 %---------------------------------------------------
	    % -------- End of Cosntruction of Hamiltonian ------
	    %=============================================================
	    */
	return ierr;
}

PetscErrorCode cHamiltonianMatrix::measurement(){

	return ierr;
}

PetscErrorCode cHamiltonianMatrix::hamiltonianConstruction(){
	ierr = MatCreate(PETSC_COMM_WORLD,&Hpolaron);CHKERRQ(ierr);
	  ierr = MatSetType(Hpolaron,MATMPIAIJ);CHKERRQ(ierr);
	  ierr = MatSetSizes(Hpolaron,PETSC_DECIDE,PETSC_DECIDE,DIM,DIM);CHKERRQ(ierr);
	  // TODO: should be able to set the symmetric/hermitian option and
	  // only do upper-right triangle part of matrix construction .
	  // and perform corresponding operations thereon.
	  // ierr = MatSetOption(Hpolaron,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
	  // ierr = MatSetOption(Hpolaron,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);

	  // TODO: what is the estimate of the pre-allocation?
	  // -- number of nonzeros per row in DIAGONAL portion of local submatrix
	  // (same value is used for all local rows) ? I put dim temporarily here.
	  // number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix
	  // (same value is used for all local rows) ?  I put dim temporarily here..
	  // More details at http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatMPIAIJSetPreallocation.html
	  ierr = MatMPIAIJSetPreallocation(Hpolaron,DIM,NULL,DIM,NULL);CHKERRQ(ierr);
	  ierr = MatSeqAIJSetPreallocation(Hpolaron,DIM,NULL);CHKERRQ(ierr);

	  ierr = MatGetOwnershipRange(Hpolaron,&rstart,&rend);CHKERRQ(ierr);
	  ierr = MatGetLocalSize(Hpolaron,&nlocal, NULL);CHKERRQ(ierr);

	  ierr = assemblance();CHKERRQ(ierr);

	  ierr = MatAssemblyBegin(Hpolaron,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Hpolaron,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

//	  ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_DENSE  );CHKERRQ(ierr);
////      ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_MATLAB  );CHKERRQ(ierr);
//	  ierr = MatView(Hpolaron,	PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr);

	  return ierr;
}

PetscErrorCode cHamiltonianMatrix::assemblance(){
	/*
	    Assemble matrix.
	*/
	int nonzeros; // TODO: check if nonzeros < __MAXNOZEROS__ is true.
	  int spin_flag;
	  PetscReal _val_;
	  for (ROW=rstart; ROW<rend; ROW++) {
		  nonzeros = 0;

//		  cout << "row is " << ROW << endl;

		  block(ROW,_jdim1,_jdim2);
//		  cout << _jdim1 << '\t' << _jdim2 << endl;
		  // random potential and interaction terms are diagonal.
		  col[nonzeros] = ROW;value[nonzeros] = compute_diag();nonzeros++;
		  // hopping terms are off-diagonal.
		  // --------- hopping of spin down -------------
		  spin_flag = 2; // to be consistent with convention in construct_basis member function.
		  compute_hopping(nonzeros,spin_flag);
		  // ----------- hopping of spin up ---------------
		  spin_flag = 1;
		  compute_hopping(nonzeros,spin_flag);
		 if (nonzeros > __MAXNOZEROS__){
	        	cerr << "nonzeros on a row " <<  nonzeros << " is larger than the pre-allocated range of"
	        	<<  __MAXNOZEROS__ <<" const arrays. Try increasing the max number in polaron.h" << endl;
	        	exit(1);
	        } else {
	        	ierr   = MatSetValues(Hpolaron,1,&ROW,nonzeros,col,value,INSERT_VALUES);CHKERRQ(ierr);
	        }
	  }
	return ierr;
}

int cHamiltonianMatrix::myindex(gsl_vector * site, const int _N_){
//	for (int var = 0; var < site->size; ++var) {
//		cout << gsl_vector_get(site,var) << endl;
//	}
//	cout << "total size of site is " << site->size << endl;
//	cout << "_N_ value is " << _N_ << endl;
	int p = 0;
	if (_N_ == 1) {
		p = gsl_vector_get(site,0);
	} else {
		int loop = 1;
		for (int i = 1; i <= gsl_vector_get(site,loop-1)-1; ++i) {
			p += gsl_sf_choose(L-i,_N_-loop);
		}
		loop ++;
		do {
			for (int i = gsl_vector_get(site,loop-2)+1; i < gsl_vector_get(site,loop-1)-1; ++i) {
				p += gsl_sf_choose(L-i,_N_-loop);
			}
			loop ++;
		} while (loop<_N_);
		p += gsl_vector_get(site,_N_-1)-gsl_vector_get(site,_N_-2); // % that's fine!
	}

//    cout << "q is " << p << endl;


	return p;
}

void cHamiltonianMatrix::compute_hopping(int &nonzeros,const int spin_flag){
	int _N, _dim, q, jdim;
	spin_flag_condition(_N,_dim,spin_flag);

	gsl_matrix *basis;
	if (spin_flag == 1) { // spin up
		basis = basis1; // (value of pointer is copied)
		jdim = _jdim1;
	} else if (spin_flag == 2) { // spin down
		basis = basis2; // (value of pointer is copied)
		jdim = _jdim2;
	} else {
		cerr << "Wrong spin flag configuration. Check construct_basis call in fock member function." << endl;
		exit(1);
	}
	gsl_vector * site = gsl_vector_alloc(_N); // TODO: why in construct basis member function, the alloc size is _N+1?

	for (int jpar = 1; jpar <= _N; ++jpar) {
	// ------------- Start: Boundary term -------------
		if (boundary==0){ // #Periodic Boundary Conditions
			// -------- boundary terms (at site L) -----------------
			if (gsl_matrix_get(basis,jpar-1,jdim)==L  && gsl_matrix_get(basis,0,jdim)!=1) {
				if (_N==1) {
					q = 1;
				} else {
					/*
					 *      nnn=2:N;
                            site=basis1(:,jdim1);
                            site(nnn)=site(nnn-1);
                            site(1)=1;
					 */
					for (int nnn = 1; nnn < _N; ++nnn) {
						gsl_vector_set(site,nnn,gsl_matrix_get(basis,nnn-1,jdim));
					}
					gsl_vector_set(site,0,1);
					q = myindex(site,_N);
				}
				compute_col_index(nonzeros,spin_flag,q);nonzeros++;
			}
			// -------- boundary terms (at site 1) -----------------
			if (gsl_matrix_get(basis,jpar-1,jdim)==1  && gsl_matrix_get(basis,_N-1,jdim)!=L) {
				if (_N==1) {
					q = L;
				} else {
				/*
				 * nnn=2:N;
					site=basis1(:,jdim1);
					site(nnn-1)=site(nnn);
					site(N)=L;
				 */
					for (int nnn = 1; nnn < _N; ++nnn) {
						gsl_vector_set(site,nnn-1,gsl_matrix_get(basis,nnn,jdim));
					}
					gsl_vector_set(site,_N-1,L);
					q = myindex(site,_N);
				}
				compute_col_index(nonzeros,spin_flag,q);nonzeros++;
			}
		} //else { // TODO: Do we consider boundary==1 case or other?
	// ------------- End: Boundary term -------------.
			if ( (jpar< _N && gsl_matrix_get(basis,jpar-1,jdim)+1<gsl_matrix_get(basis,jpar,jdim)) || (jpar==_N && gsl_matrix_get(basis,jpar-1,jdim)<L) ) {
/*
 *                     site=basis1(:,jdim1);
                    site(jpar)=site(jpar)+1;
 */
//				cout << "i am the first condition" << endl;
//				cout << "jpar = " << jpar << " jdim = " << jdim << endl;

				for (int nnn = 0; nnn < _N; ++nnn) {
					gsl_vector_set(site,nnn,gsl_matrix_get(basis,nnn,jdim));
				}
				gsl_vector_set(site,jpar-1,gsl_vector_get(site,jpar-1)+1);

//				cout << "site is" << endl;
//				for (int var = 0; var < _N; ++var) {
//					cout << gsl_vector_get(site,var) << '\t';
//				}
//				cout << endl;

                q=myindex(site,_N);
                compute_col_index(nonzeros,spin_flag,q);nonzeros++;
			}
			if ((jpar==1 && gsl_matrix_get(basis,jpar-1,jdim)>1) || (jpar>1 && gsl_matrix_get(basis,jpar-1,jdim)-1>gsl_matrix_get(basis,jpar-2,jdim))) {
			/*
			 *                     site=basis1(:,jdim1);
                    site(jpar)=site(jpar)-1;
			 */
//				cout << "i am the second condition" << endl;
//				cout << "jpar = " << jpar << " jdim = " << jdim << endl;

				for (int nnn = 0; nnn < _N; ++nnn) {
					gsl_vector_set(site,nnn,gsl_matrix_get(basis,nnn,jdim));
				}
				gsl_vector_set(site,jpar-1,gsl_vector_get(site,jpar-1)-1);

//				cout << "site is" << endl;
//				for (int var = 0; var < _N; ++var) {
//					cout << gsl_vector_get(site,var) << '\t';
//				}
//				cout << endl;

				q=myindex(site,_N);
				compute_col_index(nonzeros,spin_flag,q);nonzeros++;
			}
		//}
	}
	gsl_vector_free (site);
}

void cHamiltonianMatrix::compute_col_index(const int nonzeros, const int spin_flag, const int q){
	value[nonzeros] = -1.0; // -t term. hopping value, t, is taken as the unit.
	if (spin_flag == 1) { // spin up
		col[nonzeros] = (_jdim2)*dim+q-1;
//		cout << "spin up" << endl;
	} else if (spin_flag == 2) { // spin down
		col[nonzeros] = (q-1)*dim+_jdim1;
//		cout << "spin down" << endl;
	} else {
		cerr << "Wrong spin flag configuration. Check construct_basis call in fock member function." << endl;
		exit(1);
	}
//	cout << "col is " << col[nonzeros] << endl; // << " and value is " << value[nonzeros]
}

void cHamiltonianMatrix::spin_flag_condition(int & _N, int & _dim, const int spin_flag){
	if (spin_flag == 1) { // spin up
		_N = N;
		_dim = dim;
	} else if (spin_flag == 2) { // spin down
		_N = N2;
		_dim = dim2;
	} else {
		cerr << "Wrong spin flag configuration. Check construct_basis call in fock member function." << endl;
		exit(1);
	}
}

PetscReal cHamiltonianMatrix::compute_diag(){
	PetscReal sumAA=0, tempR;
	if (judge==1) {
		for (int jpar = 1; jpar <= N; ++jpar) {
			tempR = gsl_matrix_get(basis1,jpar-1,_jdim1); // cout << tempR << endl;
			tempR = tempR-1; // potential fail of the index 0 or 1 bug. Temporary fix.
			sumAA += gsl_vector_get(randV,int(round(tempR)));
		}
//	cout << sumAA << endl;
	}
	for (int jpar2 = 1; jpar2 <= N2; ++jpar2) {
		sumAA += gsl_vector_get(randV,gsl_matrix_get(basis2,jpar2-1,_jdim2)-1);
	}
//	cout << sumAA << endl;
	for (int jpar = 1; jpar <= N; ++jpar) {
		for (int jpar2 = 1; jpar2 <= N2; ++jpar2) {
			if (gsl_matrix_get(basis1,jpar-1,_jdim1)==gsl_matrix_get(basis2,jpar2-1,_jdim2)) {
				sumAA += U;
			}
		}
	}
//	cout << sumAA << endl;
	return sumAA;
}

void cHamiltonianMatrix::block(int irow,int &jdim1, int &jdim2){
	jdim2 = floor(irow/dim);
	jdim1 = irow - jdim2*dim;
}

void cHamiltonianMatrix::construct_basis(gsl_matrix *basis, const int spin_flag){
	int _N, _dim, loop, sumAA, len, num, sum2, index_loop, index_p, index_j;

	spin_flag_condition(_N,_dim,spin_flag);

	gsl_vector * site = gsl_vector_alloc(_N+1);
	for (int p = 1; p <= _dim; ++p) {
		index_p = p-1; // MATLAB vector and matrix index starts from 1 while C/C++ index starts from 0. Be careful about the 0 and 1 bug.
		loop = 1; index_loop = loop-1;
		if (_N==1) {
			loop++; index_loop = loop-1;
			gsl_vector_set(site,index_loop,p);
		} else {
			sumAA = 0;
			do {
				len = L-gsl_vector_get(site,index_loop);
				num = _N-loop;
				for (int jcar = len-1; jcar >= num; jcar=jcar-1) {
					sum2 = sumAA;
					sumAA += gsl_sf_choose(jcar,num);
					if (sumAA >= p) {
						loop++; index_loop = loop-1;
						sumAA=sum2;
						gsl_vector_set(site, index_loop, gsl_vector_get(site,index_loop-1)+len-jcar );
						break;
					}
				}
			} while (loop < _N);
			gsl_vector_set(site, index_loop+1, gsl_vector_get(site,index_loop)+p-sum2);
		}
		for (int j = 1; j <= _N; ++j) {
			index_j = j-1;
			gsl_matrix_set(basis, index_j, index_p, gsl_vector_get(site,index_j+1));
		}
	}
	gsl_vector_free (site);
	// Debugging purpose: output.
//	if (rank == 0){ // only root CPU prints the result to screen.
//		for (int i = 0; i < _N; ++i) {
//			for (int j = 0; j < _dim; ++j) {
//				cout << gsl_matrix_get(basis,i,j) << '\t';
//			}
//			cout << "EOL" << endl;
//		}
//	}
}

PetscErrorCode cHamiltonianMatrix::initial_state(){
	Mat H0;
	PetscScalar kr,ki;
	Vec            xr,xi;          /* real and imag part of eigenvectors */
	ierr = MatCreate(PETSC_COMM_WORLD,&H0);CHKERRQ(ierr);
	ierr = MatSetType(H0,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatSetSizes(H0,PETSC_DECIDE,PETSC_DECIDE,L,L);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(H0,2,NULL,2,NULL);CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(H0,2,NULL);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(H0,&rstart,&rend);CHKERRQ(ierr);
	int nonzeros;
	for (ROW = rstart; ROW < rend; ++ROW) {
		nonzeros = 0;
		if (ROW==0) {
			col[nonzeros] = 1; value[nonzeros] = -1.0; nonzeros++;
		} else if (ROW==L-1) {
			col[nonzeros] = L-2; value[nonzeros] = -1.0; nonzeros++;
		} else {
			col[nonzeros] = ROW-1; value[nonzeros] = -1.0; nonzeros++;
			col[nonzeros] = ROW+1; value[nonzeros] = -1.0; nonzeros++;
		}
		if (nonzeros > 2) {
			cerr << "single particle H0 is not constructed properly" << endl;exit(1);
		} else {
			ierr   = MatSetValues(H0,1,&ROW,nonzeros,col,value,INSERT_VALUES);CHKERRQ(ierr);
		}
	}
	ierr = MatAssemblyBegin(H0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(H0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatCreateVecs(H0,NULL,&xr);CHKERRQ(ierr);
//	ierr = MatCreateVecs(H0,NULL,&xi);CHKERRQ(ierr);
	ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_DENSE  );CHKERRQ(ierr);
	//      ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,	PETSC_VIEWER_ASCII_MATLAB  );CHKERRQ(ierr);
	ierr = MatView(H0,	PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr);

	ierr = EPSSetOperators(eps,H0,NULL);CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);
	ierr = EPSSetDimensions(eps,L,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = EPSSolve(eps);CHKERRQ(ierr);
	ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
//	if (nconv>0) { // The printout tells us, in this case of matrix, the eigenvalue is in ascending order and eigenvectors are normalized (well, at least close to being normalized).
// TODO: if it is not the case, use GSL sorting of eigenvalues and eigenvectors, https://www.gnu.org/software/gsl/manual/html_node/Sorting-Eigenvalues-and-Eigenvectors.html#Sorting-Eigenvalues-and-Eigenvectors
	//		for (int i=0;i<nconv;i++) {
//			ierr = EPSGetEigenpair(eps,i,&kr,NULL,xr,NULL);CHKERRQ(ierr);
//			if (rank==0) {
//				cout << i << "th eigenvalue is " << PetscRealPart(kr) << endl;
//				ierr = VecView(xr,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//				ierr = VecView(xi,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//			}
//		}
//	}
	// %initial state---zero temperature.
	ierr = VecCreate(PETSC_COMM_WORLD,&X1);CHKERRQ(ierr);
	ierr = VecSetSizes(X1,nlocal,DIM);CHKERRQ(ierr);
	ierr = VecSetFromOptions(X1);CHKERRQ(ierr);
	ierr = VecDuplicate(X1,&X2);CHKERRQ(ierr);

	gsl_vector * site2 = gsl_vector_alloc(N2);
	gsl_vector_set(site2,0,position);
	if (N2>1) for (int i = 1; i < N2; ++i) 	gsl_vector_set(site2,i,0);
	int p2 = myindex(site2,N2); //%index of impurity particle, also the position on the lattice
	double val_det;
	gsl_matrix_complex *slater = gsl_matrix_complex_alloc(N,N);
	gsl_complex val;
	for (int jdim = 0; jdim < dim; ++jdim) {
		for (int j1 = 0; j1 < N; ++j1) { // row
			for (int j2 = 0; j2 < N; ++j2) { // col
				ierr = EPSGetEigenpair(eps,j2,&kr,NULL,xr,NULL);CHKERRQ(ierr);
				VecGetValues(xr,1,) // VecGetValues can only  get values on the same processor ... Let me abandon this and move to gsl instead...
				gsl_matrix_get(basis1,j1-1,jdim),j2

				gsl_matrix_complex_set(slater,j1,j2,val);
			}
		}

		ierr = VecSetValue(X1,(p2-1)*dim+jdim,&val_det, INSERT_VALUES);
	}


	ierr = MatDestroy(&H0);CHKERRQ(ierr);
	ierr = VecDestroy(&xr);CHKERRQ(ierr);
//	ierr = VecDestroy(&xi);CHKERRQ(ierr);
	gsl_vector_free(site2);
	gsl_matrix_free(slater);
	return ierr;
}

PetscErrorCode cHamiltonianMatrix::destruction(){

//	  cout << "Object is being deleted" << endl;
	  gsl_matrix_free(basis1);
	  gsl_matrix_free(basis2);
	  gsl_vector_free(randV);
	  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	  ierr = MatDestroy(&Hpolaron);CHKERRQ(ierr);
	  ierr = VecDestroy(&X1);CHKERRQ(ierr);
	  ierr = VecDestroy(&X2);CHKERRQ(ierr);
//	  cout << " do i have a seg fault?" << endl;
	return ierr;
}

void cHamiltonianMatrix::randomPotential(gsl_vector* randV){

	const gsl_rng_type * T;
	  gsl_rng * r;

	  gsl_rng_env_setup();

	  T = gsl_rng_default;
	  r = gsl_rng_alloc (T);
	  // gsl_rng_env_setup() is getting a generator type and seed from environment variables.
	  // http://stackoverflow.com/questions/9768519/gsl-uniform-random-number-generator
	  // dev/random solution is very time consuming compared to the gettimeofday(),
	  // the gettimeofday() solution, might be better its level of accuracy is enough:
	  unsigned long int seed_number;
	  if (rank == 0) {
		  seed_number = random_seed();
//		  cout << "Before broadcasting, seed number is " << seed_number << endl;
	  }
	  MPI_Bcast(&seed_number, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD); //replace MPI_INT by MPI_UNSIGNED_LONG if int is not long enough.
//	  cout << "After broadcasting, rank " << rank << " has seed number " << seed_number << endl;
	    gsl_rng_set(r,seed_number);
//	    gsl_rng_set(r,random_seed_devrandom());

	  for (int i = 0; i < L; i++)
	    {
		  gsl_vector_set(randV, i, -W/2+W*gsl_rng_uniform (r));
	    }

	  gsl_rng_free (r);
}

unsigned long int cHamiltonianMatrix:: random_seed()
{
  struct timeval tv;
  gettimeofday(&tv,0);
  return (tv.tv_sec + tv.tv_usec);
}

unsigned long int cHamiltonianMatrix:: random_seed_devrandom()
{

  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;

  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
  } else {
    fread(&seed,sizeof(seed),1,devrandom);
    fclose(devrandom);
  }

  return(seed);

}
