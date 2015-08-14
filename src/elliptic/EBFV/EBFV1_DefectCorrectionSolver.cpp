/*
 * EBFV1__DefectCorrectionSolver.cpp
 *
 *  Created on: 25/08/2012
 *      Author: rogsoares
 */

#include <time.h>
#include "EBFV/EBFV1_elliptic.h"

namespace PRS           // PRS: Petroleum Reservoir Simulator
{

	double EBFV1_elliptic::solveIteratively(){
//		PetscLogDouble flops1;
//		PetscGetFlops(&flops1);
//		//	cout << "Mflops: " << flops1/1.0e6 << endl;
//		/*
//		 * Iterative stop criteria: relative tolerance between two consecutive solutions
//		 * Should not be confused with krylov relative tolerance.
//		 */
//		double res = 1e+8, res_old, tol = pSimPar->rtol2();
//		double norm1, norm2;
//
//		/*
//		 * Output file to monitoring solver behavior using defect correction (DC) approach
//		 */
//	#ifdef USE_DEFECT_CORRECTION_MONITOR
//		ofstream fid;
//		static int step = 0;
//		if (!P_pid()){
//			char filename[256]; sprintf(filename,"Solver_monitor__DC__step_%d.txt",++step);
//			static bool key = true;
//			if (key){
//				key = false;
//				//	fid.open(filename);
//				//if (!P_pid()) cout << "\niter    iter.    residuum      elapsed_time  total.\n";
//			}
//			//fid << "\n\nStep: " << step++ << endl;
//		}
//	#endif // USE_DEFECT_CORRECTION_MONITOR
//
//		/*
//		 * u_old: 	solution vector at k-1 iteration
//		 * v_tmp: 	temporary vector
//		 * resvec: 	difference between solution k and k-1
//		 * rhs:		right-hand-side vector (boundary condition of first and second order)
//		 */
//		//Vec u_old, v_tmp, resvec, rhs;
//
//		int ndom = pSimPar->getNumDomains();	// number fo domains
//		PetscInt its;							// counts kryolv iterations
//		PetscInt dcits = 0;							// counts DC iterations
//
//		// Create vectors
//		if (DF_key){
//			ierr = VecCreate(PETSC_COMM_WORLD,&pVec->u_old);                         CHKERRQ(ierr);
//			ierr = VecSetSizes(pVec->u_old,PETSC_DECIDE,pMData->getNum_GF_Nodes() ); CHKERRQ(ierr);
//			ierr = VecSetFromOptions(pVec->u_old);                                   CHKERRQ(ierr);
//
//			ierr = VecDuplicate(pVec->u_old,&pVec->resvec);								   CHKERRQ(ierr);
//			ierr = VecDuplicate(pVec->u_old,&pVec->v_tmp);								   CHKERRQ(ierr);
//			ierr = VecDuplicate(pVec->u_old,&pVec->rhs);								   CHKERRQ(ierr);
//
//			ierr = VecZeroEntries(pVec->u_old);										CHKERRQ(ierr);
//			ierr = VecZeroEntries(pVec->resvec);                                      CHKERRQ(ierr);
//			ierr = VecZeroEntries(pVec->v_tmp);                                      CHKERRQ(ierr);
//			ierr = VecZeroEntries(pVec->rhs);                                      CHKERRQ(ierr);
//			DF_key = false;
//		}
//		else{
//			ierr = VecZeroEntries(pVec->u_old);                                      CHKERRQ(ierr);
//			ierr = VecZeroEntries(pVec->v_tmp);                                      CHKERRQ(ierr);
//			ierr = VecZeroEntries(pVec->rhs);                                      CHKERRQ(ierr);
//
//		}
//		// u_new = u_old = 0
//		ierr = VecZeroEntries(output);CHKERRQ(ierr);
//
//		/*----------------------------------------------------------------------
//		 * "Defect-correction" like procedure to compute pressure field:
//		 *
//		 * 		RHS = RHS - (E * F) * u^n
//		 * 		G * u^(n+1) = RHS
//		 *----------------------------------------------------------------------*/
//		PetscTruth flg = PETSC_FALSE;
//		bool key = true;
//		double sum = .0;
//		int I=0;
//		while ( res > tol ){
//			double startt = MPI_Wtime();
//			PetscReal nrm;
//
//			ierr = VecZeroEntries(pVec->v_tmp);CHKERRQ(ierr);
//			for (int i=0; i<ndom; i++){
//				ierr = MatMultAdd(EF_multiDom[i],pVec->u_old,pVec->v_tmp,pVec->v_tmp);CHKERRQ(ierr);
//			}
//			double t1 = MPI_Wtime();
//
//			// w = alpha x + y.
//			// rhs = RHS - v_tmp
//			// -> rhs(w) = alpha(-1)v_tmp(x) + RHS(y)
//			// w,alpha,x,y -> rhs,-1,v_tmp,RHS
//			ierr = VecWAXPY(pVec->rhs, /* w */
//					-1.0, /* alpha */
//					pVec->v_tmp, /* x */
//					matvec_struct->RHS);/* y */
//			CHKERRQ(ierr);
//			double t2 = MPI_Wtime();
//
//			// G * u^(n+1) = RHS
//			KSP_solver(matvec_struct->G,
//					matvec_struct->G,
//					pVec->rhs,
//					output,
//					pSimPar,
//					flg,
//					KSPCG,pSimPar->getPCType(),its);CHKERRQ(ierr);
//			double t3 = MPI_Wtime();
//
//			flg = PETSC_TRUE;
//			// resvec = u_new - u_old.
//			ierr = VecWAXPY(pVec->resvec,-1.0,pVec->u_old,output);CHKERRQ(ierr);
//			ierr = VecNorm(pVec->resvec,NORM_2,&norm1);CHKERRQ(ierr);
//			ierr = VecNorm(output,NORM_2,&norm2);CHKERRQ(ierr);
//			double t4 = MPI_Wtime();
//
//			res_old = res;
//			res = (double)norm1/norm2;
//
//			if (res > res_old && !P_pid()){
//				cout << "WARNING: solver is diverging!\n";
//				static int div_status = 0; div_status++;
//				//if (div_status > 10)
//				//throw Exception(__LINE__,__FILE__,"Solver failured for not converging.");
//			}
//
//			// u_old = u_new
//			ierr = VecCopy(output,pVec->u_old);CHKERRQ(ierr);
//
//			double endt = MPI_Wtime();
//			double elapsed_time = endt-startt;
//			sum += elapsed_time;
//
//
//	#ifdef USE_DEFECT_CORRECTION_MONITOR
//			if (!P_pid()){
//				//	cout << ++dcits << "  " << its << " " << res << "  " << elapsed_time << "  "  << sum << "\n";
//			}
//	#endif // USE_DEFECT_CORRECTION_MONITOR
//		}
//		//		VecView(output,PETSC_VIEWER_STDOUT_WORLD);
//		//		PetscReal val;
//		//		VecNorm(output,NORM_2,&val);
//		//		PetscPrintf(PETSC_COMM_WORLD,"norm: %f",val);
//		//		STOP();
//
//		PetscLogDouble flops2;
//		PetscGetFlops(&flops2);
//		//cout << "Mflops: " << flops2/1.0e6 << endl;
//		//cout << "Mflops: " << (flops2-flops1)/1.0e6 << endl;
//
//
//		// free memory
//		for (int i=0; i<ndom; i++){
//			ierr = MatDestroy(EF_multiDom[i]);CHKERRQ(ierr);
//		}
//		delete [] EF_multiDom;
//		return sum;
		return 0;
	}
}
