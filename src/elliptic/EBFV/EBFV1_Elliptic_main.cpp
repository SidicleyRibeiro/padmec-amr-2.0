#include <time.h>
#include "EBFV/EBFV1_elliptic.h"

namespace PRS{
	EBFV1_elliptic::EBFV1_elliptic(){
	}
	
	EBFV1_elliptic::EBFV1_elliptic(pMesh mesh, PhysicPropData *ppd, SimulatorParameters *sp, GeomData *gcd, MeshData *md){
		matvec_struct = new Data_struct;
		pPPData = ppd;
		pGCData = gcd;
		pSimPar = sp;
		pMData = md;
		theMesh = mesh;
		one_eighth = (double)(1.0/8.0);
		firstVTK = true;

		// Auxiliary vector to assembly distributed matrix system of equations
		pMData->rowsToImport(mesh,matvec_struct->nrows,matvec_struct->rows);
	}
	
	EBFV1_elliptic::~EBFV1_elliptic(){
		if (matvec_struct){
			delete[] matvec_struct->rows;
			matvec_struct->rows = 0;
			delete matvec_struct;
			matvec_struct = 0;
		}
	}
	
	// solves system of equation for pressure field
	double EBFV1_elliptic::solver(pMesh theMesh){
		if (pSimPar->userRequiresAdaptation()){
			matvec_struct = new Data_struct;
		}
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tIN\n";
		#endif

		assembly_EFG_RHS(theMesh);
		setMatrixFreeOperation(theMesh);
		updatePressure(theMesh);
		calculatePressureGradient();
		freeMemory();

		// print VTK only for: Sw = S_initial and p = p(Sw_initial). VTK zero
		if (firstVTK){
			if (!pSimPar->timeToPrintVTK()){
				pSimPar->setTimetoPrintVTK();	// forces VTK to be printed
				pSimPar->printOutVTK(theMesh,pPPData,0,pSimPar,pGCData,exportSolutionToVTK);
				firstVTK = false;
			}
		}

		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tOUT\n";
		#endif
		return 0;
	}
	
	double EBFV1_elliptic::updatePressure(pMesh theMesh){
		CPU_Profile::Start();
		
		PetscScalar *sol, val;
		PetscInt i,m,n,row,col=0;
		PetscInt numGN = pMData->getNum_GNodes();
		Mat mSol,mLSol;
		
		// create a column matrix to receive output vector values (mSol)
		MatCreateAIJ(PETSC_COMM_WORLD,numGN,1,numGN,1,0,PETSC_NULL,0,PETSC_NULL,&mSol);
		MatSetUp(mSol);
		MatSetOption(mSol, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		
		int nLIDs, *IDs_ptr;
		pMData->getRemoteIDs(nLIDs,&IDs_ptr);
		
		// transference process: from vector to column matrix
		VecGetArray(output,&sol);

		MatSetValues(mSol,matvec_struct->nrows,matvec_struct->rows,1,&col,sol,INSERT_VALUES);
		VecRestoreArray(output,&sol);
		MatAssemblyBegin(mSol,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(mSol,MAT_FINAL_ASSEMBLY);
		//VecView(output,PETSC_VIEWER_STDOUT_WORLD); //
		
		// remote values cannot be gotten from remote matrix positions. transfer (via MatGetSubMatrixRaw) necessary remote values 
		// to each process to a second column matrix (mLSol).
		IS rows, cols;
		ISCreateGeneral(PETSC_COMM_WORLD,nLIDs,IDs_ptr,PETSC_COPY_VALUES,&rows);
		ISCreateGeneral(PETSC_COMM_WORLD,1,&col,PETSC_COPY_VALUES,&cols);
		MatGetSubMatrix(mSol,rows,cols,MAT_INITIAL_MATRIX,&mLSol);
		MatDestroy(&mSol);
		MatGetOwnershipRange(mLSol,&m,&n);
		
		// loop over IDs_ptr[i]
		row = m;
		for(i=0; i<nLIDs;i++){
			int ID = pMData->get_PETScToApp_Ordering(IDs_ptr[i]+1);
			MatGetValues(mLSol,1,&row,1,&col,&val);
			pPPData->setPressure(ID-1,val);
			row++;
		}
		MatDestroy(&mLSol);
		
		static bool key = true;
		if (key){
			int nnodes = M_numVertices(theMesh);
			for(i=0;i<nnodes;i++){
				int ID = pMData->get_AppToPETSc_Ordering(i+1);
				if ( pMData->getDirichletValue(ID,&val) ){
					pPPData->setPressure(i,val);
				}
			}
			key = false;
		}

		CPU_Profile::End("updatePressure");
		return 0;
	}
	
	double EBFV1_elliptic::freeMemory(){
		CPU_Profile::Start();

		if (!pSimPar->useDefectCorrection()){
			MatDestroy(&matrix);
		}
		MatDestroy(&matvec_struct->G);
		for(int i=0; i<pSimPar->getNumDomains(); i++){
			MatDestroy(&matvec_struct->E[i]);
			MatDestroy(&matvec_struct->F[i]);
		}
		
		/// free vectors
		VecDestroy(&matvec_struct->RHS);
		VecDestroy(&matvec_struct->z);
		VecDestroy(&output);
		delete[] matvec_struct->F;
		delete[] matvec_struct->E;
		matvec_struct->F = 0;
		matvec_struct->E = 0;
		if (pSimPar->userRequiresAdaptation()){
			delete[] matvec_struct->rows; matvec_struct->rows = 0;
			delete matvec_struct; matvec_struct = 0;
		}

		CPU_Profile::End("freeMemory");
		return 0;
	}
}
