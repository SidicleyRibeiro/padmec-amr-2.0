/*
 * EBFV1_AssemblyMatVec.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: rogerio
 */


#include "EBFV/EBFV1_elliptic.h"

namespace PRS{

	/*
	Based on a vertex centered and edge based data structured finite volume formulation, elliptic equation is calculated through a linear
	system of equation where the stiffness matrix depends on three other matrices:

	Ax=b

	where,

	A = E*F + G
	For more details about physical/numerical meaning of all these matrices see: http://www.repositorio.ufpe.br/handle/123456789/5336
	*/

	double EBFV1_elliptic::assembly_EFG_RHS(pMesh mesh){
		CPU_Profile::Start();
		int np = pMData->getNum_GNodes();
		int numGF = pMData->getNum_GF_Nodes();
		int dim = pGCData->getMeshDim();
		int ndom = pSimPar->getNumDomains();

		if (!dim || !np || !numGF || !ndom){
			char msg[256]; sprintf(msg,"dim = %d, np = %d, numGF = %d, ndom = %d.  NULL parameters found!\n",dim,np,numGF,ndom);
			throw Exception(__LINE__,__FILE__,msg);
		}
		matvec_struct->ndom = ndom;
		matvec_struct->F_nrows = np*dim;
		matvec_struct->F_ncols = np;

		// these matrices take all mesh vertices and will be delete very soon.
		Mat G_tmp, E[ndom], F[ndom];

		// Create matrix G.
		double Cij[3];
		int i,j,nedges, dom, idx0, idx1,idx0_global, idx1_global, id0, id1, dom_flag;
		MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,np,np,80,PETSC_NULL,80,PETSC_NULL,&G_tmp);
		for (dom=0; dom<ndom; dom++){
			MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,np,np*dim,100,PETSC_NULL,100,PETSC_NULL,&E[dom]);
			MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,np*dim,np,100,PETSC_NULL,100,PETSC_NULL,&F[dom]);
			nedges = pGCData->getNumEdgesPerDomain(dom);
			dom_flag = pGCData->getDomFlag(dom);
			for (i=0; i<nedges; i++){
				pGCData->getCij(dom,i,Cij);
				pGCData->getEdge(dom,i,idx0,idx1,idx0_global,idx1_global);
				pGCData->getID(dom,idx0,idx1,id0,id1);
				divergence_E(E[dom],Cij,i,dom,dom_flag,idx0_global,idx1_global,id0,id1,dim);
				divergence_G(G_tmp,Cij,i,dom,dom_flag,idx0_global,idx1_global,id0,id1,dim);
				gradient_F_edges(F[dom],Cij,dom,idx0,idx1,id0,id1,dim);
			}
			gradient_F_bdry(F[dom],dom);
			assemblyMatrix(F[dom]);
			assemblyMatrix(E[dom]);
		}
		assemblyMatrix(G_tmp);

		// Get from matrix G its contribution to RHS. matvec_struct->G correspond to all free nodes
		set_SOE(mesh,G_tmp,matvec_struct->G,true,matvec_struct->RHS,true,true);

		// Get from matrices E and F their contribution to RHS. E and F must be handled domain by domain.
		Mat EF, tmp;
		Vec EF_rhs;
		for (i=0; i<ndom; i++){
			MatMatMult(E[i],F[i],MAT_INITIAL_MATRIX,1.0,&EF);
			if (pSimPar->useDefectCorrection()){
				//set_SOE(mesh,EF,EF_multiDom[i],true,EF_rhs,true,false);
			}
			else{
				set_SOE(mesh,EF,tmp,false,EF_rhs,true,false);		// EF matrix is destroyed inside set_SOE function
			}
			VecAXPY(matvec_struct->RHS,1.0,EF_rhs);
			VecDestroy(&EF_rhs);
		}

		// Create E and F matrices related to free nodes only. Note that they were been created using all mesh nodes.
		/* TRANSFER VALUES FROM TEMPORARY E MATRIX TO matvec_struct->E */
		matvec_struct->F = new Mat[ndom];
		matvec_struct->E = new Mat[ndom];

		// Gets all E[i] columns. Same for all processors
		int nrows = matvec_struct->nrows;
		int *rows = matvec_struct->rows;

		static bool cvfm = true;
		if (!pSimPar->userRequiresAdaptation()){
			if (cvfm){
				pMData->createVectorsForMatrixF(F[0]);
				cvfm = false;
			}
		}
		else{
			pMData->createVectorsForMatrixF(F[0]);
		}

		IS rows_F, cols_F, rows_E, cols_E;
		ISCreateGeneral(PETSC_COMM_WORLD,pMData->get_F_nrows(),pMData->get_F_rows_ptr(),PETSC_COPY_VALUES,&rows_F);
		ISCreateGeneral(PETSC_COMM_WORLD,numGF,pMData->get_F_cols_ptr(),PETSC_COPY_VALUES,&cols_F);
		ISCreateGeneral(PETSC_COMM_WORLD,nrows,rows,PETSC_COPY_VALUES,&rows_E);
		ISCreateGeneral(PETSC_COMM_WORLD,np*dim,pMData->get_pos_ptr(),PETSC_COPY_VALUES,&cols_E);

		for (i=0; i<ndom; i++){
			MatGetSubMatrix(F[i],rows_F,cols_F,MAT_INITIAL_MATRIX,&matvec_struct->F[i]);
			MatGetSubMatrix(E[i],rows_E,cols_E,MAT_INITIAL_MATRIX,&matvec_struct->E[i]);
			MatDestroy(&E[i]);
			MatDestroy(&F[i]);
		}

		ISDestroy(&rows_F);
		ISDestroy(&cols_F);
		ISDestroy(&rows_E);
		ISDestroy(&cols_E);

		// Create the output vector
		VecCreate(PETSC_COMM_WORLD,&output);
		VecSetSizes(output,PETSC_DECIDE,numGF);
		VecSetFromOptions(output);
		CPU_Profile::End("MatricesAssembly");
		return 0;
	}

	int EBFV1_elliptic::assemblyMatrix(Mat A){
		MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
		return 0;
	}

	int EBFV1_elliptic::set_SOE(pMesh mesh, Mat A, Mat &LHS, bool assemblyLHS, Vec &RHS,bool assemblyRHS, bool includeWell){
		int i=0, j=0;
		int numGF = pMData->getNum_GF_Nodes();	// set global free nodes
		int numGP = pMData->getNum_GP_Nodes();	// set global prescribed (dirichlet) nodes

		// -------------------------------------------------------------------------
		// 					ASSEMBLY MATRIX LHS (STIFFNESS MATRIX)
		// -------------------------------------------------------------------------
		int nrows = matvec_struct->nrows;
		int *rows = matvec_struct->rows;

		IS rowsToImport;
		ISCreateGeneral(PETSC_COMM_WORLD,nrows,rows,PETSC_COPY_VALUES,&rowsToImport);

		IS colsToImport;
		ISCreateGeneral(PETSC_COMM_WORLD,numGF,pMData->get_idxFreecols_ptr(),PETSC_COPY_VALUES,&colsToImport);

#ifdef _SEEKFORBUGS_
		if (!nrows){
			throw Exception(__LINE__,__FILE__,"nrows NULL!");
		}
		if (!rows || !pMData->get_idxFreecols_ptr()){
			throw Exception(__LINE__,__FILE__,"rows NULL!");
		}
#endif

		if (assemblyLHS){
			MatGetSubMatrix(A,rowsToImport,colsToImport,MAT_INITIAL_MATRIX,&LHS);
		}

		// --------------------------- ASSEMBLY RHS VECTOR -------------------------
		// matrix formed by columns associated to prescribed nodes.
		// -------------------------------------------------------------------------
		if (assemblyRHS){
			Mat rhs;
			int m,n;

			IS colsToImport_dirichlet;
			ISCreateGeneral(PETSC_COMM_WORLD,numGP,pMData->get_idxn_ptr(),PETSC_COPY_VALUES,&colsToImport_dirichlet);
			MatGetSubMatrix(A,rowsToImport,colsToImport_dirichlet,MAT_INITIAL_MATRIX,&rhs);
			MatDestroy(&A);
			ISDestroy(&colsToImport_dirichlet);

			VecCreate(PETSC_COMM_WORLD,&RHS);
			VecSetSizes(RHS,PETSC_DECIDE,numGF);
			VecSetFromOptions(RHS);

			// fill RHS vector with values from rhs matrix. entries from this matrix are multiplied by their corresponding prescribed values
			MatGetOwnershipRange(rhs,&m,&n);
			for (i=m; i<n; i++){
				double sum = .0;
				MIter prescribedIter = pMData->dirichletBegin();
				for (j=0; j<numGP; j++, prescribedIter++){
					int col = j;
					double val=.0;
					MatGetValues(rhs,1,&i,1,&col,&val);
					//printf("sum_old = %f  ",sum);
					sum += -val*prescribedIter->second;
					//printf("val: %f  prescVal: %f     sum = %f\n",val,prescribedIter->second,sum);
				}
				VecSetValue(RHS,i,sum,ADD_VALUES);
			}

			if ( includeWell ){
				wellsContributionToRHS(mesh,RHS);
			}
			VecAssemblyBegin(RHS);
			VecAssemblyEnd(RHS);
			MatDestroy(&rhs);
		}
		rows = 0;
		ISDestroy(&rowsToImport);
		ISDestroy(&colsToImport);
		return 0;
	}


	// NOTE: mapNodesOnWells should not be available in the way it appears here. Some function should be implemented inside SimulatorParameters
	// to provide only the data required.
	int EBFV1_elliptic::wellsContributionToRHS(pMesh mesh, Vec &RHS){
		int node_ID, row;
		double Vt, Vi, Qi, Qt;

		if (!pSimPar->mapNodesOnWell.size()){
			throw Exception(__LINE__,__FILE__,"No wells found!");
		}

		// for each well flag
		map<int,set<int> >::iterator mit = pSimPar->mapNodesOnWell.begin();
		for (; mit!=pSimPar->mapNodesOnWell.end(); mit++){
			int well_flag = mit->first;
			// source/sink term
			Qt = pSimPar->getFlowrateValue(well_flag);
			if ( fabs(Qt)<=1e-7 ){
				throw Exception(__LINE__,__FILE__,"Flow rate NULL!");
			}

			//cout << "well-flag: " << well_flag << endl;

			// get all flagged node IDs for that well
			if (!mit->second.size()){
				throw Exception(__LINE__,__FILE__,"No wells found!");
			}
			SIter sit = mit->second.begin();
			for (; sit!=mit->second.end(); sit++){
				Vt = pSimPar->getWellVolume(well_flag);
				#ifdef _SEEKFORBUGS_
					if ( Vt<1e-12 ){
						char msg[256]; sprintf(msg,"Well with null volume V = %.6f. Vertex (%d)",Vt,node_ID);
						throw Exception(__LINE__,__FILE__,msg);
					}
				#endif //_SEEKFORBUGS_

				Vi = .0;
				node_ID = *sit;
				for (SIter_const dom=pSimPar->setDomain_begin(); dom!=pSimPar->setDomain_end(); dom++){
					pVertex node = (mEntity*)mesh->getVertex( node_ID );
					Vi += pGCData->getVolume(node,*dom);
				}

				// for node i, Q is a fraction of total well flow rate
				Qi = Qt*(Vi/Vt);

				// FPArray ('F'ree 'P'rescribed array) maps node id: node_ID -> row
				// row: position in Petsc GlobalMatrix/RHSVBector where node must be assembled
				node_ID = pMData->get_AppToPETSc_Ordering(node_ID);
				row = pMData->FPArray(node_ID-1);

				/*
				 * Do not include well flux on nodes with prescribed pressure
				 */
				if (pSimPar->isNodeFree(well_flag)){
// 										cout << "---------------------------------------------\n";
// 										cout << setprecision(8);
// 										cout << "Node = " << node_ID << endl;
// 										cout << "Qt = " << Qt << endl;
// 										cout << "Qi = " << Qi << endl;
// 										cout << "Vt = " << Vt << endl;
// 										cout << "Vi = " << Vi << endl;
// 									//	cout << "row = " << row << endl;
// 										cout << "---------------------------------------------\n";
					VecSetValue(RHS,row,Qi,ADD_VALUES);

				}
				//printf("well contribution to RHS: node %d -> row %d\n",node_ID,row);
				//printf("[%d] - Qi = %f on row: %d \n",P_pid(),Qi,row);
				//
				////				throw Exception(__LINE__,__FILE__,"Favor consertar gambiarra!");
				//				}
			}
		}
		/*
		 * This is a way to use the simulator to evaluate elliptic equation without screw-up
		 * the input data procedure.
		 * Treating source/sink terms:
		 */
	#ifdef CRUMPTON_EXAMPLE
		double x,y,coord[3];
		double srcsnk, vol1, vol2, f_xy1, f_xy2;

		int nrows = matvec_struct->nrows;	// number of free nodes
		int *rows = matvec_struct->rows;	// free nodes indices
		/*
		 *  1 - Take only free nodes
		 *  2 - This a specific problem with only two sub-domains flagged as 3300 and 3301.
		 *  3 - A node located on boundary domains has two control-volumes, each one associated to a sub-domain
		 *      Source/sink term must be applied to both
		 *  4 - When a node's volume is entirely located in a sub-domain, vol1 or vol2 will be equal to zero.
		 */
		for (int i=0; i<nrows; i++){
			pVertex vertex = mesh->getVertex(rows[i]+1);
			if (vertex){
				V_coord(vertex,coord); x = coord[0]; y = coord[1];
				vol1 = pGCData->getVolume(vertex,3300); 			// it supposes coord_x <= 0
				vol2 = pGCData->getVolume(vertex,3301); 			// it supposes coord_x  > 0
				f_xy1 = -(2.*sin(y) + cos(y))*ALPHA*x - sin(y);		// source/sink for coord_x <= 0
				f_xy2 = 2.*ALPHA*exp(x)*cos(y);						// source/sink for coord_x  > 0
				srcsnk = vol1*f_xy1 + vol2*f_xy2;					// source/sink term
				node_ID = pMData->get_AppToPETSc_Ordering(rows[i]+1);
				row = pMData->FPArray(node_ID-1);
				VecSetValue(RHS,row,-srcsnk,ADD_VALUES);
			}
		}
	#endif
		return 0;
	}
}



