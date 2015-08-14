//#include "EBFV1__pre-processors.h"
//
//typedef std::set<int> IntSet;
//
//int unifyVolumesAmongProcessors(pMesh theMesh, const set<int> &setOfDomains, const string &whatUnify, GeomData *pGCData){
//	if (P_size() == 1) return 0;
//
//	PetscErrorCode ierr;
//
//	PetscTruth flg;
//	PetscOptionsHasName(PETSC_NULL,"-refine",&flg);
//
//
//	/*
//	 * unifiedCV_matrix: sparse distributed matrix where each row and column cor-
//	 * responds, respectively, to a node ID and a domain. It unifies Control vo-
//	 * lumes of nodes located on the partition boundaries for all domains.
//	 */
//	Mat unifiedCV_matrix;
//	const int nrows = M_getMaxVId(theMesh);
//	const int ncols = (int)setOfDomains.size();
//
//	// parallel matrix to unify nodalCV from nodes with remote copies
//	ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,nrows,ncols,0,PETSC_NULL,0,PETSC_NULL,&unifiedCV_matrix);CHKERRQ(ierr);
//
//	/*
//	 * First part:
//	 */
//	int i, row, col = 0;
//	pEntity node;
//
//	IntSet remoteNodeIds;
//	IntSet::iterator SIter;
//	std::list<IntSet> domainList;
//	//printf("rank %d has %d domains ncols: %d\n",P_pid(),(int)setOfDomains.size(),ncols);
//	for (SIter = setOfDomains.begin(); SIter != setOfDomains.end(); SIter++){
//		const int dom = *SIter;	// domain's flag
//
//		// loop over nodes to get those with remote copies
//		VIter vit = M_vertexIter(theMesh);
//		while ( (node = VIter_next(vit)) ){
//			if ( M_numRemoteCopies(theMesh,node) && pGCData->nodeBelongToDomain(node,dom) ){
//				row = EN_id(node) - 1;						// get node ID
//				remoteNodeIds.insert(row);				// store node ID
//				double val = pGCData->getVolume(node,dom);	// get volume associated to domain dom
//				ierr = MatSetValue(unifiedCV_matrix,row,col,val,ADD_VALUES); CHKERRQ(ierr);
//			}
//		}
//		VIter_delete(vit);
//
//		if (!(int)remoteNodeIds.size()){
//			remoteNodeIds.insert(0);
//			//ierr = MatSetValue(unifiedCV_matrix,0,col,.0,ADD_VALUES); CHKERRQ(ierr);
//		}
//
//		domainList.push_back(remoteNodeIds);
//		//printf("%d %d numRemoteNodes: %d\n",P_pid(),dom,(int)remoteNodeIds.size());
//		remoteNodeIds.clear();
//		col++;												// next domain
//	}
//	ierr = MatAssemblyBegin(unifiedCV_matrix,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//	ierr = MatAssemblyEnd(unifiedCV_matrix,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//	//ierr = MatView(unifiedCV_matrix,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); //throw 1;
//
//	/*
//	 * Second part:
//	 */
//	col = 0; // each column -> domain
//	SIter = setOfDomains.begin();
//	std::list<IntSet>::iterator LIter = domainList.begin();
//	for (; LIter != domainList.end(); LIter++){
//		const int dom = *SIter;	// domain's flag
//		remoteNodeIds = *LIter; // node list (each set corresponds to a domain)
//
//		// transfer all remote nodes from a set container to an array
//		int numRemoteNodes = ( (int)remoteNodeIds.size() )?( (int)remoteNodeIds.size() ):1;
//		int *pNodeIDs = new int[numRemoteNodes];
//
//		if (!(int)remoteNodeIds.size()) pNodeIDs[0] = 0;
//
//		i = 0;
//		IntSet::iterator RNIDs_Iter = remoteNodeIds.begin();
//		for (; RNIDs_Iter != remoteNodeIds.end(); RNIDs_Iter++) pNodeIDs[i++] = *RNIDs_Iter;
//		remoteNodeIds.clear();
//
//		// indices used for petsc function
//		IS isrow, iscol;
//		ierr = ISCreateGeneral(PETSC_COMM_WORLD,numRemoteNodes,pNodeIDs,&isrow); CHKERRQ(ierr);
//		ierr = ISCreateGeneral(PETSC_COMM_WORLD,1,&col,&iscol); CHKERRQ(ierr);
//
//		PetscInt size;
//		ISGetSize(isrow,&size);
//		//printf("%d %d numRemoteNodes: %d  ISROWsize: %d\n",P_pid(),col,numRemoteNodes,size);
//
//		/*
//		 * Unified CV must have local access. We use MatGetSubMatrix to gather all
//		 * nodes on partition boundary required by each process.
//		 */
//		Mat newUnifiedCV_matrix; // matrix with all requested nodes
//		ierr = MatGetSubMatrix(unifiedCV_matrix, isrow, iscol, PETSC_DECIDE, MAT_INITIAL_MATRIX,&newUnifiedCV_matrix);CHKERRQ(ierr);
//		ierr = ISDestroy(isrow); CHKERRQ(ierr);
//		ierr = ISDestroy(iscol); CHKERRQ(ierr);
//
//		PetscInt m, n;
//		ierr = MatGetOwnershipRange(newUnifiedCV_matrix,&m,&n); CHKERRQ(ierr);
//		const int numLocalRows = n - m; 			// number of rows in newUnifiedCV_matrix owned to each process
//		double *nodalCV = new double[numLocalRows];
//		int *rowIndices = new int[numLocalRows];
//		for (i=0; i<numLocalRows; i++) rowIndices[i] = m + i;
//
//		/*
//		 * Get from newUnifiedCV_matrix unified nodal CV for each domain
//		 */
//		i = 0;
//		int ncol = 0;
//		ierr = MatGetValues(newUnifiedCV_matrix,numLocalRows,rowIndices,1,&ncol,nodalCV);CHKERRQ(ierr);
//		VIter vit = M_vertexIter(theMesh);
//		while ( (node = VIter_next(vit)) ){
//			// node must be located on partition boundary and belong to domain dom
//			if ( M_numRemoteCopies(theMesh,node) && pGCData->nodeBelongToDomain(node,dom) )
//				pGCData->setVolume(node,dom,nodalCV[i++]);
//		}
//		VIter_delete(vit);
//
//		ierr = MatDestroy(newUnifiedCV_matrix); CHKERRQ(ierr);
//		delete[] rowIndices; rowIndices = 0;
//		delete[] nodalCV; nodalCV = 0;
//		delete[] pNodeIDs; pNodeIDs = 0;
//		SIter++;
//		col++;
//	}
//	ierr = MatDestroy(unifiedCV_matrix); CHKERRQ(ierr);
//	return 0;
//}
