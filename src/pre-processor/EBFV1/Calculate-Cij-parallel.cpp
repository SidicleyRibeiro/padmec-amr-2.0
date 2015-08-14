///**************************************************************************************
//            UNIVERSIDADE FEDERAL DE PERNAMBUCO
//Author:     Rogerio Soares da Silva
//Date:       2006-2011
//Licence:    This program is free software; you can redistribute it and/or modify it.
//Warranty:   None! Use it at your own risk!
//Visit:      www.padmec.org
//**************************************************************************************/
//
//
//#include "EBFV1__pre-processors.h"
//
//
//int unifyCijAmongProcessors(EdgeInfo* edgeinfo, pMesh theMesh);
///* ----------------------------------------------------------------------------------------
//	In parallel Cij is calculated partially on each partition, but they need the same value
//	as they were calculated in serial. Cij(rank p, dom k) = Cij(rank q, dom k)
//	Cij(rank p, dom k) = Cij(rank p, dom k) + Cij(rank q, dom k)
//	Cij(rank q, dom k) = Cij(rank q, dom k) + Cij(rank p, dom k)
//	Only edges on partion boundary are treated here. Do nothing in serial.
//------------------------------------------------------------------------------------------- */
//
//
///*
// * Collect all infomation about edges with remote copies to be unified
// */
//int unifyCijAmongProcessors(pMesh theMesh, const set<int> &setOfDomains, GeomData *pGCData){
//	if (P_size()==1) return 1;
//
//	PetscTruth flg;
//	PetscOptionsHasName(PETSC_NULL,"-refine",&flg);
//	EdgeInfo* edgeinfo;
//	pEntity edge;
//	vector<double> Cij(3);
//
//	set<int>::iterator SIter;
//	for (SIter = setOfDomains.begin(); SIter != setOfDomains.end(); SIter++){
//		const int dom = *SIter;		// domain flag
//		edgeinfo = new EdgeInfo;	// used to store edge node's IDs and Cij coordenates
//
//		// loop over all edges: get only those with remote copies. get edge node's IDs and Cij
//		EIter eit = M_edgeIter(theMesh);
//		while ( (edge  = EIter_next(eit)) ){
//			if (pGCData->edgeBelongToDomain(edge,dom) && M_numRemoteCopies(theMesh,edge) > 0 ){
//				int id0 = EN_id(edge->get(0,0))-1;
//				int id1 = EN_id(edge->get(0,1))-1;
//				if (id1 < id0) swap(id0,id1);
//
//				// GET IDS
//				IdPair idpair(id0,id1);
//				edgeinfo->setEdgeIds(idpair);
//
//				// GET Cij
//				double val[3] = {.0,.0,.0};
//				pGCData->getCij(edge,dom,Cij); std::copy(Cij.begin(),Cij.end(),val);
//				edgeinfo->setRemoteEdgeCij(edgeinfo->xlist,val[0]);
//				edgeinfo->setRemoteEdgeCij(edgeinfo->ylist,val[1]);
//				edgeinfo->setRemoteEdgeCij(edgeinfo->zlist,val[2]);
//			}
//		}
//		edgeinfo->mapRows();
//
//		// UNIFY CIJ FOR EDGES WITH REMOTE COPIES
//		// ---------------------------------------------------------------------
//		// unify Cij_x
//
//		edgeinfo->setWhichList(edgeinfo->xlist);
//		unifyCijAmongProcessors(edgeinfo,theMesh);
//
//		// unify Cij_y
//		edgeinfo->setWhichList(edgeinfo->ylist);
//		unifyCijAmongProcessors(edgeinfo,theMesh);
//
//		// unify Cij_z
//		edgeinfo->setWhichList(edgeinfo->zlist);
//		unifyCijAmongProcessors(edgeinfo,theMesh);
//
//		// update Cij on mesh data structure
//		list<double>::iterator Cij_xIter = edgeinfo->Cij_beginIterator(edgeinfo->xlist);
//		list<double>::iterator Cij_yIter = edgeinfo->Cij_beginIterator(edgeinfo->ylist);
//		list<double>::iterator Cij_zIter = edgeinfo->Cij_beginIterator(edgeinfo->zlist);
//
//		// loop over all edges once more to update only those with remote copies
//		eit = M_edgeIter(theMesh);
//		while ( (edge  = EIter_next(eit)) ){
//			// only edge with remote copies
//			if (pGCData->edgeBelongToDomain(edge,dom) && M_numRemoteCopies(theMesh,edge) > 0){
//				Cij[0] = *Cij_xIter; Cij_xIter++;
//				Cij[1] = *Cij_yIter; Cij_yIter++;
//				Cij[2] = *Cij_zIter; Cij_zIter++;
//				pGCData->setCij(edge,dom,Cij);
//			}
//		}
//		EIter_delete(eit);
//
//		// free memory
//		edgeinfo->deleteEdgeIds();
//		edgeinfo->deleteMapRows();
//		edgeinfo->deleteCij_list(edgeinfo->xlist);
//		edgeinfo->deleteCij_list(edgeinfo->ylist);
//		edgeinfo->deleteCij_list(edgeinfo->zlist);
//		delete edgeinfo; edgeinfo = 0;
//	}
//}
//
//int unifyCijAmongProcessors(EdgeInfo* edgeinfo, pMesh theMesh){
//	PetscErrorCode ierr;
//	double val;
//	int i;
//	int numNodesGlobal = M_getMaxVId(theMesh);
//
//	// matrix to store Cij coord for all processes
//	Mat exchangeCij; // matrix to store nodes from all processes
//	ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD, PETSC_DECIDE,PETSC_DECIDE,numNodesGlobal-1,numNodesGlobal,0,PETSC_NULL, 0,PETSC_NULL, &exchangeCij); CHKERRQ(ierr);
//
//	// transfer all Cij coordenates to a parallel matrix. duplicated Cij are added
//	// maybe CijCoords can be inserted into exchangeCij without a loop.
//	const list<IdPair> idpairList = edgeinfo->getRemoteEdgeIds();
//
//	list<IdPair>::iterator Id_Iter = edgeinfo->IdPair_beginIterator();
//	list<double>::iterator Cij_Iter = edgeinfo->Cij_beginIterator();
//
//	// populate exchangeCij matrix
//	for (; Id_Iter != edgeinfo->IdPair_endIterator(); Id_Iter++, Cij_Iter++){
//		ierr = MatSetValue(exchangeCij,Id_Iter->getId0(),Id_Iter->getId1(),*Cij_Iter,ADD_VALUES); CHKERRQ(ierr);
//	}
//
//	ierr = MatAssemblyBegin(exchangeCij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//	ierr = MatAssemblyEnd(exchangeCij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//
//	// now, all nodes, regardless their location(processor), are updated.
//	Mat unifiedCij;     							// matrix with all requested nodes
//	PetscInt nrows = edgeinfo->numRowstoImport();   // how rows rank p must require from parallel matrix
//	PetscInt ncols = numNodesGlobal;
//
//	int* rows = new int[nrows];                     // rows to be imported
//	edgeinfo->rowstoImport(rows);
//
//	int* cols = new int[ncols];
//	for (i=0; i<numNodesGlobal; i++) cols[i] = i;
//
//	// coefficients will now be accessible locally
//	MatGetSubMatrixRaw(exchangeCij,nrows,rows,ncols,cols,PETSC_DECIDE,MAT_INITIAL_MATRIX,&unifiedCij); CHKERRQ(ierr);
//
//	ierr = MatDestroy(exchangeCij); CHKERRQ(ierr);
//	delete[] rows;
//	delete[] cols;
//
//	PetscInt m, n, row, col;
//	ierr = MatGetOwnershipRange(unifiedCij,&m,&n);  // range of rows from parallel matrix owned to rank p
//
//	// loop to update Cij coefficients list
//	edgeinfo->deleteCij_list();
//	for (Id_Iter = edgeinfo->IdPair_beginIterator(); Id_Iter != edgeinfo->IdPair_endIterator(); Id_Iter++, Cij_Iter++){
//		row = edgeinfo->findMappedRow(Id_Iter->getId0()) + m;
//		col = Id_Iter->getId1();
//		ierr = MatGetValues(unifiedCij,1,&row,1,&col,&val);CHKERRQ(ierr);
//		edgeinfo->setRemoteEdgeCij(val);
//	}
//	ierr = MatDestroy(unifiedCij); CHKERRQ(ierr);
//}
