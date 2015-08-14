/*
 * setCorrectNumberOfRemoteCopies.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: rogsoares
 *
 *  For a distributed mesh representing a multidomain media, some edges with
 *  remote copies are located over domains' boundary. The number of remote copies
 *  of an edge should be recalculated if its remote copies belong to other domains.
 *  For example, if an edge between two domains has a remote copy due to mesh
 *  partition, it must not be treated as it had a remote copy. It will mass matrix
 *  assembly. Then, setCorrectNumberOfRemoteCopies function was built to avoid this
 *  kind of error.
 */
#include "EBFV1__pre-processors.h"

ofstream fid;

PetscErrorCode ierr;
typedef std::list<pEntity> EList;
typedef std::set<int> setint;
typedef setint::iterator setintIter;
typedef std::map<int,setint> mapset;
typedef mapset::iterator mapsetIter;

int excludeCopiesOtherDomains(PP_Parameters*, int, EList &);
int setNewNumberOfRemoteCopies(Mat&, PP_Parameters*, EList &, int);
int getsubMatrix(Mat&, PP_Parameters*, const EList &, Mat&);
int createExchangeMatrix(PP_Parameters*, const EList &, Mat&);
void getEdgesWithRemoteCopies(PP_Parameters*, int, EList &);
void STOP();

int setCorrectNumberOfRemoteCopies(PP_Parameters* pPPP){
    #ifdef _PREPROCESSOR_DEBUG_
		printf(" setCorrectNumberOfRemoteCopies....\n");
	#endif
	if (pPPP->setOfDomain.size()==1 && P_size()==1) return 0;
	// loop over all domains
	EList edgeList;
	for (set<int>::iterator SIter = pPPP->setOfDomain.begin(); SIter != pPPP->setOfDomain.end(); SIter++){
		int dom = *SIter;		// domain flag

		// Loop over all edges belonging to domain dom. Take all edges with remote copies
		getEdgesWithRemoteCopies(pPPP,dom,edgeList);
		excludeCopiesOtherDomains(pPPP,dom,edgeList);
	}
	return 0;
}

int excludeCopiesOtherDomains(PP_Parameters* pPPP, int dom, EList &edgeList){
	/*
	 * That's the plan:
	 *
	 *  Let's create a a distributed matrix n by n, where n is the largest node global ID.
	 *  Each rank will put for each remote edge its contribution. At the end, will know
	 *  how many remote copies an edge will have for domain 'dom'.
	 *  Though, an initially edge with remote copy will not have remote copies any more, or the number
	 *  could be less than before.
	 */

	// loop over edgeList and create exchange matrix
	Mat exchangeMatrix;
	createExchangeMatrix(pPPP,edgeList,exchangeMatrix);

	// each partition take a sub-matrix from exchangeMatrix
	// the number of rows of the sub-matrix is the same of the number of partition edges
	// with remote copies which these numbers could be different from now on.
	Mat subMatrix;
	getsubMatrix(exchangeMatrix,pPPP,edgeList,subMatrix);
	setNewNumberOfRemoteCopies(subMatrix,pPPP,edgeList,dom);

	// Clean memory for next domain
	ierr = MatDestroy(exchangeMatrix); CHKERRQ(ierr);
	ierr = MatDestroy(subMatrix); CHKERRQ(ierr);
	edgeList.clear();
	return 0;
}

int setNewNumberOfRemoteCopies(Mat &subMatrix, PP_Parameters* pPPP, EList &edgeList, int dom){
	mapset mapRows_SetCols;
	for (EList::const_iterator iter = edgeList.begin(); iter != edgeList.end(); iter++){
		int id0 = EN_id( (*iter)->get(0,0) ); // ID_0
		int id1 = EN_id( (*iter)->get(0,1) ); // ID_1
		if (id1 < id0) std::swap(id0,id1);

		mapsetIter msIter = mapRows_SetCols.find(id0);
		if (msIter != mapRows_SetCols.end())
			msIter->second.insert(id1);
		else{
			setint setCols;
			setCols.insert(id1);
			mapRows_SetCols[id0] = setCols;
		}
	}

	double nrc;
	//bool hasRC = false;
	PetscInt m, n, row, col;
	ierr = MatGetOwnershipRange(subMatrix,&m,&n);CHKERRQ(ierr);
	row = m;
	// loop over "rows"
	for (mapsetIter msIter = mapRows_SetCols.begin(); msIter != mapRows_SetCols.end(); msIter++){
		// loop over "columns"
		int id0 = msIter->first;
		for (setintIter sIter = msIter->second.begin(); sIter != msIter->second.end(); sIter++){
			int id1 = *sIter;
			col = id1;
			ierr = MatGetValues(subMatrix,1,&row,1,&col,&nrc);CHKERRQ(ierr);

			mVertex* v0 = pPPP->theMesh->getVertex(id0);
			mVertex* v1 = pPPP->theMesh->getVertex(id1);
			mEntity* edge = (mEntity*)pPPP->theMesh->getEdge(v0,v1);

			// new number of remote copies = nrc - 1;
			// if nrc=1, it means that for domain dom there are one one edge, then it must not
			// be treated having remote copies. In this case, nrc would be zero.
			nrc = nrc - 1.0;
			pPPP->pGCData->setNumRC(edge,dom,nrc);
			#ifdef _PREPROCESSOR_DEBUG_
				printf("edge[%d %d] nrc_old: %d  nrc_new: %f for domain %d\n",
					id0,id1,pPPP->pGCData->getNumRemoteCopies(pPPP->theMesh,edge,pPPP->flg),nrc-1.,dom);
			#endif
		}
		row++;
	}
	return 0;
}

int getsubMatrix(Mat &exchangeMatrix, PP_Parameters* pPPP, const EList &edgeList, Mat &subMatrix){
	if (!exchangeMatrix){
		printf("NULL matrix. Exiting....\n"); throw 1;
	}
	int numNG = M_getMaxVId(pPPP->theMesh) + 1;
	std::set<int> setRows;
	for (EList::const_iterator iter = edgeList.begin(); iter != edgeList.end(); iter++){
		int id0 = EN_id( (*iter)->get(0,0) ); // ID_0
		int id1 = EN_id( (*iter)->get(0,1) ); // ID_1
		if (id1 < id0) std::swap(id0,id1);
		setRows.insert(id0);
	}

	PetscInt nrows = (int)setRows.size();
	PetscInt ncols = numNG;
	int rows[nrows];
	int cols[ncols];
	int i = 0;
	for (std::set<int>::iterator iter=setRows.begin(); iter!=setRows.end(); iter++) rows[i++] = *iter;
	for (i=0; i<numNG; i++) cols[i] = i;
	setRows.clear();

	MatGetSubMatrixRaw(exchangeMatrix,nrows,rows,ncols,cols,PETSC_DECIDE,MAT_INITIAL_MATRIX,&subMatrix); CHKERRQ(ierr);
	return 0;
}

int createExchangeMatrix(PP_Parameters* pPPP, const EList &edgeList, Mat &exchangeMatrix){
	//printf("createExchangeMatrix\n");
	int numNG = M_getMaxVId(pPPP->theMesh) + 1;
	ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,numNG,numNG,0,PETSC_NULL,0,PETSC_NULL, &exchangeMatrix); CHKERRQ(ierr);
	for (EList::const_iterator iter = edgeList.begin(); iter != edgeList.end(); iter++){
		int id0 = EN_id( (*iter)->get(0,0) ); // ID_0
		int id1 = EN_id( (*iter)->get(0,1) ); // ID_1
		if (id1 < id0) std::swap(id0,id1);
	//	fid << id0 << " " << id1 << endl;
		ierr = MatSetValue(exchangeMatrix,id0,id1,1.0,ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(exchangeMatrix,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(exchangeMatrix,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	return 0;
}

void getEdgesWithRemoteCopies(PP_Parameters* pPPP, int dom, EList &edgeList){
	pEntity edge;
	EIter eit = M_edgeIter(pPPP->theMesh);
	while ( (edge  = EIter_next(eit)) ){
		if (pPPP->pGCData->edgeBelongToDomain(edge,dom) && M_numRemoteCopies(pPPP->theMesh,edge) > 0 )
			edgeList.push_back(edge);
	}
	EIter_delete(eit);
}

//void STOP(){
//	MPI_Barrier(MPI_COMM_WORLD); throw 1;
//}

