/*
 * ErrorAnalysis_auxiliar.cpp
 *
 *  Created on: 29/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"

void getMaxMinData(pMesh theMesh, pEntity elem, ErrorAnalysis* pEA, int &maxDepth, int &maxRefLevel, int &minRefLevel);

void ErrorAnalysis::deletePointers(){
	if (pStore_h_new){
		delete[] pStore_h_new; pStore_h_new = 0;
	}
	else{
		throw Exception(__LINE__,__FILE__,"Nao era pra ter entrado aqui!");
	}
}

void ErrorAnalysis::resetNumFacesAroundNode(pMesh theMesh){
	pEntity node;
	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		EN_attachDataInt(node,MD_lookupMeshDataId( "numfaces" ),0);			// number of faces around a node
	}
	VIter_delete(vit);
}

void ErrorAnalysis::calculate_height_ratio(pMesh theMesh){
	double h_new_node, h_old_node;
	pEntity node;
	VIter vit = M_vertexIter(theMesh);
	while ( node = VIter_next(vit) ){
		EN_getDataDbl(node,MD_lookupMeshDataId( "h_old_node" ),&h_old_node);
		EN_getDataDbl(node,MD_lookupMeshDataId( "elem_height" ),&h_new_node);
		EN_attachDataDbl(node,MD_lookupMeshDataId( "height_ratio" ),(double)(h_new_node/h_old_node));

	}
	VIter_delete(vit);
}

void ErrorAnalysis::countNumFaceAroundNode(pMesh theMesh, bool singularity){
	int n;
	pEntity node, face;
	int nelement = 0;

	FIter fit = M_faceIter (theMesh);
	while ( (face = FIter_next(fit)) ){
		for (int i=0; i<3; i++){
			node = (pEntity)face->get(0,i);
			EN_getDataInt(node,MD_lookupMeshDataId( "numfaces" ),&n);
			EN_attachDataInt(node,MD_lookupMeshDataId( "numfaces" ),++n);
		}
	}
	FIter_delete(fit);

	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		EN_getDataInt(node,MD_lookupMeshDataId( "numfaces" ),&n);
		if (!n){
			throw Exception(__LINE__,__FILE__,"Num faces = 0!\n");
		}
	}
	VIter_delete(vit);
}

// set element height per node to zero for all mesh nodes
void ErrorAnalysis::resetElemHeightPerNode(pMesh theMesh){
	pEntity node;

	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		if (!isSingular(node)){
			EN_attachDataInt(node,MD_lookupMeshDataId( "elem_height" ),.0);		// number of faces around a node
		}
	}
	VIter_delete(vit);
}

void ErrorAnalysis::resetAllElementsAsSingular(pMesh theMesh){
	pEntity node, face;
	FIter fit = M_faceIter (theMesh);
	while ( (face = FIter_next(fit)) ){
		resetElementAsSingular(face);
	}
	FIter_delete(fit);

	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		resetElementAsSingular(node);
	}
	VIter_delete(vit);
}

void ErrorAnalysis::setElementsAsSingular(std::list<pEntity> &singularElemList){
	// the following step allows to compute element height per node for all nodes for the first scalar field
	// if singular elements were flagged before it will be difficult to identify who is who
	std::list<pEntity>::iterator iter = singularElemList.begin();
	for(;iter!=singularElemList.end();iter++){
		pEntity face = *iter;
		setElementAsSingular(face);
		for (int i=0; i<3; i++){
			pEntity node = (pEntity)face->get(0,i);
			setElementAsSingular(node);
		}
	}
	cout << "Num. singular elements: " << singularElemList.size() << endl;
	//		int c=0;
	//		vit = M_vertexIter(theMesh);						// weighting element height per node
	//		while ( node = VIter_next(vit) ){
	//			if (!isSingular(node)){
	//				c++;
	//			}
	//		}
	//		VIter_delete(vit);
	//		cout << "Num. singular vertices: " << M_numVertices(theMesh) - c << endl;
	setNumElements_Singularity((int)singularElemList.size());
	singularElemList.clear();
}

void ErrorAnalysis::store_h_new(pMesh theMesh){
	int i = 0;
	double height;
	pEntity node;
	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		if (!isSingular(node)){
			EN_getDataDbl(node,MD_lookupMeshDataId( "elem_height" ),&height);
			pStore_h_new[i] = std::min(pStore_h_new[i],height);
		}
		i++;
	}
	VIter_delete(vit);
//	cout << "M_numVertices(theMesh): " << M_numVertices(theMesh) << endl;
//	cout << "store_h_new: num. nodes not singular: " << c << endl;
}

void ErrorAnalysis::update_h_new(pMesh theMesh){
	int i = 0;
	pEntity node;
	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		//cout << EN_id(node)<< "\t" << pStore_h_new[i] << "\th_min = " << get_h_min() << "\tsingular? " << isSingular(node) << endl;
		EN_attachDataDbl(node,MD_lookupMeshDataId( "elem_height" ),pStore_h_new[i]);
		if (pStore_h_new[i]<.0){
			throw Exception(__LINE__,__FILE__,"Nengative h_new!!!");
		}
		if (pStore_h_new[i] < get_h_min()){
			char msg[256]; sprintf(msg,"h_new too small: %f\th_min = %f",pStore_h_new[i],get_h_min());
			throw Exception(__LINE__,__FILE__,msg);
		}
		if (pStore_h_new[i] > 1e+9){
			char msg[256]; sprintf(msg,"h_new too big: %f\th_min = %f",pStore_h_new[i],get_h_min());
			throw Exception(__LINE__,__FILE__,msg);
		}
		i++;
	}
	VIter_delete(vit);
}

void ErrorAnalysis::initializeParameters(pMesh theMesh){
	resetNumFacesAroundNode(theMesh);
	countNumFaceAroundNode(theMesh,false);
	calculate_CharacteristicDimensionLength(theMesh);
	this->SGN = .0;
	this->SGN_sing = .0;
	this->averageError = .0;
	this->averageError_Singularity = .0;
	this->globalError = .0;
	this->globalError_Singularity = .0;

	static bool EA_key = true;												// Use elements heigths to define a minimum allowed element heigh for remesh.
	if (EA_key){															// This must be done for the very first time only.
		set_h_min(theMesh);
		EA_key = false;
	}

	int i = 0;
	int nnodes = M_numVertices(theMesh);
	pStore_h_new = 0;														// pointer to NULL
	pStore_h_new = new double[nnodes];						// create array of size number of mesh nodes


	pEntity node, face, tetra;
	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		resetElementAsSingular(node);										// set node as not belonging to a singular element
		pStore_h_new[i] = 1e+10;
		i++;
	}
	VIter_delete(vit);
	
	int n;
	if (theMesh->getDim()==2){
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face)){
				setElementError(face,.0);									// set element error as null
				setLevelOfRefinement(face,0);								// force element level to assume a real value and not a fictitious one like -1000.
				resetElementAsSingular(face);								// set element as not singular
			}
		}
		FIter_delete(fit);
	}
	else{
		RIter rit = M_regionIter (theMesh);
		while ( (tetra = RIter_next(rit)) ){
			if(!theMesh->getRefinementDepth(tetra)){
				setElementError(tetra,.0);
				setLevelOfRefinement(tetra,0);								// force element level to assume a real value and not a fictitious one like -1000.
				resetElementAsSingular(tetra);
			}
		}
		RIter_delete(rit);
	}
}

// calculates an average height for elements before adaptation. that's for Madlib use.
// for each node mesh node: ratio = avg_hnew/avg_hold
// void ErrorAnalysis::calculate_hold_pernode(pMesh theMesh){
// 	if (theMesh->getDim()==2){
// 		pEntity face;
// 		std::map<int,int> facesAroundVertices;
// 		calculateNumFacesAroundVertices(theMesh,facesAroundVertices);
// 		VIter vit = M_vertexIter(theMesh);
// 		while ( (node = VIter_next(vit)) ){
// 			EN_attachDataInt(node,MD_lookupMeshDataId( "hold_pernode" ),0);
// 		}
// 		VIter_delete(vit);
// 		
// 		double h_old;
// 		FIter fit = M_faceIter(m);
// 		while ( (ent = FIter_next(fit)) ){
// 			for (int i=0; i<3; i++){
// 				v = (pEntity)ent->get(0,i);
// 				EN_getDataInt(v,MD_lookupMeshDataId( "hold_pernode" ),&h_old);
// 				
// 				EN_attachDataInt(v,MD_lookupMeshDataId( "hold_pernode" ),h_old);
// 			}
// 		}
// 		FIter_delete(fit);
// 		
// 		vit = M_vertexIter(m);
// 		while ( (ent = VIter_next(vit)) ){
// 			EN_getDataInt(ent,MD_lookupMeshDataId( "hold_pernode" ),&n);
// 			facesAroundVertices[EN_id(ent)] = n;
// 		}
// 		VIter_delete(vit);
// 	}
// }

double ErrorAnalysis::calculate_ErrorSum(pMesh theMesh, bool excludingSingularities=false){
	double error_sum = .0;
	double error;
	if (theMesh->getDim()==2){
		pEntity face;

		if (!excludingSingularities){
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face)){
				error = getElementError(face);
				error_sum -= error*error;
			}
		}
		FIter_delete(fit);
		}
		else{
			FIter fit = M_faceIter (theMesh);
			while ( (face = FIter_next(fit)) ){
				if(!theMesh->getRefinementDepth(face) && !this->isSingular(face)){
					error = getElementError(face);
					error_sum -= error*error;
				}
			}
			FIter_delete(fit);
		}
	}
	else{
		pEntity tetra;
		RIter rit = M_regionIter (theMesh);
		while ( (tetra = RIter_next(rit)) ){
			if(!theMesh->getRefinementDepth(tetra)){
				error_sum += this->getElementError(tetra);
			}
		}
		RIter_delete(rit);
	}
	return error_sum;
}

void ErrorAnalysis::countElements(pMesh theMesh, bool excludingSingularities=false){
	int numElem = 0;
	int numElem_exSing = 0;
	if (theMesh->getDim()==2){
		pEntity face;
		if (!excludingSingularities){
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face)){
				numElem++;
			}
		}
		FIter_delete(fit);
		}
		else{
			FIter fit = M_faceIter (theMesh);
			while ( (face = FIter_next(fit)) ){
				if(!theMesh->getRefinementDepth(face)){
					if (!this->isSingular(face)){
						numElem_exSing++;
					}
				}
			}
			FIter_delete(fit);
		}
	}
	else{
		// 		pEntity tetra;
		// 		RIter rit = M_regionIter (theMesh);
		// 		while ( (tetra = RIter_next(rit)) ){
		// 			if(!theMesh->getRefinementDepth(tetra)){
		// 				if (excludingSingularities && !this->isSingular(tetra)){
		// 					numElem_exSing++;
		// 				}
		// 				else{
		// 					numElem++;
		// 				}
		// 			}
		// 			RIter_delete(rit);
		// 		}
	}
	if (excludingSingularities){
		setNumElements_Singularity(numElem_exSing);
		//cout << " ############## numElem_exSing = " << numElem_exSing << endl;
	}
	else{
		setNumElements(numElem);
	}
	
}

void ErrorAnalysis::calculate_MaxMinErrorData(pMesh theMesh){
	int maxDepth = 0;
	int maxRefLevel = -10;
	int minRefLevel = 10;
	pEntity elem;
	if (theMesh->getDim()==2){
		FIter fit = M_faceIter (theMesh);
		while ( (elem = FIter_next(fit)) ){
			getMaxMinData(theMesh,elem,this,maxDepth,maxRefLevel,minRefLevel);
		}
		FIter_delete(fit);
	}
	else{
		RIter rit = M_regionIter (theMesh);
		while ( (elem = RIter_next(rit)) ){
			getMaxMinData(theMesh,elem,this,maxDepth,maxRefLevel,minRefLevel);
		}
		RIter_delete(rit);
	}
	this->setMaxDepth(maxDepth);
	this->setMaxRefinementFlag(maxRefLevel);
	this->setMinRefinementFlag(minRefLevel);
}

void ErrorAnalysis::monitoring(FIELD field,pMesh theMesh, double tol1, double tol2){
	static int i = 0;
	static bool openfile = true;
	if (openfile){
		fid.open("Error_Analysis_Monitor.csv");
		fid << "---------------------------------------------------------------------------------------------------------------------------------------------------------\n"
			   "  Property    tol1       tol2      GError       SGradNorm     GError_sing    SGradNorm_sing      MError        MError_sing        Adapt      MeshElements\n"
			   "---------------------------------------------------------------------------------------------------------------------------------------------------------\n";
		openfile = false;
	}

	string str;
	switch (field){
	case PRESSURE:
		str = "     P";
		break;
	case SATURATION:
		str = "    Sw";
		break;
	}
	string	str1 = (adapt)?"Y":"N";
	fid << setprecision(2) << scientific;
	fid << str << "      "
	<< tol1 << "   " 
	<< tol2 << "  "
	<< setprecision(5)
	<< getGlobalError() << "   "
	<< getSmoothedGradNorm() << "     "
	<< getGlobalError_Singularity() << "     " 
	<< getSmoothedGradNorm_Singularity() << "     "
	<< getAverageError() << "     "
	<< getAverageError_Singularity() << "          "
	<< str1 << "           "
	<< M_numFaces(theMesh) << endl;
}

void getMaxMinData(pMesh theMesh, pEntity elem, ErrorAnalysis* pEA, int &maxDepth, int &maxRefLevel, int &minRefLevel){
	maxDepth = std::max(maxDepth,theMesh->getRefinementDepth(elem));
	maxRefLevel = std::max(maxRefLevel,pEA->getLevelOfRefinement(elem));
	minRefLevel = std::min(minRefLevel,pEA->getLevelOfRefinement(elem));
}

bool ErrorAnalysis::checkMaximumNumberOfSubdivision(pMesh theMesh, const int &maxNumberOfSubdivision){
	pEntity elem;
	int level = -1000;		/// initialize variable to get maximum number of refinement
	if (theMesh->getDim()==2){
		FIter fit = M_faceIter (theMesh);
		while ( (elem = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(elem)){
				level = std::max(theMesh->getRefinementLevel(elem),level);
			}
		}
		FIter_delete(fit);
	}
	else{
		RIter rit = M_regionIter (theMesh);
		while ( (elem = RIter_next(rit)) ){
			if(!theMesh->getRefinementDepth(elem)){
				level = std::max(theMesh->getRefinementLevel(elem),level);
			}
		}
		RIter_delete(rit);
	}
	return (level==maxNumberOfSubdivision);
}

void ErrorAnalysis::getRefUnrefElementsList(pMesh theMesh, std::list<pEntity> &elementList, std::set<pEntity>& nodesBGMesh){
	#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: getRefUnrefElementsList\tIN\n";
	#endif
	
	// avoid messing things
	if (elementList.size()){
		elementList.clear();
		nodesBGMesh.clear();
	}
	
	pEntity elem;
	std::list<double> heightList;
	double param1 = _pSimPar->Remeshing_param1();
	double param2 = _pSimPar->Remeshing_param2();
	printf("param1: %.5f\tparam2: %.5f\n",param1,param2);
	double h_old, h_new;
	if (!theMesh){
		throw Exception(__LINE__,__FILE__,"NULL Mesh!");
	}
	
	int dim = theMesh->getDim();
	if (dim==2){
// 		double hmin = get_h_min();
// 		avgError = getAverageError();
		FIter fit = M_faceIter (theMesh);
		while ( (elem = FIter_next(fit)) ){
// 			error = getElementError(elem);
 			h_old = getElement_CDL(elem);
// 			h_new = h_old;
			EN_getDataDbl(elem,MD_lookupMeshDataId( "h_new" ),&h_new);
// 			if (fabs(error) > 1e-8){
// 				h_new = h_old*(avgError/error);
// 			}
			EN_attachDataInt(elem,MD_lookupMeshDataId( "elem_to_remove" ),0);
			if ( h_new < param1*h_old){
				elementList.push_back(elem);
				//heightList.push_back(h_new);
				EN_attachDataInt(elem,MD_lookupMeshDataId( "elem_to_remove" ),1);
			}
			if (h_new > param2*h_old){
				elementList.push_back(elem);
				//heightList.push_back(h_new);
				EN_attachDataInt(elem,MD_lookupMeshDataId( "elem_to_remove" ),1);
			}
		}
		FIter_delete(fit);
	}
	else if (dim==3){
		throw Exception(__LINE__,__FILE__,"Under construction!");
	}
	else{
		throw Exception(__LINE__,__FILE__,"Mesh dimension unknown.");
	}
	
	// calculate average element height per node. Each element node have the average new height of all elements choose to be (un)refined
// 	int i;
// 	double height;
// 	pEntity v;
// 	
// 	cout << "Checking mesh\n";
// 	cout << "Elements: " << M_numFaces(theMesh) << endl;
// 	cout << "elementList: " << elementList.size() << endl;
// 	cout << "heightList: " << heightList.size() << endl;
// 	
// 	if (elementList.size() != heightList.size()){
// 		throw Exception(__LINE__,__FILE__,"ElementsList and HeightList sizes don't match!\n");
// 	}
// 	checkMesh(theMesh);
// 	
// 	theMesh->modifyState(0,1);
// 	theMesh->modifyState(0,2);
// 	theMesh->modifyState(2,1);
// 	theMesh->modifyState(2,0);
// 	theMesh->modifyState(1,2);
// 	
// 	/// set to zero all elem_height
// 	pEntity ent;
// 	VIter vit = M_vertexIter(theMesh);
// 	while ( (ent = VIter_next(vit)) ){
// 		EN_attachDataDbl(ent,MD_lookupMeshDataId( "elem_height" ),.0);
// 	}
// 	VIter_delete(vit);
// 	
// 	
// 	int k=0;
// 	std::list<pEntity>::iterator iter1 = elementList.begin();
// 	std::list<double>::iterator iter2 = heightList.begin();
// 	for ( ;iter1 != elementList.end(); iter1++, iter2++){
// 		for (i=0; i<dim+1; i++){
// 			v = (pEntity)(*iter1)->get(0,i);
// 			#ifdef _SEEKFORBUGS_
// 			if (!v){
// 				throw Exception(__LINE__,__FILE__,"Null pointer!!!\n");
// 			}
// 			#endif
// 			nodesBGMesh.insert(v);
// 			height = .0;
// 			EN_getDataDbl(v,MD_lookupMeshDataId( "elem_height" ),&height);
// 			height += *iter2;
// 			EN_attachDataDbl(v,MD_lookupMeshDataId( "elem_height" ),height);
// 		}
// 	}
// 	// last step: do not allow h_new be excessively small. Let's lut a limit: if h_new < h_min then h_new = h_min;(useless!)
// 	std::map<int,int> facesAroundVertices;
// 	calculateNumFacesAroundVertices(theMesh,facesAroundVertices);
// 	int n;
// 	double h_min = get_h_min();
// 	std::set<pEntity>::iterator iter3 = nodesBGMesh.begin();
// 	k = 0;
// 	for ( ;iter3 != nodesBGMesh.end(); iter3++){
// 		v = (pEntity)(*iter3);
// 		height = .0;
// 		n = facesAroundVertices[EN_id(v)];
// 		EN_getDataDbl(v,MD_lookupMeshDataId( "elem_height" ),&height);
// 		height /= n;
// 		if (!n){
// 			cout << "WARNING: division by ZERO!\n";
// 			throw Exception(__LINE__,__FILE__,"WARNING: division by ZERO!");
// 		}
// 		
// 		// do not allow that element with height less than the minimum be removed
// // 		if ( height < h_min ){
// // 			height = h_min;
// // 		}
// 		EN_attachDataDbl(v,MD_lookupMeshDataId( "elem_height" ),height);
// 	}
//	cout << "\n\n##################################################\n";
//	cout << "Number of elements to be removed: " << heightList.size() << endl;
//	cout << "##################################################\n\n";
	//facesAroundVertices.clear();
	heightList.clear();
	#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: getRefUnrefElementsList\tOUT\n";
	#endif
}

void ErrorAnalysis::set_h_min(pMesh theMesh){
	static bool h_min_set = true;
	if ( h_min_set ){
		pEntity elem;
		h_min = 1e+10;
		FIter fit = M_faceIter (theMesh);
		while ( (elem = FIter_next(fit)) ){
			h_min = std::min(h_min, getElement_CDL(elem) );
		}
		FIter_delete(fit);
		h_min /= _pSimPar->Remeshing_param3();
		h_min_set = false;
	}
}
