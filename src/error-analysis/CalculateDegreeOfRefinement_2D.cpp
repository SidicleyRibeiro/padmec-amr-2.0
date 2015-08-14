/*
 * DefineDegreeOfRefinement.cpp
 *
 *  Created on: 01/03/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"

int getLevelToRefine(const double& ratio);
int getLevelToUnrefine(const double& ratio);
void ErrorAnalysis_2D::calculate_DegreeOfRefinement(pMesh theMesh, SimulatorParameters* pSimPar, bool singularity){

	int maxNSub = pSimPar->getMax2D();							// maximum number of subdivisions
	int numSubStep = pSimPar->getNumSubdivision_perStep();		// number of subdivisions per step
	switch( pSimPar->getRefStrategy() ){
		case H_REFINEMENT:
			calculate_DegreeOfRefinement(theMesh,maxNSub,numSubStep,singularity);
			break;
		case ADAPTIVE_REMESHING:
			calculate_NewElementHeight(theMesh,singularity);
			weightNewElementHeight(theMesh,singularity,singularElemList);
			store_h_new(theMesh);
			setElementsAsSingular(singularElemList);
			break;
		case RH_REFINEMENT:
			throw Exception(__LINE__,__FILE__,"Under construction!");
			break;
		default:
			throw Exception(__LINE__,__FILE__,"Unknown adaptation strategy.");
	}
													}
													
// Calculate new element height (h_new) for remeshing. Assigns h_new variable to each element for the entire mesh
// singularity=true means that elements set as singular must not be used
void ErrorAnalysis_2D::calculate_NewElementHeight(pMesh theMesh, bool singularity){
	pEntity face;
	double h_new, h_old, error, h_new_debug;
	double avgError = (singularity)?getAverageError_Singularity():getAverageError();
	
	static int counter = 0;
	FIter fit = M_faceIter (theMesh);
	while (face = FIter_next(fit)){
		if (!isSingular(face)){							// calculate h_new over elements flagged as singular
			error = getElementError(face);				// get element error
			h_old = getElement_CDL(face);				// element characteristic length
			h_new = h_old;								// define a new characteristic element dimension (h_new1 = d1 na tese)
			if (fabs(error) > 1e-8){
				h_new = h_old*(avgError/error);
			}
			EN_attachDataDbl(face,MD_lookupMeshDataId( "h_new" ),h_new);
		}
	}
	FIter_delete(fit);
}

// assign new element height to mesh nodes weigthing by number of sharing elements around node
void ErrorAnalysis_2D::weightNewElementHeight(pMesh theMesh, bool singularity, std::list<pEntity> &singularElemList){
	resetElemHeightPerNode(theMesh);
	int n;
	double h_new, height;
	pEntity node, face;
	double h_min = get_h_min();								// lower limit for element height
	// cumulative element height per node					// singularity = false means that all nodes must have an accumulated height
	FIter fit = M_faceIter (theMesh);						// singularity = true means that all nodes except those flagged as singular
	while (face = FIter_next(fit)){							// must have an accumulated height
		if (!isSingular(face)){
			EN_getDataDbl(face,MD_lookupMeshDataId( "h_new" ),&h_new);
			if (h_new<h_min){
				h_new = h_min;
				singularElemList.push_back(face);
				EN_attachDataDbl(face,MD_lookupMeshDataId( "h_new" ),h_new);
			}
			for (int i=0; i<3; i++){						// accumulative height per node
				node = (pEntity)face->get(0,i);
				if (!isSingular(node)){
					EN_getDataDbl(node,MD_lookupMeshDataId( "elem_height" ),&height);
					EN_attachDataDbl(node,MD_lookupMeshDataId( "elem_height" ),height+h_new);
				}
			}
		}
	}
	FIter_delete(fit);
	
	VIter vit = M_vertexIter(theMesh);						// weighting element height per node
	while ( node = VIter_next(vit) ){
		if (!isSingular(node)){
			EN_getDataInt(node,MD_lookupMeshDataId( "numfaces" ),&n);
			if (!n){
				throw Exception(__LINE__,__FILE__,"Num faces = 0!");
			}

			EN_getDataDbl(node,MD_lookupMeshDataId( "elem_height" ),&height);
			if (height==.0){
				throw Exception(__LINE__,__FILE__,"height = 0!");
			}
			height = (double)(height/n);
			if (height<h_min){
				height = h_min;
			}
			EN_attachDataDbl(node,MD_lookupMeshDataId( "elem_height" ),height);
		}
	}
	VIter_delete(vit);
}

// Define degree of refinement for all elements
void ErrorAnalysis_2D::calculate_DegreeOfRefinement(pMesh theMesh, int maxNSub, int numSubStep, bool singularity){
	
	#ifdef __ERROR_ANALYSIS_DEBUG__
	bool EADebugging = false;
	#endif
	
	double avgError = (singularity)?getAverageError_Singularity():getAverageError();
	int numSigularElements = 0;
	int refunref_level;
	FIter fit = M_faceIter (theMesh);
	while (pEntity face = FIter_next(fit)){
		if (!theMesh->getRefinementDepth(face)){
			double error = getElementError(face); // get element error
			double h_old = getElement_CDL(face);  // element characteristic length
			double h_new = h_old*(avgError/error);		// define a new characteristic element dimension (h_new1 = d1 na tese)
			//eq. 4.27 e 3.2.2 (p�g. 32) - tese de Filipe
			// The following code lines were taken from the error analysis written by Felipe Araujo and update by Bruno Luna in (mar/2006).
			double ratio = h_new/h_old;
			
			if(ratio >= 1.01){		// ----------------------------------- Unrefine
				refunref_level = getLevelToUnrefine(ratio);
			}
			else{					// ----------------------------------- Refine
				refunref_level = getLevelToRefine(ratio);
				if (refunref_level == 4){
					setElementAsSingular(face);
					numSigularElements++;
				}
			}
			
			/* Do not allow root element to be refined more than the maximum number of refinement set by user. If root element's depth 
			 * + leave's level to refine is greater than maximum number of subdivisions, then correct refunref_level. If it happens, 
			 * error analysis MUST NOT be performed once more, because tolerance will not be satisfied due the presence of the singularities.
			 */
			int rootDepth = theMesh->getRefinementDepth(face->root());
			if ( rootDepth + refunref_level > maxNSub ){
				refunref_level = maxNSub - rootDepth;
			}
			else{
				// Before set the new element level or refinement, check if the previous one is less or greater the new level. 
				// Greater level MUST always prevail
				refunref_level = std::max(refunref_level,getLevelOfRefinement(face));
			}
			
			// stops element to be refined completely saving memory.
			if (refunref_level > numSubStep){
				refunref_level = numSubStep;
			}
			setLevelOfRefinement(face,refunref_level);
			
			#ifdef __ERROR_ANALYSIS_DEBUG__
			if (refunref_level){
				EADebugging = true;
			}
			#endif
		}
	}
	FIter_delete(fit);
	
	
	#ifdef __ERROR_ANALYSIS_DEBUG__
	if (!EADebugging){
		//throw Exception(__LINE__,__FILE__,"Global error is greater than tolerance but any element has been flagged to be refined.\n");
	}
	#endif
}

void ErrorAnalysis_2D::calculate_CharacteristicDimensionLength(pMesh theMesh){
	mVertex* v;
	double h_old, h_old_node, edIJ[3][3];
	Trellis_Util::mPoint pts[3];
	pEntity node, face;
	
	// calculate h_old per element
	FIter fit = M_faceIter (theMesh);
	while (face = FIter_next(fit)){
		if (!theMesh->getRefinementDepth(face)){
			for (int i=0; i<3; i++){
				v = (mVertex*)(face->get(0,i));
				pts[i] = v->point();
			}
			for (int i=0; i<3; i++) {
				edIJ[0][i] = pts[0](i)-pts[1](i);
				edIJ[1][i] = pts[0](i)-pts[2](i);
				edIJ[2][i] = pts[1](i)-pts[2](i);
			}
			double dimLen = .0;
			for (int i=0; i<3; i++) {
				dimLen += sqrt(pow(edIJ[i][0],2)+pow(edIJ[i][1],2));
			}
			double CDL = (double)(dimLen/3.0);
			setElement_CDL(face,CDL);
		}
	}
	FIter_delete(fit);

	// calculate h_old per node
	int n, i;
	VIter vit = M_vertexIter(theMesh);										// weighting element height per node
	while ( node = VIter_next(vit) ){
		EN_attachDataDbl(node,MD_lookupMeshDataId( "h_old_node" ),.0);
	}
	VIter_delete(vit);

	fit = M_faceIter (theMesh);
	while (pEntity face = FIter_next(fit)){
		h_old = getElement_CDL(face);
		EN_getDataDbl(face,MD_lookupMeshDataId( "h_old_node" ),&h_old_node);
		for (i=0; i<3; i++){
			node = (pEntity)face->get(0,i);
			EN_attachDataDbl(node,MD_lookupMeshDataId( "h_old_node" ),h_old_node+h_old);
		}
	}
	FIter_delete(fit);

	vit = M_vertexIter(theMesh);
	while ( node = VIter_next(vit) ){
		EN_getDataInt(node,MD_lookupMeshDataId( "numfaces" ),&n);
		EN_getDataDbl(node,MD_lookupMeshDataId( "h_old_node" ),&h_old_node);
		EN_attachDataDbl(node,MD_lookupMeshDataId( "h_old_node" ),(double)(h_old_node/n));
	}
	VIter_delete(vit);
}

int getLevelToRefine(const double& ratio){
	if ( ratio >= 1 && ratio < 1.2 ){
		return 0;
	}
	else if (ratio<1 && ratio>=0.5){
		return 1;
	}
	else if (ratio<0.5 && ratio>=0.25){
		return 2;
	}
	else if (ratio<0.25 && ratio>=0.125){
		return 3;
	}
	else{
		return 4;
	}
}

int getLevelToUnrefine(const double& ratio){
	if(ratio > 4){
		return -3;
	}
	else if (ratio <= 4){
		return -2;
	}
	else if (ratio > 2.5){
		return -1;
	}
}
