/*
 * CalculateSmoothGradientNorm.cpp
 *
 *  Created on: 24/04/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"


/*
 * Calculate smooth gradient Norm (SGN):
 * 	SNG = sum((|gradSw_I|+|gradSw_J|+|gradSw_K|)/3.0)^2
 * 	where sum is for each triangle element with nodes I, J, K.
 */
void ErrorAnalysis_2D::calculate_SmoothedGradientNorm(pMesh theMesh, SimulatorParameters *pSimPar, GeomData* pGCData, FuncPointer_GetGradient getGradient, FIELD field){
	double grad_0[3], grad_1[3], grad_2[3], SGN, dot_0, dot_1, dot_2;
	int dom_counter, idx, idx0, idx1, idx2, idx0_global, idx1_global, idx2_global;

	pEntity face;
	dom_counter = 0;
	for (SIter_const dom=pSimPar->setDomain_begin(); dom!=pSimPar->setDomain_end();dom++){
		SGN = .0;
		idx = 0;
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face)){
				pGCData->getFace(dom_counter,idx,idx0,idx1,idx2,idx0_global,idx1_global,idx2_global);
				getGradient(field,dom_counter,idx0,idx0_global,grad_0);
				getGradient(field,dom_counter,idx1,idx1_global,grad_1);
				getGradient(field,dom_counter,idx2,idx2_global,grad_2);
				dot_0 = inner_product(grad_0,grad_0,2);
				dot_1 = inner_product(grad_1,grad_1,2);
				dot_2 = inner_product(grad_2,grad_2,2);
				SGN += (dot_0 + dot_1 + dot_2)/3.0;
				idx++;
			}
		}
		FIter_delete(fit);
		// Allow only one loop over elements to calculate error based on saturation field.
		// Saturation gradient is not per domain like pressure.
		if (field == SATURATION){
			break;
		}
		dom_counter++;
	}
	setSmoothedGradNorm(sqrt(SGN));
}

void ErrorAnalysis_2D::calculate_SmoothedGradientNorm_Singularity(pMesh theMesh, SimulatorParameters *pSimPar, GeomData* pGCData, FuncPointer_GetGradient getGradient, FIELD field){
	double grad_0[3], grad_1[3], grad_2[3], SGN, dot_0, dot_1, dot_2, gradDiff;
	int idx, idx0, idx1, idx2, idx0_global, idx1_global, idx2_global, i, dom_counter ;
	pEntity face;

	gradDiff  = .0;
	dom_counter = 0;
	for (SIter_const dom=pSimPar->setDomain_begin(); dom!=pSimPar->setDomain_end();dom++){
		SGN = .0;
		idx = 0;
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face) && !this->isSingular(face)){
				pGCData->getFace(dom_counter,idx,idx0,idx1,idx2,idx0_global,idx1_global,idx2_global);
				getGradient(field,dom_counter,idx0,idx0_global,grad_0);
				getGradient(field,dom_counter,idx1,idx1_global,grad_1);
				getGradient(field,dom_counter,idx2,idx2_global,grad_2);
				dot_0 = inner_product(grad_0,grad_0,2);
				dot_1 = inner_product(grad_1,grad_1,2);
				dot_2 = inner_product(grad_2,grad_2,2);
				SGN += (dot_0 + dot_1 + dot_2)/3.0;
				idx++;
			}
		}
		FIter_delete(fit);
		// Allow only one loop over elements to calculate error based on saturation field. Saturation gradient is not per domain like pressure.
		if (field == SATURATION){
			break;
		}
		dom_counter++;
	}
	setSmoothedGradNorm_Singularity(sqrt(SGN));
}
