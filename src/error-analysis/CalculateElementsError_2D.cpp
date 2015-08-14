/*
 * CalculateElementsError.cpp
 *
 *  Created on: 24/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"

void ErrorAnalysis_2D::calculate_ElementsError(pMesh theMesh, SimulatorParameters *pSimPar,
		                                       GeomData* pGCData, FuncPointer_GetGradient getGradient, FIELD field){

	int dim = theMesh->getDim();
	double v0_coords[3], v1_coords[3], v2_coords[3];	// face's nodes coords array
	double vector_0[3], vector_1[3], vector_2[3];		// face's edge vectors
	double grad_0[3], grad_1[3], grad_2[3];				// face's node gradients
	int idx, idx0, idx1, idx2, idx0_global, idx1_global, idx2_global;

	int dom_counter = 0;
	for (SIter_const dom=pSimPar->setDomain_begin(); dom!=pSimPar->setDomain_end();dom++){		// loop over domains to calculate error element
		idx = 0;
		FIter fit = M_faceIter (theMesh);														// loop over elements
		while (pEntity face = FIter_next(fit)){
			int flag = getFaceFlag(face);
			if ( !theMesh->getRefinementDepth(face) && flag==*dom){
				pGCData->getFace(dom_counter,idx,idx0,idx1,idx2,idx0_global,idx1_global,idx2_global);
				V_coord(face->get(0,0),v0_coords);
				V_coord(face->get(0,1),v1_coords);
				V_coord(face->get(0,2),v2_coords);
				makeVector(v0_coords,v1_coords,vector_0);
				makeVector(v1_coords,v2_coords,vector_1);
				makeVector(v2_coords,v0_coords,vector_2);
				getGradient(field,dom_counter,idx0,idx0_global,grad_0);
				getGradient(field,dom_counter,idx1,idx1_global,grad_1);
				getGradient(field,dom_counter,idx2,idx2_global,grad_2);
				//cout << grad_0[0] << " " << grad_0[1] << ",  " << grad_1[0] << " " << grad_1[1] << ",  " << grad_2[0] << " " << grad_2[1] << endl;
				double delta_grad0[3] = {grad_1[0]-grad_0[0], grad_1[1]-grad_0[1], grad_1[2]-grad_0[2]};
				double delta_grad1[3] = {grad_2[0]-grad_1[0], grad_2[1]-grad_1[1], grad_2[2]-grad_1[2]};
				double delta_grad2[3] = {grad_0[0]-grad_2[0], grad_0[1]-grad_2[1], grad_0[2]-grad_2[2]};
				double error = (fabs(dot(delta_grad0,vector_0,2)) + fabs(dot(delta_grad1,vector_1,2)) + fabs(dot(delta_grad2,vector_2,2)))/3.0;
				this->setElementError(face, sqrt(error));
				idx++;
			}
		}
		FIter_delete(fit);
		// Allow only one loop over elements to calculate error based on saturation field. Saturation gradient is not per domain like pressure.
		dom_counter = (field == SATURATION)?0:dom_counter++;
	}
}
