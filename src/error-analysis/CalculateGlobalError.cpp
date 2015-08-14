/*
 * CalculateGlobalError.cpp
 *
 *  Created on: 28/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis.h"

void ErrorAnalysis::calculate_GlobalError(pMesh theMesh){
	double error_sum = calculate_ErrorSum(theMesh,false);
	double smooth_gradNorm = getSmoothedGradNorm();
	#ifdef __ERROR_ANALYSIS_DEBUG__
	if (fabs(error_sum) < 1e-8 || fabs(smooth_gradNorm) < 1e-8){
		char msg[256]; sprintf(msg,"Null data: error_sum = %.8f\t smooth_gradNorm = %.8f.",error_sum,smooth_gradNorm);
		throw Exception(__LINE__,__FILE__,msg);
	}
	#endif
	double global_error = sqrt(fabs(error_sum/smooth_gradNorm));
	setGlobalError(global_error);
}

void ErrorAnalysis::calculate_GlobalError_Singularity(pMesh theMesh){
	double error_sum = calculate_ErrorSum(theMesh,true);
	double smooth_gradNorm = getSmoothedGradNorm_Singularity();
	setGlobalError_Singularity(sqrt(fabs(error_sum/smooth_gradNorm)));
}

// Define an average error. It's an distributed error over the mesh. * 	//eq. 4.24 e 3.2.1 (p�g 32) do algoritmo - tese de Filipe
void ErrorAnalysis::calculate_AvgError(pMesh theMesh, double tol, bool singularity){
	double SGN = (singularity)?getSmoothedGradNorm_Singularity():getSmoothedGradNorm();
	countElements(theMesh,singularity);
	double numElements = (double)(singularity)?getNumElements_Singularity():getNumElements();
	double avgError = (numElements)?(tol*( sqrt(SGN/numElements))):.0;
	if (singularity){
		setAverageError_Singularity(avgError);
	}
	else{
		setAverageError(avgError);
	}
}
