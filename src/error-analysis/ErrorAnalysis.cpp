#include "ErrorAnalysis.h"

/*
 * Steps for error analysis (taken from Felipe Araujo Dissertation, p47):
 * 
 * step 1: Gradient Recovery (skip this step because gradients have already been calculated.)
 * step 2: Calculate L2 norms of all nodal gradients and elements errors.
 *		 The error analysis works always favorable to security, it means that when evaluating element error related to two or more variables (e.g pressure and saturation)
 *		 the highest error will be considered.
 * 
 * 
 * step 3: Calculate global error for solution Neta
 * 
 *         if (global_error > tolerance)
 * 
 * step 4: calculate mean error based on Neta_1 for whole mesh
 * step 5: calculate d1 (new element dimension) parameter for all mesh elements. Defines the degree of refinement for all elements.
 * 
 * 
 * step 6: if (d1 < d_min) d1 = d_min
 * step 7: recalculate mean error based on the new Neta_2 for whole domain excluding the singularities
 * 
 * Singularities must be understood as regions where gradients are too high which leads the error analysis to set those elements with a very high refinement level
 * 
 * step 8: calculate d2 (new element dimension) parameter for all mesh elements excluding elements with singularities.
 * step 9: if (d2 < d_min) d2 = d_min
 * step 10: d_new = min(d1,d2)
 *          In fact, step 10 is inside step 9, because as long the new element degree
 */

bool calculate_ErrorAnalysis(ErrorAnalysis *pEA, pMesh theMesh, SimulatorParameters *pSimPar, GeomData* pGCData, FuncPointer_GetGradient getGradient){

#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Error Analysis\tIN\n";
#endif

	pEA->_pSimPar = pSimPar;			// for remeshing use
	pEA->adapt = false;					// Says if it will be necessary perform a mesh adaptation
	pEA->initializeParameters(theMesh);	// Every new error analysis, error and level of refinement of all mesh elements 
	// are set to avoid mistakes from previous analysis.

	// Before perform the error analysis, check if some original element number of subdivisions 
	// reached the maximum number of subdivision allowed. If true, return 0
	if ( pEA->checkMaximumNumberOfSubdivision(theMesh,pSimPar->getMax2D()) ){
#ifdef __ERROR_ANALYSIS_DEBUG__
		cout << "Mesh has reached maximum number of element refinement allowed by user.";
#endif
		return pEA->adapt;
	}

	double tol1 = pSimPar->getToleranceForAllElements();
	double tol2 = pSimPar->getToleranceForAllElements_excludingSingularities();
	FIELD field = PRESSURE;
	for (int i=0; i<2; i++){
		pEA->resetAllElementsAsSingular(theMesh);									// do not preserve any element flagged as singular from previous analysis
		pEA->calculate_ElementsError(theMesh,pSimPar,pGCData,getGradient,field);
		pEA->calculate_SmoothedGradientNorm(theMesh,pSimPar,pGCData,getGradient,field);
		pEA->calculate_GlobalError(theMesh);
		if ( pEA->getGlobalError() > tol1 ){
			pEA->adapt = true;
			pEA->calculate_AvgError(theMesh,tol1,false);
			pEA->calculate_DegreeOfRefinement(theMesh,pSimPar,false);
			// if exist some element flagged as singular....
			if (pEA->getNumElements_Singularity()){
				// Repeat calculations for all mesh elements excluding those flagged as singular
				pEA->calculate_SmoothedGradientNorm_Singularity(theMesh,pSimPar,pGCData,getGradient,field);
				pEA->calculate_GlobalError_Singularity(theMesh);
				if (pEA->getGlobalError_Singularity()>tol2){
					pEA->calculate_AvgError(theMesh,tol2,true);
					pEA->calculate_DegreeOfRefinement(theMesh,pSimPar,true);
				}
			}
		}
		else{
			cout << "Adaptation not necessary for field " << i << "\n--------------------------\n";
		}
		pEA->monitoring(field,theMesh,tol1,tol2);
		field = SATURATION;
	}

	// after all fields had been analized, set the correct values for h_new
	if (pEA->adapt){
		pEA->update_h_new(theMesh);
		pEA->calculate_height_ratio(theMesh);
	}
	pEA->deletePointers();
	//STOP();

#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Error Analysis\tOUT\n";
#endif
	return pEA->adapt;
}

