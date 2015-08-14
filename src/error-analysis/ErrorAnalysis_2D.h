/*
 * ErrorAnalysis-2D.h
 *
 *  Created on: 24/04/2012
 *      Author: rogsoares
 */

#ifndef ERRORANALYSIS_2D_H_
#define ERRORANALYSIS_2D_H_

#include "ErrorAnalysis.h"


/*! \class ErrorAnalysis_2D
 *  \brief Used to perform error analysis in 2-D unstructured (triangles) meshes. General analysis functions are implemented in base class.
 */
class ErrorAnalysis_2D : public ErrorAnalysis{

public:

	ErrorAnalysis_2D(){}
	~ErrorAnalysis_2D(){}

private:
	/*! \brief: Realize
	 * \param theMesh FMDB mesh object
	 * \param GetPFuncGrad function pointer to get nodal gradient
	 */
	void calculate_ElementsError(pMesh theMesh, SimulatorParameters* pSimPar, GeomData*, FuncPointer_GetGradient,FIELD);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_SmoothedGradientNorm(pMesh, SimulatorParameters *pSimPar, GeomData*, FuncPointer_GetGradient, FIELD);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_SmoothedGradientNorm_Singularity(pMesh, SimulatorParameters*, GeomData*, FuncPointer_GetGradient, FIELD);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_CharacteristicDimensionLength(pMesh);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_DegreeOfRefinement(pMesh, SimulatorParameters *, bool);
	void calculate_DegreeOfRefinement(pMesh, int, int, bool);
	void calculate_NewElementHeight(pMesh, bool);
	void weightNewElementHeight(pMesh,bool,std::list<pEntity> &);
};



#endif /* ERRORANALYSIS_2D_H_ */
