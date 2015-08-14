/*
 * ADAPTATIVE MESH REFINEMENT CLASS
 *
 * AMR.h
 *
 *  Created on: 24/01/2013
 *      Author: Prof. Rogério Soares
 */

#ifndef AMR_H_
#define AMR_H_

#include "ErrorAnalysis.h"
#include "exportVTK.h"


class AMR{
public:
	AMR(){}
	virtual ~AMR(){}

	virtual void rodar(ErrorAnalysis *pErrorAnalysis, pMesh & theMesh)=0;
};


#endif /* AMR_H_ */
