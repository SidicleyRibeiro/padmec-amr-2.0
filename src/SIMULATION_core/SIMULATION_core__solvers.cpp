/*
 * SIMULATION_core__solvers.cpp
 *
 *  Created on: 27/08/2012
 *      Author: rogsoares
 */

#include "SIMULATION_core.h"
//#include <sstream>

namespace PRS{

	int SIMULATION_core::solver(){
		double timeStep;
		if (simFlag==STEADY_STATE){
			bool adapt;
			int count = 0;
			do{
				pElliptic_eq->solver(theMesh);
				pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,pGCData,exportSolutionToVTK);
				#ifndef NOADAPTATION
					adapt = adaptation();
				#endif
				count++;
			}while (adapt);
			cout<< "Loops :"<< count <<endl;
		}
		else if (simFlag==TRANSIENT){
			while ( !pSimPar->finishSimulation() ){
				pElliptic_eq->solver(theMesh);
				pHyperbolic_eq->solver(theMesh,timeStep);
				#ifndef NOADAPTATION
					adaptation();
				#endif
			}
		}
		return 0;
	}
}
