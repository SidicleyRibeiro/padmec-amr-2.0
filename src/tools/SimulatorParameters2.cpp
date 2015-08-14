/*
 * SimulatorParameters2.cpp
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#include "SimulatorParameters.h"

namespace PRS{

	void SimulatorParameters::defineExactSolution(){
		bc_external_definition = true;

		// put all command line strings into one string object
		string casename;
		for(int i=2; i<__argc; i++){
			casename.append(__argv[i]);
		}

		if ( casename.find("benchmark3d_case1") != std::string::npos ){
			exact_solution = Benchmark3D_case1__ES;
		}
		else{
			bc_external_definition = false;
		}
	}
}
