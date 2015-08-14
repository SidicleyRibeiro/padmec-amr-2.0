/*
 * Boundary_Conditions.cpp
 *
 *  Created on: Jan 29, 2015
 *      Author: rogerio
 */

#include "Boundary_Conditions.h"

double Benchmark3D_case1__ES(double x, double y, double z){
	return 1.0 + sin(pi*x) * sin(pi*(y+1/2)) * sin(pi*(z+1/3));
}


