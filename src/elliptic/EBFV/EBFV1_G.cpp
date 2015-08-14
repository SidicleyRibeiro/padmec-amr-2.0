#include "EBFV/EBFV1_elliptic.h"

// =============================================================================
// Matrix G is equivalent to the matrix of discretization of elliptic terms of
// DFVC method in a dual Voronoi mesh. G is assembled with a loop over all edges.
// In parallel, if an edge is on the partition boundary, Gij must be divided by
// (NRC) Number of Remote Copies. Each partition will add Gij/NRC to the global
// matrix E and the resulting Gij contribution will be the same as it was assem-
// bled in serial.
// =============================================================================

namespace PRS{
//	int EBFV1_elliptic::divergence_G(Mat G, pEntity edge, const int &dom, int dim, double *Cij){
int EBFV1_elliptic::divergence_G(Mat G, double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim){
		int i, j;
		double versor[3], length, sign = 1.0;
		if (id0 > id1){
			std::swap(id0,id1);
			sign = -sign;
		}

		// index for global Fg
		id0 = pMData->get_AppToPETSc_Ordering(id0);
		id1 = pMData->get_AppToPETSc_Ordering(id1);

		// permeability tensor
		const double *Permeability = pSimPar->getPermeability(dom_flag);

		double Sw_I, Sw_J, MobI, MobJ, MobIJ;
		pPPData->getSaturation(idx0_global,Sw_I);
		pPPData->getSaturation(idx1_global,Sw_J);
		MobI = pPPData->getTotalMobility(Sw_I);
		MobJ = pPPData->getTotalMobility(Sw_J);
		MobIJ = 0.5*(MobI + MobJ);

		pGCData->getLength(dom,edge,length);
		pGCData->getVersor(dom,edge,versor);
		for (i=0; i<dim; i++){
			versor[i] *= sign;
		}

		// ####################################################################
		// what is calculated below: Gij = -(K*Mob_IJ)*versor/(length)*Cij*|1 -1|
		//											                    |-1 1|
		// ####################################################################
		// KL = Permeability * versor
		double KL[3] = {.0, .0, .0};

		int pos = 0;
		if ( pSimPar->is_K_Isotropic() ){
			for (i=0; i<dim; i++){
				KL[i] = Permeability[pos]*versor[i];
				pos += dim+1;
			}
		}
		else{
			for (i=0; i<dim; i++)
				for (j=0; j<dim; j++) KL[i] += Permeability[dim*i+j]*versor[j];
		}

		double aux = .0;
		for (i=0; i<dim; i++) aux += Cij[i]*KL[i];
		aux /= -length;
		double nrc = 1.0;
		//			// if an edge is on partition boundary, the result from all processors
		//			//  contribution for this edge must be equal as if was ran in serial
		//			// Gij (serial) = Gij (parallel) if edge is on partition boundary
		aux /= nrc;

		double Gij_I[2] = { MobIJ*aux,-MobIJ*aux};
		double Gij_J[2] = {-MobIJ*aux, MobIJ*aux};
		int rows[2] = {id0-1,id1-1};	// where Gij must be assembled in global G
		int cols[2] = {id0-1,id1-1};	// where Gij must be assembled in global G

		MatSetValues(G,1,&rows[0],2,cols,Gij_I,ADD_VALUES);
		MatSetValues(G,1,&rows[1],2,cols,Gij_J,ADD_VALUES);
		return 0;
	}
}
