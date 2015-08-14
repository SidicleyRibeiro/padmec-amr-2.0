#include "EBFV/EBFV1_elliptic.h"

namespace PRS{

	int EBFV1_elliptic::gradient_F_bdry(Mat F, int dom){
		double Dij[3], Fij_column1[4], Fij_column2[4], volumeI, volumeJ, volumeK;
		double v0[3], v1[3], v2[3];
		int ndom, nedges, nfaces, id0, id1, id2, idx_0, idx_1, idx_2, idx0_global, idx1_global, idx2_global, i, j, pos1, pos2, pos3;

		double sum[3]={.0,.0,.0};
		if (pGCData->getMeshDim()==2){
			nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
			for (j = 0; j<nedges; j++){
				pGCData->getBdryEdge(dom,j,idx_0,idx_1);
				pGCData->getBdryID(dom,idx_0,idx_1,id0,id1);
				pGCData->getBdryVolume(dom,idx_0,idx_1,volumeI,volumeJ);
				pGCData->getDij(dom,j,Dij);
				Fij_column1 [0] = 5.*Dij[0]/(6.*volumeI);
				Fij_column1 [1] = 5.*Dij[1]/(6.*volumeI);
				Fij_column1 [2] = Dij[0]/(6.*volumeJ);
				Fij_column1 [3] = Dij[1]/(6.*volumeJ);
				Fij_column2 [0] = Dij[0]/(6.*volumeI);
				Fij_column2 [1] = Dij[1]/(6.*volumeI);
				Fij_column2 [2] = 5.*Dij[0]/(6.*volumeJ);
				Fij_column2 [3] = 5.*Dij[1]/(6.*volumeJ);
				id0 = pMData->get_AppToPETSc_Ordering(id0);
				id1 = pMData->get_AppToPETSc_Ordering(id1);
				pos1 = 2*(id0-1);
				pos2 = 2*(id1-1);
				int idxm[4] = {pos1,pos1+1,pos2,pos2+1};
				int idxn[2] = {id0-1,id1-1};
				MatSetValues(F,4,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
				MatSetValues(F,4,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
			}
		}
		else{
			nfaces = pGCData->getNumBdryFacesPerDomain(dom);
			//cout << "Faces: " << nfaces << endl;
			for (j = 0; j<nfaces; j++){
				pGCData->getBdryFace(dom,j,idx_0,idx_1,idx_2,idx0_global,idx1_global,idx2_global);
				pGCData->getBdryID(dom,idx_0,idx_1,idx_2,id0,id1,id2);
				pGCData->getBdryVolume(dom,idx_0,idx_1,idx_2,volumeI,volumeJ,volumeK);
				pGCData->getDij(dom,j,Dij);

//				sum[0] += Dij[0];
//				sum[1] += Dij[1];
//				sum[2] += Dij[2];
//				cout << setprecision(8) << "Dij_x: " << Dij[0] << endl;


				double tmp[3] = {1./(8.*volumeI), 1./(8.*volumeJ), 1./(8.*volumeK)};
				double aux[3][3] = {{6.*tmp[0],tmp[0],tmp[0]},{tmp[1],6.*tmp[1],tmp[1]},{tmp[2],tmp[2],6.*tmp[2]}};

				// fill edge matrix
				double Fij_column1[9], Fij_column2[9], Fij_column3[9];
				for (i=0; i<3; i++){
					Fij_column1[3*i] =   aux[i][0]*Dij[0];
					Fij_column1[3*i+1] = aux[i][0]*Dij[1];
					Fij_column1[3*i+2] = aux[i][0]*Dij[2];

					Fij_column2[3*i] = aux[i][1]*Dij[0];
					Fij_column2[3*i+1] = aux[i][1]*Dij[1];
					Fij_column2[3*i+2] = aux[i][1]*Dij[2];

					Fij_column3[3*i] = aux[i][2]*Dij[0];
					Fij_column3[3*i+1] = aux[i][2]*Dij[1];
					Fij_column3[3*i+2] = aux[i][2]*Dij[2];
				}

				// index for global Fg
				// where edge matrix must be assembled.
				id0 = pMData->get_AppToPETSc_Ordering(id0);
				id1 = pMData->get_AppToPETSc_Ordering(id1);
				id2 = pMData->get_AppToPETSc_Ordering(id2);

				int pos1 = 3*(id0-1);
				int pos2 = 3*(id1-1);
				int pos3 = 3*(id2-1);

				//cout << "id: " << id0 << "  " << id1 << "  " << id2 << endl;

				int idxm[9] = {pos1,pos1+1,pos1+2, pos2,pos2+1,pos2+2, pos3,pos3+1,pos3+2};
				int idxn[3] = {id0-1, id1-1, id2-1};

				MatSetValues(F,9,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
				MatSetValues(F,9,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
				MatSetValues(F,9,idxm,1,&idxn[2],Fij_column3,ADD_VALUES);

				//break;

//				double C[3] = {one_eighth*volumeI, one_eighth*volumeJ, one_eighth*volumeK };
//				for (i=0;i<3;i++){
//					v0[i] = C[0]*Dij[i];
//					v1[i] = C[1]*Dij[i];
//					v2[i] = C[2]*Dij[i];
//				}
//				double Fij[27] = {6.0*v0[0], 6.0*v0[1], 6.0*v0[2],     v0[0],     v0[1],     v0[2],     v0[0],     v0[1],     v0[2],
//						              v1[0],     v1[1],     v1[2], 6.0*v1[0], 6.0*v1[1], 6.0*v1[2],     v1[0],     v1[1],     v1[2],
//						              v2[0],     v2[1],     v2[2],     v2[0],     v2[1],     v2[2], 6.0*v2[0], 6.0*v2[1], 6.0*v2[2]};
//
//				pos1 = 3*(pMData->get_AppToPETSc_Ordering(id0)-1);
//				pos2 = 3*(pMData->get_AppToPETSc_Ordering(id1)-1);
//				pos3 = 3*(pMData->get_AppToPETSc_Ordering(id2)-1);
//
//				int idxm[9] = {pos1,pos1+1,pos1+2, pos2,pos2+1,pos2+2, pos3,pos3+1,pos3+2};
//				int idxn[3] = {id0-1, id1-1, id2-1};
//				MatSetValues(F,9,idxm,3,idxn,Fij,ADD_VALUES);
			}
//			cout << setprecision(8) << scientific;
//			cout << "Domain: " << dom << "\t" << sum[0] << " " << sum[1] << " " << sum[2] << "\n";

//			if ( fabs(sum[0])>1e-10 || fabs(sum[0])>1e-10 || fabs(sum[0])>1e-10){
//				throw Exception(__LINE__,__FILE__,"nao fecha!");
//			}
		}
		return 0;
	}
}
