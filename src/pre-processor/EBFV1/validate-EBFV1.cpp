///**
// * **********************************************************************************
// * Implemented by: Rogerio Soares
// * Date: 2011
// * Updated: 03/Set/2012
// *
// *
// *
// * NAO É CONFIÁVEL ESTE PROGRAMA!!!!!!!!!!!!!1
// ***********************************************************************************
// *
// * Validates Cij and Dij coefficients calculation computing nodal control volumes for 2-D and 3-D meshes.
// * Cij and Dij summation for each node should be null or less the 1e-15 to guarrantee the correct calculation.
// *
// */
//
//#include "EBFV1__pre-processors.h"
//
//void validate_EBFV1(GeomData *pGCData, pMesh theMesh, std::set<int> &setOfDomain){
//	cout<< "validate_EBFV1"<<endl;
//	if (!setOfDomain.size()){
//		throw Exception(__LINE__,__FILE__,"Number of domains NULL?\n");
//	}
//	// loop over all domains
//	for (std::set<int>::iterator iter=setOfDomain.begin(); iter!=setOfDomain.end(); iter++){
//		//double cvSum = .0;								// control volumes summation
//// 		mapSUM mapNodeCVsum; 							// map surface control volume for all mesh node IDs
////
// 		int flag_domain = *iter;						// domain iterator
//// 		initializeCoeffSummation(pGCData,theMesh,mapNodeCVsum);
//// 		getEdgesContribution(pGCData,theMesh,mapNodeCVsum,flag_domain);
//// 		getBoundaryEdgescontribution(pGCData,theMesh,mapNodeCVsum,flag_domain);
//// 		//		getFacesContribution(pGCData,theMesh,mapNodeCVsum,flag_domain);
//// 		//		getRemoteNodesContribution(pGCData,theMesh,mapNodeCVsum,flag_domain);
//// 		calculateCVSummation(pGCData,theMesh,mapNodeCVsum,flag_domain);
//		double total_volume = .0;
//		pEntity node;
//		VIter vit = M_vertexIter(theMesh);
//		while ( node = VIter_next(vit) ){
//			//cout << EN_id(node) << "\t" << pGCData->getVolume(node,flag_domain) << endl;
//			total_volume += pGCData->getVolume(node,flag_domain);
//		}
//		VIter_delete(vit);
//		cout << "total_volume = " << total_volume << endl;
//
//
//		theMesh->modifyState(0,1);
//		std::vector<double> Cij(3), Dij(3), sum_Cij(3,.0);
//		vit = M_vertexIter(theMesh);
//		while ( node = VIter_next(vit) ){
//			int ID = EN_id(node);
//			double sumx = .0;
//			double sumy = .0;
//			for(int i=0; i<V_numEdges(node); i++){
//				pEdge edge = (pEntity)node->get(1,i);
//				int id0 = EN_id(edge->get(0,0));
//				int id1 = EN_id(edge->get(0,1));
//				pGCData->getCij(edge,flag_domain,Cij);
//				if ( id0 != ID ){
//					if (id0 < ID){
//						sumx += -Cij[0];
//						sumy += -Cij[1];
//					}
//					else{
//						sumx += Cij[0];
//						sumy += Cij[1];
//					}
//				}
//				if (id1 != ID){
//					if (id1 < ID){
//						sumx += -Cij[0];
//						sumy += -Cij[1];
//					}
//					else{
//						sumx += Cij[0];
//						sumy += Cij[1];
//					}
//				}
//
//				if (E_numFaces(edge)==1){
//					pGCData->getDij(edge,flag_domain,Dij);
//					sumx += Dij[0];
//					sumy += Dij[1];
//				}
//			}
//			//cout << ID-1 << "\tsumx = " << sumx << "\tsumy = " << sumy  << endl;
//		}
//		VIter_delete(vit);
//
//
//		std::vector<double> sum_I(3,.0);
//		//cout << "sum_I[0] = " << sum_I[0] << "\tsum_I[1] = " << sum_I[1] << endl;
//		pEntity edge;
//		EIter eit = M_edgeIter(theMesh);
//		while (edge  = EIter_next(eit) ){
//			if ( E_numFaces(edge)==1 ){
//				pGCData->getDij(edge,flag_domain,Dij);
//				sum_I[0] += Dij[0];
//				sum_I[1] += Dij[1];
//				//cout << "sum_I[0] = " << sum_I[0] << "\tsum_I[1] = " << sum_I[1] << endl;
//			}
//		}
//		EIter_delete(eit);
//		cout << "sum_I[0] = " << sum_I[0] << "\tsum_I[1] = " << sum_I[1] << endl;
//		//STOP();
//	}
//}
//
//void initializeCoeffSummation(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum){
//	std::vector<double> vec(3,.0);
//	VIter vit = M_vertexIter(theMesh);
//	while ( pEntity node = VIter_next(vit) ){
//		mapNodeCVsum[ EN_id(node) ] = vec;			// initialize surface control volumes
//	}
//	VIter_delete(vit);
//}
//
//void getEdgesContribution(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum, int flag_domain){
//	std::vector<double> Cij(3), sum_I(3), sum_J(3);
//	EIter eit = M_edgeIter(theMesh);
//	while ( pEntity edge = EIter_next(eit) ){
//		if (!theMesh->getRefinementDepth(edge) && pGCData->edgeBelongToDomain(edge,flag_domain) ){
//			int ID_I = EN_id(edge->get(0,0));		// edge ID I
//			int ID_J = EN_id(edge->get(0,1));		// edge ID J
//			pGCData->getCij(edge,flag_domain,Cij);	// Cij coefficient
//			//			cout << "Cij = " << Cij[0] << "\t" << Cij[1] << "\t" << Cij[2] << "\n";
//			sum_I = mapNodeCVsum[ ID_I ];			// get SC volume summation node ID_I
//			sum_J = mapNodeCVsum[ ID_J ];			// get SC volume summation node ID_I
//			double sinal = (ID_I < ID_J)?1.0:-1.0;	// verify edge orientation
//			for (int i=0; i<3; i++){				// perform summation
//				sum_I[i] += sinal*Cij[i];
//				sum_J[i] += -sinal*Cij[i];
//			}
//			mapNodeCVsum[ ID_I ] = sum_I;			// update SC volume summation node ID_I
//			mapNodeCVsum[ ID_J ] = sum_J;			// update SC volume summation node ID_J
//		}
//	}
//	EIter_delete(eit);
//}
//
//void getBoundaryEdgescontribution(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum, int flag_domain){
//	if (theMesh->getDim()==3){
//		return;
//	}
//
//	std::vector<double> Dij(3), sum_I(3), sum_J(3);
//	EIter eit = M_edgeIter(theMesh);
//	while ( pEntity edge = EIter_next(eit) ){
//		if (!theMesh->getRefinementDepth(edge)){
//			// take only boundary faces
//			if ( pGCData->belongsToBoundary(edge) )
//				if ( pGCData->edgeBelongToDomain(edge,3300) ){
//					pGCData->getDij(edge,flag_domain,Dij);
//					int ID_I = EN_id(edge->get(0,0));		// edge ID I
//					int ID_J = EN_id(edge->get(0,1));		// edge ID J
//					sum_I = mapNodeCVsum[ ID_I ];			// get SC volume summation node ID_I
//					sum_J = mapNodeCVsum[ ID_J ];			// get SC volume summation node ID_I
//					for (int i=0; i<3; i++){				// perform summation
//						sum_I[i] += Dij[i];
//						sum_J[i] += Dij[i];
//					}
//					mapNodeCVsum[ ID_I ] = sum_I;			// update SC volume summation node ID_I
//					mapNodeCVsum[ ID_J ] = sum_J;			// update SC volume summation node ID_J
//				}
//		}
//	}
//	EIter_delete(eit);
//}
//
//void calculateCVSummation(GeomData *pGCData, pMesh theMesh, mapSUM &mapNodeCVsum, int flag_domain){
//	// every time coefficient validation is involked a new text file is created with all out put
//	static int CVcalling = 0;
//	char filename[256];
//	sprintf(filename,"logs/CoefficientValidation-%d.txt",++CVcalling);
//	ofstream fid;
//	fid.open(filename);
//	fid << "Summation of all surface vectors for each control volume: x - y - z\n";
//	fid << "columns: 1 - ID, 2 - x, 3 - y, 4 - z\n\n";
//	fid << scientific << setprecision(8);
//
//	bool PASSED = true;							// check validation
//	const double accuracy = 1e-10;				// accuracy
//	std::vector<double> sum_I(3), sum(3);		// normal surface vector must be 0.0 for all coordinates (x,y,z)
//
//// 	VIter vit = M_vertexIter(theMesh);
//// 	while ( pEntity node = VIter_next(vit) ){
//// 		if ( pGCData->nodeBelongToDomain(node,flag_domain) ){	// pass only node over domain dom
//// 			int ID_I = EN_id(node);						// get node ID
//// 			sum_I = mapNodeCVsum[ ID_I ];				// get x,y,z summations components
//// 			fid << ID_I << "  " << sum_I[0] << " " << sum_I[1] << " " << sum_I[2] << endl;
//// 			for (int i=0; i<3; i++){
//// 				sum[i] += sum_I[i];
//// 				if (fabs(sum[i] > accuracy)){
//// 					fid << "WARNNING: node [" << EN_id(node) << "] failed while checking Cij/Dij control volume surface.\n";
//// 					fid << scientific << setprecision(8) << sum_I[0] << " " << sum_I[1] << " " << sum_I[2] << endl;
//// 					PASSED = false;
//// 				}
//// 			}
//// 		}
//// 	}
//// 	VIter_delete(vit);
//
//	pEntity node;
//	std::vector<double> sumI(3);
//	double volume  =.0;
//	VIter vit = M_vertexIter(theMesh);
//	while ( (node = VIter_next(vit)) ){
//		volume += pGCData->getVolume(node,3300);
//	}
//	VIter_delete(vit);
//
//	if (volume > 1.0){
//		char msg[256];
//		sprintf(msg,"Volume = %.8f",volume);
//		throw Exception(__LINE__,__FILE__,msg);
//	}
//
//	//	cout << "---------------------------------------------------------------------------------------\n";
//	//	cout << "Results:";
//	//	cout << "Total summation: " << sum[0] << " " << sum[1] << " " << sum[2] << "\n\n";
//	//	cout << "---------------------------------------------------------------------------------------\n";
//	fid << "---------------------------------------------------------------------------------------\n";
//	fid << "Results:";
//	fid << "Total summation: " << sum[0] << " " << sum[1] << " " << sum[2] << "\n\n";
//	fid << "---------------------------------------------------------------------------------------\n";
//
//	mapNodeCVsum.clear();
//
//	if (PASSED){
//		//printf("Geometric Coefficient Validation PASSED for domain %d\n",flag_domain);
//		fid << "Geometric Coefficient Validation PASSED for domain " << flag_domain << endl;
//	}
//	else{
//		fid << "Coefficient validation FAILED!!!!\n";
//		throw Exception(__LINE__,__FILE__,"Coefficient validation FAILED!!!!\n");
//	}
//}
//
//
////void getFacesContribution(GeomData *pGCData, pMesh theMesh, std::set<int> &setOfDomain, int flag_domain){
////	if (theMesh->getDim()!=3) return;
////	printf("\n\n");
////	pEntity face;
////	std::vector<double> Dij(3), sumI(3), sumJ(3), sumK(3);
////	FIter fit = M_faceIter(theMesh);
////	while ( (face = FIter_next(fit)) ){
////		int flag = GEN_tag(face->getClassification());
////		// take only boundary faces
////		std::set<int>::iterator iter = setOfDomain.find(flag);
////		if ( iter==setOfDomain.end() ){
////			pEntity node[3] = {(pEntity)face->get(0,0), (pEntity)face->get(0,1), (pEntity)face->get(0,2)};
////			if (pGCData->getDij(face,flag_domain,Dij)){
////				pGCData->getSumIJ(node[0],sumI);
////				pGCData->getSumIJ(node[1],sumJ);
////				pGCData->getSumIJ(node[2],sumK);
////				for (int i=0; i<3; i++){
////					sumI[i] += Dij[i];
////					sumJ[i] += Dij[i];
////					sumK[i] += Dij[i];
////				}
////				pGCData->setSumIJ(node[0],sumI);
////				pGCData->setSumIJ(node[1],sumJ);
////				pGCData->setSumIJ(node[2],sumK);
////			}
////		}
////	}
////	FIter_delete(fit);
////}
//
////void getRemoteNodesContribution(GeomData *pGCData, pMesh theMesh, std::set<int> &setOfDomain, int flag_domain){
////	pEntity node;
////	std::vector<double> sumI(3);
////	VIter vit = M_vertexIter(theMesh);
////	while ( (node = VIter_next(vit)) ){
////		// take only nodes with remote copies
////		if (M_numRemoteCopies(theMesh,node)){
////			if (pGCData->getVolume(node,flag_domain) != .0){
////				pGCData->getSumIJ(node,sumI);
////			}
////		}
////	}
////	VIter_delete(vit);
////}
