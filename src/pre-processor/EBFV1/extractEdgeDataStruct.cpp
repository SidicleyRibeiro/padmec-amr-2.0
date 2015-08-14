///*
// * extractEdgeDataStruct.cpp
// *
// *  Created on: Sep 19, 2014
// *      Author: rogerio
// */
//
//#include "GeomData.h"
//
//typedef std::set<int> theSet;
//typedef std::set<int>::iterator theSetIter;
//
//
//struct ElemIndices{
//	theSet myset;
//	int idx;
//};
//
//
//void extractEDS(const char* filename){
//
//	int i,j;
//	int numElem;
//	int* ID,
//	double** coords;
//	int *triangles;
//	int **tetras;
//	int sortedIDs[4];
//
//	readGmshFile(filename,ID,coords,triangles,tetras,numElem);
//
//
//	theSet tmpset;
//	std::map<int,theSet> edges_Map;
//	for (i=0;i<numElem; i++){
//
//		// edges:
//		//			id0 - id1
//		//			id0 - id2
//		//			id0 - id3
//		//			id1 - id2
//		//			id1 - id3
//		//			id2 - id3
//
//		tmpset.insert(tetra[i][0]);
//		tmpset.insert(tetra[i][1]);
//		tmpset.insert(tetra[i][2]);
//		tmpset.insert(tetra[i][3]);
//
//		j = 0;
//		for(theSetIter iter=tmpset.begin(); iter!=tmpset.end(); iter++){
//			sortedIDs[j++] = *iter;
//		}
//
//		theSet set1 = edges_Map[ sortedIDs[0] ];
//		set1.insert( sortedIDs[1] );
//		set1.insert( sortedIDs[2] );
//		set1.insert( sortedIDs[3] );
//		edges_Map[ sortedIDs[0] ] = set1;
//
//		theSet set2 = edges_Map[ sortedIDs[1] ];
//		set2.insert( sortedIDs[2] );
//		set2.insert( sortedIDs[3] );
//		edges_Map[ sortedIDs[1] ] = set2;
//
//		theSet set3 = edges_Map[ sortedIDs[2] ];
//		set3.insert( sortedIDs[3] );
//		edges_Map[ sortedIDs[2] ] = set3;
//	}
//
//	// calculate the number of edges
//	int nedges = 0;
//	for (std::map<int,theSet>::iterator iter = edges_Map.begin(); iter!=edges_Map.end(); iter++){
//		nedges += (int)(*iter)->second.size();
//	}
//
//	// allocate an indexed structure to store the edges
//	Matrix<int> edges;
//	edges.allocateMemory(nedges,2);
//
//	// Transfer edges from edges_Map to the indexed structure
//	int IDs[2], row = 0;
//	for (std::map<int,theSet>::iterator iter = edges_Map.begin(); iter!=edges_Map.end(); iter++){
//		IDs[0] = (*iter)->first;
//		for (theSetIter iters = (*iter)->second.begin(); iters!=(*iter)->second.end(); iters++){
//			IDs[1] = *iters;
//			edges.setRow(row,IDs);
//			row++;
//		}
//	}
//
//	//clean memory
//	for (std::map<int,theSet>::iterator iter = edges_Map.begin(); iter!=edges_Map.end(); iter++){
//		(*iter)->second.clear();
//	}
//	edges_Map.clear();
//
//
//
//}
