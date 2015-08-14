#ifndef EXPORT_VTK_H_
#define EXPORT_VTK_H_

#include "PhysicPropData.h"
#include "ErrorAnalysis.h"
#include "GeomData.h"


void exportSolutionToVTK(pMesh,void*,void*,void*,void*,string);
void printVTKHeader(ofstream &, pMesh);
void printVerticesCoordenates(ofstream &, pMesh);
void printElementConnectivities(ofstream&, pMesh, int, int);
void printElementConnectivities(ofstream&, pEntity, int);
void printCellTypeList(ofstream &, int, int);
void printPressure(ofstream &, pMesh, PRS::PhysicPropData*);
void printSaturation(ofstream &, pMesh, PRS::PhysicPropData*);
void printSaturation(ofstream &, pMesh, PRS::PhysicPropData*, PRS::GeomData*);

void printDegreeOfRefinement(ofstream &, pMesh, ErrorAnalysis*);
void printCharacteristicLentgh(ofstream &fid, pMesh theMesh);
void print_hNew(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis);
void print_hOld(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis);
void print_ElementsToBeRemoved(ofstream &fid, pMesh theMesh);
void print_hnew_hold_percentual(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis);

void printElementError(ofstream &, pMesh, ErrorAnalysis*);
void printCharac_Lenth(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis);
void print_Swgrad(ofstream &fid, pMesh theMesh, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);
void print_pw_GradientNorm(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);
void print_Sw_GradientNorm(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);
void print_Sw_GradientNorm2(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);
void print_SingularElements(ofstream& fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis);

// New VTK source file
// --------------------------------------------------------------
void print_headers(ofstream&,int);
void print_coodinates(ofstream&,const GeomData*);
void print_connectivities(ofstream&, const GeomData*);
void print_celltype(ofstream&,int,int);
void print_pressure(ofstream&, const PhysicPropData*);
void print_saturation(ofstream&, const PhysicPropData*);


#endif /*EXPORT_VTK_H_*/
