#include "auxiliar.h"

void failOpeningFile(string filename, int line, const char* cppfile){
	char msg[256]; sprintf(msg,"File '%s' could not be opened or do not exist.\n ",filename.c_str());
	throw Exception(line,cppfile,msg);
}

double F_area(const double *ptn1, const double *ptn2, const double *ptn3){
	double xyz[3][3], n[3];
	for (int i=0;i<3;i++){
		xyz[0][i] = ptn1[i];
		xyz[1][i] = ptn2[i];
		xyz[2][i] = ptn3[i];
	}
	double a[3] = { xyz[1][0]-xyz[0][0], xyz[1][1]-xyz[0][1], xyz[1][2]-xyz[0][2] } ;
	double b[3] = { xyz[2][0]-xyz[0][0], xyz[2][1]-xyz[0][1], xyz[2][2]-xyz[0][2] } ;
	computeCrossProduct(a,b,n);
	return .5*norm(n);
}

double norm(const double *n){
	return sqrt( pow(n[0],2) + pow(n[1],2) + pow(n[2],2) );
}

void getEdgeVector(pEdge edge, dblarray &edgVec){
	double coord1[3];
	double coord2[3];

	V_coord(edge->get(0,0),coord1);			// extract coordenates from ondes
	V_coord(edge->get(0,1),coord2);			// extract coordenates from ondes

	double signal = (EN_id(edge->get(0,0)) < EN_id(edge->get(0,1)))?1.0:-1.0;
	for (int i=0; i<3; i++) edgVec[i] = signal*(coord2[i] - coord1[i]);
}

double E_length(pEdge edge){
	dblarray v(3,.0);
	getEdgeVector(edge,v);
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// multiply a tensor K by a vector vec
void multTensorVector2D(const double *tensor, const double *vec, double *mult){
	mult[0] = tensor[0]*vec[0]+tensor[1]*vec[1];
	mult[1] = tensor[2]*vec[0]+tensor[3]*vec[1];
}

// converts a string read from file to a double
double strToDouble(string &str){
	return strtod(getSubString(str), 0);
}

// converts a string read from file to an interger
int strToInteger(string &str){
	return atoi(getSubString(str));
}

// filters string the numeric part - internal use
const char* getSubString(string &str){
	string::size_type loc = str.find( "=", 0 );
	string numberstr = str.substr(loc+1, str.size()-loc);
	return numberstr.c_str();
}

// Make a vector given to nodes
void makeVector(const double *A, const double *B, double *v){
	for (int i=0; i<3; i++) v[i] = B[i] - A[i];
}

bool isEdgeExternal(int flag){
	return ( flag>=1000 && flag<=1499 );
}

bool isEdgeInternal(int flag){
	return ( flag>=1500 && flag<=1999 );
}

bool isFaceExternal(int flag){
	return ( flag>=2000 && flag<=2499 );
}

bool isFaceInternal(int flag){
	return ( flag>=2500 && flag<=2999 );
}


double* getIdentityMatrix(const int &dim){
	int i,j,k = 0;
	double *Identity = new double[dim*dim];
	for (i=0; i<dim; i++)
		for (j=0; j<dim; j++)
			Identity[k++] = (i==j)?1.0:.0;
	return Identity;
}

int getVertexFlag(pVertex node){
	return getEntityFlag(0,node);
}

int getEdgeFlag(pEdge edge){
	return getEntityFlag(1,edge);
}

int getFaceFlag(pFace face){
	return getEntityFlag(2,face);
}

int getTetraFlag(pRegion tetra){
	return getEntityFlag(3,tetra);
}

int getEntityFlag(int i, pEntity ent){
	pGEntity gEntity;
	switch (i)
	{
	case 0:
		gEntity = (pGEntity)V_whatIn(ent);
		break;
	case 1:
		gEntity = (pGEntity)E_whatIn(ent);
		break;
	case 2:
		gEntity = (pGEntity)F_whatIn(ent);
		break;
	case 3:
		gEntity = (pGEntity)R_whatIn(ent);
		break;
	}
	return (gEntity)?GEN_tag(gEntity):0;
}

double getSmallestEdgeLength(pMesh theMesh){
	pEntity edge;
	double delta_x = 1e30; // infinity delta_x

	// loop over all edges to get the smallest one
	EIter eit = M_edgeIter(theMesh);
	while ( (edge = EIter_next(eit)) ){
		double length = E_length(edge);
		if (delta_x>length)	delta_x = length;
	}
	EIter_delete(eit);
	//printf("smallest length = %f\n",delta_x);
	return P_getMinDbl(delta_x);
}

void replaceAllOccurencesOnString(string &theString, string::size_type size,
		string seekFor, string replaceFor){
	string::size_type pos = 0;
	while (true){
		pos = theString.find(seekFor,pos);
		if (pos != string::npos)
			theString.replace(pos,size,replaceFor);
		else
			break;
	}
}

void getIJnodes(pEdge edge, std::vector<pEntity> &vertex){
	M_GetVertices(edge,vertex);
	if (EN_id(vertex[0])>EN_id(vertex[1])) std::swap(vertex[0],vertex[1]);
}

void F_getEdges(pMesh theMesh, pEntity face, std::vector<pEntity> &edges){
	std::vector<pVertex> vertex;
	M_GetVertices(face,vertex);
	int IDs[3]={EN_id(vertex[0]),EN_id(vertex[1]),EN_id(vertex[2])};
	edges[0] = (pEntity)theMesh->getEdge(theMesh->getVertex(IDs[0]),
			theMesh->getVertex(IDs[1]));
	edges[1] = (pEntity)theMesh->getEdge(theMesh->getVertex(IDs[0]),
			theMesh->getVertex(IDs[2]));
	edges[2] = (pEntity)theMesh->getEdge(theMesh->getVertex(IDs[1]),
			theMesh->getVertex(IDs[2]));

	vertex.clear();
}

VPoint E_getVertexPoint(pEntity edge, int i){
	return ((mVertex*)edge->get(0,i))->point();
}


void E_getVerticesCoord(pEntity edge, double *I, double *J){
	VPoint point_I = ((mVertex*)edge->get(0,0))->point();
	VPoint point_J = ((mVertex*)edge->get(0,1))->point();
	for(int i=0; i<3; i++){
		I[i] = point_I(i);
		J[i] = point_J(i);
	}
}

void E_IJvector(double *I, double *J, double *IJ){
	for (int i=0; i<3; i++) IJ[i] = J[i]-I[i];
}

void E_vertices(pEntity edge, std::vector<pEntity> & vertex){
	vertex[0] = (pEntity)edge->get(0,0);
	vertex[1] = (pEntity)edge->get(0,1);
}

void printSimulationHeader(){
	if (!P_pid())
		std::cout<< "\n\n\t\t\t=====================================================================================\n"
		"\t\t\t\t\tUNIVERSIDADE FEDERAL DE PERNAMBUCO\n"
		"\t\t\t\t\tCENTRO ACADEMICO DO AGRESTE - CAA\n"
		"\t\t\t\t\tNUCLEO DE TECNOLOGIA - NT\n"
		"\t\t\t\t\tPADMEC - DEMEC - CNPq\n"
		"\t\t\t\t\t2008-2014\n"
		"\n\t\t\t\t\tPRS: Parallel Reservoir Simulator\n"
		"\t\t\t\t\tA Finite Volume Scheme for the Two-Phase Oil-Water Flow in 2-D and\n\t\t\t\t\t3-D heterogeneous/anisotropic porous media\n"
		"\n\t\t\t\t\tAuthor: Rogerio Soares da Silva.\n"
		"\t\t\t\t\tNO WARRANTY. IT MEANS: USE IT AT YOUR OWN RISK\n"
		"\t\t\t=====================================================================================\n\n";
}


double ratio(double a, double b);

//double dot(const double* a, const double *b){
//	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
//}

double inner_product(const dblarray &a1, const dblarray &a2){
	return std::inner_product(a1.begin(),a1.end(),a2.begin(),.0);
}

double inner_product(const double* ptr1, const double* ptr2, int n){
	double IP = .0;
	for(int i=0;i<n;i++){
		IP += ptr1[i]*ptr2[i];
	}
	return IP;
}

double ratio(double a, double b){
	return (double)(a/b);
}

void unitary_vector(const dblarray &a1, dblarray &a2){
	std::copy(a1.begin(),a1.end(),a2.begin());
	dblarray norm(a1.size(),sqrt( std::inner_product(a1.begin(),a1.end(),a1.begin(),.0) ));
	std::transform(a2.begin(),a2.end(),norm.begin(),a2.begin(),ratio);
}

void unitary_vector(dblarray &a1){
	dblarray norm(a1.size(),sqrt( std::inner_product(a1.begin(),a1.end(),a1.begin(),.0) ));
	std::transform(a1.begin(),a1.end(),norm.begin(),a1.begin(),ratio);
}

double norm2(const dblarray &a1){
	return sqrt( std::inner_product(a1.begin(),a1.end(),a1.begin(),.0) );
}

double dot(const double *u, const double *v, int dim){
	double _dot = .0;
	for (int i=0;i<dim;i++){
		_dot += u[i]*v[i];
	}
	return _dot;
}

double norm_L2(const double *p, int dim){
	return sqrt( dot(p,p,dim) );
}

void setBarrier(){
	MPI_Barrier(MPI_COMM_WORLD);
	throw 1;
}

void checklinepassing(int line, const char* filepath, std::string somemessage){
	static ofstream fid1;
	static bool openonce = true;
	static int count = 0;

	if (openonce){
		openonce=false;
		char filename[256]; sprintf(filename,"Debbugging-for-rank__%d.txt",P_pid());
		fid1.open(filename);
	}

	fid1 << "File	: " << filepath << std::endl;
	fid1 << "Line	: " << line << std::endl;
	fid1 << "Counter	: " << ++count << std::endl;
	fid1 << "Msg		: " << somemessage << std::endl;
	fid1 << std::endl << std::endl;
}

void convertSecToTime(double t, double *h, double *m, double *s){
	double frac;
	frac = modf(t/3600.,h);
	frac = modf(frac*60.,m);
	frac = modf(frac*60.,s);
}


void LogFiles(double timeStep, double assemblyT, double solverT, double gradT, int KSPiter, double hyperbolicCPU,LOG_FILES LG, string path,
		      bool restart, int last_step, double cumulativeSTime_Restart, double CPUTime_Restart){

	static int step_counter = (restart)?last_step:0;
	static double cumulativeSTime = (restart)?cumulativeSTime_Restart:0;
	static double cumulativeCPU = (restart)?CPUTime_Restart:0; // cumulative CPU time.
	static ofstream fid;

	switch (LG){
	case OPENLG:{
		// print results on file
		if (!P_pid()){
			char filename[256]; sprintf(filename,"%s_simulation-monitor-%d.csv",path.c_str(),P_size());
			if (restart){
				fid.open(filename,ios_base::app);
			}
			else{
				fid.open(filename);
				fid << "#step time-step PVI cumulative-TS assembly-CPU solver-CPU grad-CPU KSP-iter hyperbolic-CPU cumulative-CPU\n";
			}
			std::cout << "\n--------------------------------------------------\nStart Simulation\n--------------------------------------------------\n\n";
		}
	}
	break;
	case UPDATELG:{
		// rank 0 is in charge to print the output CPU-time
		if (!P_pid()){
			cumulativeSTime += timeStep;
			cumulativeCPU += assemblyT + solverT + gradT + hyperbolicCPU;
			fid << scientific << setprecision(8);
			fid << ++step_counter << " " << timeStep << " " << (int)(100*cumulativeSTime/0.07) << " " << cumulativeSTime << " " << assemblyT << " " <<  solverT
				<< " " << gradT << " " << KSPiter << " " << hyperbolicCPU << " " << cumulativeCPU << endl;
		}
	}
	break;
	case CLOSELG:{
		if (!P_pid()){
			fid.close();
		}
	}
	}
}

void STOP(){
	MPI_Barrier(MPI_COMM_WORLD); throw 1;
}

int printMatrixToFile(Mat& m,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	MatView(m,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

int printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

void getEdgesChildren(pEdge edge, mEdge** edgesChildren){
	std::list<mEntity*> edgeLeaves;
	std::list<mEntity*>::iterator iter;
	edge->getLeaves(edgeLeaves);
	iter = edgeLeaves.begin();
	edgesChildren[0] = (mEdge*)(*iter);
	iter++;
	edgesChildren[1] = (mEdge*)(*iter);
	edgeLeaves.clear();
}

void getFacesAroundFace(pFace face, pFace *faces, int &numFaces){
	numFaces = 0;
	// If a face has one or two boundary edge(s), this face will have two or one neighbor(s).
	for (int i=0; i<3; i++){
		pEdge edge = face->get(1,i);
		// if true, take face's neighbor
		if (E_numFaces(edge) == 2){
			if (edge->get(2,0)==face)
				faces[numFaces] = edge->get(2,1);
			else
				faces[numFaces] = edge->get(2,0);
			numFaces++;
		}
	}
}

int EN_getFlag(pEntity ent){
	if (!ent->getClassification()){
		pEntity root = ent->root();
		if (!root){
			throw Exception(__LINE__,__FILE__,"Null root?!\n");
		}
		return (!root->getClassification())?0:GEN_tag( root->getClassification() );
	}
	else{
		return GEN_tag( ent->getClassification() );
	}
}

void getTriVerticesIDs(pFace face, int* IDs){
	IDs[0] = EN_id( face->get(0,0) );
	IDs[1] = EN_id( face->get(0,1) );
	IDs[2] = EN_id( face->get(0,2) );
}

void printFaceIDs(pEntity face){
	int IDs[3];
	getTriVerticesIDs(face,IDs);
	cout << "Face IDs: " << IDs[0] << " " << IDs[1] << " "<< IDs[2] << "\n";
}

void getNeighboursFace(pEntity face, std::vector<pEntity> &neighbousFace, int &size){
	if (!face){
		throw Exception(__LINE__,__FILE__,"Null pointer!.");
	}
	int k = 0;
	for (int i=0; i<3; i++){
		pEdge edge = (pEdge)face->get(1,i);
		if (E_numFaces(edge)==2){
			if (edge->get(2,0) != face){
				neighbousFace[k++] = edge->get(2,0);
			}
			else{
				neighbousFace[k++] = edge->get(2,1);
			}
		}
	}
	//	for (int i=0; i<k; i++){
	//		printFaceIDs(neighbousFace[i]);
	//	}
	//	if (k==2){
	//		cout << "Face is on boundary\n";
	//	}
	size = k;
}

void makeMeshCopy2(pMesh m, PADMEC_mesh *pm, void(*pGetPressure)(int,double&), void(*pGetSaturation)(int,double&)){
	pm->numVertices = M_numVertices(m);
	pm->numElements = M_numFaces(m);
	pm->ID = new int[M_numVertices(m)];
	pm->coords = new double[3*M_numVertices(m)];
	pm->elements = new int[3*M_numFaces(m)];
	pm->field1 = new double[M_numVertices(m)];
	pm->field2= new double[M_numVertices(m)];

	int i;
	int k = 0;
	int iter = 0;
	double coord[3];
	pEntity ent;
	VIter vit = M_vertexIter(m);
	while ( (ent = VIter_next(vit)) ){
		pm->ID[iter] = EN_id(ent);
		V_coord(ent,coord);
		for(i=0;i<3;i++){
			pm->coords[k] = coord[i];
			k++;
		}
		pGetPressure(iter,pm->field1[iter]);
		pGetSaturation(iter,pm->field2[iter]);
		iter++;
	}
	VIter_delete(vit);
	//STOP();


	k = 0;
	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		for(i=0;i<3;i++){
			pm->elements[k] = EN_id((mVertex*)ent->get(0,i));
			k++;
		}
	}
	FIter_delete(fit);
}

void makeMeshCopy2(PADMEC_mesh* pm,pMesh m, void(*pSetPressure)(int,double), void(*pSetSaturation)(int,double)){
	int i,j, k = 0;
	double x,y,z;
	for (i=0; i<pm->numVertices; i++){
		x = pm->coords[k++];
		y = pm->coords[k++];
		z = pm->coords[k++];
		mVertex* v = m->createVertex(pm->ID[i],x,y,z,0);
		pSetPressure(i,pm->field1[i]);
		pSetSaturation(i,pm->field2[i]);
	}

	k = 0;
	mVertex* vertices[3];
	for (i=0; i<pm->numElements; i++){
		vertices[0] = m->getVertex( pm->elements[k++] );
		vertices[1] = m->getVertex( pm->elements[k++] );
		vertices[2] = m->getVertex( pm->elements[k++] );
		m->createFaceWithVertices(vertices[0],vertices[1],vertices[2],0);
	}
	m->modifyState(0,2,0);
	m->modifyState(0,2);
	m->modifyState(2,0);
}

void makeMeshCopy(pMesh m1, pMesh m2, void(*pSetPressure)(pEntity,double), double(*pGetPressure)(pEntity), void(*pSetSaturation)(pEntity,double), double(*pGetSaturation)(pEntity)){
	mEntity* ent;
	mVertex* v;
	//mVertex* vertices[4];
	pGEntity pGEnt;
	double coord[3];
	//int i,j,k,ID[4];
	double val;

	VIter vit = M_vertexIter(m1);
	while ( (ent = VIter_next(vit)) ){
		pGEnt = 0;
		if (ent->getClassification()){
			pGEnt = ent->getClassification();
		}
		V_coord(ent,coord);
		mVertex* v = m2->createVertex(EN_id(ent),coord[0],coord[1],coord[2],pGEnt);
		val = pGetPressure(ent);
		pSetPressure((pEntity)v,val);
		val = pGetSaturation(ent);
		pSetSaturation((pEntity)v,val);

		//cout << pGetPressure(v) << "\t" << pGetSaturation(v) << endl;
	}
	VIter_delete(vit);
	//STOP();

	EIter eit = M_edgeIter(m1);
	while ( (ent = EIter_next(eit)) ){
		pGEnt = 0;
		if (ent->getClassification()){
			pGEnt = ent->getClassification();
		}
		m2->createEdge((mVertex*)ent->get(0,0), (mVertex*)ent->get(0,1),pGEnt);
	}
	EIter_delete(eit);

	if (m1->getDim()==2){
		FIter fit = M_faceIter(m1);
		while ( (ent = FIter_next(fit)) ){
			pGEnt = 0;
			if (ent->getClassification()){
				pGEnt = ent->getClassification();
			}
			m2->createFaceWithVertices((mVertex*)ent->get(0,0), (mVertex*)ent->get(0,1),(mVertex*)ent->get(0,2),pGEnt);
		}
		FIter_delete(fit);
	}
	else{
	}

#ifdef _SEEKFORBUGS2_
	/*	pEntity ent1, ent2;
	int IDs1[4],IDs2[4];

	VIter vit1 = M_vertexIter(m1);
	VIter vit2 = M_vertexIter(m2);
	while ( (ent1 = VIter_next(vit1)) ){
		ent2 = VIter_next(vit2);
		printf("%d [%d]\t%d [%d]\n",EN_id(ent1),GEN_tag(ent1->getClassification()),EN_id(ent2),GEN_tag(ent2->getClassification()));
	}
	VIter_delete(vit1);
	VIter_delete(vit2);
	//exit(1);
	EIter eit1 = M_edgeIter(m1);
	EIter eit2 = M_edgeIter(m2);
	while ( (ent1 = EIter_next(eit1)) ){
		ent2 = EIter_next(eit2);
		printf("%d %d [%d]\t %d %d [%d]\n",
				EN_id((mVertex*)ent1->get(0,0)),EN_id((mVertex*)ent1->get(0,1)),GEN_tag(ent1->getClassification()),
				EN_id((mVertex*)ent2->get(0,0)),EN_id((mVertex*)ent2->get(0,1)),GEN_tag(ent2->getClassification()));
	}
	EIter_delete(eit1);
	EIter_delete(eit2);
	//exit(1);

	// let's check if copy was made successfully
	if (m1->getDim()==2){
		FIter fit1 = M_faceIter(m1);
		FIter fit2 = M_faceIter(m2);
		while ( (ent1 = FIter_next(fit1)) ){
			ent2 = FIter_next(fit2);
			getTriVerticesIDs(ent1,IDs1);
			getTriVerticesIDs(ent2,IDs2);
			printf("%d %d %d [%d]\t%d %d %d [%d]\n",
					IDs1[0],IDs1[1],IDs1[2],GEN_tag(ent1->getClassification()),
					IDs2[0],IDs2[1],IDs2[2],GEN_tag(ent2->getClassification()));
			//			printf("%d %d %d [%d]\t%d %d %d [%d]\n",
			//								IDs1[0],IDs1[1],IDs1[2],0,IDs2[0],IDs2[1],IDs2[2],0);
		}
		FIter_delete(fit1);
		FIter_delete(fit2);


		// let's delete m1 and see if m2 still alive
		FIter fit = M_faceIter(m1);
		while ( (ent = FIter_next(fit)) ){
			m1->DEL(ent);
		}
		FIter_delete(fit);

		EIter eit = M_edgeIter(m1);
		while ( (ent = EIter_next(eit)) ){
			m1->DEL(ent);
		}
		EIter_delete(eit);

		VIter vit = M_vertexIter(m1);
		while ( (ent = VIter_next(vit)) ){
			m1->DEL(ent);
		}
		VIter_delete(vit);

		cout << "m1     m2\n";
		cout << "V: " << M_numVertices(m1);
		cout << "\tV: " << M_numVertices(m2) << endl;
		cout << "E: " << M_numEdges(m1);
		cout << "\tE: " << M_numEdges(m2) << endl;
		cout << "F: " << M_numFaces(m1);
		cout << "\tF: " << M_numFaces(m2) << endl;


		throw Exception(__LINE__,__FILE__,"Finish him!\n");


	}
	else{
		throw Exception(__LINE__,__FILE__,"To be implemented soon...\n");
		//		RIter rit = M_regionIter(m1);
		//		while ( (ent = RIter_next(rit)) ){
		//			v[0] = (mVertex*)ent->get(0,0);
		//			v[1] = (mVertex*)ent->get(0,1);
		//			v[2] = (mVertex*)ent->get(0,2);
		//			v[3] = (mVertex*)ent->get(0,3);
		//			pGEnt = ent->getClassification();
		//			m2->createTetWithVertices(v[0],v[1],v[2],v[3],pGEnt);
		//		}
		//		RIter_delete(rit);
	}*/
#endif
}
void deleteMesh(PADMEC_mesh* pm){
	delete[] pm->ID; pm->ID = 0;
	delete[] pm->coords; pm->coords = 0;
	delete[] pm->elements; pm->elements = 0;
	delete[] pm->field1; pm->field1 = 0;
	delete[] pm->field2; pm->field2 = 0;
}

void deleteMesh(pMesh m){
	// 	cout << "V: " << M_numVertices(m) << endl;
	// 	cout << "E: " << M_numEdges(m) << endl;
	// 	cout << "F: " << M_numFaces(m) << endl;
	// 	cout << "R: " << M_numRegions(m) << endl;
	pEntity ent;
	cout << "Deleting mesh: \n";
	if (m->getDim()==3){
		///cout << "Deleting tetrahedras\n";
		RIter rit = M_regionIter(m);
		while ( (ent = RIter_next(rit)) ){
			m->DEL(ent);
		}
		RIter_delete(rit);
	}

	//	cout << "Deleting triangles\n";
	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		m->DEL(ent);
	}
	FIter_delete(fit);

	//	cout << "Deleting edges\n";
	EIter eit = M_edgeIter(m);
	while ( (ent = EIter_next(eit)) ){
		m->DEL(ent);
	}
	EIter_delete(eit);

	VIter vit = M_vertexIter(m);
	while ( (ent = VIter_next(vit)) ){
		m->DEL(ent);
	}
	VIter_delete(vit);
	// 	cout << "V: " << M_numVertices(m) << endl;
	// 	cout << "E: " << M_numEdges(m) << endl;
	// 	cout << "F: " << M_numFaces(m) << endl;
	// 	cout << "R: " << M_numRegions(m) << endl;

}


/*
 * 
 * seta as entidade das malhas com os flags necessarios para uma nova simulação (condições de contorno). 
 * isso deve vir direto do adaptador. Mas, enquanto isso nao é resolvido por saulo, vamos "adaptando" as coisas.
 * Rogério: 23/05/2013
 * 
 * ---------------------------------------------------------------------------------------------------------------
 * ATENÇÃO: ESTAMOS CONSIDERANDO UM FIVE-SPOT HOMOGENEO, OU SEJA, O MESMO FLAG 3300 PARA TODOS OS TRIANGULOS,
 *          FLAG 2000 PARA AS ARESTAS DO CONTORNO E 10/51 PARA OS POÇOS !!!!!!!
 * * -------------------------------------------------------------------------------------------------------------
 */
void PADMEC_GAMBIARRA(pMesh m){
	cout << "PADMEC_GAMBIARRA   IN\n";
	pEntity ent;

	m->modifyState(2,1);
	m->modifyState(1,2);

	// insere flag nos elementos
	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		ent->classify( m->getGEntity(3300,2) );
	}
	FIter_delete(fit);

	/// first, set to all nodes the triangles' flags
	VIter vit = M_vertexIter(m);
	while ( (ent = VIter_next(vit)) ){
		ent->classify( m->getGEntity(3300,0) );
	}
	VIter_delete(vit);

	/// second, set to all nodes on boundary edges the boundary flag
	mVertex *v1, *v2, *p1, *p2, *injectionWell, *productionWell;
	double coord1[3], coord2[3];
	double xmax = 0.0, xmin = 1e3;
	double ymax = 0.0, ymin = 1e3;

	EIter eit = M_edgeIter(m);
	while ( (ent = EIter_next(eit)) ){
		ent->classify( m->getGEntity(3300,1) );
		if ( E_numFaces(ent)==1 ){
			v1 = (mVertex*)ent->get(0,0);
			v2 = (mVertex*)ent->get(0,1);
			v1->classify( m->getGEntity(2000,0) );
			v2->classify( m->getGEntity(2000,0) );
			ent->classify( m->getGEntity(2000,1) );
			V_coord(v1,coord1);
			V_coord(v2,coord2);

			/// look for dirichlet and neumann nodes
			if (coord1[0] < xmin){
				xmin = coord1[0];
				injectionWell = v1;
			}
			if (coord1[1] < ymin){
				ymin = coord1[1];
			}
			if (coord1[0] > xmax){
				xmax = coord1[0];
			}
			if (coord1[1] > ymax){
				ymax = coord1[1];
			}
			if (coord2[0] < xmin){
				xmin = coord2[0];
			}
			if (coord2[1] < ymin){
				ymin = coord2[1];
			}
			if (coord2[0] > xmax){
				xmax = coord2[0];
			}
			if (coord2[1] > ymax){
				ymax = coord2[1];
			}
			if (coord1[0] >= xmax && coord1[1] >= ymax){
				productionWell = v1;
			}
			if (coord2[0] >= xmax && coord2[1] >= ymax){
				productionWell = v2;
			}
			if (coord1[0] <= xmin && coord1[1] <= ymin){
				injectionWell = v1;
			}
			if (coord2[0] <= xmin && coord2[1] <= ymin){
				injectionWell = v2;
			}
			if (coord1[0] <= xmin && coord1[1] >= ymax){
				p1 = v1;
			}
			if (coord2[0] <= xmin && coord2[1] >= ymax){
				p1 = v2;
			}
			if (coord1[0] >= xmax && coord1[1] <= ymin){
				p2 = v1;
			}
			if (coord2[0] >= xmax && coord2[1] <= ymin){
				p2 = v2;
			}
		}
	}
	EIter_delete(eit);

	/// set dirichlet and neumann nodes flags
	injectionWell->classify( m->getGEntity(10,0) );
	productionWell->classify( m->getGEntity(51,0) );
	p1->classify( m->getGEntity(1100,0) );
	p2->classify( m->getGEntity(1100,0) );
	cout << "PADMEC_GAMBIARRA   OUT\n";
}


void checkNumEdges(pMesh theMesh, int line, char* file){
	if (!M_numEdges(theMesh)){
		throw Exception(line,file,"Number of edges is 0!");
	}
}

void checkMesh(pMesh m){
	ofstream fid;
	fid.open("checkingMesh.txt");
	fid << "Nodes: " << M_numVertices(m) << endl;
	pEntity ent;
	double coords[3];
	VIter vit = M_vertexIter(m);
	while ( (ent = VIter_next(vit)) ){
		fid << EN_id(ent) << "\t";
		V_coord(ent,coords);
		for(int i=0;i<3;i++){
			fid << coords[i] << "\t";
		}
		fid << endl;
	}
	VIter_delete(vit);
	fid << endl;
	fid << endl;
	fid << "Elements: " << M_numFaces(m) << endl;
	int IDs[3];

	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		getTriVerticesIDs(ent,IDs);
		fid << "Face IDs: " << IDs[0] << " " << IDs[1] << " "<< IDs[2] << "\n";
	}
	FIter_delete(fit);
	fid.close();
}

void calculateNumFacesAroundVertices(pMesh m, std::map<int,int> &facesAroundVertices){
	pEntity ent, v;
	int n;


	VIter vit = M_vertexIter(m);
	while ( (ent = VIter_next(vit)) ){
		EN_attachDataInt(ent,MD_lookupMeshDataId( "numfaces" ),0);
	}
	VIter_delete(vit);

	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		for (int i=0; i<3; i++){
			v = (pEntity)ent->get(0,i);
			EN_getDataInt(v,MD_lookupMeshDataId( "numfaces" ),&n);
			EN_attachDataInt(v,MD_lookupMeshDataId( "numfaces" ),++n);
		}
	}
	FIter_delete(fit);

	vit = M_vertexIter(m);
	while ( (ent = VIter_next(vit)) ){
		EN_getDataInt(ent,MD_lookupMeshDataId( "numfaces" ),&n);
		facesAroundVertices[EN_id(ent)] = n;
	}
	VIter_delete(vit);
}

void getdomains(pMesh m, int &n, int* domList){
	set<int> setdomains;
	pEntity ent;
	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		setdomains.insert(GEN_tag(ent->getClassification()));
	}
	FIter_delete(fit);
	cout << "Number of domains: " << setdomains.size() << endl;
	n = (int)setdomains.size();
	domList = new int[n];
	int i = 0;
	set<int>::iterator iter = setdomains.begin();
	for(;iter!=setdomains.end();iter++){
		domList[i] = *iter;
		i++;
	}
}

void readmesh(pMesh m,char* filename){
	cout << "Lendo malha... ";
	ifstream fid;
	fid.open(filename);
	int NbNod;
	int iNod;
	double x,y,z;
	char line[256];
	fid.getline (line,256);
	fid >> NbNod;
	cout << "\nNodes: " << NbNod;
	for(int i=0;i<NbNod;i++){
		fid >> iNod >> x >> y >> z;
		m->createVertex(iNod,x,y,z,0);
	}
	fid.getline (line,256);
	fid.getline (line,256);
	fid.getline (line,256);
	fid >> NbNod;
	cout << "\tElements: " << NbNod << endl;
	int face_id = 0;
	for (int i=0; i<NbNod; i++){
		int iNbNod,iTyp,iGrp,iElm,iNbSub,id;
		mVertex *nod[100];
		fid >> iElm >> iTyp >> iGrp >> iNbSub >> iNbNod;

		//cout << "iTyp: " << iTyp << "\tiGrp: " << iGrp << endl;
		for(int i=0;i<iNbNod;i++){
			fid >> id;
			nod[i] = m->getVertex(id);
		}

		mEntity *theEntity = 0;
		switch(iTyp){
		case 2 :
			theEntity = m->createFaceWithVertices(nod[0],nod[1],nod[2],m->getGEntity(iGrp,2));
			EN_setID((pEntity)theEntity,++face_id);
			break;
		case 1 :
			theEntity = m->createEdge(nod[0],nod[1],m->getGEntity(iGrp,1));
			break;
		case 15 :
			mVertex *v = nod[0];
			v->classify(m->getGEntity(iGrp,0));
			break;
		}
	}
	cout << "OK!\n";
	pVertex node;

	int flag;
	pEntity edge, ent;
	EIter eit = M_edgeIter(m);
	while ( (edge = EIter_next(eit)) ){
		flag = getEdgeFlag(edge);
		for (int i=0; i<2; i++){
			pVertex v = edge->get(0,i);
			if ( !v->getClassification() ){
				v->classify( m->getGEntity(flag,0) );
			}
		}
	}
	EIter_delete(eit);

	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		flag = GEN_tag(ent->getClassification());
		for (int i=0; i<3; i++){
			pVertex v = ent->get(0,i);
			if ( !v->getClassification() ){
				v->classify( m->getGEntity(flag,0) );
			}
		}
	}
	FIter_delete(fit);

	//	VIter vit = M_vertexIter(m);
	//	while ( (node = VIter_next(vit)) ){
	//		cout << "Node " << EN_id(node) << " has classification ";
	//		if ( !node->getClassification() ){
	//			cout << "NULL\n";
	//		}
	//		else{
	//			cout << GEN_tag(node->getClassification()) << endl;
	//		}
	//	}
	//	VIter_delete(vit);
}
