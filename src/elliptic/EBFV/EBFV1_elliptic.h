#ifndef _EBFV1_ELLIPTIC_EQUATION_H_
#define _EBFV1_ELLIPTIC_EQUATION_H_

#include "Elliptic_equation.h"
#include "CPU_Profiling.h"

namespace PRS{           // PRS: Petroleum Reservoir Simulator
	/**
	 * For EBFV1 elliptic formulation, a matrix-free procedure is applied to solve pressure field. A pointer-function is used to tell 
	 * Petsc how to make a matrix-vector product and all matrices and vector are passed to this function as a structure.
	 */
	struct Data_struct{
		Mat *E;
		Mat G;
		Mat *F;
		Vec RHS;
		int ndom;
		int nrows;
		int *rows;

		Vec z;
		int F_nrows;
		int F_ncols;
	};

	struct Vectors{
		Vec u_old;
		Vec v_tmp;
		Vec resvec;
		Vec rhs;
		bool key;
	};


	int MatMultUser(Mat mat, Vec u, Vec y);

	class EBFV1_elliptic : public Elliptic_equation{
	public:

		EBFV1_elliptic();
		EBFV1_elliptic(pMesh, PhysicPropData *, SimulatorParameters *, GeomData *, MeshData *);
		~EBFV1_elliptic();
		double solver(pMesh);

//		void getCPUtime(double &assemblyT, double &solverT, double &gradientT, int &KSPiter){
//			assemblyT = _assemblyT;
//			solverT = _solverT;
//			gradientT = _gradientT;
//			KSPiter = _KSPiter;
//		}


	private:

		bool DF_key;
		dblarray Lij;
		Vectors* pVec;
		Mat *EF_multiDom;

		// After fill the matrix, it must be assembled (PETSC)
		int assemblyMatrix(Mat);
		
		// Fill matrices E,G and F for a specific edge and domain. Theses function are called inside a loop of domains
		int divergence_E(Mat E, double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim);
		int divergence_G(Mat G, double *Cij, int edge, int dom, int dom_flag, int idx0_global, int idx1_global, int id0, int id1, int dim);
		int gradient_F_edges(Mat F, double *Cij, int dom, int idx0, int idx1, int id0, int id1, int dim);
		int gradient_F_bdry(Mat, int);
		int F_bdryFaces(pMesh, Mat, const int&);
		int F_bdryEdges(int,Mat);

		 // Solves (EF + G)u = q as showed below:
		 // 			G * u^k+1 = q - EF * u^k		 
		double solveIteratively();

		// Associate to mesh nodes new pressure values computed
		double pressureGradient(pMesh theMesh);
		double updatePressure(pMesh theMesh);
		void calculatePressureGradient();
		void resetPressureGradient();
		void calc_p_grad_1(int,int);
		void calc_p_grad_2(int,int);
		void calc_p_grad_3(int);

		// PETSc matrix for matrix-free procedure
		Mat matrix;

		// Pointer to be used inside matrix-free procedure
		Data_struct *matvec_struct;
		double setMatrixFreeOperation(pMesh);

		/* Group set of calls to compute and assembly matrices E,F and G*/
		double assembly_EFG_RHS(pMesh);

		/* MatMultUser is called several time by Petsc during KSP solver until convergence be reached. */
		static int MatMultUser(Mat mat, Vec u, Vec y){
			void *ctx;
			MatShellGetContext(mat,&ctx);
			Data_struct* mats = (Data_struct*)(ctx);

			// the desired multiplication => A*u = y, where A = EF+G
			// step 1: y = G*u;
			MatMult(mats->G,u,y);

			// step 2:  y = y + (E*F*u)_dom1 + (E*F*u)_dom2 + ... + (E*F*u)_domN
			for (int i=0; i<mats->ndom; i++){
				// z = [F]*u
				VecZeroEntries(mats->z);
				MatMult(mats->F[i],u,mats->z);
				// y = y + ([E]*z)_dom_i
				MatMultAdd(mats->E[i],mats->z,y,y);
			}
			return 0;
		}

		/*set well contribution to right hand side*/
		int wellsContributionToRHS(pMesh theMesh, Vec&);

		/* Matrices E,F and G are assembled for all nodes.
		 * After that, sub matrices are extracted to form the system of equations.
		 * It is not seem to be the wisest way to assembly the final system of equation
		 * but EBFV1 formulation is a little bit complex to make something efficient.*/
		int set_SOE(pMesh theMesh, Mat, Mat&, bool, Vec&, bool, bool);

		/*Free memory from matrices related to Data_struct*/
		double freeMemory();

		double one_eighth;

		// CPU time and monitoring
		int _KSPiter;
		double _assemblyT;
		double _solverT;
		double _gradientT;

		bool firstVTK;


		PhysicPropData *pPPData;
		SimulatorParameters *pSimPar;
		GeomData *pGCData;
		MeshData *pMData;
		pMesh theMesh;
	};
}
#endif
