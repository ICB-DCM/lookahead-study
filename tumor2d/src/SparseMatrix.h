#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

#include "VoronoiDiagramExtended.h"

#define DEBUG 0

#define COMPARE_WITH_EXACT_SOLUTION 0

typedef struct _SparseMatrix SparseMatrix;
typedef struct _SparseVector SparseVector;

extern SparseMatrix *sA;
extern SparseMatrix *M;
extern SparseMatrix *J;
extern float *x;
extern float *b;
extern float *v0;
extern double *v0d;
extern double *v1d;
extern double *v2d;
extern double *v3d;
extern double *v4d;
extern double *v5d;
extern double *v6d;
extern float *v1;
extern float *v2;
extern float *v3;
extern float *v4;
extern float *v5;
extern float *v6;
extern float *v7;
extern float *v8;

// voronoi cell type
struct _SparseMatrix {
	
	float **A; // values
	int  **JA;  // column index of values	
	int *sizeA;
	int *maxSizeA;
	//int *IA;  // indeces of first element of each row
	//int sizeIA;
	
	
	int dimI, dimJ; // dimensions
	
	
	// methodes
	static SparseMatrix *newSparseMatrix( int, int);
	static void deleteSparseMatrix(SparseMatrix *matrix);
	
	float get( int i, int j);

	void set( int i, int j, float value);
	void setLast( int i, int j, float value);

	void add( int i, int j, float value);

	void resetRow( int i);

};

struct _SparseVector {
	
	float *v; // values
	int sizev;
	int *Iv;  // row index of values	
	
	int dim; // dimensions
	
	
	// methodes
	void newSparseVector( int);
	void deleteSparseVector();
	
	float get( int i);
	float set( int i);
};
	
	
void sparseMatrixVectorProduct( SparseMatrix *sA, float *b, float *x);
void sparseMatrixVectorProduct( SparseMatrix *sA, double *b, double *x);
void vectorSparseMatrixProduct( float *b, SparseMatrix *sA, float *x);
void vectorSparseMatrixProduct( double *b, SparseMatrix *sA, double *x);

double UpdateSystemImplicitSparse( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference);
double UpdateSystemNewtonSparse( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference);
double UpdateSystemNonLinearCGSparse( VoronoiDiagram *voronoiDiagram, double time, double end_time, double timeStep);
double UpdateSystemNewtonCGSparse( VoronoiDiagram *voronoiDiagram, double time, double end_time, double timeStep);
double UpdateGrowthFactorsNonLinearCGSparse( VoronoiDiagram *voronoiDiagram, double time, double end_time, double timeStep);

void setupMatrixImplicitSparse( VoronoiDiagram *voronoiDiagram, double timeStep, char molecule);
void JacobiPreconditioner( SparseMatrix *A, SparseMatrix *M);
void SuccessiveOverRelaxationPreconditioner( SparseMatrix *A, SparseMatrix *M);

void SolveExplicit( SparseMatrix *A, float *b, float *x);
int SolveBiCGSTAB( SparseMatrix *A, float *b, float *x, int maxit, float minerr);
int SolveBiCGSTAB( SparseMatrix *A, double *b, double *x, int maxit, double minerr);
void SolveGaussSeidelMethod( SparseMatrix *A, float *b, float *x0, int maxit, float minerr);
void SolveGaussSeidelMethod( SparseMatrix *A, double *b, double *x0, int maxit, double minerr);
void SolveGaussSeidelMethod( SparseMatrix *A, double *b, double *x0, int maxit, double minerr);
void SolveSuccessiveOverRelaxationMethod( SparseMatrix *A, double *b, double *x0, int maxit, double minerr, double omega);
void SolveJacobiMethod( SparseMatrix *A, double *b, double *x0, int maxit, double minerr);
float PreconditionedConjugateGradientSparse( SparseMatrix *A, SparseMatrix *M, float *b, float *x, int maxit, float minerr);
double PreconditionedConjugateGradientSparse( SparseMatrix *A, SparseMatrix *M, double *b, double *x, int maxit, float minerr);
void ConjugateGradientSparse( SparseMatrix *A, float *b, float *x, int maxit, float minerr);
void ConjugateGradientSparse( SparseMatrix *A, double *b, double *x, int maxit, double minerr);

void PreconditionJacobi( SparseMatrix *A, double *b);
void PreconditionSOR( SparseMatrix *A, double *b, double omega);
void SolveDirectly( SparseMatrix *Ai, double *bi, double *x, VoronoiDiagram *vd);

float GetLactateProductionRate( VoronoiCell *cell, float glucose, float oxygen);
float GetGlucoseConsumptionRate( VoronoiCell *cell, float glucose, float oxygen);
float GetOxygenConsumptionRate( VoronoiCell *cell, float glucose, float oxygen);
float OxygenConsumptionRate_PartialDerivativeOxygen( VoronoiCell *cell, float glucose, float oxygen);
float OxygenConsumptionRate_PartialDerivativeGlucose( VoronoiCell *cell, float glucose, float oxygen);
float GlucoseConsumptionRate_PartialDerivativeGlucose( VoronoiCell *cell, float glucose, float oxygen);
float GlucoseConsumptionRate_PartialDerivativeOxygen( VoronoiCell *cell, float glucose, float oxygen);


#endif
