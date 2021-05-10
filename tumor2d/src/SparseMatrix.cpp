#include "SparseMatrix.h"
#include "Mathematix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <omp.h>


#include "Substrate.h"
#include "Agent.h"
#define DX 13.

#define ROW_EXTENSION_SIZE 10


#define NO_NEGATIVES
//#define INIT_OXY 0.07
//#define INIT_GLU 0.8

// symmetric properties => CG
// FREE == LIQUID
#define FREE_DIFFUSION_SPEEDUP 30.//30.//15.
// FREE == TISSUE
//#define FREE_DIFFUSION_SPEEDUP 1.
#define PRECONDITIONED

// assymmetric properties => Newton
//#define FREE_IS_BOUNDARY
#define BORDER_IS_BOUNDARY


//#define CRANK_NICHOLSON
#define STEADY_STATE

#define PARTLY_LINEARIZED

#define DONT_ACTUALIZE_MATRIX 0
#define ACTUALIZE_MATRIX      1

#define DONT_ACTUALIZE_VECTOR 0
#define ACTUALIZE_VECTOR      1


#define CG_CRITICAL_ERROR		1e-5//0.5
#define NEWTON_CRITICAL_ERROR	0.5

#define COUPLING


//NOTES
// STEADY_STATE best
// PRECONDITIONED works fine
// lysis is bad
// migration is bad
// value reset to INIT_OXY and INIT_GLU worked fine


SparseMatrix *sA;
SparseMatrix *M;
SparseMatrix *J;

extern float *b;//[20];
extern float *x;//[20];
extern double *exactSolution;
 float *v0;
 float *v1;
 float *v2;
 float *v3;
 float *v4;
 float *v5;
 float *v6;
 float *v7;
 float *v8;
 double *v0d;
 double *v1d;
 double *v2d;
 double *v3d;
 double *v4d;
 double *v5d;
 double *v6d;

//typedef struct _SparseMatrix SparseMatrix;



SparseMatrix* SparseMatrix::newSparseMatrix( int dimI, int dimJ)
{
	fprintf( stderr, "SparseMatrix::newSparseMatrix()\n");

	SparseMatrix *newMatrix = (SparseMatrix*) malloc( sizeof(SparseMatrix));

	newMatrix->A        = (float**) malloc( dimI * sizeof(float*)); // values
	if( newMatrix->A == NULL) memoryAllocationError((char*) "newMatrix->A");
	newMatrix->JA       = (int**)   malloc( dimI * sizeof(int*));  // column index of values	
	if( newMatrix->A == NULL) memoryAllocationError((char*) "newMatrix->JA");
	newMatrix->sizeA    = (int*)    malloc( dimI * sizeof(int));
	if( newMatrix->A == NULL) memoryAllocationError((char*) "newMatrix->sizeA");
	newMatrix->maxSizeA = (int*)    malloc( dimI * sizeof(int));
	if( newMatrix->A == NULL) memoryAllocationError((char*) "newMatrix->maxSizeA");

	for( int i=0; i<dimI; i++){
		newMatrix->A[i]        = (float*) malloc( ROW_EXTENSION_SIZE * sizeof(float)); // values
		if( newMatrix->A[i] == NULL) memoryAllocationError((char*) "newMatrix->A[i]");
		newMatrix->JA[i]       = (int*)   malloc( ROW_EXTENSION_SIZE * sizeof(int));  // column index of values	
		if( newMatrix->JA[i] == NULL) memoryAllocationError((char*) "newMatrix->JA[i]");
		newMatrix->sizeA[i]    = 0;
		newMatrix->maxSizeA[i] = ROW_EXTENSION_SIZE;
	}

	newMatrix->dimI = dimI;
	newMatrix->dimJ = dimJ;
	
	return newMatrix;
}


void SparseMatrix::deleteSparseMatrix( SparseMatrix* matrix)
{
	for( int i=0; i<matrix->dimI; i++){
		free( matrix->A[i]);
		free( matrix->JA[i]);
	}

	free( matrix->A);
	free( matrix->JA);
	free( matrix->sizeA);
	free( matrix->maxSizeA);
	
	free( matrix);
}

float SparseMatrix::get( int i, int j)
{
	for( int jj=0; jj<this->sizeA[i]; jj++){
		if( JA[i][jj] == j){
			return A[i][jj];	
		}else if( JA[i][jj] > j){
			return 0.;
		}
	}
	
	return 0.;
}



void SparseMatrix::resetRow( int i)
{
	this->sizeA[i] = 0;
}


void SparseMatrix::setLast( int i, int j, float value)
{
	if( value == 0.) return;
	
	int k = this->sizeA[i];
	if( k > 0 && this->JA[i][k-1] >= j){
		fprintf( stderr, "ERROR in SparseMatrix::setLast(): Element cannot be inserted at last position!! \n(i,j) = (%i,%i), sizeA(i) = %i \n", i, j, k);
		exit( 0);	
	}
	
	// allocate memory
	if( this->sizeA[i] == this->maxSizeA[i]){
		this->maxSizeA[i] += ROW_EXTENSION_SIZE;
		this->JA[i] = (int*)   realloc( this->JA[i], this->maxSizeA[i] * sizeof(int));
		if( this->JA[i] == NULL) memoryAllocationError((char*) "this->JA[i]");
			this->A[i]  = (float*) realloc( this->A[i],  this->maxSizeA[i] * sizeof(float));
		if( this->A[i] == NULL) memoryAllocationError((char*) "this->A[i]");
	}
	
	// add element
	this->JA[i][k] = j;
	this->A[i][k] = value;
	this->sizeA[i] ++;		
	
	return;
	
}


void SparseMatrix::set( int i, int j, float value)
{
	if( value == 0.) return;
	
	int ji = 0;
	int jj = this->sizeA[i]-1;
	int k  = 0;
	// look for position
	if( this->sizeA[i]>0){
		ji = 0;
		int x = j;

		do{
			k = (ji + jj) / 2;
			if( x > this->JA[i][k])
				ji = k + 1;
			else 
				jj = k - 1;
		}while( (this->JA[i][k] != x) && (ji <= jj));
		// jj > ji => 
	}

	if( this->sizeA[i]==0 || (this->JA[i][k] != j && jj<ji)){
		if( jj<ji)
			k = ji;
		// realloc memory
		if( this->sizeA[i] == this->maxSizeA[i]){
			this->maxSizeA[i] += ROW_EXTENSION_SIZE;
			this->JA[i] = (int*)   realloc( this->JA[i], this->maxSizeA[i] * sizeof(int));
			if( this->JA[i] == NULL) memoryAllocationError((char*) "this->JA[i]");
			this->A[i]  = (float*) realloc( this->A[i],  this->maxSizeA[i] * sizeof(float));
			if( this->A[i] == NULL) memoryAllocationError((char*) "this->A[i]");
		}
		// shift following elements
		for( int jjj=this->sizeA[i]; jjj>k; jjj--){
			this->JA[i][jjj] = this->JA[i][jjj-1];
			this->A[i][jjj] = this->A[i][jjj-1];
		}
		this->JA[i][k] = j;
		this->A[i][k] = value;
		this->sizeA[i] ++;		

	}else{
			A[i][k] = value;
	}
	return;
}


void SparseMatrix::add( int i, int j, float value)
{
	if( value == 0.) return;

	int ji = 0;
	int jj = this->sizeA[i]-1;
	int k  = 0;
	// look for position
	if( this->sizeA[i]>0){
		ji = 0;
		int x = j;

		do{
			k = (ji + jj) / 2;
			if( x > this->JA[i][k])
				ji = k + 1;
			else
				jj = k - 1;
		}while( (this->JA[i][k] != x) && (ji <= jj));
	}


	if( this->sizeA[i]==0 || (this->JA[i][k] != j && jj<ji)){
		if( jj<ji)
			k = ji;
		// realloc memory
		if( this->sizeA[i] == this->maxSizeA[i]){
			this->maxSizeA[i] += ROW_EXTENSION_SIZE;
			this->JA[i] = (int*)   realloc( this->JA[i], this->maxSizeA[i] * sizeof(int));
			if( this->JA[i] == NULL) memoryAllocationError((char*) "this->JA[i]");
			this->A[i]  = (float*) realloc( this->A[i],  this->maxSizeA[i] * sizeof(float));
			if( this->A[i] == NULL) memoryAllocationError((char*) "this->A[i]");
		}
		// shift following elements
		for( int jjj=this->sizeA[i]; jjj>k; jjj--){
			this->JA[i][jjj] = this->JA[i][jjj-1];
			this->A[i][jjj] = this->A[i][jjj-1];
		}
		this->JA[i][k] = j;
		this->A[i][k] = value; // set first time
		this->sizeA[i] ++;

	}else{
			A[i][k] += value; // add to old value
	}
	return;
}



void setupMatrixImplicitSparse( VoronoiDiagram *voronoiDiagram, double timeStep, char molecule)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	
	float r = 0.;
	switch( molecule){
		case 'G':
			r = Glucose_Diffusion * timeStep/(DX*DX);
			break;
		case 'O':
			r = Oxygen_Diffusion * timeStep/(DX*DX);
			break;
	}
	
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		// actual element
		int m = i*di + ii*dii + iii*diii;

		// element inside domain
		if( TRUE
		#ifdef BORDER_IS_BOUNDARY
		    && i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		#endif
		#ifdef FREE_IS_BOUNDARY
		    && voronoiDiagram->voronoiCells[m]->getState() != FREE
		#endif
		)
		{
			//matrix
			
			// DIRICHLET
			
		#ifdef BORDER_IS_BOUNDARY
			#ifdef FREE_DIFFUSION_SPEEDUP
			int a = 0;
			
			if( voronoiDiagram->voronoiCells[m]->getState() == FREE && voronoiDiagram->voronoiCells[m-diii]->getState() == FREE)
			{sA->set(m, m-diii, -r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-diii, -r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[m]->getState() == FREE && voronoiDiagram->voronoiCells[m-dii]->getState() == FREE)
			{sA->set(m, m-dii, -r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-dii, -r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[m]->getState() == FREE && voronoiDiagram->voronoiCells[m-di]->getState() == FREE)
			{sA->set(m, m-di, -r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-di, -r); 

			#ifndef FREE_DIFFUSION_SPEEDUP
			switch( molecule){
				case 'G': sA->set(m, m, 1 + 6*r); break;
				case 'O': sA->set(m, m, 1 + 6*r);  break;
			}
			#endif

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[m]->getState() == FREE && voronoiDiagram->voronoiCells[m+di]->getState() == FREE)
			{sA->set(m, m+di, -r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+di, -r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[m]->getState() == FREE && voronoiDiagram->voronoiCells[m+dii]->getState() == FREE)
			{sA->set(m, m+dii, -r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+dii, -r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[m]->getState() == FREE && voronoiDiagram->voronoiCells[m+diii]->getState() == FREE)
			{sA->set(m, m+diii, -r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+diii, -r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			switch( molecule){
				case 'G': sA->set(m, m, 1 + (6-a)*r + a*r*FREE_DIFFUSION_SPEEDUP); break;
				case 'O': sA->set(m, m, 1 + (6-a)*r + a*r*FREE_DIFFUSION_SPEEDUP);  break;
			}
			#endif
		#else	

			// NEUMANN

			float diagValue = 1;	
	
			if(iii>0){ sA->set(m, m-diii, -r); diagValue += r;}
			if(ii>0){ sA->set(m, m-dii, -r); diagValue += r;}
			if( i>0){ sA->set(m, m-di, -r); diagValue += r;}
			if(i<voronoiDiagram->xN[0]-1){ sA->set(m, m+di, -r); diagValue += r;}
			if(ii<voronoiDiagram->xN[1]-1){ sA->set(m, m+dii, -r); diagValue += r;}
			if(iii<voronoiDiagram->xN[2]-1){ sA->set(m, m+diii, -r); diagValue += r;}
			sA->set(m, m, diagValue);	
		#endif
			// vector
			switch( molecule){
				case 'G':
					b[m] = voronoiDiagram->voronoiCells[m]->glucose * ( 1. - timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m]));
					break;
				case 'O':
					b[m] = voronoiDiagram->voronoiCells[m]->oxygen * ( 1. - timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m]));
					break;
			}

		}
		
		// border condition
		else{
			//matrix
			sA->set(m, m, 1);

			// vector
			switch( molecule){
				case 'G':
					b[m] = voronoiDiagram->voronoiCells[m]->glucose;
					break;
				case 'O':
					b[m] = voronoiDiagram->voronoiCells[m]->oxygen;
					break;
			}			
		}
		
		
	}
}


void setupMatrixImplicitCrankNicholsonSparse( VoronoiDiagram *voronoiDiagram, double timeStep, char molecule)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	
	//float a = 0.5;
	//float b = 1. - a;
	//float c = 1.
	
	float r = 0.;
	//float r2;
	 
	switch( molecule){
		case 'G':
			r = Glucose_Diffusion * timeStep/(DX*DX);
			break;
		case 'O':
			r = Oxygen_Diffusion * timeStep/(DX*DX);
			break;
	}
	//r2 = r * 10.;
	
	//int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		// actual element
		int m = i*di + ii*dii + iii*diii;

		// element inside domain
		if( i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		#ifdef FREE_IS_BOUNDARY
		    && voronoiDiagram->voronoiCells[m]->getState() != FREE
		#endif
		){
			//matrix
			sA->set(m, m-diii, -r); 
			sA->set(m, m-dii, -r); 
			sA->set(m, m-di, -r); 
			sA->set(m, m, 2 + 6*r);			
			sA->set(m, m+di, -r); 
			sA->set(m, m+dii, -r); 
			sA->set(m, m+diii, -r); 
			
			// vector
			float c = 0.5;
			float nc = 0.;
			switch( molecule){
				case 'G':
					b[m] = voronoiDiagram->voronoiCells[m]->glucose * ( 2. - 6.*r - timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m]))
					       + voronoiDiagram->voronoiCells[m-diii]->glucose*r 
					       + voronoiDiagram->voronoiCells[m-dii]->glucose*r 
					       + voronoiDiagram->voronoiCells[m-di]->glucose*r 
					       + voronoiDiagram->voronoiCells[m+di]->glucose*r 
					       + voronoiDiagram->voronoiCells[m+dii]->glucose*r 
					       + voronoiDiagram->voronoiCells[m+diii]->glucose*r 
					       + nc*c*(- voronoiDiagram->voronoiCells[m-diii]->glucose*timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m-diii])
					       - voronoiDiagram->voronoiCells[m-dii]->glucose*timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m-dii])
					       - voronoiDiagram->voronoiCells[m-di]->glucose*timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m-di])
					       - voronoiDiagram->voronoiCells[m-di]->glucose*timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m+di])
					       - voronoiDiagram->voronoiCells[m-dii]->glucose*timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m+dii])
					       - voronoiDiagram->voronoiCells[m-diii]->glucose*timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m+diii]))
					       ;
					break;
				case 'O':
					b[m] = voronoiDiagram->voronoiCells[m]->oxygen * ( 2. - 2.*r - timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m]))
					       + voronoiDiagram->voronoiCells[m-diii]->oxygen*r 
					       + voronoiDiagram->voronoiCells[m-dii]->oxygen*r 
					       + voronoiDiagram->voronoiCells[m-di]->oxygen*r 
					       + voronoiDiagram->voronoiCells[m+di]->oxygen*r 
					       + voronoiDiagram->voronoiCells[m+dii]->oxygen*r 
					       + voronoiDiagram->voronoiCells[m+diii]->oxygen*r 
					       + nc*c*(- voronoiDiagram->voronoiCells[m-diii]->oxygen*timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m-diii])
					       - voronoiDiagram->voronoiCells[m-dii]->oxygen*timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m-dii])
					       - voronoiDiagram->voronoiCells[m-di]->oxygen*timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m-di])
					       - voronoiDiagram->voronoiCells[m-di]->oxygen*timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m+di])
					       - voronoiDiagram->voronoiCells[m-dii]->oxygen*timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m+dii])
					       - voronoiDiagram->voronoiCells[m-diii]->oxygen*timeStep * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m+diii]))
					       ;
					break;
			}
		}
		
		// border condition
		else{
			//matrix
			sA->set(m, m, 1);

			// vector
			switch( molecule){
				case 'G':
					b[m] = voronoiDiagram->voronoiCells[m]->glucose;
					break;
				case 'O':
					b[m] = voronoiDiagram->voronoiCells[m]->oxygen;
					break;
			}			
		}
		
		
	}
}

void setupMatrixImplicitSteadyState( VoronoiDiagram *voronoiDiagram, char molecule)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];


	float r;

	switch( molecule){
		case 'G':
			r = Glucose_Diffusion /(DX*DX);
			break;
		case 'O':
			r = Oxygen_Diffusion /(DX*DX);
			break;
	}
	//r2 = r * 10.;
	
	//int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	int mi=0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		if( mi==0) {molecule ='G'; r = Glucose_Diffusion /(DX*DX);}
		if( mi==1) {molecule ='O'; r = Oxygen_Diffusion /(DX*DX); }
		
		// actual element of the matrix
		int m = i*di + ii*dii + iii*diii  + mi*d;

		// actual element of the voronoi diagram
		int v = i*di + ii*dii + iii*diii;

		// reset row: set all entries to zero
		sA->resetRow( m);

		// element inside domain
		if( 
			voronoiDiagram->voronoiCells[v]->getState() != VESSEL //TRUE
		#ifdef BORDER_IS_BOUNDARY
		    && i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		#endif
		#ifdef FREE_IS_BOUNDARY
		    //&& voronoiDiagram->voronoiCells[v]->getState() != FREE
		    && voronoiDiagram->voronoiCells[v]->agent != NULL
		#endif
		)
		{
			//matrix
			
			// DIRICHLET
			
		#ifdef BORDER_IS_BOUNDARY
			#ifdef FREE_DIFFUSION_SPEEDUP
			int a = 0;
			
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v-diii]->getState() == FREE)
			{sA->set(m, m-diii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-diii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v-dii]->getState() == FREE)
			{sA->set(m, m-dii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-dii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v-di]->getState() == FREE)
			{sA->set(m, m-di, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-di, r); 

			#ifndef FREE_DIFFUSION_SPEEDUP
			switch( molecule){
				case 'G': sA->set(m, m, - 6*r - GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v])); break;
				case 'O': sA->set(m, m, - 6*r - GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));  break;
			}
			#endif

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v+di]->getState() == FREE)
			{sA->set(m, m+di, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+di, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v+dii]->getState() == FREE)
			{sA->set(m, m+dii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+dii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v+diii]->getState() == FREE)
			{sA->set(m, m+diii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+diii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			switch( molecule){
				case 'G': sA->set(m, m, - (6-a)*r - a*r*FREE_DIFFUSION_SPEEDUP - GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v])); break;
				case 'O': sA->set(m, m, - (6-a)*r - a*r*FREE_DIFFUSION_SPEEDUP - GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));  break;
			}
			#endif
		#else	

			// NEUMANN

			float diagValue = 0.;	
			switch( molecule){
				case 'G': sA->set(m, m, diagValue = -GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v])); break;
				case 'O': sA->set(m, m, diagValue = -GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));  break;
			}
	
			if(iii>0){ sA->set(m, m-diii, r); diagValue -= r;}
			if(ii>0){ sA->set(m, m-dii, r); diagValue -= r;}
			if( i>0){ sA->set(m, m-di, r); diagValue -= r;}
			if(i<voronoiDiagram->xN[0]-1){ sA->set(m, m+di, r); diagValue -= r;}
			if(ii<voronoiDiagram->xN[1]-1){ sA->set(m, m+dii, r); diagValue -= r;}
			if(iii<voronoiDiagram->xN[2]-1){ sA->set(m, m+diii, r); diagValue -= r;}
			sA->set(m, m, diagValue);	
		#endif
			// vector
			switch( molecule){
				case 'G':
					b[m] = 0.;
					break;
				case 'O':
					b[m] = 0.;
					break;
			}

		}
		
		// border condition
		else{
			//matrix
			sA->set(m, m, 1);

			// vector
			switch( molecule){
				case 'G':
					b[m] = voronoiDiagram->voronoiCells[v]->glucose;
					break;
				case 'O':
					b[m] = voronoiDiagram->voronoiCells[v]->oxygen;
					break;
			}			
		}
		

	}
}

void setupGrowthFactorMatrixImplicit( VoronoiDiagram *voronoiDiagram, SparseMatrix *A, float *b, char actualize_vector, char actualize_matrix, float timeStep, float spaceStep)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	float r = GrowthFactors_Diffusion * timeStep/(spaceStep*spaceStep);
	
	int m = 0;
	int v = 0;

	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		
		// reset row: set all entries to zero
		if( actualize_matrix)
		sA->resetRow( m);

		// element inside domain
		if( 
		    voronoiDiagram->voronoiCells[v]->getState() != NECROTIC &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			if( actualize_matrix){
			//matrix
#ifdef FREE_DIFFUSION_SPEEDUP
			float a = 0.;
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-diii]->getState()==FREE)
			{sA->set(m, m-diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-diii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-dii]->getState()==FREE)
			{sA->set(m, m-dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-di]->getState()==FREE)
			{sA->set(m, m-di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+di]->getState()==FREE)
			{sA->set(m, m+di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+dii]->getState()==FREE)
			{sA->set(m, m+dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+diii]->getState()==FREE)
			{sA->set(m, m+diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+diii, 1);

			sA->set(m, m, - (a*FREE_DIFFUSION_SPEEDUP + 6-a));	
#else
			sA->set(m, m-diii, 1);
			sA->set(m, m-dii, 1);
			sA->set(m, m-di, 1);

			sA->set(m, m, - 6 - 1/r);	

			sA->set(m, m+di, 1);
			sA->set(m, m+dii, 1);
			sA->set(m, m+diii, 1);
#endif
			}

			// vector
			if( actualize_vector)
			b[m] = -voronoiDiagram->voronoiCells[v]->growthfactors/r;
		
		}
		
		// border condition
		else{
			//matrix
			if( actualize_matrix)
			sA->set(m, m, 1);

			// vector
			if( actualize_vector){
				if( voronoiDiagram->voronoiCells[v]->getState() == NECROTIC)
					b[m] = 1; //mM
				else
					b[m] = voronoiDiagram->voronoiCells[v]->growthfactors;
			}
		}
		m++;
		v = (v + 1) % d;
	}
}





void setupTestMatrixSteadyStateNewton( VoronoiDiagram *voronoiDiagram, SparseMatrix * sA, float *x, char actualize_vector, char actualize_matrix)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	//int nx = voronoiDiagram->xN[0],
	//    ny = voronoiDiagram->xN[1],
	//    nz = voronoiDiagram->xN[2];
	
	float r,
	      //rg = 1, //Glucose_Diffusion /(DX*DX), //1
	      //ro = 10; //Oxygen_Diffusion /(DX*DX);
	      rg = Glucose_Diffusion /(DX*DX), //1
	      ro = Oxygen_Diffusion /(DX*DX);
	//float r = 1.;
	
	int mi=0;
	int m = 0;
	int v = 0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		


		// reset row: set all entries to zero
		if( actualize_matrix)
		sA->resetRow( m);

		// element inside domain
		if( 
			//voronoiDiagram->voronoiCells[v]->getState() != FREE &&
		    voronoiDiagram->voronoiCells[v]->getState() != VESSEL &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			float consumption, diagValue = 0., linkingValue = 0.;
			
			// consumption
			consumption = ( mi==0 ? GiveMeTheGlucoseRate(voronoiDiagram->voronoiCells[v], x[m], x[m+d])
			                      : GiveMeTheOxygenRate( voronoiDiagram->voronoiCells[v], x[m-d], x[m]));
			if( voronoiDiagram->voronoiCells[v]->getState() != COMPARTMENT)
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;
			
			// partial derivative of consumption
			// I. Glucose, II. Oxygen
			switch( mi){
				case 0: diagValue    = GlucoseConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				        linkingValue = GlucoseConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				        //J->set(m, m+d, linkingValue/r);
				        break;
				case 1: diagValue    = OxygenConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        linkingValue = OxygenConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        //J->set(m, m-d, linkingValue/r);
				        break;
			}
			
			r = (mi==0 ? rg : ro);
			if( actualize_matrix){
			//matrix
						
			if( mi == 1) sA->setLast(m, m-d, -linkingValue/r);
			
#ifdef FREE_DIFFUSION_SPEEDUP
			float a = 0.;
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-diii]->getState()==FREE)
			{sA->setLast(m, m-diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-diii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-dii]->getState()==FREE)
			{sA->setLast(m, m-dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-di]->getState()==FREE)
			{sA->setLast(m, m-di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+di]->getState()==FREE)
			{sA->setLast(m, m+di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+dii]->getState()==FREE)
			{sA->setLast(m, m+dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+diii]->getState()==FREE)
			{sA->setLast(m, m+diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+diii, 1);

			float consumption = ( mi==0 ? GiveMeTheGlucoseRate(voronoiDiagram->voronoiCells[v])
				                        : GiveMeTheOxygenRate( voronoiDiagram->voronoiCells[v]));
			if( voronoiDiagram->voronoiCells[v]->getState() != COMPARTMENT)
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;
			sA->set(m, m, - (a*FREE_DIFFUSION_SPEEDUP + 6-a) - consumption/r);	
#else
			sA->setLast(m, m-diii, 1);
			sA->setLast(m, m-dii, 1);
			sA->setLast(m, m-di, 1);

			// partly linearized
			sA->setLast(m, m, - 6 - diagValue/r);	

			sA->setLast(m, m+di, 1);
			sA->setLast(m, m+dii, 1);
			sA->setLast(m, m+diii, 1);
#endif

			if( mi == 0) sA->setLast(m, m+d, -linkingValue/r);

			}

			// vector
			// partly linearized
			if( actualize_vector)
			b[m] = consumption/r - diagValue/r*x[m] - linkingValue/r*( mi==0 ? x[m+d] : x[m+d]);
		
		}
		
		// border condition
		else{
			//matrix
			if( actualize_matrix)
			sA->setLast(m, m, 1);

			// vector
			if( actualize_vector)
			b[m] = (mi==0 ? voronoiDiagram->voronoiCells[v]->glucose
			              : voronoiDiagram->voronoiCells[v]->oxygen);
		}
		m++;
		v = (v + 1) % d;
	}
}









/*void setupBloodVesselMatrixSteadyState( VoronoiDiagram *voronoiDiagram, SparseMatrix * sA, float *x, char actualize_vector, char actualize_matrix)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	float r,
	      rbv = Bloodvessel_Diffusion /(DX*DX);
	
	int mi=0;
	int m = 0;
	int v = 0;
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		// reset row: set all entries to zero
		if( actualize_matrix)
		sA->resetRow( m);

		// element inside domain
		if( 
		    voronoiDiagram->voronoiCells[v]->getState() != VESSEL &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			float consumption, diagValue, linkingValue;
			
			// consumption
			consumption = ( mi==0 ? GiveMeTheGlucoseRate(voronoiDiagram->voronoiCells[v], x[m], x[m+d])
			                      : GiveMeTheOxygenRate( voronoiDiagram->voronoiCells[v], x[m-d], x[m]));
			if( voronoiDiagram->voronoiCells[v]->getState() != COMPARTMENT){
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;
			
			// partial derivative of consumption
			// I. Glucose, II. Oxygen
			
			r = (mi==0 ? rg : ro);
			if( actualize_matrix){
			//matrix
						
			sA->setLast(m, m-diii, 1);
			sA->setLast(m, m-dii, 1);
			sA->setLast(m, m-di, 1);

			// partly linearized
			sA->setLast(m, m, - 6 - diagValue/r);	

			sA->setLast(m, m+di, 1);
			sA->setLast(m, m+dii, 1);
			sA->setLast(m, m+diii, 1);

			if( mi == 0) sA->setLast(m, m+d, -linkingValue/r);

			}

			// vector
			// partly linearized
			if( actualize_vector)
			b[m] = consumption/r - diagValue/r*x[m] - linkingValue/r*( mi==0 ? x[m+d] : x[m+d]);
		
		}
		
		// border condition
		else{
			//matrix
			if( actualize_matrix)
			sA->setLast(m, m, 1);

			// vector
			if( actualize_vector)
			b[m] = (mi==0 ? voronoiDiagram->voronoiCells[v]->glucose
			              : voronoiDiagram->voronoiCells[v]->oxygen);
		}
		m++;
		v = (v + 1) % d;
	}
}*/




void setupTestMatrixSteadyState( VoronoiDiagram *voronoiDiagram, SparseMatrix *sA, float* x, char actualize_vector, char actualize_matrix)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]
#if DIMENSIONS >= 2
	        *voronoiDiagram->xN[1]
#endif
#if DIMENSIONS >= 3
*voronoiDiagram->xN[2]
#endif
	        ;


	//int nx = voronoiDiagram->xN[0],
	//    ny = voronoiDiagram->xN[1],
	//    nz = voronoiDiagram->xN[2];

	float r,
	      //rg = 1, //Glucose_Diffusion /(DX*DX), //1
	      //ro = 10; //Oxygen_Diffusion /(DX*DX);
	      rg = Glucose_Diffusion /(DX*DX), //1
	      ro = Oxygen_Diffusion /(DX*DX);
	//float r = 1.;

	int mi=0;
	int m = 0;
	int v = 0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
#if DIMENSIONS >= 3
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
#endif
#if DIMENSIONS >= 2
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
#endif
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{

		// reset row: set all entries to zero
		if( actualize_matrix)
		sA->resetRow( m);

		// element inside domain
		if(
			//voronoiDiagram->voronoiCells[v]->getState() != FREE &&
		    voronoiDiagram->voronoiCells[v]->getState() != VESSEL &&
		    i>0 && i<voronoiDiagram->xN[0]-1
#if DIMENSIONS >= 2
		    && ii>0 && ii<voronoiDiagram->xN[1]-1
#endif
#if DIMENSIONS >= 3
		    && iii>0 && iii<voronoiDiagram->xN[2]-1
#endif
		)
		{

			float consumption = ( mi==0 ? GiveMeTheGlucoseRate(voronoiDiagram->voronoiCells[v], x[m], x[m+d])
				                        : GiveMeTheOxygenRate( voronoiDiagram->voronoiCells[v], x[m-d], x[m]));
			if( voronoiDiagram->voronoiCells[v]->getState() != COMPARTMENT){
				if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE)
					consumption = 0.;
			}

			r = (mi==0 ? rg : ro);
			if( actualize_matrix){
			//matrix
#ifdef FREE_DIFFUSION_SPEEDUP
			float a = 0.;

#if DIMENSIONS == 3
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-diii]->getState()==FREE)
			{sA->setLast(m, m-diii, FREE_DIFFUSION_SPEEDUP); a+=FREE_DIFFUSION_SPEEDUP;} else{ sA->setLast(m, m-diii, 1); a+=1;}
#endif
#if DIMENSIONS >= 2
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-dii]->getState()==FREE)
			{sA->setLast(m, m-dii,  FREE_DIFFUSION_SPEEDUP); a+=FREE_DIFFUSION_SPEEDUP;} else{ sA->setLast(m, m-dii, 1); a+=1;}
#endif
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-di]->getState()==FREE)
			{sA->setLast(m, m-di,   FREE_DIFFUSION_SPEEDUP); a+=FREE_DIFFUSION_SPEEDUP;} else{ sA->setLast(m, m-di, 1); a+=1;}

			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+di]->getState()==FREE)
			{sA->setLast(m, m+di,   FREE_DIFFUSION_SPEEDUP); a+=FREE_DIFFUSION_SPEEDUP;} else{ sA->setLast(m, m+di, 1); a+=1;}
#if DIMENSIONS >= 2
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+dii]->getState()==FREE)
			{sA->setLast(m, m+dii,  FREE_DIFFUSION_SPEEDUP); a+=FREE_DIFFUSION_SPEEDUP;} else{ sA->setLast(m, m+dii, 1); a+=1;}
#endif
#if DIMENSIONS == 3
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+diii]->getState()==FREE)
			{sA->setLast(m, m+diii, FREE_DIFFUSION_SPEEDUP); a+=FREE_DIFFUSION_SPEEDUP;} else{ sA->setLast(m, m+diii, 1); a+=1;}
#endif

			// partly linearized
#ifdef PARTLY_LINEARIZED

			float diag = -a - consumption/(r*x[m]);
			if( isnan( -a - consumption/(r*x[m]))){
				diag = -a;
			}
			sA->set(m, m, diag);


#else

			// fully linearized
			sA->setLast(m, m, - a);
#endif


#else
	#if DIMENSIONS == 3
			sA->setLast(m, m-diii, 1);
	#endif
	#if DIMENSIONS >= 2
			sA->setLast(m, m-dii, 1);
	#endif
			sA->setLast(m, m-di, 1);

			// partly linearized
	#ifdef PARTLY_LINEARIZED

			float diag = -2*DIMENSIONS - consumption/(r*x[m]);
			if( isnan( - 2*DIMENSIONS - consumption/(r*x[m]))){
				diag = -2.*DIMENSIONS;
			}
			sA->setLast(m, m, diag);


	#else

			// fully linearized
			sA->setLast(m, m, - 6);
	#endif

			sA->setLast(m, m+di, 1);
	#if DIMENSIONS >= 2
			sA->setLast(m, m+dii, 1);
	#endif
	#if DIMENSIONS == 3
			sA->setLast(m, m+diii, 1);
	#endif
#endif
		}

			// vector
			// partly linearized
	#ifdef PARTLY_LINEARIZED
			if( actualize_vector)
			b[m] = 0.;

	#else
			// fully linearized
			float consumption = ( mi==0 ? GiveMeTheGlucoseRate(voronoiDiagram->voronoiCells[v])
				                        : GiveMeTheOxygenRate( voronoiDiagram->voronoiCells[v]));
			if( voronoiDiagram->voronoiCells[v]->getState() != COMPARTMENT)
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;

			if( actualize_vector)
			b[m] = ( mi==0 ? voronoiDiagram->voronoiCells[v]->glucose:voronoiDiagram->voronoiCells[v]->oxygen)*consumption/r;
	#endif
		}

		// border condition
		else{
			//fprintf( stderr, "BORDER!!!\n");
			//matrix
			if( actualize_matrix)
			sA->setLast(m, m, 1);

			// vector
			if( actualize_vector)
			b[m] = (mi==0 ? voronoiDiagram->voronoiCells[v]->glucose
			              : voronoiDiagram->voronoiCells[v]->oxygen);
		}
		m++;
		v = (v + 1) % d;
	}

	for(int i=0; i<sA->dimI; i++){
		float aii = sA->get( i,i);
		for( int jj=0; jj<sA->sizeA[i]; jj++){
			//int j=sA->JA[i][jj];
			sA->A[i][jj] /= aii;
		}
		b[i] /= aii;
	}
}




void setupJacobiMatrixSteadyState( VoronoiDiagram *voronoiDiagram, SparseMatrix *J, float* x)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	//int nx = voronoiDiagram->xN[0],
	//    ny = voronoiDiagram->xN[1],
	//    nz = voronoiDiagram->xN[2];
	
	float r,
	      //rg = 1, //Glucose_Diffusion /(DX*DX), //1
	      //ro = 10; //Oxygen_Diffusion /(DX*DX);
	      rg = Glucose_Diffusion /(DX*DX), //1
	      ro = Oxygen_Diffusion /(DX*DX);
	//float r = 1.;
	
	int mi=0;
	int m = 0;
	int v = 0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		
		// actual element of the matrix
		//int m = i*di + ii*dii + iii*diii  + mi*d;

		// actual element of the voronoi diagram
		//int v = i*di + ii*dii + iii*diii;
		//fprintf( stderr, "m:%i, v:%i\n ", m, v);

		// reset row: set all entries to zero
		//if( actualize_matrix)
		J->resetRow( m);

		// element inside domain
		if( 
			//voronoiDiagram->voronoiCells[v]->getState() != FREE &&
		    voronoiDiagram->voronoiCells[v]->getState() != VESSEL &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			//float consumption =  (mi==0 ? GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v]) 
			//                            : GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));
			//int radius = 3;
			//float consumption;
			/*float oxy = voronoiDiagram->voronoiCells[v]->oxygen; // 0.07
			float glu = voronoiDiagram->voronoiCells[v]->glucose;
			float min = 0.001,
			      max = 0.09,
			      kgg = 0.01,
			      kgo = 0.01;
			float consumption = ( mi==0 ? glu/(kgg+glu)*(max - (max-min) * oxy/(kgo+oxy))
			                            : oxy/(kgg+oxy)*(max - (max-min) * glu/(kgo+glu)));*/
			
			/////////////////////
			//consumption
			float diagValue = 0.;
			float linkingValue = 0.;
			// I. Glucose, II. Oxygen
			switch( mi){
				case 0: diagValue    = GlucoseConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				        linkingValue = GlucoseConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				        //J->set(m, m+d, linkingValue/r);
				        break;
				case 1: diagValue    = OxygenConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        linkingValue = OxygenConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        //J->set(m, m-d, linkingValue/r);
				        break;
			}
			/////////
			//fprintf( stderr, "test1\n");
					
			r = (mi==0 ? rg : ro);		
			//if( !init_only_vector)
			{

			if( mi == 1) J->setLast(m, m-d, -linkingValue/r);
			//fprintf( stderr, "test2\n");

			//matrix
#ifdef FREE_DIFFUSION_SPEEDUP
			float a = 0.;
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-diii]->getState()==FREE)
			{J->setLast(m, m-diii, FREE_DIFFUSION_SPEEDUP); a++;} else J->setLast(m, m-diii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-dii]->getState()==FREE)
			{J->setLast(m, m-dii, FREE_DIFFUSION_SPEEDUP); a++;} else J->setLast(m, m-dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-di]->getState()==FREE)
			{J->setLast(m, m-di, FREE_DIFFUSION_SPEEDUP); a++;} else J->setLast(m, m-di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+di]->getState()==FREE)
			{J->setLast(m, m+di, FREE_DIFFUSION_SPEEDUP); a++;} else J->setLast(m, m+di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+dii]->getState()==FREE)
			{J->setLast(m, m+dii, FREE_DIFFUSION_SPEEDUP); a++;} else J->setLast(m, m+dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+diii]->getState()==FREE)
			{J->setLast(m, m+diii, FREE_DIFFUSION_SPEEDUP); a++;} else J->setLast(m, m+diii, 1);

			float consumption = diagValue;
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;
			J->set(m, m, - (a*FREE_DIFFUSION_SPEEDUP + 6-a) - consumption/r);	
#else
			//fprintf( stderr, "test3\n");
			J->setLast(m, m-diii, 1);
			J->setLast(m, m-dii, 1);
			J->setLast(m, m-di, 1);

			// partly linearized
	#ifdef PARTLY_LINEARIZED
			float consumption = diagValue;
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;
			J->setLast(m, m, - 6 - consumption/r);	
	#else			
			
			// fully linearized
			J->setLast(m, m, - 6);	
	#endif			

			//fprintf( stderr, "test4\n");
			J->setLast(m, m+di, 1);
			J->setLast(m, m+dii, 1);
			J->setLast(m, m+diii, 1);
#endif

			if( mi == 0) J->setLast(m, m+d, -linkingValue/r);
			
			//fprintf( stderr, "test\n");

			}

		}
		
		// border condition
		else{
			//matrix
			//if( !init_only_vector)
			J->setLast(m, m, 1);
		}
		m++;
		v = (v + 1) % d;
	}
}



void setupTestMatrixImplicit( VoronoiDiagram *voronoiDiagram, char init_only_vector, float timeStep, float spaceStep)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	int nx = voronoiDiagram->xN[0],
	    ny = voronoiDiagram->xN[1],
	    nz = voronoiDiagram->xN[2];
	
	float r,
	      rg = 1, //Glucose_Diffusion /(DX*DX), //1
	      ro = 10; //Oxygen_Diffusion /(DX*DX);
	//float r = 1.;
	
	int mi=0;
	int m = 0;
	int v = 0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		
		// actual element of the matrix
		//int m = i*di + ii*dii + iii*diii  + mi*d;

		// actual element of the voronoi diagram
		//int v = i*di + ii*dii + iii*diii;
		//fprintf( stderr, "m:%i, v:%i\n ", m, v);

		// reset row: set all entries to zero
		if( !init_only_vector)
		sA->resetRow( m);

		// element inside domain
		if( 
			//voronoiDiagram->voronoiCells[v]->getState() != FREE &&
		    voronoiDiagram->voronoiCells[v]->getState() != VESSEL &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			//float consumption =  (mi==0 ? GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v]) 
			//                            : GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));
			//int radius = 3;
			float oxy = voronoiDiagram->voronoiCells[v]->oxygen; // 0.07
			float glu = voronoiDiagram->voronoiCells[v]->glucose;
			float min = 0.001,
			      max = 0.09,
			      kgg = 0.01,
			      kgo = 0.01;
			
			if( !init_only_vector){
			r = (mi==0 ? rg : ro) * timeStep / (spaceStep*spaceStep);
			//matrix
			sA->set(m, m-diii, r);
			sA->set(m, m-dii, r);
			sA->set(m, m-di, r);
			float consumption = ( mi==0 ? 1/(kgg+glu)*(max - (max-min) * oxy/(kgo+oxy))
			                            : 1/(kgg+oxy)*(max - (max-min) * glu/(kgo+glu)));
			//sA->set(m, m, - ( voronoiDiagram->voronoiCells[v]->getState() != FREE ? consumption : 0.) - 6*r);	
			////sA->set(m, m, - (abs(i-nx/2)<3 && abs(ii-ny/2)<3 && abs(iii-nz/2)<3 ? consumption : 0.) - 6*r);	
		
			sA->set(m, m, - 6*r - 1/*spaceStep*spaceStep / timeStep*/ - (abs(i-nx/2)<3 && abs(ii-ny/2)<3 && abs(iii-nz/2)<3 ? consumption : 0.) * timeStep/*spaceStep*spaceStep*/);	
			sA->set(m, m+di, r);
			sA->set(m, m+dii, r);
			sA->set(m, m+diii, r);
			}

			// vector
			//b[m] = 0.;
		
			//float consumption = ( mi==0 ? glu/(kgg+glu)*(max - (max-min) * oxy/(kgo+oxy))
			//                            : oxy/(kgg+oxy)*(max - (max-min) * glu/(kgo+glu)));
			b[m] = - (mi==0 ? glu : oxy) /** spaceStep*spaceStep / timeStep*/
			       //+ (abs(i-nx/2)<radius && abs(ii-ny/2)<radius && abs(iii-nz/2)<radius ? consumption : 0.) * spaceStep*spaceStep
			       ;
			//b[m] = (voronoiDiagram->voronoiCells[v]->getState() != FREE ? consumption : 0.);
		}
		
		// border condition
		else{
			//matrix
			if( !init_only_vector)
			sA->set(m, m, 1);

			// vector
			b[m] = (mi==0 ? voronoiDiagram->voronoiCells[v]->glucose
			              : voronoiDiagram->voronoiCells[v]->oxygen);
		}
		m++;
		v = (v + 1) % d;
	}
}


void setupJacobiMatrixCrankNicholson( VoronoiDiagram *voronoiDiagram, SparseMatrix *J, float* x, float timeStep, float spaceStep, float beta)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	//int nx = voronoiDiagram->xN[0],
	//    ny = voronoiDiagram->xN[1],
	//    nz = voronoiDiagram->xN[2];
	
	float r,
	      //rg = 1, //Glucose_Diffusion /(DX*DX), //1
	      //ro = 10; //Oxygen_Diffusion /(DX*DX);
	      rg = Glucose_Diffusion /(DX*DX), //1
	      ro = Oxygen_Diffusion /(DX*DX);
	//float r = 1.;
	
	int mi=0;
	int m = 0;
	int v = 0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		
		// actual element of the matrix
		//int m = i*di + ii*dii + iii*diii  + mi*d;

		// actual element of the voronoi diagram
		//int v = i*di + ii*dii + iii*diii;
		//fprintf( stderr, "m:%i, v:%i\n ", m, v);

		// reset row: set all entries to zero
		//if( !init_only_vector)//{
		J->resetRow( m); //printf("reset row %i\n", m);}

		// element inside domain
		if( 
			//voronoiDiagram->voronoiCells[v]->getState() != FREE &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			//float consumption =  (mi==0 ? GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v]) 
			//                            : GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));
			//int radius = 3;
			/*float oxy = voronoiDiagram->voronoiCells[v]->oxygen; // 0.07
			float glu = voronoiDiagram->voronoiCells[v]->glucose;
			float min = 0.001,
			      max = 0.09,
			      kgg = 0.01,
			      kgo = 0.01;
			*/
			float consumption = 0.;
			//if( !init_only_vector)
			{
			
			/////////////////////
			//consumption
			float diagValue = 0.;
			float linkingValue = 0.;
			// I. Glucose, II. Oxygen
			switch( mi){
				case 0: diagValue    = GlucoseConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				        linkingValue = GlucoseConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				        //J->set(m, m+d, linkingValue/r);
				        break;
				case 1: diagValue    = OxygenConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        linkingValue = OxygenConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        //J->set(m, m-d, linkingValue/r);
				        break;
			}
			/////////
							
			// diffusion
			r = (mi==0 ? rg : ro) * timeStep / (spaceStep*spaceStep);

			// consumption
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) 
				consumption = 0.;
			else
				consumption = diagValue;


			//printf("test1\n");
			if( mi == 1) J->setLast(m, m-d, - linkingValue * timeStep/r);
			
			//printf("test2\n");
			
#ifdef FREE_DIFFUSION_SPEEDUP
			float a = 0.;

			//matrix
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-diii]->getState()==FREE)
			{J->setLast(m, m-diii, FREE_DIFFUSION_SPEEDUP); a++;} else J->set(m, m-diii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-dii]->getState()==FREE)
			{J->setLast(m, m-dii, FREE_DIFFUSION_SPEEDUP); a++;} else J->set(m, m-dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-di]->getState()==FREE)
			{J->setLast(m, m-di, FREE_DIFFUSION_SPEEDUP); a++;} else J->set(m, m-di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+di]->getState()==FREE)
			{J->setLast(m, m+di, FREE_DIFFUSION_SPEEDUP); a++;} else J->set(m, m+di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+dii]->getState()==FREE)
			{J->setLast(m, m+dii, FREE_DIFFUSION_SPEEDUP); a++;} else J->set(m, m+dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+diii]->getState()==FREE)
			{J->setLast(m, m+diii, FREE_DIFFUSION_SPEEDUP); a++;} else J->set(m, m+diii, 1);
			
			J->set(m, m, - (a*FREE_DIFFUSION_SPEEDUP + 6-a) - 1/r/beta - consumption/r * timeStep);	
#else
			J->setLast(m, m-diii, 1);
			J->setLast(m, m-dii, 1);
			J->setLast(m, m-di, 1);
			//float consumption = ( mi==0 ? 1/(kgg+glu)*(max - (max-min) * oxy/(kgo+oxy))
			//                            : 1/(kgg+oxy)*(max - (max-min) * glu/(kgo+glu)));
			//int tempState = voronoiDiagram->voronoiCells[v]->getState();
			J->setLast(m, m, - 6 - 1/r/beta - consumption * timeStep/r);	
		
			J->setLast(m, m+di, 1);
			J->setLast(m, m+dii, 1);
			J->setLast(m, m+diii, 1);
#endif

			//printf("test3\n");
			if( mi == 0) J->setLast(m, m+d, - linkingValue * timeStep/r);
			
			}
		}
		
		// border condition
		else{
			

			J->setLast(m, m, 1);
			
		}
		m++;
		v = (v + 1) % d;
	}
}


void setupTestMatrixCrankNicholson( VoronoiDiagram *voronoiDiagram, char actualize_vector, char actualize_matrix, float timeStep, float spaceStep, float beta)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	//int nx = voronoiDiagram->xN[0],
	//    ny = voronoiDiagram->xN[1],
	//    nz = voronoiDiagram->xN[2];
	
	float r = 0.,
	      //rg = 1, //Glucose_Diffusion /(DX*DX), //1
	      //ro = 10; //Oxygen_Diffusion /(DX*DX);
	      rg = Glucose_Diffusion, //1
	      ro = Oxygen_Diffusion;
	//float r = 1.;
	
	int mi=0;
	int m = 0;
	int v = 0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		
		// actual element of the matrix
		//int m = i*di + ii*dii + iii*diii  + mi*d;

		// actual element of the voronoi diagram
		//int v = i*di + ii*dii + iii*diii;
		//fprintf( stderr, "m:%i, v:%i\n ", m, v);

		// reset row: set all entries to zero
		if( actualize_matrix)
		sA->resetRow( m);

		// element inside domain
		if( 
			//voronoiDiagram->voronoiCells[v]->getState() != FREE &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			//float consumption =  (mi==0 ? GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v]) 
			//                            : GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));
			//int radius = 3;
			/*float oxy = voronoiDiagram->voronoiCells[v]->oxygen; // 0.07
			float glu = voronoiDiagram->voronoiCells[v]->glucose;
			float min = 0.001,
			      max = 0.09,
			      kgg = 0.01,
			      kgo = 0.01;
			*/
			float consumption = 0.;
			if( actualize_matrix)
			{
				
			// diffusion
			r = (mi==0 ? rg : ro) * timeStep / (spaceStep*spaceStep);

			// consumption
			consumption = ( mi==0 ? GiveMeTheGlucoseRate(voronoiDiagram->voronoiCells[v])
			                      : GiveMeTheOxygenRate( voronoiDiagram->voronoiCells[v]));
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;
			//consumption = (abs(i-nx/2)<3 && abs(ii-ny/2)<3 && abs(iii-nz/2)<3 ? consumption : 0.);

#ifdef FREE_DIFFUSION_SPEEDUP
			float a = 0.;

			//matrix
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-diii]->getState()==FREE)
			{sA->setLast(m, m-diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-diii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-dii]->getState()==FREE)
			{sA->setLast(m, m-dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-di]->getState()==FREE)
			{sA->setLast(m, m-di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+di]->getState()==FREE)
			{sA->setLast(m, m+di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+dii]->getState()==FREE)
			{sA->setLast(m, m+dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+diii]->getState()==FREE)
			{sA->setLast(m, m+diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+diii, 1);
			
			sA->set(m, m, - (a*FREE_DIFFUSION_SPEEDUP + 6-a) - 1/r/beta - consumption/r * timeStep);	
#else
			sA->setLast(m, m-diii, 1);
			sA->setLast(m, m-dii, 1);
			sA->setLast(m, m-di, 1);
			//float consumption = ( mi==0 ? 1/(kgg+glu)*(max - (max-min) * oxy/(kgo+oxy))
			//                            : 1/(kgg+oxy)*(max - (max-min) * glu/(kgo+glu)));
			//int tempState = voronoiDiagram->voronoiCells[v]->getState();
			sA->setLast(m, m, - 6 - 1/r/beta - consumption/r * timeStep);	
		
			sA->setLast(m, m+di, 1);
			sA->setLast(m, m+dii, 1);
			sA->setLast(m, m+diii, 1);
#endif

			
			}

			// vector
			//b[m] = 0.;
		
			//float consumption = ( mi==0 ? glu/(kgg+glu)*(max - (max-min) * oxy/(kgo+oxy))
			//                            : oxy/(kgg+oxy)*(max - (max-min) * glu/(kgo+glu)));
			//b[m] = - (mi==0 ? glu : oxy);
			
			if( actualize_vector)
			/*b[m] = (mi==0 ? -voronoiDiagram->voronoiCells[v]->glucose/r/beta 
			                -(1/beta-1)*(voronoiDiagram->voronoiCells[v-di]->glucose+voronoiDiagram->voronoiCells[v-dii]->glucose+voronoiDiagram->voronoiCells[v-diii]->glucose
			                           +voronoiDiagram->voronoiCells[v+di]->glucose+voronoiDiagram->voronoiCells[v+dii]->glucose+voronoiDiagram->voronoiCells[v+diii]->glucose
			                           -6*voronoiDiagram->voronoiCells[v]->glucose - consumption*voronoiDiagram->voronoiCells[v]->glucose/r * timeStep) 
			              : -voronoiDiagram->voronoiCells[v]->oxygen/r/beta 
			                -(1/beta-1)*(voronoiDiagram->voronoiCells[v-di]->oxygen+voronoiDiagram->voronoiCells[v-dii]->oxygen+voronoiDiagram->voronoiCells[v-diii]->oxygen
			                            +voronoiDiagram->voronoiCells[v+di]->oxygen+voronoiDiagram->voronoiCells[v+dii]->oxygen+voronoiDiagram->voronoiCells[v+diii]->oxygen
			                            -6*voronoiDiagram->voronoiCells[v]->oxygen - consumption*voronoiDiagram->voronoiCells[v]->oxygen/r * timeStep));
			*/
			b[m] = (mi==0 ? (1-1/beta)*(voronoiDiagram->voronoiCells[v-di]->glucose
			                           +voronoiDiagram->voronoiCells[v-dii]->glucose
			                           +voronoiDiagram->voronoiCells[v-diii]->glucose
			                           -6*voronoiDiagram->voronoiCells[v]->glucose
			                           +voronoiDiagram->voronoiCells[v+diii]->glucose
			                           +voronoiDiagram->voronoiCells[v+dii]->glucose
			                           +voronoiDiagram->voronoiCells[v+di]->glucose) // diffusion 
			               -(1-1/beta)*timeStep/r*consumption*voronoiDiagram->voronoiCells[v]->glucose // consumption
			               -voronoiDiagram->voronoiCells[v]->glucose/r/beta
			                
			              : (1-1/beta)*(voronoiDiagram->voronoiCells[v-di]->oxygen
			                           +voronoiDiagram->voronoiCells[v-dii]->oxygen
			                           +voronoiDiagram->voronoiCells[v-diii]->oxygen
			                           -6*voronoiDiagram->voronoiCells[v]->oxygen
			                           +voronoiDiagram->voronoiCells[v+diii]->oxygen
			                           +voronoiDiagram->voronoiCells[v+dii]->oxygen
			                           +voronoiDiagram->voronoiCells[v+di]->oxygen) // diffusion 
			               -(1-1/beta)*timeStep/r*consumption*voronoiDiagram->voronoiCells[v]->oxygen // consumption
			               -voronoiDiagram->voronoiCells[v]->oxygen/r/beta
			        );
			//b[m] = (voronoiDiagram->voronoiCells[v]->getState() != FREE ? consumption : 0.);
		}
		
		// border condition
		else{
			//matrix
			if( actualize_matrix)
			sA->setLast(m, m, 1);

			// vector
			if( actualize_vector)
			b[m] = (mi==0 ? voronoiDiagram->voronoiCells[v]->glucose
			              : voronoiDiagram->voronoiCells[v]->oxygen);
		}
		m++;
		v = (v + 1) % d;
	}
}


void setupTestMatrixCrankNicholsonNewton( VoronoiDiagram *voronoiDiagram, SparseMatrix *sA, float* x, char actualize_vector, char actualize_matrix, float timeStep, float spaceStep, float beta)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d    = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	//int nx = voronoiDiagram->xN[0],
	//    ny = voronoiDiagram->xN[1],
	//    nz = voronoiDiagram->xN[2];
	
	float r = 0.,
	      //rg = 1, //Glucose_Diffusion /(DX*DX), //1
	      //ro = 10; //Oxygen_Diffusion /(DX*DX);
	      rg = Glucose_Diffusion, //1
	      ro = Oxygen_Diffusion;
	//float r = 1.;
	
	int mi=0;
	int m = 0;
	int v = 0;
#ifdef COUPLING
	for( mi=0; mi<2; mi++)
#endif
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		
		// actual element of the matrix
		//int m = i*di + ii*dii + iii*diii  + mi*d;

		// actual element of the voronoi diagram
		//int v = i*di + ii*dii + iii*diii;
		//fprintf( stderr, "m:%i, v:%i\n ", m, v);

		// reset row: set all entries to zero
		if( actualize_matrix)
		sA->resetRow( m);

		// element inside domain
		if( 
			//voronoiDiagram->voronoiCells[v]->getState() != FREE &&
		    i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		)
		{
			//float consumption =  (mi==0 ? GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v]) 
			//                            : GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));
			//int radius = 3;
			/*float oxy = voronoiDiagram->voronoiCells[v]->oxygen; // 0.07
			float glu = voronoiDiagram->voronoiCells[v]->glucose;
			float min = 0.001,
			      max = 0.09,
			      kgg = 0.01,
			      kgo = 0.01;
			*/
			float consumption = 0.;
			float diagValue = 0.;
			//float linkingValue = 0.;

			if( actualize_matrix)
			{
				
			// diffusion
			r = (mi==0 ? rg : ro) * timeStep / (spaceStep*spaceStep);

			// consumption rate
			consumption = ( mi==0 ? GiveMeTheGlucoseRate(voronoiDiagram->voronoiCells[v], x[m], x[m+d])
			                      : GiveMeTheOxygenRate( voronoiDiagram->voronoiCells[v], x[m-d], x[m]));
			if( voronoiDiagram->voronoiCells[v]->getState() != ACTIVE && voronoiDiagram->voronoiCells[v]->getState() != NONACTIVE) consumption = 0.;
			
			// partial derivative of consumption
			// I. Glucose, II. Oxygen
			switch( mi){
				case 0: diagValue    = GlucoseConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				        //linkingValue = GlucoseConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m], x[m+d]);
				        //J->set(m, m+d, linkingValue/r);
				        break;
				case 1: diagValue    = OxygenConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        //linkingValue = OxygenConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				        //J->set(m, m-d, linkingValue/r);
				        break;
			}


#ifdef FREE_DIFFUSION_SPEEDUP
			float a = 0.;

			//matrix
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-diii]->getState()==FREE)
			{sA->setLast(m, m-diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-diii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-dii]->getState()==FREE)
			{sA->setLast(m, m-dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v-di]->getState()==FREE)
			{sA->setLast(m, m-di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m-di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+di]->getState()==FREE)
			{sA->setLast(m, m+di, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+di, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+dii]->getState()==FREE)
			{sA->setLast(m, m+dii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+dii, 1);
			if(voronoiDiagram->voronoiCells[v]->getState()==FREE && voronoiDiagram->voronoiCells[v+diii]->getState()==FREE)
			{sA->setLast(m, m+diii, FREE_DIFFUSION_SPEEDUP); a++;} else sA->set(m, m+diii, 1);
			
			sA->set(m, m, - (a*FREE_DIFFUSION_SPEEDUP + 6-a) - 1/r/beta - consumption/r * timeStep);	
#else
			sA->setLast(m, m-diii, 1);
			sA->setLast(m, m-dii, 1);
			sA->setLast(m, m-di, 1);
			
			//sA->setLast(m, m, - 6 - 1/r/beta - consumption/r * timeStep);	
			sA->setLast(m, m, - 6 - 1/r/beta - diagValue/r/beta * timeStep);	
		
			sA->setLast(m, m+di, 1);
			sA->setLast(m, m+dii, 1);
			sA->setLast(m, m+diii, 1);
#endif

			
			}

			if( actualize_vector)
			b[m] = (mi==0 ? (1-1/beta)*(voronoiDiagram->voronoiCells[v-di]->glucose
			                           +voronoiDiagram->voronoiCells[v-dii]->glucose
			                           +voronoiDiagram->voronoiCells[v-diii]->glucose
			                           -6*voronoiDiagram->voronoiCells[v]->glucose
			                           +voronoiDiagram->voronoiCells[v+diii]->glucose
			                           +voronoiDiagram->voronoiCells[v+dii]->glucose
			                           +voronoiDiagram->voronoiCells[v+di]->glucose) // diffusion 
			               //-(1-1/beta)*timeStep/r*consumption*voronoiDiagram->voronoiCells[v]->glucose // consumption
			               +consumption/r/beta - diagValue/r/beta*voronoiDiagram->voronoiCells[v]->glucose
			               -voronoiDiagram->voronoiCells[v]->glucose/r/beta
			                
			              : (1-1/beta)*(voronoiDiagram->voronoiCells[v-di]->oxygen
			                           +voronoiDiagram->voronoiCells[v-dii]->oxygen
			                           +voronoiDiagram->voronoiCells[v-diii]->oxygen
			                           -6*voronoiDiagram->voronoiCells[v]->oxygen
			                           +voronoiDiagram->voronoiCells[v+diii]->oxygen
			                           +voronoiDiagram->voronoiCells[v+dii]->oxygen
			                           +voronoiDiagram->voronoiCells[v+di]->oxygen) // diffusion 
			               //-(1-1/beta)*timeStep/r*consumption*voronoiDiagram->voronoiCells[v]->oxygen // consumption
			               +consumption/r/beta - diagValue/r/beta*voronoiDiagram->voronoiCells[v]->oxygen
			               -voronoiDiagram->voronoiCells[v]->oxygen/r/beta
			        );
		}
		
		// border condition
		else{
			//matrix
			if( actualize_matrix)
			sA->setLast(m, m, 1);

			// vector
			if( actualize_vector)
			b[m] = (mi==0 ? voronoiDiagram->voronoiCells[v]->glucose
			              : voronoiDiagram->voronoiCells[v]->oxygen);
		}
		m++;
		v = (v + 1) % d;
	}
}


void setupJacobianMatrixImplicitSteadyState( VoronoiDiagram *voronoiDiagram,  SparseMatrix *J, float *x, char molecule)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	int d = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	//float a = 0.5;
	//float b = 1. - a;
	//float c = 1.
	
	//fprintf( stderr, "Build Steady State Matrix: setupJacobianMatrixImplicitSteadyState()\n ");
	
	float r = 0.;
	//float r2;
	 
	switch( molecule){
		case 'G':
			r = Glucose_Diffusion /(DX*DX);
			break;
		case 'O':
			r = Oxygen_Diffusion /(DX*DX);
			break;
	}
	//r2 = r * 10.;
	
	//int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	for( int mi=0; mi<2; mi++)
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		if( mi==0) molecule ='G';
		if( mi==1) {molecule ='O'; r = Oxygen_Diffusion /(DX*DX); }
		
		// actual element
		// actual element of the matrix
		int m = mi*d   + i*di + ii*dii + iii*diii;

		// actual element of the voronoi diagram
		int v = i*di + ii*dii + iii*diii;
		//fprintf( stderr, "m:%i, v:%i\n ", m, v);

		// reset row: set all entries to zero
		sA->resetRow( m);

		// element inside domain
		if( 
			voronoiDiagram->voronoiCells[v]->getState() != VESSEL //TRUE
		#ifdef BORDER_IS_BOUNDARY
		    && i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1 
		#endif
		#ifdef FREE_IS_BOUNDARY
		    //&& voronoiDiagram->voronoiCells[v]->getState() != FREE
		    && voronoiDiagram->voronoiCells[v]->agent != NULL
		#endif
		)
		{
			//matrix
			
			// DIRICHLET
			
		/*#ifdef BORDER_IS_BOUNDARY
			#ifdef FREE_DIFFUSION_SPEEDUP
			int a = 0;
			
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v-diii]->getState() == FREE)
			{sA->set(m, m-diii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-diii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v-dii]->getState() == FREE)
			{sA->set(m, m-dii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-dii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v-di]->getState() == FREE)
			{sA->set(m, m-di, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m-di, r); 

			#ifndef FREE_DIFFUSION_SPEEDUP
			switch( molecule){
				case 'G': sA->set(m, m, - 6*r - GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v])); break;
				case 'O': sA->set(m, m, - 6*r - GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));  break;
			}
			#endif

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v+di]->getState() == FREE)
			{sA->set(m, m+di, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+di, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v+dii]->getState() == FREE)
			{sA->set(m, m+dii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+dii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			if( voronoiDiagram->voronoiCells[v]->getState() == FREE && voronoiDiagram->voronoiCells[v+diii]->getState() == FREE)
			{sA->set(m, m+diii, r * FREE_DIFFUSION_SPEEDUP); a++;}
			else
			#endif 
			sA->set(m, m+diii, r); 

			#ifdef FREE_DIFFUSION_SPEEDUP
			switch( molecule){
				case 'G': sA->set(m, m, - (6-a)*r - a*r*FREE_DIFFUSION_SPEEDUP - GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v])); break;
				case 'O': sA->set(m, m, - (6-a)*r - a*r*FREE_DIFFUSION_SPEEDUP - GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));  break;
			}
			#endif
		#else	

			// NEUMANN

			float diagValue = 0.;	
			switch( molecule){
				case 'G': sA->set(m, m, diagValue = -GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[v])); break;
				case 'O': sA->set(m, m, diagValue = -GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[v]));  break;
			}
	
			if(iii>0){ sA->set(m, m-diii, r); diagValue -= r;}
			if(ii>0){ sA->set(m, m-dii, r); diagValue -= r;}
			if( i>0){ sA->set(m, m-di, r); diagValue -= r;}
			if(i<voronoiDiagram->xN[0]-1){ sA->set(m, m+di, r); diagValue -= r;}
			if(ii<voronoiDiagram->xN[1]-1){ sA->set(m, m+dii, r); diagValue -= r;}
			if(iii<voronoiDiagram->xN[2]-1){ sA->set(m, m+diii, r); diagValue -= r;}
			sA->set(m, m, diagValue);	
		#endif*/
					//consumption
			float diagValue = 0.;
			float linkingValue = 0.;
			// I. Glucose, II. Oxygen
			switch( molecule){
				case 'G': diagValue    = - GlucoseConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[v], x[m], x[m+d]); 
				          linkingValue = - GlucoseConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[i], x[m], x[m+d]); 
				          J->set(m, m+d, linkingValue);
				          break;
				case 'O': diagValue    = - OxygenConsumptionRate_PartialDerivativeOxygen(  voronoiDiagram->voronoiCells[v], x[m-d], x[m]);
				          linkingValue = - OxygenConsumptionRate_PartialDerivativeGlucose( voronoiDiagram->voronoiCells[i], x[m-d], x[m]);
				          J->set(m, m-d, linkingValue);
				          break;
			}

			//matrixˀ
		#ifdef BORDER_IS_BOUNDARY
			// DIRICHLET
			diagValue += -6*r;
			J->set(m, m-diii, r); 
			J->set(m, m-dii, r); 
			J->set(m, m-di, r); 
			J->set(m, m, diagValue);					
			J->set(m, m+di, r); 
			J->set(m, m+dii, r); 
			J->set(m, m+diii, r); 
		#else	
			// NEUMANN
			if(iii>0){ J->set(m, m-diii, r); diagValue -= r;}
			if(ii>0){ J->set(m, m-dii, r); diagValue -= r;}
			if( i>0){ J->set(m, m-di, r); diagValue -= r;}
			if(i<voronoiDiagram->xN[0]-1){ J->set(m, m+di, r); diagValue -= r;}
			if(ii<voronoiDiagram->xN[1]-1){ J->set(m, m+dii, r); diagValue -= r;}
			if(iii<voronoiDiagram->xN[2]-1){ J->set(m, m+diii, r); diagValue -= r;}
			J->set(m, m, diagValue);	
		#endif

		}
		
		// border condition
		else{
			//matrix'
			J->set(m, m, 0);
		}
		
		/*fprintf( stderr, "A[%i]: ", m);
		for( int jj=0; jj<J->sizeA[m]; jj++){
			int j = J->JA[m][jj];
			fprintf( stderr, "%e ", J->A[m][jj]);
		}
		fprintf( stderr, "\n");
*/
		
	}
	//		fprintf( stderr, "...finished\n ");

}




///////////// OXYGEN CONSUMPTION  /////////////////////////

float GetOxygenConsumptionRate( VoronoiCell *cell, float glucose, float oxygen)
{
	if( cell->getState() != COMPARTMENT)
	if( cell->getState() != NONACTIVE && cell->getState() != ACTIVE)
		return 0.;


	// consumption per cell
	float k1 = Agent::NICK_O_CRITICAL_OXY;
	float k2 = (Agent::NICK_O_MAX - (Agent::NICK_O_MAX - Agent::NICK_O_MIN) * ( glucose / (glucose + Agent::NICK_O_CRITICAL_GLU))) * MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR;

	float oxygenUptake = k2 * oxygen / (oxygen + k1);

	// consumption per voronoi cell
	oxygenUptake *= 2/3.;

	if( cell->getState() == COMPARTMENT)
		oxygenUptake *= (double)GetAgent(cell)->cellCount/(double)GetAgent(cell)->maxCellCount;

	return oxygenUptake;
}

// OXYGEN CONSUMPTION: df( g, o)/do
float OxygenConsumptionRate_PartialDerivativeOxygen( VoronoiCell *cell, float glucose, float oxygen)
{
	if( cell->getState() != COMPARTMENT)
	if( cell->getState() != NONACTIVE && cell->getState() != ACTIVE)
		return 0.;
	

	// consumption per cell
	float k1 = Agent::NICK_O_CRITICAL_OXY;
	float k2 = Agent::NICK_O_MAX;
	float k3 = Agent::NICK_O_MIN;
	float k4 = Agent::NICK_O_CRITICAL_GLU;
	
	float oxygenUptake = k1 / pow(oxygen+k1, 2.) * (k2 - (k2-k3) * glucose/(glucose + k4));

	// consumption per voronoi cell
	oxygenUptake *= 2/3.;

	if( cell->getState() == COMPARTMENT)
		oxygenUptake *= (double)GetAgent(cell)->cellCount/(double)GetAgent(cell)->maxCellCount;

	return oxygenUptake;
}

// OXYGEN CONSUMPTION: df( g, o)/dg
float OxygenConsumptionRate_PartialDerivativeGlucose( VoronoiCell *cell, float glucose, float oxygen)
{
	if( cell->getState() != COMPARTMENT)
	if( cell->getState() != NONACTIVE && cell->getState() != ACTIVE)
		return 0.;
	

	// consumption per cell
	float k1 = Agent::NICK_O_CRITICAL_OXY;
	float k2 = Agent::NICK_O_MAX;
	float k3 = Agent::NICK_O_MIN;
	float k4 = Agent::NICK_O_CRITICAL_GLU;
	
	float oxygenUptake = - oxygen/(oxygen+k1) * (k2-k3) * k4 / pow(glucose + k4, 2.);

	// consumption per voronoi cell
	oxygenUptake *= 2/3.;

	if( cell->getState() == COMPARTMENT)
		oxygenUptake *= (double)GetAgent(cell)->cellCount/(double)GetAgent(cell)->maxCellCount;

	return oxygenUptake;
}



/////////////// GLUCOSE CONSUMPTION //////////////////////////

float GetGlucoseConsumptionRate( VoronoiCell *cell, float glucose, float oxygen)
{
	if( cell->getState() != COMPARTMENT)
	if( cell->getState() != NONACTIVE && cell->getState() != ACTIVE)
		return 0.;


	// consumption per cell
	float k1 = Agent::NICK_G_CRITICAL_GLU;
	float k2 = (Agent::NICK_G_MAX - (Agent::NICK_G_MAX - Agent::NICK_G_MIN) * ( oxygen / (oxygen + Agent::NICK_G_CRITICAL_OXY))) * MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR;

	float glucoseUptake = k2 * glucose / (glucose + k1);

	// consumption per voronoi cell
	glucoseUptake *= 2/3.;

	if( cell->getState() == COMPARTMENT)
		glucoseUptake *= (double)GetAgent(cell)->cellCount/(double)GetAgent(cell)->maxCellCount;

	return glucoseUptake;
}

float GlucoseConsumptionRate_PartialDerivativeGlucose( VoronoiCell *cell, float glucose, float oxygen)
{
	if( cell->getState() != COMPARTMENT)
	if( cell->getState() != NONACTIVE && cell->getState() != ACTIVE)
		return 0.;
	

	// consumption per cell
	float k1 = Agent::NICK_G_CRITICAL_GLU;
	float k2 = Agent::NICK_G_MAX;
	float k3 = Agent::NICK_G_MIN;
	float k4 = Agent::NICK_G_CRITICAL_OXY;
	
	float glucoseUptake = k1 / pow(glucose+k1, 2.) * (k2 - (k2-k3) * oxygen/(oxygen + k4));

	// consumption per voronoi cell
	glucoseUptake *= 2/3.;

	if( cell->getState() == COMPARTMENT)
		glucoseUptake *= (double)GetAgent(cell)->cellCount/(double)GetAgent(cell)->maxCellCount;

	return glucoseUptake;
}

float GlucoseConsumptionRate_PartialDerivativeOxygen( VoronoiCell *cell, float glucose, float oxygen)
{
	if( cell->getState() != COMPARTMENT)
	if( cell->getState() != NONACTIVE && cell->getState() != ACTIVE)
		return 0.;
	

	// consumption per cell
	float k1 = Agent::NICK_G_CRITICAL_GLU;
	float k2 = Agent::NICK_G_MAX;
	float k3 = Agent::NICK_G_MIN;
	float k4 = Agent::NICK_G_CRITICAL_OXY;
	
	float glucoseUptake = - glucose/(glucose+k1) * (k2-k3) * k4 / pow(oxygen + k4, 2.);

	// consumption per voronoi cell
	glucoseUptake *= 2/3.;

	if( cell->getState() == COMPARTMENT)
		glucoseUptake *= (double)GetAgent(cell)->cellCount/(double)GetAgent(cell)->maxCellCount;

	return glucoseUptake;
}


float GetLactateProductionRate( VoronoiCell *cell, float glucose, float oxygen)
{
	float cg = GetGlucoseConsumptionRate( cell, glucose, oxygen);
	float co = GetOxygenConsumptionRate(  cell, glucose, oxygen);

	return 2*MAX( 0, cg - co/6.);
}



double UpdateSystemNonLinearCGSparse( VoronoiDiagram *voronoiDiagram, double time, double end_time, double timeStep)
{
	//double time;
	
	int N = voronoiDiagram->xN[0]
#if DIMENSIONS >= 2
	        *voronoiDiagram->xN[1]
#endif
#if DIMENSIONS >= 3
*voronoiDiagram->xN[2]
#endif
	        ;

	float max_err, mean_err;

	int i;

	float *temp = v8;

#ifndef COUPLING
	sA->dimI = N;
	sA->dimJ = N;
#else
	sA->dimI = 2*N;
	sA->dimJ = 2*N;
#endif
	// timestep loop
	// IMPLICIT
	//for( ; time+timeStep <= end_time; time += timeStep)
	// STEADY STATE
	time = end_time;
	{
		fprintf( stderr, "time = %lf / %lf (dt = %lf) ...\n", time, end_time, timeStep);

#pragma omp parallel for
		for(int m=0; m<N; m++){
			temp[m]= voronoiDiagram->voronoiCells[m]->glucose;
			x[m]   = voronoiDiagram->voronoiCells[m]->glucose;
#ifdef COUPLING
			temp[m+N] = x[m+N] = voronoiDiagram->voronoiCells[m]->oxygen;
#endif
		}


		// iteration loop (non-linearity)
		i = 0;
		do{
//			fprintf( stderr, "Outer Iteration %i...\n", i+1);

			// set x 

			//fprintf( stderr, "initialize vector x\n");
			/*for(int m=0; m<N; m++){
				temp[m]= voronoiDiagram->voronoiCells[m]->glucose;
				x[m]   = voronoiDiagram->voronoiCells[m]->glucose;
#ifdef COUPLING
				temp[m+N] = x[m+N] = voronoiDiagram->voronoiCells[m]->oxygen;
#endif
			}*/
#ifdef COUPLING
			N *= 2; 
#endif
			// set implicit matrix
			//setupMatrixImplicitSteadyState( voronoiDiagram, 'G');
		
	#ifdef PARTLY_LINEARIZED		
			// STEADY STATE (partly linearized)
			setupTestMatrixSteadyState( voronoiDiagram, sA, x, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX); time = end_time;
			//setupTestMatrixSteadyStateNewton( voronoiDiagram, sA, x, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX); time = end_time;
			//setupTestMatrixImplicit( voronoiDiagram, 0, timeStep , 1.);

	#else
			// STEADY STATE (fully linearized)
			setupTestMatrixSteadyState( voronoiDiagram, ACTUALIZE_VECTOR, (i==0 ? ACTUALIZE_MATRIX : DONT_ACTUALIZE_MATRIX); time = end_time;
	#endif
	
			// TIME EVOLVING
			//float beta = 1.;
			//setupTestMatrixCrankNicholson( voronoiDiagram, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep , DX, beta);
			//setupTestMatrixCrankNicholsonNewton( voronoiDiagram, sA, x, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep , DX, beta);
			
			/*float beta = 1.;
			float pre_beta = 1.;
			if( time < timeStep){
				int it=1;
				//int it=4;
				for( int ii=1; ii<it; ii++){
					setupTestMatrixCrankNicholson( voronoiDiagram, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep/(double)it , 1., pre_beta); // 0 Explicit, 0.5 Crank-Nicholson, 1 Implicit
					SolveBiCGSTAB( sA, b, x, 1000, CG_CRITICAL_ERROR); // WORKS
					fprintf( stderr, "Implicit! %lf\n", timeStep/(double)it);
					for(int m=0; m<N/2; m++){
						voronoiDiagram->voronoiCells[m]->glucose = (x[m]  >=0.?x[m]:0.);
						voronoiDiagram->voronoiCells[m]->oxygen  = (x[m+N/2]>=0.?x[m+N/2]:0.);
					}
				}
				setupTestMatrixCrankNicholson( voronoiDiagram, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep/(double)it , 1., pre_beta); // 0 Explicit, 0.5 Crank-Nicholson, 1 Implicit
				//setupTestMatrixCrankNicholson( voronoiDiagram, 0, timeStep/2. , 1., 1.); // 0 Explicit, 0.5 Crank-Nicholson, 1 Implicit
			}else
				setupTestMatrixCrankNicholson( voronoiDiagram, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep , 1., beta); // 0 Explicit, 0.5 Crank-Nicholson, 1 Implicit
		*/
			// solve linear system
			//JacobiPreconditioner( sA, M);
			//SuccessiveOverRelaxationPreconditioner( sA, M);
			//PreconditionedConjugateGradientSparse( sA, M, b, x, 10000, N*1e-4);
			//ConjugateGradientSparse( sA, b, x, 200, CG_CRITICAL_ERROR);
			//SolveExplicit( sA, b, x);
			SolveBiCGSTAB( sA, b, x, 1000, /*N*/CG_CRITICAL_ERROR); // WORKS
			//SolveJacobiMethod( sA, b, x, 200, 1e-20); // WORKS FINE!!!
			//SolveGaussSeidelMethod( sA, b, x, 200, N*1e-4); // WORKS

			// change
			/*vectorDifference( temp, x, temp, N);
			old_change = change;
			change = sqrt( dotProduct( temp, temp, N));

			// error
			abs_b = sqrt( dotProduct( b, b, N));*/
		
		
//#define EXPLICIT
#ifndef EXPLICIT
			//sparseMatrixVectorProduct( sA, x, temp);
			//vectorDifference( temp, b, temp, N);
			for(int m=0; m<N; m++)
				if(x[m]<0.)x[m]=0.;
			vectorDifference( temp, x, temp, N);
			
			
			// max error
			//old_max_err = max_err;
			max_err = 0.; 
			mean_err = 0.;
			for(int m=0; m<N; m++){
				mean_err += fabs(temp[m] /* x[m]*/);
				if( max_err < fabs(temp[m] /* x[m]*/))
					max_err = fabs(temp[m] /* x[m]*/);
				temp[m]=x[m];
			}
			mean_err /= N;
			//if(!i)
				//old_max_err = max_err;
#endif

#ifdef COUPLING
			N /= 2; 
#endif
//			fprintf( stderr, "...finished after %i Iterations -> max err.: %.3e, mean err.: %.3e\n", i, max_err, mean_err);
			//fprintf( stdout, "%i %lf %lf %lf\n", i, time, max_err, mean_err);
			/*for(int m=0; m<N; m++){
				voronoiDiagram->voronoiCells[m]->glucose = (x[m]  >=0.?x[m]:0.);
#ifdef COUPLING
				voronoiDiagram->voronoiCells[m]->oxygen  = (x[m+N]>=0.?x[m+N]:0.);
#endif
			}*/
			i++;

		}while( /*err/abs_b > 5e-5 &&*/ /*old_change != change && change > 1e-10 && */ 
	        i < 20 &&
	        max_err > CG_CRITICAL_ERROR * 10
	        //false
	        /*&& old_max_err >= max_err*/
	      );

		// end of iteration loop
		
		for(int m=0; m<N; m++){
			voronoiDiagram->voronoiCells[m]->glucose = (x[m]  >=0.?x[m]:0.);
#ifdef COUPLING
			voronoiDiagram->voronoiCells[m]->oxygen  = (x[m+N]>=0.?x[m+N]:0.);
#endif
		}


	
		// OUTPUT of ERROR
		//char filename[512];
		//FILE *fp;
/*#ifdef PARTLY_LINEARIZED
		sprintf( filename, "error_semiLinSteadyState.dat");
#else
		sprintf( filename, "error_fullLinSteadystate.dat");
#endif
*/		/*sprintf( filename, "error_CrankNicholson.dat");
		if(( fp = fopen( filename, ( time==0. ? "w" : "a")))==NULL){
			fprintf( stderr, "Error opening file %s for writing!\n", filename);
			exit(0);
		}
	
		fprintf( fp, "%lf %e %e %i\n", time+timeStep, max_err, mean_err, i);
		fclose( fp);	
	*/
	}
	
	// end of time loop
	
//	fprintf( stderr, "\n========================\n");

	return end_time - time;
}

double UpdateGrowthFactorsNonLinearCGSparse( VoronoiDiagram *voronoiDiagram, double time, double end_time, double timeStep)
{
	int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	float max_err, mean_err;
	int i;
	
	float *temp = v8;

	sA->dimI = N;
	sA->dimJ = N;

	// timestep loop
	// IMPLICIT
	for( ; time+timeStep <= end_time; time += timeStep)
	// STEADY STATE
	
	{
		fprintf( stderr, "time = %lf / %lf (dt = %lf) ...\n", time, end_time, timeStep);


		for(int m=0; m<N; m++){
			x[m]   = voronoiDiagram->voronoiCells[m]->growthfactors;
		}


		// iteration loop (non-linearity)
		i = 0;
		do{
			fprintf( stderr, "Outer Iteration %i...\n", i+1);

			// SET MATRIX
		
	/*#ifdef PARTLY_LINEARIZED		
			// STEADY STATE (partly linearized)
			setupGrowthFactorMatrixSteadyState( voronoiDiagram, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX); time = end_time;

	#else
			// STEADY STATE (fully linearized)
			setupGrowthFactorMatrixSteadyState( voronoiDiagram, ACTUALIZE_VECTOR, (i==0 ? ACTUALIZE_MATRIX : DONT_ACTUALIZE_MATRIX); time = end_time;
	#endif*/
	
			// TIME EVOLVING
			//setupTestMatrixCrankNicholson( voronoiDiagram, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep , 1., beta); // 0 Explicit, 0.5 Crank-Nicholson, 1 Implicit
			setupGrowthFactorMatrixImplicit( voronoiDiagram, sA, b, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep, DX);

			// SOLVE SYSTEM
			SolveBiCGSTAB( sA, b, x, 1000, /*N*/CG_CRITICAL_ERROR); // WORKS


		
			// ERROR ESTIMATION
			sparseMatrixVectorProduct( sA, x, temp);
			vectorDifference( temp, b, temp, N);

			//old_max_err = max_err;
			max_err = 0.; 
			mean_err = 0.;
			for(int m=0; m<N; m++){
				mean_err += fabs(temp[m] /* x[m]*/);
				if( max_err < fabs(temp[m] /* x[m]*/))
					max_err = fabs(temp[m] /* x[m]*/);
			}
			mean_err /= N;
			//if(!i)
			//	old_max_err = max_err;

			fprintf( stderr, "...finished after %i Iterations -> max err.: %.3e, mean err.: %.3e\n", i, max_err, mean_err);

			i++;

		}while(  
	        i < 20 &&
	        max_err > CG_CRITICAL_ERROR * 10
	      );

		// end of iteration loop
		
		for(int m=0; m<N; m++){
			voronoiDiagram->voronoiCells[m]->growthfactors = (x[m]  >=0.?x[m]:0.);
		}
	}
	
	// end of time loop
	
	fprintf( stderr, "\n========================\n");

	return end_time - time;
}

double UpdateSystemNewtonSparse( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference)
{
	//double time;
	
	int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	float f[2*N];
	float h[2*N];
	//SparseMatrix *J = SparseMatrix::newSparseMatrix( N, N);
	float err, cg_err;

	// set x
	//fprintf( stderr, "initialize vector x\n");
	for(int m=0; m<N; m++){
		x[m] = voronoiDiagram->voronoiCells[m]->glucose;
		x[m+N] = voronoiDiagram->voronoiCells[m]->oxygen;

		h[m]=0.;
		h[m+N]=0.;
	}

	// timestep loop
//	for( time = 0; time+timeStep <= timeDifference; time += timeStep){
//	for(int i=0; i<10; i++){
	int i = 0;
	do{
		fprintf( stderr, "Outer Iteration %i...\n", i+1);

		//// GLUCOSE ////////////////////

		// set A and b
		//fprintf( stderr, "setupMatrixImplicitSteadyState()\n");
		setupMatrixImplicitSteadyState( voronoiDiagram, 'G');

		// newton loop
		int it = 0;
		//do{
//			fprintf( stderr, "Newton Iteration %i\n", (int)((time + timeStep)/timeStep));
			
			// set f = b - A*x
			//fprintf( stderr, "sparseMatrixVectorProduct()\n");
			sparseMatrixVectorProduct( sA, x, f);
			//fprintf( stderr, "vectorDifference()\n");
			vectorDifference( b, f, f, 2*N);

			// set J
			//fprintf( stderr, "setupJacobianMatrixImplicitSteadyState()\n");
			setupJacobianMatrixImplicitSteadyState( voronoiDiagram, J, x, 'G');

			// solve h*J = f
			//ConjugateGradientSparse( J, f, h, 2*N, 2);
			//fprintf( stderr, "JacobiPreconditioner()\n");
			//SuccessiveOverRelaxationPreconditioner( J, M);
			JacobiPreconditioner( J, M);
			//fprintf( stderr, "PreconditionedConjugateGradientSparse()\n");
			cg_err = PreconditionedConjugateGradientSparse( J, M, f, h, 2*N, 100);

			// update x = x + h
			//fprintf( stderr, "vectorSum()\n");
			vectorSum( x, h, x, 2*N);

			err = sqrt( dotProduct( h, h, 2*N));
			it++;


		//}while( err > 1e-5 && it < 100);
		//fprintf( stderr, "...finished after %i Newton Iterations -> error: %e\n", i, err);
		for(int m=0; m<N; m++){
			voronoiDiagram->voronoiCells[m]->glucose = (x[m]  >=0.?x[m]:0.);
			voronoiDiagram->voronoiCells[m]->oxygen  = (x[m+N]>=0.?x[m+N]:0.);
		}




		//// OXYGEN ////////////////////
/*
		// set A and b
		setupMatrixImplicitSteadyState( voronoiDiagram, 'O');

		// newton loop
		it = 0;
		do{
//			fprintf( stderr, "Newton Iteration %i\n", (int)((time + timeStep)/timeStep));
			
			// set f = b - A*x
			sparseMatrixVectorProduct( sA, x, f, N);
			vectorDifference( b, f, f, N);
		
			// set J
			setupJacobianMatrixImplicitSteadyState( voronoiDiagram, J, sA, b, x, 'O');
		
			// solve h*J = f
			ConjugateGradientSparse( J, f, h, N, 2);

			// update x = x + h
			vectorSum( x, h, x, N);

			err = sqrt( dotProduct( h, h, N));
			it++;

			fprintf( stderr, "\r\t\t\t\t\t  Newton: error = %e \b", err);

			//fprintf( stderr, "Newton Iteration %i -> error: %e\n", it, err);
		}while( err > 1e-5 && it < 10);
		//fprintf( stderr, "...finished after %i Newton Iterations -> error: %e\n", i, err);
		for(int m=0; m<N; m++)
			voronoiDiagram->voronoiCells[m]->oxygen = x[m];

		fprintf( stderr, "\rTotal: error = %e \b", err);
*/
		i++;
	}while( err > NEWTON_CRITICAL_ERROR /*&& i < 100*/ || cg_err > CG_CRITICAL_ERROR);

	fprintf( stderr, "\n========================\n");

	return 0.;
}

double UpdateSystemNewtonCGSparse( VoronoiDiagram *voronoiDiagram, double time, double end_time, double timeStep)
{
	//double time;

	int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	//float err;
	float max_err, mean_err;
	//float old_max_err = 1e-5;
	//float change = 0., old_change;
	//float abs_b;
	int i;

	//float *temp = v0;
	//float *temp = v8;
	float *f = v8;
	float *h = v7;

#ifndef COUPLING
	sA->dimI = N;
	sA->dimJ = N;
#endif


		// set x
		for(int m=0; m<N; m++){
			//temp[m]   =
			x[m]   = voronoiDiagram->voronoiCells[m]->glucose;
			//temp[m+N] =
			x[m+N] = voronoiDiagram->voronoiCells[m]->oxygen;
		}
		N *= 2;

	// timestep loop
	// IMPLICIT
	//for( ; time+timeStep <= end_time; time += timeStep)
	// STEADY STATE
	time = end_time;
	{
		fprintf( stderr, "time = %lf / %lf (dt = %lf) ...\n", time, end_time, timeStep);

		//float beta = 1.;

		// build linear system:
		// set A & b
		//setupTestMatrixCrankNicholson( voronoiDiagram, ACTUALIZE_VECTOR, ACTUALIZE_MATRIX, timeStep , 1., beta); // 0 Explicit, 0.5 Crank-Nicholson, 1 Implicit
		//setupTestMatrixSteadyState( voronoiDiagram, 0);

		// set f = b - A*x
		//sparseMatrixVectorProduct( sA, x, f);
		//vectorDifference( b, f, f, N);
		//vectorDifference( f, b, f, N);

		// iteration loop (non-linearity)
		i = 0;
		do{
			fprintf( stderr, "Outer Iteration %i...\n", i+1);

			// IMPLICIT
			//setupTestMatrixCrankNicholson( voronoiDiagram, ( i==0 ? ACTUALIZE_VECTOR : DONT_ACTUALIZE_VECTOR) , ACTUALIZE_MATRIX, timeStep , 1., (time<=1.?1.:beta));
			// STEADY STATE
			setupTestMatrixSteadyState( voronoiDiagram, sA, x, ( i==0 ? ACTUALIZE_VECTOR : DONT_ACTUALIZE_VECTOR) , ACTUALIZE_MATRIX);

			// set f = b - A*x
			sparseMatrixVectorProduct( sA, x, f);
			vectorDifference( b, f, f, N);
			//vectorDifference( f, b, f, N);


			// set J (Jacobi-Matrix)
			// IMPLICIT
			//setupJacobiMatrixCrankNicholson( voronoiDiagram, J, x, timeStep , 1., (time<=1.?1.:beta));
			// STEADY STATE
			setupJacobiMatrixSteadyState( voronoiDiagram, J, x);

			// solve h*J = f
			SolveBiCGSTAB( J, f, h, 1000, CG_CRITICAL_ERROR); // WORKS

			// solve x*J = f = 0
			//for( int m=0; m<N; m++) f[m] = 0.;
			//SolveBiCGSTAB( J, f, x, 1000, CG_CRITICAL_ERROR/10); // WORKS

			// update x = x + h
			vectorSum( x, h, x, N);

			// update A
			N /= 2;
			for(int m=0; m<N; m++){
				voronoiDiagram->voronoiCells[m]->glucose = (x[m]  >=0.?x[m]:0.);
				voronoiDiagram->voronoiCells[m]->oxygen  = (x[m+N]>=0.?x[m+N]:0.);
			}
			N *= 2;


			//sparseMatrixVectorProduct( J, x, f);

			// max error
			///old_max_err = max_err;
			max_err = 0.;
			mean_err = 0.;
			for(int m=0; m<N; m++){
				mean_err += fabs(h[m] /* x[m]*/);
				if( max_err < fabs(h[m] /* x[m]*/))
					max_err = fabs(h[m] /* x[m]*/);
			}
			mean_err /= N;
			//if(!i)
				//old_max_err = max_err;


			fprintf( stderr, "...finished after %i Iterations -> max err.: %.3e, mean err.: %.3e\n", i, max_err, mean_err);
			i++;

		}while(
		        i < 20 &&
		        max_err > CG_CRITICAL_ERROR*10
		);

		// end of iteration loop

		N /= 2;
		for(int m=0; m<N; m++){
			voronoiDiagram->voronoiCells[m]->glucose = (x[m]  >=0.?x[m]:0.);
			voronoiDiagram->voronoiCells[m]->oxygen  = (x[m+N]>=0.?x[m+N]:0.);
		}
		N *= 2;
		
	}

	// end of time loop

	fprintf( stderr, "\n========================\n");

	return end_time - time;
}

double UpdateSystemImplicitSparse( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference)
{

	double time;
	
	int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	//int itG, itO;
	int passedTime = clock();
	fprintf( stderr, "ConjugateGradientSparse... \n");
	
	for( time = 0; time+timeStep <= timeDifference; time += timeStep){		
	//do{
		// glucose  ///////////////
		
		//int passedTime = clock();
		//fprintf( stderr, "Set Sparse Matrix... \n");
		#ifdef CRANK_NICHOLSON
		setupMatrixImplicitCrankNicholsonSparse( voronoiDiagram, timeStep, 'G');
		#else 
			#ifdef STEADY_STATE
		setupMatrixImplicitSteadyState( voronoiDiagram, 'G');
			#else
		setupMatrixImplicitSparse( voronoiDiagram, timeStep, 'G');
			#endif 
		#endif
		//fprintf( stderr, "...finished ( %li clocks, %.3lf sec)\n", (clock() - passedTime), (float)(clock() - passedTime)/CLOCKS_PER_SEC);

		for(int m=0; m<N; m++)
		#ifndef INIT_GLU
			x[m] = voronoiDiagram->voronoiCells[m]->glucose;
		#else
			x[m] = INIT_GLU;
		#endif
		
		passedTime = clock();
		//fprintf( stderr, "ConjugateGradientSparse... \n");
		#ifdef PRECONDITIONED
		JacobiPreconditioner( sA, M);
		PreconditionedConjugateGradientSparse( sA, M, b, x, N, 1000);
		#else
		/*itG = */ConjugateGradientSparse( sA, b, x, N, 1000);
		#endif
		//fprintf( stderr, "...finished ( %li clocks, %.3lf sec)\n", (clock() - passedTime), (float)(clock() - passedTime)/CLOCKS_PER_SEC);

		for(int m=0; m<N; m++)
		#ifdef NO_NEGATIVES
			voronoiDiagram->voronoiCells[m]->glucose = (x[m]>=0.?x[m]:0.);
		#else
			voronoiDiagram->voronoiCells[m]->glucose = x[m]; 
		#endif



		//   oxygen   //////////////////
		
		//passedTime = clock();
		//fprintf( stderr, "Set Sparse Matrix... \n");
		#ifdef CRANK_NICHOLSON
		setupMatrixImplicitCrankNicholsonSparse( voronoiDiagram, timeStep, 'O');
		#else
			#ifdef STEADY_STATE
		setupMatrixImplicitSteadyState( voronoiDiagram, 'O');
			#else
		setupMatrixImplicitSparse( voronoiDiagram, timeStep, 'O');
			#endif
		#endif
		//fprintf( stderr, "...finished ( %li clocks, %.3lf sec)\n", (clock() - passedTime), (float)(clock() - passedTime)/CLOCKS_PER_SEC);

		for(int m=0; m<N; m++)
		#ifndef INIT_OXY
			x[m] = voronoiDiagram->voronoiCells[m]->oxygen;
		#else
			x[m] = INIT_OXY;
		#endif
		
		//passedTime = clock();
		//fprintf( stderr, "ConjugateGradientSparse... \n");
		#ifdef PRECONDITIONED
		JacobiPreconditioner( sA, M);
		PreconditionedConjugateGradientSparse( sA, M, b, x, N, 1000);
		#else
		/*itO = */ConjugateGradientSparse( sA, b, x, N, 1000);
		#endif
		//fprintf( stderr, "...finished ( %li clocks, %.3lf sec)\n", (clock() - passedTime), (float)(clock() - passedTime)/CLOCKS_PER_SEC);

		for(int m=0; m<N; m++)
		#ifdef NO_NEGATIVES
			voronoiDiagram->voronoiCells[m]->oxygen = (x[m]>=0.?x[m]:0.);
		#else
			voronoiDiagram->voronoiCells[m]->oxygen = x[m]; 
		#endif
	}
	//}while( itG > 10 || itO > 10);

	fprintf( stderr, "...finished ( %li clocks, %.3lf sec)\n", (clock() - passedTime), (float)(clock() - passedTime)/CLOCKS_PER_SEC);


	return timeDifference - time;
	
}



void SolveJacobiMethod( SparseMatrix *A, float *b, float *x0, int maxit, float minerr)
{
	float *x1 = v0;
	float *temp;
	float err = minerr;

	for( int m=0; m<maxit && err >= minerr; m++){
		for( int i=0; i<A->dimI; i++){
			x1[i] = 0.;
			for( int jj=0; jj<A->sizeA[i]; jj++){
				int j=A->JA[i][jj];
				if( j!=i)
					x1[i] += A->A[i][jj] * x0[j];
			}
			x1[i] = (b[i]-x1[i]) / A->get(i,i);
			// max error
			err = ( err>fabs(x1[i]-x0[i]) ? err : fabs(x1[i]-x0[i]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;
	}
}

void SolveJacobiMethod( SparseMatrix *A, double *b, double *x0, int maxit, double minerr)
{
	char filename[512];
	FILE *fp;

	sprintf( filename, "convergenceJacobi_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();

	double *x1 = v0d;
	double *temp;
	double err = minerr;

	for( int m=0; m<maxit && err >= minerr; m++){
		for( int i=0; i<A->dimI; i++){
			x1[i] = 0.;
			for( int jj=0; jj<A->sizeA[i]; jj++){
				int j=A->JA[i][jj];
				if( j!=i)
					x1[i] += A->A[i][jj] * x0[j];
			}
			x1[i] = (b[i]-x1[i]) / A->get(i,i);
			if( isnan(x1[i]))
				exit( 0);
			// max error
			err = ( err>fabs(x1[i]-x0[i]) ? err : fabs(x1[i]-x0[i]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;

		// maximal error
		double myError = 0;
		double myError2 = 0;
		for(int n=0; n<A->dimI; n++){
			//if( myError < fabs(r[m]))
			//	myError = fabs(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( myError2 < fabs(exactSolution[n]-x0[n]))
				myError2 = fabs(exactSolution[n]-x0[n]);
#endif
		}

		//fprintf( fp, "%i %e %e %e \n", m, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);
	}
	fclose(fp);
}

void SolveGaussSeidelMethod( SparseMatrix *A, float *b, float *x0, int maxit, float minerr)
{
	float *x1 = v0;
	float *temp;
	float err = minerr;

	for( int m=0; m<maxit && err >= minerr; m++){
		for( int k=0; k<A->dimI; k++){
			x0[k] = x1[k];
			
			x1[k] = 0.;
			for( int ii=0; ii<A->sizeA[k]; ii++){
				int i=A->JA[k][ii];
				if( i<k)
					x1[k] += A->A[k][ii] * x1[i];
				if( i>k)
					x1[k] += A->A[k][ii] * x0[i];
			}
			
			x1[k] = (b[k]-x1[k]) / A->get(k,k);
			err = ( err>fabs(x1[k]-x0[k]) ? err : fabs(x1[k]-x0[k]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;
	}
}

void SolveGaussSeidelMethod( SparseMatrix *A, double *b, double *x0, int maxit, double minerr)
{
	char filename[512];
	FILE *fp;

	sprintf( filename, "convergenceGaussSeidel_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();

	double *x1 = v0d;
	double *temp;
	double err = minerr;

	for( int m=0; m<maxit && err >= minerr; m++){
		for( int k=0; k<A->dimI; k++){
			//x0[k] = x1[k];

			x1[k] = 0.;
			for( int ii=0; ii<A->sizeA[k]; ii++){
				int i=A->JA[k][ii];
				if( i<k)
					x1[k] += A->A[k][ii] * x1[i];
				if( i>k)
					x1[k] += A->A[k][ii] * x0[i];
			}

			x1[k] = (b[k]-x1[k]) / A->get(k,k);
			err = ( err>fabs(x1[k]-x0[k]) ? err : fabs(x1[k]-x0[k]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;

		// maximal error
		double myError = 0;
		double myError2 = 0;
		for(int n=0; n<A->dimI; n++){
			//if( myError < fabs(r[m]))
			//	myError = fabs(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( myError2 < fabs(exactSolution[n]-x0[n]))
				myError2 = fabs(exactSolution[n]-x0[n]);
#endif
		}

		//fprintf( fp, "%i %e %e %e \n", m, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);
	}

	fclose(fp);
}

void SolveSuccessiveOverRelaxationMethod( SparseMatrix *A, double *b, double *x0, int maxit, double minerr, double omega)
{
	char filename[512];
	FILE *fp;

	sprintf( filename, "convergenceSOR_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();

	double *x1 = v0d;
	double *temp;
	double err = minerr;

	double _omega = 1.-omega;

	for( int m=0; m<maxit && err >= minerr; m++){
		for( int i=0; i<A->dimI; i++){
			//x0[k] = x1[k];
			double Aii = 1.;
			x1[i] = 0.;
			for( int jj=0; jj<A->sizeA[i]; jj++){
				int j=A->JA[i][jj];
				if( j<i)
					x1[i] += A->A[i][jj] * x1[j];
				else if( j==i)
					Aii = A->A[i][jj];
				else if( j>i)
					x1[i] += A->A[i][jj] * x0[j];
			}

			//x1[k] = (1.-omega)*x0[k] + (b[k]-x1[k]) * omega / A->get(k,k);
			x1[i] = _omega*x0[i] + (b[i]-x1[i]) * omega / Aii;
			err = ( err>fabs(x1[i]-x0[i]) ? err : fabs(x1[i]-x0[i]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;

		// maximal error
		double myError = 0;
		double myError2 = 0;
		for(int n=0; n<A->dimI; n++){
			//if( myError < fabs(r[m]))
			//	myError = fabs(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( myError2 < fabs(exactSolution[n]-x0[n]))
				myError2 = fabs(exactSolution[n]-x0[n]);
#endif
		}

		fprintf( fp, "%i %e %e %e \n", m, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);
	}

	fclose(fp);
}

void ConjugateGradientSparse( SparseMatrix *A, float *b, float *x, int maxit, float minerr)
// If A is symmetric positive matrix
{
	int N = A->dimI;
	
	// vectors
	float r[N];
	float w[N];
	float z[N];
	
	// scalars
	float alpha;
	float beta;
	
	float temp[N];
	
	//float err = 1e-10;
	float myError = 0.;
	
	// initialization of residual vector: r
	sparseMatrixVectorProduct( A, x, temp);
	vectorDifference( b, temp, r, N);
	
	// w
	vectorScale( r, -1, w, N);
	
	// z
	sparseMatrixVectorProduct( A, w, z);
	
	// alpha
	alpha = dotProduct( r, w, N) / dotProduct( w, z, N);
	
	// beta
	beta = 0.;
	
	// x
	vectorScale( w, alpha, temp, N);
	vectorSum( x, temp, x, N);
	
	for( int i=0; i<maxit; i++){
		//fprintf( stderr, "Inner CG-Iteration %i\n", i+1);

		// r = r - alpha*z
		vectorScale( z, alpha, temp, N);
		vectorDifference( r, temp, r, N);
		
		/*float errSqr = dotProduct( r, r, N);
		if( isnan(errSqr)){
			fprintf( stderr, "ConjugateGradientSparse: ERROR is nan!\n");
			exit( 0);
		}
		myError = sqrt( errSqr);*/

		//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);

		myError = 0;
		for(int m=0; m<N; m++)
			if( myError < fabs(r[m]))
				myError = fabs(r[m]);

		if( myError < minerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, err, i+1);
			return;// i;
		}
		
		// B = (r'*z)/(w'*z);
		beta = dotProduct( r, z, N) / dotProduct( w, z, N);
		//if( beta<0.) beta = 0.;
		
		// w = -r + B*w;
		vectorScale( w, beta, w, N);
		vectorDifference( w, r, w, N);
		
		// z = A*w;
		sparseMatrixVectorProduct( A, w, z);
		
		// a = (r'*w)/(w'*z);
		alpha = dotProduct( r, w, N) / dotProduct( w, z, N);
		
		// x = x + a*w;
		vectorScale( w, alpha, temp, N);
		vectorSum( x, temp, x, N);
	}
	
	//fprintf( stderr, "error is still %e (=sqrt(%e)) < %e\n", myError, dotProduct( r, r, N), err);
	
	//return iterations;
}


void ConjugateGradientSparse( SparseMatrix *A, double *b, double *x, int maxit, double minerr)
// If A is symmetric positive matrix
{
	char filename[512];
	FILE *fp;

	sprintf( filename, "convergencePreCondCG_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();

	int N = A->dimI;

	// vectors
	double *r = v0d;
	double *w = v1d;
	double *z = v2d;

	// scalars
	double alpha;
	double beta;

	double *temp = v3d;

	//float err = 1e-10;
	double myError = 0.;
	//fprintf( stderr, "1\n");
	// initialization of residual vector: r
	sparseMatrixVectorProduct( A, x, temp);
	vectorDifference( b, temp, r, N);

	// w
	vectorScale( r, -1, w, N);
	//fprintf( stderr, "2\n");
	// z
	sparseMatrixVectorProduct( A, w, z);

	// alpha
	alpha = dotProduct( r, w, N) / dotProduct( w, z, N);

	// beta
	beta = 0.;

	// x
	vectorScale( w, alpha, temp, N);
	vectorSum( x, temp, x, N);

	for( int i=0; i<maxit; i++){
		//fprintf( stderr, "Inner CG-Iteration %i\n", i+1);

		// r = r - alpha*z
		vectorScale( z, alpha, temp, N);
		vectorDifference( r, temp, r, N);

		/*float errSqr = dotProduct( r, r, N);
		if( isnan(errSqr)){
			fprintf( stderr, "ConjugateGradientSparse: ERROR is nan!\n");
			exit( 0);
		}
		myError = sqrt( errSqr);*/

		//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);

		myError = 0;
		double myError2 = 0;
		for(int m=0; m<N; m++){
			if( myError < fabs(r[m]))
				myError = fabs(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( myError2 < fabs(exactSolution[m]-x[m]))
				myError2 = fabs(exactSolution[m]-x[m]);
#endif
		}

		fprintf( fp, "%i %e %e %e\n", i, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);

		if( myError < minerr){
			fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, minerr, i+1);
			return;// i;
		}

		// B = (r'*z)/(w'*z);
		beta = dotProduct( r, z, N) / dotProduct( w, z, N);
		//if( beta<0.) beta = 0.;

		// w = -r + B*w;
		vectorScale( w, beta, w, N);
		vectorDifference( w, r, w, N);
		//fprintf( stderr, "3\n");
		// z = A*w;
		sparseMatrixVectorProduct( A, w, z);

		// a = (r'*w)/(w'*z);
		alpha = dotProduct( r, w, N) / dotProduct( w, z, N);

		// x = x + a*w;
		vectorScale( w, alpha, temp, N);
		vectorSum( x, temp, x, N);
	}

	fprintf( stderr, "error is still %e (=sqrt(%e)) < %e\n", myError, dotProduct( r, r, N), minerr);

	fclose(fp);
	//return iterations;
}


/*float MaximalError( SparseMatrix *A, float *b, float *x)
{
	float max_error = 0.;
	for( int m=0; m<A->dimI; m++){
		if( )
	}
}*/

// ************ BiCGSTAB ************ //

// x = x - (omega*w + alpha*p)
void calculate_x( float *x, float omega, float *w, float alpha, float *p, int N)
{
	for( int i=0; i<N; i++)
		x[i] -=	(omega*w[i] + alpha*p[i]);
}

// r = r - (omega*v + alpha*s )
void calculate_r( float *r, float omega, float *v, float alpha, float *s, int N)
{
	for( int i=0; i<N; i++)
		r[i] -=	(omega*v[i] + alpha*s[i]);
}

// p = r - beta*(p - omega*s)
void calculate_p( float *x, float *y, float alpha, float *a, float beta, float *b, int N)
{
	for( int i=0; i<N; i++)
		x[i] = y[i] + alpha*a[i] + beta*b[i];
}

// p = r - beta*(p - omega*s)
void calculate_p1( float *x, float *y, float alpha, float *a, int N)
{
	for( int i=0; i<N; i++)
		x[i] = y[i] + alpha*a[i];
}

void SolveExplicit( SparseMatrix *A, float *b, float *x)
{
#pragma omp parallel for
	for( int i=0; i<A->dimI; i++)
		x[i] = b[i] / A->A[i][0];
}


// x = x - (omega*w + alpha*p)
void calculate_x( double *x, double omega, double *w, double alpha, double *p, int N)
{
#pragma omp parallel for
	for( int i=0; i<N; i++)
		x[i] -=	(omega*w[i] + alpha*p[i]);
}

// r = r - (omega*v + alpha*s )
void calculate_r( double *r, double omega, double *v, double alpha, double *s, int N)
{
#pragma omp parallel for
	for( int i=0; i<N; i++)
		r[i] -=	(omega*v[i] + alpha*s[i]);
}

// p = r - beta*(p - omega*s)
void calculate_p( double *x, double *y, double alpha, double *a, double beta, double *b, int N)
{
#pragma omp parallel for
	for( int i=0; i<N; i++)
		x[i] = y[i] + alpha*a[i] + beta*b[i];
}

// p = r - beta*(p - omega*s)
void calculate_p1( double *x, double *y, double alpha, double *a, int N)
{
#pragma omp parallel for
	for( int i=0; i<N; i++)
		x[i] = y[i] + alpha*a[i];
}


int SolveBiCGSTAB( SparseMatrix *A, float *b, float *x, int maxit, float minerr)
// If A is symmetric positive matrix
{
	int N = A->dimI;
	//fprintf( stderr, "test\n");
	// vectors
/*	float r0[N];
	float r[N];

	float p[N];

	float s[N];

	float v[N];
	float w[N];

	float temp[N];*/
	float *r0 = v0;
	float *r  = v1;

	float *p  = v2;

	float *s  = v3;

	float *v  = v4;
	float *w  = v5;

	float *temp = v6;
	//fprintf( stderr, "test end\n");
	
	// scalars
	float alpha;
	float beta;
	float omega;
	float rho;
	float sigma;
	
	
	// initialization of residual vector: r
	//fprintf( stderr, "r = Ax\n");
	sparseMatrixVectorProduct( A, x, temp);
	//fprintf( stderr, "...finished\n");
	vectorDifference( temp, b, r, N);
	/*for(int m=0; m<N; m++)
		if( isnan(r[m])){
			fprintf(stderr, "ERROR: r[%i] = %lf -> b[%i]=%lf x[%i]=%lf\n", m, r[m], m, b[m], m, x[m]);
			for( int jj=0; jj<A->sizeA[m]; jj++)
				fprintf(stderr, "A[%i][%i]=%lf\n", m, A->JA[m][jj], A->A[m][jj]);
			exit(0);
		}else
			fprintf(stderr, ", r[%i] = %.3lf", m, r[m]);
	 */

	// y, p, v
	for( int i=0; i<N; i++){
		r0[i] = r[i];
		p[i]  = r[i];
		//v[i] = 0.;
	}
	
	// rho
	rho = dotProduct( r, r0, N);
	
	// omega
	//omega = 1;
	
	// alpha
	//alpha = 1;

	int i=0;
	for( ; i<maxit; i++){
		//fprintf( stderr, "Inner CG-Iteration %i\n", i+1);

		// v = Ap
		//fprintf( stderr, "s = Ap\n");
		sparseMatrixVectorProduct( A, p, s);
		//fprintf( stderr, "...finished\n");
		/*for(int m=0; m<N; m++)
			if( isnan(s[m])){
				fprintf(stderr, "ERROR: s[%i] = %lf\n", m, s[m]); exit(0);
			}else
				fprintf(stderr, ", s[%i] = %.3lf", m, s[m]);
		 */

		// sigma
		//beta = rho;
		sigma = dotProduct( s, r0, N);
		if(sigma == 0){
			fprintf( stderr, "sigma = %lf, it = %i\n", sigma, i+1);
			return i;
		}
		
		// alpha 
		alpha = rho / sigma;
		if(isnan(alpha))
			alpha=0;
		
		// w = r - alpha*s
		/*vectorScale( s, alpha, temp, N);
		vectorDifference( r, temp, w, N);*/
		calculate_p1( w, r, -alpha, s, N);
		/*for(int m=0; m<N; m++)
			if( isnan(w[m])){
				fprintf(stderr, "ERROR: w[%i] = %lf\n", m, w[m]); exit(0);
			}else
				fprintf(stderr, ", w[%i] = %.3lf", m, w[m]);
		*/

		// v = Ap
		//fprintf( stderr, "v = Aw\n");
		sparseMatrixVectorProduct( A, w, v);
		//fprintf( stderr, "...finished\n");
		/*for(int m=0; m<N; m++)
			if( isnan(v[m])){
				fprintf(stderr, "ERROR: v[%i] = %lf\n", m, v[m]); exit(0);
			}else
				fprintf(stderr, ", v[%i] = %.3lf", m, v[m]);
		 */

		// omega
		omega = dotProduct( v, w, N) / dotProduct( v, v, N);
		if(isnan(omega))
			omega=0;
		//fprintf(stderr, "v*w = %lf, v*v = %lf\n", dotProduct( v, w, N), dotProduct( v, v, N));

		// x = x - (omega*w + alpha*p)
		/*vectorScale( w, omega, w, N);
		vectorScale( p, alpha, temp, N);
		vectorSum( w, temp, temp, N);
		vectorDifference( x, temp, x, N);*/
		calculate_x( x, omega, w, alpha, p, N);
		for(int m=0; m<N; m++)
			if( isnan(x[m])){
				fprintf(stderr, "ERROR: x[%i] = %lf, omega = %lf, alpha = %lf\n", m, x[m], omega, alpha); exit(0);
			}


		// r = r - (omega*v + alpha*s )
		/*vectorScale( v, omega, v, N);
		vectorScale( s, alpha, temp, N);
		vectorSum( v, temp, temp, N);
		vectorDifference( r, temp, r, N);*/
		calculate_r( r, omega, v, alpha, s, N);

		// rho
		beta = rho;
		rho = dotProduct( r, r0, N);

		// quadratic error
		//float myError = sqrt( dotProduct( r, r, N));
		
		// maximal error
		float myError = 0;
		for(int m=0; m<N; m++){
			if( myError < fabs(r[m]))
				myError = fabs(r[m]);
			/*if( isnan(r[m])){
				fprintf(stderr, "ERROR: r[%i] = %lf\n", m, r[m]); exit(0);
			}*/
		}
		
//		fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);
		if( myError < minerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, err, i+1);
			return i+1;
		}
		
		// new search direction
		// beta
		beta = rho / beta * alpha / omega;

		// p = r + beta*(p - omega*s)
		/*vectorScale( s, -omega, s, N);
		vectorSum( p, s, temp, N);
		vectorScale( temp, beta, temp, N);
		vectorSum( r, temp, p, N);*/
		calculate_p( p, r, beta, p, -beta*omega, s, N);

		
	}
	
	//fprintf( stderr, "error is still %e (=sqrt(%e)) < %e\n", myError, dotProduct( r, r, N), err);
	
	return i;
}

int SolveBiCGSTAB( SparseMatrix *A, double *b, double *x, int maxit, double minerr)
// If A is symmetric positive matrix
{
	char filename[512];
	FILE *fp;

	sprintf( filename, "convergenceBiCGSTAB_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();

	int N = A->dimI;
	//fprintf( stderr, "test\n");
	// vectors
/*	float r0[N];
	float r[N];

	float p[N];

	float s[N];

	float v[N];
	float w[N];

	float temp[N];*/
	double *r0 = v0d;
	double *r  = v1d;

	double *p  = v2d;

	double *s  = v3d;

	double *v  = v4d;
	double *w  = v5d;

	double *temp = v6d;
	//fprintf( stderr, "test end\n");

	// scalars
	double alpha;
	double beta;
	double omega;
	double rho;
	double sigma;


	// initialization of residual vector: r
	//fprintf( stderr, "r = Ax\n");
	sparseMatrixVectorProduct( A, x, temp);
	//fprintf( stderr, "...finished\n");
	vectorDifference( temp, b, r, N);

	// y, p, v
	for( int i=0; i<N; i++){
		r0[i] = r[i];
		p[i]  = r[i];
		//v[i] = 0.;
	}

	// rho
	rho = dotProduct( r, r0, N);

	// omega
	//omega = 1;

	// alpha
	//alpha = 1;

	int i=0;
	for( ; i<maxit; i++){
		//fprintf( stderr, "Inner CG-Iteration %i\n", i+1);

		// v = Ap
		//fprintf( stderr, "s = Ap\n");
		sparseMatrixVectorProduct( A, p, s);
		//fprintf( stderr, "...finished\n");

		// sigma
		//beta = rho;
		sigma = dotProduct( s, r0, N);
		if(sigma == 0){
			fprintf( stderr, "sigma = %lf, it = %i\n", sigma, i+1);
			return i;
		}

		// alpha
		alpha = rho / sigma;

		// w = r - alpha*s
		/*vectorScale( s, alpha, temp, N);
		vectorDifference( r, temp, w, N);*/
		calculate_p1( w, r, -alpha, s, N);

		// v = Ap
		//fprintf( stderr, "v = Aw\n");
		sparseMatrixVectorProduct( A, w, v);
		//fprintf( stderr, "...finished\n");

		// omega
		omega = dotProduct( v, w, N) / dotProduct( v, v, N);

		// x = x - (omega*w + alpha*p)
		/*vectorScale( w, omega, w, N);
		vectorScale( p, alpha, temp, N);
		vectorSum( w, temp, temp, N);
		vectorDifference( x, temp, x, N);*/
		calculate_x( x, omega, w, alpha, p, N);

		// r = r - (omega*v + alpha*s )
		/*vectorScale( v, omega, v, N);
		vectorScale( s, alpha, temp, N);
		vectorSum( v, temp, temp, N);
		vectorDifference( r, temp, r, N);*/
		calculate_r( r, omega, v, alpha, s, N);

		// rho
		beta = rho;
		rho = dotProduct( r, r0, N);

		// quadratic error
		//float myError = sqrt( dotProduct( r, r, N));

		// maximal error
		double myError = 0;
		double myError2 = 0;
		for(int m=0; m<N; m++){
			if( myError < fabs(r[m]))
				myError = fabs(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( myError2 < fabs(exactSolution[m]-x[m]))
				myError2 = fabs(exactSolution[m]-x[m]);
#endif
		}

		//fprintf( fp, "%i %e %e %e \n", i, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);

		//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);
		if( myError < minerr){
			fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, minerr, i+1);
			return i+1;
		}

		// new search direction
		// beta
		beta = rho / beta * alpha / omega;

		// p = r + beta*(p - omega*s)
		/*vectorScale( s, -omega, s, N);
		vectorSum( p, s, temp, N);
		vectorScale( temp, beta, temp, N);
		vectorSum( r, temp, p, N);*/
		calculate_p( p, r, beta, p, -beta*omega, s, N);


	}

	//fprintf( stderr, "error is still %e (=sqrt(%e)) < %e\n", myError, dotProduct( r, r, N), err);

	fclose(fp);
	return i;
	//return iterations;
}


/*void SolveBiCGSTAB( SparseMatrix *A, float *b, float *x, int maxit, float minerr)
// If A is symmetric positive matrix
{
	int N = A->dimI;
	
	// vectors
	float r[N];
	float p[N];
	float z[N];
	float y[N];
	float s[N];
	float v[N];
	float t[N];
	
	// scalars
	float alpha;
	float beta;
	float omega;
	float rho;
	
	float temp[N];
	
	//float err = 1e-10;
	float myError = 0.;
	
	// initialization of residual vector: r
	fprintf( stderr, "r = Ax\n");
	sparseMatrixVectorProduct( A, x, temp);
	fprintf( stderr, "...finished\n");
	vectorDifference( temp, b, r, N);
	
	// y, p, v
	for( int i=0; i<N; i++){
		y[i] = r[i];
		p[i] = 0.;
		v[i] = 0.;
	}
	
	// rho
	rho = 1;
	
	// omega
	omega = 1;
	
	// alpha
	alpha = 1;
	
	
	for( int i=0; i<maxit; i++){
		//fprintf( stderr, "Inner CG-Iteration %i\n", i+1);

		// rho
		beta = rho;
		rho = dotProduct( y, r, N);
		
		// beta
		fprintf( stderr, "beta (%lf) = rho_n (%lf) / rho_n-1 (%lf) * alpha (%lf) / omega_n (%lf) \n", 
		rho * alpha / (beta * omega), rho, beta, alpha, omega);
		beta = rho * alpha / (beta * omega);
		
		// p
		vectorScale( v, omega, temp, N);
		vectorSum( p, temp, temp, N);
		vectorScale( temp, beta, temp, N);
		vectorSum( r, temp, p, N);
		
		// v = Ap
		fprintf( stderr, "v = Ap\n");
		sparseMatrixVectorProduct( A, p, v);
		fprintf( stderr, "...finished\n");
		
		// ALPHA
		alpha = dotProduct( y, v, N);
		fprintf( stderr, "alpha (%lf)\n", alpha);
		
		// s
		vectorScale( v, alpha, temp, N);
		vectorDifference( r, temp, s, N);
		
		// t
		fprintf( stderr, "t = As\n");
		sparseMatrixVectorProduct( A, s, t);
		fprintf( stderr, "...finished\n");
		
		// omega
		omega = dotProduct( t, s, N) / dotProduct( t, t, N);
		fprintf( stderr, "omega (%lf) = <t,s> (%lf) / <t,t> (%lf)\n", omega, dotProduct( t, s, N) / dotProduct( t, t, N));
		// r
		vectorScale( t, omega, temp, N);
		vectorDifference( s, temp, r, N);
		
		// x
		vectorScale( s, omega, s, N);
		vectorScale( p, alpha, temp, N);
		vectorSum( s, temp, temp, N);
		vectorSum( x, temp, x, N);
		
		//if( sqrt( dotProduct( r, r, N)) < err){
		if( sqrt( dotProduct( r, r, N)) < minerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, err, i+1);
			return;// i;
		}
		
	}
	
	//fprintf( stderr, "error is still %e (=sqrt(%e)) < %e\n", myError, dotProduct( r, r, N), err);
	
	//return iterations;
}*/


void JacobiPreconditioner( SparseMatrix *A, SparseMatrix *M)
{
#pragma omp parallel for
	for( int i=0; i<A->dimI; i++){
		M->set( i,i, 1. / A->get(i,i));
		//M->set( i,i, A->get(i,i));
	}
}


void SuccessiveOverRelaxationPreconditioner( SparseMatrix *A, SparseMatrix *M)
{
	// (ixj) = (ixn) x (nxj)
	double omega = .01;
	
	// row
	for( int i=0; i<A->dimI; i++){
		//double factor = 0.;
		//column
		//for( int j=0; j<A->dimI; j++){
		for( int jj=0; jj<A->sizeA[i]; jj++){
			int j = A->JA[i][jj];
			double value = 0.;
			double ab = 0.;
			//double b = 0.;
			double c = 0.;
			
			// matrix product
			for( int pp=0; pp<A->sizeA[i]; pp++){
				int p = A->JA[i][pp];
				ab=0.;
				c=0.;
				// A*B
				if( i>p){
					// A*B = L * w/(2-w)D^-1
					ab = A->A[i][pp] * omega/(2.-omega) / A->get(p,p);
				}
				if( p==i){
					// A*B = D/w * w/(2-w)D^-1
					ab = A->A[i][pp]/(2.-omega) / A->get(p,p);
				}
				
				// C
				//if( ab!=0.){
				if( p<j){
					// C = U
					c = A->get(p,j) * omega/(2.-omega) / A->get(p,p);
				}
				if( p==j){
					// C = D/w
					c = A->get(p,p)/omega;
				}
				//}
				
				value += ab	* c;
			}
			
			if(value!=0.)
				M->set( i,j, value);
		}
		
		//M->set( i,i, A->get(i,i));
	}
//	fprintf( stderr, "SuccessiveOverRelaxationPreconditioner still not implimented\n");
//	exit( 0);
}


/*void SuccessiveOverRelaxationPreconditioner( SparseMatrix *A, SparseMatrix *M)
{
	// (ixj) = (ixn) x (nxj)
	double omega = .1;
	
	for( int i=0; i<A->dimI; i++){
		double factor = 0.;
		for( int jj=0; jj<A->sizeA[i]; jj++){
			int j = A->JA[i][jj];
			if( j<i){
				// L * w/(2-w)D^-1
				M->set( i,j, A->A[i][jj] * omega/(2.-omega) / A->get(j,j));
			}
			if( j==i){
				// D/w * w/(2-w)D^-1
				M->set( i,j, A->A[i][jj]/(2.-omega) / A->get(j,j));
			}
		}
		
		//M->set( i,i, A->get(i,i));
	}
//	fprintf( stderr, "SuccessiveOverRelaxationPreconditioner still not implimented\n");
//	exit( 0);
}*/


float PreconditionedConjugateGradientSparse( SparseMatrix *A, SparseMatrix *M, float *b, float *x, int maxit, float minerr)
{
	// vectors
	//float r[N];
	//float z[N];
	//float p[N];
	int N = A->dimI;
	
	float *r = v0;
	float *z = v1;
	float *p = v2;
	
	// scalars
	float alpha;
	float beta;
	float rho;
	
	//float temp[N];
	float *temp = v3;
	
	//float err = 1e-4;
	float myError = 0.;
	
	// initialization: 
	
	// r = b - A*x
	/*for( int i=0; i<N; i++)
		if( isnan(x[i])){
			fprintf( stderr, "x[%i] = %lf\n", i, x[i]);
			exit(0);
		}*/
	//fprintf( stderr, "r = b - A*x\n");
	sparseMatrixVectorProduct( A, x, temp);
	vectorDifference( b, temp, r, N);
	
	// z = M*r
	//fprintf( stderr, "z = M*r\n");
	sparseMatrixVectorProduct( M, r, z);

	// p = z
	vectorScale( z, 1, p, N);
	
	
	for( int i=0; i<maxit; i++){
		// rho = r*z
		rho = dotProduct( r, z, N);
		
		// alpha = r*z / p*A*p
		vectorSparseMatrixProduct( p, A, temp);
		alpha = rho / dotProduct( temp, p, N);
	
		// x = x + alpha*p
		vectorScale( p, alpha, temp, N);
		vectorSum( x, temp, x, N);
	
		// r = r - alpha*Ap
		sparseMatrixVectorProduct( A, p, temp);
		vectorScale( temp, alpha, temp, N);
		vectorDifference( r, temp, r, N);
		
		myError = sqrt( dotProduct( r, r, N));
		//if( sqrt( dotProduct( r, r, N)) < err){
		fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);

		if( myError < minerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, err, i+1);
			//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);
			return myError;
		}else if( myError > 1.){
			//fprintf( stderr, "WARNING: ConjugateGradientSparse reached error %e > 1 after %i iterations\n", myError, i+1);	
		}
		
		// z = M*r
		sparseMatrixVectorProduct( M, r, z);
		
		// beta = r*z/rho
		beta = dotProduct( r, z, N) / rho;
		
		// p = z + beta*p
		vectorScale( p, beta, p, N);
		vectorSum( z, p, p, N);
	}
	
	fprintf( stderr, "error is still %e < %e\n", myError, minerr);
	return myError;
}





double PreconditionedConjugateGradientSparse( SparseMatrix *A, SparseMatrix *M, double *b, double *x, int maxit, float minerr)
{
	char filename[512];
	FILE *fp;

	sprintf( filename, "convergencePreconditionedCG_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();

	// vectors
	//float r[N];
	//float z[N];
	//float p[N];
	int N = A->dimI;

	double *r = v0d;
	double *z = v1d;
	double *p = v2d;

	// scalars
	double alpha;
	double beta;
	double rho;

	//float temp[N];
	double *temp = v3d;

	//float err = 1e-4;
	double myError = 0.;

	// initialization:

	// r = b - A*x
	/*for( int i=0; i<N; i++)
		if( isnan(x[i])){
			fprintf( stderr, "x[%i] = %lf\n", i, x[i]);
			exit(0);
		}*/
	//fprintf( stderr, "r = b - A*x\n");
	sparseMatrixVectorProduct( A, x, temp);
	vectorDifference( b, temp, r, N);

	// z = M*r
	//fprintf( stderr, "z = M*r\n");
	sparseMatrixVectorProduct( M, r, z);
	//SolveSuccessiveOverRelaxationMethod( A, r, z, 50, 1e-6, 1.5);
	//SolveJacobiMethod( A, r, z, 50, 1e-6);

	// p = z
	vectorScale( z, 1, p, N);


	for( int i=0; i<maxit; i++){
		// rho = r*z
		rho = dotProduct( r, z, N);

		// alpha = r*z / p*A*p
		vectorSparseMatrixProduct( p, A, temp);
		alpha = rho / dotProduct( temp, p, N);

		// x = x + alpha*p
		vectorScale( p, alpha, temp, N);
		vectorSum( x, temp, x, N);

		// r = r - alpha*Ap
		sparseMatrixVectorProduct( A, p, temp);
		vectorScale( temp, alpha, temp, N);
		vectorDifference( r, temp, r, N);

		myError = sqrt( dotProduct( r, r, N));
		//if( sqrt( dotProduct( r, r, N)) < err){
		fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);

		myError = 0;
		double myError2 = 0;
		for(int m=0; m<N; m++){
			if( myError < fabs(r[m]))
				myError = fabs(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( myError2 < fabs(exactSolution[m]-x[m]))
				myError2 = fabs(exactSolution[m]-x[m]);
#endif
		}

		//fprintf( fp, "%i %e %e %e\n", i, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);

		if( myError < minerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, err, i+1);
			//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);
			return myError;
		}else if( myError > 1.){
			//fprintf( stderr, "WARNING: ConjugateGradientSparse reached error %e > 1 after %i iterations\n", myError, i+1);
		}

		// z = M*r
		sparseMatrixVectorProduct( M, r, z);
		//SolveSuccessiveOverRelaxationMethod( A, r, z, 50, 1e-6, 1.5);
		//SolveJacobiMethod( A, r, z, 50, 1e-6);

		// beta = r*z/rho
		beta = dotProduct( r, z, N) / rho;

		// p = z + beta*p
		vectorScale( p, beta, p, N);
		vectorSum( z, p, p, N);
	}
	fclose(fp);
	fprintf( stderr, "error is still %e < %e\n", myError, minerr);
	return myError;
}





/*void PreconditionedConjugateGradientCoupledSparse( SparseMatrix *A, SparseMatrix *M, float *b, float *x, int iterations)
{
	// vectors
	float r[N];
	float z[N];
	float p[N];
	
	// scalars
	float alpha;
	float beta;
	float rho;
	
	float temp[N];
	
	//float err = 1e-4;
	float myError = 0.;
	
	// initialization: 
	
	// r = b - A*x
	sparseMatrixVectorProduct( A, x, temp, N);
	vectorDifference( b, temp, r, N);
	
	// z = M*r
	sparseMatrixVectorProduct( M, r, z, N);

	// p = z
	vectorScale( z, 1, p, N);
	
	
	for( int i=0; i<iterations; i++){
		// rho = r*z
		rho = dotProduct( r, z, N);
		
		// alpha = r*z / p*A*p
		vectorsparseMatrixProduct( p, A, temp, N);
		alpha = rho / dotProduct( temp, p, N);
	
		// x = x + alpha*p
		vectorScale( p, alpha, temp, N);
		vectorSum( x, temp, x, N);
	
		// r = r - alpha*Ap
		sparseMatrixVectorProduct( A, p, temp, N);
		vectorScale( temp, alpha, temp, N);
		vectorDifference( r, temp, r, N);
		
		myError = sqrt( dotProduct( r, r, N));
		//if( sqrt( dotProduct( r, r, N)) < err){
		//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);

		if( myError < CG_CRITICAL_ERROR){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, err, i+1);
			return;
		}else if( myError > 1.){
			//fprintf( stderr, "WARNING: ConjugateGradientSparse reached error %e > 1 after %i iterations\n", myError, i+1);	
		}
		
		// z = M*r
		sparseMatrixVectorProduct( M, r, z, N);
		
		// beta = r*z/rho;
		beta = dotProduct( r, z, N) / rho;
		
		// p = z + beta*p
		vectorScale( p, beta, p, N);
		vectorSum( z, p, p, N);
	}
	
//	fprintf( stderr, "error is still %e < %e\n", myError, err);
}*/



void sparseMatrixVectorProduct( SparseMatrix *sA, float *b, float *x)
{
	//x = A*b
	// (n x m) * (m x p) = (n x p)
	// (i x j) * (j x 1) = (i x 1)
	for( int i=0; i<sA->dimI; i++){
		x[i] = 0.;
		for( int jj=0; jj<sA->sizeA[i]; jj++){
			int j = sA->JA[i][jj];
			x[i] += sA->A[i][jj] * b[j];
#if DEBUG > 0
			if( isnan(x[i])){
				fprintf( stderr, "nan occures in sparseMatrixVectorProduct\n");
				//fprintf( stderr, "vector b[%i] = %lf\n", j, b[j]);
				x[i] = 0.;
				for( int jj=0; jj<sA->sizeA[i]; jj++){
					j = sA->JA[i][jj];
					x[i] += sA->A[i][jj] * b[j];
					fprintf( stderr, "matrix A[%i][%i] = %lf, vector b[%i] = %lf, x[%i] = %lf\n", i, j, sA->A[i][jj], j, b[j], i, x[i]);
				}
				exit( 0);
			}
#endif
		}	
	}
}
 
void sparseMatrixVectorProduct( SparseMatrix *sA, double *b, double *x)
{
	//x = A*b
	// (n x m) * (m x p) = (n x p)
	// (i x j) * (j x 1) = (i x 1)
#pragma omp parallel for //shared(sA,b,x) private(i)
	for( int i=0; i<sA->dimI; i++){
		x[i] = 0.;
		for( int jj=0; jj<sA->sizeA[i]; jj++){
			int j = sA->JA[i][jj];
			x[i] += sA->A[i][jj] * b[j];
#if DEBUG > 0
			if( isnan(x[i])){
				fprintf( stderr, "nan occures in sparseMatrixVectorProduct\n");
				//fprintf( stderr, "vector b[%i] = %lf\n", j, b[j]);
				x[i] = 0.;
				for( int jj=0; jj<sA->sizeA[i]; jj++){
					j = sA->JA[i][jj];
					x[i] += sA->A[i][jj] * b[j];
					fprintf( stderr, "matrix A[%i][%i] = %lf, vector b[%i] = %lf, x[%i] = %lf\n", i, j, sA->A[i][jj], j, b[j], i, x[i]);
				}
				exit( 0);
			}
#endif
		}
	}
}

// NEW!!!
void vectorSparseMatrixProduct( float *b, SparseMatrix *sA, float *x)
{
	//x = b*A
	// (n x m) * (m x p) = (n x p)
	// (1 x i) * (i x j) = (1 x j)
	//int jj = 0;
	for( int j=0; j<sA->dimI; j++)
		x[j] = 0.;

	for( int i=0; i<sA->dimI; i++){
		for( int jj=0; jj<sA->sizeA[i]; jj++){
			int j = sA->JA[i][jj];
			x[j] += b[i] * sA->A[i][jj];
#if DEBUG > 0
			if( isnan(x[j])){
				fprintf( stderr, "nan occures in vectorsparseMatrixProduct\n");
				fprintf( stderr, "vector b[%i] = %lf\n", i, b[i]);
				for( int jj=0; jj<sA->sizeA[i]; jj++){
					j = sA->JA[i][jj];
					fprintf( stderr, "matrix A[%i][%i] = %lf\n", i, j, sA->A[i][jj]);
				}
				exit( 0);
			}
#endif
		}	
	}
}
 
void vectorSparseMatrixProduct( double *b, SparseMatrix *sA, double *x)
{
	//x = b*A
	// (n x m) * (m x p) = (n x p)
	// (1 x i) * (i x j) = (1 x j)
	//int jj = 0;
	for( int j=0; j<sA->dimI; j++)
		x[j] = 0.;

	for( int i=0; i<sA->dimI; i++){
		for( int jj=0; jj<sA->sizeA[i]; jj++){
			int j = sA->JA[i][jj];
			x[j] += b[i] * sA->A[i][jj];
#if DEBUG > 0
			if( isnan(x[j])){
				fprintf( stderr, "nan occures in vectorsparseMatrixProduct\n");
				fprintf( stderr, "vector b[%i] = %lf\n", i, b[i]);
				for( int jj=0; jj<sA->sizeA[i]; jj++){
					j = sA->JA[i][jj];
					fprintf( stderr, "matrix A[%i][%i] = %lf\n", i, j, sA->A[i][jj]);
				}
				exit( 0);
			}
#endif
		}
	}
}

void LUdecomposition( SparseMatrix *sA, SparseMatrix *L, SparseMatrix *U, float *b)
{

}





void PreconditionJacobi( SparseMatrix *A, double *b)
{
#pragma omp parallel for
	for( int i=0; i<A->dimI; i++){
		double Aii = A->get(i,i);
		for( int jj=0; jj<A->sizeA[i]; jj++){
			//int j = A->JA[i][jj];
			A->A[i][jj] /= Aii;
		}
		b[i] /= Aii;
	}
}

void PreconditionSOR( SparseMatrix *A, double *b, double omega)
{
	double _omega = 1.-omega;
	for( int i=0; i<A->dimI; i++){
		double Aii = 1.;
		double sigma = 0.;
		for( int jj=0; jj<A->sizeA[i]; jj++){
			int j=A->JA[i][jj];
			if( j!=i)
				sigma += A->A[i][jj];
			else
				Aii = A->A[i][jj];
		}

		for( int jj=0; jj<A->sizeA[i]; jj++){
			A->A[i][jj] *= (_omega*sigma + (b[i]-sigma) * omega / Aii);
		}
		b[i] *= (_omega*sigma + (b[i]-sigma) * omega / Aii);
	}
}
