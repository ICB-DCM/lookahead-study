#include <stdio.h>
#include "Substrate.h"
#include "Mathematix.h"
#include "finiteDifferences.h"
//#include <boost/lambda/lambda.hpp>

#define DX 13.

#define EXPLICIT        0
#define IMPLICIT        1
#define CRANK_NICHOLSON 2
#define ADI             3

//#define DIFFUSION_METHODE           EXPLICIT
//#define DIFFUSION_METHODE           IMPLICIT
//#define DIFFUSION_METHODE           CRANK_NICHOLSON
#define DIFFUSION_METHODE           ADI
#define DIFFUSION_REACTION_SEPARAT FALSE
//#define DIFFUSION_REACTION_SEPARAT TRUE

//#define DT 0.01


float **A;
float *b;//[20];
float *x;//[20];
int ADI_direction = 1;

void initMatrix( VoronoiDiagram *voronoiDiagram/*, int imax, int iimax*/)
{
	int matrixDim = voronoiDiagram->xN[0];

	for( int dim=1; dim<DIMENSIONS; dim++)
		if(voronoiDiagram->xN[dim]>matrixDim)
			matrixDim = voronoiDiagram->xN[dim];

	A = newMatrix( matrixDim, matrixDim);
	
	b = (float*) malloc( matrixDim * sizeof(float));
	x = (float*) malloc( matrixDim * sizeof(float));
/*	for( int i=0; i<imax; i++)
		for( int ii=0; ii<iimax; ii++){
			A[i][ii] = voronoiDiagram->voronoiCells[i + ii*20];
		}
*/}
double factor = 1.;

/*void setupMatrixOxygenX( VoronoiDiagram *voronoiDiagram, int iset, double timeStep)
{
	double r = factor*Oxygen_Diffusion * timeStep/(DX*DX);
	
	for( int m=0; m<20; m++)
		for( int n=0; n<20; n++)
			A[m][n] = 0;

	for( int m=0; m<20; m++){
		A[m][m]   = 1;

		if( m>0){
			A[m][m-1] = -r;
			A[m][m]  += r;
		}
		if( m<19){
			A[m][m]  += r;
			A[m][m+1] = -r;
		}
		
		b[m] = voronoiDiagram->voronoiCells[iset + m*20]->oxygen - timeStep * voronoiDiagram->voronoiCells[iset + m*20]->oxygen * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[iset + m*20]);
	}
}

void setupMatrixOxygenY( VoronoiDiagram *voronoiDiagram, int iiset, double timeStep)
{
	double r = factor*Oxygen_Diffusion * timeStep/(DX*DX);
	
	for( int m=0; m<20; m++)
		for( int n=0; n<20; n++)
			A[m][n] = 0;

	for( int m=0; m<20; m++){
		A[m][m]   = 1;

		if( m>0){
			A[m][m-1] = -r;
			A[m][m]  += r;
		}
		if( m<19){
			A[m][m]  += r;
			A[m][m+1] = -r;
		}
		
		b[m] = voronoiDiagram->voronoiCells[m + iiset*20]->oxygen - timeStep * voronoiDiagram->voronoiCells[m+iiset*20]->oxygen * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[m+iiset*20]);
	}
}*/

void setupMatrix( VoronoiDiagram *voronoiDiagram, int base_index, double timeStep, char molecule, char direction)
{
	double r;
	//double consumption;
	
	if( molecule == 'o'){
		r = factor*Oxygen_Diffusion * timeStep/(DX*DX);
	}else{
		r = factor*Glucose_Diffusion * timeStep/(DX*DX);	
	}
	
	for( int m=0; m<voronoiDiagram->xN[(int)direction]/*20*/; m++)
		for( int n=0; n<voronoiDiagram->xN[(int)direction]/*20*/; n++)
			A[m][n] = 0;

	for( int m=0; m<voronoiDiagram->xN[(int)direction]/*20*/; m++){
		if( m>0 && m<voronoiDiagram->xN[(int)direction]-1){
			A[m][m-1] = -r;
			A[m][m]   =  1 + 2*r;
			A[m][m+1] = -r;
		}
		else{
			A[m][m]   = 1;
		}
		
		int di = (int)pow(voronoiDiagram->xN[(int)direction],direction);
		int index = base_index + (int)(m*di);
		//fprintf(stderr, "m=%i, dir=%i,  %i + %i => %i (< %i)\n", m, direction, base_index, (int)(m*pow(20,direction)), index, (int)pow(20,3));

		if( molecule == 'o')
			b[m] = voronoiDiagram->voronoiCells[index]->oxygen
#if DIFFUSION_REACTION_SEPARAT == FALSE
			     - timeStep/DIMENSIONS * voronoiDiagram->voronoiCells[index]->oxygen * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[index])
#endif
			    ;
		else
			b[m] = voronoiDiagram->voronoiCells[index]->glucose
#if DIFFUSION_REACTION_SEPARAT == FALSE
			     - timeStep/DIMENSIONS * voronoiDiagram->voronoiCells[index]->glucose * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[index])
#endif
			    ;
		/*for( int n=0; n<20; n++)
			fprintf( stderr, "%lf ", A[m][n]);
		fprintf( stderr, "| %lf\n", b[m]);*/
	}
}

void setupMatrixADI( VoronoiDiagram *voronoiDiagram, int base_index, double timeStep, char molecule, char direction)
{
	double r;
	//double consumption;
	
	if( molecule == 'o'){
		r = factor*Oxygen_Diffusion * timeStep/(DX*DX);
	}else{
		r = factor*Glucose_Diffusion * timeStep/(DX*DX);	
	}
	
	for( int m=0; m<voronoiDiagram->xN[(int)direction]/*20*/; m++)
		for( int n=0; n<voronoiDiagram->xN[(int)direction]/*20*/; n++)
			A[m][n] = 0;

	for( int m=0; m<voronoiDiagram->xN[(int)direction]/*20*/; m++){
		if( m>0 && m<voronoiDiagram->xN[(int)direction]-1){
			A[m][m-1] = -r;
			A[m][m]   =  1 + 2*r;
			A[m][m+1] = -r;
		}
		else{
			A[m][m]   = 1;
		}
		/*if( m>0){
			A[m][m-1] = -r;
			A[m][m]  += r;
		}
		if( m<voronoiDiagram->xN[(int)direction]-1){
			A[m][m]  += r;
			A[m][m+1] = -r;
		}*/
		
		int di = (int)pow(voronoiDiagram->xN[(int)direction],direction);
		int index = base_index + (int)(m*di);
		//fprintf(stderr, "m=%i, dir=%i,  %i + %i => %i (< %i)\n", m, direction, base_index, (int)(m*pow(20,direction)), index, (int)pow(20,3));

		if( molecule == 'o')
			b[m] = voronoiDiagram->voronoiCells[index]->oxygen
#if DIFFUSION_REACTION_SEPARAT == FALSE
			     - timeStep/DIMENSIONS * voronoiDiagram->voronoiCells[index]->oxygen * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[index])
#endif
			    ;
		else
			b[m] = voronoiDiagram->voronoiCells[index]->glucose
#if DIFFUSION_REACTION_SEPARAT == FALSE
			     - timeStep/DIMENSIONS * voronoiDiagram->voronoiCells[index]->glucose * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[index])
#endif
			    ;
		/*for( int n=0; n<20; n++)
			fprintf( stderr, "%lf ", A[m][n]);
		fprintf( stderr, "| %lf\n", b[m]);*/
	}
}

void setupMatrixCrankNicholson( VoronoiDiagram *voronoiDiagram, int base_index, double timeStep, char molecule, char direction)
{
	double r;
	//double consumption;
	
	if( molecule == 'o'){
		r = Oxygen_Diffusion * timeStep/(DX*DX);
	}else{
		r = Glucose_Diffusion * timeStep/(DX*DX);	
	}
	
	for( int m=0; m<20; m++)
		for( int n=0; n<20; n++)
			A[m][n] = 0;

	for( int m=0; m<20; m++){
		//A[m][m]   = 2;

		if( m>0){
			A[m][m-1] = -r;
			A[m][m]  += 1+r;
		}
		if( m<19){
			A[m][m]  += 1+r;
			A[m][m+1] = -r;
		}
		
		int index, index_minus, index_plus;
		
		index = base_index + (int)(m*pow(20,direction));
		index_plus  = base_index + (int)((m+1)*pow(20,direction));
		index_minus = base_index + (int)((m-1)*pow(20,direction));
		
		
		double factor, factor_minus, factor_plus;
		factor       = 1-r;
		factor_plus  = 0;
		factor_minus = 0;

#if DIFFUSION_REACTION_SEPARAT == FALSE
		if( molecule == 'o'){
			factor       -= GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[index]);
			if( m<19)factor_plus  = 1 - GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[index_plus]);
			if( m>0) factor_minus = 1 - GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[index_minus]);
		}else{
			factor       -= GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[index]);
			if( m<19)factor_plus  = 1 - GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[index_plus]);
			if( m>0) factor_minus = 1 - GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[index_minus]);
		}
#endif
		b[m] = 0.;
		if( molecule == 'o'){
			if( m>0){
				b[m] += timeStep/DIMENSIONS * factor * voronoiDiagram->voronoiCells[index]->oxygen + timeStep/DIMENSIONS * factor_minus * voronoiDiagram->voronoiCells[index_minus]->oxygen;
			}
			if( m<19){
				b[m] += timeStep/DIMENSIONS * factor * voronoiDiagram->voronoiCells[index]->oxygen + timeStep/DIMENSIONS * factor_plus * voronoiDiagram->voronoiCells[index_plus]->oxygen;
			}
		}else{
			if( m>0){
				b[m] += timeStep/DIMENSIONS * factor * voronoiDiagram->voronoiCells[index]->glucose + timeStep/DIMENSIONS * factor_minus * voronoiDiagram->voronoiCells[index_minus]->glucose;
			}
			if( m<19){
				b[m] += timeStep/DIMENSIONS * factor * voronoiDiagram->voronoiCells[index]->glucose +timeStep/DIMENSIONS *  factor_plus * voronoiDiagram->voronoiCells[index_plus]->glucose;
			}
		}
		/*for( int n=0; n<20; n++)
			fprintf( stderr, "%lf ", A[m][n]);
		fprintf( stderr, "| %lf\n", b[m]);*/
	}
}

double secondDerivationOxygen( VoronoiCell * cell, VoronoiDiagram *voronoiDiagram)
{
	double temp = 0.;

/*	temp = - ((double)cell->countNeighborCells) * cell->oxygen;

	for( int i=0; i<cell->countNeighborCells; i++)
		temp += cell->neighborCells[i]->oxygen;
		
	return temp/(DX*DX) * 4./(double)cell->countNeighborCells;
*/
	int i   = 0, 
	    ii  = 0,
	    iii = 0;
	
	int index = cell->index;
	
#if DIMENSIONS > 2
	iii = index / (20*20);
	index = index - iii*20*20;
#endif
#if DIMENSIONS > 1
	ii = index / 20;
	index = index - ii*20;
#endif
	i = index;
	//fprintf( stderr, "%i: (%i,%i,%i) => %i\n", cell->index, i, ii, iii, i + 20*ii + 20*20*iii);
	
	if( i>0)
		temp += voronoiDiagram->voronoiCells[i-1 + 20*ii + 20*20*iii]->oxygen - cell->oxygen;
	if( i<19)
		temp += voronoiDiagram->voronoiCells[i+1 + 20*ii + 20*20*iii]->oxygen - cell->oxygen;
#if DIMENSIONS > 1
	if( ii>0)
		temp += voronoiDiagram->voronoiCells[i + 20*(ii-1) + 20*20*iii]->oxygen - cell->oxygen;
	if( ii<19)
		temp += voronoiDiagram->voronoiCells[i + 20*(ii+1) + 20*20*iii]->oxygen - cell->oxygen;
#endif
#if DIMENSIONS > 2
	if( iii>0)
		temp += voronoiDiagram->voronoiCells[i + 20*ii + 20*20*(iii-1)]->oxygen - cell->oxygen;
	if( iii<19)
		temp += voronoiDiagram->voronoiCells[i + 20*ii + 20*20*(iii+1)]->oxygen - cell->oxygen;
#endif
	
	return temp/(DX*DX);
}

double secondDerivationGlucose( VoronoiCell * cell, VoronoiDiagram *voronoiDiagram)
{
	double temp = 0.;

/*	temp = - ((double)cell->countNeighborCells) * cell->glucose;

	for( int i=0; i<cell->countNeighborCells; i++)
		temp += cell->neighborCells[i]->glucose;
		
	return temp/(DX*DX) * 4./(double)cell->countNeighborCells;
*/
	int i   = 0, 
	    ii  = 0,
	    iii = 0;
	
	int index = cell->index;
	
#if DIMENSIONS > 2
	iii = index / (20*20);
	index = index - iii*20*20;
#endif
#if DIMENSIONS > 1
	ii = index / 20;
	index = index - ii*20;
#endif
	i = index;
	//fprintf( stderr, "%i: (%i,%i,%i) => %i\n", cell->index, i, ii, iii, i + 20*ii + 20*20*iii);
	
	if( i>0)
		temp += voronoiDiagram->voronoiCells[i-1 + 20*ii + 20*20*iii]->glucose - cell->glucose;
	if( i<19)
		temp += voronoiDiagram->voronoiCells[i+1 + 20*ii + 20*20*iii]->glucose - cell->glucose;
#if DIMENSIONS > 1
	if( ii>0)
		temp += voronoiDiagram->voronoiCells[i + 20*(ii-1) + 20*20*iii]->glucose - cell->glucose;
	if( ii<19)
		temp += voronoiDiagram->voronoiCells[i + 20*(ii+1) + 20*20*iii]->glucose - cell->glucose;
#endif
#if DIMENSIONS > 2
	if( iii>0)
		temp += voronoiDiagram->voronoiCells[i + 20*ii + 20*20*(iii-1)]->glucose - cell->glucose;
	if( iii<19)
		temp += voronoiDiagram->voronoiCells[i + 20*ii + 20*20*(iii+1)]->glucose - cell->glucose;
#endif
	
	return temp/(DX*DX);
}



/*#if (DIFFUSION_METHODE != ADI)

double UpdateSystem( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference){
	
	double time;
	
	int direction = (int)((myRand()*DIMENSIONS)/2.);
	
	for( time = 0; time+timeStep <= timeDifference; time += timeStep){
		
#if DIFFUSION_METHODE == EXPLICIT

		//fprintf( stderr, "EXPLICIT\n");
		// Finite Differences: explicit method	
		for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
			
			voronoiDiagram->voronoiCells[i]->doxygen  = timeStep * (Oxygen_Diffusion  * secondDerivationOxygen(  voronoiDiagram->voronoiCells[i], voronoiDiagram) // diffusion
	#if DIFFUSION_REACTION_SEPARAT == FALSE
			                                          - voronoiDiagram->voronoiCells[i]->oxygen  * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[i]) // reaction
	#endif
			                                          );
			voronoiDiagram->voronoiCells[i]->dglucose = timeStep * (Glucose_Diffusion * secondDerivationGlucose( voronoiDiagram->voronoiCells[i], voronoiDiagram) // diffusion
	#if DIFFUSION_REACTION_SEPARAT == FALSE
			                                          - voronoiDiagram->voronoiCells[i]->glucose * GiveMeTheGlucoseRate( voronoiDiagram->voronoiCells[i]) // reaction
	#endif
			                                          );			                                          
		}

#endif
		
#if (DIFFUSION_METHODE == IMPLICIT) || (DIFFUSION_METHODE == CRANK_NICHOLSON) || (DIFFUSION_METHODE == ADI)
		// Finite Differences: alternating direction implicit (ADI)
		double weight = 1.;
		// diffusion in x-direction
		for( int iset=0; iset<20; iset++){
			//setupMatrixOxygenX( voronoiDiagram, iset, timeStep);
	#if   DIFFUSION_METHODE == CRANK_NICHOLSON
			setupMatrixCrankNicholson( voronoiDiagram, iset, timeStep, 'o', 0);
			//fprintf( stderr, "CRANK_NICHOLSON\n");
	#elif DIFFUSION_METHODE == ADI
			setupMatrixADI( voronoiDiagram, iset, timeStep, 'o', direction);
	#else
			//fprintf( stderr, "ADI\n");
			setupMatrix( voronoiDiagram, iset, 0.5*timeStep, 'o', 0);
	#endif
			solveLinearSystem( A, b, x, 20);
			for( int ii=0; ii<20; ii++)
				voronoiDiagram->voronoiCells[iset + ii*20]->doxygen = weight*(b[ii] - voronoiDiagram->voronoiCells[iset + ii*20]->oxygen);

	#if DIFFUSION_METHODE == CRANK_NICHOLSON
			setupMatrixCrankNicholson( voronoiDiagram, iset, timeStep, 'g', 0);
	#elif DIFFUSION_METHODE == ADI
			setupMatrixADI( voronoiDiagram, iset, timeStep, 'g', direction);
	#else
			setupMatrix( voronoiDiagram, iset, 0.5*timeStep, 'g', 0);
	#endif
			solveLinearSystem( A, b, x, 20);
			for( int ii=0; ii<20; ii++)
				voronoiDiagram->voronoiCells[iset + ii*20]->dglucose = weight*(b[ii] - voronoiDiagram->voronoiCells[iset + ii*20]->glucose);
		}
	#if DIFFUSION_METHODE == ADI
		for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
			if( voronoiDiagram->voronoiCells[i]->position[0]>voronoiDiagram->xMin[0]+1
			 && voronoiDiagram->voronoiCells[i]->position[0]<voronoiDiagram->xMax[0]-1 
			 && voronoiDiagram->voronoiCells[i]->position[1]>voronoiDiagram->xMin[1]+1
			 && voronoiDiagram->voronoiCells[i]->position[1]<voronoiDiagram->xMax[1]-1
		#if DIMENSIONS > 2
			 && voronoiDiagram->voronoiCells[i]->position[2]>voronoiDiagram->xMin[2]+1
			 && voronoiDiagram->voronoiCells[i]->position[2]<voronoiDiagram->xMax[2]-1 
		#endif
			){	
				voronoiDiagram->voronoiCells[i]->oxygen  += voronoiDiagram->voronoiCells[i]->doxygen;
				if( voronoiDiagram->voronoiCells[i]->oxygen < 0.) voronoiDiagram->voronoiCells[i]->oxygen = 0.;
				voronoiDiagram->voronoiCells[i]->glucose += voronoiDiagram->voronoiCells[i]->dglucose;
				if( voronoiDiagram->voronoiCells[i]->glucose < 0.) voronoiDiagram->voronoiCells[i]->glucose = 0.;
			}
		}
		direction = (direction+1)%DIMENSIONS;
	#endif

		// diffusion in y-direction
		for( int iiset=0; iiset<20; iiset++){
//			setupMatrixOxygenY( voronoiDiagram, iiset, timeStep);
	#if DIFFUSION_METHODE == CRANK_NICHOLSON
			setupMatrixCrankNicholson( voronoiDiagram, iiset, timeStep, 'o', 1);
	#elif DIFFUSION_METHODE == ADI
			setupMatrixADI( voronoiDiagram, iiset, timeStep, 'o', direction);
	#else
			setupMatrix( voronoiDiagram, iiset, 0.5*timeStep, 'o', 1);
	#endif
			solveLinearSystem( A, b, x, 20);
			for( int i=0; i<20; i++){
				voronoiDiagram->voronoiCells[i+iiset*20]->doxygen += weight*(b[i] - voronoiDiagram->voronoiCells[i+iiset*20]->oxygen);
			}

	#if DIFFUSION_METHODE == CRANK_NICHOLSON
			setupMatrixCrankNicholson( voronoiDiagram, iiset, timeStep, 'g', 1);
	#elif DIFFUSION_METHODE == ADI
			setupMatrixADI( voronoiDiagram, iiset, timeStep, 'g', direction);
	#else
			setupMatrix( voronoiDiagram, iiset, 0.5*timeStep, 'g', 1);
	#endif
			solveLinearSystem( A, b, x, 20);
			for( int i=0; i<20; i++){
				voronoiDiagram->voronoiCells[i+iiset*20]->dglucose += weight*(b[i] - voronoiDiagram->voronoiCells[i+iiset*20]->glucose);
			}
		}


#endif

		// Update
		for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
			if( voronoiDiagram->voronoiCells[i]->position[0]>voronoiDiagram->xMin[0]+1
			 && voronoiDiagram->voronoiCells[i]->position[0]<voronoiDiagram->xMax[0]-1 
			 && voronoiDiagram->voronoiCells[i]->position[1]>voronoiDiagram->xMin[1]+1
			 && voronoiDiagram->voronoiCells[i]->position[1]<voronoiDiagram->xMax[1]-1
#if DIMENSIONS > 2
			 && voronoiDiagram->voronoiCells[i]->position[2]>voronoiDiagram->xMin[2]+1
			 && voronoiDiagram->voronoiCells[i]->position[2]<voronoiDiagram->xMax[2]-1 
#endif
			){	
				voronoiDiagram->voronoiCells[i]->oxygen  += voronoiDiagram->voronoiCells[i]->doxygen
#if DIFFUSION_REACTION_SEPARAT == TRUE
				                                          - timeStep * voronoiDiagram->voronoiCells[i]->oxygen  * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[i])
#endif
				                                          ;
				if( voronoiDiagram->voronoiCells[i]->oxygen < 0.) voronoiDiagram->voronoiCells[i]->oxygen = 0.;
				voronoiDiagram->voronoiCells[i]->glucose += voronoiDiagram->voronoiCells[i]->dglucose
#if DIFFUSION_REACTION_SEPARAT == TRUE
				                                          - timeStep * voronoiDiagram->voronoiCells[i]->glucose * GiveMeTheGlucoseRate( voronoiDiagram->voronoiCells[i])
#endif
				                                          ;
				if( voronoiDiagram->voronoiCells[i]->glucose < 0.) voronoiDiagram->voronoiCells[i]->glucose = 0.;
			}
		}
	}

	return timeDifference - time;
}

#else*/

double UpdateSystem( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference){
	
	double time;
	
	//int direction = (int)(myRand()*DIMENSIONS+0.5);
	int direction = ADI_direction; ADI_direction = (ADI_direction+1)%DIMENSIONS;
	
	
	for( time = 0; time+timeStep <= timeDifference; time += timeStep){

#if DIFFUSION_METHODE == EXPLICIT

		//fprintf( stderr, "EXPLICIT\n");
		// Finite Differences: explicit method	
		for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
			
			voronoiDiagram->voronoiCells[i]->doxygen  = timeStep * (Oxygen_Diffusion  * secondDerivationOxygen(  voronoiDiagram->voronoiCells[i], voronoiDiagram) // diffusion
	#if DIFFUSION_REACTION_SEPARAT == FALSE
			                                          - voronoiDiagram->voronoiCells[i]->oxygen  * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[i]) // reaction
	#endif
			                                          );
			voronoiDiagram->voronoiCells[i]->dglucose = timeStep * (Glucose_Diffusion * secondDerivationGlucose( voronoiDiagram->voronoiCells[i], voronoiDiagram) // diffusion
	#if DIFFUSION_REACTION_SEPARAT == FALSE
			                                          - voronoiDiagram->voronoiCells[i]->glucose * GiveMeTheGlucoseRate( voronoiDiagram->voronoiCells[i]) // reaction
	#endif
			                                          );			                                          
		}

#else


	for( int d=0; d<DIMENSIONS; d++){
	//fprintf(stderr, "dir=%i\n", direction);
		// Finite Differences: alternating direction implicit (ADI)
		double weight = 1.;
		
		//int i = 0;
		//int ii = 0;
		//int iii = 0;
		
		for( int iset=0; iset<voronoiDiagram->xN[(d==0?1:0)]/*20*/; iset++){
			for( int iiset=0; iiset<voronoiDiagram->xN[(d<=1?2:1)]/*20*/; iiset++){

				// base index
				int index=0;
				int index_count = 0;
				for( int dd=0; dd<DIMENSIONS; dd++){
					if( dd!=direction){
						if(index_count==0)
							index += (int)( iset*pow( voronoiDiagram->xN[(d==0?1:0)]/*20*/, dd));
						else
							index += (int)( iiset*pow( voronoiDiagram->xN[(d<=1?2:1)]/*20*/, dd));
						index_count++;
					}
				}

				//fprintf(stderr, "setupMatrixADI(o): %i <= (%i,%i)\n", index, iset, iiset);

	#if   DIFFUSION_METHODE == CRANK_NICHOLSON
				setupMatrixCrankNicholson( voronoiDiagram, index, timeStep, 'o', direction);
				//fprintf( stderr, "CRANK_NICHOLSON\n");
	#elif DIFFUSION_METHODE == ADI
				setupMatrixADI( voronoiDiagram, index, timeStep, 'o', direction);
	#elif DIFFUSION_METHODE == IMPLICIT
				//fprintf( stderr, "ADI\n");
				setupMatrix( voronoiDiagram, index, timeStep, 'o', direction);
	#endif
				//setupMatrixADI( voronoiDiagram, index, timeStep, 'o', direction);
				//fprintf(stderr, "solve\n");
				//solveLinearSystemTridiagonalMatrix( A, b, x, voronoiDiagram->xN[d]);
				solveLinearSystem( A, b, x, voronoiDiagram->xN[d]/*20*/);
				for( int i=0; i<voronoiDiagram->xN[d]/*20*/; i++)
	#if DIFFUSION_METHODE == IMPLICIT
					if( d!=0)
						voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->doxygen += weight*(b[i] - voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->oxygen);
					else
	#endif
						voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->doxygen = weight*(b[i] - voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->oxygen);

				//fprintf(stderr, "setupMatrixADI(g): %i <= (%i,%i)\n", index, iset, iiset);

	#if   DIFFUSION_METHODE == CRANK_NICHOLSON
				setupMatrixCrankNicholson( voronoiDiagram, index, timeStep, 'g', direction);
				//fprintf( stderr, "CRANK_NICHOLSON\n");
	#elif DIFFUSION_METHODE == ADI
				setupMatrixADI( voronoiDiagram, index, timeStep, 'g', direction);
	#elif DIFFUSION_METHODE == IMPLICIT
				//fprintf( stderr, "Implicit\n");
				setupMatrix( voronoiDiagram, index, timeStep, 'g', direction);
	#endif
				//setupMatrixADI( voronoiDiagram, index, timeStep, 'g', direction);
				solveLinearSystem( A, b, x, voronoiDiagram->xN[d]/*20*/);
				//solveLinearSystemTridiagonalMatrix( A, b, x, voronoiDiagram->xN[d]);
				for( int i=0; i<voronoiDiagram->xN[d]/*20*/; i++)
	#if DIFFUSION_METHODE == IMPLICIT
					if( d!=0)
						//add differences
						voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->dglucose += weight*(b[i] - voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->glucose);
					else
	#endif
						// set differences
						voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->dglucose = weight*(b[i] - voronoiDiagram->voronoiCells[(int)(index+i*pow( voronoiDiagram->xN[d]/*20*/, direction))]->glucose);
			}
		}
#endif
		

#if DIFFUSION_METHODE == IMPLICIT
		if(d+1==DIMENSIONS)
#endif
		// update half time step
		{//fprintf(stderr, "update half time step\n");
		for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
			if( voronoiDiagram->voronoiCells[i]->position[0]>voronoiDiagram->xMin[0]+1
			 && voronoiDiagram->voronoiCells[i]->position[0]<voronoiDiagram->xMax[0]-1 
			 && voronoiDiagram->voronoiCells[i]->position[1]>voronoiDiagram->xMin[1]+1
			 && voronoiDiagram->voronoiCells[i]->position[1]<voronoiDiagram->xMax[1]-1
		#if DIMENSIONS > 2
			 && voronoiDiagram->voronoiCells[i]->position[2]>voronoiDiagram->xMin[2]+1
			 && voronoiDiagram->voronoiCells[i]->position[2]<voronoiDiagram->xMax[2]-1 
		#endif
			){	
				voronoiDiagram->voronoiCells[i]->oxygen  += voronoiDiagram->voronoiCells[i]->doxygen;
				if( voronoiDiagram->voronoiCells[i]->oxygen < 0.) voronoiDiagram->voronoiCells[i]->oxygen = 0.;
				voronoiDiagram->voronoiCells[i]->glucose += voronoiDiagram->voronoiCells[i]->dglucose;
				if( voronoiDiagram->voronoiCells[i]->glucose < 0.) voronoiDiagram->voronoiCells[i]->glucose = 0.;
			}
		}}

		// change direction
		direction = (direction+1)%DIMENSIONS;
	}

#if DIFFUSION_REACTION_SEPARAT == TRUE
	// reaction
	for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
		if( voronoiDiagram->voronoiCells[i]->position[0]>voronoiDiagram->xMin[0]+1
		 && voronoiDiagram->voronoiCells[i]->position[0]<voronoiDiagram->xMax[0]-1 
		 && voronoiDiagram->voronoiCells[i]->position[1]>voronoiDiagram->xMin[1]+1
		 && voronoiDiagram->voronoiCells[i]->position[1]<voronoiDiagram->xMax[1]-1
	#if DIMENSIONS > 2
		 && voronoiDiagram->voronoiCells[i]->position[2]>voronoiDiagram->xMin[2]+1
		 && voronoiDiagram->voronoiCells[i]->position[2]<voronoiDiagram->xMax[2]-1 
	#endif
		){	
			voronoiDiagram->voronoiCells[i]->oxygen  += - timeStep * voronoiDiagram->voronoiCells[i]->oxygen  * GiveMeTheOxygenRate(  voronoiDiagram->voronoiCells[i]);
			if( voronoiDiagram->voronoiCells[i]->oxygen < 0.) voronoiDiagram->voronoiCells[i]->oxygen = 0.;
			voronoiDiagram->voronoiCells[i]->glucose += - timeStep * voronoiDiagram->voronoiCells[i]->glucose * GiveMeTheGlucoseRate( voronoiDiagram->voronoiCells[i]);
			if( voronoiDiagram->voronoiCells[i]->glucose < 0.) voronoiDiagram->voronoiCells[i]->glucose = 0.;
		}
	}
#endif
#if DIFFUSION_METHODE != EXPLICIT
	}
#endif
	return timeDifference - time;
}

//#endif

//////////////////////////////////////////////////////////////////////////////////////

void initMatrixImplicit( VoronoiDiagram *voronoiDiagram)
{
	int matrixDim = 1;
	
	for( int d=0; d<DIMENSIONS; d++)
		matrixDim *= voronoiDiagram->xN[d];

	b = (float*) malloc( matrixDim * sizeof(float));

	x = (float*) malloc( matrixDim * sizeof(float));
		
	A = newMatrix( matrixDim, matrixDim);
}


void setupMatrixImplicit( VoronoiDiagram *voronoiDiagram, double timeStep, char molecule)
{
	int di   = 1;
	int dii  = voronoiDiagram->xN[0];
	int diii = voronoiDiagram->xN[0]*voronoiDiagram->xN[1];
	
	double r = Glucose_Diffusion * timeStep/(DX*DX);
	
	int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	
	for( int iii=0; iii<voronoiDiagram->xN[2]; iii++)
	for( int ii=0; ii<voronoiDiagram->xN[1]; ii++)
	for( int i=0; i<voronoiDiagram->xN[0]; i++)
	{
		// actual element
		int m = i*di + ii*dii + iii*diii;
		//fprintf( stderr, "%i ", m);

		// init matrix row
		for( int n=0; n<N; n++)
			A[m][n] = 0.;			
		
		//matrix
		A[m][m] = 1;
		
		// x
		/*if( i>0){
			A[m][m] = r;
			A[m][m-di] = -r; 
		}
		if( i<voronoiDiagram->xN[0]-1){
			A[m][m] = r;
			A[m][m+di] = -r; 
		}
		
		// y
		if( ii>0){
			A[m][m] = r;
			A[m][m-dii] = -r; 
		}
		if( ii<voronoiDiagram->xN[1]-1){
			A[m][m] = r;
			A[m][m+dii] = -r; 
		}
		
		// z
		if( iii>0){
			A[m][m] = r;
			A[m][m-diii] = -r; 
		}
		if( iii<voronoiDiagram->xN[2]-1){
			A[m][m] = r;
			A[m][m+diii] = -r; 
		}*/
		if( i>0 && i<voronoiDiagram->xN[0]-1 && ii>0 && ii<voronoiDiagram->xN[1]-1 && iii>0 && iii<voronoiDiagram->xN[2]-1)
		{
			A[m][m] += 6*r;			
			A[m][m-di] = -r; 
			A[m][m+di] = -r; 
			A[m][m-dii] = -r; 
			A[m][m+dii] = -r; 
			A[m][m-diii] = -r; 
			A[m][m+diii] = -r; 
		}
		
		// vector
		b[m] = voronoiDiagram->voronoiCells[m]->glucose * ( 1. - timeStep * GiveMeTheGlucoseRate(  voronoiDiagram->voronoiCells[m]));
		
	}
}


double UpdateSystemImplicit( VoronoiDiagram *voronoiDiagram, double timeStep, double timeDifference)
{

	double time;
	
	int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	//float **B = newMatrix( N, N);
	
	for( time = 0; time+timeStep <= timeDifference; time += timeStep){
		//fprintf( stderr, "%i. iteration:\nSetup Matrix\n", (int)(time / timeStep + 0.5));
		int passedTime = clock();
		fprintf( stderr, "Set Matrix... \n");
		setupMatrixImplicit( voronoiDiagram, timeStep, 'G');
		fprintf( stderr, "...finished ( %li clocks, %.3lf sec)\n", (clock() - passedTime), (float)(clock() - passedTime)/CLOCKS_PER_SEC);

		/*for(int m=0; m<N; m++){
			for(int n=0; n<N; n++)
				fprintf( stderr, "%6.0lf ", A[m][n]);
			fprintf( stderr, "\n");
		}*/
		/*for(int m=0; m<N; m++){
			for(int n=0; n<N; n++)
				fprintf( stdout, "%i %i %lf \n", m, n, (A[m][n]!=0.?1.:0.));
			//fprintf( stderr, "\n");
		}*/
		
		/*fprintf( stderr, "Solve Matrix\n");
		solveLinearSystemB( A, b, x, N, B);
		
		fprintf( stderr, "Actualize Values\n");
		for(int m=0; m<N; m++)
			voronoiDiagram->voronoiCells[m]->glucose = b[m];
		*/	
		for(int m=0; m<N; m++)
			x[m] = voronoiDiagram->voronoiCells[m]->glucose;

		passedTime = clock();
		fprintf( stderr, "ConjugateGradient... \n");
		ConjugateGradient( A, b, x, N, 1);
		fprintf( stderr, "...finished ( %li clocks, %.3lf sec)\n", (clock() - passedTime), (float)(clock() - passedTime)/CLOCKS_PER_SEC);
		//fprintf( stderr, "...finished ( %lisec)\n", (clock() - passedTime)/CLOCKS_PER_SEC);

		for(int m=0; m<N; m++)
			voronoiDiagram->voronoiCells[m]->glucose = x[m];
	}

	return timeDifference - time;
	
}

void ConjugateGradient_old( float **A, float *b, float *x, int N, int iterations)
{
	float r[N];
	float r2[N];
	
	float p[N];
	float q[N];

	//float beta[N];
	//float alpha[N];
	float beta;
	float alpha;
	
	float temp[N];
	
	
	// initialization of residual vector: r_0
	matrixVectorProduct( A, x, temp, N);
	vectorDifference( b, temp, r, N);
	
	// initialization of p_0
	vectorCopy( r, p, N);


	
	for( int i=0; i<iterations; i++){
		// q_k
		matrixVectorProduct( A, p, q, N);
		
		alpha = dotProduct( r, r, N);
		
		if(i>0){
			beta = alpha / dotProduct( r2, r2, N);
			
			vectorScale( p, beta, p, N);
			vectorSum( r, p, p, N);
		}
		
		alpha /= dotProduct( p, q, N);

		// next p
		vectorScale( p, alpha, temp, N);
		vectorSum( p, temp, p, N);
		
		// next x
		vectorScale( x, alpha, temp, N);
		vectorSum( x, temp, x, N);
	}
}

void ConjugateGradient3( float **A, float *b, float *x, int N, int iterations)
{
	float r[N];
	float r2[N];
	
	float p[N];
	float q[N];

	//float beta[N];
	//float alpha[N];
	float beta;
	float alpha;
	
	float temp[N];
	
	
	// initialization of residual vector: r_0
	matrixVectorProduct( A, x, temp, N);
	vectorDifference( b, temp, r, N);
	
	// initialization of p_0
	vectorCopy( r, p, N);


	
	for( int i=0; i<iterations; i++){
		// q_k
		matrixVectorProduct( A, p, q, N);
		
		alpha = dotProduct( r, r, N);
		
		if(i>0){
			beta = alpha / dotProduct( r2, r2, N);
			
			vectorScale( p, beta, p, N);
			vectorSum( r, p, p, N);
		}
		
		alpha /= dotProduct( p, q, N);

		// next x
		vectorScale( p, alpha, temp, N);
		vectorSum( x, temp, x, N);
		
		// old r
		vectorCopy( r, r2, N);
		
		// next r
		vectorScale( q, -alpha, temp, N);
		vectorSum( r, temp, r, N);
		
		//float *pointer ;
	}
}
void ConjugateGradient( float **A, float *b, float *x, int N, int iterations)
{
	// vectors
	float r[N];
	float w[N];
	float z[N];
	
	// scalars
	float alpha;
	float beta;
	
	float temp[N];
	
	
	// initialization of residual vector: r
	matrixVectorProduct( A, x, temp, N);
	vectorDifference( b, temp, r, N);
	
	// w
	vectorScale( r, -1, w, N);
	
	// z
	matrixVectorProduct( A, w, z, N);
	
	// alpha
	alpha = dotProduct( r, w, N) / dotProduct( w, z, N);
	
	// beta
	beta = 0.;
	
	// x
	vectorScale( w, alpha, temp, N);
	vectorSum( x, temp, x, N);
	
	for( int i=0; i<iterations; i++){
		vectorScale( z, alpha, temp, N);
		vectorDifference( r, temp, r, N);
		
		if( sqrt( dotProduct( r, r, N)) < 1e-10)
			return;
		
		// B = (r'*z)/(w'*z);
		beta = dotProduct( r, z, N) / dotProduct( w, z, N);
		
		// w = -r + B*w;
		vectorScale( w, beta, w, N);
		vectorDifference( w, r, w, N);
		
		// z = A*w;
		matrixVectorProduct( A, w, z, N);
		
		// a = (r'*w)/(w'*z);
		alpha = dotProduct( r, w, N) / dotProduct( w, z, N);
		
		// x = x + a*w;
		vectorScale( w, alpha, temp, N);
		vectorSum( x, temp, x, N);
	}
}
