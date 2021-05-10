/*
 * Interpolation.cpp
 *
 *  Created on: 22.02.2010
 *      Author: jagiella
 */

#include <math.h>
#include <stdio.h>

#include "VoronoiDiagramExtended.h"
#include "Interpolation.h"


#define COLORS 3
//#define DIM_MAX		100

double getHeight( double *p0, double *p1, double *p2, double *p3){
	double vector1[3], vector2[3], normal[3];

	vector1[0] = p2[0] - p1[0];
	vector1[1] = p2[1] - p1[1];
	vector1[2] = p2[2] - p1[2];

	vector2[0] = p3[0] - p1[0];
	vector2[1] = p3[1] - p1[1];
	vector2[2] = p3[2] - p1[2];

	// normal vector of plane p1-p2-p3
	normal[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
	normal[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
	normal[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];

	// normalize normal vector
	double length = sqrt( pow(normal[0],2.) + pow(normal[1],2.) + pow(normal[2],2.));
	if( length==0.){
		return 0.;
	}
	normal[0] /= length;
	normal[1] /= length;
	normal[2] /= length;

	// distance: D = n*(v-pi), i=1,2 or 3
	return fabs( normal[0]*(p0[0]-p1[0]) + normal[1]*(p0[1]-p1[1]) + normal[2]*(p0[2]-p1[2]));
}

double getDistance( double *p1, double *p2){
	return sqrt( pow(p1[0]-p2[0], 2.) + pow(p1[1]-p2[1], 2.) + pow(p1[2]-p2[2], 2.));
}

double getSurface( double *p1, double *p2, double *p3){
	double a = getDistance( p1, p2);
	double b = getDistance( p2, p3);
	double c = getDistance( p3, p1);
	double s = (a+b+c)/2.;
	return sqrt( s*(s-a)*(s-b)*(s-c));
}

double getVolume( double *p0, double *p1, double *p2, double *p3)
{
	return 1./3. * getSurface( p1, p2, p3) * getHeight( p0, p1, p2, p3);
}


float linearInterpolation( float x, double A, double B, float cA, float cB)
{
	//  A-----B

	return cA + (x-A)*(cB-cA)/(B-A);
}

float bilinearInterpolation( float x, float y, double *A, double *D, float cA, float cB, float cC, float cD)
{
	//  C-----D
	//  |     |
	//  |     |
	//  A-----B

	float xRange = D[0]-A[0];
	float yRange = D[1]-A[1];

	return (cA*(D[0]-x)*(D[1]-y) +
			cB*(x-A[0])*(D[1]-y) +
			cC*(D[0]-x)*(y-A[1]) +
			cD*(x-A[0])*(y-A[1]))/xRange/yRange;
}

float trilinearInterpolation( float x, float y, float z, double *A, double *H,
		float cA, float cB, float cC, float cD, float cE, float cF, float cG, float cH)
{
	//     G-----H
	//  C-----D  |
	//  |  |  |  |
	//  |  E  |  F
	//  A-----B

	double D[3] = {H[0], H[1], A[2]};
	double E[3] = {A[0], A[1], H[2]};

	return linearInterpolation(
			z,
			A[2],
			H[2],
			bilinearInterpolation( x, y, A, D, cA, cB, cC, cD),
			bilinearInterpolation( x, y, E, H, cE, cF, cG, cH));
}


float getTetrahedronVolumeAboveThreshold( double *A, double *B, double *C, double *D, float cA, float cB, float cC, float cD, float threshold)
{
	// sorting: big < small
	if( cA < cB){
		float  c = cA; cA = cB; cB = c;
		double *t = A;  A =  B;  B = t;
	}
	if( cB < cC){
		float  c = cB; cB = cC; cC = c;
		double *t = B;  B =  C;  C = t;
	}
	if( cC < cD){
		float  c = cC; cC = cD; cD = c;
		double *t = C;  C =  D;  D = t;
	}

	if( cA < cB){
		float  c = cA; cA = cB; cB = c;
		double *t = A;  A =  B;  B = t;
	}
	if( cB < cC){
		float  c = cB; cB = cC; cC = c;
		double *t = B;  B =  C;  C = t;
	}

	if( cA < cB){
		float  c = cA; cA = cB; cB = c;
		double *t = A;  A = B;   B = t;
	}
	//fprintf( stderr, "%f > %f > %f > %f\n", cA, cB, cC, cD);

	if( cA < threshold){
		//fprintf( stderr, "all below threshold\n");
		return 0.;
	}
	if( cB < threshold){
		//fprintf( stderr, "threshold between A and (B, C, D)\n");
		double i;

		i = (threshold - cA)/(cB - cA);
		double tB[3] = {A[0] + i*(B[0]-A[0]), A[1] + i*(B[1]-A[1]), A[2] + i*(B[2]-A[2])};

		i = (threshold - cA)/(cC - cA);
		double tC[3] = {A[0] + i*(C[0]-A[0]), A[1] + i*(C[1]-A[1]), A[2] + i*(C[2]-A[2])};

		i = (threshold - cA)/(cD - cA);
		double tD[3] = {A[0] + i*(D[0]-A[0]), A[1] + i*(D[1]-A[1]), A[2] + i*(D[2]-A[2])};

		//fprintf( stderr, "=> %f\n", getVolume( A, tB, tC, tD));
		return getVolume( A, tB, tC, tD);
	}
	if( cC < threshold){
		//fprintf( stderr, "threshold between (A, B) and (C, D)\n"); //return 0.;
		double i;

		i = (threshold - cA)/(cC - cA); // A->C
		double tA1[3] = {A[0] + i*(C[0]-A[0]), A[1] + i*(C[1]-A[1]), A[2] + i*(C[2]-A[2])};
		i = (threshold - cA)/(cD - cA); // A->D
		double tA2[3] = {A[0] + i*(D[0]-A[0]), A[1] + i*(D[1]-A[1]), A[2] + i*(D[2]-A[2])};

		i = (threshold - cB)/(cC - cB); // B->C
		double tB1[3] = {B[0] + i*(C[0]-B[0]), B[1] + i*(C[1]-B[1]), B[2] + i*(C[2]-B[2])};
		i = (threshold - cB)/(cD - cB); // B->D
		double tB2[3] = {B[0] + i*(D[0]-B[0]), B[1] + i*(D[1]-B[1]), B[2] + i*(D[2]-B[2])};

		//fprintf( stderr, "=> %f\n", getVolume( A, tA1, tA2, B) + getVolume( tA1, tA2, tB2, B) + getVolume( tA1, tB1, tB2, B));
		return getVolume( A, tA1, tA2, B) + getVolume( tA1, tA2, tB2, B) + getVolume( tA1, tB1, tB2, B);
	}
	if( cD < threshold){
		//fprintf( stderr, "threshold between (A, B, C) and D\n"); //return 0.;
		double i;

		i = (threshold - cD)/(cA - cD);
		double tA[3] = {D[0] + i*(A[0]-D[0]), D[1] + i*(A[1]-D[1]), D[2] + i*(A[2]-D[2])};
		//fprintf( stderr, "%lf | ", i);
		i = (threshold - cD)/(cB - cD);
		double tB[3] = {D[0] + i*(B[0]-D[0]), D[1] + i*(B[1]-D[1]), D[2] + i*(B[2]-D[2])};
		//fprintf( stderr, "%lf | ", i);
		i = (threshold - cD)/(cC - cD);
		double tC[3] = {D[0] + i*(C[0]-D[0]), D[1] + i*(C[1]-D[1]), D[2] + i*(C[2]-D[2])};
		//fprintf( stderr, "%lf | ", i);

		//fprintf( stderr, "=> %f = %f - %f\n", getVolume( A, B, C, D) - getVolume( A, tB, tC, tD), getVolume( A, B, C, D), getVolume( A, tB, tC, tD));
		return getVolume( A, B, C, D) - getVolume( A, tA, tB, tC);
	}
	//fprintf( stderr, "all above threshold\n");
	//fprintf( stderr, "=> %f\n", getVolume( A, B, C, D));
	return getVolume( A, B, C, D);

	/*if( cA >= threshold && cB >= threshold && cC >= threshold){
		// all above threshold
		fprintf( stderr, "all above threshold\n"); return getArea( A, B, C);
	}

	if( cA < threshold && cB < threshold && cC < threshold){
		// all below threshold
		fprintf( stderr, "all below threshold\n"); return 0.;
	}

	if( cA >= threshold && cB < threshold){
		fprintf( stderr, "threshold between A and B\n");
	}
	if( cA < threshold && cB >= threshold){
		fprintf( stderr, "threshold between B and A\n");
	}
	if( cA >= threshold && cC < threshold){
		fprintf( stderr, "threshold between A and C\n");
	}
	if( cA < threshold && cC >= threshold){
		fprintf( stderr, "threshold between C and A\n");
	}
	if( cB >= threshold && cC < threshold){
		fprintf( stderr, "threshold between B and C\n");
	}
	if( cB < threshold && cC >= threshold){
		fprintf( stderr, "threshold between C and B\n");
	}*/

}

VoronoiDiagram *Interpolation::voronoiDiagram = NULL;

#define VALUE( vc) \
		(vc->oxygen * vc->glucose)
//	(vc->oxygen)
//	(vc->glucose)

float Interpolation::getRectangleAreaAboveThreshold3D( int index, char mode, float threshold)
{
	int countDim = (int) (pow( voronoiDiagram->countVoronoiCells, 1./DIMENSIONS) + 0.5);
	//fprintf( stderr, "countDim %i\n", countDim);
	// spatial index
	int ix = (int) (countDim * voronoiDiagram->voronoiCells[index]->position[0]),
		iy = (int) (countDim * voronoiDiagram->voronoiCells[index]->position[1]),
		iz = (int) (countDim * voronoiDiagram->voronoiCells[index]->position[2]);
	//fprintf( stderr, "point (%lf, %lf)->index %i ==> %i\n", voronoiDiagram->voronoiCells[index]->position[0], voronoiDiagram->voronoiCells[index]->position[1], voronoiDiagram->voronoiCells[index]->index, iy*countDim + ix);

	double A[3] = { 0., 0., 0.},
		   D[2] = { 1., 1.},
		   H[3] = { 1., 1., 1.};

	float cA, cB, cC, cD, cE, cF, cG, cH;
	cA = VALUE( voronoiDiagram->voronoiCells[index]);
	//fprintf( stderr, "cA: %f\n", cA);
	//int count = 0;
	float area = 0.;
	//float steps = 100.;
	//float stepSize = 1./steps;
	for( int dx=-1; dx<=1; dx+=2)
	for( int dy=-1; dy<=1; dy+=2)
	for( int dz=-1; dz<=1; dz+=2)
		if( ix+dx>=0 && ix+dx<countDim && iy+dy>=0 && iy+dy<countDim && iz+dz>=0 && iz+dz<countDim){
			//     G-----H
			//  C-----D  |
			//  |  |  |  |
			//  |  E  |  F
			//  A-----B
			cB = VALUE( voronoiDiagram->voronoiCells[index + dx]);
			cC = VALUE( voronoiDiagram->voronoiCells[index + dy*countDim]);
			cD = VALUE( voronoiDiagram->voronoiCells[index + dx + dy*countDim]);

			cE = VALUE( voronoiDiagram->voronoiCells[index + dz*countDim*countDim]);
			cF = VALUE( voronoiDiagram->voronoiCells[index + dx + dz*countDim*countDim]);
			cG = VALUE( voronoiDiagram->voronoiCells[index + dy*countDim + dz*countDim*countDim]);
			cH = VALUE( voronoiDiagram->voronoiCells[index + dx + dy*countDim + dz*countDim*countDim]);


			// numerical discretization (bilinear interpol.)
//			for( float x=0.5*stepSize; x<=0.5; x+=stepSize)
//			for( float y=0.5*stepSize; y<=0.5; y+=stepSize)
//			for( float z=0.5*stepSize; z<=0.5; z+=stepSize)
//				count += ( trilinearInterpolation( x, y, z, A, H, cA, cB, cC, cD, cE, cF, cG, cH) >= threshold? 1 : 0);

			// exact caculation (triangular interpol.)
			cH = trilinearInterpolation( 0.5, 0.5, 0.5, A, H, cA, cB, cC, cD, cE, cF, cG, cH);

			cD = bilinearInterpolation( 0.5, 0.5, A, D, cA, cB, cC, cD);
			cF = bilinearInterpolation( 0.5, 0.5, A, D, cA, cB, cE, cF);
			cG = bilinearInterpolation( 0.5, 0.5, A, D, cA, cE, cC, cG);

			cB = (VALUE( voronoiDiagram->voronoiCells[index + dx])+cA)/2.;
			cC = (VALUE( voronoiDiagram->voronoiCells[index + dy*countDim])+cA)/2.;
			cE = (VALUE( voronoiDiagram->voronoiCells[index + dz*countDim*countDim])+cA)/2.;

			double tA[3] = { 0., 0., 0.},
				   tB[3] = { .5, 0., 0.},
				   tC[3] = { 0., .5, 0.},
				   tD[3] = { .5, .5, 0.},
				   tE[3] = { 0., 0., .5},
				   tF[3] = { .5, 0., .5},
				   tG[3] = { 0., .5, .5},
				   tH[3] = { .5, .5, .5};

			//  C
			//  | \           .
			//  | ,E\         .
			//  A-----B
			area += getTetrahedronVolumeAboveThreshold( tA, tB, tC, tE, cA, cB, cC, cE, threshold);

			//    _ _ ---H
			//  C-----D''/
			//    \   | /
			//      \ |/
			//        B
			area += getTetrahedronVolumeAboveThreshold( tB, tC, tD, tH, cB, cC, cD, cH, threshold);

			//    _G-----H
			//  C  |   /
			//    \| /
			//     E
			//
			area += getTetrahedronVolumeAboveThreshold( tC, tE, tG, tH, cC, cE, cG, cH, threshold);

			//           H
			//         //|
			//       / / |
			//     E  / _F
			//       \B-
			area += getTetrahedronVolumeAboveThreshold( tB, tE, tF, tH, cB, cE, cF, cH, threshold);

			//    _ _ - -H
			//  C      //
			//   \ \ / /
			//    E, \/
			//      '-B
			area += getTetrahedronVolumeAboveThreshold( tB, tC, tE, tH, cB, cC, cE, cH, threshold);
		}else{
			if(cA>=threshold){
//				count += (int)(steps*steps*steps/8);
				area += pow(0.5, DIMENSIONS);
			}
		}

	//fprintf( stderr, "%lf <=> %lf\n", (float)count/(steps*steps*steps), area);

	return area;
}

float Interpolation::getDiscritizedVolumeAboveThreshold3D( int index, char mode, float threshold)
{
	int countDim = (int) (pow( voronoiDiagram->countVoronoiCells, 1./DIMENSIONS) + 0.5);
	//fprintf( stderr, "countDim %i\n", countDim);
	// spatial index
	int ix = (int) (voronoiDiagram->voronoiCells[index]->position[0]),
		iy = (int) (voronoiDiagram->voronoiCells[index]->position[1]),
		iz = (int) (voronoiDiagram->voronoiCells[index]->position[2]);
	//fprintf( stderr, "point (%lf, %lf)->index %i ==> %i\n", voronoiDiagram->voronoiCells[index]->position[0], voronoiDiagram->voronoiCells[index]->position[1], voronoiDiagram->voronoiCells[index]->index, iy*countDim + ix);
	//fprintf( stderr, "Interpolation (countDim=%i, index=%i, ix:%i, iy:%i, iz:%i)\n", countDim, index, ix, iy, iz);

	double A[3] = { 0., 0., 0.},
	//	   D[2] = { 1., 1.},
		   H[3] = { 1., 1., 1.};

	float cA, cB, cC, cD, cE, cF, cG, cH;
	cA = VALUE( voronoiDiagram->voronoiCells[index]);
	//fprintf( stderr, "cA: %f\n", cA);
	int count = 0;
	//float area = 0.;
	float steps = 10.;
	float stepSize = 1./steps;
	for( int dx=-1; dx<=1; dx+=2)
	for( int dy=-1; dy<=1; dy+=2)
	for( int dz=-1; dz<=1; dz+=2)
		if( ix+dx>=0 && ix+dx<countDim && iy+dy>=0 && iy+dy<countDim && iz+dz>=0 && iz+dz<countDim){
			//     G-----H
			//  C-----D  |
			//  |  |  |  |
			//  |  E  |  F
			//  A-----B
			//fprintf( stderr, " -> (%i, %i, %i, %i, %i, %i, %i)\n",
			//		index + dx, index + dy*countDim, index + dx + dy*countDim, index + dz*countDim*countDim, index + dx + dz*countDim*countDim, index + dy*countDim + dz*countDim*countDim, index + dx + dy*countDim + dz*countDim*countDim);
			cB = VALUE( voronoiDiagram->voronoiCells[index + dx*countDim*countDim]);
			cC = VALUE( voronoiDiagram->voronoiCells[index + dy*countDim]);
			cD = VALUE( voronoiDiagram->voronoiCells[index + dx*countDim*countDim + dy*countDim]);

			cE = VALUE( voronoiDiagram->voronoiCells[index + dz]);
			cF = VALUE( voronoiDiagram->voronoiCells[index + dx*countDim*countDim + dz]);
			cG = VALUE( voronoiDiagram->voronoiCells[index + dy*countDim + dz]);
			cH = VALUE( voronoiDiagram->voronoiCells[index + dx*countDim*countDim + dy*countDim + dz]);


			// numerical discretization (bilinear interpol.)
			for( float x=0.5*stepSize; x<=0.5; x+=stepSize)
			for( float y=0.5*stepSize; y<=0.5; y+=stepSize)
			for( float z=0.5*stepSize; z<=0.5; z+=stepSize)
				count += ( trilinearInterpolation( x, y, z, A, H, cA, cB, cC, cD, cE, cF, cG, cH) >= threshold? 1 : 0);

			// exact caculation (triangular interpol.)
/*			cH = trilinearInterpolation( 0.5, 0.5, 0.5, A, H, cA, cB, cC, cD, cE, cF, cG, cH);

			cD = bilinearInterpolation( 0.5, 0.5, A, D, cA, cB, cC, cD);
			cF = bilinearInterpolation( 0.5, 0.5, A, D, cA, cB, cE, cF);
			cG = bilinearInterpolation( 0.5, 0.5, A, D, cA, cE, cC, cG);

			cB = (VALUE( voronoiDiagram->voronoiCells[index + dx])+cA)/2.;
			cC = (VALUE( voronoiDiagram->voronoiCells[index + dy*countDim])+cA)/2.;
			cE = (VALUE( voronoiDiagram->voronoiCells[index + dz*countDim*countDim])+cA)/2.;

			double tA[3] = { 0., 0., 0.},
				   tB[3] = { .5, 0., 0.},
				   tC[3] = { 0., .5, 0.},
				   tD[3] = { .5, .5, 0.},
				   tE[3] = { 0., 0., .5},
				   tF[3] = { .5, 0., .5},
				   tG[3] = { 0., .5, .5},
				   tH[3] = { .5, .5, .5};

			//  C
			//  | \
			//  | ,E\
			//  A-----B
			area += getTetrahedronVolumeAboveThreshold( tA, tB, tC, tE, cA, cB, cC, cE, threshold);

			//    _ _ ---H
			//  C-----D''/
			//    \   | /
			//      \ |/
			//        B
			area += getTetrahedronVolumeAboveThreshold( tB, tC, tD, tH, cB, cC, cD, cH, threshold);

			//    _G-----H
			//  C  |   /
			//    \| /
			//     E
			//
			area += getTetrahedronVolumeAboveThreshold( tC, tE, tG, tH, cC, cE, cG, cH, threshold);

			//           H
			//         //|
			//       / / |
			//     E  / _F
			//       \B-
			area += getTetrahedronVolumeAboveThreshold( tB, tE, tF, tH, cB, cE, cF, cH, threshold);

			//    _ _ - -H
			//  C      //
			//   \ \ / /
			//    E, \/
			//      '-B
			area += getTetrahedronVolumeAboveThreshold( tB, tC, tE, tH, cB, cC, cE, cH, threshold);
*/
		}else{
			if(cA>=threshold){
				count += (int)(steps*steps*steps/8);
//				area += pow(0.5, DIMENSIONS);
			}
		}

	//fprintf( stderr, "%lf <=> %lf\n", (float)count/(steps*steps*steps), area);

	//fprintf( stderr, "...finished Interpolation\n");

	return (float)count/(steps*steps*steps);
}


float Interpolation::getDiscritizedVolumeAboveThreshold2D( int index, char mode, float threshold)
{
	int countDim = (int) (pow( voronoiDiagram->countVoronoiCells, 1./DIMENSIONS) + 0.5);
	//fprintf( stderr, "countDim %i\n", countDim);
	// spatial index
	int ix = (int) (voronoiDiagram->voronoiCells[index]->position[0]),
		iy = (int) (voronoiDiagram->voronoiCells[index]->position[1]);
	//fprintf( stderr, "point (%lf, %lf)->index %i ==> %i\n", voronoiDiagram->voronoiCells[index]->position[0], voronoiDiagram->voronoiCells[index]->position[1], voronoiDiagram->voronoiCells[index]->index, iy*countDim + ix);
	//fprintf( stderr, "Interpolation (countDim=%i, index=%i, ix:%i, iy:%i)\n", countDim, index, ix, iy);

	double A[2] = { 0., 0.},
		   D[2] = { 1., 1.};

	float cA, cB, cC, cD;
	cA = VALUE( voronoiDiagram->voronoiCells[index]);
	//fprintf( stderr, "cA: %f\n", cA);
	int count = 0;
	//float area = 0.;
	float steps = 10.;
	float stepSize = 1./steps;
	for( int dx=-1; dx<=1; dx+=2)
	for( int dy=-1; dy<=1; dy+=2)
		if( ix+dx>=0 && ix+dx<countDim && iy+dy>=0 && iy+dy<countDim){
			//  C-----D
			//  |     |
			//  |     |
			//  A-----B
			//fprintf( stderr, " -> (%i, %i, %i, %i) -> ix:%i, iy:%i -> dx:%i, dy:%i\n",
			//		index, index + dx*countDim, index + dy, index + dx*countDim + dy, ix, iy, dx, dy);
			cB = VALUE( voronoiDiagram->voronoiCells[index + dx]);
			cC = VALUE( voronoiDiagram->voronoiCells[index + dy*countDim]);
			cD = VALUE( voronoiDiagram->voronoiCells[index + dx + dy*countDim]);


			// numerical discretization (bilinear interpol.)
			for( float x=0.5*stepSize; x<=0.5; x+=stepSize)
			for( float y=0.5*stepSize; y<=0.5; y+=stepSize)
				count += ( bilinearInterpolation( x, y, A, D, cA, cB, cC, cD) >= threshold? 1 : 0);

			// exact caculation (triangular interpol.)
/*			cH = trilinearInterpolation( 0.5, 0.5, 0.5, A, H, cA, cB, cC, cD, cE, cF, cG, cH);

			cD = bilinearInterpolation( 0.5, 0.5, A, D, cA, cB, cC, cD);
			cF = bilinearInterpolation( 0.5, 0.5, A, D, cA, cB, cE, cF);
			cG = bilinearInterpolation( 0.5, 0.5, A, D, cA, cE, cC, cG);

			cB = (VALUE( voronoiDiagram->voronoiCells[index + dx])+cA)/2.;
			cC = (VALUE( voronoiDiagram->voronoiCells[index + dy*countDim])+cA)/2.;
			cE = (VALUE( voronoiDiagram->voronoiCells[index + dz*countDim*countDim])+cA)/2.;

			double tA[3] = { 0., 0., 0.},
				   tB[3] = { .5, 0., 0.},
				   tC[3] = { 0., .5, 0.},
				   tD[3] = { .5, .5, 0.},
				   tE[3] = { 0., 0., .5},
				   tF[3] = { .5, 0., .5},
				   tG[3] = { 0., .5, .5},
				   tH[3] = { .5, .5, .5};

			//  C
			//  | \
			//  | ,E\
			//  A-----B
			area += getTetrahedronVolumeAboveThreshold( tA, tB, tC, tE, cA, cB, cC, cE, threshold);

			//    _ _ ---H
			//  C-----D''/
			//    \   | /
			//      \ |/
			//        B
			area += getTetrahedronVolumeAboveThreshold( tB, tC, tD, tH, cB, cC, cD, cH, threshold);

			//    _G-----H
			//  C  |   /
			//    \| /
			//     E
			//
			area += getTetrahedronVolumeAboveThreshold( tC, tE, tG, tH, cC, cE, cG, cH, threshold);

			//           H
			//         //|
			//       / / |
			//     E  / _F
			//       \B-
			area += getTetrahedronVolumeAboveThreshold( tB, tE, tF, tH, cB, cE, cF, cH, threshold);

			//    _ _ - -H
			//  C      //
			//   \ \ / /
			//    E, \/
			//      '-B
			area += getTetrahedronVolumeAboveThreshold( tB, tC, tE, tH, cB, cC, cE, cH, threshold);
*/
		}else{
			if(cA>=threshold){
				count += (int)(steps*steps/4);
//				area += pow(0.5, DIMENSIONS);
			}
		}

	//fprintf( stderr, "%lf <=> %lf\n", (float)count/(steps*steps*steps), area);

	//fprintf( stderr, "...finished Interpolation\n");

	return (float)count/(steps*steps);
}


float Interpolation::getDiscritizedVolumeAboveThreshold1D( int index, char mode, float threshold)
{
	int countDim = (int) (pow( voronoiDiagram->countVoronoiCells, 1./DIMENSIONS) + 0.5);
	//fprintf( stderr, "countDim %i\n", countDim);
	// spatial index
	int ix = (int) (voronoiDiagram->voronoiCells[index]->position[0]);
	//fprintf( stderr, "point (%lf, %lf)->index %i ==> %i\n", voronoiDiagram->voronoiCells[index]->position[0], voronoiDiagram->voronoiCells[index]->position[1], voronoiDiagram->voronoiCells[index]->index, iy*countDim + ix);
	//fprintf( stderr, "Interpolation (countDim=%i, index=%i, ix:%i, iy:%i, iz:%i)\n", countDim, index, ix, iy, iz);

	//fprintf( stderr, "ix = %i\n", ix);

	double A = 0.,
		   B = 1.;

	float cA, cB;
	cA = VALUE( voronoiDiagram->voronoiCells[index]);
	//fprintf( stderr, "cA: %f\n", cA);
	int count = 0;
	//float area = 0.;
	float steps = 100.;
	float stepSize = 1./steps;
	for( int dx=-1; dx<=1; dx+=2)
		if( ix+dx>=0 && ix+dx<countDim){
			//fprintf( stderr, "test\n");
			//  A-----B
			//fprintf( stderr, " -> (%i, %i, %i, %i, %i, %i, %i)\n",
			//		index + dx, index + dy*countDim, index + dx + dy*countDim, index + dz*countDim*countDim, index + dx + dz*countDim*countDim, index + dy*countDim + dz*countDim*countDim, index + dx + dy*countDim + dz*countDim*countDim);
			cB = VALUE( voronoiDiagram->voronoiCells[index + dx]);
			//fprintf( stderr, "Interval (%i, %i) -> (%lf, %lf)\n", ix, ix+dx, cA, cB);

			// numerical discretization (bilinear interpol.)
			//fprintf( stderr, "intervall: [%lf, %lf]\n", 0.5*stepSize, 0.5);
			for( float x=0.5*stepSize; x<=0.5; x+=stepSize)
				count += ( linearInterpolation( x, A, B, cA, cB) >= threshold? 1 : 0);
		}else{
			if(cA>=threshold){
				count += (int)(steps/2);
//				area += pow(0.5, DIMENSIONS);
			}
		}

	//fprintf( stderr, "%lf <=> %lf\n", (float)count/(steps*steps*steps), area);

	//fprintf( stderr, "...finished Interpolation\n");

	return (float)count/steps;
}

#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

float Interpolation::getDiscritizedVolumeAboveThreshold1DCubicSpline( int index, char mode, float threshold)
{
	int countDim = (int) (pow( voronoiDiagram->countVoronoiCells, 1./DIMENSIONS) + 0.5);
	//fprintf( stderr, "countDim %i\n", countDim);
	// spatial index
	int ix = (int) (voronoiDiagram->voronoiCells[index]->position[0]);
	//fprintf( stderr, "point (%lf, %lf)->index %i ==> %i\n", voronoiDiagram->voronoiCells[index]->position[0], voronoiDiagram->voronoiCells[index]->position[1], voronoiDiagram->voronoiCells[index]->index, iy*countDim + ix);
	//fprintf( stderr, "Interpolation (countDim=%i, index=%i, ix:%i, iy:%i, iz:%i)\n", countDim, index, ix, iy, iz);

	//fprintf( stderr, "ix = %i\n", ix);

	// set range
	int intervalSize = 3;
	int min = MAX( ix - intervalSize, 0);
	int max = MIN( ix + intervalSize, countDim-1);
	int N = (max-min)+1; // number of points

	double X[N], Y[N], a[N], b[N], c[N], d[N];
	//for( int i=0; i<N; i++)
	for( int i=0; i<N; i++)
	{
		X[i] = i+min+0.5;
		Y[i] = VALUE( voronoiDiagram->voronoiCells[i+min]);
		//fprintf( stderr, "X=%lf, Y=%lf\n", X[i], Y[i]);
	}

	//fprintf( stderr, "ix=%i, min=%i, max=%i, N=%i\n", ix, min, max, N);

	cubicSplines1Dcubic( N, X, Y, a, b, c, d);
	//cubicSplines1Dconstraint( N, X, Y, a, b, c, d);

	int count = 0;
	int countSteps = 0;
	float steps = 100.;
	float stepSize = 1./steps;
	// for both splines
	/*for( int i=ix-1; i<=ix; i++)
		if(i>=0 && i+1<countDim){
			for( float x=(i+0.5) + 0.5*stepSize; x<=0.5; x+=stepSize)
				count += ( a[i]*pow(x - X[i],3.) + b[i]*pow(x - X[i],2.) + c[i]*(x - X[i]) + d[i] >= threshold ? 1 : 0);
		}
*/
	for( float x=ix + 0.5*stepSize; x< ix+1.; x+=stepSize){
		int i = (x<ix+0.5 ? ix-1-min: ix-min);
		//count += ( a[i]*pow(x - X[i],3.) + b[i]*pow(x - X[i],2.) + c[i]*(x - X[i]) + d[i] >= threshold ? 1 : 0);
		//fprintf( stdout, "%lf %lf\n", x, a[i]*pow(x - X[i],3.) + b[i]*pow(x - X[i],2.) + c[i]*(x - X[i]) + d[i]);
		count += ( d[i]*pow(x,3.) + c[i]*pow(x,2.) + b[i]*(x) + a[i] >= threshold ? 1 : 0);
		//fprintf( stdout, "%lf %lf\n", x, d[i]*pow(x,3.) + c[i]*pow(x,2.) + b[i]*(x) + a[i]);
		countSteps++;
	}

	/*for( int dx=-1; dx<=1; dx+=2)
		if( ix+dx>=0 && ix+dx<countDim){
			//fprintf( stderr, "test\n");
			//  A-----B
			//fprintf( stderr, " -> (%i, %i, %i, %i, %i, %i, %i)\n",
			//		index + dx, index + dy*countDim, index + dx + dy*countDim, index + dz*countDim*countDim, index + dx + dz*countDim*countDim, index + dy*countDim + dz*countDim*countDim, index + dx + dy*countDim + dz*countDim*countDim);
			cB = VALUE( voronoiDiagram->voronoiCells[index + dx]);
			//fprintf( stderr, "Interval (%i, %i) -> (%lf, %lf)\n", ix, ix+dx, cA, cB);

			// numerical discretization (bilinear interpol.)
			//fprintf( stderr, "intervall: [%lf, %lf]\n", 0.5*stepSize, 0.5);
			for( float x=0.5*stepSize; x<=0.5; x+=stepSize)
				count += ( linearInterpolation( x, A, B, cA, cB) >= threshold? 1 : 0);
		}else{
			if(cA>=threshold){
				count += (int)(steps/2);
//				area += pow(0.5, DIMENSIONS);
			}
		}*/

	//fprintf( stderr, "%lf <=> %lf\n", (float)count/(steps*steps*steps), area);

	//fprintf( stderr, "...finished Interpolation\n");

	fprintf( stderr, "count=%i, countSteps=%i, steps=%lf\n", count, countSteps, steps);

	return (float)count/steps;
}

void TridiagonalSolve (const double *a, const double *b, double *c, double *d, double *x, int n){

	/* Modify the coefficients. */
	c[0] /= b[0];	/* Division by zero risk. */
	d[0] /= b[0];	/* Division by zero would imply a singular matrix. */
	for (int i = 1; i < n; i++){
		double id = 1 / (b[i] - c[i-1] * a[i]);  /* Division by zero risk. */
		c[i] *= id;	                         /* Last value calculated is redundant. */
		d[i] = (d[i] - d[i-1] * a[i]) * id;
	}

	/* Now back substitute. */
	x[n - 1] = d[n - 1];
	for (int i = n - 2; i >= 0; i--)
		x[i] = d[i] - c[i] * x[i + 1];
}


void cubicSplines1Dcubic( int N, double *X, double *Y,
		double *a, double *b, double *c, double *d)
{
	// cubic spline interpolation: implicit
	//double a[N];
	//double b[N];
	//double c[N];
	//double d[N];
	double x[N];
	int M = N-1;

	// solve lin. system to get b[i]
	b[0] = 1; c[0] = 0; d[0] = 0;
	for( int i=1; i<M; i++){
		a[i] = X[i] - X[i-1];
		b[i] = 2*(X[i+1] - X[i-1]);
		c[i] = X[i+1] - X[i];
		// wikipedia
		//d[i] = 6.*(Y[i+1] - Y[i])/(X[i+1] - X[i]) - 6.*(Y[i] - Y[i-1])/(X[i] - X[i-1]);
		// nick
		d[i] = 3.*(Y[i+1] - Y[i])/(X[i+1] - X[i]) - 3.*(Y[i] - Y[i-1])/(X[i] - X[i-1]);
	}
	a[M] = 0; b[M] = 1; d[M] = 0;
	TridiagonalSolve (a, b, c, d, x, N);

	// calculate all a[i], c[i], d[i]
	for( int i=0; i<M+1; i++)
		// b[i]
		b[i] = x[i];

	for( int i=0; i<M+1; i++){
		// d[i]
		d[i] = Y[i];

		// c[i]
		c[i] = (Y[i+1] - Y[i])/(X[i+1] - X[i]) - (b[i+1]+2.*b[i])/3.*(X[i+1] - X[i]);

		// a[i]
		a[i] = (b[i+1]-b[i])/3./(X[i+1] - X[i]);
	}

	// output splines
	FILE *fp = fopen( "s3.dat", "w+");
	for(int i=0; i<M; i++)
		for( double x = X[i]; x <= X[i+1]; x+=(X[i+1]-X[i])/10.){

			//printf( "%10.3lf %10.3lf\n", X[i], Y[i]);
			//fprintf( fp, "%10.3lf %10.3lf\n", x, a[i]*pow(x - X[i],3.) + b[i]*pow(x - X[i],2.) + c[i]*(x - X[i]) + d[i]);
			//fprintf( fp, "%10.3lf %10.3lf\n", x, a[i]*pow(x - X[i],3.) + b[i]*pow(x - X[i],2.) + c[i]*(x - X[i]) + d[i]);
			/*fprintf( fp, "%10.3lf %10.3lf\n", X,
				(x[i+1]*pow(X-X[i],3)+x[i]*pow(X[i+1]-X,3)) / 6. / (X[i+1]-X[i])
				+(Y[i+1] / (X[i+1]-X[i]) - (X[i+1]-X[i]) / 6. * x[i+1]) * (X-f[i  ][0])
				-(f[i  ][1] / (X[i+1]-X[i]) - (X[i+1]-X[i]) / 6. * x[i  ]) * (X-X[i+1])
			);*/
	}
	fclose( fp);

}

#define SIGN(x) (x<0 ? -1 : 1)
void cubicSplines1Dconstraint( int N, double *X, double *Y,
		double *a, double *b, double *c, double *d)
{
	// cubic spline interpolation: implicit
	/*double a[N];
	double b[N];
	double c[N];
	double d[N];
	double x[N];*/
	int M = N-1;

	// calculate all a[i], c[i], d[i]

	// first segment
	{
		int i=0;

		double f_11 = 2./((X[i+2]-X[i+1])/(Y[i+2]-Y[i+1]) + (X[i+1]-X[i])/(Y[i+1]-Y[i]));
		if(SIGN((Y[i+2] - Y[i+1])/(X[i+2] - X[i+1])) != SIGN((Y[i+1] - Y[i])/(X[i+1] - X[i])))
			f_11 = 0.;
		double f_10 = 3./2.*(Y[i+1]-Y[i])/(X[i+1]-X[i]) - f_11/2.;
		double f__10 = -2.*(f_11 + 2.*f_10)/(X[i+1] - X[i]) + 6.*(Y[i+1] - Y[i])/ pow(X[i+1] - X[i],2.);
		double f__11 =  2.*(2.*f_11 + f_10)/(X[i+1] - X[i]) - 6.*(Y[i+1] - Y[i])/ pow(X[i+1] - X[i],2.);

		//fprintf( stderr, "f'1(X[i+1])=%lf, f'1(X[i])=%lf, f''1(X[i])=%lf, f''1(X[i+1])=%lf\n", f_11, f_10, f__10, f__11);

		d[i] = 1./6. * (f__11 - f__10)/(X[i+1] - X[i]);
		c[i] = 1./2. * (X[i+1]*f__10 - X[i]*f__11)/(X[i+1] - X[i]);
		b[i] = ((Y[i+1] - Y[i]) - c[i]*(X[i+1]*X[i+1] - X[i]*X[i])-d[i]*(X[i+1]*X[i+1]*X[i+1] - X[i]*X[i]*X[i]))/(X[i+1]-X[i]);
		a[i] = Y[i] - b[i]*X[i] - c[i]*X[i]*X[i] - d[i]*X[i]*X[i]*X[i];
		//fprintf( stderr, "a[i]=%lf, b[i]=%lf, c[i]=%lf, d[i]=%lf\n", a[i], b[i], c[i], d[i]);
	}

	for( int i=1; i<M-1; i++){

		double f_11 = 2./((X[i+2] - X[i+1])/(Y[i+2] - Y[i+1]) + (X[i+1] - X[i])/(Y[i+1] - Y[i]));
		if(SIGN((Y[i+2] - Y[i+1])/(X[i+2] - X[i+1])) != SIGN((Y[i+1] - Y[i])/(X[i+1] - X[i])))
			f_11 = 0.;//fprintf( stderr, "SIGN CHANGE\n");}
		double f_10 = 2./((X[i+1] - X[i])/(Y[i+1] - Y[i]) + (X[i] - X[i-1])/(Y[i] - Y[i-1]));
		if(SIGN((Y[i+1] - Y[i])/(X[i+1] - X[i])) != SIGN((Y[i] - Y[i-1])/(X[i] - X[i-1])))
			f_10 = 0.;//fprintf( stderr, "SIGN CHANGE\n");}
		double f__10 = -2.*(   f_11 + 2.*f_10)/(X[i+1] - X[i]) + 6.*(Y[i+1] - Y[i])/ pow(X[i+1] - X[i],2);
		double f__11 =  2.*(2.*f_11 +    f_10)/(X[i+1] - X[i]) - 6.*(Y[i+1] - Y[i])/ pow(X[i+1] - X[i],2);

		d[i] = 1/6. * (f__11 - f__10)/(X[i+1] - X[i]);
		c[i] = 1/2. * (X[i+1]*f__10 - X[i]*f__11)/(X[i+1] - X[i]);
		b[i] = ((Y[i+1] - Y[i]) - c[i]*(X[i+1]*X[i+1] - X[i]*X[i])-d[i]*(X[i+1]*X[i+1]*X[i+1] - X[i]*X[i]*X[i]))/(X[i+1]-X[i]);
		a[i] = Y[i] - b[i]*X[i] - c[i]*X[i]*X[i] - d[i]*X[i]*X[i]*X[i];
	}

	// last segment
	{
		int i=M-1;

		double f_10 = 2./((X[i+1]-X[i])/(Y[i+1]-Y[i]) + (X[i]-X[i-1])/(Y[i]-Y[i-1]));
		if(SIGN((Y[i+1] - Y[i])/(X[i+1] - X[i])) != SIGN((Y[i] - Y[i-1])/(X[i] - X[i-1])))
			f_10 = 0.;//fprintf( stderr, "SIGN CHANGE\n");}
		double f_11 = 3./2.*(Y[i+1]-Y[i])/(X[i+1]-X[i]) - f_10/2.;
		double f__10 = -2.*(f_11 + 2.*f_10)/(X[i+1] - X[i]) + 6.*(Y[i+1] - Y[i])/ pow(X[i+1] - X[i],2.);
		double f__11 =  2.*(2.*f_11 + f_10)/(X[i+1] - X[i]) - 6.*(Y[i+1] - Y[i])/ pow(X[i+1] - X[i],2.);

		//fprintf( stderr, "f'1(X[i+1])=%lf, f'1(X[i])=%lf, f''1(X[i])=%lf, f''1(X[i+1])=%lf\n", f_11, f_10, f__10, f__11);

		d[i] = 1./6. * (f__11 - f__10)/(X[i+1] - X[i]);
		c[i] = 1./2. * (X[i+1]*f__10 - X[i]*f__11)/(X[i+1] - X[i]);
		b[i] = ((Y[i+1] - Y[i]) - c[i]*(X[i+1]*X[i+1] - X[i]*X[i])-d[i]*(X[i+1]*X[i+1]*X[i+1] - X[i]*X[i]*X[i]))/(X[i+1]-X[i]);
		a[i] = Y[i] - b[i]*X[i] - c[i]*X[i]*X[i] - d[i]*X[i]*X[i]*X[i];
		//fprintf( stderr, "a[i]=%lf, b[i]=%lf, c[i]=%lf, d[i]=%lf\n", a[i], b[i], c[i], d[i]);
	}

	// output splines
	FILE *fp = fopen( "s3cs.dat", "w+");
	for(int i=0; i<M; i++)
		for( double x = X[i]; x <= X[i+1]; x+=(X[i+1]-X[i])/10.){

			//printf( "%10.3lf %10.3lf\n", X[i], Y[i]);
			//fprintf( fp, "%10.3lf %10.3lf\n", x, a[i]*pow(x - X[i],3.) + b[i]*pow(x - X[i],2.) + c[i]*(x - X[i]) + d[i]);
			//fprintf( fp, "%10.3lf %10.3lf\n", x, d[i]*pow(x,3.) + c[i]*pow(x,2.) + b[i]*(x) + a[i]);
			/*fprintf( fp, "%10.3lf %10.3lf\n", X,
				(x[i+1]*pow(X-X[i],3)+x[i]*pow(X[i+1]-X,3)) / 6. / (X[i+1]-X[i])
				+(Y[i+1] / (X[i+1]-X[i]) - (X[i+1]-X[i]) / 6. * x[i+1]) * (X-f[i  ][0])
				-(f[i  ][1] / (X[i+1]-X[i]) - (X[i+1]-X[i]) / 6. * x[i  ]) * (X-X[i+1])
			);*/
	}
	fclose( fp);

}

