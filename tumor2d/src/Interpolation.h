/*
 * Interpolation.h
 *
 *  Created on: 22.02.2010
 *      Author: jagiella
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "VoronoiDiagramExtended.h"

class Interpolation
{
public:
	static VoronoiDiagram *voronoiDiagram;
	static float getRectangleAreaAboveThreshold3D( int index, char mode, float threshold);
	static float getDiscritizedVolumeAboveThreshold1D( int index, char mode, float threshold);
	static float getDiscritizedVolumeAboveThreshold2D( int index, char mode, float threshold);
	static float getDiscritizedVolumeAboveThreshold3D( int index, char mode, float threshold);
	static float getDiscritizedVolumeAboveThreshold1DCubicSpline( int index, char mode, float threshold);
};


double getVolume( double *p0, double *p1, double *p2, double *p3);
void cubicSplines1Dcubic( int N, double *X, double *Y, double *a, double *b, double *c, double *d);
void cubicSplines1Dconstraint( int N, double *X, double *Y, double *a, double *b, double *c, double *d);
#endif /* INTERPOLATION_H_ */
