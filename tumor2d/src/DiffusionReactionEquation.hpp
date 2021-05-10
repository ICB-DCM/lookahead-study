/*
 * DiffusionReactionEquation.hpp
 *
 *  Created on: Jul 22, 2011
 *      Author: jagiella
 */

#ifndef DIFFUSIONREACTIONEQUATION_HPP_
#define DIFFUSIONREACTIONEQUATION_HPP_


#include "VoronoiDiagramExtended.h"
#include "SparseMatrix.h"


class DiffusionReactionEquation{
public:
	// SCHEMES
	static const char EXPLICIT		= 1;
	static const char IMPLICIT		= 2;
	static const char STEADYSTATE	= 3;
	static const char CRANKNICHOLSON= 4;

	// BOUNDARY CONDITION
	static const char NEUMANN		= 1;
	static const char DIRICHLET		= 2;

private:
	VoronoiDiagram *lattice;

	float *temp;
	float *x;
	float *b;
	SparseMatrix *A;
	bool systemUpToDate;

	char reactionScheme;
	float *reaction;
	float *reactionDerivative;
	char diffusionScheme;
	float *diffusion;
	char timeScheme;
	float timeLeft;
	float timeStep;
	float spaceStep;
	char boundaryCondition;

	int maxIterations;
	float maxError;
	float maxStepSize;

public:
	DiffusionReactionEquation( VoronoiDiagram *vd, char timeScheme, char reactionScheme, char diffusionScheme);
	~DiffusionReactionEquation();

	void initSystem();
	void initSystem( float *x);
	void initSystem( SparseMatrix* A, float *b, float *x);

	void setTimeStep( float);
	void setSpaceStep(float);
	void setBoundaryCondition( char bc);

	float getTimeStep();
	float getSpaceStep();

	void initValues();
	void initReaction();
	void initDiffusion();
	void initValues(    float *values);
	void initReaction(  float *reaction);
	void initReactionDerivative(  float *reaction);
	void initDiffusion( float *diffusion);

	void setValues(    int index, float values);
	void setReaction(  int index, float &reaction);
	void setDiffusion( int index, float &diffusion);
	float *getValues();
	float *getReaction();
	float *getDiffusion();

	void assembleSystem();

	void update( float duration);
	void update();
};

#endif /* DIFFUSIONREACTIONEQUATION_HPP_ */
