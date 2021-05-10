/*
 * DiffusionReaction.cpp
 *
 *  Created on: Jul 20, 2011
 *      Author: jagiella
 */

#include "DiffusionReactionEquation.hpp"
#include "SparseMatrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

DiffusionReactionEquation::DiffusionReactionEquation(VoronoiDiagram *vd,
		char timeScheme, char reactionScheme, char diffusionScheme) {
	// TYPE OF SCHEME
	this->timeScheme 		= timeScheme;
	this->reactionScheme 	= reactionScheme;//EXPLICIT;
	this->diffusionScheme 	= diffusionScheme;
	this->boundaryCondition = DiffusionReactionEquation::DIRICHLET;

	// SPATIAL DISCRITIZATION
	this->lattice = vd;

	// LINEAR SYSTEM
	this->b = 0;
	this->x = 0;
	this->A = 0;
	this->temp = 0;

	// PARAMETERS
	this->reaction = 0;
	this->diffusion = 0;
	this->timeStep = 1.;
	this->timeLeft = 0.;
	this->spaceStep = 1;

	maxIterations = 1000;
	maxError    = 1e-4;
	maxStepSize = 1e-8;

	this->systemUpToDate = false;
}

void DiffusionReactionEquation::update(float duration) {
	this->timeLeft += duration;

	while (this->timeLeft >= this->timeStep) {
		update();

		if (timeScheme == DiffusionReactionEquation::STEADYSTATE)
			this->timeLeft = 0;
		else
			this->timeLeft -= this->timeStep;
	}
}

void DiffusionReactionEquation::update() {
	if (!systemUpToDate)
		assembleSystem();

	if (timeScheme == DiffusionReactionEquation::STEADYSTATE){
		int it, it_total = 0;
		if( temp==0)
			temp = (float*)malloc( sizeof(float) * this->lattice->countVoronoiCells);
		float max_err = 0;
		do{
			for(int i=0; i<this->lattice->countVoronoiCells; i++)
				temp[i] = x[i];

			// solve system
			it = SolveBiCGSTAB(A, b, x, maxIterations, maxError);
			it_total += it;

			// calculate max error / norm L^inf
			max_err = 0;
			for(int i=0; i<this->lattice->countVoronoiCells; i++)
				if( fabs( temp[i] - x[i]) > max_err){
				     max_err = fabs( temp[i] - x[i]);
				}


		}while( it > 1 && it_total <= maxIterations && max_err > maxStepSize);
	}
	else
		SolveBiCGSTAB(A, b, x, maxIterations, maxError);
}

void DiffusionReactionEquation::assembleSystem() {
	int di = 1;
#if DIMENSIONS >=2
	int dii = lattice->xN[0];
#endif
#if DIMENSIONS >= 3
	int diii = lattice->xN[0] * lattice->xN[1];
#endif

	for (int i = 0; i < lattice->countVoronoiCells; i++) {
		// reset column
		A->resetRow(i);
		b[i] = 0.;

		if (	this->boundaryCondition == DiffusionReactionEquation::NEUMANN
				||(
				   (int) lattice->voronoiCells[i]->position[0] > 0
				&& (int) lattice->voronoiCells[i]->position[0] < lattice->xN[0] - 1
#if DIMENSIONS >= 2
				&& (int) lattice->voronoiCells[i]->position[1] > 0
				&& (int) lattice->voronoiCells[i]->position[1] < lattice->xN[1] - 1
#endif
#if DIMENSIONS >= 3
				&& (int) lattice->voronoiCells[i]->position[2] > 0
				&& (int) lattice->voronoiCells[i]->position[2] < lattice->xN[2] - 1
#endif
						)) {

			// reaction
			double diagonal = 0.;
			double vector = 0.;
			switch (reactionScheme) {
			case DiffusionReactionEquation::EXPLICIT:
				// ... = ... - r_n
				vector -= reaction[i];
				break;
			case DiffusionReactionEquation::IMPLICIT:
				// ... + r_n+1 = ...
				//if( x[i]>0)
				//	diagonal += reaction[i] / x[i];
				//else
				//	vector -= reaction[i];
				diagonal += reactionDerivative[i];
				vector -= reaction[i];
				break;
			case DiffusionReactionEquation::CRANKNICHOLSON:
				diagonal += 0.5 * reaction[i] / x[i];
				vector -= 0.5 * reaction[i];
				break;
			default:
				fprintf(stderr, "ERROR: The used reaction scheme (%i) is not implemented!\n", (int)reactionScheme);
				exit(0);
			}

			// time derivative
			switch (timeScheme) {
			case DiffusionReactionEquation::STEADYSTATE:
				// (C_n+1 - C_n) / dt = 0
				break;
			default:
				// (C_n+1 - C_n) / dt = ...
				diagonal += 1. / timeStep;
				vector += x[i] / timeStep;
				break;
			}

			// diffusion
			double r = 1. / (spaceStep * spaceStep);
			switch (diffusionScheme) {
			case DiffusionReactionEquation::IMPLICIT:{
				/*diagonal += DIMENSIONS * 2. * r * this->diffusion[i];
				A->setLast(i, i - diii, -r / diagonal * this->diffusion[i - diii]);
				A->setLast(i, i - dii,	-r / diagonal * this->diffusion[i - dii]);
				A->setLast(i, i - di,   -r / diagonal * this->diffusion[i - di]);
				A->setLast(i, i, 1.);
				A->setLast(i, i + di,   -r / diagonal * this->diffusion[i + di]);
				A->setLast(i, i + dii,  -r / diagonal * this->diffusion[i + dii]);
				A->setLast(i, i + diii, -r / diagonal * this->diffusion[i + diii]);
				vector /= diagonal;*/

				//diagonal += DIMENSIONS * 2. * r * this->diffusion[i];
				float diffusionSum=0.;
#if DIMENSIONS >= 3
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-3] > 0){
					A->setLast(i, i - diii, -r * (this->diffusion[i - diii]+this->diffusion[i])*0.5);
					diffusionSum+=(this->diffusion[i - diii]+this->diffusion[i])*0.5;
				}
#endif
#if DIMENSIONS >= 2
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-2] > 0){
					A->setLast(i, i - dii,	-r * (this->diffusion[i - dii]+this->diffusion[i])*0.5);
					diffusionSum+=(this->diffusion[i - dii]+this->diffusion[i])*0.5;
				}
#endif
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-1] > 0){
					A->setLast(i, i - di,   -r * (this->diffusion[i - di]+this->diffusion[i])*0.5);
					diffusionSum+=(this->diffusion[i - di]+this->diffusion[i])*0.5;
				}
				//A->setLast(i, i, diagonal);
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-1] < lattice->xN[DIMENSIONS-1] - 1){
					A->setLast(i, i + di,   -r * (this->diffusion[i + di]+this->diffusion[i])*0.5);
					diffusionSum+=(this->diffusion[i + di]+this->diffusion[i])*0.5;
				}
#if DIMENSIONS >= 2
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-2] < lattice->xN[DIMENSIONS-2] - 1){
					A->setLast(i, i + dii,  -r * (this->diffusion[i + dii]+this->diffusion[i])*0.5);
					diffusionSum+=(this->diffusion[i + dii]+this->diffusion[i])*0.5;
				}
#endif
#if DIMENSIONS >= 3
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-3] < lattice->xN[DIMENSIONS-3] - 1){
					A->setLast(i, i + diii, -r * (this->diffusion[i + diii]+this->diffusion[i])*0.5);
					diffusionSum+=(this->diffusion[i + diii]+this->diffusion[i])*0.5;
				}
#endif
				diagonal += r * diffusionSum;
				A->set(i, i, diagonal);

				b[i] = vector;
				break;
			}
			case DiffusionReactionEquation::EXPLICIT:{

				int a=0;
#if DIMENSIONS >= 3
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-3] > 0){
					vector +=  r * this->diffusion[i - diii] * x[i - diii];
					a++;
				}
#endif
#if DIMENSIONS >= 2
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-2] > 0){
					vector +=  r * this->diffusion[i - dii]  * x[i - dii];
					a++;
				}
#endif
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-1] > 0){
					vector +=  r * this->diffusion[i - di]   * x[i - di];
					a++;
				}

				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-1] < lattice->xN[DIMENSIONS-1] - 1){
					vector +=  r * this->diffusion[i + di]   * x[i + di];
					a++;
				}
#if DIMENSIONS >= 2
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-2] < lattice->xN[DIMENSIONS-2] - 1){
					vector +=  r * this->diffusion[i + dii]  * x[i + dii];
					a++;
				}
#endif
#if DIMENSIONS >= 3
				if((int) lattice->voronoiCells[i]->position[DIMENSIONS-3] < lattice->xN[DIMENSIONS-3] - 1){
					vector +=  r * this->diffusion[i + diii] * x[i + diii];
					a++;
				}
#endif
				vector += -r * this->diffusion[i]        * x[i] * a;

				//A->setLast(i, i, diagonal);
				//b[i] = vector;
				A->setLast(i, i, 1.);
				b[i] = vector/diagonal;
			}
				break;
			default:
				fprintf(stderr, "ERROR: The used diffusion scheme (%i) is not implemented!\n", (int)diffusionScheme);
				exit(0);
			}
		}

		else {
			A->setLast(i, i, 1.);
			b[i] = x[i];
		}
	}

	for(int i=0; i<A->dimI; i++){
		float aii = A->get( i,i);
		for( int jj=0; jj<A->sizeA[i]; jj++){
			//int j=sA->JA[i][jj];
			A->A[i][jj] /= aii;
		}
		b[i] /= aii;
	}

}

void DiffusionReactionEquation::initReaction(float *reaction) {
	this->reaction = reaction;
}

void DiffusionReactionEquation::initReactionDerivative(float *reaction) {
	this->reactionDerivative = reaction;
}

void DiffusionReactionEquation::initDiffusion(float *diffusion) {
	this->diffusion = diffusion;
}

void DiffusionReactionEquation::initSystem(SparseMatrix* A, float *b,
		float *x) {
	this->A = A;
	this->b = b;
	this->x = x;
}

void DiffusionReactionEquation::initSystem(float *x) {
	int N = this->lattice->countVoronoiCells;
	this->A = SparseMatrix::newSparseMatrix(N, N);
	this->b = (float*) malloc(sizeof(float) * N);
	this->x = x;
}

void DiffusionReactionEquation::initSystem() {
	int N = this->lattice->countVoronoiCells;
	this->A = SparseMatrix::newSparseMatrix(N, N);
	this->b = (float*) malloc(sizeof(float) * N);
	this->x = (float*) malloc(sizeof(float) * N);
}

void DiffusionReactionEquation::setValues( int index, float value)
{
	x[index] = value;
}

void DiffusionReactionEquation::setTimeStep( float timeStep)
{
	this->timeStep = timeStep;
}

void DiffusionReactionEquation::setBoundaryCondition( char bc)
{
	this->boundaryCondition = bc;
}

void DiffusionReactionEquation::setSpaceStep(float spaceStep)
{
	this->spaceStep = spaceStep;
}

float DiffusionReactionEquation::getTimeStep()
{
	return this->timeStep;
}

float DiffusionReactionEquation::getSpaceStep()
{
	return this->spaceStep;
}
