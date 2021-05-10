/****************************************************************************
 * Includes                                                                 *
 ****************************************************************************/
#include "Montecarlo.h"

#include <stdio.h>            // Standard-Ein- und Ausgabe
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <vector>


#include "VoronoiDiagramExtended.h"
#include "SquareLattice2D.h"
#include "HexagonalLattice2D.h"
#include "Agent.h"
#include "CellProcesses.h"
#include "Mathematix.h"
#include "Substrate.h"
#include "VesselNetwork.h"
#include "finiteDifferences.h"
#include "SparseMatrix.h"
#include "Interpolation.h"

#include "DiffusionReactionEquation.hpp"

#include "matrix.hpp"
#include "statistics.hpp"
#include "AsciiIO.h"
#include <time.h>

#define MAX_CELLS 200000
#define PROB_REENTERING_CELL_CYCLE( DELTA_L) ( ReentranceProbabilityFunction == HEAVYSIDE ? (DELTA_L <= ReentranceProbabilityReferenceLength ? 1 : 0) : exp(-DELTA_L/ReentranceProbabilityReferenceLength))

double TempTime = 0;

// TEXT STYLE
#define RESET		0
#define BRIGHT 		1
#define DIM			2
#define UNDERLINE 	3
#define BLINK		4
#define REVERSE		7
#define HIDDEN		8

// TEXT COLORS
#define BLACK 		0
#define RED			1
#define GREEN		2
#define YELLOW		3
#define BLUE		4
#define MAGENTA		5
#define CYAN		6
#define	WHITE		7



int Case = 1;

typedef struct _Statisticx {
	int *count_cells;
	int *count_cell_volume;
	int *count_expanded_cells;
	double *count_divided_cells;
	int *count_vessel_cells;
	int *count_necrotic_cells;
	int *count_necrotic_volume;
	int *count_inner_empty_volume;
	int *count_lysed_cells;
	int *count_lysed_volume;
} Statistix;


void updateHistogram( AgentList *agentArray, VoronoiDiagram *vd, int* histogram, int* histogramDividing, int* histogramNecrotic, int* histogramFree, double* histogramECM)
{
	//Time = EndTime;

	// GET ALL BORDER CELLS
	int countBorderCells = 0;
	int maxBorderCells = 100000;

	VoronoiCell *borderCells[100000];
	//fprintf(stderr, "[Estimate Border]\n");

	for (int a = 0; countBorderCells!=maxBorderCells && a < agentArray->countActiveAgents; a++)
		for (int l = 0; countBorderCells!=maxBorderCells && l < agentArray->agents[a]->countLocations; l++){
			//bool unoccupiedNeighbor = false;
			for( int n=0; countBorderCells!=maxBorderCells && n<agentArray->agents[a]->location[l]->countNeighborCells /*&& !unoccupiedNeighbor*/; n++)
				if( !agentArray->agents[a]->location[l]->neighborCells[n]->agent
						|| agentArray->agents[a]->location[l]->neighborCells[n]->index < 0
						|| agentArray->agents[a]->location[l]->isDomainBorder(vd)){

					//unoccupiedNeighbor=true;
					//if( countBorderCells==100000)
					//	exit(0);
					borderCells[countBorderCells++] = agentArray->agents[a]->location[l]->neighborCells[n];
				}
		}
	//fprintf(stderr, "[Border Estimated]\n");

	for (int a = 0; a < agentArray->countActiveAgents; a++)
	//if( agentArray->agents[a]->state != FREE)
			{
		//fprintf(stderr, "\r[Agent %i]: State                 \b", a );
		// DIVIDING?
		bool dividing = false;
		bool necrotic = false;
		bool free = false;
		for (int l = 0;
				l < agentArray->agents[a]->countLocations;
				l++) {
			if (agentArray->agents[a]->state == FREE)
				free = true;
			if (agentArray->agents[a]->state == NECROTIC)
				necrotic = true;
			if (agentArray->agents[a]->state == ACTIVE
					&& agentArray->agents[a]->divide == 1
					&& !IsQuiescent(agentArray->agents[a]))
				dividing = true;
		}

		// distance to border
		VoronoiCell *border = 0;

		//fprintf(stderr, "\r[Agent %i]: Dist to Border                 \b", a );
		if(countBorderCells){
			border = borderCells[0];
			double distClosestBorder = border->getDistanceTo( agentArray->agents[a]->location[0]);
			for( int b=1; b<countBorderCells; b++){
				double dist = borderCells[b]->getDistanceTo( agentArray->agents[a]->location[0]);
				if( distClosestBorder>dist){
					border = borderCells[b];
					distClosestBorder=dist;
				}
			}
		}else{
			//borderCells[]
		}

		double dist = 1000000;
		if (border) {
			dist = border->getDistanceTo(
					agentArray->agents[a]->location[0])
					* SPATIAL_UNIT;
		} else {
			//fprintf(stderr, "No border found!\n");
			// shortest distance to domain border
			for( int i=0; i<DIMENSIONS; i++){
				dist = fmin( dist, agentArray->agents[a]->location[0]->position[i]);
				dist = fmin( dist, 100 /*pow(agentArray->countAgents, 1./DIMENSIONS)*/ - agentArray->agents[a]->location[0]->position[i]);
			}
			dist *= SPATIAL_UNIT;
		}

		for (int i = (dist - AGENT_DIAMETER > 0. ? (int) (dist - AGENT_DIAMETER) : 0);
				i < (int) (dist); i++) {
			histogram[i]++;
			if (dividing)
				histogramDividing[i]++;
			if (necrotic)
				histogramNecrotic[i]++;
			if (free)
				histogramFree[i]++;
			histogramECM[i] +=
					agentArray->agents[a]->location[0]->ecm;
		}

	}
	//fprintf(stderr, "[finished]\n");

}

double GetDistanceToClosestFreeNeighbor(VoronoiDiagram *voronoiDiagram,
		VoronoiCell *voronoiCell) {
	int dist = 0;
	int min_dist = 0;
	for (int d = 0; d < DIMENSIONS; d++)
		min_dist += voronoiDiagram->xN[d];

	if (VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD) {
		VoronoiCell *closestVC;
		if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
			closestVC =
					voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
							voronoiCell, VoronoiCell::extendedNeighborDepth,
							VoronoiCell::symbolicExtendedNeighborhoodSize,
							VoronoiCell::symbolicExtendedNeighborhood);
		else
			closestVC =
					voronoiDiagram->searchClosestFreeVoronoiCell(voronoiCell,
							VoronoiCell::extendedNeighborDepth,
							VoronoiCell::symbolicExtendedNeighborhoodSize,
							VoronoiCell::symbolicExtendedNeighborhood);
		if (closestVC)
			min_dist = closestVC->getDistanceSquareTo(voronoiCell);
	} else if (VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD) {
		VoronoiCell *closestVC;
		if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
			closestVC =
					voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
							voronoiCell, VoronoiCell::extendedNeighborDepth);
		else
			closestVC =
					voronoiDiagram->searchClosestFreeVoronoiCell(voronoiCell,
							VoronoiCell::extendedNeighborDepth);

		if (closestVC)
			min_dist = closestVC->getDistanceSquareTo(voronoiCell);
	} else {
		if (voronoiCell->countFreeNeighborCells)
			for (int i = 0; i < voronoiCell->countNeighborCells; i++) {
				dist = 0;
				for (int d = 0; d < DIMENSIONS; d++)
					dist +=
							pow(
									voronoiCell->position[d]
											- voronoiCell->neighborCells[i]->position[d],
									2);
				if (min_dist > dist)
					min_dist = dist;
			}

		if (voronoiCell->countFreeExtendedNeighborCells)
			for (int i = 0; i < voronoiCell->countExtendedNeighborCells; i++) {
				dist = 0;
				for (int d = 0; d < DIMENSIONS; d++)
					dist +=
							pow(
									voronoiCell->position[d]
											- voronoiCell->extendedNeighborhood[i]->position[d],
									2);
				if (min_dist > dist)
					min_dist = dist;
			}
	}

	return sqrt(min_dist);
}

void refineSurrounding(VoronoiDiagram *voronoiDiagram, VoronoiCell *cell,
		ActionTree *actionList, AgentList *agentArray, int scale){
	int index; // = floor(cell->position[0])*voronoiDiagram->xN[1] + floor(cell->position[1]);
	int x = (int) floor(cell->position[0]);
	int y = (int) floor(cell->position[1]);
#if DIMENSIONS == 3
	int z = (int)floor(cell->position[2]);
	//int dz = 1;
	int dy = voronoiDiagram->xN[2];
	int dx = voronoiDiagram->xN[1] * dy;
	int minz = MAX(0,z-1);
	int maxz = MIN(voronoiDiagram->xN[2]-1,z+1);
#elif DIMENSIONS == 2
	int dy = 1;
	int dx = voronoiDiagram->xN[1];
#endif
	int n = 0;
	//for( int ix=MAX(0,x-1); ix<MIN(voronoiDiagram->xN[0],x+2); ix++)
	//for( int iy=MAX(0,y-1); iy<MIN(voronoiDiagram->xN[1],y+2); iy++)
	int minx = MAX(0,x-1);
	int miny = MAX(0,y-1);
	int maxx = MIN(voronoiDiagram->xN[0]-1,x+1);
	int maxy = MIN(voronoiDiagram->xN[1]-1,y+1);

	for (int ix = minx; ix <= maxx; ix++)
		for (int iy = miny; iy <= maxy; iy++)
#if DIMENSIONS == 2
				{
			index = ix * dx + iy;
#elif DIMENSIONS == 3
			for( int iz=minz; iz<=maxz; iz++)
			//for( int iz=MAX(0,z-1); iz<MIN(voronoiDiagram->xN[2],z+2); iz++)
			{
				index = ix*dx +iy*dy +iz;
#endif
			//{
			//index = ix*dx +iy*dy;
			if (voronoiDiagram->voronoiCells[index]->refined == false
					&& voronoiDiagram->voronoiCells[index]->isFree()) {

				//long passedTime = clock();
				VoronoiCell* comp = voronoiDiagram->voronoiCells[index];
				voronoiDiagram->refine(voronoiDiagram->voronoiCells[index],
						scale, actionList);
				//fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);

				/*int n=0;
				 for( ; n<GetVoronoiCell( selected_action->originalCell)->countNeighborCells; )
				 {
				 if( GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->index >= 0 &&
				 GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->refined == false
				 && GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->isFree()){
				 //printTriangulation( voronoiDiagram, "beforeRef.eps", 1);
				 passedTime = clock();
				 VoronoiCell* comp = GetVoronoiCell( selected_action->originalCell)->neighborCells[n];
				 voronoiDiagram->refine( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], CountCellsPerVoronoiCell, actionList);
				 fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);
				 //printTriangulation( voronoiDiagram, "afterRef.eps", 1);
				 //voronoiDiagram->coarsen( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], 10);
				 n=0;*/

				//TEST
				if (comp->agent == NULL)
				//exit(0);
				{
					//fprintf(stderr, "Add Compartment Agent\n");

					Agent *compAgent = agentArray->activateAgent();
					compAgent->attach(comp);
					compAgent->state = COMPARTMENT; //NONACTIVE;
					compAgent->cellCount = 0;
					compAgent->maxCellCount = scale * scale
#if DIMENSIONS == 3
					*scale
#endif
					;
					compAgent->cellCount = compAgent->maxCellCount;
					compAgent->countFree = compAgent->maxCellCount;
					compAgent->countActive = 0;
					compAgent->countNonactive = 0;
					//for( int n=0; n<comp->countNeighborCells; n++)
					//	comp->neighborCells[n]->countFreeNeighborCells--;
				}
				//TEST END
			} else {
				n++;
			}
#if DIMENSIONS >= 2
		}
#endif
}

void refineNeighborhood(VoronoiDiagram *voronoiDiagram, VoronoiCell *cell,
		ActionTree *actionList, AgentList *agentArray, int scale) {
	int n = 0;
	for (; n < cell->countNeighborCells;) {
		if (cell->neighborCells[n]->index >= 0
				&& cell->neighborCells[n]->refined == false
				&& cell->neighborCells[n]->isFree()) {
			//printTriangulation( voronoiDiagram, "beforeRef.eps", 1);
			//long passedTime = clock();
			VoronoiCell* comp = cell->neighborCells[n];
			voronoiDiagram->refine(cell->neighborCells[n],
					CountCellsPerVoronoiCell, actionList);
			//fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);
			//printTriangulation( voronoiDiagram, "afterRef.eps", 1);
			//voronoiDiagram->coarsen( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], 10);
			n = 0;

			//TEST
			if (comp->agent == NULL)
			//exit(0);
			{
				//fprintf(stderr, "Add Compartment Agent\n");

				Agent *compAgent = agentArray->activateAgent();
				compAgent->attach(comp);
				compAgent->state = COMPARTMENT; //NONACTIVE;
				compAgent->cellCount = 0;
				compAgent->maxCellCount = scale * scale
#if DIMENSIONS == 3
				*scale
#endif
				;
				compAgent->cellCount = compAgent->maxCellCount;
				compAgent->countFree = compAgent->maxCellCount;
				compAgent->countActive = 0;
				compAgent->countNonactive = 0;
				//for( int n=0; n<comp->countNeighborCells; n++)
				//	comp->neighborCells[n]->countFreeNeighborCells--;
			}
			//TEST END
		} else {
			n++;
		}
	}
}




void performMigration(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);

void performChemotacticMigration(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);

VoronoiCell* performGrowth(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);
void performDivision(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);
VoronoiCell* performGrowthAndDivision(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);


double getAvgGlucose(AgentList *agentList) {
	double avg = 0.;
	double N = 0;

	for (int i = 0; i < agentList->countActiveAgents; i++) {
		if (agentList->agents[i]->growingTumorCellCount > 0){
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++){
				if (agentList->agents[i]->state == COMPARTMENT) {
					avg +=
							agentList->agents[i]->location[ii]->glucose
									* agentList->agents[i]->growingTumorCellCount;
					N += agentList->agents[i]->growingTumorCellCount;
				} else {
					avg += agentList->agents[i]->location[ii]->glucose;
					N++;
				}
			}
		}
	}

	return (N > 0. ? avg / N : 0.);
}

double getAvgOxygen(AgentList *agentList) {
	double avg = 0.;
	double N = 0;

	for (int i = 0; i < agentList->countActiveAgents; i++) {
		if (agentList->agents[i]->growingTumorCellCount > 0){
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++){
				if (agentList->agents[i]->state == COMPARTMENT) {
					if (N == 0
							|| avg
									> agentList->agents[i]->location[ii]->glucose
											* agentList->agents[i]->location[ii]->oxygen) {
						avg =
								agentList->agents[i]->location[ii]->glucose
										* agentList->agents[i]->location[ii]->oxygen;
						//			if( N==0 || avg > agentList->agents[i]->location[ii]->glucose){
						//				avg = agentList->agents[i]->location[ii]->glucose;
						N = 1;
					}
					//			avg += agentList->agents[i]->location[ii]->oxygen * agentList->agents[i]->growingTumorCellCount;
					//			N += agentList->agents[i]->growingTumorCellCount;
				} else {
					if (N == 0
							|| avg
									> agentList->agents[i]->location[ii]->glucose
											* agentList->agents[i]->location[ii]->oxygen) {
						avg =
								agentList->agents[i]->location[ii]->glucose
										* agentList->agents[i]->location[ii]->oxygen;
						//			if( N==0 || avg > agentList->agents[i]->location[ii]->glucose){
						//				avg = agentList->agents[i]->location[ii]->glucose;
						N = 1;
					}
					//			avg += agentList->agents[i]->location[ii]->oxygen;
					//			N ++;
				}
			}
		}
	}

	return (N > 0. ? avg / N : 0.);
}

VoronoiCell **getOuterBorder(AgentList *agentList, int &borderSize) {
	int max_queueSize = agentList->countActiveAgents
			* (MAX_SUBCELLULAR_COMPONENTS) * 2;
	int queueSize = 0;
	queueSize = 0;
	VoronoiCell **queue = (VoronoiCell**) calloc(max_queueSize,
			sizeof(VoronoiCell*));

	//fprintf( stderr, "ENTER!\n");

	//fprintf( stderr, "agentList->countActiveAgents: %i\n", agentList->countActiveAgents);

	// find first border location
	for (int i = 0; !queueSize && i < agentList->countActiveAgents; i++) {
		//	fprintf( stderr, "agent: %i\n", agentList->agents[i]->index);
		if (agentList->agents[i]->state != FREE)
			for (int ii = 0;
					!queueSize && ii < agentList->agents[i]->countLocations;
					ii++) {
				//		fprintf( stderr, "-->location: %i\n", agentList->agents[i]->location[ii]->index);
				for (int iii = 0;
						!queueSize
								&& iii
										< agentList->agents[i]->location[ii]->countNeighborCells;
						iii++) {
					//			fprintf( stderr, "---->neighbor: %i\n", agentList->agents[i]->location[ii]->neighborCells[iii]->index);
					if (GetAgent( agentList->agents[i]->location[ii]->neighborCells[iii])
							== NULL) {
						if (max_queueSize == queueSize) {
							fprintf(
									stderr,
									"WARNING: realloc queue: %i->%i\n",
									max_queueSize,
									max_queueSize + agentList->countActiveAgents);
							max_queueSize += agentList->countActiveAgents;
							queue =
									(VoronoiCell**) realloc(
											queue,
											max_queueSize
													* sizeof(VoronoiCell *));
						}
						queue[queueSize++] =
								agentList->agents[i]->location[ii]->neighborCells[iii];
						//				fprintf( stderr, "    >ADD TO QUEUE!\n");
					}
				}
			}
	}
	//fprintf( stderr, "AHA?!\n");
	//fprintf( stderr, "queueSize: %i, queue[0]: %i\n", queueSize, queue[0]->index);
	if (queueSize == 0)
		return NULL;
	//fprintf( stderr, "AHA!!!!\n");

	// collect all border points
	int max_borderSize = agentList->countActiveAgents
			* (MAX_SUBCELLULAR_COMPONENTS) * 2;
	borderSize = 0;
	VoronoiCell **border = (VoronoiCell**) malloc(
			max_borderSize * sizeof(VoronoiCell*));
	do {
		// take next element from queue
		VoronoiCell *actual = queue[--queueSize];

		// test if containes empty neighbors
		//fprintf( stderr, "-->next: %i\n", actual->index);
		char occupiedNeighborFound = FALSE;
		for (int iii = 0;
				!occupiedNeighborFound && iii < actual->countNeighborCells;
				iii++) {
			//	fprintf( stderr, "---->neighbor: %i\n", actual->neighborCells[iii]->index);
			// test if containes empty neighbors
			if (GetAgent( actual->neighborCells[iii]) != NULL)
				occupiedNeighborFound = TRUE;
		}

		if (occupiedNeighborFound) {
			// add neighbors to queue
			for (int iii = 0; iii < actual->countNeighborCells; iii++) {
				if (GetAgent( actual->neighborCells[iii]) == NULL
						&& !isElementOf(border, borderSize,
								actual->neighborCells[iii])
						&& !isElementOf(queue, queueSize,
								actual->neighborCells[iii])) {

					if (max_queueSize == queueSize) {
						fprintf(stderr, "WARNING: realloc queue: %i->%i\n",
								max_queueSize,
								max_queueSize + agentList->countActiveAgents);
						max_queueSize += agentList->countActiveAgents;
						queue =
								(VoronoiCell**) realloc(queue,
										max_queueSize * sizeof(VoronoiCell *));
					}

					queue[queueSize++] = actual->neighborCells[iii];
					//	fprintf( stderr, "    >ADD TO QUEUE!\n");
				}
			}

			// add to border
			//if( max_borderSize==borderSize)
			//	exit( 0);
			if (max_borderSize == borderSize) {
				fprintf(stderr, "WARNING: realloc border: %i->%i\n",
						max_borderSize,
						max_borderSize + agentList->countActiveAgents);
				max_borderSize += agentList->countActiveAgents;
				border =
						(VoronoiCell**) realloc(border,
								max_borderSize * sizeof(VoronoiCell *));
			}
			border[borderSize++] = actual;
		}
	} while (queueSize > 0);

	free(queue);

	//fprintf( stderr, "border points:");		
	for (int i = 0; i < borderSize; i++) {
		//fprintf( stderr, " %i", border[i]->index);		
	}
	//fprintf( stderr, "\n");
	//fprintf( stderr, "LEAVE!\n");

	return border;

	free(border);
	borderSize = 0;
	return (VoronoiCell**) malloc(1 * sizeof(VoronoiCell*));
}

VoronoiCell **getOuterBorder2(AgentList *agentList, int &borderSize) {
	int max_borderSize = agentList->countActiveAgents
			* (MAX_SUBCELLULAR_COMPONENTS) * 2;
	borderSize = 0;
	VoronoiCell **border = (VoronoiCell**) malloc(
			max_borderSize * sizeof(VoronoiCell*));

	//fprintf( stderr, "ENTER!\n");

	//fprintf( stderr, "agentList->countActiveAgents: %i\n", agentList->countActiveAgents);

	// find first border location
	for (int i = 0; i < agentList->countActiveAgents; i++) {
		//	fprintf( stderr, "agent: %i\n", agentList->agents[i]->index);
		if (agentList->agents[i]->state == ACTIVE)
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++) {
				//		fprintf( stderr, "-->location: %i\n", agentList->agents[i]->location[ii]->index);
				for (int iii = 0;
						iii
								< agentList->agents[i]->location[ii]->countNeighborCells;
						iii++) {
					//			fprintf( stderr, "---->neighbor: %i\n", agentList->agents[i]->location[ii]->neighborCells[iii]->index);
					if (GetAgent( agentList->agents[i]->location[ii]->neighborCells[iii])
							== NULL
							&& !isElementOf(
									border,
									borderSize,
									agentList->agents[i]->location[ii]->neighborCells[iii])) {
						if (max_borderSize == borderSize) {
							fprintf(
									stderr,
									"WARNING: realloc queue: %i->%i\n",
									max_borderSize,
									max_borderSize
											+ agentList->countActiveAgents);
							max_borderSize += agentList->countActiveAgents;
							border =
									(VoronoiCell**) realloc(
											border,
											max_borderSize
													* sizeof(VoronoiCell *));
						}
						border[borderSize++] =
								agentList->agents[i]->location[ii]->neighborCells[iii];
						//				fprintf( stderr, "    >ADD TO QUEUE!\n");
					}
				}
			}
	}
	//fprintf( stderr, "AHA?!\n");
	//fprintf( stderr, "queueSize: %i, queue[0]: %i\n", queueSize, queue[0]->index);
	//if( queueSize == 0)
	//	return NULL;
	//fprintf( stderr, "AHA!!!!\n");

	return border;
}

double getGyrationRadius(AgentList *agentList) {
	// center of mass
	//fprintf( stderr, "center of mass:");
	double centerMass[DIMENSIONS];
	int count=agentList->countActiveAgents;
	for (int i = 0; i < DIMENSIONS; i++) {
		centerMass[i] = 0.;
		count=0;
		for (int ii = 0; ii < agentList->countActiveAgents; ii++)
			if( agentList->agents[ii]->state != FREE){
			centerMass[i] += agentList->agents[ii]->location[0]->position[i];
			count++;
		}
		centerMass[i] /= (double) count;
		//	fprintf( stderr, " %lf", centerMass[i]);
	}
	//fprintf( stderr, "\n");

	// radius of gyration
	double gyrRadius = 0.;
	for (int ii = 0; ii < agentList->countActiveAgents; ii++)
		if( agentList->agents[ii]->state != FREE){
		for (int i = 0; i < DIMENSIONS; i++) {
			gyrRadius += pow(agentList->agents[ii]->location[0]->position[i] - centerMass[i], 2.);
		}
	}
	gyrRadius /= (double) count;
	//fprintf( stderr, "radius of gyration: %lf\n", gyrRadius);

	return gyrRadius;
}

double getGyrationRadiusOfBorder(AgentList *agentList, VoronoiDiagram *vd) {
	// collect border points
	int max_borderSize = agentList->countActiveAgents * MAX_SUBCELLULAR_COMPONENTS;
	int borderSize = 0;
	//VoronoiCell *border[max_borderSize];
	VoronoiCell **border = (VoronoiCell**) malloc(sizeof(VoronoiCell*) * max_borderSize);
	//fprintf(stderr, "Call\n");
	for (int i = 0; i < agentList->countActiveAgents; i++) {
		//	fprintf( stderr, "agent: %i\n", agentList->agents[i]->index);
		if (agentList->agents[i]->state != VESSEL)
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++) {
				//		fprintf( stderr, "-->location: %i\n", agentList->agents[i]->location[ii]->index);

				char foundEmptyNeighbor = FALSE;
				// check for border
				if( agentList->agents[i]->location[ii]->isDomainBorder(vd)){
					if (borderSize == max_borderSize) {
						max_borderSize += agentList->countActiveAgents;
						border = (VoronoiCell**) realloc( border, max_borderSize * sizeof(VoronoiCell *));
					}
					//fprintf(stderr, "Found\n");
					border[borderSize++] = agentList->agents[i]->location[ii];

					foundEmptyNeighbor = TRUE;
				}

				// check for empty neighbors
				for (int iii = 0; !foundEmptyNeighbor && iii < agentList->agents[i]->location[ii]->countNeighborCells; iii++) {
					//			fprintf( stderr, "---->neighbor: %i\n", agentList->agents[i]->location[ii]->neighborCells[iii]->index);
					if ( GetAgent( agentList->agents[i]->location[ii]->neighborCells[iii]) == NULL) {
						if (borderSize == max_borderSize) {
							max_borderSize += agentList->countActiveAgents;
							border = (VoronoiCell**) realloc( border, max_borderSize * sizeof(VoronoiCell *));
						}

						border[borderSize++] = agentList->agents[i]->location[ii];
						foundEmptyNeighbor = TRUE;
						//				fprintf( stderr, "    >ADD TO QUEUE!\n");
						if (agentList->agents[i]->location[ii] == NULL) {
							fprintf(
									stderr,
									"agentList->agents[%i]->location[%i / %i] == NULL\n",
									i, ii,
									agentList->agents[i]->countLocations);
							exit(0);
						}
						if (agentList->agents[i]->location[ii]->agent == NULL) {
							fprintf(
									stderr,
									"agentList->agents[%i]->location[%i]->agent == NULL\n",
									i, ii);
							exit(0);
						}
						if (agentList->agents[i]->location[ii]->agent
								!= agentList->agents[i]) {
							fprintf(
									stderr,
									"agentList->agents[%i]->location[%i](%i)->agent(%i) != agentList->agents[%i](%i)\n",
									i,
									ii,
									agentList->agents[i]->location[ii]->index,
									GetAgent(agentList->agents[i]->location[ii])->index,
									i, agentList->agents[i]->index);
							exit(0);
						}
						/*if( GetAgent( border[ii]) == NULL){
						 fprintf( stderr, "border element %i == NULL => borderSize=%i\n", ii, borderSize);
						 exit(0);
						 }*/
					}
				}
			}
	}

	if (borderSize > max_borderSize) {
		fprintf(stderr, "border size to big: %i > %i\n", borderSize,
				max_borderSize);
		exit(0);
	}

	// center of mass
	//fprintf( stderr, "center of mass:");		
	double centerMass[DIMENSIONS];
	for (int i = 0; i < DIMENSIONS; i++) {
		centerMass[i] = 0.;
		for (int ii = 0; ii < borderSize; ii++) {
			if (GetAgent( border[ii]) == NULL) {
				fprintf(stderr, "border element %i == NULL => borderSize=%i\n",
						ii, borderSize);
				exit(0);
			}
			//fprintf( stderr, "border %i == (%lf,%lf,%lf)\n", ii, border[ii]->position[0], border[ii]->position[1], border[ii]->position[2]);
			centerMass[i] += border[ii]->position[i];
		}
		centerMass[i] /= (double) borderSize;
		//	fprintf( stderr, " %lf", centerMass[i]);
	}
	//fprintf( stderr, "\n");		

	// radius of gyration
	double gyrRadius = 0.;
	for (int ii = 0; ii < borderSize; ii++) {
		for (int i = 0; i < DIMENSIONS; i++) {
			gyrRadius += pow(border[ii]->position[i] - centerMass[i], 2.);
		}
	}
	gyrRadius /= (double) borderSize;
	//fprintf( stderr, "radius of gyration: %lf\n", gyrRadius);		

	free(border);

	return gyrRadius;
	return 0.;
}

Action *selectActionFixTimeStep(AgentList *agentArray, double &Time) {

	int i;
	// get oldest generation of dividing cells
	Agent *growingAgents[2 * (int) CellDivisionDepth];
	int countGrowingAgents = 0;
	int oldestGeneration = -1;
	int youngestGeneration = -1;
	for (i = 0; i < agentArray->countActiveAgents; i++) {
		if (agentArray->agents[i]->actions[INDEX_GROWTH]->top != NULL) {
			if (youngestGeneration < 0
					|| agentArray->agents[i]->generation > youngestGeneration) {
				youngestGeneration = agentArray->agents[i]->generation;
			}
			if (oldestGeneration < 0
					|| agentArray->agents[i]->generation < oldestGeneration) {
				oldestGeneration = agentArray->agents[i]->generation;
				countGrowingAgents = 0;
				//fprintf( stderr, "*");
			}
			//fprintf( stderr, "%i ", agentArray->agents[i]->generation);
			if (agentArray->agents[i]->generation == oldestGeneration)
				growingAgents[countGrowingAgents++] = agentArray->agents[i];
		}
	}
	//fprintf( stderr, " (%i x %i, %i)\n", countGrowingAgents, oldestGeneration, youngestGeneration);

	// chose random cell to divide
	int random = (int) (myRand() * countGrowingAgents);
	Action *selected_action = growingAgents[random]->actions[INDEX_GROWTH];

	// update time
	if (oldestGeneration == youngestGeneration)
		Time += 1. / selected_action->rate;

	return selected_action;
}

Action *selectActionVariablTimeStep(AgentList *agentArray, double &Time) {

	int i;
	// get oldest generation of dividing cells
	Agent *growingAgents[2 * (int) CellDivisionDepth];
	int countGrowingAgents = 0;
	int oldestGeneration = -1;
	int youngestGeneration = -1;
	for (i = 0; i < agentArray->countActiveAgents; i++) {
		if (agentArray->agents[i]->actions[INDEX_GROWTH]->top != NULL) {
			if (youngestGeneration < 0
					|| agentArray->agents[i]->generation > youngestGeneration) {
				youngestGeneration = agentArray->agents[i]->generation;
			}
			if (oldestGeneration < 0
					|| agentArray->agents[i]->generation < oldestGeneration) {
				oldestGeneration = agentArray->agents[i]->generation;
				//countGrowingAgents = 0;
				//fprintf( stderr, "*");
			}
			//fprintf( stderr, "%i ", agentArray->agents[i]->generation);
			//if( agentArray->agents[i]->generation==oldestGeneration)
			growingAgents[countGrowingAgents++] = agentArray->agents[i];
		}
	}
	//fprintf( stderr, " (%i x %i, %i)\n", countGrowingAgents, oldestGeneration, youngestGeneration);

	// sum rates
	double rateSum = 0.;
	for (i = 0; i < countGrowingAgents; i++)
		rateSum += growingAgents[i]->actions[INDEX_GROWTH]->rate;

	// chose random cell to divide
	double randomRateSum = myRand() * rateSum;
	double tempRateSum = 0.;
	//Action *selected_action = growingAgents[random]->actions[INDEX_GROWTH];
	for (i = 0; i < countGrowingAgents; i++) {
		tempRateSum += growingAgents[i]->actions[INDEX_GROWTH]->rate;
		if (tempRateSum >= randomRateSum) {

			// update time
			Time += 1. / rateSum;
			return growingAgents[i]->actions[INDEX_GROWTH];
		}
	}

	return NULL;
}



void statisticsOfSpheroid(VoronoiDiagram* voronoiDiagram,
		double** global_oxygen_concentration,
		double** global_glucose_concentration,
		double** global_active_cell_density,
		double** global_nonactive_cell_density,
		double** global_necrotic_cell_density, double radius,
		int radiusIntervals, int timeIndex,
		char gc_tmpout[], char ecm_tmpout[] ,char prolif_tmpout[]) {

	int i, k, kk;
	double radiusStep = radius / (double) radiusIntervals;
	//double timeStep =   time / (double) timeIntervals;

	/*double** global_oxygen_concentration = ( double** ) calloc ( EndTime_int, sizeof( double*));
	 for(k = 0;k < EndTime_int; k++){
	 global_oxygen_concentration[k] = ( double* ) calloc ( EndRadius_int, sizeof( double));
	 for(kk = 0;kk < EndRadius_int; kk++){
	 global_oxygen_concentration[k][kk] = 0.0;
	 }
	 }*/

	// center of mass
	double centerOfMass[DIMENSIONS];
	int countFoundCells = 0;
	for (k = 0; k < DIMENSIONS; k++)
		centerOfMass[k] = 0.0;

	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		if (voronoiDiagram->voronoiCells[i] != NULL) {
			//if( !voronoiDiagram->voronoiCells[i]->isFree()){
			for (k = 0; k < DIMENSIONS; k++)
				centerOfMass[k] += voronoiDiagram->voronoiCells[i]->position[k];
			countFoundCells++;
		}
	}
	if (countFoundCells)
		for (k = 0; k < DIMENSIONS; k++) {
			centerOfMass[k] /= (double) countFoundCells;
			//fprintf( stderr, "%lf\n", centerOfMass[k]);
		}
	else
		for (k = 0; k < DIMENSIONS; k++) {
			centerOfMass[k] =
					(voronoiDiagram->xMax[k] + voronoiDiagram->xMin[k]) / 2.;
			//fprintf( stderr, "%lf\n", centerOfMass[k]);
		}
	//exit( 0);

	// scan radial oxygen concentration
	//k = timeIndex;
	double distanceToMassCenter;
	int countCellsInInterval[radiusIntervals];
	for (i = 0; i < radiusIntervals; i++)
		countCellsInInterval[i] = 0;
	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		//if( GetAgent(voronoiDiagram->voronoiCells[i])->state != FREE){
		// determine distance to mass center
		distanceToMassCenter = 0.;
		for (k = 0; k < DIMENSIONS; k++)
			distanceToMassCenter +=
					pow(
							centerOfMass[k]
									- voronoiDiagram->voronoiCells[i]->position[k],
							2);
		distanceToMassCenter = sqrt(distanceToMassCenter);

		//
		kk = (int) (distanceToMassCenter / radiusStep);
		if (kk < radiusIntervals)
			countCellsInInterval[kk]++;
		//}
	}

	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		//if( GetAgent(voronoiDiagram->voronoiCells[i])->state != FREE){
		// determine distance to mass center
		distanceToMassCenter = 0.;
		for (k = 0; k < DIMENSIONS; k++)
			distanceToMassCenter +=
					pow(
							centerOfMass[k]
									- voronoiDiagram->voronoiCells[i]->position[k],
							2.);
		distanceToMassCenter = sqrt(distanceToMassCenter);

		//
		kk = (int) (distanceToMassCenter / radiusStep);
		if (kk >= radiusIntervals) {
			fprintf(stderr, "Index too big: %i (< %i)!!! <= %lf\n", kk,
					radiusIntervals, distanceToMassCenter);

			exit(0);
		}
		if (kk < radiusIntervals) {
			global_oxygen_concentration[timeIndex][kk] +=
					voronoiDiagram->voronoiCells[i]->oxygen
							/ (double) countCellsInInterval[kk];
			global_glucose_concentration[timeIndex][kk] +=
					voronoiDiagram->voronoiCells[i]->glucose
							/ (double) countCellsInInterval[kk];
		}
		if (voronoiDiagram->voronoiCells[i]->getState() == ACTIVE)
			global_active_cell_density[timeIndex][kk] +=
					1. / (double) countCellsInInterval[kk];
		if (voronoiDiagram->voronoiCells[i]->getState() == NONACTIVE)
			global_nonactive_cell_density[timeIndex][kk] +=
					1. / (double) countCellsInInterval[kk];
		if (voronoiDiagram->voronoiCells[i]->getState() == NECROTIC)
			global_necrotic_cell_density[timeIndex][kk] +=
					1. / (double) countCellsInInterval[kk];

		//}
	}
}

double InitialOxygenConcentration = 1.;
double InitialGlucoseConcentration = 1.;

double montecarlo(double InitialRadius_in, double InitialQuiescentFraction_in, double MaxCellDivisionRate_in, double DivisionDepth_in,
					double ECMThresholdQuiescence_in, double ECMProductionRate_in, double ECMDegradationRate_in,
					double EndTimeIn, double OutputRate, double profileTime, int profileDepth, int rand_seed,
					std::vector<double> &gc_out, std::vector<double> &ecm_out, std::vector<double> &prolif_out)
{
    double InitialRadius = InitialRadius_in;
    double InitialQuiescentFraction = InitialQuiescentFraction_in;
    MaxCellDivisionRate = MaxCellDivisionRate_in;
    double DivisionDepth = DivisionDepth_in;
	double ECMThresholdQuiescence = ECMThresholdQuiescence_in;
	double ECMProductionRate = ECMProductionRate_in;
	double ECMDegradationRate = ECMDegradationRate_in;

    ReentranceProbabilityReferenceLength = DivisionDepth;

	// Time
	double Time, Last_Time, EndEndTime, EndTime;
	int AdaptEndTime = FALSE;
	EndEndTime = EndTime = EndTimeIn;
	
	// Variables introduced by dennis
	bool writeFullSimulation = false;
	bool writeRawComparison = false;
	bool writeRawSimulation = true;
	
	// Related to potential outputs
	int nGrowthCurveTime;
	double* growthCurveOut_time;
	double* growthCurveOut_value;
	int nProfileDepth;
	double* ecm_profile;
	double* prolif_profile;
	
	//OutputRate = atof(optarg);
	nGrowthCurveTime = (int)(EndEndTime / OutputRate);	
	//fprintf(stdout,"\n N Time points: %i\n",nGrowthCurveTime );
	
	// Variables introduced by dennis that will later be provided by the user
	
	int histogram[10000];
	int histogramDividing[10000];
	int histogramNecrotic[10000];
	int histogramFree[10000];
	double histogramECM[10000];

	// AVERAGES & STD DERIVATION
	int histogramSquaresDividing[10000];
	int histogramSquaresNecrotic[10000];
	int histogramSquaresDividingFraction[10000];
	int histogramSquaresNecroticFraction[10000];
	int histogramSquaresFree[10000];
	double histogramSquaresECM[10000];

	int histogramSumDividing[10000];
	int histogramSumNecrotic[10000];
	int histogramSumDividingFraction[10000];
	int histogramSumNecroticFraction[10000];
	int histogramSumFree[10000];
	double histogramSumECM[10000];
	// -

	for (int i = 0; i < 10000; i++) {
		histogram[i] = 0;

		histogramDividing[i] = 0;
		histogramNecrotic[i] = 0;
		histogramFree[i] = 0;
		histogramECM[i] = 0.;

		histogramSquaresDividing[i] = 0;
		histogramSquaresNecrotic[i] = 0;
		histogramSquaresDividingFraction[i] = 0;
		histogramSquaresNecroticFraction[i] = 0;
		histogramSquaresFree[i] = 0;
		histogramSquaresECM[i] = 0.;

		histogramSumDividing[i] = 0;
		histogramSumNecrotic[i] = 0;
		histogramSumDividingFraction[i] = 0;
		histogramSumNecroticFraction[i] = 0;
		histogramSumFree[i] = 0;
		histogramSumECM[i] = 0.;
	}

	int i, k, l;
	//int j;
	unsigned short seed_start = 0;
	int c;
	//Grid *voronoiGrid;    // Grid structure
#if( DETERMINE_COARSE_GRAINED_BEHAVIOR)
	Grid *coarseVoronoiGrid; // Grid structure for the coarse grained model
#endif
	ActionTree *actionList; // probability list structure
	Action *selected_action;

	// File and Directory Variables
	char parameter_file[FILENAMESIZE];
	char dirname[FILENAMESIZE];
	char outfilename[FILENAMESIZE];
	char parameterfilename[FILENAMESIZE];
	//char init_filename[ FILENAMESIZE ]; 

	char filename_fineGrid[FILENAMESIZE] = "";
	char filename_coarseGrid[FILENAMESIZE];

	char pnm_suffix[FILENAMESIZE + 10], pnmstring[FILENAMESIZE];
	parameter_file[0] = '\0';

	// Filepointer
	FILE *fp;
	FILE *fp_param;
	char boundaryCondition[FILENAMESIZE] = "";

	// Statistic Variables
	int count_inoculated_cells = 0;
	double count_divided_cells = 0;
	int count_expanded_cells = 0;
	int count_cells = 0;
	int count_cell_volume = 0;
	int count_vessel_cells = 0;
	int count_necrotic_cells = 0;
	int count_necrotic_volume = 0;
	int count_inner_empty_volume = 0;
	int count_lysed_cells = 0;
	int count_lysed_volume = 0;
	double gyrRadius = 0.;
	double prob_X_equal_k = 0.;
	double prob_X_smaller_k = 0.;
	double prob_X_bigger_k = 0.;

	Statistix myStatistics;

	myStatistics.count_cells = &count_cells;
	myStatistics.count_cell_volume = &count_cell_volume;
	myStatistics.count_expanded_cells = &count_expanded_cells;
	myStatistics.count_divided_cells = &count_divided_cells;
	myStatistics.count_vessel_cells = &count_vessel_cells;
	myStatistics.count_necrotic_cells = &count_necrotic_cells;
	myStatistics.count_necrotic_volume = &count_necrotic_volume;
	myStatistics.count_inner_empty_volume = &count_inner_empty_volume;
	myStatistics.count_lysed_cells = &count_lysed_cells;
	myStatistics.count_lysed_volume = &count_lysed_volume;


	int Last_NumberOfCells = 0;
	double Last_gyrRadius = 0;


	// Scaling Variables
	int EndTime_int = 0;
	//double OutputRate = 1.;

	// Averaging Arrays
	double *global_cells;
	double *global_volume;
	double *global_inner_empty_volume;
	double *global_necrotic_cells;
	double *global_necrotic_volume;
	double *global_vessel_cells;
	double *global_prob_X_equal_k;
	double *global_prob_X_smaller_k;
	double *global_prob_X_bigger_k;

	// Time & Timer Variables
	time_t timer;
	time_t temp_t;
	char local_time[256];
	struct tm *loctime = (struct tm *) malloc(sizeof(struct tm));
	long passedTime;
	long benchmarkTime = 0;

	int resolution = 10;

	// VESSEL STUFF
	double DistanceToInitialCell = 10.0;
	int CountVessel = 1;
	double BranchingProbability = 0.;
	int BranchingLength = 20;

	int Averages = 0;
	int NumberOfCellsPerSite = 2;

	bool NoRadialProfiles = true;
	bool NoSliceOutput = true;
	bool PovrayOutput = false;

	int     RadialProfilesCount = 1;
	double *RadialProfilesTime  =  (double*)malloc(sizeof(double));
	RadialProfilesTime[0] = profileTime;
	

	double WilliamsCodeParameter = 0.;

	// OUTPUTS
	int OutputAnimation = FALSE;
	char CustomDirectoryName = FALSE;
	int rotations = 0;
	double setoffX = 30, setoffY = 30, setoffZ = 30;
	double startAngle = -0.5 * PI;

	double customTimestep = 0.1;
	double customDiffusionCoefficientFactor = 1.;

	double custom_Waste_Diffusion = Waste_Diffusion;

	//double InitialRadius = 1;
	//double InitialQuiescentFraction = 0.;

	// DATA
	comparison_t data_growthcurve = create_comparison();
	comparison_t data_KI67 = create_comparison();
	comparison_t data_ECM  = create_comparison();

	double maxRadius  = 0;
	double maxEpsilon = DBL_MAX;
	double cumEpsilon = 0;
	int    idxEpsilon = 0;
	double measurement_error = 0;
	double measurement_error_gc   = 0;
	double measurement_error_ki67 = 0;
	double measurement_error_ecm  = 0;


	// PARAMETERS SET TO DEFAULT VALUES BY DENNIS
	ReentranceProbabilityFunction = EXPONENTIAL;
	M_gro = 9;
	CellDivisionDepth = 10.0;
	

	// Parameters set to default values by Dennis:
	Averages = 1;

	CustomDirectoryName = TRUE;
	seed_start =  rand_seed;
	
	// Parameters used for fitting
	VoronoiCell::ECM_THRESHOLD_QUIESCENCE = ECMThresholdQuiescence;
	VoronoiCell::ECM_PRODUCTION_RATE = ECMProductionRate;
	VoronoiCell::ECM_DEGRADATION_RATE = ECMDegradationRate;
	
	srand(0);
	//customDiffusionCoefficientFactor /= pow( (double)CountCellsPerVoronoiCell,2./3.);
	Oxygen_Diffusion /= pow((double) CountCellsPerVoronoiCell, 2. / DIMENSIONS);
	Glucose_Diffusion /= pow((double) CountCellsPerVoronoiCell, 2. / DIMENSIONS);


	// ADDITIVE MEASUREMENT ERROR
	measurement_error_gc = 0;
	for( int j=0; j<data_growthcurve.dim; j++)
		measurement_error_gc = max<double>( measurement_error_gc, data_growthcurve.m[j]);
	measurement_error_gc *= measurement_error;

	measurement_error_ki67 = 0;
	for( int j=0; j<data_KI67.dim; j++)
		measurement_error_ki67 = max<double>( measurement_error_ki67, data_KI67.m[j]);
	measurement_error_ki67 *= measurement_error;

	measurement_error_ecm = 0;
	for( int j=0; j<data_ECM.dim; j++)
		measurement_error_ecm = max<double>( measurement_error_ecm, data_ECM.m[j]);
	measurement_error_ecm *= measurement_error;


	// FILE AND DIRECTORY NAMES
	if (false && !CustomDirectoryName) {
		sprintf(dirname, "C%iD%iM%im%ik%in%i", Case, DIMENSIONS, M_gro + 1,
				M_div + 1, (int) CellDivisionDepth, CountCellsPerVoronoiCell);
		mkdir(dirname MODUS);
	}



	VoronoiDiagram * voronoiDiagram;
	VoronoiDiagram * voronoiDiagramTmp;
	int countPoints[DIMENSIONS];

	if( strlen(filename_fineGrid)){

		//fprintf(stderr, "Reading Voronoi Lattice from File... \n");
		voronoiDiagram= new VoronoiDiagram(filename_fineGrid);


		for (i = 0; i < DIMENSIONS; i++)
			countPoints[i] = (int) (pow(voronoiDiagram->countVoronoiCells, 1. / DIMENSIONS)	+ 0.5);
	}else{
		for (i = 0; i < DIMENSIONS; i++)
			countPoints[i] = 100;

		char latticeFilename[512], testFilename[512];
		int  latticeSize = 100;
		sprintf(latticeFilename, "staticdata/%ipow%i", latticeSize, DIMENSIONS);


		voronoiDiagramTmp = VoronoiDiagram::newVoronoiDiagram( latticeSize, latticeSize, (DIMENSIONS == 2 ? 1 : latticeSize));
		voronoiDiagram = new VoronoiDiagram(latticeFilename, voronoiDiagramTmp);
	}
	voronoiDiagram->boundaryThickness = 1.;

	
	voronoiDiagram->setDomain();
	AgentList* agentArray = AgentList::newAgentList(
			voronoiDiagram->countVoronoiCells);

	// set extended neighborhood
	VoronoiCell::extendedNeighborDepth = CellDivisionDepth;
	if (!(VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD
			|| VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD)
			&& CellDivisionDepth > 1) {

		// VORONOI DIAGRAM
		if (!voronoiDiagram->readExtendedNeighborhoodToFile(filename_fineGrid,
				(int) CellDivisionDepth)) {
			//fprintf( stderr, "Couldn't read Extended Neighborhood from file -> Calculate my self!\n");
			double innerDomainRadius = voronoiDiagram->xMax[0]
					- voronoiDiagram->xMin[0];
			for (i = 1; i < DIMENSIONS; i++) {
				if (innerDomainRadius
						> voronoiDiagram->xMax[i] - voronoiDiagram->xMin[i])
					innerDomainRadius =
							voronoiDiagram->xMax[i] - voronoiDiagram->xMin[i];
			}
			innerDomainRadius /= 2.;
			fprintf(
					stderr,
					"innerDomainRadius = %lf, CellDivisionDepth = %lf, CentralCell = %p\n",
					innerDomainRadius, CellDivisionDepth,
					getCentralCell(voronoiDiagram));
			voronoiDiagram->NEWsetExtendedNeighborhood((int) CellDivisionDepth);

		}


	}

	// initialize diffusion/kinetics

	Substrate substrate(countPoints, Case, voronoiDiagram, agentArray,
			InitialOxygenConcentration, InitialGlucoseConcentration,
			WilliamsCodeParameter, dirname, boundaryCondition, customTimestep);

	DiffusionReactionEquation *lactateDynamics = 0;
	float *lactate = 0;
	float *lactateProduction = 0;
	float *lactateDiffusion = 0;

	if (VoronoiCell::USE_LACTATE) {
		fprintf(stderr, "Initialize Lactate Dynamics\n");
		lactateDynamics =
				new DiffusionReactionEquation(voronoiDiagram,
						DiffusionReactionEquation::STEADYSTATE, //time
						DiffusionReactionEquation::EXPLICIT, //reaction
						DiffusionReactionEquation::IMPLICIT //diffusion
						);
		lactateDynamics->setSpaceStep(SPATIAL_UNIT);
		lactate =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		lactateProduction =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		lactateDiffusion =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);

		lactateDynamics->initSystem(lactate);
		lactateDynamics->initReaction(lactateProduction);
		lactateDynamics->initDiffusion(lactateDiffusion);

		for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			lactate[i] = 0.;
			lactateDiffusion[i] = Lactate_Diffusion;
		}
	}

	DiffusionReactionEquation *wasteDynamics = 0;
	float *waste = 0;
	float *wasteProduction = 0;
	float *wasteUptake = 0;
	float *wasteDiffusion = 0;

	if (VoronoiCell::USE_WASTE) {
		fprintf(stderr, "Initialize Waste Dynamics\n");
		wasteDynamics =
				new DiffusionReactionEquation(voronoiDiagram,
						DiffusionReactionEquation::IMPLICIT, //time
						DiffusionReactionEquation::IMPLICIT, //reaction
						DiffusionReactionEquation::IMPLICIT //diffusion
						);
		wasteDynamics->setSpaceStep(SPATIAL_UNIT);
		waste =	(float*) malloc(sizeof(float) * voronoiDiagram->countVoronoiCells);
		wasteProduction = (float*) malloc(sizeof(float) * voronoiDiagram->countVoronoiCells);
		wasteUptake = (float*) malloc(sizeof(float) * voronoiDiagram->countVoronoiCells);
		wasteDiffusion = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);

		wasteDynamics->initSystem(waste);
		wasteDynamics->initReaction(wasteProduction);
		wasteDynamics->initReactionDerivative(wasteUptake);
		wasteDynamics->initDiffusion(wasteDiffusion);

		for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			waste[i] = 0.;
			wasteDiffusion[i] = custom_Waste_Diffusion;
		}
	}

	DiffusionReactionEquation *morphogenDynamics = 0;
	float *morphogen = 0;
	float *morphogenProduction = 0;
	float *morphogenDecay = 0;
	float *morphogenDiffusion = 0;

	if (VoronoiCell::USE_MORPHOGEN) {
		fprintf(stderr, "Initialize Morphogen Dynamics\n");
		morphogenDynamics =
				new DiffusionReactionEquation(voronoiDiagram,
						DiffusionReactionEquation::IMPLICIT, //time
						DiffusionReactionEquation::IMPLICIT,//EXPLICIT, //reaction
						DiffusionReactionEquation::IMPLICIT //diffusion
						);
		morphogen =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		morphogenProduction =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		morphogenDecay =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		morphogenDiffusion =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);

		morphogenDynamics->initSystem(morphogen);
		morphogenDynamics->initReaction(morphogenProduction);
		morphogenDynamics->initReactionDerivative(morphogenDecay);
		morphogenDynamics->initDiffusion(morphogenDiffusion);
		morphogenDynamics->setTimeStep(customTimestep);

		for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			morphogen[i] = 0.;
			morphogenDiffusion[i] = Lactate_Diffusion;
			morphogenProduction[i] = 0.;//-300;
		}
	}

	double *dECM = 0;
	if (VoronoiCell::ECM_THRESHOLD_QUIESCENCE != 0) {

		dECM = (double*)malloc( voronoiDiagram->countVoronoiCells*sizeof(double));
		for( int v=0; v<voronoiDiagram->countVoronoiCells; v++)
			dECM[v]=0;
	}




	// INITIALIZE STATISTICAL ARRAYS

	// calculate number of discrete measure points 
	EndTime_int = nrOfSteps( BeginningTime, EndTime, OutputRate);

	// allocate memory for arrays
	global_cells = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_cells[k] = 0.0;

	global_necrotic_cells = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_necrotic_cells[k] = 0.0;

	global_necrotic_volume = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_necrotic_volume[k] = 0.0;

	global_vessel_cells = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_vessel_cells[k] = 0.0;

	global_volume = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_volume[k] = 0.0;

	global_inner_empty_volume =
			(double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_inner_empty_volume[k] = 0.0;

	double *global_gyrRadius = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_gyrRadius[k] = 0.0;

	double *global_gyrRadiusSquare = (double*) calloc(EndTime_int + 1,
			sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_gyrRadiusSquare[k] = 0.0;

	double *global_avgOxygen = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_avgOxygen[k] = 0.0;

	double *global_avgGlucose = (double*) calloc(EndTime_int + 1,
			sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_avgGlucose[k] = 0.0;


	// WRITE ALL PARAMETERS TO FILE
	sprintf(parameterfilename, "%s%cparameter.dat", dirname, SEPERATOR);
	if(writeFullSimulation){
		if ((fp_param = fopen(parameterfilename, "w")) == NULL) {
			fprintf(stderr, "Error opening file %s for writing!\n",
					parameterfilename);
		}
		fprintf(fp_param, "Command Line:\n");
		
		//fprintf(fp_param, "\n");

		fprintf(fp_param, "******************************** \n");
		fprintf(fp_param, "\nDefined Parameters: \n");
		fprintf(fp_param, "#define MIN_SPECIFIC_DEATH_RATE %lf\n",
				MIN_SPECIFIC_DEATH_RATE);
		fprintf(fp_param, "#define MAX_SPECIFIC_DEATH_RATE %lf\n",
				MAX_SPECIFIC_DEATH_RATE);
		fprintf(fp_param, "#define MONOD_PARAMETER_DEATH %lf\n",
				MONOD_PARAMETER_DEATH);
		fprintf(fp_param, "#define MAX_SPECIFIC_GROWTH_RATE %lf\n",
				MAX_SPECIFIC_GROWTH_RATE);
		fprintf(fp_param, "#define MONOD_PARAMETER_GROWTH %lf\n",
				MONOD_PARAMETER_GROWTH);
		fprintf(fp_param, "#define CELL_GROWTH_YIELD_OXYGEN %lf\n",
		CELL_GROWTH_YIELD_OXYGEN);
		fprintf(fp_param, "#define CELL_GROWTH_YIELD_GLUCOSE %lf\n",
		CELL_GROWTH_YIELD_GLUCOSE);
		fprintf(fp_param, "#define MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN %lf\n",
		MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN);
		fprintf(fp_param, "#define MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE %lf\n",
		MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE);
		fprintf(fp_param, "#define MONOD_KINETICS %i\n", MONOD_KINETICS);
		fprintf(fp_param, "#define USE_DIVISION %i\n", USE_DIVISION);
		fprintf(fp_param, "#define USE_MIGRATION %i\n", USE_MIGRATION);
		fprintf(fp_param, "#define USE_APOPTOSIS %i\n", USE_APOPTOSIS);
		fprintf(fp_param, "#define USE_NECROSIS %i\n", USE_NECROSIS);
		fprintf(fp_param, "#define USE_LYSIS %i\n", USE_LYSIS);
		fprintf(fp_param, "#define MONOD_KINETICS %i\n", MONOD_KINETICS);
		fprintf(fp_param, "#define MULTISCALE %i\n", MULTISCALE);
		fprintf(fp_param, "#define AGENT_DIAMETER %lf\n", AGENT_DIAMETER);
		fprintf(fp_param, "#define SPATIAL_UNIT %lf\n", SPATIAL_UNIT);

		fclose(fp_param);
	}
	/** STOCHASTIC APPROACHE *************************************************/

	double timeActualization = 0, timeSelect = 0, timeExecution = 0;

#define USE_SPARSE_MATRIX

	double timeStep = customTimestep; //0.1;
	double timeDifference = 0.;
	if (timeStep > timeDifference)
		timeDifference = timeStep;

	int N = 2 * voronoiDiagram->xN[0]
#if DIMENSIONS >= 2
			* voronoiDiagram->xN[1]
#endif
#if DIMENSIONS >= 3
	* voronoiDiagram->xN[2]
#endif
	;
	if (Case > 1) {
		sA = SparseMatrix::newSparseMatrix(N, N);
		b = (float *) malloc(N * sizeof(float));
		x = (float *) malloc(N * sizeof(float));
		v0 = (float *) malloc(N * sizeof(float));
		v1 = (float *) malloc(N * sizeof(float));
		v2 = (float *) malloc(N * sizeof(float));
		v3 = (float *) malloc(N * sizeof(float));
		v4 = (float *) malloc(N * sizeof(float));
		v5 = (float *) malloc(N * sizeof(float));
		v6 = (float *) malloc(N * sizeof(float));
		v7 = (float *) malloc(N * sizeof(float));
		v8 = (float *) malloc(N * sizeof(float));
	}

	// Interpolation
	Interpolation::voronoiDiagram = voronoiDiagram;

	for (k = 0; k < Averages; k++) {
		// INITIALIZE RANDOM GENERATOR
		srand(seed_start + k);

		Last_Time = Time = TempTime = 0;
		gyrRadius = 0;
		count_inoculated_cells = 0;
		count_divided_cells = 0;
		count_expanded_cells = 0;
		count_cells = 0;
		count_cell_volume = 0;
		count_vessel_cells = 0;
		count_necrotic_cells = 0;
		count_necrotic_volume = 0;
		count_inner_empty_volume = 0;
		count_lysed_cells = 0;
		count_lysed_volume = 0;
		prob_X_equal_k = 0;
		prob_X_smaller_k = 0;
		prob_X_bigger_k = 0;
		double global_oxygen_concentration = InitialOxygenConcentration;
		double global_glucose_concentration = InitialGlucoseConcentration;

		timer = time(NULL);


		// INITIALIZE PROBABILITY LIST
		actionList = ActionTree::newActionTree();

		// SET INITIAL VESSEL NETWORK
		VoronoiCell *centralCell = getCentralCell(voronoiDiagram);
		centralCell =
				centralCell->neighborCells[(int) (centralCell->countNeighborCells
						* myRand())];


		if( VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD){
		// Init Symbolic Extended Neighborhood
			// extended neighborship
			VoronoiCell ***extendedNeighborhood;
			int extendedNeighborhoodSize[(int) CellDivisionDepth + 1];
			int extendedNeighborhoodMax[(int) CellDivisionDepth + 1];
			extendedNeighborhood =
					(VoronoiCell ***) malloc(
							(int) (CellDivisionDepth + 1)
									* sizeof(VoronoiCell **));

			// 1st gen
			extendedNeighborhood[0] =
					(VoronoiCell **) malloc(1 * sizeof(VoronoiCell *));
			extendedNeighborhood[0][0] = centralCell;
			extendedNeighborhoodSize[0] = 1;

			// 2nd gen
			extendedNeighborhood[1] =
					(VoronoiCell **) malloc(
							centralCell->countNeighborCells
									* sizeof(VoronoiCell *));
			for (int n = 0; n < centralCell->countNeighborCells; n++)
				extendedNeighborhood[1][n] = centralCell->neighborCells[n];
			extendedNeighborhoodSize[1] = centralCell->countNeighborCells;

			// following gen
			for (int g = 2; g <= (int) CellDivisionDepth; g++) {
				// alloc memory
				extendedNeighborhood[g] =
						(VoronoiCell **) malloc(1000 * sizeof(VoronoiCell *));
				extendedNeighborhoodSize[g] = 0;
				extendedNeighborhoodMax[g] = 1000;

				//fprintf( stderr, "gen:%i\n", g);

				// look for all points in old generation
				for (int i = 0; i < extendedNeighborhoodSize[g - 1]; i++)
					// look at surrounding points
					for (int n = 0;
							n
									< extendedNeighborhood[g - 1][i]->countNeighborCells;
							n++) {
						// if not added already
						//fprintf( stderr, "TEST cell:%i\n", extendedNeighborhood[g-1][i]->neighborCells[n]->index);
						if (!isElementOf(
								extendedNeighborhood[g - 2],
								extendedNeighborhoodSize[g - 2],
								extendedNeighborhood[g - 1][i]->neighborCells[n])
								&& !isElementOf(
										extendedNeighborhood[g - 1],
										extendedNeighborhoodSize[g - 1],
										extendedNeighborhood[g - 1][i]->neighborCells[n])
								&& !isElementOf(
										extendedNeighborhood[g],
										extendedNeighborhoodSize[g],
										extendedNeighborhood[g - 1][i]->neighborCells[n])) {
							// realloc
							if (extendedNeighborhoodSize[g]
									== extendedNeighborhoodMax[g]) {
								extendedNeighborhoodMax[g] += 1000;
								extendedNeighborhood[g] =
										(VoronoiCell **) realloc(
												extendedNeighborhood[g],
												extendedNeighborhoodMax[g]
														* sizeof(VoronoiCell *));
							}
							// add
							extendedNeighborhood[g][extendedNeighborhoodSize[g]] =
									extendedNeighborhood[g - 1][i]->neighborCells[n];
							extendedNeighborhoodSize[g]++;
						}
					}
			}

			// symbolic extended (with direct neighbors)
			VoronoiCell::symbolicExtendedNeighborhood =
					(int**) malloc((int) (CellDivisionDepth) * sizeof(int *));
			VoronoiCell::symbolicExtendedNeighborhoodSize =
					(int*) malloc((int) (CellDivisionDepth) * sizeof(int));
			for (int i = 0; i < (int) CellDivisionDepth; i++) {
				//fprintf(stderr, "generation %i: ", i);
				VoronoiCell::symbolicExtendedNeighborhoodSize[i] =
						extendedNeighborhoodSize[i + 1];
				VoronoiCell::symbolicExtendedNeighborhood[i] =
						(int*) malloc(
								VoronoiCell::symbolicExtendedNeighborhoodSize[i]
										* sizeof(int));
				for (int j = 0;
						j < VoronoiCell::symbolicExtendedNeighborhoodSize[i];
						j++) {
					/*fprintf(
							stderr,
							" %i,",
							extendedNeighborhood[i + 1][j]->index
									- centralCell->index);*/
					VoronoiCell::symbolicExtendedNeighborhood[i][j] =
							extendedNeighborhood[i + 1][j]->index
									- centralCell->index;
				}
				//fprintf(stderr, "\n");
				free(extendedNeighborhood[i + 1]);
			}
			free(extendedNeighborhood[0]);
			//free( extendedNeighborhood[1]);
			free(extendedNeighborhood);
		}

		double DistanceToInitialCellArray[DIMENSIONS];
		//DistanceToInitialCellArray[0] = DistanceToInitialCell;
		//for( i=1; i<DIMENSIONS; i++)
		for (i = 0; i < DIMENSIONS; i++)
			DistanceToInitialCellArray[i] = DistanceToInitialCell; // default value

		switch (Case) {
		case 1:
		case 2:
			break;
		case 3:
			break;
		case 4:
		case 5:
			fprintf(stderr, "Initialize Vessel Network...\n");
			//@ author Emmanuel
			#if(DIMENSIONS == 3)
			SetRegularInitialVesselNetwork(voronoiDiagram, agentArray,
					actionList, centralCell, DistanceToInitialCellArray,
					BranchingProbability, BranchingLength);
			fprintf(stderr, "...finished\n");
			#endif
			break;
		}


		substrate.Initialization(Case);

		for (int f = 0; f < voronoiDiagram->countFramePoints; f++) {
			agentArray->activateAgent()->attach(voronoiDiagram->framePoints[f]);
			GetAgent(voronoiDiagram->framePoints[f])->state = FRAME;
		}


		//fprintf(stderr, "Init cells\n");
		switch( 0){
		case 0:{
			// START WITH BALL OF CELLS
			//double InitialRadius = 10;
			for (int v = 0; v < voronoiDiagram->countVoronoiCells; v++) {
				double dist = centralCell->getDistanceTo(
						voronoiDiagram->voronoiCells[v]);
				if (dist < InitialRadius) {

					// AGENT
					agentArray->activateAgent()->attach(
							voronoiDiagram->voronoiCells[v]);
					GetAgent( voronoiDiagram->voronoiCells[v])->growingTumorCellCount = 1;
					GetAgent( voronoiDiagram->voronoiCells[v])->dividingTumorCellCount = 0;
					GetAgent( voronoiDiagram->voronoiCells[v])->cellCount = 1;
					if( CountCellsPerVoronoiCell == 1){
						GetAgent( voronoiDiagram->voronoiCells[v])->state = ACTIVE; //COMPARTMENT;
						update_surrounding_added_cell(actionList, voronoiDiagram->voronoiCells[v], voronoiDiagram);
					}else{
						GetAgent( voronoiDiagram->voronoiCells[v])->state = COMPARTMENT;
						GetAgent( voronoiDiagram->voronoiCells[v])->maxCellCount = CountCellsPerVoronoiCell;
						GetAgent( voronoiDiagram->voronoiCells[v])->actualize( actionList);
					}


					// STATS
					count_inoculated_cells++;
					count_cells++;
					count_cell_volume++;
				}
			}
		}	break;

		case 1:{
			// START with a DENSITY OF CELLS in a SPHERE
			double center[3]={50,50,50};
			for (int v = 0; v < voronoiDiagram->countVoronoiCells; v++)
			if(voronoiDiagram->voronoiCells[v]->getDistanceTo(center)<50.){
				if (0.1 > myRand()) {

					// AGENT
					agentArray->activateAgent()->attach(
							voronoiDiagram->voronoiCells[v]);
					GetAgent( voronoiDiagram->voronoiCells[v])->growingTumorCellCount =
							1;
					GetAgent( voronoiDiagram->voronoiCells[v])->dividingTumorCellCount =
							0;
					GetAgent( voronoiDiagram->voronoiCells[v])->cellCount = 1;
					GetAgent( voronoiDiagram->voronoiCells[v])->state = ACTIVE; //COMPARTMENT;
					update_surrounding_added_cell(actionList,
							voronoiDiagram->voronoiCells[v], voronoiDiagram);

					// STATS
					count_inoculated_cells++;
					count_cells++;
					count_cell_volume++;
				}
			}
		}	break;
		}

		// SET 3/4 of the CELLS QUIESCENT
		for (int a = 0; a < agentArray->countActiveAgents; a++) {
			agentArray->agents[a]->divide = //1;
					(myRand() < (1.-InitialQuiescentFraction) * PROB_REENTERING_CELL_CYCLE( SPATIAL_UNIT * GetDistanceToClosestFreeNeighbor( voronoiDiagram, agentArray->agents[a]->location[0])) ? 1: 0);
		}


		// INITIALIZE STATISTICAL VARIABLES
		timer = time(NULL);
		Last_Time = Time = 0.0;
		//gyrRadius = getGyrationRadius(agentArray);
		Last_gyrRadius = gyrRadius = getGyrationRadiusOfBorder(agentArray, voronoiDiagram);

		l = (int) Time;
		global_cells[l] += count_cells;
		global_volume[l] += count_cell_volume;
		global_necrotic_cells[l] += count_necrotic_cells;
		global_necrotic_volume[l] += count_necrotic_volume;
		global_inner_empty_volume[l] += count_inner_empty_volume;
		global_vessel_cells[l] += count_vessel_cells;
		global_gyrRadius[l] += sqrt(gyrRadius);
		global_gyrRadiusSquare[l] += gyrRadius;
		global_avgGlucose[l] += getAvgGlucose(agentArray);
		global_avgOxygen[l] += getAvgOxygen(agentArray);


		actionList->addAction( newAction( GAP, 1./OutputRate));
		//actionList->addAction( newAction( GAP, 1./TIME_STEP));

		/** SIMULATION **/

		for (i = 0; i < voronoiDiagram->countVoronoiCells; i++)
			voronoiDiagram->voronoiCells[i]->ecm = 0.;

		Last_Time = -1.;

		int flag_actualizeProcessRates = TRUE;
		FILE *fp_dbg;
		do {
			//fp_dbg = fopen( "outdir/dbg.txt", "a+");	fprintf(fp_dbg ," Start of time loop \n");	fclose(fp_dbg );
			//printf( "count_cells: %i == %i (active agents)\n", count_cells, agentArray->countActiveAgents);
			// INITIALIZE STATISTICAL VARIABLES
			Last_Time = Time;
			Last_NumberOfCells =
					count_inoculated_cells + (int) count_divided_cells
							- count_necrotic_cells - count_lysed_cells;

			if (flag_actualizeProcessRates == TRUE) {
				//fprintf( stderr, "actualizeProcessRates\n");
				actionList->actualizeAllRates(voronoiDiagram);
				flag_actualizeProcessRates = FALSE;
				//actionList->getDepth( actionList->root);
			}
			//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Select action \n");fclose(fp_dbg );
			
			// SELECT ACTION
			if (actionList->size != 0 && actionList->rateSum != 0.) {

				bool nullEvent = false;

				//actionList->getDepth( actionList->root);
				//fprintf(stderr, "actualizeProcessRates();\n");
				passedTime = clock();
				timeActualization +=
						(double) (clock() - passedTime)
								/ (double) CLOCKS_PER_SEC;

				//printf("selectAction();\n");
				passedTime = clock();
				selected_action = actionList->selectAction(&Time);
				timeSelect +=
						(double) (clock() - passedTime)
								/ (double) CLOCKS_PER_SEC;

				//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Execute action \n");fclose(fp_dbg );
				// EXECUTE ACTION
				//printf("EXECUTE ACTION\n");
				if (selected_action->rate != 0 && Time <= EndTime) {
					passedTime = clock();
					//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Starting case switch....\n");fclose(fp_dbg );fp_dbg = fopen( "outdir/dbg.txt", "a+");
					switch (selected_action->type) {

					// PROCESSES ON THE MICRO CARRIER SURFACE

					case MIGRATION:
						//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," MIGRATION\n");
						//fprintf(stderr, "MIGRATION\n");
						if (0.5 < myRand()
								&& GetVoronoiCell(selected_action->originalCell)->isDomainBorder(
										voronoiDiagram)) {
							selected_action->originalCell->state = FREE;
							selected_action->originalCell->actualize(
									actionList);
							selected_action->originalCell->detach();
							fprintf(stderr,
									"WARNING: Agent migrates out of domain!\n");
						} else{


							if( Agent::USE_GRAVITY && Time < 0.2/3600 /* 0.1 sec */){
							//performChemotacticMigration
								//performFreeMigration
								performMigration(voronoiDiagram, agentArray,
										actionList, selected_action,
										&myStatistics, Time, EndTime);

								// UPDATE MIGRATION RATE
								double rho_cell = 1035; // kg/m^3
								double rho_liquid = 992.2; // kg/m^3
								double Cd = 0.1;
								double g = 9.81; // m/s^2
								double d = 20 * pow(10, -6); // m
								double k1 = g * (1 - rho_liquid / rho_cell);
								double k2 = Cd * 6 / 8 / d * rho_liquid
										/ rho_cell;
								double velocity = sqrt(k1 / k2)
										* tanh(Time * 60 * 60 * sqrt(k1 * k2)); // m/s
#define pi 3.14159265
								double lattice_constant = d
										* pow(pi / 6., 1. / 3.); // m
								CellMigrationRate =
										velocity / (pi / 4. * lattice_constant)
												* 60 * 60; // h^-1
								//fprintf(stderr, "Velocity: %e [m/s] => Hopping Rate: %e [mu m/h]\n", velocity, CellMigrationRate);
								flag_actualizeProcessRates = TRUE;
							}else
							{
								if( Agent::USE_GRAVITY){
									alpha_ref = 10000;
									CellMigrationRate = 1.17;
									Agent::USE_GRAVITY = false;
									flag_actualizeProcessRates = TRUE;
								}

								//performChemotacticMigration(
								//performFreeMigration
								performMigration(
									voronoiDiagram, agentArray,
									actionList, selected_action,
									&myStatistics, Time, EndTime);

							}
						}
						//diffuse_cell(actionList, selected_action, voronoiDiagram);
						break;

					case GROWTH: {
						//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," GROWTH\n");fclose(fp_dbg );
						if (selected_action->originalCell->state == COMPARTMENT) {
							//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," A: Compartment start\n");fclose(fp_dbg );
							int whichCell =
									(int) myRandE(
											(double) selected_action->originalCell->growingTumorCellCount);
							int sumCells = 0;
							int randM = -1;
							for (int m = 0; m <= M_gro && randM < 0; m++) {
								if (selected_action->internalStateM[m] > 0) {
									sumCells +=
											selected_action->internalStateM[m];
									if (sumCells > whichCell) {
										randM = m;
									}
								}
							}

							if (randM == M_gro) {
								VoronoiCell *destination =
										performGrowthAndDivision(voronoiDiagram,
												agentArray, actionList,
												selected_action, &myStatistics,
												Time, EndTime);
								if (destination != NULL) {
									if (GetVoronoiCell(selected_action->originalCell)->isDomainBorder(
											voronoiDiagram)) {
										EndTime = Time;
									}

									selected_action->internalStateM[randM]--;
									selected_action->internalStateM[0]++;
									GetAgent(destination)->actions[INDEX_GROWTH]->internalStateM[0]++;
								}
							} else {
								selected_action->internalStateM[randM]--;
								selected_action->internalStateM[randM + 1]++;
							}
							//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," A: Compartment end\n");fclose(fp_dbg );
						}else {
							//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B: Something else start\n");fclose(fp_dbg );
							if( selected_action->originalCell->location[0]->waste > Agent::WASTE_THRESHOLD_SLOWED_GROWTH
									|| selected_action->originalCell->getOxygen() < 0.07)
								selected_action->originalCell->intoxication += 1./selected_action->rate;
							//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B: Prior to next decision...\n");fclose(fp_dbg );
							if(selected_action->originalCell->divide){ // <-- DIVISION / GROWTH
								if (selected_action->internalState == M_gro) {
									//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B1: 1 In the light of the fact\n");fclose(fp_dbg );
									VoronoiCell *destination = performGrowth(
											voronoiDiagram, agentArray, actionList,
											selected_action, &myStatistics, Time,
											EndTime);
									//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B1: 2  On the lone and levelo\n");fclose(fp_dbg );
									if (M_div < 0 && destination != NULL/*&& selected_action->originalCell->countLocations >= 2*MIN_SUBCELLULAR_COMPONENTS*/) {
										performDivision(
												voronoiDiagram,
												agentArray,
												actionList,
												GetAgent(destination)->actions[INDEX_DIVISION],
												&myStatistics, Time, EndTime);
										//if (count_cells == MAX_CELLS)
										//	EndTime = Time;
										if (count_cells == voronoiDiagram->countVoronoiCells)
											Time = EndTime;
										if (AdaptEndTime
												&& GetAgent(destination)
														!= selected_action->originalCell)
											if (GetVoronoiCell(selected_action->originalCell)->isDomainBorder(
													voronoiDiagram)) {
												Time = EndTime;
												//EndTime = Time;
											}
									}
									//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B1: 3  Sand stretch far away\n");fclose(fp_dbg );
									selected_action->internalState = 0;
									//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B1: 4  In the heat of the action\n");fclose(fp_dbg );
								}else {
									//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B2: hu?\n");fclose(fp_dbg );
									selected_action->internalState++;
									//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B2: ha!\n");fclose(fp_dbg );
								}
							//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," B: Something else end\n");fclose(fp_dbg );
							}else{ /// <-- QUIESCENCE
								//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," C: Something else start\n");fclose(fp_dbg );
								VoronoiCell *closestVC = voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
										selected_action->originalCell->location[0],
										VoronoiCell::extendedNeighborDepth,
										VoronoiCell::symbolicExtendedNeighborhoodSize,
										VoronoiCell::symbolicExtendedNeighborhood);


								if (closestVC && selected_action->originalCell->countFreeNeighbors()){
									double min_dist =	closestVC->getDistanceSquareTo(
											selected_action->originalCell->location[0]);

									double prob_reentering_cell_cycle =
													PROB_REENTERING_CELL_CYCLE( sqrt((double)min_dist) * SPATIAL_UNIT);

									// INTOXICATION -> QUIESCENCE
									if( selected_action->originalCell->intoxication > Agent::WASTE_INTOXICATED_CELL_CYCLES/MaxCellDivisionRate)
										prob_reentering_cell_cycle=0;

									if (
													( myRand() <= prob_reentering_cell_cycle
													&& VoronoiCell::ECM_THRESHOLD_QUIESCENCE <= selected_action->originalCell->location[0]->ecm))
										selected_action->originalCell->divide = 1.;

								}
								//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," C: Something else end\n");fclose(fp_dbg );
							} /// <--
						}
					}
					//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," GROWTH END\n");fclose(fp_dbg );
					break;


					case DIVISION:
						//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," DIVISION\n");fclose(fp_dbg );fp_dbg = fopen( "outdir/dbg.txt", "a+");						
						//fprintf(stderr, "DIVISION?\n");
						if (selected_action->internalState == M_div) {
							selected_action->internalState = 0;
							performDivision(voronoiDiagram, agentArray,
									actionList, selected_action, &myStatistics,
									Time, EndTime);
						} else {
							selected_action->internalState++;
						}

						if (AdaptEndTime)
							if (!(
									GetPosition( selected_action->originalCell)[0]
									        > voronoiDiagram->xMin[0]
											        + voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[0]
											< voronoiDiagram->xMax[0]
													- voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[1]
											> voronoiDiagram->xMin[1]
													+ voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[1]
											< voronoiDiagram->xMax[1]
#if DIMENSIONS == 3
													- voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[2]
											> voronoiDiagram->xMin[2]
													+ voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[2]
											< voronoiDiagram->xMax[2]
													- voronoiDiagram->boundaryThickness
#endif
							)
						) {
								Time = EndTime;
								//EndTime = Time;
							}
						break;

					case NECROSIS:
						//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ,"NECROSIS\n");fclose(fp_dbg );
						if (selected_action->internalState == M_nec) {
							if (selected_action->originalCell->state
									== COMPARTMENT) {
								//fprintf(stderr, "NECROSIS\n");//exit( 0);

								// CHOSE DYING CELL
								int whichCell =
										(int) myRandE(
												(double) selected_action->originalCell->growingTumorCellCount);
								int sumCells = 0;
								int randM = -1;
								for (int m = 0; m <= M_gro && randM < 0; m++) {
									if (selected_action->internalStateM[m] > 0) {
										sumCells +=
												selected_action->internalStateM[m];
										if (sumCells > whichCell) {
											randM = m;
										}
									}
								}

								// CHANGE STATE
								selected_action->originalCell->growingTumorCellCount--;
								selected_action->originalCell->actions[INDEX_GROWTH]->internalStateM[randM]--;
								selected_action->originalCell->necroticCellCount++;

								// ACTUALIZE ACTIONS OF CELL: NECROSIS & GROWTH
								actionList->actualizeRate(
										selected_action->originalCell->actions[INDEX_GROWTH],
										selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
								actionList->actualizeRate(
										selected_action->originalCell->actions[INDEX_NECROSIS],
										selected_action->originalCell->actions[INDEX_NECROSIS]->getActualRate());
								selected_action->originalCell->actualize(
										actionList);

								// STATISTICS
								count_necrotic_cells++;
								count_necrotic_volume +=
										selected_action->originalCell->countLocations;
								count_cells--;
								count_cell_volume -=
										selected_action->originalCell->countLocations;
							}

							else {
								//fprintf(stderr, "NECROSIS\n");

								// CHANGE STATE
								selected_action->originalCell->state = NECROTIC;

								// ACTUALIZE ACTIONS OF CELL
								if (selected_action->originalCell->growingTumorCellCount
										> 0)
									selected_action->originalCell->growingTumorCellCount--;
								else
									selected_action->originalCell->dividingTumorCellCount--;
								selected_action->originalCell->necroticCellCount++;
								selected_action->originalCell->actualize(
										actionList);

								// STATISTICS
								count_necrotic_cells++;
								count_necrotic_volume +=
										selected_action->originalCell->countLocations;
								count_cells--;
								count_cell_volume -=
										selected_action->originalCell->countLocations;
							}

							//fprintf( stderr, "INFO: new cell state: %s, cell count: %i\n", cellTypeToString(selected_action->originalCell->state), selected_action->originalCell->cellCount);
						} else {
							selected_action->internalState++;
						}

						break;
						case LYSIS:
						//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," LYSIS\n");fclose(fp_dbg );
						count_lysed_volume +=
								selected_action->originalCell->countLocations;
						count_necrotic_volume -=
								selected_action->originalCell->countLocations;
						count_inner_empty_volume +=
								selected_action->originalCell->countLocations;
#if MULTISCALE
						if( selected_action->originalCell->cellCount-- == selected_action->originalCell->maxCellCount) {
							update_surrounding_deleted_cell(actionList, selected_action->originalCell, voronoiDiagram);
							if( selected_action->originalCell->state == NONACTIVE)
							addDivisionAction( actionList, selected_action->originalCell, voronoiDiagram);
							selected_action->originalCell->state = FREE;
						}
						if( selected_action->originalCell->cellCount == 0) {
							update_deleted_cell( actionList, selected_action->originalCell, voronoiDiagram);
						}
#else

						//fprintf( stderr, "INFO: do single based lysis\n");
						if (selected_action->originalCell->state != COMPARTMENT) {

							//update_surrounding_deleted_cell(actionList, selected_action->originalCell, voronoiDiagram);
							update_surrounding_removed_cell(actionList,
									selected_action->originalCell,
									voronoiDiagram);
							//deleteAction( actionList, selected_action);
							selected_action->originalCell->state = FREE;
							selected_action->originalCell->actualize(
									actionList);
						}
						//selected_action->originalCell->detach();
						//agentArray->deactivateAgent( selected_action->originalCell);
#endif

						selected_action->originalCell->necroticCellCount--;
						selected_action->originalCell->cellCount--;
						count_necrotic_cells--;
						count_lysed_cells++;

						//flag_actualizeProcessRates = TRUE;

						break;

					case GAP:
						//fprintf(stderr, "\n GAP \n");
						//Time = Last_Time + 1. / selected_action->rate;
						break;

					default:
						fprintf(
								stderr,
								"\nError: Selected action (Type: %i, %s) is unknown!\n",
								selected_action->type,
								actionTypeToString(selected_action->type));
						exit(0);
					}

					timeExecution +=
							(double) (clock() - passedTime)
									/ (double) CLOCKS_PER_SEC;

				} else {
					if (selected_action->rate == 0.) {
						//	fprintf( stderr, "selected_action->rate = %lf\n", selected_action->rate);
						Time = Last_Time + TIME_STEP * 1.;
					}
					if (Time > EndTime) {
						//	fprintf( stderr, "Time(%lf)>EndTime(%lf)\n", Time, EndTime);
						// correct time
						Time = EndTime;
					}
				}
				//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Execute action done\n");fclose(fp_dbg );
			} else {
				//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Empty list\n");fclose(fp_dbg );
				// Empty Probability List
				//Last_Time = Time;
				if (actionList->rateSum == 0.) {
					Time += TIME_STEP * 10.;
				}
				if (actionList->size == 0) {
					Time = EndTime;
				}
			}
			//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Line 2583 \n");fclose(fp_dbg );
			if (Time > EndTime) {
				// correct time
				Time = EndTime;
			}
			//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Line 2589 \n");fclose(fp_dbg )
			
			if (Time != EndTime && Case != 1 /*&& count_inoculated_cells + (int)count_divided_cells - count_necrotic_cells - count_lysed_cells!=0*/) {
				timeDifference += Time - Last_Time;

				// UPDATE MORPHOGEN
				if (timeDifference > customTimestep)
				if (VoronoiCell::USE_MORPHOGEN) {
					morphogenDynamics->setBoundaryCondition( DiffusionReactionEquation::DIRICHLET);
					for (int i = 0; i < voronoiDiagram->countVoronoiCells;
							i++) {

						// GENERATION
						if (GetAgent(voronoiDiagram->voronoiCells[i])
								&& voronoiDiagram->voronoiCells[i]->getState() != FREE
								//&& voronoiDiagram->voronoiCells[i]->getState() == NECROTIC
										)
							morphogenProduction[i] = -300;

						// DECAY
						morphogenDecay[i] = 100;

					}

					// update morphogen
					morphogenDynamics->update(timeDifference);

					// hard copy
					for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++)
						voronoiDiagram->voronoiCells[i]->morphogen =
								morphogen[i];
				}

				// UPDATE GLUCOSE & OXYGEN & LACTATE

				if (timeDifference > customTimestep) {

					// UPDATE GLUCOSE & OXYGEN

					// CG for linear system
					passedTime = clock();
					if (Case == 5) {
						// Diffusion Reaction of Growth Factors
						//substrate.Produce_Growthfactors();
						UpdateGrowthFactorsNonLinearCGSparse(voronoiDiagram,
								Time - timeDifference, Time, customTimestep);

						// Update Vessel Network
						if (substrate.vasculature.Update_Network(timeDifference,
								actionList, voronoiDiagram)) {
							// Update Pressure in Vessels
							substrate.vasculature.GiveMeTheBlood();
							// Reset Initial Oxygen Concentration in Vessels
							substrate.Refill_Oxygen(Case);
							// Reset Initial Glucose Concentration in Vessels
							substrate.Refill_Glucose(Case);
						}
					}

					// Diffusion Reaction of Glucose & Oxygen
					timeDifference =
							UpdateSystemNonLinearCGSparse(voronoiDiagram,
									Time - timeDifference, Time, customTimestep);

					//pharmaco->evolveSystem( voronoiDiagram, timeDifference, customTimestep);
					//timeDifference = UpdateSystemNewtonCGSparse( voronoiDiagram, Time - timeDifference, Time, customTimestep);

					fprintf(stderr,
							"...finished ( %lisec) -> time rest = %lf\n",
							(clock() - passedTime) / CLOCKS_PER_SEC,
							timeDifference);
					benchmarkTime += (clock() - passedTime);
					//	exit(0);
					// Newton for non-linear system
					//timeDifference = UpdateSystemNewtonSparse( voronoiDiagram, timeDifference/10, timeDifference);

					//timeDifference = UpdateSystemImplicitSparse( voronoiDiagram, customTimestep, timeDifference);
					flag_actualizeProcessRates = TRUE;
					//timeDifference = 0.;

					// UPDATE LACTATE
					if (VoronoiCell::USE_LACTATE) {
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++) {
							lactateProduction[i] =
									-GetLactateProductionRate(
											voronoiDiagram->voronoiCells[i],
											voronoiDiagram->voronoiCells[i]->glucose,
											voronoiDiagram->voronoiCells[i]->oxygen);
							if (voronoiDiagram->voronoiCells[i]->isFree())
								lactateDiffusion[i] = 30 * Lactate_Diffusion;
							else
								lactateDiffusion[i] = Lactate_Diffusion;
						}

						// update lactate
						lactateDynamics->update();

						// hard copy
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++)
							voronoiDiagram->voronoiCells[i]->lactate =
									lactate[i];

					}

					// UPDATE WASTE
					if (VoronoiCell::USE_WASTE) {
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++) {

							// waste production
							wasteProduction[i] = (voronoiDiagram->voronoiCells[i]->getState() == NECROTIC ? -10 : 0);
							// waste uptake
							if (!voronoiDiagram->voronoiCells[i]->isFree() && voronoiDiagram->voronoiCells[i]->getState() != NECROTIC){
								wasteUptake[i] = Agent::WASTE_UPTAKE;
								GetAgent( voronoiDiagram->voronoiCells[i])->waste += Agent::WASTE_UPTAKE * waste[i] * customTimestep;
							}else{
								wasteUptake[i] = 0;
							}


							if (voronoiDiagram->voronoiCells[i]->isFree())
								wasteDiffusion[i] = 30 * custom_Waste_Diffusion;
							else{
								wasteDiffusion[i] = custom_Waste_Diffusion;
							}
						}

						// update waste
						wasteDynamics->update();

						// hard copy
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++)
							voronoiDiagram->voronoiCells[i]->waste = waste[i];

					}

				}

			}
			//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg ," Line 2731, pre cont model\n");fclose(fp_dbg );
			// CONTINUUM MODEL (ECM)
			if (VoronoiCell::ECM_THRESHOLD_QUIESCENCE != 0) {
				double TimeStep = 0.01;
				for (; TempTime <= Time; TempTime += TimeStep) {

					for (int v = 0; v < voronoiDiagram->countVoronoiCells;
							v++) {
						// DEGRADATION
						dECM[v] =
								-(TimeStep) * VoronoiCell::ECM_DEGRADATION_RATE
										* voronoiDiagram->voronoiCells[v]->ecm;

						// DIFFUSION
						for (int n = 0;
								n
										< voronoiDiagram->voronoiCells[v]->countNeighborCells;
								n++) {
						}
					}

					for (int a = 0; a < agentArray->countActiveAgents; a++) {
						for (int l = 0;
								l < agentArray->agents[a]->countLocations; l++)
							if (agentArray->agents[a]->state == ACTIVE
									|| agentArray->agents[a]->state == NONACTIVE) {

								// PRODUCTION
								dECM[agentArray->agents[a]->location[l]->index] +=
										(TimeStep) * VoronoiCell::ECM_PRODUCTION_RATE / 2.
								;

								for (int n = 0;
										n
												< agentArray->agents[a]->location[l]->countNeighborCells;
										n++) {
									dECM[agentArray->agents[a]->location[l]->neighborCells[n]->index] +=
											(TimeStep) * VoronoiCell::ECM_PRODUCTION_RATE
													/ 2.
													/ (double) agentArray->agents[a]->location[l]->countNeighborCells;
								}

							}
					}

					for (int v = 0; v < voronoiDiagram->countVoronoiCells;
							v++) {
						// DEGRADATION
						voronoiDiagram->voronoiCells[v]->ecm += dECM[v];
					}
				}
			}
			
			//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg,"Output Start\n");fclose(fp_dbg );
			// OUTPUT
			if( RadialProfilesCount && inbound<double>( RadialProfilesTime, floor(Last_Time/24.), floor(Time/24.), RadialProfilesCount, 1) ) {

				// UPDATE
				if( ceil(Last_Time/1) != ceil(Time/1)){
					updateHistogram( agentArray, voronoiDiagram,	histogram, histogramDividing, histogramNecrotic, histogramFree, histogramECM);
				}

				for( int day=(int)ceil(Last_Time/24.); day < (int)ceil(Time/24.); day++ )
				if(	inbound<double>( RadialProfilesTime, day-1, day-1, RadialProfilesCount, 1))
				{
					// Dennis ECM
					if(writeRawSimulation){
						//FILE *fp_raw = fopen( ecm_tmpout, "w");
						double ecm_out_val;
						for( int j=0; j<profileDepth; j++){
							if( histogram[j])
								ecm_out_val = (double) histogramECM[j]	/ (double) histogram[j];
							else
								ecm_out_val = 0;
							ecm_out_val = max( 0.,  ecm_out_val + measurement_error_ecm*normrnd());
							//fprintf( fp_raw, "%i %e \n", j, ecm_out_val);
							ecm_out.push_back(ecm_out_val);
						}
						//fclose(fp_raw);
					}
					
					// Dennis proliferation
					if(writeRawSimulation){
						//FILE *fp_raw = fopen(prolif_tmpout, "w");
						double prolif_out_val;
						for( int j=0; j<profileDepth; j++){
							prolif_out_val = ( histogram[j] - histogramFree[j]>0 ? (double) histogramDividing[j] / (double) (histogram[j] - histogramFree[j]) : 0);
							prolif_out_val = max( 0.,  prolif_out_val + measurement_error_ki67*normrnd());
						//	fprintf( fp_raw, "%i %e\n", j, prolif_out_val);
							prolif_out.push_back(prolif_out_val);
						}
						//fclose(fp_raw);
					}					

					// RESET
					for( int i=0; i<10000; i++){
						histogram[i]=0;
						histogramDividing[i]=0;
						histogramNecrotic[i]=0;
						histogramFree[i]=0;
						histogramECM[i]=0;
					}
				}
			}

			// actual timestep
			if (indexOfTime( Last_Time, BeginningTime, OutputRate)
					< indexOfTime( Time, BeginningTime, OutputRate)) {
				if(writeRawComparison && data_growthcurve.dim){
					// passed data points
					for( ; data_growthcurve.x[idxEpsilon] < Time && idxEpsilon<data_growthcurve.dim; idxEpsilon++){
						double radius = ( isnan(gyrRadius) ? sqrt(DIMENSIONS*50*50) : sqrt(gyrRadius) ) * AGENT_DIAMETER * 0.8;

						data_growthcurve.m[idxEpsilon] = max( 0.,  data_growthcurve.m[idxEpsilon] + measurement_error_gc*normrnd());

						cumEpsilon += 0.5 * pow( (data_growthcurve.m[idxEpsilon] - radius)/data_growthcurve.s[idxEpsilon], 2)  / data_growthcurve.dim;

						sprintf(outfilename, "%s/raw_gc.dat", dirname);
						FILE *fp_raw = fopen( outfilename, "a+");
						fprintf( fp_raw, "%e %e %e %e\n", data_growthcurve.x[idxEpsilon], radius, data_growthcurve.m[idxEpsilon], data_growthcurve.s[idxEpsilon]);
						fclose(fp_raw);

						if( isnan(cumEpsilon)){
							fprintf(stderr, "{GC}!\n"); exit(0);
						}

						if( cumEpsilon > maxEpsilon)
							return cumEpsilon;
					}
				}

				Last_gyrRadius = gyrRadius;
				gyrRadius = getGyrationRadiusOfBorder(agentArray, voronoiDiagram);				
				
				if(writeRawSimulation){
					int l;									
					double radius = ( isnan(gyrRadius) ? sqrt(DIMENSIONS*50*50) : sqrt(gyrRadius) ) * AGENT_DIAMETER * 0.8;
					//FILE *fp_raw = fopen( gc_tmpout, "a+");
					for (l = indexOfTime( Last_Time, BeginningTime, OutputRate) + 1;
						 l < indexOfTime( Time, BeginningTime, OutputRate);
						 l++) {
						//fprintf(fp_raw, "%e %e\n", timeOfIndex ( l, BeginningTime, OutputRate),radius);
						gc_out.push_back(radius);
					}
					l = indexOfTime( Time, BeginningTime, OutputRate);
					gc_out.push_back(radius);					
					//fprintf(fp_raw, "%e %e\n", timeOfIndex ( l, BeginningTime, OutputRate), radius );						
					//fclose(fp_raw);
				}
				
				
				Last_gyrRadius = gyrRadius;
				gyrRadius = getGyrationRadiusOfBorder(agentArray, voronoiDiagram);


				// fill passed time steps
				for (l = indexOfTime( Last_Time, BeginningTime, OutputRate) + 1;
						l < indexOfTime( Time, BeginningTime, OutputRate);
						l++) {

					// update arrays
					global_cells[l] += Last_NumberOfCells; //count_cells;
					global_volume[l] += count_cell_volume;
					global_necrotic_cells[l] += count_necrotic_cells;
					global_necrotic_volume[l] += count_necrotic_volume;
					global_inner_empty_volume[l] += count_inner_empty_volume;
					global_vessel_cells[l] += count_vessel_cells;
					global_gyrRadius[l] += sqrt(Last_gyrRadius);
					global_gyrRadiusSquare[l] += Last_gyrRadius;
					global_avgGlucose[l] += getAvgGlucose(agentArray);
					global_avgOxygen[l] += getAvgOxygen(agentArray);
				}

				// calculate array position for actual time   
				l = indexOfTime( Time, BeginningTime, OutputRate);

				// update arrays
				global_cells[l] += count_cells;
				global_volume[l] += count_cell_volume;
				global_necrotic_cells[l] += count_necrotic_cells;
				global_necrotic_volume[l] += count_necrotic_volume;
				global_inner_empty_volume[l] += count_inner_empty_volume;
				global_vessel_cells[l] += count_vessel_cells;
				global_gyrRadius[l] += sqrt(gyrRadius);
				global_gyrRadiusSquare[l] += gyrRadius;
				global_avgGlucose[l] += getAvgGlucose(agentArray);
				global_avgOxygen[l] += getAvgOxygen(agentArray);

			}
			//fp_dbg = fopen( "outdir/dbg.txt", "a+");fprintf(fp_dbg,"Output End\n");fclose(fp_dbg);
		} while (Time < EndTime && actionList->size != 0);
		/** END OF SIMULATION **/

		// FREE MEMORY
		// DESTROY ACTION LIST
		actionList->destroyActionTree();

		//fprintf(stderr,"DETACH, DEACTIVATE AND REINIT AGENTS\n");
		// DETACH, DEACTIVATE AND REINIT AGENTS
		while (agentArray->countActiveAgents != 0) {
			agentArray->agents[0]->detach();
			agentArray->deactivateAgent(agentArray->agents[0]);
		}
		for (i = 0; i < agentArray->countAgents; i++)
			if (agentArray->agents[i]->actionsInitialized)
				for (int ii = 0; ii <= INDEX_GROWTH; ii++)
					agentArray->agents[i]->actions[ii]->internalState = 0;

		// REINIT VORONOI CELLS
		for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			voronoiDiagram->voronoiCells[i]->countFreeNeighborCells =
					voronoiDiagram->voronoiCells[i]->countNeighborCells;
			voronoiDiagram->voronoiCells[i]->countFreeExtendedNeighborCells =
					voronoiDiagram->voronoiCells[i]->countExtendedNeighborCells;
		}

	} // END AVERAGES


	/** END OF STOCHASTIC APPROACHE ******************************************/
	/* Get the current time. */
	temp_t = time(NULL);
	loctime = localtime_r(&temp_t, loctime);
	strftime(local_time, 256, "%A, %B %d. %k:%M:%S ", loctime);
	free(loctime);

	// ######## CLEAN STUFF UP ###########
	// Agents	
	for (i = 0; i < agentArray->countAgents; i++) {
		if (agentArray->agents[i]->actionsInitialized) {
			for (k = 0; k < INDEX_MIGRATION + USE_MIGRATION; k++)
				free(agentArray->agents[i]->actions[k]);
			free(agentArray->agents[i]->actions);
			agentArray->agents[i]->actionsInitialized = FALSE;
		}
		if (agentArray->agents[i]->location != NULL)
			free(agentArray->agents[i]->location);
		free(agentArray->agents[i]);
	}
	free(agentArray->agents);
	free(agentArray);

	// Voronoi Grids
	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		free(voronoiDiagram->voronoiCells[i]->extendedNeighborhood);
	}

	// Statistical Variables
	free(global_cells);
	free(global_necrotic_cells);
	free(global_necrotic_volume);
	free(global_vessel_cells);
	free(global_volume);
	free(global_inner_empty_volume);
	free(global_gyrRadius);
	free(global_gyrRadiusSquare);
	free(global_avgOxygen);
	free(global_avgGlucose);


	// Free Continuum Part
	if (Case > 1) {
		sA->deleteSparseMatrix(sA);
		free(b);
		free(x);
		free(v0);
		free(v1);
		free(v2);
		free(v3);
		free(v4);
		free(v5);
		free(v6);
		free(v7);
		free(v8);
	}
	return cumEpsilon;
}
/****************************************************************************/

double DESubstrate(double substrate, double cells) {
	double cellGrowthRate = MaxCellDivisionRate, Ks = 0.2, // (mmol / L)
			ms = 0.0056, // (mmol / L * h * 10^5 cells)
			Ys = 1.37; // ()

	double apparentGrowthRate = cellGrowthRate * substrate / (Ks + substrate),
			qs = apparentGrowthRate / Ys + ms;

	return -qs * (double) cells; // / 100000;
	//return -( (double)substrate * cellGrowthRate / (Ks + (double)substrate) / Ys + ms ) * (double)cells;
}
/****************************************************************************/

double rungeKutta(double(*equation)(double, double), double xStart, double xEnd,
		double yStart, double z, int steps) {
	double stepSize = (xEnd - xStart) / (double) steps;
	int i;

	for (i = 0; i < steps; i++) {
		yStart = rungeKuttaStep(equation, xStart, xStart + stepSize, yStart, z);
	}

	return yStart;
}
/****************************************************************************/

double rungeKuttaStep(double(*equation)(double, double), double xStart,
		double xEnd, double yStart, double z) {
	double stepSize = xEnd - xStart;
	double yA, yB, yC, dyStart, dyA, dyB, dyC;

	dyStart = equation(yStart, z);

	yA = yStart + stepSize / 2. * dyStart;
	dyA = equation(xStart + stepSize / 2., yA);

	yB = yStart + stepSize / 2. * dyA;
	dyB = equation(xStart + stepSize / 2., yB);

	yC = yStart + stepSize * dyB;
	dyC = equation(xStart + stepSize, yC);

	return yStart + stepSize / 6. * (dyStart + 2. * (dyA + dyB) + dyC);
}
/****************************************************************************/


#define INTRA_AGENT_ENERGY 4
#define INTER_AGENT_ENERGY 2 //0.5
////////////////////////////////////////
double getEnergyDifferenceCadherin(VoronoiCell *oldPosition,
		VoronoiCell *newPosition) {
	double energyBilance = 0;

	// loss of contact
	for (int iii = 0; iii < oldPosition->countNeighborCells; iii++) {
		if (oldPosition->neighborCells[iii]->getState() != FREE) {
			if (GetAgent( oldPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance -= INTRA_AGENT_ENERGY;
			else
				energyBilance -= INTER_AGENT_ENERGY;
		}
	}

	// gain new contact
	for (int iii = 0; iii < newPosition->countNeighborCells; iii++) {
		if (newPosition->neighborCells[iii]->getState() != FREE
				&& newPosition->neighborCells[iii] != oldPosition) {
			if (GetAgent( newPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance += INTRA_AGENT_ENERGY;
			else
				energyBilance += INTER_AGENT_ENERGY;
		}
	}

	return energyBilance;
}
//////////////////////////////////////

double getEnergyRemove(VoronoiCell *oldPosition) {
	double energy = 0;

	for (int iii = 0; iii < oldPosition->countNeighborCells; iii++) {
		if (oldPosition->neighborCells[iii]->getState() != FREE) {
			if (GetAgent( oldPosition->neighborCells[iii])
					== GetAgent( oldPosition)) {
				// cell deformation
				energy -= INTRA_AGENT_ENERGY;
			} else {
				// loss of cell-cell-contact
				if (GetAgent( oldPosition)->getGlucose()
						* GetAgent( oldPosition)->getOxygen()
						> THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)
					energy -= INTER_AGENT_ENERGY;
			}
		}
	}

	return energy;
}
//////////////////////////////////////

double getEnergyPlaced(VoronoiCell *newPosition, VoronoiCell *oldPosition) {
	double energyBilance = 0;

	//fprintf( stderr, "getEnergyPlaced()\n");

	for (int iii = 0; iii < newPosition->countNeighborCells; iii++) {
		if (newPosition->neighborCells[iii]->getState() != FREE
				&& newPosition->neighborCells[iii] != oldPosition) {
			if (GetAgent( newPosition->neighborCells[iii])
					== GetAgent( oldPosition)) {
				// cell relaxation
				energyBilance += INTRA_AGENT_ENERGY;
			} else {
				// gain new contact
				//if( GetAgent( oldPosition)->getGlucose() * GetAgent( oldPosition)->getOxygen() > THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)
				energyBilance += INTER_AGENT_ENERGY;
			}
		}
	}

	return energyBilance;
}
//////////////////////////////////////

double getEnergyDifference(VoronoiCell *oldPosition, VoronoiCell *newPosition) {
	double energyBilance = 0;

	// loss of contact
	for (int iii = 0; iii < oldPosition->countNeighborCells; iii++) {
		if (oldPosition->neighborCells[iii]->getState() != FREE) {
			if (GetAgent( oldPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance -= INTRA_AGENT_ENERGY;
			else
				energyBilance -= INTER_AGENT_ENERGY;
		}
	}

	// gain new contact
	for (int iii = 0; iii < newPosition->countNeighborCells; iii++) {
		if (newPosition->neighborCells[iii]->getState() != FREE
				&& newPosition->neighborCells[iii] != oldPosition) {
			if (GetAgent( newPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance += INTRA_AGENT_ENERGY;
			else
				energyBilance += INTER_AGENT_ENERGY;
		}
	}

	return energyBilance;
}

double getConnections(VoronoiCell *start, VoronoiCell *ignore, Agent *agent) {
	int ii;
	int stackLength = 1;
	int stackPosition = 0;
	VoronoiCell *actualLocation; // = agentArray->agents[i]->location[0];
	VoronoiCell *locationStack[100];
	locationStack[0] = start;

	//fprintf( stderr, "getConnections()\n");

	do {
		actualLocation = locationStack[stackPosition];
		stackPosition++;
		//int iii;
		//fprintf( stderr, "getConnections(): stack position %i\n", stackPosition);
		for (ii = 0; ii < actualLocation->countNeighborCells; ii++) {
			// neighbor is occupied by same neighbor?
			if (actualLocation->neighborCells[ii]->getState() != FREE
					&& actualLocation->neighborCells[ii] != ignore
					&& agent == GetAgent( actualLocation->neighborCells[ii])) {
				// is already in stack?
				int iv;
				for (iv = 0;
						iv < stackLength
								&& locationStack[iv]
										!= actualLocation->neighborCells[ii];
						iv++)
					;
				if (iv == stackLength) {
					locationStack[stackLength++] =
							actualLocation->neighborCells[ii];
				}
			}
		}
	} while (stackLength != stackPosition);

	//fprintf( stderr, "getConnections() end\n");
	return stackLength;
}


void performMigration(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	VoronoiCell *oldLocation = NULL, *newLocation = NULL;

	// ENERGY-MINIMIZING MIGRATION

	// count candidates
	int countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0;
				ii
						< selected_action->originalCell->location[i]->countNeighborCells;
				ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState()
					== FREE)
				countFreeNeighborSites++;

	// calculate energies
	double energies[countFreeNeighborSites];
	double energySum = 0;
	countFreeNeighborSites = 0;
	double center[2] = {50,50};
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0;
				ii
						< selected_action->originalCell->location[i]->countNeighborCells;
				ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE &&
					selected_action->originalCell->location[i]->neighborCells[ii]->getDistanceTo(center)<50.) {
				// each contact to other cell
				double energyDiff =
						getEnergyDifference(
								selected_action->originalCell->location[i],
								selected_action->originalCell->location[i]->neighborCells[ii]);

				if(Agent::USE_GRAVITY){
					double dx = selected_action->originalCell->location[i]->neighborCells[ii]->position[0] - selected_action->originalCell->location[i]->position[0];
					double dy = selected_action->originalCell->location[i]->neighborCells[ii]->position[1] - selected_action->originalCell->location[i]->position[1];
#if DIMENSIONS == 3
					double dz = selected_action->originalCell->location[i]->neighborCells[ii]->position[2] - selected_action->originalCell->location[i]->position[2];
#else
					double dz=0;
#endif
					double vx=0, vy=-1, vz=0;

					double gravitation = exp( - acos(vy*dy / sqrt(vx*vx+vy*vy+vz*vz) / sqrt(dx*dx+dy*dy+dz*dz)) / alpha_ref);

					energies[countFreeNeighborSites] = gravitation;
				}else{
					// NOTOX
					double gravitation = 0;
					double temperature = 1e-2;
					energies[countFreeNeighborSites] =
							exp( (-gravitation*10 + energyDiff) / temperature);

				}

				if (getConnections(
						selected_action->originalCell->location[i]->neighborCells[ii],
						selected_action->originalCell->location[i],
						selected_action->originalCell)
						< getConnections(
								selected_action->originalCell->location[i],
								selected_action->originalCell->location[i]->neighborCells[ii],
								selected_action->originalCell))
					energies[countFreeNeighborSites] = 0;
				
				energySum += energies[countFreeNeighborSites];
				countFreeNeighborSites++;
			}

	// chose winner
	double randSum = myRand() * energySum;
	energySum = 0;
	char foundWinner = FALSE;
	countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0;
				ii
						< selected_action->originalCell->location[i]->countNeighborCells;
				ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE  &&
					selected_action->originalCell->location[i]->neighborCells[ii]->getDistanceTo(center)<50.) {
				// each contact to other cell
				energySum += energies[countFreeNeighborSites];
				if (!foundWinner && energySum > randSum) {
					oldLocation = selected_action->originalCell->location[i];
					newLocation =
							selected_action->originalCell->location[i]->neighborCells[ii];
					foundWinner = TRUE;
				}

				countFreeNeighborSites++;
			}

	// ENERGY-MINIMIZING MIGRATION INCLUDING CASE OF NO MIGRATION
	//fprintf( stderr, "Migrate: %lf\n", choiceProb/totalProbSum);
	// Perform Migration
	if (newLocation != NULL) {

		selected_action->originalCell->detach(oldLocation);
		update_surrounding_reduced_cell(actionList, oldLocation,
				voronoiDiagram);

		// exchange FREE and migrating cell
		if (newLocation->agent != NULL) {
			Agent *freeAgent = GetAgent( newLocation);
			freeAgent->detach(newLocation);
			freeAgent->attach(oldLocation);
		}
		// place FREE agent and migrating cell
		else {
			Agent* oldAgent = agentArray->activateAgent();
			oldAgent->attach(oldLocation);
			oldAgent->state = FREE;
		}

		selected_action->originalCell->attach(newLocation);
		update_surrounding_expanded_cell(actionList, newLocation,
				voronoiDiagram);

		// exchange concentrations
		double temp_oxy = newLocation->oxygen;
		double temp_glu = newLocation->glucose;
		newLocation->oxygen = oldLocation->oxygen;
		newLocation->glucose = oldLocation->glucose;
		oldLocation->oxygen = temp_oxy;
		oldLocation->glucose = temp_glu;
	}

	//fprintf( stderr, "END MIGRATE\n");
}
/****************************************************************************/

void performChemotacticMigration(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	//fprintf( stderr, "MIGRATE\n");

	VoronoiCell *oldLocation = NULL, *newLocation = NULL;

	// ENERGY-MINIMIZING MIGRATION

	//LUNGSYS
	//if (selected_action->originalCell->location[0]->morphogen < 1e-9)
	//	return;

	// count candidates
	int countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0; ii < selected_action->originalCell->location[i]->countNeighborCells; ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState()== FREE)
				countFreeNeighborSites++;

	// calculate energies
	double energies[countFreeNeighborSites];
	double energySum = 0;
	countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0; ii	< selected_action->originalCell->location[i]->countNeighborCells; ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState()== FREE) {
				// each contact to other cell
				double temperature = 1e-5; // NOTOX

				energies[countFreeNeighborSites] =
						exp( ((double)selected_action->originalCell->location[i]->neighborCells[ii]->morphogen
							- (double)selected_action->originalCell->location[i]->morphogen) / temperature);
				if(isinf( energies[countFreeNeighborSites])){
					printf("ERROR: infinit energy! Decreasing the temperature might help!\n");
					energies[countFreeNeighborSites] = (double)FLT_MAX;
				}
				if(energies[countFreeNeighborSites]==0){
					printf("ERROR: zero energy! Decreasing the temperature might help!\n");
					energies[countFreeNeighborSites] = 1./(double)FLT_MAX;
				}

				energySum += energies[countFreeNeighborSites];
				countFreeNeighborSites++;
			}

	//LUNGSYS
	double randSum = myRand() * energySum;
	energySum = 0;
	bool foundWinner = false;
	countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0; ii < selected_action->originalCell->location[i]->countNeighborCells; ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE) {
				// each contact to other cell
				energySum += energies[countFreeNeighborSites];
				if (!foundWinner && energySum >= randSum) {
					oldLocation = selected_action->originalCell->location[i];
					newLocation =
							selected_action->originalCell->location[i]->neighborCells[ii];
					foundWinner = true;
				}

				countFreeNeighborSites++;
			}

	if( !foundWinner && countFreeNeighborSites){
		printf("ERROR: energySum=%e, randSum=%e, countFreeNeighborSites=%i\n",energySum,randSum, countFreeNeighborSites);
		exit(0);
	}
	// ENERGY-MINIMIZING MIGRATION INCLUDING CASE OF NO MIGRATION
	// Perform Migration
	if (newLocation != NULL) {

		selected_action->originalCell->detach(oldLocation);
		update_surrounding_reduced_cell(actionList, oldLocation,
				voronoiDiagram);

		// exchange FREE and migrating cell
		if (newLocation->agent != NULL) {
			Agent *freeAgent = GetAgent( newLocation);
			freeAgent->detach(newLocation);
			freeAgent->attach(oldLocation);
		}
		// place FREE agent and migrating cell
		else {
			Agent* oldAgent = agentArray->activateAgent();
			oldAgent->attach(oldLocation);
			oldAgent->state = FREE;
		}

		selected_action->originalCell->attach(newLocation);
		update_surrounding_expanded_cell(actionList, newLocation,
				voronoiDiagram);

		// exchange concentrations
		double temp_oxy = newLocation->oxygen;
		double temp_glu = newLocation->glucose;
		newLocation->oxygen = oldLocation->oxygen;
		newLocation->glucose = oldLocation->glucose;
		oldLocation->oxygen = temp_oxy;
		oldLocation->glucose = temp_glu;
	}

	//fprintf( stderr, "END MIGRATE\n");
}
/****************************************************************************/

VoronoiCell* performGrowthAndDivision(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	if (selected_action->originalCell->state == COMPARTMENT) {
		Agent *daughterCell = NULL;
		//printf("%i %lf\n", selected_action->originalCell->growingTumorCellCount, selected_action->rate);

		// VOLUME EXPANSION
		if (selected_action->originalCell->cellCount
				< selected_action->originalCell->maxCellCount) {
			// into same compartment	
			//fprintf( stderr, "GROW cell %i (%i=%ig+%id/%i)\n", selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->growingTumorCellCount, selected_action->originalCell->dividingTumorCellCount, selected_action->originalCell->maxCellCount);	
			daughterCell = selected_action->originalCell;

			daughterCell->cellCount++;
			daughterCell->growingTumorCellCount++;
			daughterCell->actualize(actionList);
		} else {
			int i;

			// count free neighbors
			int countFreeNeighbors = 0;
			for (i = 0;
					i
							< GetVoronoiCell(selected_action->originalCell)->countNeighborCells;
					i++) {
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])
						== NULL) {
					Agent *newAgent = agentArray->activateAgent();
					newAgent->state = COMPARTMENT;
					newAgent->attach(
							GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
					//fprintf( stderr, "ACTIVATE cell %i\n", newAgent->index);
				}

				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount)
					countFreeNeighbors++;
			}
			if (countFreeNeighbors == 0 /*&& countSlightlyFreeNeighbors == 0*/) {
				fprintf(
						stderr,
						"ERROR in performGrowthAndDivision(): countFreeNeighbors = %i\n",
						countFreeNeighbors);
				//exit(0);
				return NULL;
			}

			// chose free neighbor
			//double rand = myRandIE( 0., (double)countFreeNeighbors);
			//int whichFreeNeighbor = (int)(myRand() * (double)countFreeNeighbors);
			int whichFreeNeighbor = (int) (myRandE((double) countFreeNeighbors));
			/*if( whichFreeNeighbor >= countFreeNeighbors){
			 fprintf( stderr, "ERROR in performGrowthAndDivision(): rand(%lf) * countFreeNeighbors(%i) =: whichFreeNeighbor(%i) >= countFreeNeighbors(%i)\n", rand, countFreeNeighbors, whichFreeNeighbor, countFreeNeighbors);
			 exit( 0);
			 }*/
			//fprintf( stderr, "which neighbor %i/%i\n", whichFreeNeighbor, countFreeNeighbors);
			countFreeNeighbors = 0;
			//for( i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells && countFreeNeighbors!=whichFreeNeighbor; i++)
			for (i = 0;
					i
							< GetVoronoiCell(selected_action->originalCell)->countNeighborCells
							&& daughterCell == NULL; i++) {
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount) {
					if (countFreeNeighbors == whichFreeNeighbor)
						daughterCell =
								GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
					countFreeNeighbors++;
				}
			}
			//selected_action->originalCell->cellCount--;
			//selected_action->originalCell->growingTumorCellCount--;

			//fprintf( stderr, "GROW cell %i (%i/%i) to %i (%i/%i)\n", 
			//	selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->maxCellCount, 
			//	daughterCell->index, daughterCell->cellCount, daughterCell->maxCellCount);	

			if (daughterCell == NULL) {
				fprintf(
						stderr,
						"ERROR in performGrowthAndDivision(): daughterCell==NULL\n");
				fprintf(stderr,
						"whichFreeNeighbor(%i) >= countFreeNeighbors(%i) == \n",
						whichFreeNeighbor, countFreeNeighbors);
				exit(0);
			}

			daughterCell->cellCount++;
			daughterCell->growingTumorCellCount++;
			daughterCell->actualize(actionList);

			selected_action->originalCell->actualize(actionList);
		}
		//fprintf( stderr, "DIVISION: %i -> %i\n", selected_action->originalCell->index, daughterCell->index);

		myStatistics->count_expanded_cells[0]++;
		myStatistics->count_divided_cells[0]++;
		myStatistics->count_cell_volume[0]++;
		myStatistics->count_cells[0]++;
		//actionList->actualizeAllRates( voronoiDiagram);
		//actionList->actualizeRate( selected_action, selected_action->getActualRate());
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_GROWTH],
				selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
		actionList->actualizeRate(daughterCell->actions[INDEX_GROWTH],
				daughterCell->actions[INDEX_GROWTH]->getActualRate());
#if USE_MIGRATION
		actionList->actualizeRate(daughterCell->actions[INDEX_MIGRATION],
				daughterCell->actions[INDEX_MIGRATION]->getActualRate());
#endif
		update_surrounding_expanded_compartment(actionList,
				GetVoronoiCell(daughterCell), voronoiDiagram);

		return GetVoronoiCell(daughterCell);
	} else {
		fprintf(
				stderr,
				"NOT YET IMPLIMENTED in performGrowthAndDivision(): case selected_action->originalCell->state != COMPARTMENT\n");
		exit(0);
	}
}

VoronoiCell* performGrowth(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {

	FILE* fp_dbg;
	//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ," Pre case \n");fclose(fp_dbg );			
	if (selected_action->originalCell->state == COMPARTMENT) {
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Case1: start \n");fclose(fp_dbg );							
		Agent *daughterCell = NULL;
		//printf("%i %lf\n", selected_action->originalCell->growingTumorCellCount, selected_action->rate);


		if (selected_action->originalCell->cellCount
				< selected_action->originalCell->maxCellCount) {
			// into same compartment	
			//fprintf( stderr, "GROW cell %i (%i=%ig+%id/%i)\n", selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->growingTumorCellCount, selected_action->originalCell->dividingTumorCellCount, selected_action->originalCell->maxCellCount);	
			selected_action->originalCell->cellCount++;
			selected_action->originalCell->growingTumorCellCount--;
			selected_action->originalCell->dividingTumorCellCount++;
			selected_action->originalCell->actualize(actionList);
			daughterCell = selected_action->originalCell;
		} else {
			int i;
			// count free neighbors
			int countFreeNeighbors = 0;
			int countSlightlyFreeNeighbors = 0;
			for (i = 0;
					i
							< GetVoronoiCell(selected_action->originalCell)->countNeighborCells;
					i++) {
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])
						== NULL) {
					Agent *newAgent = agentArray->activateAgent();
					newAgent->state = COMPARTMENT;
					newAgent->attach(
							GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
					//fprintf( stderr, "ACTIVATE cell %i\n", newAgent->index);
				}

				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						+ 1
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount)
					countFreeNeighbors++;
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount)
					countSlightlyFreeNeighbors++;
			}
			if (countFreeNeighbors == 0 /*&& countSlightlyFreeNeighbors == 0*/)
				return NULL;
			if (countFreeNeighbors == 0) {
				int whichSlightlyFreeNeighbor = (int) (myRand()
						* (double) countSlightlyFreeNeighbors);
				countSlightlyFreeNeighbors = 0;
				//for( i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells && countFreeNeighbors!=whichFreeNeighbor; i++)
				for (i = 0;
						i
								< GetVoronoiCell(selected_action->originalCell)->countNeighborCells
								&& daughterCell == NULL; i++) {
					if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
							< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount) {
						if (countSlightlyFreeNeighbors
								== whichSlightlyFreeNeighbor) {
							daughterCell =
									GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
							daughterCell->maxCellCount++;
						}
						countSlightlyFreeNeighbors++;
					}
				}

			} else {
				//if( countFreeNeighbors == 0){
				//	return NULL;
				//}

				// chose free neighbor
				int whichFreeNeighbor = (int) (myRand()
						* (double) countFreeNeighbors);
				//fprintf( stderr, "which neighbor %i/%i\n", whichFreeNeighbor, countFreeNeighbors);
				countFreeNeighbors = 0;
				//for( i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells && countFreeNeighbors!=whichFreeNeighbor; i++)
				for (i = 0;
						i
								< GetVoronoiCell(selected_action->originalCell)->countNeighborCells
								&& daughterCell == NULL; i++) {
					if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
							+ 1
							< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount) {
						if (countFreeNeighbors == whichFreeNeighbor)
							daughterCell =
									GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
						countFreeNeighbors++;
					}
				}
			}
			selected_action->originalCell->cellCount--;
			selected_action->originalCell->growingTumorCellCount--;

			daughterCell->cellCount += 2;
			daughterCell->dividingTumorCellCount++;

			selected_action->originalCell->actualize(actionList);
			daughterCell->actualize(actionList);
		}
		myStatistics->count_expanded_cells[0]++;
		myStatistics->count_cell_volume[0]++;
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_GROWTH],
				selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
		actionList->actualizeRate(daughterCell->actions[INDEX_DIVISION],
				daughterCell->actions[INDEX_DIVISION]->getActualRate());
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Case1: end\n");fclose(fp_dbg );
		return GetVoronoiCell(daughterCell);
	}

	
	//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Case1: not chosen \n");fclose(fp_dbg );	
	int i;
	int shifted = FALSE;

	//fprintf(stderr,"GROWTH\n");
	VoronoiCell *newLocation = NULL;
	VoronoiCell *sourceLocation;
	VoronoiCell *occupiedLocation = growCell(voronoiDiagram,selected_action->originalCell, shifted, &sourceLocation);
	if(occupiedLocation != NULL) {
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"ext harmless line\n");fclose(fp_dbg );	
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Init ext nb size:%i\n", occupiedLocation->countExtendedNeighborCells);fclose(fp_dbg );		
	}else{
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"occupied location is null....\n");fclose(fp_dbg );	
	}
	
	int shiftPathLength;
	VoronoiCell **_shiftPath = NULL;

	// SHIFT
	if (shifted) {
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Case2: start \n");fclose(fp_dbg );
		_shiftPath =
				getShiftPath(selected_action->originalCell, &occupiedLocation,
						sourceLocation, shiftPathLength);
		if (_shiftPath[shiftPathLength - 1]->agent != NULL) {
			if (GetAgent( _shiftPath[shiftPathLength-1])->countLocations == 1)
				agentArray->deactivateAgent(
						GetAgent( _shiftPath[shiftPathLength-1]));
			GetAgent( _shiftPath[shiftPathLength-1])->detach(
					_shiftPath[shiftPathLength - 1]);
			myStatistics->count_inner_empty_volume[0]--;
		}
		shiftPath(_shiftPath, shiftPathLength);
		newLocation = _shiftPath[0];
	} else {
		newLocation = occupiedLocation;
	}
#ifdef REFINE
	if( occupiedLocation->coarseParent==NULL) {
		fprintf( stderr, "ERROR in performDivision(): cell: %i has no coarse parent!\n", occupiedLocation->index);
		fprintf( stderr, "growing agent sits on cell %i\n", selected_action->originalCell->location[0]->index);
		printVoronoiDiagram( voronoiDiagram, "errorVD.eps", true);

	}
	if( GetAgent(occupiedLocation->coarseParent)->countFree == GetAgent(occupiedLocation->coarseParent)->maxCellCount) {
		refineNeighborhood( voronoiDiagram, occupiedLocation, actionList, agentArray, pow( GetAgent(occupiedLocation->coarseParent)->maxCellCount,1./DIMENSIONS ));
		refineSurrounding( voronoiDiagram, occupiedLocation, actionList, agentArray, pow( GetAgent(occupiedLocation->coarseParent)->maxCellCount,1./DIMENSIONS ) + 0.5);
		refineNeighborhood( voronoiDiagram, occupiedLocation, actionList, agentArray, pow( GetAgent(occupiedLocation->coarseParent)->maxCellCount,1./DIMENSIONS ) + 0.5);
	}
	GetAgent(occupiedLocation->coarseParent)->countActive++;
	GetAgent(occupiedLocation->coarseParent)->countFree--;
	fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");
#endif
	
	if (newLocation == NULL
			&& (VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD || VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD)) {
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"leaving now....\n");fclose(fp_dbg );	
		selected_action->originalCell->state = NONACTIVE;
		return NULL;
	}
	//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"I'd rather stay.....\n");fclose(fp_dbg );	

	
	if (newLocation == NULL) {
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 123: So it starts\n... ");fclose(fp_dbg );
		//|| ((newLocation->extendedNeighborCellsInitialized == FALSE
		//		|| selected_action->originalCell->location[0]->extendedNeighborCellsInitialized == FALSE) && CellDivisionDepth > 1.)) {
		fprintf(stderr, "Cell %i isn't able to divide!\n",
				selected_action->originalCell->index);
		//fprintf(stderr,"newLocation->extendedNeighborCellsInitialized == %i\n", newLocation->extendedNeighborCellsInitialized);
		fprintf(stderr, "Free Neighbors: %i/%i\n",
				selected_action->originalCell->location[0]->countFreeNeighborCells,
				selected_action->originalCell->location[0]->countNeighborCells);
		fprintf(stderr, "Free Extended Neighbors: %i/%i\n",
				selected_action->originalCell->location[0]->countFreeExtendedNeighborCells,
				selected_action->originalCell->location[0]->countExtendedNeighborCells);
		int countFreeNeighbors = 0;
		for (i = 0;
				i
						< selected_action->originalCell->location[0]->countExtendedNeighborCells;
				i++) {
			if (selected_action->originalCell->location[0]->extendedNeighborhood[i]->isFree())
				countFreeNeighbors++;
		}
		fprintf(stderr, "Found Free Extended Neighbors: %i\n",
				countFreeNeighbors);
		EndTime = (Time < EndTime ? Time : EndTime);
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 123: So it end\n... ");fclose(fp_dbg );
	} else {
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 321: So it starts...\n");fclose(fp_dbg );
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 321: 1...\n");fclose(fp_dbg );		
		if (newLocation->agent != NULL) {
			if (GetAgent( newLocation)->countLocations == 1)
				agentArray->deactivateAgent(GetAgent( newLocation));
			GetAgent( newLocation)->detach(newLocation);
			myStatistics->count_inner_empty_volume[0]--;
		}
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 321: 2...\n");fclose(fp_dbg );		
		selected_action->originalCell->attach(newLocation);
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 321: 3 going into update... ext nb size: %i\n", occupiedLocation->countExtendedNeighborCells);fclose(fp_dbg );		
		update_surrounding_expanded_cell(actionList, occupiedLocation, voronoiDiagram);
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 321: 4 coming back from update\n");fclose(fp_dbg );		
		selected_action->originalCell->actualize(actionList);
		//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Casefile 321: So it ends...\n");fclose(fp_dbg );		
	}
	//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"of screaming\n");fclose(fp_dbg );
	//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"Someone kick me out of my mind... ");fclose(fp_dbg );
	if (shifted) {

		// actualize shifted agents
		//fprintf(stderr,"Actualize shifted agents!\n");
		//fprintf( stderr, "[ %i (%i)", sourceLocation->index, GetAgent(sourceLocation)->index);
		for (i = 0; i < shiftPathLength; i++) {
			GetAgent( _shiftPath[i])->actualize(actionList);
			//	fprintf( stderr, ", %i (%i)", _shiftPath[i]->index, GetAgent(_shiftPath[i])->index);
		}
		//fprintf( stderr, "]\n");
		free(_shiftPath);
	}
	//fp_dbg = fopen( "outdir/voronoi_dbg.txt", "a+");fprintf(fp_dbg ,"I hate these thoughts I can't deny it\n");fclose(fp_dbg );
	selected_action->internalState = 0;
	myStatistics->count_expanded_cells[0]++;
	myStatistics->count_cell_volume[0]++;

	selected_action->originalCell->growingTumorCellCount--;
	selected_action->originalCell->dividingTumorCellCount++;
	//daughterCell->cellCount++;

	return newLocation;

}

/****************************************************************************/

void performDivision(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	if (selected_action->originalCell->state == COMPARTMENT) {

		//fprintf( stderr, "DIVIDE cell %i (%i/%i)\n", selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->maxCellCount);
		selected_action->originalCell->growingTumorCellCount += 2;
		selected_action->originalCell->dividingTumorCellCount--;
		selected_action->originalCell->actualize(actionList);

		myStatistics->count_divided_cells[0]++;
		myStatistics->count_cells[0]++;
		//actionList->actualizeAllRates( voronoiDiagram);
		//actionList->actualizeRate( selected_action, selected_action->getActualRate());
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_GROWTH],
				selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_DIVISION],
				selected_action->originalCell->actions[INDEX_DIVISION]->getActualRate());

		return;
	}


	/*#if SUBCELLULAR_COMPONENTS > 0
	 if( selected_action->originalCell->countLocations == 2*SUBCELLULAR_COMPONENTS){


	 #endif*/
	//sumDivTimes += Time - selected_action->originalCell->timeOfCreation;
	myStatistics->count_divided_cells[0]++;
	myStatistics->count_cells[0]++;
	//global_count_divided_cells++;
	//fprintf( stderr, "DIVISION\n");
#if SUBCELLULAR_COMPONENTS > 0
	{
		// CELL DIVISION

		//fprintf(stderr,"CELL DIVISION\n");
		// STILL WORK IN PROGRESS
		// chose sub-cellular components forming 2nd daughter cell
		int i = selected_action->originalCell->countLocations - 1; //(int)(selected_action->originalCell->countLocations*myRand());
		VoronoiCell * daughterCellLocation =
				selected_action->originalCell->location[i];

		// adapt 1st daughter cell
		selected_action->originalCell->detach(daughterCellLocation);
		//selected_action->originalCell->location[i] = selected_action->originalCell->location[--(selected_action->originalCell->countLocations)];

		// create 2nd daughter cell
		Agent* daughterCell = agentArray->activateAgent();
		daughterCell->attach(daughterCellLocation);
		daughterCell->state = ACTIVE;
		initCellActions(daughterCell);

//		GetAgent(daughterCellLocation->coarseParent)->countActive++;
//		GetAgent(daughterCellLocation->coarseParent)->countFree--;
		/*fprintf( stderr, "daughterCell: active:%i, nonact:%i, free:%i, max:%i\n",
		 GetAgent(daughterCellLocation->coarseParent)->countActive,
		 GetAgent(daughterCellLocation->coarseParent)->countNonactive,
		 GetAgent(daughterCellLocation->coarseParent)->countFree,
		 GetAgent(daughterCellLocation->coarseParent)->maxCellCount);*/

		// Get Distance to closest free Voronoi Cell
		int dist = 0;
		int min_dist = 0;
		for (int d = 0; d < DIMENSIONS; d++)
			min_dist += voronoiDiagram->xN[d];

		if (VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD) {
			VoronoiCell *closestVC;
			if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
				closestVC =
						voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth,
								VoronoiCell::symbolicExtendedNeighborhoodSize,
								VoronoiCell::symbolicExtendedNeighborhood);
			else
				closestVC =
						voronoiDiagram->searchClosestFreeVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth,
								VoronoiCell::symbolicExtendedNeighborhoodSize,
								VoronoiCell::symbolicExtendedNeighborhood);
			if (closestVC)
				min_dist =
						closestVC->getDistanceSquareTo(
								daughterCell->location[0]);
		} else if (VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD) {
			VoronoiCell *closestVC;
			if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
				closestVC =
						voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth);
			else
				closestVC =
						voronoiDiagram->searchClosestFreeVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth);

			if (closestVC)
				min_dist =
						closestVC->getDistanceSquareTo(
								daughterCell->location[0]);
		} else {
			if (daughterCell->location[0]->countFreeNeighborCells)
				for (int i = 0;
						i < daughterCell->location[0]->countNeighborCells;
						i++) {
					dist = 0;
					for (int d = 0; d < DIMENSIONS; d++)
						dist +=
								pow(
										daughterCell->location[0]->position[d]
												- daughterCell->location[0]->neighborCells[i]->position[d],
										2);
					if (min_dist > dist)
						min_dist = dist;
				}

			if (daughterCell->location[0]->countFreeExtendedNeighborCells)
				for (int i = 0;
						i
								< daughterCell->location[0]->countExtendedNeighborCells;
						i++) {
					dist = 0;
					for (int d = 0; d < DIMENSIONS; d++)
						dist +=
								pow(
										daughterCell->location[0]->position[d]
												- daughterCell->location[0]->extendedNeighborhood[i]->position[d],
										2);
					if (min_dist > dist)
						min_dist = dist;
				}
		}

		double prob_reentering_cell_cycle =
				PROB_REENTERING_CELL_CYCLE( sqrt((double)min_dist) * SPATIAL_UNIT);

		// INTOXICATION -> QUIESCENCE
		if( selected_action->originalCell->intoxication > Agent::WASTE_INTOXICATED_CELL_CYCLES/MaxCellDivisionRate)
		//if( selected_action->originalCell->intoxication > Agent::WASTE_INTOXICATED_CELL_CYCLES)
			prob_reentering_cell_cycle=0;

		//if(agentArray->countActiveAgents >100 && selected_action->originalCell->location[0]->ecm <  )
		//	prob_reentering_cell_cycle = 0.;

		//double threshold = 0.003;

		if ( /*agentArray->countActiveAgents <100 ||*/
				( myRand() <= prob_reentering_cell_cycle
				&& VoronoiCell::ECM_THRESHOLD_QUIESCENCE <= daughterCell->location[0]->ecm))
			daughterCell->divide = 1.;
		else
			daughterCell->divide = 0.;

		if ( /*agentArray->countActiveAgents <100 ||*/
				( myRand() <= prob_reentering_cell_cycle
				&& VoronoiCell::ECM_THRESHOLD_QUIESCENCE <= selected_action->originalCell->location[0]->ecm))
			selected_action->originalCell->divide = 1.;
		else
			selected_action->originalCell->divide = 0.;

		//actionList->print();
		//fprintf( stderr, "daughterCell->actualize( actionList);\n ");
		//daughterCell->tumorCellCount++;
		daughterCell->actualize(actionList);
		//actionList->getDepth( actionList->root);
		//actionList->print();
		//fprintf( stderr, "selected_action->originalCell->actualize( actionList);\n ");
		selected_action->originalCell->actualize(actionList);

		selected_action->originalCell->generation++;
		daughterCell->generation = selected_action->originalCell->generation;
		daughterCell->waste = selected_action->originalCell->waste;
		daughterCell->intoxication = selected_action->originalCell->intoxication;

		selected_action->originalCell->growingTumorCellCount++;
		selected_action->originalCell->dividingTumorCellCount--;
		daughterCell->growingTumorCellCount++;
		daughterCell->cellCount++;

		//actionList->print();
		//actionList->getDepth( actionList->root);
		//fprintf( stderr, "after selected_action->originalCell->actualize( actionList);\n ");
		//fprintf(stderr,"CELL DIVISION (a:%i -> a:%i)\n", selected_action->originalCell->index, daughterCell->index);

		//update_surrounding_added_cell( actionList, daughterCell, voronoiDiagram);

		/*fprintf( stderr, "CELL DIVISION: ");
		 fprintf( stderr, "Agent %i ---[division to]---> agent %i { ",
		 selected_action->originalCell->index, selected_action->originalCell->index);
		 int ii;
		 for(ii=0; ii<selected_action->originalCell->countLocations; ii++){
		 fprintf( stderr, "%i |", selected_action->originalCell->location[ii]->index);
		 }
		 fprintf( stderr, "\b} && agent %i { ",
		 daughterCell->index);

		 for(ii=0; ii<daughterCell->countLocations; ii++){
		 fprintf( stderr, "%i |", daughterCell->location[ii]->index);
		 }
		 fprintf( stderr, "\b}\n");*/

		//flag_actualizeProcessRates = TRUE;
	}

#else // SUBCELLULAR_COMPONENTS == 0
	/*VoronoiCell * daughterCellLocation = growCell( selected_action->originalCell);

	 // create 2nd daughter cell
	 Agent* daughterCell = agentArray->activateAgent();
	 daughterCell->attach( daughterCellLocation);
	 update_surrounding_added_cell( actionList, daughterCell, voronoiDiagram);
	 daughterCell->tumorCellCount++;*/

	// I. CELL EXPANSION
	/*VoronoiCell * daughterCellLocation = growCell( selected_action->originalCell);
	 //selected_action->originalCell->attach( daughterCellLocation);
	 //update_surrounding_expanded_cell( actionList, daughterCellLocation, voronoiDiagram);
	 //selected_action->originalCell->actualize( actionList);

	 // II. CELL DIVISION
	 //int i = (int)(SUBCELLULAR_COMPONENTS*myRand());
	 //VoronoiCell * daughterCellLocation = selected_action->originalCell->location[i];

	 // adapt 1st daughter cell
	 //selected_action->originalCell->detach( daughterCellLocation);


	 // create 2nd daughter cell
	 Agent* daughterCell = agentArray->activateAgent();
	 daughterCell->attach( daughterCellLocation);
	 daughterCell->state = ACTIVE;
	 initCellActions( daughterCell);
	 update_surrounding_expanded_cell( actionList, daughterCellLocation, voronoiDiagram);
	 daughterCell->actualize( actionList);
	 selected_action->originalCell->actualize( actionList);

	 daughterCell->tumorCellCount++;
	 selected_action->internalState=0;
	 */

#endif // SUBCELLULAR_COMPONENTS > 0
	/*#if SUBCELLULAR_COMPONENTS > 0
	 }else{
	 // CELL EXPANSION
	 fprintf(stderr,"EXPANSION\n");
	 exit( 0);
	 //VoronoiCell * newLocation = growCell( selected_action->originalCell);
	 selected_action->originalCell->attach( newLocation);
	 update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);
	 //update_expanded_cell( actionList, newLocation, voronoiDiagram);
	 selected_action->originalCell->actualize( actionList);
	 / *fprintf( stderr,"CELL EXPANSION: agent %i ---[expandes to]---> location %i ==> locations {",
	 selected_action->originalCell->index, newLocation->index);
	 for(i=0; i<selected_action->originalCell->countLocations; i++){
	 fprintf( stderr, "%i |", selected_action->originalCell->location[i]->index);
	 }
	 fprintf( stderr, "\b}\n");* /

	 //agentArray->printActiveAgents();
	 }
	 #endif*/

}


