#include <math.h>
#include <stdio.h>
#include "Mathematix.h"
#include "VesselNetwork.h"
#include "VoronoiDiagramExtended.h"
#include "CellProcesses.h"
//#include "../couple.h"


//double DistanceToInitialCell = 10.;
//int CountVessel = 1;
//double BranchingProbability = 0.00;
//int BranchingLength = 20;


/*****************************************************************************/



void SetInitialVesselNetwork( VoronoiDiagram *voronoiDiagram, AgentList* agentArray, ActionTree *probTree, VoronoiCell *initialCell, double distanceToInitialCell, int countVessel, double branchingProbability, int branchingLength){
	
	int i;

	// place vessels around initial cell
	double angleBetweenVessels = 2. * PI/(double)countVessel;
	double tempAngle = 0.;
/*
	for( i=0; i<countVessel; i++){
		SetVessel( voronoiDiagram, probTree,
		           GetClosestGridPoint( voronoiDiagram, (xminDOMAIN+xmaxDOMAIN)/2. + cos(tempAngle)*INIT_DISTANCE_TUMOR_VESSEL, 
		                                             yminDOMAIN, 
		                                             (zminDOMAIN+zmaxDOMAIN)/2. + sin(tempAngle)*INIT_DISTANCE_TUMOR_VESSEL), 
		           GetClosestGridPoint( voronoiDiagram, (xminDOMAIN+xmaxDOMAIN)/2. + cos(tempAngle)*INIT_DISTANCE_TUMOR_VESSEL, 
		                                             ymaxDOMAIN, 
		                                             (zminDOMAIN+zmaxDOMAIN)/2. + sin(tempAngle)*INIT_DISTANCE_TUMOR_VESSEL),
		           branchingProbability, branchingLength);
		tempAngle += angleBetweenVessels;*/
	
	if( !voronoiDiagram->domainSet)
		voronoiDiagram->setDomain();
		//setDomain( voronoiDiagram);
	
	for( i=0; i<countVessel; i++){
		SetVessel( voronoiDiagram, agentArray, probTree,
		           GetClosestGridPoint( voronoiDiagram, (voronoiDiagram->xMin[0]+voronoiDiagram->xMax[0])/2. + cos(tempAngle)*distanceToInitialCell, 
		                                             voronoiDiagram->xMin[1], 
		                                             (voronoiDiagram->xMin[2]+voronoiDiagram->xMax[2])/2. + sin(tempAngle)*distanceToInitialCell), 
		           GetClosestGridPoint( voronoiDiagram, (voronoiDiagram->xMin[0]+voronoiDiagram->xMax[0])/2. + cos(tempAngle)*distanceToInitialCell, 
		                                             voronoiDiagram->xMax[1], 
		                                             (voronoiDiagram->xMin[2]+voronoiDiagram->xMax[2])/2. + sin(tempAngle)*distanceToInitialCell),
		           branchingProbability, branchingLength);
		tempAngle += angleBetweenVessels;

	}
}
/*****************************************************************************/


void SetRegularInitialVesselNetwork( VoronoiDiagram *voronoiDiagram, AgentList* agentArray, ActionTree *probTree, VoronoiCell *initialCell, double distanceBetweenVessels[DIMENSIONS], double branchingProbability, int branchingLength){
	
	//int i, ii;

	// place vessels around initial cell
	//double angleBetweenVessels = 2. * PI/(double)countVessel;
	//double tempAngle = 0.;
	double vesselOffSet[DIMENSIONS];
	
	for( int i=0; i<DIMENSIONS; i++)
		{
		//this provoque a bug in my angiogenesis function that I am not able to understand.
		//if( i==0)
		//vesselOffSet[i] = (voronoiDiagram->xMax[i]-voronoiDiagram->xMin[i])*0.5+distanceBetweenVessels[i] * myRand();
		//else
		vesselOffSet[i] = distanceBetweenVessels[i] * myRand();
		
		//vesselOffSet[i] = 0.;
		#ifdef __myDEBUG__
		vesselOffSet[i] = 0.;
		#endif
		}
	if( !voronoiDiagram->domainSet)
		voronoiDiagram->setDomain();
	

	
	for( double x = voronoiDiagram->xMin[0]+vesselOffSet[0]; x < voronoiDiagram->xMax[0]; x+=distanceBetweenVessels[0])
	for( double y = voronoiDiagram->xMin[1]+vesselOffSet[1]; y < voronoiDiagram->xMax[1]; y+=distanceBetweenVessels[1])
	for( double z = voronoiDiagram->xMin[2]+vesselOffSet[2]; z < voronoiDiagram->xMax[2]; z+=distanceBetweenVessels[2]){
		SetVessel( voronoiDiagram, agentArray, probTree,
		           GetClosestGridPoint( voronoiDiagram, x,                           y, z), 
		           GetClosestGridPoint( voronoiDiagram, x+distanceBetweenVessels[0], y, z),
		           branchingProbability, branchingLength);
		SetVessel( voronoiDiagram, agentArray, probTree,
		           GetClosestGridPoint( voronoiDiagram, x, y                          , z), 
		           GetClosestGridPoint( voronoiDiagram, x, y+distanceBetweenVessels[1], z),
		           branchingProbability, branchingLength);
		SetVessel( voronoiDiagram, agentArray, probTree,
		           GetClosestGridPoint( voronoiDiagram, x, y, z                          ), 
		           GetClosestGridPoint( voronoiDiagram, x, y, z+distanceBetweenVessels[2]),
		           branchingProbability, branchingLength);

	}

}
/*****************************************************************************/


VoronoiCell *GetClosestGridPoint( VoronoiDiagram *voronoiDiagram, double x, double y, double z){

	int i;
	double tempDist,
	       minDist = sqrt( pow( voronoiDiagram->voronoiCells[0]->position[0] - x, 2.) + pow( voronoiDiagram->voronoiCells[0]->position[1] - y, 2.) + pow( voronoiDiagram->voronoiCells[0]->position[2] - z, 2.)); 
        VoronoiCell *startPoint,
	     *nextPoint = voronoiDiagram->voronoiCells[0];

	do{
		startPoint = nextPoint;
		for( i=0; i<startPoint->countNeighborCells; i++){
			tempDist = sqrt( pow( startPoint->neighborCells[i]->position[0] - x, 2.) + pow( startPoint->neighborCells[i]->position[1] - y, 2.) + pow( startPoint->neighborCells[i]->position[2] - z, 2.));             
			if( minDist > tempDist){
				minDist = tempDist;
 				nextPoint = startPoint->neighborCells[i];
			}
		}
	}while( startPoint != nextPoint);

	return startPoint;
}
/*****************************************************************************/



void SetVessel( VoronoiDiagram *voronoiDiagram, AgentList* agentArray, ActionTree *probTree, VoronoiCell *startCell, VoronoiCell *endCell, double branchingProbability, int branchingLength){

		// set vessel cells
		//fprintf( stderr, "INFO: set vessel cells\n");
		int i;
		VoronoiCell *nextCell = startCell;
		double minDist, tempDist;
		do{
			// Main Vessel
			startCell = nextCell;
			//fprintf( stderr, "INFO: Turn cell %i into a vessel segment!\n", startCell->nr);
			if( startCell->getState() == FREE)
				addedVesselAndUpdateSurrounding( probTree, agentArray, startCell, voronoiDiagram);
			//else
				//fprintf( stderr, "INFO: cell %i is already declared as %s!\n", startCell->index, cellTypeToString(startCell->getState()));
			
			minDist = sqrt( pow( startCell->position[0] - endCell->position[0], 2.) + pow( startCell->position[1] - endCell->position[1], 2.) + pow( startCell->position[2] - endCell->position[2], 2.));
			// look for best neighbor
			for( i=0; i<startCell->countNeighborCells; i++){
				//if( startCell->neighborCells[i]->getState() != VESSEL)
				{
					tempDist = sqrt( pow( startCell->neighborCells[i]->position[0] - endCell->position[0], 2.) 
					               + pow( startCell->neighborCells[i]->position[1] - endCell->position[1], 2.) 
					               + pow( startCell->neighborCells[i]->position[2] - endCell->position[2], 2.)); 
					double neighborDist = sqrt( pow( startCell->neighborCells[i]->position[0] - startCell->position[0], 2.) 
						                  + pow( startCell->neighborCells[i]->position[1] - startCell->position[1], 2.) 
						                  + pow( startCell->neighborCells[i]->position[2] - startCell->position[2], 2.)); 
					            
					if( minDist > tempDist && neighborDist < 2.){
						minDist = tempDist;
						nextCell = startCell->neighborCells[i];
					}
				}
			}

			// Branching / Sub-vessels
			if(myRand() < branchingProbability){
				//printf("INFO:  Branching!\n");
				MakeRandomBranch( voronoiDiagram, agentArray, probTree, startCell, branchingLength);
				//printf("INFO:  Branching ended!\n");
			}
		}while( nextCell != startCell ); 
}
/*****************************************************************************/


void MakeRandomBranch( VoronoiDiagram *voronoiDiagram, AgentList* agentArray, ActionTree *probTree, VoronoiCell *startCell, int branchingLength){

	int i, ii;

	if( branchingLength != 0){
		// chose next vessel segment
		int countPossibleCandidates = 0;
		VoronoiCell *possibleCandidates[MAX_NNS];
		for( i=0; i<startCell->countNeighborCells; i++){ // for each candidate
			VoronoiCell *nextCell = startCell->neighborCells[i];	
			if( nextCell->getState()!=VESSEL){ // check that is not already a vessel
				int valid = 1;
				for( ii=0; ii<nextCell->countNeighborCells; ii++){ // for all neighbors of candidate
					if( nextCell->neighborCells[ii]->getState()==VESSEL && nextCell->neighborCells[ii]!=startCell) // check that no contact to vessels
						valid = 0;
				}
				if( valid){
					possibleCandidates[ countPossibleCandidates] = nextCell;
					countPossibleCandidates++;
				}
			}
		}

		i = (int)(myRand() * (double)countPossibleCandidates);
		if( countPossibleCandidates!=0){
		//double probToBeChoosen = 1./(double)countPossibleCandidates;
		//for( i=0; i<countPossibleCandidates; i++){
			//if( myRand()<=probToBeChoosen){
				// next segment choosen
				VoronoiCell *nextCell = possibleCandidates[ i];
				//fprintf( stderr, "INFO: Turn cell %i into a vessel segment!\n", nextCell->nr);
				addedVesselAndUpdateSurrounding( probTree, agentArray, nextCell, voronoiDiagram);
				
				//fprintf( stderr, "INFO: Turn cell %i into a vessel segment!\n", startCell->nr);
				//printf( "INFO: Subbranching: %lf < %lf < %lf, %lf < %lf < %lf, %lf < %lf < %lf\n", 
				//        xminDOMAIN, nextCell->position[0], xmaxDOMAIN, yminDOMAIN, nextCell->position[1], ymaxDOMAIN, zminDOMAIN, nextCell->position[2], zmaxDOMAIN);
				double BoundaryLayer = 1.;
				if( nextCell->position[0] > voronoiDiagram->xMin[0]+BoundaryLayer && nextCell->position[0] < voronoiDiagram->xMax[0]-BoundaryLayer &&
				    nextCell->position[1] > voronoiDiagram->xMin[1]+BoundaryLayer && nextCell->position[1] < voronoiDiagram->xMax[1]-BoundaryLayer &&
				    nextCell->position[2] > voronoiDiagram->xMin[2]+BoundaryLayer && nextCell->position[2] < voronoiDiagram->xMax[2]-BoundaryLayer){

					// continue this branche
					//printf("INFO: Subbranching\n");
					MakeRandomBranch( voronoiDiagram, agentArray, probTree, nextCell, branchingLength-1);
				
				}
				return;
		//	}
		}
	}

	return;		
}

