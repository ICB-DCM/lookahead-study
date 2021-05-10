#ifndef __VESSEL_NETWORK_H
#define __VESSEL_NETWORK_H

// Vascularization
#define GROWTHFACTOR_THRESHOLD 0.1


//#include "tools_coarse.h"
#include "VoronoiDiagramExtended.h"
#include "CellProcesses.h"


void SetInitialVesselNetwork( VoronoiDiagram *voronoiDiagram3D, AgentList* agentArray, ActionTree *probTree, VoronoiCell *initialCell, double distanceToInitialCell, int countVessel, double branchingProbability, int branchingLength);
void SetRegularInitialVesselNetwork( VoronoiDiagram *voronoiDiagram, AgentList* agentArray, ActionTree *probTree, VoronoiCell *initialCell, double distanceBetweenVessels[DIMENSIONS], double branchingProbability, int branchingLength);

void SetVessel( VoronoiDiagram *voronoiDiagram3D, AgentList* agentArray, ActionTree *probTree, VoronoiCell *startCell, VoronoiCell *endCell, double branchingProbability, int branchingLength);
void MakeRandomBranch( VoronoiDiagram *voronoiDiagram3D, AgentList* agentArray, ActionTree *probTree, VoronoiCell *startCell, int branchingLength);
VoronoiCell * GetClosestGridPoint( VoronoiDiagram *voronoiDiagram3D, double x, double y, double z);


#endif

