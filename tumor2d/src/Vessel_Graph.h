#ifndef __VESSEL_GRAPH_H
#define __VESSEL_GRAPH_H

//INCLUDES


#include "Vessel_Unit.h"
#include "VoronoiDiagramExtended.h"
#include "Angiogenesis.h"
#include "Agent.h"
#include <math.h>

//DEFINES
//GLOBAL VARIABLES

class Vessel_Unit;
//CLASSES
class Vessel_Graph
	{
	public:
	Vessel_Graph();
	
	int index;
	double x,y,z;
	//VoronoiCell *p_voronoiCell;
	Agent *p_agent;
	//double growth_factor;

	//to compute the pressure with a trilinear function
	double d_max, d_min, d;


	double pressure;
	double viscosity;
	double radius;
	double outflow; //the total flow is zero so outflow = -inflow. We take the outflow positive.
	bool circulation;
	bool InDomain;
	int CountBranches;
	
	Vessel_Graph **branch; //branch[i] is the pointer of the node
	Vessel_Unit  **edge;   //edge[i] is the pointer of the edge between the node and its branch[i]
	void PressureComputation();
	bool DefineInDomain();
	int FindBranchNumber( Vessel_Graph *b);
	}; 

//FUNCTION PROTOTYPES

#endif
