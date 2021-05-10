#ifndef __VESSEL_UNIT_H
#define __VESSEL_UNIT_H
//INCLUDES

#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "Vessel_Graph.h"

//DEFINES
//GLOBAL VARIABLES
//CLASSES

class Vessel_Graph;

class Vessel_Unit
	{
	public:
	Vessel_Unit()
		{
		radius = 1;
		velocity = 1;
		viscosity = 1;
		flow = 1;
		flowproportion = 1;
		shear = 1;
		circulation = true;
		}
		
	int index;
	bool circulation;
	Vessel_Graph *node1;
	Vessel_Graph *node2;//the nodes in the graph which are connected by this edge.
	double theta,phi,length;
	
	double radius;
	double velocity;
	double viscosity; // could be a fonction.
	double flow; //The positive blood flow (in the direction of the circulation)
	double flowproportion; //the percentage of the flow from a node to the other.(depends of the direction of the flow)
	double shear; //The shear stress on slice's wall
	//double growth_factor;
	//double oxygen;
	
	//Carefull, the computation of values may depends of the other values, it has to be compute in a certain order.
	void Compute_length();
	void Compute_theta();
	void Compute_phi();
	void Compute_radius();
	void Compute_viscosity();
	void Compute_flow();
	void Compute_velocity();
	void Compute_shearstress();
	};
//FUNCTION PROTOTYPES	


	
#endif

