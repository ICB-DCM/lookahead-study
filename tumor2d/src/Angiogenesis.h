#ifndef __ANGIOGENESIS_H
#define __ANGIOGENESIS_H

#include <cstdlib>
#include <iostream>

#include "Agent.h"
#include "CellProcesses.h"
#include "Vessel_Unit.h"
#include "Vessel_Graph.h"

//#define COMMENTS
//practical stuffs
#define SQR2 1.4142135623
#define SQR3 1.7320508075
#define LN2  0.693147181
#define NONINDOMAIN -7
#define Pi 3.141592653589793238462643383279502
#define MICRO_EPSILON (1e-37)

//All Parameters


//Oxygen stuff: Inlet is the quantity of oxygen at the beginning of a vessel
#define BOUNDARYOXYGEN 0.0
#define OxyEPSILON 0.5
#define m_oxygen 0.0056
#define Y_oxygen 1.37
#define k_oxygen 0.2
#define mu_max_oxygen 0.047

//glucose

#define m_glucose 0.0056
#define Y_glucose 1.37
#define k_glucose 0.2
#define mu_max_glucose 0.047
//Gf

//Vessel stuff
#define INIT_VISCOSITY 0.1
#define INIT_RADIUS 10 //10
#define MAXRADIUS (3.5 * INIT_RADIUS)

//Pressure stuff
#define InletPRESSURE 	100.0
#define OutletPRESSURE 	-100.0
#define WallPRESSURE 0.0
#define PRESSURE_EPSILON 0.1


//Parameters for the vessels to change
#define MINOXYGEN 0.01 //0.01
#define MINSHEAR 0.01

#define PROBASPROUT 0.90
#define PROBAFUSION 0.20
#define PROBADILATE 0.90
#define PROBAREGRESSION 0.9
#define PROBACOLLAPSE 0.9


#define EC_PROLIFERATION_TIME 40
#define MAXSPROUT 8 // in terms of voronoicells
#define SURROUNDING_PERCENTAGE 0.8

extern double ANGIOGENESIS_GROWTHFACTOR_THRESHOLD;

extern double xminDOMAIN;
extern double yminDOMAIN;
extern double zminDOMAIN;
extern double xmaxDOMAIN;
extern double ymaxDOMAIN;
extern double zmaxDOMAIN;
extern double BoundaryLayer;


using namespace std;



class Vessel_Graph;
class Vessel_Unit;

extern bool DoWeKnowEachOther(Vessel_Graph *N, Vessel_Graph *M);
extern bool ProbaToSprout(double dt);
extern bool ProbaToCollapse(double dt);
extern bool ProbaToRegress();
extern bool SurroundingNode(Vessel_Graph * Vg);
extern bool SurroundingSlice(Vessel_Unit * Vu);

class Angiogenesis
	{
	public:
	Angiogenesis();
	~Angiogenesis();
	

	void Initialize(AgentList *agentList);
	bool Update_Network(double dt, ActionTree *p_list, VoronoiDiagram* voronoiDiagram);	
	void GiveMeTheBlood();
	void Update_Graph(Agent *agentmoved);



	AgentList* angioAgentList;
//	int VoronoiCell_Number;
//	VoronoiCell **VoronoiCell_List;

	Vessel_Graph **Node;
	int nbrNodes; //the number of nodes
	int maxNodes;

	Vessel_Unit** Slice;
	int nbrSlices; //the number of slices
	int maxSlices;
	


	private:


	void KillSlice(Vessel_Unit *X);
	void MakeOneEdge(Vessel_Graph *a, Vessel_Graph *b);
	void KillNode(Vessel_Graph *N);
	Vessel_Graph* MakeOneNode();
	bool DoTheSprout(Vessel_Graph *N, ActionTree *p_list, VoronoiDiagram* voronoiDiagram);
	bool DoesExistTheSlice(Vessel_Graph* A, Vessel_Graph* B);
	void DoTheFusion(Vessel_Graph *N, Vessel_Graph *M);
	bool DoTheDilation(Vessel_Graph *N);
	bool DoTheVesselCollapse(double dt);
	bool KillConnection(Vessel_Unit *Vu);
	void Holocaust();
	Vessel_Graph* FindMyNode(Agent *A);
	
	void GiveMeThePressure(double epsilon);
	void GiveMeTheFlow();
	void GiveMeTheSpeed();
	void GiveMeTheShearStress();

	void VasculatureCreation();

	//debug
	};
#endif
