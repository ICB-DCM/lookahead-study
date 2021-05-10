#ifndef __CELL_PROCESSES_H
#define __CELL_PROCESSES_H


// INCLUDES

#include "Mathematix.h"
#include "VoronoiDiagramExtended.h"
#include "Agent.h"


// GLOBAL VARIABLES

extern double MaxCellDivisionRate;
extern double CellDivisionDepth;
extern int    M_div;
extern int    M_gro;
extern int    M_nec;

extern double CellMigrationRate;
extern double CellApoptosisRate;
extern double CellNecrosisRate;
extern double CellLysisRate;
extern double VesselGrowthRate;

extern double alpha_ref;
extern char ReentranceProbabilityFunction;
extern double ReentranceProbabilityReferenceLength;

#define HEAVYSIDE   0
#define EXPONENTIAL 1


#define SCHALLER false
#define JAGIELLA true

// DEFINES

#define USE_ACTION_TREE	TRUE
//#define USE_ACTION_TREE	FALSE

// Cell Death
#define MIN_SPECIFIC_DEATH_RATE 0.00011375
#define MAX_SPECIFIC_DEATH_RATE 0.00583333
#define MONOD_PARAMETER_DEATH 0.001
//#define MONOD_PARAMETER_DEATH 0.02

// Cell Growth
//#define MAX_SPECIFIC_GROWTH_RATE 0.023/0.82842712474619 // m=2;    (m*( pow( 2., 1./m) - 1.)) //0.014583333
//#define MAX_SPECIFIC_GROWTH_RATE (1./22.) // m=10;   (m*( pow( 2., 1./m) - 1.)) //0.014583333
//#define MAX_SPECIFIC_GROWTH_RATE 0.023/0.717734625 // m=10;   (m*( pow( 2., 1./m) - 1.)) //0.014583333
#define MAX_SPECIFIC_GROWTH_RATE (0.023/0.717734625)
//#define MAX_SPECIFIC_GROWTH_RATE (1./24.)//(0.023/0.717734625)
//0.693147180559945 // m=inf
//#define MAX_SPECIFIC_GROWTH_RATE 0.023 // m=10;   (m*( pow( 2., 1./m) - 1.)) //0.014583333

// SCHALLER
//#define MAX_SPECIFIC_GROWTH_RATE (3600/54000) // m=10;   (m*( pow( 2., 1./m) - 1.)) //0.014583333
//#define MAX_SPECIFIC_GROWTH_RATE (0.06666) // m=10;   (m*( pow( 2., 1./m) - 1.)) //0.014583333

#define MONOD_PARAMETER_GROWTH 0.01
//#define MONOD_PARAMETER_GROWTH 0.02

// Monod Kinetics
#define CELL_GROWTH_YIELD_OXYGEN  1.37 * 100000
#define CELL_GROWTH_YIELD_GLUCOSE 1.37 * 100000
#define MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN  0.0056 / 100000
#define MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE 0.0056 / 100000

// RIEGER
#define R_CELL_GROWTH_THRESHOLD		0.01
#define R_MAX_SPECIFIC_GROWTH_RATE	0.1
#define R_CELL_DEATH_THRESHOLD		0.001
#define R_MAX_SPECIFIC_DEATH_RATE 	0.01

//PARAMETERS micrometers, hours, mMol
//micrometers^2.hours^-1 


#define INIT_GROWTHFACTORS 1.
#define INIT_LACTATE 0.

//koo=0.030752, kgo=0.100326, omin=10.171515, omax=22.800151, kgg=0.068049, kog=0.030789, gmin=14.727716, gmax=53.672035

/*

//Glucose Uptake
#define NICK_G_CRITICAL_GLU 0.068049//0.068049
#define NICK_G_CRITICAL_OXY 0.030789
#define NICK_G_MIN (14.727716e-17)
#define NICK_G_MAX (53.672035e-17)


//Oxygen Uptake
#define NICK_O_CRITICAL_OXY 0.030752//0.030752
#define NICK_O_CRITICAL_GLU 0.100326
#define NICK_O_MIN (10.171515e-17)
#define NICK_O_MAX (22.800151e-17)

*/

//Schaller 6 300 000 //1750 * 3600

//Othmer 6552000 //6552000 // 1.82*3.6*1000000


//Schaller 378 000 //105 * 3600
//Othmer 396000  // 396000 // 1.10*3.6*100000



//#define Othmer_OxygenUptake1  77.0055120    //ref: 77.0055120  = 1.0642 * 2.01 * 3.6 * 10
// #define Othmer_OxygenUptake2  43.5621672    //ref: 43.5621672  = 6.0202 * 2.01 * 3.6
// #define Othmer_OxygenUptake3  0.55
// #define Othmer_CriticalOxygenConcentration  0.00464 //4.64 * 0.001

// #define Othmer_GlucoseUptake1  77.0055120   //ref: 77.0055120  = 1.0642 * 2.01 * 3.6 * 10 
// #define Othmer_GlucoseUptake2  12.9372444   //ref: 12.9372444  = 1.7879 * 2.01 * 3.6
// #define Othmer_GlucoseUptake3  0.04
// #define Othmer_CriticalGlucoseConcentration  0.04 //4.0 * 0.01

//#define INITIAL_OXYGEN 0.2
//#define INITIAL_GLUCOSE 25 //25	

#define MILLIMOLAR_TO_MOLAR 0.0001
#define CELL_VOLUME 1
//3e-9
#define CUBECENTIMETER_TO_CELL (1/CELL_VOLUME)
#define CUBICMICROMETER_TO_LITER 1e-15
#define SEC_TO_HOUR (1./3600)
#define TO_MILLI 1000.
#define MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR (TO_MILLI/SEC_TO_HOUR/(AGENT_VOLUME*CUBICMICROMETER_TO_LITER))

#define OXYGEN_OPTIMAL 0.28
#define OXYGEN_THRESHOLD 0.02
#define GLUCOSE_OPTIMAL 5.5
#define GLUCOSE_THRESHOLD 0.06

#define OXYGEN_CONSUMPTION_QUIESCENT (50 * CELL_VOLUME)
#define GLUCOSE_CONSUMPTION_QUIESCENT (80 * CELL_VOLUME)
#define LACTATE_PRODUCTION_QUIESCENT (110 * CELL_VOLUME)

#define GROWTHFACTORS_PRODUCTION (1.0 * CELL_VOLUME)
#define INHIBITORYFACTORS_PRODUCTION (2.0 * CELL_VOLUME)

#define OXYGEN_CONSUMPTION_PROLIFERATING (108 * CELL_VOLUME)
#define GLUCOSE_CONSUMPTION_PROLIFERATING (162 * CELL_VOLUME)
#define LACTATE_PRODUCTION_PROLIFERATING (2 * GLUCOSE_CONSUMPTION_PROLIFERATING)
//#define LACTATE_PRODUCTION_PROLIFERATING (240 * CELL_VOLUME)
//#define LACTATE_PRODUCTION_PROLIFERATING (300 * CELL_VOLUME)

// NEW KINETICS


#define THRESHOLD_QUIESCENCE_GLUCOSE_OXYGEN	0.025//0.00025

#define THRESHOLD_QUIESCENCE_GLUCOSE	0.
#define THRESHOLD_QUIESCENCE_OXYGEN		0.//01
#define THRESHOLD_QUIESCENCE_ATP		0.//400

#define THRESHOLD_NECROSIS_GLUCOSE_OXYGEN	0.025//0.00025
#define THRESHOLD_NECROSIS_GLUCOSE		0.
#define THRESHOLD_NECROSIS_OXYGEN		0.//02
#define THRESHOLD_NECROSIS_ATP			0.//400

#define Schaller_OxygenUptake   (MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR * 20e-18)
//  20.0 * 10^-18 * 1000 * 3600 / 3.10^-12
#define Schaller_GlucoseUptake  (MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR * 95e-18)
//  95.0 * 10^-18 * 1000 * 3600 / 3.10^-12

// JAGIELLA-WEENS
//#define MONOD_PARAMETER_GROWTH_OXYGEN		0.00007
//#define MONOD_PARAMETER_GROWTH_OXYGEN_N		4

//#define MONOD_PARAMETER_GROWTH_GLUCOSE		30.
//#define MONOD_PARAMETER_GROWTH_GLUCOSE_N	3
//ng=10,kg=7,no=4,ko=0.0007
//ng=3,kg=30,no=4,ko=0.00007

#define MONOD_PARAMETER_GROWTH_OXYGEN		0.000007
#define MONOD_PARAMETER_GROWTH_OXYGEN_N		5

#define MONOD_PARAMETER_GROWTH_GLUCOSE		30
#define MONOD_PARAMETER_GROWTH_GLUCOSE_N	3
//ng=3,kg=30,no=6,ko=0.000007

// FREYER
#define THRESHOLD_GLUCOSE	0.02
#define THRESHOLD_OXYGEN	0.06
#define THRESHOLD_LACTATE	8.


// CONSTANTS
#define CUBICMICROMETER_TO_MILLILITER (CUBICMICROMETER_TO_CUBICMETER * CUBICMETER_TO_LITER * LITER_TO_MILLILITER)    //1 000 000 000 000.
#define CUBICMETER_TO_LITER 1000.
#define CUBICMICROMETER_TO_CUBICMETER 0.000000000000000001
#define LITER_TO_MILLILITER 1000.
//#define DIMENSIONS 3

// action types
#define DIVISION                0    // state of Action is a division
#define MIGRATION               1    // state of Action is a diffusion
#define APOPTOSIS               2    // state of Action is a apoptosis
#define NECROSIS                3    // state of Action is a necrosis
#define LYSIS					4    // state of Action is a lyse of necrotic cells in suspension
#define VESSEL_GROWTH           5    // state of Action is a diffusion of a cell in suspension
#define GROWTH                  6    // state of Action is a division
#define GAP						7

#define USE_DIVISION 		TRUE
#define USE_MIGRATION 		TRUE
#define USE_APOPTOSIS 		FALSE
#define USE_NECROSIS 		TRUE
#define USE_LYSIS 			TRUE
#define USE_VESSEL_GROWTH 	FALSE
#define USE_GROWTH 			TRUE

#define INDEX_DIVISION 		0
#define INDEX_APOPTOSIS 	INDEX_DIVISION      + USE_DIVISION	//INDEX_MIGRATION + USE_MIGRATION
#define INDEX_NECROSIS 		INDEX_APOPTOSIS     + USE_APOPTOSIS
#define INDEX_LYSIS 		INDEX_NECROSIS      + USE_NECROSIS
#define INDEX_VESSEL_GROWTH INDEX_LYSIS         + USE_LYSIS
#define INDEX_GROWTH 		INDEX_VESSEL_GROWTH + USE_VESSEL_GROWTH	
#define INDEX_MIGRATION 	INDEX_GROWTH        + USE_GROWTH


// diffusion types
#define TYPE_OF_MIGRATION	FREEMIGRATION

#define FREEMIGRATION		0
#define BORDERMIGRATION		1
#define ENERGYMINIMIZING	2
#define ENERGY			3

#define MONOD_KINETICS 5 // 0:no kinetics, 1:monod kinetics, 2:RIEGER, 3:JAGIELLA-WEENS, 4:Schaller-Meyer-Herrmann, 5:Freyer 

// COMPARTMENTS
#define SUBCELLULAR_COMPONENTS 1
#define MIN_SUBCELLULAR_COMPONENTS 1
#define MAX_SUBCELLULAR_COMPONENTS 2//4
//#define INITIAL_GROWTH_DECISION	FALSE

#define MULTISCALE	0
#define CELLS_PER_SITE NumberOfCellsPerSite//3
#define IMPORT_COARSE_GRAINED_BEHAVIOR 0
#if( MULTISCALE)
	#define DETERMINE_COARSE_GRAINED_BEHAVIOR 0
#else
	#define DETERMINE_COARSE_GRAINED_BEHAVIOR 0
#endif





// DATA TYPES

typedef struct _Action     Action;
typedef struct _ActionList ActionList;

// DATA STRUCTURES

struct _Action{

	double actualizeRate();
	double getActualRate();
	int active();
	void print();

	// action specific attributes
	char    type;
	float rate;    
	char    internalState;      // for increasing the internal state to a maximum M
	int    internalStateM[100];      // for increasing the internal state to a maximum M
	Agent *originalCell; 
	Agent *destinationCell;

	// action tree attributes
	Action*	top;     // next action in probability list
	
	Action*	next;     // next action in probability list
	double	rateSumNext;
	int	sizeNext;

	Action*	prev;     // previous action in probability list
	double	rateSumPrev;
	int	sizePrev;
};



struct _ActionList{

  double sum_of_prob;
  int    length;                  // count of possible actions
  Action *head;    // first element in probability list
};



typedef struct _ActionTree ActionTree;
//typedef struct _ActionTreeElement ActionTreeElement;
struct _ActionTree{

	// methodes
	static ActionTree *newActionTree();
	void   destroyActionTree();
	void   addAction( Action* action);
	void   deleteAction( Action* action);
	void   deleteAllActions( Agent* cell);
	void   actualizeRate( Action* action, double newRate);
	void   actualizeAllRates( VoronoiDiagram* voronoiDiagram);
	Action*selectAction( double *time);
	int    getDepth( Action* action);
	void   print();
	
	// attributes
	double	rateSum;
	int	size;                  // count of possible actions
	Action*	root;    // first element in probability list
};

/*struct _ActionTreeElement{

	double probabilitySum;
	Action *action;
	ActionTreeElement *left;    // first element in probability list
	ActionTreeElement *right;    // first element in probability list
};*/



// PROBABILITY LIST AND ELEMENTS

// memory allocation and initialisation of probability list
ActionList* newActionList();                            // Constructor
void        destroyActionList( ActionList *probList);   // Destructor
void        initCellActions( Agent* cell);

Action* newAction( int type, double rate);

void addDivisionAction( ActionList *actionList, Agent* cell, VoronoiDiagram* voronoiDiagram);
void addNecrosisAction( ActionList *actionList, Agent* cell, VoronoiDiagram* voronoiDiagram);

void addAction( ActionList *probList, Action* cellAction);
void deleteAction( ActionList *probList, Action* cellAction);

void addedVesselAndUpdateSurrounding( ActionTree *p_list, AgentList* agentArray, VoronoiCell* addedVessel,  VoronoiDiagram* voronoiDiagram);
void removeVesselAndUpdateSurrounding( ActionTree *actionList, AgentList* agentArray, VoronoiCell* deletedVessel, VoronoiDiagram* voronoiDiagram );
void update_surrounding_added_cell  ( ActionTree *p_list, VoronoiCell* addedCell,    VoronoiDiagram* voronoiDiagram);
void update_surrounding_expanded_cell(ActionTree *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram );
void update_surrounding_expanded_compartment(ActionTree *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram );
void update_expanded_cell(ActionTree *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram );
void update_surrounding_deleted_cell( ActionTree *p_list, Agent* deletedCell,  VoronoiDiagram* voronoiDiagram);
void update_surrounding_reduced_compartment( ActionTree *actionList, Agent* deleted_cell, VoronoiDiagram* voronoiDiagram );
void update_surrounding_reduced_cell( ActionTree *actionList, VoronoiCell* deleted_cell, VoronoiDiagram* voronoiDiagram );
void update_deleted_cell            ( ActionTree *p_list, Agent* deletedCell,  VoronoiDiagram* voronoiDiagram);
void update_surrounding_removed_cell( ActionTree *actionList, Agent* deleted_cell, VoronoiDiagram* voronoiDiagram );
//void update_surrounding_necrotic_cell(ActionList *p_list, VoronoiCell* necroticCell, VoronoiCell* v_VoronoiCell);

// selecting action randomly but probability weighted 
Action* selectAction( ActionList *actionList, double* time);
void actualizeProcessRates( VoronoiDiagram* voronoiDiagram, ActionList *actionList, double ***Oxygen_Concentration, double ***GrosseFactors_Concentration, double ***Glucose_Concentration);

// 
VoronoiCell* growCell( VoronoiDiagram *voronoiDiagram, Agent* cell, int &shift, VoronoiCell **sourceLocation);
Agent* divideCell( ActionList *actionList, Action* action, VoronoiDiagram* voronoiDiagram);
VoronoiCell* getChildCell( Agent *parent_cell, VoronoiDiagram *v_grid, int *shift_flag, ActionList *actionList);
VoronoiCell** getShiftPath( Agent *expandingAgent, VoronoiCell **end_point, VoronoiCell *sourceLocation, int &pathLength);
void shiftPath( VoronoiCell** path, int path_length);
void initShiftingNeighborhoods( VoronoiDiagram *voronoiDiagram, double radius);


// ToString
const char *actionTypeToString( int type);
void printActionList( ActionList *p_list);

// Kinetics
bool IsQuiescent( Agent * agent);

double GetGrowthRate( Agent * agent);
double GetDivisionRate( Agent * agent);
double GetDeathRate( Agent * agent);
double GetConsumptionRateOxygen( VoronoiCell * cell);
double GetConsumptionRateGlucose( VoronoiCell * cell);
double GetConsumptionRateGrowthFactors( Agent * agent);
double GetDecayRateGrowthFactors( Agent * agent);


double GiveMeTheOxygenRate(VoronoiCell *thecell);
double GiveMeTheGlucoseRate(VoronoiCell *thecell);
double GiveMeTheOxygenRate(VoronoiCell *thecell, float glu, float oxy);
double GiveMeTheGlucoseRate(VoronoiCell *thecell, float glu, float oxy);
double GiveMeTheOxygenRate(Agent *thecell);
double GiveMeTheGlucoseRate(Agent *thecell);
double GiveMeTheATP(VoronoiCell *thecell);
double GiveMeTheATP(Agent *thecell);

#endif

