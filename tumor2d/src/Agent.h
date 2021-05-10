#ifndef __AGENT_H
#define __AGENT_H


// DEFINES

// cell states
#define FREE             0    
#define ACTIVE           1
#define NONACTIVE        2
#define NECROTIC         3
#define VESSEL           4
#define COMPARTMENT      5
#define QUIESCENT        6
#define FRAME            7



// INCLUDES

#include <stdlib.h>
#include <math.h>
#include "VoronoiDiagramExtended.h"
//#inclu4de "CellProcesses.h"


// DEFINES

//#define AGENT_VOLUME		2000 
//#define AGENT_DIAMETER	15.6 // volume = 2000

// EMT6/Ro cells [Freyer, Sutherland]
//#define AGENT_VOLUME	(1200) // mym^3
//#define AGENT_DIAMETER	13.1844153 // volume = 1200

// SK-MES-1 cells
//#define AGENT_DIAMETER	15.1 // volume = 1200
#define AGENT_DIAMETER		16.8// volume = 1200
#define NUCLEUS_DIAMETER	10// volume = 1200
#define AGENT_VOLUME	(1800) // mym^3
//#define AGENT_VOLUME	(AGENT_DIAMETER*AGENT_DIAMETER*AGENT_DIAMETER * M_PI/6.) // mym^3

//#define AGENT_VOLUME	(1500) // mym^3
//#define AGENT_DIAMETER	14.3 // volume = 1500

//#define AGENT_VOLUME	(3000) // mym^3
//#define AGENT_DIAMETER	18.  // volume = 3000
#define DIAMETER_TO_CUBE_LENGTH( diameter) (diameter * 0.805995977008235)
#define CUBE_LENGTH_TO_DIAMETER( length)   (length / 0.805995977008235)
#define SPATIAL_UNIT	DIAMETER_TO_CUBE_LENGTH( AGENT_DIAMETER)

#define GetNeighborAgent( agent, j) \
		( GetAgent( GetVoronoiCell( agent)->neighborCells[j]) )

#define GetNeighborCell( agent, j) \
		( GetVoronoiCell( agent)->neighborCells[j] )

#define GetExtendedNeighborCell( agent, j) \
		( GetVoronoiCell( agent)->extendedNeighborhood[j] )


#define GetExtendedNeighborAgent( agent, j) \
		( GetAgent( GetVoronoiCell( agent)->extendedNeighborhood[j]) )


#define GetVoronoiCell( agent) \
		( ( VoronoiCell*) agent->location[0] )


#define GetIndex( agent) \
		( agent->index )


#define GetPosition( agent) \
		( (( VoronoiCell*) agent->location[0])->position)
//		( agent->position)

#define AddLocation( agent, latticeSite) \
		{ agent->location = (VoronoiCell**) realloc( agent->location, ++(agent->countLocations) * sizeof(VoronoiCell*)); \
		  agent->location[agent->countLocations-1] = latticeSite; }

#define GetAgent( voronoiCell) \
		( ( Agent*) voronoiCell->agent )


#define AddAgent( voronoiCell) \
		{ voronoiCell->data = (Agent*) malloc( sizeof(Agent)); }


#define RemoveAgent( voronoiCell) \
		( voronoiCell->data = NULL )


#define DestroyAgent( voronoiCell) \
		( free( voronoiCell->data) )


// CLASSES

typedef struct _Action     Action;
typedef struct _ActionList     ActionList;
typedef struct _ActionTree     ActionTree;
//typedef struct _Agent     Agent;
typedef struct _AgentList AgentList;

// Cell Data datatype

class Agent{
public:
	//Glucose Uptake
	static double NICK_G_CRITICAL_GLU;
	static double NICK_G_CRITICAL_OXY;
	static double NICK_G_MIN;
	static double NICK_G_MAX;


	//Oxygen Uptake
	static double NICK_O_CRITICAL_OXY;
	static double NICK_O_CRITICAL_GLU;
	static double NICK_O_MIN;
	static double NICK_O_MAX;

	static bool USE_GRAVITY;

	static float WASTE_UPTAKE;
	static float WASTE_THRESHOLD_QUIESCENCE;
	//static float WASTE_THRESHOLD_QUIESCENCE_INTRACELLULAR;
	static float WASTE_THRESHOLD_SLOWED_GROWTH;
	static int   WASTE_INTOXICATED_CELL_CYCLES;

	static float ReentranceRate;


	Agent();
	Agent( VoronoiCell * correspondingVoronoiCell);

	void initAgent();
	void actualize( ActionList *actionList);
	void actualize( ActionTree *actionTree);
	void validate();
	void attach( VoronoiCell * correspondingVoronoiCell);
	void detach();
	void detach( VoronoiCell * correspondingVoronoiCell);
	void print();

	int isConnected();

	int isFree();
	int isGrowing();
	int isDividing();
	int isDying();
	int isLysing();
	
	int countFreeNeighbors();

	double getOxygen();
	double getGlucose();

	
	// index
	int index;
	
	// location
	VoronoiCell** location;
	int     countLocations;

	// cell state
	char	state;				// is there a cell 0 - FREE, 1 - ACTIVE, 2 - NONACTIVE, 3 - NECROTIC, 4 - VESSEL

	// multiscale
	int	maxCellCount;
	int	cellCount;
	int	growingTumorCellCount;
	int	dividingTumorCellCount;
	int	necroticCellCount;

	int generation;

	// histogram
	int countActive;
	int countNonactive;
	int countFree;


	// cell kinetics
	//int	mx,my,mz;         // integer position in the different diffusion-computation meshes
	//void *	pTissueData;
	//double	glucose;
	//double	oxygen;
	//double	growthfactors;

	// neighborships
	//int	countDirectNeighbors;		// number of free neighbours
	//int	countDirectFreeNeighbors;	// number of free neighbours
	//int	countReachableNeighbors;	// number of points within a radius k_r
	//int	countReachableFreeNeighbors;
	//Agent** reachableNeighbors;	// list with points within a radius k_r

	// cell dynamics
	Action** actions;		// each cell knows its possible actions
	char 	actionsInitialized;
	//double	timeOfCreation;

	float divide;

	float waste;
	float intoxication;
};


struct _AgentList{
	
	// methodes
	static AgentList* newAgentList( int countNewAgents);
	Agent* activateAgent();
	//void deactivateAgent( Agent* agent);
	void deactivateAgent( Agent* agent);
	void deleteAgentList();
	void print();
	void printActiveAgents();
	
	// attributes
	Agent ** agents;
	int countAgents;
	int countActiveAgents;	
};

// GLOBAL VARIABLES

extern int CountCellsPerVoronoiCell;



// FUNCTION PROTOTYPES


//AgentArray * newAgentArray( VoronoiDiagram*);
Agent ** initAgents( VoronoiDiagram* voronoiDiagram);
const char *cellTypeToString( int type);

#endif

