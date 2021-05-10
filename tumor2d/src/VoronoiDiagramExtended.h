#ifndef __VORONOI_DIAGRAM_EXTENDED_H
#define __VORONOI_DIAGRAM_EXTENDED_H

//#include "CellProcesses.h"

// DEFINES
#define FILENAMESIZE 512
#define READBUFFERSIZE 20000

#define MAX_TRIANGLES 6727210
#define MAX_NNS 30

#define INIT_FACES_ARRAY_SIZE 10

#define FACES_ARRAY_EXTENSION_SIZE 10
#define TETRAHEDRA_ARRAY_EXTENSION_SIZE 10
#define POINTS_ARRAY_EXTENSION_SIZE 100

#ifndef DIMENSIONS
#define DIMENSIONS 3
#endif

#define NR_FACE_POINTS        DIMENSIONS
#define NR_TETRAHEDRON_POINTS (DIMENSIONS + 1)

#define TRACK_DELAUNAY_REGION_NEIGHBORHOOD 1

#define LIMITED_NEIGHBOR_DISTANCE	3

#define POINT_POSITION_PRECISION	(1e-9)

#define USE_AGENT

//#define USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD true
//#define USE_DYNAMIC_EXTENDED_NEIGHBORHOOD true

// TYPES
typedef double POSITION_T;
typedef double DISTANCE_T;

// DATA STRUCTURES
//typedef struct _Action Action;
typedef struct _GridPoint GridPoint;
//typedef struct _VoronoiCell VoronoiCell;
class VoronoiCell;
typedef struct _Tetrahedron Tetrahedron;
//typedef struct _VoronoiDiagram VoronoiDiagram;
class VoronoiDiagram;
typedef struct _SphericalCoordinates SphericalCoordinates;
typedef struct _ActionTree ActionTree;
typedef struct _AgentList AgentList;

struct _SphericalCoordinates{
	double phi;		// 0 <= phi   <  2*PI
	double theta;	// 0 <= theta <= PI
};

// voronoi cell type
struct _GridPoint {
	
	// methods
	//static newGridPoint();
	
	// index
	int index;
	
	// spatial position
	double position[DIMENSIONS];

	// neighbor points
	GridPoint **neighborPoints;
	int countNeighborPoints;
};

// voronoi cell type
class VoronoiCell {
public:
	// methodes
	VoronoiCell();
	VoronoiCell( double x, double y, double z);
	~VoronoiCell();

	void actualizeExtendedNeighborhood( VoronoiDiagram *voronoiDiagram, int base_index, int generations, int* nr_kr_neighbors_gen, int radius);
	void actualizeFreeNeighbors();
	void checkFreeNeighbors();
	int isFree();
	int getState();
	double getDistanceTo( double[DIMENSIONS]);
	double getDistanceTo( VoronoiCell* cell);
	double getDistanceSquareTo( VoronoiCell* cell);
	bool isDomainBorder( VoronoiDiagram *vd);

	int index;

	// spatial position
	POSITION_T position[DIMENSIONS];

	// neighborhood
	char neighborCellsInitialized;
	int countNeighborCells;
	int countFreeNeighborCells;
	VoronoiCell **neighborCells;

	// neighborhood
	/*static const char STATIC   = 0;
	static const char SYMBOLIC = 1;
	static const char DYNAMIC  = 2;
	static char extendedNeighborhoodType;
	 */
	static bool USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD;
	static bool USE_DYNAMIC_EXTENDED_NEIGHBORHOOD;
	static bool SHIFT_TO_UNOCCUPIED;

	static bool   USE_LACTATE;
	static double LACTATE_THRESHOLD_QUIESCENCE;
	static double LACTATE_THRESHOLD_DEATH;

	static bool   USE_WASTE;
	static double WASTE_THRESHOLD_QUIESCENCE;
	static double WASTE_THRESHOLD_DEATH;

	static bool   USE_MORPHOGEN;

	static double ECM_THRESHOLD_QUIESCENCE;
	static double ECM_PRODUCTION_RATE;
	static double ECM_DEGRADATION_RATE;


	static double ATP_THRESHOLD_QUIESCENCE;
	static double ATP_THRESHOLD_DEATH;
	static double GLUCOSE_THRESHOLD_QUIESCENCE;
	static double GLUCOSE_THRESHOLD_DEATH;
	static double OXYGEN_THRESHOLD_QUIESCENCE;
	static double OXYGEN_THRESHOLD_DEATH;
	static double GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_QUIESCENCE;
	static double GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_DEATH;

	static int extendedNeighborDepth;

	static int **symbolicExtendedNeighborhood;
	static int *symbolicExtendedNeighborhoodSize;

	char extendedNeighborCellsInitialized;
	int countExtendedNeighborCells;
	int countFreeExtendedNeighborCells;
	VoronoiCell **extendedNeighborhood;

	// meta data
	void*    agent;
	bool	refined;
	int refinement;
	VoronoiCell *coarseParent;

	//Action** actions;		// each cell knows its possible actions
	//int 	 actionsInitialized;
	double	glucose;
	double	oxygen;
	double	growthfactors;
	double  morphogen;
	double	lactate;
	double	waste;
	double  ecm;

	double  dglucose;
	double  doxygen;

	double	CPT11in;
	double	CPT11out;
	SphericalCoordinates sphericalCoordinates;
};

// tetrahedron type
struct _Tetrahedron {

	// methodes
	int addTetrahedronNeighbor( Tetrahedron* neighborTet);

	int index;

	// vertices
	VoronoiCell *vertices[DIMENSIONS+1];

	// neighbor tetrahedra
	int countNeighborTetrahedra;
	Tetrahedron *neighborTetrahedra[DIMENSIONS+1];

	// circumsphere
	char circumsphereInitialized;
	double circumsphere[DIMENSIONS+1];
	double radius;
//	double numericError;
};


// voronoi diagram type
class VoronoiDiagram {
public:
	// methodes
	VoronoiDiagram( char* filename);
	VoronoiDiagram( char* filename, VoronoiDiagram * DummyDiagram);

	void setVoronoiGrid();
	void setDomain();
	void triangulate();
	void sortTetrahedra();
	void setConvexHullNeighborhood();
	void NEWsetExtendedNeighborhood( int radius);
	void setExtendedNeighborhoodAroundVoronoiCell( int radius, int explorationRadius, VoronoiCell *explorationCenter);
	void setExtendedNeighborhoodWithinSphere( int radius, double sphereRadius, double *sphereCenter);
	void setExtendedNeighborhood( double radius);
	void setExtendedNeighborhood( int cells);
	
	static VoronoiCell *searchClosestVoronoiCell( VoronoiCell *explorationCenter, double targetPosition[DIMENSIONS]);
	VoronoiCell *searchClosestFreeVoronoiCell( VoronoiCell *explorationCenter, int explorationRadius);
	VoronoiCell *searchClosestFreeVoronoiCell( VoronoiCell *explorationCenter, int maxDepth, int *symbolicExtendedNeighborhoodSize, int **symbolicExtendedNeighborhood);
	VoronoiCell *searchClosestUnoccupiedVoronoiCellNAIVE( VoronoiCell *explorationCenter);
	VoronoiCell *searchClosestUnoccupiedVoronoiCell( VoronoiCell *explorationCenter, int explorationRadius);
	VoronoiCell *searchClosestUnoccupiedVoronoiCell( VoronoiCell *explorationCenter, int explorationRadius, int *symbolicExtendedNeighborhoodSize, int **symbolicExtendedNeighborhood);
	//static VoronoiCell *searchForVoronoiCell( int countIgnoredPoints, VoronoiCell** ignoredPoints, int countSourcePoints, VoronoiCell** sourcePoints, int &countExploredPoints, VoronoiCell** exploredPoints);
	static VoronoiCell *searchForVoronoiCell( int countIgnoredPoints, VoronoiCell** ignoredPoints, int countSourcePoints, VoronoiCell** sourcePoints, int &countExploredPoints, VoronoiCell** exploredPoints);
	
	static VoronoiDiagram* newVoronoiDiagram();
	static VoronoiDiagram* newVoronoiDiagram( int x, int y, int z);
	static VoronoiDiagram* newVoronoiDiagramFromFile( char* filename);
	void deleteVoronoiDiagram();
	void setPointsOnSphericalSurface( double radius, int N, double variance, double center[DIMENSIONS]);

	void getConvexHull( double thresholdDistance);
	bool PointInsideConvexHull( double *p);
	void setFramePoints();

	int readExtendedNeighborhoodToFile( char* filename, int radius);

	void refine(  VoronoiCell *vc, int scale, ActionTree *actionTree);
	VoronoiCell* coarsen( VoronoiCell *vc, int scale, ActionTree *actionTree, AgentList *agentList);
	bool isHomogen( VoronoiCell *vc, int scale, int &type);

	bool tetrahedronContainsFramePoints( Tetrahedron* tet);
	void writeToFile( char *filename);

	// voronoi cell
	int countVoronoiCells;
	int maxVoronoiCells;
	VoronoiCell **voronoiCells;
	
	// frame
	int countFramePoints;
	VoronoiCell **framePoints;

	// tetrahedra
	int countTetrahedra;
	int maxTetrahedra;
	Tetrahedron **tetrahedra;

	// voronoi grid point
	char voronoiGridSet;
	GridPoint *voronoiGridPoints;

	// link between voronoi grid and voronoi cells
	GridPoint ***voronoiCellToVoronoiGridPoints;
	int        *voronoiCellCountVoronoiGridPoints;

	// faces
	int countFaces;
	int maxFaces;
	VoronoiCell ***faces;
	
	// convex faces
	int countConvexFaces;
	int maxConvexFaces;
	VoronoiCell ***convexFaces;

	// domain
	char domainSet;
	float boundaryThickness;
	float xMin[DIMENSIONS];
	float xMax[DIMENSIONS];
	int xN[DIMENSIONS];
	
	
};





// FUNCTION PROTOTYPES

void getCircumCenter( Tetrahedron *tet, double center[DIMENSIONS]);

// functions for voronoi diagram
VoronoiDiagram* newVoronoiDiagram();

void deleteVoronoiDiagram(  VoronoiDiagram* voronoiDiagram);
VoronoiCell* findNearestCell3D( VoronoiDiagram* voronoiDiagram, VoronoiCell* newVoronoiCell);
void insertVoronoiCell( VoronoiDiagram* voronoiDiagram, VoronoiCell* newVoronoiCell);
void removeVoronoiCell( VoronoiDiagram* voronoiDiagram, VoronoiCell* removedVoronoiCell);
void printPoints( const char * title, VoronoiCell ** allPoints, int api );
void setTetrahedronNeighbors( VoronoiDiagram* voronoiDiagram);
void printTetrahedronNeighbors( Tetrahedron* tetrahedron);
void printTetrahedraNeighbors( VoronoiDiagram* voronoiDiagram);
VoronoiCell * getCentralCell( VoronoiDiagram *voronoiDiagram);
void setDomain( VoronoiDiagram *voronoiDiagram);
double getDistanceOfPointToCircumsphere( VoronoiCell* point, Tetrahedron* tet);
double getCircumsphereRadius( Tetrahedron* tet);
int isElementOf( VoronoiCell **pointList, int nrOfPoints, VoronoiCell *point);

// functions for voronoi cell

// functions for tetrahedron
Tetrahedron* newTetrahedron();

// testing
void checkDelaunayCondition( VoronoiDiagram *voronoiDiagram, VoronoiCell **newCell, int countNewCells);




#endif
