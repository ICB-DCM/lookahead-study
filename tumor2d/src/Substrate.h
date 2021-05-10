#ifndef __SUBSTRATE_H
#define __SUBSTRATE_H



#define TIME_STEP .1 //.04

//#define Oxygen_Diffusion (6300000) //(6300000)
//#define Glucose_Diffusion (378000) // (378000)
extern double Oxygen_Diffusion; //(6300000)
extern double Glucose_Diffusion; // (378000)
//#define Lactate_Diffusion (212400)
extern double Lactate_Diffusion;// (212400)
extern double Waste_Diffusion;// (212400)
//#define Lactate_Diffusion 8484.00
#define GrowthFactors_Diffusion 100

//#define QUALITY_CONTROL
//#define __myOMP__
//#define __TestSpeed__

#ifdef __WithGUI__
#include <gui/window.h>
class Window;
#endif



#ifdef __MultiMade__
#include <time.h>
#endif

#ifdef __myOMP__
#include <omp.h>
#endif

#include <ctime>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Montecarlo.h"

//#include "Agent.h"
#include "CellProcesses.h"
#include "Molecule.h"
#include "Angiogenesis.h"

//COMMENTS TO USE THE CLASS
//Substrate(..,..,)
//DomainSize is an edge of a cube means that the volume of the cube is DomainSize to the power 3
//dx is the size of the small cube where the agent will placed randomy in the vorono� grid.
//type = 0 means that the domain is a cube type = 1 means that the domain is a torus
//allagents is an array on pointers of agent
//Init_Oxy is the initial concentration of oxygen or at the boundary (vessels, free sites)
//Init_Glu the same but for the glucose


	//possibilities for each molecule:
		//CONSUMPTION
		//PRODUCTION
		//DIFFUSION
		//DECAY
		//REFILL
		//BOUNDARY CONDITIONS


using namespace std;


class DataRecord
	{
	public:
	DataRecord(){}

	void Init()
		{
		for(int i = 0; i < MaxRecord; i++)
			{
			time[i] = 0.;
			living_cells[i] = 0;
			vessels[i] = 0;
			radius_of_gyration[i] = 0;
			calls[i] = 0;
			}
		}
	
	static const int MaxRecord = 3000;
	static const double TimeFrequency;

	double time[MaxRecord];
	int living_cells[MaxRecord];
	int vessels[MaxRecord];
	double radius_of_gyration[MaxRecord];
	int calls[MaxRecord];

	void addvalue(double t, int a, int b, double c)
		{
		int mm = (int)(t * TimeFrequency);
		if( mm < MaxRecord )
			{
			calls[mm]++;
			time[mm] = t;
			vessels[mm] = ((vessels[mm] * (calls[mm] - 1)) + a ) / calls[mm];
			living_cells[mm] = ((living_cells[mm] * (calls[mm] - 1)) + b ) / calls[mm];
			radius_of_gyration[mm] = ((radius_of_gyration[mm] * (calls[mm] - 1)) + c ) / calls[mm];
			}


		}	
	};

class Timer
	{
	private:
	clock_t start_clock;
	clock_t mytimer;
	
	public:
	Timer(): mytimer(0),calls(0) {}
	int calls;
	void Start() {start_clock = clock(); calls++;}
	void End(){mytimer += (clock() - start_clock);}
	
	double Mean(){return ((double)(mytimer / CLOCKS_PER_SEC)/(double)calls) ;}	
	double Total(){return ((double)(mytimer/ CLOCKS_PER_SEC)) ;}	
	};
class Substrate
	{
	public:

	#ifndef __WithGUI__
	Substrate(int* sizeofaxes, int type, VoronoiDiagram *voronoiDiagram, AgentList *agentList, double Init_Oxy, double Init_Glu, double Work_Parameter, const char *thedirname, char*, double customTimestep);
	#else
	bool gui_is_defined;
	Window *gui;
	Substrate(int* sizeofaxes, int type, VoronoiDiagram *voronoiDiagram, AgentList *agentList, double Init_Oxy, double Init_Glu, double Work_Parameter, const char *thedirname, char *boundaryCondition, double customTimestep, Window *thegui);
	void SwitchOnVisualization();
	void SwitchOffVisualization();
	bool pause;
	#endif

	Angiogenesis vasculature;

	~Substrate();
	void Initialization(int mode);
	void RebuildMatrix();
	double getdata(int WhatData);

	int axe_size[3];


	AgentList *myAgentList;
	VoronoiDiagram* myVoronoiDiagram; 
	VoronoiCell **voronoiCells;

	double spacestep;
	double timestep;

	//Used by the gui
	double Initial_Oxygen;
	double Initial_Glucose;
	double Initial_Growthfactors;
	double Initial_Lactate;

	double oxygen_dirichlet_value;
	double glucose_dirichlet_value;
	double growthfactors_dirichlet_value;
	double lactate_dirichlet_value;
	
	double TotalSimulationTime;

	void Refill_Oxygen(int type);
	void Refill_Glucose(int type);
	void Produce_Growthfactors();

	private:
	const char *outputpath;
	
	Timer timer;

	static Molecule Oxygen;
	static Molecule Glucose;
	static Molecule Growthfactors;
	static Molecule Lactate;
	
	void Initialize_Oxygen();
	void Initialize_Glucose();
	void Initialize_Growthfactors();
	void Initialize_Lactate();
	
	VoronoiCell ****TheCubes;
	
	double TimeSum;


	//Refill functions
	//void Refill_Oxygen(int type);
	//void Refill_Glucose(int type);
	void Refill_Lactate();
	
	//Production functions
	void Produce_Lactate();
	//void Produce_Growthfactors();

	
	//Gnuplot
#ifdef __With_Outputs__
	static const int dimensionview = 1;
	static const int space_output_frequency = 1;
	static const double time_output_frequency = 0.1;
	
	DataRecord myRecord;
	
	int y_viewfixpoint;
	int z_viewfixpoint;

	int max_output_number;
	int simulation_number;
	int outputnumber;
	int Output(int Number);
	void MakeGnuplotScript(int NumberOfFiles);
	int SimpleOutput();
	
	int getdata_angiogenesis(int WhatData);


#endif

	//Quality Control	
#ifdef QUALITY_CONTROL
	void QualityControl(int mode);
	void Instationnaire_q(int mode);
	void Instationnaire_qu(int mode);
	void Parabole(int mode);
	void DirichletZero(int mode);
	void TimeStepControl(int mode);
	void Output_QualityControl(double *tab, int factor,char *filename);
	static double FixRate(VoronoiCell *thecell);
	static double fixrate;
	Molecule mymolecule;
	double alpha_mymolecule;
	double mymolecule_diffusion;
	double mymolecule_dirichletvalue;
	static const int quality_output_frequency = 1;
	double *centralCellValue;
	static const double ControlDuration = 1.;

#endif

	}; 


#endif


