#ifndef __MOLECULE_H
#define __MOLECULE_H


#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include "CellProcesses.h"
#include "Agent.h"

#define PURE_DIFFUSION 0
#define CONSUMPTION_DIFFUSION 1
#define CONSUMPTION_DIFFUSION_TEST_1 2
#define CONSUMPTION_DIFFUSION_TEST_2 3


#define BORDERAREBORDER		0
#define FREEAREBORDER 		1
#define VESSELAREBORDER		2
#define NECROTICAREBORDER	3
using namespace std;

class Molecule
	{
	public:
	
	Molecule()
		{
		modulo = 0;
		}
	~Molecule();
	
	static void SetDomain(int x_size,int y_size,int z_size, double thetimestep, VoronoiCell ****thedomain);
	static void SetAgents( AgentList *theAgentArray);
	static void ChangeTimeStep(double ts){dt = ts;}
	void SetMethod(int methode, int refill, string bordercondition);
	void SetMemory();
	void SetMatrix(double parameter);
	void SetValuePointer(int mx, int my, int mz,double *thevalueaddress); 

	void (*GiveMe)(Molecule *Mo);
	double (*Rate)(VoronoiCell *thecell);

	void (*BuildMatrix)(int axe, int posx, int posy, int posz, Molecule *Mo, int bc);

	private:
	static int size[3];
	static double dt;

	int modulo;
	double alpha;
		
	double *diagonal[3];
	double *upperdiagonal[3];
	double *lowerdiagonal[3];
	
	double *second_member[3];
	double *solution[3];

	

	double ****values;

	static VoronoiCell ****Domain;
	static AgentList *agentArray;	
	string border_condition;

	int SolveTridiagonaleMatrix(const int axe);
	static void PureDiffusion3D(Molecule *Mo);
	static void PureDiffusion2D(Molecule *Mo);
	static void PureDiffusion1D(Molecule *Mo);

	static void ConsumptionWhithinDiffusion1D(Molecule *Mo);
	static void ConsumptionWhithinDiffusion2D(Molecule *Mo);
	static void ConsumptionWhithinDiffusion3D(Molecule *Mo);
	static void ConsumptionDiffusionMoche1D(Molecule *Mo);
	static void ConsumptionDiffusionMoche3D(Molecule *Mo);

	static void BuildMatrixPure(int axe, int posx, int posy, int posz, Molecule *Mo, int bc);
	static void BuildMatrixPureWihtoutNecrotic(int axe, int posx, int posy, int posz, Molecule *Mo, int bc);

	static void BuildMatrixRate(int axe, int posx, int posy, int posz, Molecule *Mo, int bc);
	static void BuildMatrixRateWihtoutFree(int axe, int posx, int posy, int posz, Molecule *Mo, int bc);
	static void BuildMatrixRateWihtoutVessel(int axe, int posx, int posy, int posz, Molecule *Mo, int bc);

	static list<Molecule*> allmolecules;
	}; 
	
#endif


