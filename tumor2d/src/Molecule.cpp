#include "Molecule.h"



int Molecule::size[3];
double Molecule::dt(1);
VoronoiCell ****Molecule::Domain;
AgentList *Molecule::agentArray;	
list<Molecule*> Molecule::allmolecules;	



void Molecule::SetAgents( AgentList *theAgentArray)
	{
	agentArray = theAgentArray;
	}

void Molecule::SetDomain(int x_size,int y_size,int z_size,double thetimestep,VoronoiCell ****thedomain)
	{
	//Size of the regular domain


	#if  VERBOSE>= 1
	cout << "Molecule said: Set the domain of computation\n";
	#endif

	size[0] = x_size;
	size[1] = y_size;
	size[2] = z_size;
	
	dt = thetimestep;

	Domain = thedomain;
	}


void Molecule::SetMethod(int method, int refill, string bordercondition)
	{

	//allmolecules.push_back(this);

	#if  VERBOSE>= 1
	cout << "Molecule said: Set the method to compute in "<< DIMENSIONS <<" dimension(s)\n";
	#endif
	

	if(refill == BORDERAREBORDER)
		{
		#if  VERBOSE>= 1
		cout << "Molecule said: The borders are the borders\n";
		#endif
		if( method == PURE_DIFFUSION )
			{BuildMatrix = &Molecule::BuildMatrixPure;}
		else
			{BuildMatrix = &Molecule::BuildMatrixRate;}
		}
	
	if(refill == FREEAREBORDER)
		{
		#if  VERBOSE>= 1
		cout << "Molecule said: The borders are the initial free sites\n";
		#endif
		BuildMatrix = &Molecule::BuildMatrixRateWihtoutFree;
		}

	if(refill == VESSELAREBORDER)
		{
		#if  VERBOSE>= 1
		cout << "Molecule said: The borders are the vessels (even the news)\n";
		#endif
		BuildMatrix = &Molecule::BuildMatrixRateWihtoutVessel;
		}
	

	if(refill == NECROTICAREBORDER)
		{
		#if  VERBOSE>= 1
		cout << "Molecule said: The borders are the necrotic cells\n";
		#endif
		BuildMatrix = &Molecule::BuildMatrixPureWihtoutNecrotic;
		}

	bool methodset = false;
	//Method
	if( method == PURE_DIFFUSION )
		{
		methodset = true;

		if(DIMENSIONS == 3)
			GiveMe = &PureDiffusion3D;

		if(DIMENSIONS == 2)
			GiveMe = &PureDiffusion2D;
		
		if(DIMENSIONS == 1)
			GiveMe = &PureDiffusion1D;
		
		
		#if  VERBOSE>= 1
		if(methodset)
		cout << "Molecule said: the method is pure diffusion\n";
		#endif
		}


	if( method == CONSUMPTION_DIFFUSION )
		{

			
		if(DIMENSIONS == 3 && JAGIELLA)
			GiveMe = &ConsumptionWhithinDiffusion3D;
		
		if(DIMENSIONS == 2 && JAGIELLA)
			GiveMe = &ConsumptionWhithinDiffusion2D;
			
		if(DIMENSIONS == 1 && JAGIELLA)
			GiveMe = &ConsumptionWhithinDiffusion1D;
			
		methodset = true;
			#if  VERBOSE>= 1
			if(JAGIELLA)
			cout << "Molecule said: the method is diffusion with consumption Jagiella's way\n";
			#endif
		
		
		methodset = true;
		if(DIMENSIONS == 3 && SCHALLER)
		GiveMe = &ConsumptionDiffusionMoche3D;

		if(DIMENSIONS == 2 && SCHALLER)
		methodset = false;
		
		if(DIMENSIONS == 1 && SCHALLER)
		GiveMe = &ConsumptionDiffusionMoche1D;

			#if  VERBOSE>= 1
			if(SCHALLER)
			cout << "Molecule said: the method is diffusion with consumption Schaller's way\n";
			#endif
		}
	
	if( method == CONSUMPTION_DIFFUSION_TEST_1 )
		{

			
		if(DIMENSIONS == 3)
			GiveMe = &ConsumptionWhithinDiffusion3D;
		
		if(DIMENSIONS == 2)
			GiveMe = &ConsumptionWhithinDiffusion2D;
			
		if(DIMENSIONS == 1)
			GiveMe = &ConsumptionWhithinDiffusion1D;
			
		methodset = true;
			#if  VERBOSE>= 1
			cout << "Molecule said: the method is diffusion with consumption Jagiella's way\n";
			#endif
		}

	if( method == CONSUMPTION_DIFFUSION_TEST_2 )
		{		
		methodset = true;
		if(DIMENSIONS == 3)
		GiveMe = &ConsumptionDiffusionMoche3D;

		if(DIMENSIONS == 2)
		methodset = false;
		
		if(DIMENSIONS == 1)
		GiveMe = &ConsumptionDiffusionMoche1D;

			#if  VERBOSE>= 1
			cout << "Molecule said: the method is diffusion with consumption Schaller's way\n";
			#endif
		}
	if(!methodset)
		{cout << "Molecule said: no method to compute.\n";exit(1);}


	border_condition = bordercondition;
		

	}

void Molecule::SetMemory()
	{
	//Intialization of the Matrix (sparse)
	diagonal[0] = new double[ size[0] ];
	upperdiagonal[0] = new double[ size[0] ];
	lowerdiagonal[0] = new double[ size[0] ];
	
	diagonal[1] = new double[ size[1] ];
	upperdiagonal[1] = new double[ size[1] ];
	lowerdiagonal[1] = new double[ size[1] ];
	
	diagonal[2] = new double[ size[2] ];
	upperdiagonal[2] = new double[ size[2] ];
	lowerdiagonal[2] = new double[ size[2] ];
	
	solution[0] = new double[ size[0] ];
	second_member[0] = new double[ size[0] ];

	solution[1] = new double[ size[1] ];
	second_member[1] = new double[ size[1] ];

	solution[2] = new double[ size[2] ];
	second_member[2] = new double[ size[2] ];

	values = new double***[ size[0] ];
	for(int i = 0; i < size[0]; i++)
		{
		values[i] = new double**[ size[1] ];
		for(int j = 0; j < size[1]; j++)
			{
			values[i][j] = new double*[ size[2] ];
			}
		}
	}

void Molecule::SetValuePointer(int mx, int my, int mz, double *thevalueaddress)
	{
	values[mx][my][mz] = thevalueaddress;
	}
	
Molecule::~Molecule()
	{
	delete[] diagonal[0];
	/*	delete[] diagonal[1];
			delete[] diagonal[2];
			
	delete[] upperdiagonal[0];
		delete[] upperdiagonal[1];
			delete[] upperdiagonal[2];
			
	delete[] lowerdiagonal[0];
		delete[] lowerdiagonal[1];
			delete[] lowerdiagonal[2];
			
	delete[] solution[0];
		delete[] solution[1];
			delete[] solution[2];
			
	delete[] second_member[0];
	delete[] second_member[1];	
	delete[] second_member[2];*/
	
	}
	
void Molecule::SetMatrix(double parameter)
	{
	alpha = parameter;	
	int i;
	/*
	for(i = 0; i < size[0]; i++)
		{
		upperdiagonal[0][i] = -alpha;
		diagonal[0][i] = 1 + 2*alpha;
		lowerdiagonal[0][i] = -alpha;
		}
		
	for(i = 0; i < size[1]; i++)
		{
		upperdiagonal[1][i] = -alpha;
		diagonal[1][i] = 1 + 2*alpha;
		lowerdiagonal[1][i] = -alpha;
		}
		
	for(i = 0; i < size[2]; i++)
		{
		upperdiagonal[2][i] = -alpha;
		diagonal[2][i] = 1 + 2*alpha;
		lowerdiagonal[2][i] = -alpha;
		}
	*/
	
	for(i = 0; i < size[0]; i++)
		{
		upperdiagonal[0][i] = 0;
		diagonal[0][i] = 1;
		lowerdiagonal[0][i] = 0;
		}
		
	for(i = 0; i < size[1]; i++)
		{
		upperdiagonal[1][i] = 0;
		diagonal[1][i] = 1;
		lowerdiagonal[1][i] = 0;
		}
		
	for(i = 0; i < size[2]; i++)
		{
		upperdiagonal[2][i] = 0;
		diagonal[2][i] = 1;
		lowerdiagonal[2][i] = 0;
		}
	
	//Border conditions
	if(border_condition == "Dirichlet")
		{

		diagonal[0][0] = 1;
		diagonal[0][ size[0]  - 1] = 1;
		upperdiagonal[0][0] = 0;
		lowerdiagonal[0][ size[0]  - 1] = 0;

		diagonal[1][0] = 1;
		diagonal[1][ size[1]  - 1] = 1;
		upperdiagonal[1][0] = 0;
		lowerdiagonal[1][ size[1]  - 1] = 0;
		
		diagonal[2][0] = 1;
		diagonal[2][ size[2]  - 1] = 1;
		upperdiagonal[2][0] = 0;
		lowerdiagonal[2][ size[2]  - 1] = 0;
		
		#if VERBOSE >= 3
		cout << "Molecule said: set the matrix to Dirichlet conditions\n";
		#endif
		}

	if(border_condition == "Neumann")
		{
		#if VERBOSE >= 3
		cout << "Molecule said: set the matrix to Neumann conditions\n";
		#endif
		diagonal[0][0] = 1 + 2*alpha;
		diagonal[0][ size[0]  - 1] = 1 + 2*alpha;
		upperdiagonal[0][0] = -2*alpha;
		lowerdiagonal[0][ size[0]  - 1] = -2*alpha;
		
		diagonal[1][0] = 1 + 2*alpha;
		diagonal[1][ size[1]  - 1] = 1 + 2*alpha;
		upperdiagonal[1][0] = -2*alpha;
		lowerdiagonal[1][ size[1]  - 1] = -2*alpha;

		diagonal[2][0] = 1 + 2*alpha;
		diagonal[2][ size[2]  - 1] = 1 + 2*alpha;
		upperdiagonal[2][0] = -2*alpha;
		lowerdiagonal[2][ size[2]  - 1] = -2*alpha;
		}
	}
////////MATHS
//GENERAL MATHS METHODS
/*
int SolveTridiagonaleMatrix(const double *R, double *U, const int N, const double *alpha, const double *beta, const double *gamma)
	{
	// Solves for a vector U of length N the tridiagonal linear set
	// MU = R, where alpha, beta and gamma are the three main diagonals of matrix
	// M(N,N), the other terms are 0. R is the right side vector.
	
	//Compute the diagonale of the matrix
	// b c . . . . . . . . . . . . .
	// a b c . . . . . . . . . . . .
	// . a b c . . . . . . . . . . .
	// . . a b c . . . . . . . . . .
	// . . . a b c . . . . . . . . .
	
	
	int j;
	
	if (beta[0] == 0.0)  return 1;
	
	double G[N];
	double B = beta[0];

	U[0] = R[0]/ B;

	for (j=1; j < N; j++) 
		{         // Decomposition and forward substitution
		G[j] = gamma[j-1] / B;
		B = beta[j] - alpha[j]*G[j];
		if (B ==0.0) return 2;         // Algorithm fails
		
		U[j]=(R[j]-alpha[j]*U[j-1]) / B;
		}
	
	for (j = N - 2; j >= 0; j--)             // Back substitution
		U[j] -= G[j+1]*U[j+1];
	
	return 0;
	}*/


//second_member[axe],solution[axe],size[axe],lowerdiagonal[axe],diagonal,upperdiagonal[axe]
int Molecule::SolveTridiagonaleMatrix(const int axe)
	{
	// Solves for a vector U of length N the tridiagonal linear set
	// MU = R, where alpha, beta and gamma are the three main diagonals of matrix
	// M(N,N), the other terms are 0. R is the right side vector.
	
	//Compute the diagonale of the matrix
	// b c . . . . . . . . . . . . .
	// a b c . . . . . . . . . . . .
	// . a b c . . . . . . . . . . .
	// . . a b c . . . . . . . . . .
	// . . . a b c . . . . . . . . .

	int j;
	
	#ifdef __myDEBUG__
	if (diagonal[axe][0] == 0.0)  return 1;
	#endif
	
	double G[size[axe]];
	double B = diagonal[axe][0];

	solution[axe][0] = second_member[axe][0]/ B;

	for (j=1; j < size[axe]; j++) 
		{         // Decomposition and forward substitution
		G[j] = upperdiagonal[axe][j-1] / B;
		B = diagonal[axe][j] - lowerdiagonal[axe][j]*G[j];
		
		#ifdef __myDEBUG__
		if (B ==0.0) return 2;         // Algorithm fails
		#endif
			
		solution[axe][j]=(second_member[axe][j]-lowerdiagonal[axe][j]*solution[axe][j-1]) / B;
		}
	
	for (j = size[axe] - 2; j >= 0; j--)             // Back substitution
		solution[axe][j] -= G[j+1]*solution[axe][j+1];
	
	return 0;
	}
//////////////////



void Molecule::PureDiffusion1D(Molecule *Mo)
	{
	
	int i;
	
	int a = 0;
	if(Mo->border_condition == "Dirichlet"){a = 1;}
	Mo->BuildMatrix(0,-1,-1,-1,Mo,a);

	for(i = 0; i < size[0]; i++)
		{
		Mo->second_member[0][i] = *(Mo->values[i][0][0]);
		}



		
	for(i = a; i < size[0] - a; i++)
		{
		*(Mo->values[i][0][0]) = Mo->solution[0][i];
		if( Mo->solution[0][i] < 0){exit(1);}
		}
	
	}
void Molecule::PureDiffusion2D(Molecule *Mo)
	{

	int i,j;

	int *x = NULL;
	int *y = NULL;

	int i_size, j_size;
	int i_axe, j_axe;
	
	if( Mo->modulo % 2 == 0 )
	{x = &i; y = &j; i_size = size[0]; j_size =size[1]; i_axe = 0; j_axe = 1; }
	else
	{x = &j; y = &i; i_size = size[1]; j_size =size[0]; i_axe = 1; j_axe = 0; }
	
	int a = 0;
	//third time Step in x
	if(Mo->border_condition == "Dirichlet"){a = 1;}	
	
	
		for(j = a; j < j_size - a; j++)
			{
			for(i = 0; i < i_size; i++)
				{
				Mo->second_member[i_axe][i] =  *(Mo->values[*x][*y][0]) ;
				}
			Mo->BuildMatrix(i_axe,i,j,0,Mo,a); 

		
			for(i = a; i < i_size-a; i++)
				{
				*(Mo->values[*x][*y][0]) = Mo->solution[i_axe][i];
				if( Mo->solution[i_axe][i] < 0)
					{
				cout << "Molecule::PureDiffusion2D said Concentration negative" << endl;
					exit(1);
					}
				}
			}
		

	//third time Step in y

		for(i = a; i < i_size - a; i++)
			{
			for(j = 0; j < j_size; j++)
				{
				Mo->second_member[j_axe][j] = *(Mo->values[*x][*y][0]);
				}
			Mo->BuildMatrix(j_axe,i,j,0,Mo,a); 
			
			if(Mo->SolveTridiagonaleMatrix(j_axe) != 0)
				cout << "Matrix Error" << endl;
			
			for(j = a; j < j_size-a; j++)
				{*(Mo->values[*x][*y][0]) = Mo->solution[j_axe][j];}
			}
	}


void Molecule::PureDiffusion3D(Molecule *Mo)	
	{

	int i,j,k;

	int *x = NULL;
	int *y = NULL;
	int *z = NULL;
	int i_size, j_size, k_size;
	int i_axe, j_axe, k_axe;
	
	switch( Mo->modulo )
		{
		case 0:
		x = &i; y = &j; z = &k; i_size = size[0]; j_size =size[1]; k_size = size[2];i_axe = 0; j_axe = 1; k_axe = 2;
		break;
		
		case 1:
		x = &j; y = &k; z = &i; i_size = size[2]; j_size =size[0]; k_size = size[1];i_axe = 2; j_axe = 0; k_axe = 1;
		break;

		case 2:
		x = &k; y = &i; z = &j; i_size = size[1]; j_size =size[2]; k_size = size[0];i_axe = 1; j_axe = 2; k_axe = 0;
		break;

		case 3:
		x = &i; y = &k; z = &j; i_size = size[0]; j_size =size[2]; k_size = size[1];i_axe = 0; j_axe = 2; k_axe = 1;
		break;

		case 4:
		x = &j; y = &i; z = &k; i_size = size[1]; j_size =size[0]; k_size = size[2];i_axe = 1; j_axe = 0; k_axe = 2;
		break;

		case 5:
		x = &k; y = &j; z = &i; i_size = size[2]; j_size =size[1]; k_size = size[0];i_axe = 2; j_axe = 1; k_axe = 0;
		break;

		default:
		x = &i; y = &j; z = &k; i_size = size[0]; j_size =size[1]; k_size = size[2];i_axe = 0; j_axe = 1; k_axe = 2;
		cerr << endl << "PureDiffusion3D says: what's that modulo?" << endl;
		exit(1);
		break;
		}

	int a = 0;
	//third time Step in x
	if(Mo->border_condition == "Dirichlet"){a = 1;}	

	

	
	for(k = a; k < k_size - a; k++)
		{
		for(j = a; j < j_size - a; j++)
			{
			for(i = 0; i < i_size; i++)
				{
				Mo->second_member[i_axe][i] = *(Mo->values[*x][*y][*z]);
				}
			Mo->BuildMatrix(i_axe,i,j,k,Mo,a); 
			
			
			if(Mo->SolveTridiagonaleMatrix(i_axe) != 0)
				cout << "Matrix Error" << endl;
		
			for(i = a; i < i_size-a; i++)
				{
				*(Mo->values[*x][*y][*z]) = Mo->solution[i_axe][i];
				if( Mo->solution[i_axe][i] < 0)
					{
				cout << "Molecule::PureDiffusion3D said Concentration negative" << endl;
					exit(1);
					}
				}
			}
		}

	//third time Step in y
	for(k = a; k < k_size - a; k++)
		{
		for(i = a; i < i_size - a; i++)
			{
			for(j = 0; j < j_size; j++)
				{
				Mo->second_member[j_axe][j] = *(Mo->values[*x][*y][*z]);
				}
			Mo->BuildMatrix(j_axe,i,j,k,Mo,a); 
			
			if(Mo->SolveTridiagonaleMatrix(j_axe) != 0)
				cout << "Matrix Error" << endl;
			
			for(j = a; j < j_size-a; j++)
				{*(Mo->values[*x][*y][*z]) = Mo->solution[j_axe][j];}
			}
		}

	//third time Step in z
	for(j = a; j < j_size - a; j++)
		{
		for(i = a; i < i_size - a; i++)
			{
			
			
			for(k = 0; k < k_size; k++)
				{
				Mo->second_member[k_axe][k] = *(Mo->values[*x][*y][*z]);	
				}

			Mo->BuildMatrix(k_axe,i,j,k,Mo,a);
		
			if(Mo->SolveTridiagonaleMatrix(k_axe) != 0)
				cout << "Matrix Error" << endl;
		
			for(k = a; k < k_size-a; k++)
				{*(Mo->values[*x][*y][*z]) = Mo->solution[k_axe][k];}
			}
		}

	}

//Schaller - DIFFUSION
void Molecule::ConsumptionDiffusionMoche1D(Molecule *Mo)
	{
	int i;
	
	int a = 0;
	if(Mo->border_condition == "Dirichlet"){a = 1;}
	Mo->BuildMatrix(0,-1,0,0,Mo,a);

	for(i = 0; i < size[0]; i++)
		{
		Mo->second_member[0][i] = (*(Mo->values[i][0][0])) - dt*Mo->Rate(Domain[i][0][0]);
		
		
		if( Mo->second_member[0][i] < 0)
			{
			#if  VERBOSE>= 7
			cout << "Concentration negative " << Mo->second_member[i_axe][i] << endl;
			#endif
			Mo->second_member[0][i] = 0.;
			}
		
		}


	if(Mo->SolveTridiagonaleMatrix(0) != 0)
		cout << "Matrix Error" << endl;
		
	for(i = a; i < size[0] - a; i++)
		{
		*(Mo->values[i][0][0]) = Mo->solution[0][i];
		if( Mo->solution[0][i] < 0){cout << "Concentration negative" << endl; exit(1);}
		}
	


	}

void Molecule::ConsumptionDiffusionMoche3D(Molecule *Mo)
	{
	
	int i,j,k;

	int *x = NULL;
	int *y = NULL;
	int *z = NULL;
	int i_size, j_size, k_size;
	int i_axe, j_axe, k_axe;
	
	switch( Mo->modulo )
		{
		case 0:
		x = &i; y = &j; z = &k; i_size = size[0]; j_size =size[1]; k_size = size[2];i_axe = 0; j_axe = 1; k_axe = 2;
		break;
		
		case 1:
		x = &j; y = &k; z = &i; i_size = size[2]; j_size =size[0]; k_size = size[1];i_axe = 2; j_axe = 0; k_axe = 1;
		break;

		case 2:
		x = &k; y = &i; z = &j; i_size = size[1]; j_size =size[2]; k_size = size[0];i_axe = 1; j_axe = 2; k_axe = 0;
		break;

		case 3:
		x = &i; y = &k; z = &j; i_size = size[0]; j_size =size[2]; k_size = size[1];i_axe = 0; j_axe = 2; k_axe = 1;
		break;

		case 4:
		x = &j; y = &i; z = &k; i_size = size[1]; j_size =size[0]; k_size = size[2];i_axe = 1; j_axe = 0; k_axe = 2;
		break;

		case 5:
		x = &k; y = &j; z = &i; i_size = size[2]; j_size =size[1]; k_size = size[0];i_axe = 2; j_axe = 1; k_axe = 0;
		break;

		default:
		x = &i; y = &j; z = &k; i_size = size[0]; j_size =size[1]; k_size = size[2];i_axe = 0; j_axe = 1; k_axe = 2;
		cerr << endl << "ConsumptionDiffusionMoche3D says: what's that modulo?" << endl;
		exit(1);
		break;
		}

	int a = 0;
	//third time Step in x
	if(Mo->border_condition == "Dirichlet")
		{
		//To avoid the computation on the border
		a = 1;
		}	
	
	for(k = a; k < k_size - a; k++)
		{
		for(j = a; j < j_size - a; j++)
			{
			for(i = 0; i < i_size; i++)
				{
				Mo->second_member[i_axe][i] = (*(Mo->values[*x][*y][*z])) - dt*Mo->Rate(Domain[*x][*y][*z]);
				if( Mo->second_member[i_axe][i] < 0)
					{
					#if  VERBOSE>= 7
					cout << "Concentration negative " << Mo->second_member[i_axe][i] << endl;
					#endif
					Mo->second_member[i_axe][i] = 0.;
					}
				}

			Mo->BuildMatrix(i_axe,i,j,k,Mo,a); 

			if(Mo->SolveTridiagonaleMatrix(i_axe) != 0)
				cout << "Matrix Error" << endl;
		
			for(i = a; i < i_size-a; i++)
				{
				*(Mo->values[*x][*y][*z]) = Mo->solution[i_axe][i];
				if( Mo->solution[i_axe][i] < 0){cout << "Concentration negative" << endl; exit(1);}
				}
			}
		}

	//third time Step in y
	for(k = a; k < k_size - a; k++)
		{
		for(i = a; i < i_size - a; i++)
			{
			for(j = 0; j < j_size; j++)
			{Mo->second_member[j_axe][j] = *(Mo->values[*x][*y][*z]);}

			Mo->BuildMatrix(j_axe,i,j,k,Mo,a); 

			if(Mo->SolveTridiagonaleMatrix(j_axe) != 0)
				cout << "Matrix Error" << endl;
			
			for(j = a; j < j_size-a; j++)
				{*(Mo->values[*x][*y][*z]) = Mo->solution[j_axe][j];}
			}
		}

	//third time Step in z
	for(j = a; j < j_size - a; j++)
		{
		for(i = a; i < i_size - a; i++)
			{
			for(k = 0; k < k_size; k++)
			{Mo->second_member[k_axe][k] = *(Mo->values[*x][*y][*z]);}

			Mo->BuildMatrix(k_axe,i,j,k,Mo,a); 

			if(Mo->SolveTridiagonaleMatrix(k_axe) != 0)
				cout << "Matrix Error" << endl;
		
			for(k = a; k < k_size-a; k++)
				{*(Mo->values[*x][*y][*z]) = Mo->solution[k_axe][k];}
			}
		}

	}



//JAGIELLA
void Molecule::ConsumptionWhithinDiffusion1D(Molecule *Mo)
	{
	int i;
	double rate;
	int a = 0;
	if(Mo->border_condition == "Dirichlet"){a = 1;}
	
	Mo->BuildMatrix(0,-1,0,0,Mo,a); 

	for(i = a; i < size[0] - a; i++)
		{
		rate = dt*Mo->Rate(Domain[i][0][0]);
		Mo->second_member[0][i] = ( *(Mo->values[i][0][0]) )/ (1 + rate);
		}

	
	if(Mo->SolveTridiagonaleMatrix(0) != 0)
		cout << "Matrix Error" << endl;
		
	for(i = a; i < size[0] - a; i++)
		{
		*(Mo->values[i][0][0]) = Mo->solution[0][i];
		if( Mo->solution[0][i] < 0){cout << "Concentration negative" << endl; exit(1);}
		}
	}

void Molecule::ConsumptionWhithinDiffusion2D(Molecule *Mo)
	{

	int i,j;

	int *x = NULL;
	int *y = NULL;

	int i_size, j_size;
	int i_axe, j_axe;
	
	if( Mo->modulo % 2 == 0 )
	{x = &i; y = &j; i_size = size[0]; j_size =size[1]; i_axe = 0; j_axe = 1; }
	else
	{x = &j; y = &i; i_size = size[1]; j_size =size[0]; i_axe = 1; j_axe = 0; }
	
	double rate;
	int a = 0;
	//third time Step in x
	if(Mo->border_condition == "Dirichlet"){a = 1;}	
	
	
		for(j = a; j < j_size - a; j++)
			{
			for(i = 0; i < i_size; i++)
				{
				rate = dt*Mo->Rate(Domain[*x][*y][0]);
				Mo->second_member[i_axe][i] = ( *(Mo->values[*x][*y][0]) )/ (1 + rate);
				}
			Mo->BuildMatrix(i_axe,i,j,0,Mo,a); 
			
			
			if(Mo->SolveTridiagonaleMatrix(i_axe) != 0)
				cout << "Matrix Error" << endl;
		
			for(i = a; i < i_size-a; i++)
				{
				*(Mo->values[*x][*y][0]) = Mo->solution[i_axe][i];
				if( Mo->solution[i_axe][i] < 0)
					{
				cout << "Molecule::ConsumptionWhithinDiffusion2D said Concentration negative" << endl;
					exit(1);
					}
				}
			}
		

	//third time Step in y

		for(i = a; i < i_size - a; i++)
			{
			for(j = 0; j < j_size; j++)
				{
				Mo->second_member[j_axe][j] = *(Mo->values[*x][*y][0]);
				}
			Mo->BuildMatrix(j_axe,i,j,0,Mo,a); 
			
			if(Mo->SolveTridiagonaleMatrix(j_axe) != 0)
				cout << "Matrix Error" << endl;
			
			for(j = a; j < j_size-a; j++)
				{*(Mo->values[*x][*y][0]) = Mo->solution[j_axe][j];}
			}
	}
	
void Molecule::ConsumptionWhithinDiffusion3D(Molecule *Mo)
	{
		

	int i,j,k;

	int *x = NULL;
	int *y = NULL;
	int *z = NULL;
	int i_size, j_size, k_size;
	int i_axe, j_axe, k_axe;
	
	switch( Mo->modulo )
		{
		case 0:
		x = &i; y = &j; z = &k; i_size = size[0]; j_size =size[1]; k_size = size[2];i_axe = 0; j_axe = 1; k_axe = 2;
		break;
		
		case 1:
		x = &j; y = &k; z = &i; i_size = size[2]; j_size =size[0]; k_size = size[1];i_axe = 2; j_axe = 0; k_axe = 1;
		break;

		case 2:
		x = &k; y = &i; z = &j; i_size = size[1]; j_size =size[2]; k_size = size[0];i_axe = 1; j_axe = 2; k_axe = 0;
		break;

		case 3:
		x = &i; y = &k; z = &j; i_size = size[0]; j_size =size[2]; k_size = size[1];i_axe = 0; j_axe = 2; k_axe = 1;
		break;

		case 4:
		x = &j; y = &i; z = &k; i_size = size[1]; j_size =size[0]; k_size = size[2];i_axe = 1; j_axe = 0; k_axe = 2;
		break;

		case 5:
		x = &k; y = &j; z = &i; i_size = size[2]; j_size =size[1]; k_size = size[0];i_axe = 2; j_axe = 1; k_axe = 0;
		break;

		default:
		x = &i; y = &j; z = &k; i_size = size[0]; j_size =size[1]; k_size = size[2];i_axe = 0; j_axe = 1; k_axe = 2;
		cerr << endl << "ConsumptionWhithinDiffusion3D says: what's that modulo?" << endl;
		exit(1);
		break;
		}




	for( i=0; i<size[0]; i++)
	for( j=0; j<size[1]; j++)
	for( k=0; k<size[2]; k++)
	if( Domain[i][j][k]->agent != NULL && Domain[i][j][k]->getState() == FREE){
		


		int l;
		for( l=0; l<Domain[i][j][k]->countNeighborCells && Domain[i][j][k]->neighborCells[l]->agent != NULL; l++);
		if( l != Domain[i][j][k]->countNeighborCells){
			if( GetAgent( Domain[i][j][k])->countLocations == 1){
				agentArray->deactivateAgent( GetAgent( Domain[i][j][k]));
			}
			GetAgent( Domain[i][j][k])->detach( Domain[i][j][k]);
			Domain[i][j][k]->oxygen  = Domain[i][j][k]->neighborCells[l]->oxygen;
			Domain[i][j][k]->glucose = Domain[i][j][k]->neighborCells[l]->glucose;
		}
	}
	
	double rate;
	int a = 0;
	//third time Step in x
	if(Mo->border_condition == "Dirichlet"){a = 1;}	

	

	
	for(k = a; k < k_size - a; k++)
		{
		for(j = a; j < j_size - a; j++)
			{
			for(i = 0; i < i_size; i++)
				{
				rate = dt*Mo->Rate(Domain[*x][*y][*z]);
				Mo->second_member[i_axe][i] = ( *(Mo->values[*x][*y][*z]) )/ (1 + rate);
				}
			Mo->BuildMatrix(i_axe,i,j,k,Mo,a); 
			
			
			if(Mo->SolveTridiagonaleMatrix(i_axe) != 0)
				cout << "Matrix Error" << endl;
		
			for(i = a; i < i_size-a; i++)
				{
				*(Mo->values[*x][*y][*z]) = Mo->solution[i_axe][i];
				if( Mo->solution[i_axe][i] < 0)
					{
				cout << "Molecule::ConsumptionWhithinDiffusion3D said Concentration negative" << endl;
					exit(1);
					}
				}
			}
		}

	//third time Step in y
	for(k = a; k < k_size - a; k++)
		{
		for(i = a; i < i_size - a; i++)
			{
			for(j = 0; j < j_size; j++)
				{
				Mo->second_member[j_axe][j] = *(Mo->values[*x][*y][*z]);
				}
			Mo->BuildMatrix(j_axe,i,j,k,Mo,a); 
			
			if(Mo->SolveTridiagonaleMatrix(j_axe) != 0)
				cout << "Matrix Error" << endl;
			
			for(j = a; j < j_size-a; j++)
				{*(Mo->values[*x][*y][*z]) = Mo->solution[j_axe][j];}
			}
		}

	//third time Step in z
	for(j = a; j < j_size - a; j++)
		{
		for(i = a; i < i_size - a; i++)
			{
			
			
			for(k = 0; k < k_size; k++)
				{
				Mo->second_member[k_axe][k] = *(Mo->values[*x][*y][*z]);	
				}

			Mo->BuildMatrix(k_axe,i,j,k,Mo,a);
		
			if(Mo->SolveTridiagonaleMatrix(k_axe) != 0)
				cout << "Matrix Error" << endl;
		
			for(k = a; k < k_size-a; k++)
				{*(Mo->values[*x][*y][*z]) = Mo->solution[k_axe][k];}
			}
		}

	}
	

	
void Molecule::BuildMatrixPure(int axe, int posx, int posy, int posz, Molecule *Mo, int bc)
	{
	int n;	
	for(n = bc; n < size[axe] - bc; n++)
		{
		Mo->diagonal[axe][n] = 1 + 2*Mo->alpha;
		Mo->lowerdiagonal[axe][n] = -Mo->alpha;
		Mo->upperdiagonal[axe][n] = -Mo->alpha;
		}
	if(bc == 0)
		{
		n = 0;
		Mo->upperdiagonal[axe][n] = -2*Mo->alpha;

		n = size[axe] - 1;	
		Mo->lowerdiagonal[axe][n] = -2*Mo->alpha;
		}
	}

void Molecule::BuildMatrixRate(int axe, int posx, int posy, int posz, Molecule *Mo, int bc)
	{
	
	int *i = &posx, *j = &posy, *k = &posz;
	int n;
	switch( axe )
		{
		case 0:
		i = &n;
		break;
		
		case 1:
		j = &n;
		break;
		
		case 2:
		k = &n;
		break;
		
		default:
		i = &n;
		break;
		}
		
	double rate;
	
	for(n = bc; n < size[axe] - bc; n++)
		{
		rate = dt*Mo->Rate(Domain[*i][*j][*k]);
		Mo->diagonal[axe][n] = 1 + 2*Mo->alpha  /  (1 + rate);
		Mo->lowerdiagonal[axe][n] = -Mo->alpha  /  (1 + rate);
		Mo->upperdiagonal[axe][n] = -Mo->alpha  /  (1 + rate);
		}
	if(bc == 0)
		{
		n = 0;
		rate = dt*Mo->Rate(Domain[*i][*j][*k]);
		Mo->upperdiagonal[axe][n] = -2*Mo->alpha /  (1 + rate);

		n = size[axe] - 1;
		rate = dt*Mo->Rate(Domain[*i][*j][*k]);
		Mo->lowerdiagonal[axe][n] = -2*Mo->alpha /  (1 + rate);
		}
	}

void Molecule::BuildMatrixRateWihtoutFree(int axe, int posx, int posy, int posz, Molecule *Mo, int bc)
	{
	//This function build the matrix as if free sites were border conditions.
	int *i = &posx, *j = &posy, *k = &posz;
	int n;
	switch( axe )
		{
		case 0:
		i = &n;
		break;
		
		case 1:
		j = &n;
		break;
		
		case 2:
		k = &n;
		break;
		
		default:
		i = &n;
		break;
		}
	double rate;
	
	for(n = bc; n < size[axe] - bc; n++)
		{

		if(Domain[*i][*j][*k]->agent != NULL )
			{
			rate = dt*Mo->Rate(Domain[*i][*j][*k]);
			Mo->diagonal[axe][n] = 1 + 2*Mo->alpha  /  (1 + rate);
			Mo->lowerdiagonal[axe][n] = -Mo->alpha  /  (1 + rate);
			Mo->upperdiagonal[axe][n] = -Mo->alpha  /  (1 + rate);

			}
		else
			{
			Mo->diagonal[axe][n] = 1 ;
			Mo->lowerdiagonal[axe][n] = 0;
			Mo->upperdiagonal[axe][n] = 0;
			}
		}
	
	if(bc == 0)
		{

		if(Domain[*i][*j][*k]->agent != NULL )
			{
			n = 0;
			rate = dt*Mo->Rate(Domain[*i][*j][*k]);
			Mo->upperdiagonal[axe][n] = -2*Mo->alpha /  (1 + rate);
			n = size[axe] - 1;
			rate = dt*Mo->Rate(Domain[*i][*j][*k]);
			Mo->lowerdiagonal[axe][n] = -2*Mo->alpha /  (1 + rate);
			}
		}
	}

void Molecule::BuildMatrixPureWihtoutNecrotic(int axe, int posx, int posy, int posz, Molecule *Mo, int bc)
	{
	//This function build the matrix as if free sites were border conditions.
	int *i = &posx, *j = &posy, *k = &posz;
	int n;
	switch( axe )
		{
		case 0:
		i = &n;
		break;
		
		case 1:
		j = &n;
		break;
		
		case 2:
		k = &n;
		break;
		
		default:
		i = &n;
		break;
		}

	
	for(n = bc; n < size[axe] - bc; n++)
		{
		if(Domain[*i][*j][*k]->getState() != NECROTIC )
			{
			Mo->diagonal[axe][n] = 1 + 2*Mo->alpha;
			Mo->lowerdiagonal[axe][n] = -Mo->alpha;
			Mo->upperdiagonal[axe][n] = -Mo->alpha;
			//cout << "i=" << *i << "\tj=" << *j << "\tk="<< *k <<endl;
			}
		else
			{
			Mo->diagonal[axe][n] = 1 ;
			Mo->lowerdiagonal[axe][n] = 0;
			Mo->upperdiagonal[axe][n] = 0;
			}
		}
	
	if(bc == 0)
		{
		if(Domain[*i][*j][*k]->getState() != NECROTIC )
			{
			n = 0;
			Mo->upperdiagonal[axe][n] = -2*Mo->alpha;
			n = size[axe] - 1;
			Mo->lowerdiagonal[axe][n] = -2*Mo->alpha;
			}
		}
	}

void Molecule::BuildMatrixRateWihtoutVessel(int axe, int posx, int posy, int posz, Molecule *Mo, int bc)
	{
	//This function build the matrix as if free sites were border conditions.
	int *i = &posx, *j = &posy, *k = &posz;
	int n;
	switch( axe )
		{
		case 0:
		i = &n;
		break;
		
		case 1:
		j = &n;
		break;
		
		case 2:
		k = &n;
		break;
		
		default:
		i = &n;
		break;
		}
	double rate;
	
	for(n = bc; n < size[axe] - bc; n++)
		{
		if(Domain[*i][*j][*k]->getState() != VESSEL )
			{
			rate = dt*Mo->Rate(Domain[*i][*j][*k]);
			Mo->diagonal[axe][n] = 1 + 2*Mo->alpha  /  (1 + rate);
			Mo->lowerdiagonal[axe][n] = -Mo->alpha  /  (1 + rate);
			Mo->upperdiagonal[axe][n] = -Mo->alpha  /  (1 + rate);
			}
		else
			{
			Mo->diagonal[axe][n] = 1 ;
			Mo->lowerdiagonal[axe][n] = 0;
			Mo->upperdiagonal[axe][n] = 0;
			}
		}

	if(bc == 0)
		{
		if(Domain[*i][*j][*k]->getState() != VESSEL )
			{
			n = 0;
			rate = dt*Mo->Rate(Domain[*i][*j][*k]);
			Mo->upperdiagonal[axe][n] = -2*Mo->alpha /  (1 + rate);
			n = size[axe] - 1;
			rate = dt*Mo->Rate(Domain[*i][*j][*k]);
			Mo->lowerdiagonal[axe][n] = -2*Mo->alpha /  (1 + rate);
			}
		}
	}



