#include "Angiogenesis.h"


#ifdef __myDEBUG__
int KillNode_counter = 0;
int KillSlice_counter = 0;
int MakeOneEdge_counter = 0;
int DoTheDilation_counter = 0;
int KillConnection_counter = 0;
int Ko_counter = 0;
int Regression_counter = 0;
int DoTheSprout_counter = 0;
#endif

Angiogenesis::Angiogenesis()
	{
	nbrNodes = 0;
	maxNodes = 0;
	Node = NULL;

	nbrSlices = 0;
	maxSlices = 0;
	Slice = NULL;
	}

Angiogenesis::~Angiogenesis()
	{
	int i;
	for(i = 0; i < nbrSlices; i++)
		{
		free(Slice[i]);
		}
	free(Slice);
	
	for(i = 0; i < nbrNodes; i++)
		{
		free(Node[i]->branch);
		free(Node[i]->edge);
		free(Node[i]);	
		}
	free(Node);
	}

void Angiogenesis::Initialize(AgentList *agentList)
	{

	angioAgentList = agentList;
	VasculatureCreation();
	GiveMeTheBlood();


	}



void Angiogenesis::KillSlice(Vessel_Unit *X)	
	{
	#ifdef __myDEBUG__
	KillSlice_counter++;
	#endif

	int I = X->index;
	Vessel_Unit *Temp;
	nbrSlices--;
	Temp = Slice[nbrSlices];
	
	Slice[nbrSlices] = X;
	Slice[nbrSlices]->index = nbrSlices;
	
	Slice[I] = Temp;
	Slice[I]->index = I;
	}

void Angiogenesis::MakeOneEdge(Vessel_Graph *a, Vessel_Graph *b)
	{
	#ifdef __myDEBUG__
	MakeOneEdge_counter++;
	#endif
	int N = nbrSlices;
	nbrSlices++;
	if(N == maxSlices)
		{
		maxSlices++;
		Slice = (Vessel_Unit **) realloc(Slice, sizeof(Vessel_Unit*) * (maxSlices));
		Slice[N] = new Vessel_Unit;	
		}
	Slice[N]->index = N;
	Slice[N]->circulation = true;
	Slice[N]->node1 = a;
	Slice[N]->node2 = b;
	Slice[N]->shear = 10*MINSHEAR; // to avoid the collapse just after birth and before the computation of the real value	
	Slice[N]->flow = 1;
	Slice[N]->velocity = 1;
	Slice[N]->Compute_length();
	Slice[N]->Compute_theta();
	Slice[N]->Compute_phi();
	
	Slice[N]->Compute_radius();
	Slice[N]->Compute_viscosity();
	
	
	int m;
	m = a->FindBranchNumber(b);
	a->edge[m] = Slice[N];
	m = b->FindBranchNumber(a);
	b->edge[m] = Slice[N];

	
	}



void Angiogenesis::KillNode(Vessel_Graph *N)
	{
	#ifdef __myDEBUG__
	KillNode_counter++;
	#endif

    	int I = N->index;
	//delete the connections with the other nodes
	int j,k;
	Vessel_Graph *temp;
	Vessel_Unit *tamp;
	for(j = 0; j < N->CountBranches; j++)
		{
		//Kill every slices between the nodes
		KillSlice(N->edge[j]);

		k = N->branch[j]->FindBranchNumber(N);
		N->branch[j]->CountBranches -= 1;

		temp = N->branch[j]->branch[N->branch[j]->CountBranches];
		N->branch[j]->branch[N->branch[j]->CountBranches] = N->branch[j]->branch[k];
		N->branch[j]->branch[k] = temp;

		tamp = N->branch[j]->edge[N->branch[j]->CountBranches];
		N->branch[j]->edge[N->branch[j]->CountBranches] = N->branch[j]->edge[k];
		N->branch[j]->edge[k] = tamp;
		}
	nbrNodes--;
	temp = Node[nbrNodes];
	Node[nbrNodes] = N;
	Node[nbrNodes]->index = nbrNodes;
	
	Node[I] = temp;
	Node[I]->index = I;
	//printf("Oh my god, They killed a node\n");
    }

Vessel_Graph* Angiogenesis::MakeOneNode()
	{
	int P = nbrNodes;
	nbrNodes++;
	if(P == maxNodes)
		{
		maxNodes++;
		Node = (Vessel_Graph **) realloc(Node, sizeof(Vessel_Graph*) * (maxNodes));
		Node[P] = new Vessel_Graph;
		Node[P]->circulation = false;
		}
		
	Node[P]->index = P;
	Node[P]->CountBranches = 0;
	Node[P]->branch = NULL;
	Node[P]->edge = NULL;
	Node[P]->pressure = 0;
	Node[P]->radius = INIT_RADIUS;
	Node[P]->viscosity = INIT_VISCOSITY;
	return Node[P];
	}


Vessel_Graph* Angiogenesis::FindMyNode(Agent *A)
	{
	for(int i = 0; i < nbrNodes; i++)
		{
		if(Node[i]->p_agent == A){return Node[i];}
		}
	return NULL;
	}

bool Angiogenesis::DoesExistTheSlice(Vessel_Graph* A, Vessel_Graph* B)
	{
	int k;
	
	for(k = 0; k < nbrSlices; k++)
		{
		if(Slice[k]->node1 == A && Slice[k]->node2 == B)
			{return true;}
		if(Slice[k]->node2 == A && Slice[k]->node1 == B)
			{return true;}
		}

	return false;
	}


void Angiogenesis::VasculatureCreation() //Initialization of the function
	{
	//Initialization of the network
	nbrNodes = 0;
	maxNodes = 0;
	Node = NULL;

	nbrSlices = 0;
	maxSlices = 0;
	Slice = NULL;
	
	int i,j;
	
	Vessel_Graph *N;
	for(i = 0; i < angioAgentList->countActiveAgents; i++)
		{

		if( angioAgentList->agents[i]->state == VESSEL)
			{

			N = MakeOneNode();
			N->p_agent = angioAgentList->agents[i];
			//N->p_voronoiCell = GetVoronoiCell( angioAgentList->agents[i]);
			N->x = GetVoronoiCell( N->p_agent )->position[0];
			N->y = GetVoronoiCell( N->p_agent )->position[1];
			N->z = GetVoronoiCell( N->p_agent )->position[2];
			}
		}
	
	for(i = 0; i < nbrNodes; i++)
		{
		for(j = 0; j < GetVoronoiCell( Node[i]->p_agent )->countNeighborCells; j++)
			{
			if( GetVoronoiCell( Node[i]->p_agent )->neighborCells[j]->getState()== VESSEL)
				{
	Node[i]->CountBranches++;
	Node[i]->branch = (Vessel_Graph**) realloc(Node[i]->branch, sizeof(Vessel_Graph*) * Node[i]->CountBranches);
	Node[i]->edge = (Vessel_Unit**) realloc(Node[i]->edge, sizeof(Vessel_Unit*) * Node[i]->CountBranches);
	Node[i]->branch[Node[i]->CountBranches - 1] = FindMyNode(GetAgent( GetVoronoiCell( Node[i]->p_agent )->neighborCells[j] ));
				}


			}
		}
	for(i = 0; i < nbrNodes; i++)
		{
		for(j = 0; j < Node[i]->CountBranches; j++)
			{
			if(!DoesExistTheSlice(Node[i],Node[i]->branch[j]))
				{MakeOneEdge(Node[i],Node[i]->branch[j]);}
			}
		}
	#if VERBOSE >= 1
	printf("Numbre of Nodes created: %i\n",nbrNodes);
	printf("Numbre of Slices created: %i\n",nbrSlices);
	#endif
	for(i = 0; i < nbrNodes; i++)
		{
		Node[i]->InDomain = Node[i]->DefineInDomain();
		}


	#ifdef __myDEBUG__
	cerr << "Check if every vessel is been added as a node in the list\n";
	for(i = 0; i < VoronoiCell_Number; i++)
		{
		if( VoronoiCell_List[i]->getState() == VESSEL)
			{
			FindMyNode(VoronoiCell_List[i]);	
			}
		}
	#endif
	}


bool DoWeKnowEachOther(Vessel_Graph *N, Vessel_Graph *M)
     {
     //Check if the nodes are already connected
     if(N == M){return true;}
     int j;
     for(j = 0; j < N->CountBranches; j++)
           {
           if(N->branch[j] == M){return true;}
           }
     return false;
     }

void Angiogenesis::DoTheFusion(Vessel_Graph *N, Vessel_Graph *M)
	{
	M->CountBranches++;
	M->branch = (Vessel_Graph**) realloc(M->branch, sizeof(Vessel_Graph*)*M->CountBranches);
	M->edge = (Vessel_Unit**) realloc(M->edge, sizeof(Vessel_Unit*)*M->CountBranches);
	M->branch[M->CountBranches - 1]  = N;
	
	N->CountBranches++;
	N->branch = (Vessel_Graph**) realloc(N->branch, sizeof(Vessel_Graph*) * N->CountBranches);
	N->edge = (Vessel_Unit**) realloc(N->edge, sizeof(Vessel_Unit*) * N->CountBranches);
	N->branch[N->CountBranches - 1] = M;
		
	MakeOneEdge(N,M);
	}

bool Angiogenesis::DoTheDilation(Vessel_Graph *N)
	{

	#ifdef __myDEBUG__
	DoTheDilation_counter++;
	#endif
	int i;



     N->radius += INIT_RADIUS/(2.0*PI);
     if( N->radius > MAXRADIUS){N->radius = MAXRADIUS;}
     for(i = 0; i < N->CountBranches ;i++)
         {
         N->edge[i]->Compute_radius();
       	 N->edge[i]->Compute_viscosity();
         }
     return true;
     }


bool Angiogenesis::DoTheSprout(Vessel_Graph *N, ActionTree *p_list, VoronoiDiagram* voronoiDiagram)
	{
	

	//RULES
	//Hypothèses de Rieger pour le sprouting
	//ok Si la place est FREE
	//Si aucun vessel aux alentours
	//ok Si inferieure à N cellules de distance
	//ok Dans la direction des VEGF

	//Hypothèses supplémentaires
		//Premier vessel rencontré forme la fin du vaisseau.

	int i,j,count;

	//Keep only the first with no vessels around except the starting point. (prevent loops and bridges size one.)
	VoronoiCell *first = NULL;
	double max_vegf = 0.;
	for(i = 0; i < GetVoronoiCell( N->p_agent )->countNeighborCells; i++)
		{
		if( GetVoronoiCell( N->p_agent )->neighborCells[i]->getState() == FREE)
			{
			count = 0;
			for(j = 0; j < GetVoronoiCell( N->p_agent )->neighborCells[i]->countNeighborCells; j++)
				{
				if(GetVoronoiCell( N->p_agent )->neighborCells[i]->neighborCells[j]->getState() == VESSEL)
					{count++;}
				}
			if(count < 2 && GetVoronoiCell( N->p_agent )->neighborCells[i]->growthfactors > max_vegf){first = GetVoronoiCell( N->p_agent )->neighborCells[i]; max_vegf = first->growthfactors;}
			}
		}

	if(first == NULL)
		return false;


	VoronoiCell*tab[MAXSPROUT+1];
	tab[0] = GetVoronoiCell( N->p_agent );
	tab[1] = first;
	bool ivefoundone;
	bool vesselend = false;
	for(i = 1; i < MAXSPROUT && !vesselend; i++)
		{
		ivefoundone = false;
		max_vegf = 0.;
		for(j = 0; j < tab[i]->countNeighborCells; j++)
			{
			if(tab[i]->neighborCells[j]->growthfactors > max_vegf && tab[i]->neighborCells[j] != tab[i-1])
				{
				//Or this is a vessel and not the father, this means, the end of the sprouting
				if(tab[i]->neighborCells[j]->getState() == VESSEL && tab[i]->neighborCells[j] != tab[0])
					{
					tab[i+1] = tab[i]->neighborCells[j];
					vesselend = true;
					ivefoundone = true;
					max_vegf = tab[i]->neighborCells[j]->growthfactors;
					}
				//or this is free
				if(!vesselend && tab[i]->neighborCells[j]->getState() == FREE)
					{
					tab[i+1] = tab[i]->neighborCells[j];
					ivefoundone = true;
					max_vegf = tab[i]->neighborCells[j]->growthfactors;
					}
				}
			}
		if(ivefoundone == false)
			{return false;}
		}
		
	if(vesselend == false)
		{return false;}


	i = 1;
	Vessel_Graph *VGtab[MAXSPROUT+1];
	int m = 0;
	while(tab[i]->getState() != VESSEL)
		{
		addedVesselAndUpdateSurrounding( p_list, angioAgentList, tab[i], voronoiDiagram);
		
		
		VGtab[m] = MakeOneNode();
		VGtab[m]->p_agent = GetAgent( tab[i] );
		VGtab[m]->x = GetVoronoiCell( VGtab[m]->p_agent )->position[0];
		VGtab[m]->y = GetVoronoiCell( VGtab[m]->p_agent )->position[1];
		VGtab[m]->z = GetVoronoiCell( VGtab[m]->p_agent )->position[2];
		VGtab[m]->radius = INIT_RADIUS;
		VGtab[m]->InDomain = VGtab[m]->DefineInDomain();


		m++;
		i++;
		}
	
	Vessel_Graph *last = FindMyNode(GetAgent( tab[i] ) );


	#if defined ( __myDEBUG__ ) && VERBOSE >= 1
	if(last == NULL)
		{
		int a_test = 0;
		while(a_test < m + 2)
			{
			cerr << "entering the check: "; 
			FindMyNode(tab[a_test]);
			cerr << "exiting the check.\n";
			a_test++;
			}

		//OUTPUT ALL THE SYSTEM THERE
		exit(1);
		}
	#endif

	//The root of the sproot
	N->CountBranches++;
	N->branch = (Vessel_Graph**) realloc(N->branch, sizeof(Vessel_Graph*) * N->CountBranches);
	N->edge = (Vessel_Unit**) realloc(N->edge, sizeof(Vessel_Unit*) * N->CountBranches);
	N->branch[N->CountBranches - 1] = VGtab[0];
	
	VGtab[0]->CountBranches = 2;
	VGtab[0]->branch = (Vessel_Graph**) malloc( sizeof(Vessel_Graph*) * 2);
	VGtab[0]->edge = (Vessel_Unit**) malloc( sizeof(Vessel_Unit*) * 2);
	VGtab[0]->branch[0] = N;
	VGtab[0]->branch[1] = VGtab[1];

	MakeOneEdge(N,VGtab[0]);

	for(i = 1; i < m-1; i++)
		{
		VGtab[i]->CountBranches = 2;
		VGtab[i]->branch = (Vessel_Graph**) malloc(sizeof(Vessel_Graph*) * 2);
		VGtab[i]->edge = (Vessel_Unit**) malloc(sizeof(Vessel_Unit*) * 2);
		VGtab[i]->branch[0] = VGtab[i-1];
		VGtab[i]->branch[1] = VGtab[i+1];
		
		MakeOneEdge(VGtab[i - 1],VGtab[i]);
		}
	

	 
	//The end of the sprout
	VGtab[m-1]->CountBranches = 2;
	VGtab[m-1]->branch = (Vessel_Graph**) malloc( sizeof(Vessel_Graph*) * 2);
	VGtab[m-1]->edge = (Vessel_Unit**) malloc( sizeof(Vessel_Unit*) * 2);
	VGtab[m-1]->branch[0] = VGtab[m-2];
	VGtab[m-1]->branch[1] = last;	

	MakeOneEdge(VGtab[m - 2],VGtab[m - 1]);

	last->CountBranches++;
	last->branch = (Vessel_Graph**) realloc(last->branch, sizeof(Vessel_Graph*) * last->CountBranches);
	last->edge = (Vessel_Unit**) realloc(last->edge, sizeof(Vessel_Unit*) * last->CountBranches);
	last->branch[last->CountBranches - 1] = VGtab[m-1];

	MakeOneEdge(VGtab[m - 1],last);


	#ifdef __myDEBUG__
	DoTheSprout_counter++;

	#endif




     return false;
     }

bool Angiogenesis::KillConnection(Vessel_Unit *Vu)
     {

	#ifdef __myDEBUG__
	KillConnection_counter++;
	#endif

     //Kill slice between the nodes
     

     int k = Vu->node1->FindBranchNumber(Vu->node2);
     Vu->node1->CountBranches--;
     Vessel_Graph *temp = Vu->node1->branch[k];
     Vu->node1->branch[k] = Vu->node1->branch[Vu->node1->CountBranches];
     Vu->node1->branch[Vu->node1->CountBranches] = temp;
     Vessel_Unit *tamp = Vu->node1->edge[k];
     Vu->node1->edge[k] = Vu->node1->edge[Vu->node1->CountBranches];
     Vu->node1->edge[Vu->node1->CountBranches] = tamp;
     
     
     k = Vu->node2->FindBranchNumber(Vu->node1);
     Vu->node2->CountBranches--;
     temp = Vu->node2->branch[k];
     Vu->node2->branch[k] = Vu->node2->branch[Vu->node2->CountBranches];
     Vu->node2->branch[Vu->node2->CountBranches] = temp;
     tamp = Vu->node2->edge[k];
     Vu->node2->edge[k] = Vu->node2->edge[Vu->node2->CountBranches];
     Vu->node2->edge[Vu->node2->CountBranches] = tamp;

	
     KillSlice(Vu);
     
     return true;
     }

bool ProbaToSprout(double dt)
	{
	if(myRand() < (dt / EC_PROLIFERATION_TIME))
		return true;

	return false;
	}

bool ProbaToCollapse(double dt)
	{
	if(myRand() < (dt / 10*EC_PROLIFERATION_TIME))
		return true;

	return false;
	}

bool ProbaToRegress()
	{
	if(myRand() < 0.5)
		return true;

	return false;
	}

bool SurroundingNode(Vessel_Graph * Vg)
	{

	int nb_tumor_cells = 0;
	//If surrounding by 80% of tumor cells - RIEGER - Pysical Review Letter
	for(int j = 0; j < GetVoronoiCell( Vg->p_agent )->countNeighborCells; j++)
		{
		if( GetVoronoiCell( Vg->p_agent )->neighborCells[j]->getState() == ACTIVE    )
			{nb_tumor_cells++;}

		if( GetVoronoiCell( Vg->p_agent )->neighborCells[j]->getState() == NONACTIVE )
			{nb_tumor_cells++;}
		}

	if( (double)nb_tumor_cells/ (double)GetVoronoiCell( Vg->p_agent )->countNeighborCells  < SURROUNDING_PERCENTAGE )
		{
		return false;
		}


	return true;
	}

bool SurroundingSlice(Vessel_Unit * Vu)
	{
	
	int j,nb_tumor_cells = 0;
	//If surrounding by 80% of tumor cells - RIEGER - Pysical Review Letter
	for(j = 0; j < GetVoronoiCell( Vu->node1->p_agent )->countNeighborCells; j++)
		{
		if( GetVoronoiCell( Vu->node1->p_agent )->neighborCells[j]->getState() == ACTIVE    )
			{nb_tumor_cells++;}

		if( GetVoronoiCell( Vu->node1->p_agent )->neighborCells[j]->getState() == NONACTIVE )
			{nb_tumor_cells++;}
		}

	if( (double)nb_tumor_cells/ (double)GetVoronoiCell( Vu->node1->p_agent )->countNeighborCells  < SURROUNDING_PERCENTAGE )
		{
		return false;
		}
	
	nb_tumor_cells = 0;
	for(j = 0; j < GetVoronoiCell( Vu->node2->p_agent )->countNeighborCells; j++)
		{
		if( GetVoronoiCell( Vu->node2->p_agent )->neighborCells[j]->getState() == ACTIVE    )
			{nb_tumor_cells++;}

		if( GetVoronoiCell( Vu->node2->p_agent )->neighborCells[j]->getState() == NONACTIVE )
			{nb_tumor_cells++;}
		}

	if( (double)nb_tumor_cells/ (double)GetVoronoiCell( Vu->node2->p_agent )->countNeighborCells  < SURROUNDING_PERCENTAGE )
		{
		return false;
		}


	return true;
	}

bool Angiogenesis::DoTheVesselCollapse(double dt)
     	{
     	bool modification = false;
     	
     	for(int i = nbrSlices - 1; i >= 0; i--)
		{
		if(Slice[i]->shear < MINSHEAR && ProbaToCollapse(dt) && SurroundingSlice(Slice[i]) )
			{

			modification = KillConnection(Slice[i]);	
			}
		}
	
     	return modification;
	}

bool Angiogenesis::Update_Network(double dt, ActionTree *p_list, VoronoiDiagram* voronoiDiagram)
	{


	bool modification = false;
	
	int i;
     	for(i = 0; i < nbrSlices; i++)
		{			
		if(Slice[i]->flow <= MICRO_EPSILON)
			{
			Slice[i]->circulation = false;
			}
		else
			{
			Slice[i]->circulation = true;
			}
		}

	for(i = nbrNodes - 1; i >= 0; i--) //The new nodes are added at the end.
		{

		if(Node[i]->CountBranches == 0)
			{

			#ifdef __myDEBUG__
			Ko_counter++;
			#endif
			removeVesselAndUpdateSurrounding( p_list, angioAgentList, GetVoronoiCell( Node[i]->p_agent ), voronoiDiagram);
			KillNode(Node[i]);
			}
		else
			{
			if(GetVoronoiCell( Node[i]->p_agent )->growthfactors > ANGIOGENESIS_GROWTHFACTOR_THRESHOLD)
				{
				if(ProbaToSprout(dt))
					{
					if(SurroundingNode(Node[i]))
               					{modification |= DoTheDilation(Node[i]);}

               				modification |= DoTheSprout(Node[i],p_list,voronoiDiagram);
                   			}
				}
			
			//else  //that condition tells that a node with Gf > Threshold can't regress
				{
				if(GetVoronoiCell( Node[i]->p_agent )->oxygen < MINOXYGEN)
					{
					if(ProbaToRegress())
               					{
               					modification = true;

						#ifdef __myDEBUG__
						Regression_counter++;
						#endif
						removeVesselAndUpdateSurrounding( p_list, angioAgentList, GetVoronoiCell( Node[i]->p_agent ), voronoiDiagram);
               					KillNode(Node[i]);
              					}
					}
				}
			}
		}

	modification |= DoTheVesselCollapse(dt);

	return modification;
	}
	
void Angiogenesis::Holocaust()
     {
     int i;
     for(i = nbrNodes - 1; i >=0 ; i--)
           {
           KillNode(Node[i]);
           }
     for(i = nbrSlices - 1; i >=0 ; i--)
           {
           KillSlice(Slice[i]);
           }
     }	


void Angiogenesis::GiveMeThePressure(double epsilon)	
	{
	int i;
	double TempPressure;

	bool Convergence = false;
	while(Convergence == false)
		{
		Convergence = true;
		for(i = 0; i < nbrNodes; i++)
			{
			TempPressure = Node[i]->pressure;
			Node[i]->PressureComputation();
				
			//Check the convergence, it has to be good for each node.
			if(fabs(Node[i]->pressure - TempPressure) > epsilon)
			{Convergence = false;}		
                	}
		}
	}
void Angiogenesis::GiveMeTheFlow()
	{
	int i,j;
	double diffPressure,out;
	for(i = 0; i < nbrSlices; i++)
		{
           	Slice[i]->Compute_flow();
           	}
	//Compute the Outflow
	for(i = 0; i < nbrNodes; i++)
           	{
		out = 0;
		for(j = 0; j < Node[i]->CountBranches; j++)
			{
			diffPressure = Node[i]->pressure - Node[i]->branch[j]->pressure;
				
			if(diffPressure > 0)
				{
				out += Node[i]->edge[j]->flow;
				}
			Node[i]->outflow = out;
			}
        	}
	//Compute the Proportion of the flow from a node to another
	for(i = 0; i < nbrNodes; i++)
           	{
		out = 0;
		for(j = 0; j < Node[i]->CountBranches; j++)
			{
			diffPressure = Node[i]->pressure - Node[i]->branch[j]->pressure;	
			if(diffPressure > 0)
				{
				Node[i]->edge[j]->flowproportion = Node[i]->edge[j]->flow / Node[i]->outflow;
				}
			}
        	}
     	}
   
void Angiogenesis::GiveMeTheSpeed()
     {
     int i;
     for(i = 0; i < nbrSlices; i++)
	{Slice[i]->Compute_velocity();}
     }
void Angiogenesis::GiveMeTheShearStress()
     {
     int i;
     for(i = 0; i < nbrSlices; i++)
	{Slice[i]->Compute_shearstress();}
     }

void Angiogenesis::GiveMeTheBlood()
	{
        ////Find the pressure, the flow, the speed and the shearstress in each node or slice
	GiveMeThePressure(PRESSURE_EPSILON);
        GiveMeTheFlow();
	GiveMeTheSpeed();
	GiveMeTheShearStress();
        }

void Angiogenesis::Update_Graph(Agent *agentmoved)
{
	Vessel_Graph *N = FindMyNode( agentmoved );

	bool belongstoneighborhood;

	for(int i = 0; i < N->CountBranches; ++i)
		{
		belongstoneighborhood = false;
		for(int j = 0; j < GetVoronoiCell( agentmoved )->countNeighborCells; j++)
			{
			if( GetVoronoiCell( N->branch[i]->p_agent ) == GetVoronoiCell( agentmoved )->neighborCells[j])
				{ belongstoneighborhood |= true;}
			}
		if( belongstoneighborhood == false)
			{
			KillConnection( N->edge[i] );
			}
		}
}
