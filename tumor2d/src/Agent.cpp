#include <stdio.h>
#include <float.h>
#include <limits.h>

#include "Agent.h"
#include "CellProcesses.h"



int CountCellsPerVoronoiCell = 1;

float Agent::ReentranceRate = 0;

double Agent::NICK_G_CRITICAL_GLU = 0.068049;
double Agent::NICK_G_CRITICAL_OXY = 0.030789;
double Agent::NICK_G_MIN = (14.727716e-17);
double Agent::NICK_G_MAX = (53.672035e-17);


//Oxygen Uptake
double Agent::NICK_O_CRITICAL_OXY = 0.030752;
double Agent::NICK_O_CRITICAL_GLU = 0.100326;
double Agent::NICK_O_MIN = (10.171515e-17);
double Agent::NICK_O_MAX = (22.800151e-17);

bool Agent::USE_GRAVITY	= false;

float Agent::WASTE_UPTAKE = 0;
float Agent::WASTE_THRESHOLD_QUIESCENCE = FLT_MAX;
//float Agent::WASTE_THRESHOLD_QUIESCENCE_INTRACELLULAR = FLT_MAX;
float Agent::WASTE_THRESHOLD_SLOWED_GROWTH = FLT_MAX;
int   Agent::WASTE_INTOXICATED_CELL_CYCLES = INT_MAX;

double Agent::getOxygen()
{
	double oxyConc = 0.;
	
	for( int i=0; i<this->countLocations; i++)
		oxyConc += this->location[i]->oxygen;

	return oxyConc / (double) this->countLocations;
}

double Agent::getGlucose()
{
	double gluConc = 0.;
	
	for( int i=0; i<this->countLocations; i++)
		gluConc += this->location[i]->glucose;

	return gluConc / (double) this->countLocations;
}

int Agent::isConnected()
{
		int ii;
		int stackLength = 1;
		int stackPosition = 0;
		VoronoiCell *actualLocation;// = agentArray->agents[i]->location[0];
		VoronoiCell *locationStack[this->countLocations];
		locationStack[0]  = this->location[0];
		
		do{
			actualLocation = locationStack[stackPosition];
			stackPosition++;
			//int iii;
			for( ii=0; ii<actualLocation->countNeighborCells; ii++)
			{
				// neighbor is occupied by same neighbor?
				if( GetAgent(actualLocation->neighborCells[ii]) == this){
					// is already in stack?
					int iv;
					for( iv=0; iv<stackLength && locationStack[iv]!=actualLocation->neighborCells[ii]; iv++);
					if( iv==stackLength){
						locationStack[stackLength++] = actualLocation->neighborCells[ii];
					}
				}
			}
		}while( stackLength != stackPosition);

	return this->countLocations == stackLength;
}


void Agent::validate()
{
	int i;
	int countFreeNeighbors = 0;

	for( i=0; i<this->countLocations; i++){
		countFreeNeighbors += this->location[i]->countFreeNeighborCells + this->location[i]->countFreeExtendedNeighborCells;
		if( this->location[i]->countFreeNeighborCells<0 || this->location[i]->countFreeExtendedNeighborCells<0){
			fprintf( stderr, "ERROR: Agent %i -> location %i (%s) has free number of neighbors below zero!\n", this->index, this->location[i]->index, cellTypeToString( this->state));
			fprintf( stderr, "%i free neighbors and %i free extended neighbors\n", this->location[i]->countFreeNeighborCells, this->location[i]->countFreeExtendedNeighborCells);
			exit( 0);
		}
	}


	// ACTUALIZE STATE
	switch( this->state){
		case FREE:
			break;

		case NONACTIVE:
			if( countFreeNeighbors>0){
				//agent has free neighbors
				fprintf( stderr, "ERROR: Agent %i (%s) should be ACTIVE!\n", this->index, cellTypeToString( this->state));
				fprintf( stderr, "%i free neighbors and extended neighbors\n", countFreeNeighbors);
				exit(0);
			}
			break;

		case ACTIVE:
			//for( i=0; i<this->countLocations && this->location[i]->countFreeNeighborCells==0 && this->location[i]->countFreeExtendedNeighborCells==0; i++) 
			//	countFreeNeighbors += this->location[i]->countFreeNeighborCells + this->location[i]->countFreeExtendedNeighborCells;
			//if( i==this->countLocations){
			if( countFreeNeighbors==0){
				//agent has no free neighbors
				fprintf( stderr, "ERROR: Agent %i (%s) should be NONACTIVE!\n", this->index, cellTypeToString( this->state));
				fprintf( stderr, "%i free neighbors and extended neighbors\n", countFreeNeighbors);
				exit(0);
			}
			break;
	}

	// ACTUALIZE ACTIONS
	switch( this->state){
		case FREE:
			if( this->actionsInitialized)
			{
#if USE_DIVISION
			// NO DIVISION
			if( this->actions[INDEX_DIVISION]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains DIVISION!\n", this->index, cellTypeToString( this->state));
			exit(0);}
#endif
#if USE_NECROSIS
			// NO NECROSIS
			if( this->actions[INDEX_NECROSIS]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains NECROSIS!\n", this->index, cellTypeToString( this->state));
			exit(0);}
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains LYSIS!\n", this->index, cellTypeToString( this->state));
			exit(0);}
#endif
			}
			break;

		case ACTIVE:
#if USE_DIVISION
			// -> DIVISION
			if( this->actions[INDEX_DIVISION]->next==NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) doesn't contain DIVISION!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
#if USE_NECROSIS
			// -> NECROSIS
			if( this->actions[INDEX_NECROSIS]->next==NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) doesn't contain NECROSIS!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains LYSIS!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
			break;

		case NONACTIVE:
#if USE_DIVISION
			// NO DIVISION
			if( this->actions[INDEX_DIVISION]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains DIVISION!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
#if USE_NECROSIS
			// -> NECROSIS
			if( this->actions[INDEX_NECROSIS]->next==NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) doesn't contain NECROSIS!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains LYSIS!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
			break;

		case NECROTIC:
#if USE_DIVISION
			// NO DIVISION
			if( this->actions[INDEX_DIVISION]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains DIVISION!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
#if USE_NECROSIS
			// NO NECROSIS
			if( this->actions[INDEX_NECROSIS]->next!=NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) contains NECROSIS!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
#if USE_LYSIS
			// -> LYSIS
			if( this->actions[INDEX_LYSIS]->next==NULL){
				fprintf( stderr, "ERROR: Agent %i (%s) doesn't contain LYSIS!\n", this->index, cellTypeToString( this->state));
exit(0);}
#endif
			break;
	}	

}

#if (USE_ACTION_TREE == TRUE)
void Agent::actualize( ActionTree *actionTree)
{
	int i;

	//fprintf( stderr, "Agent::actualize( ActionList *actionList)\n");

	//============ ACTUALIZE STATE ================//

	//fprintf( stderr, "agent was %s\n", cellTypeToString( this->state));

	switch( this->state){
		case COMPARTMENT:
			break;
			
		case FREE:
			break;

		case NONACTIVE:
			for( i=0; i<this->countLocations && this->location[i]->countFreeNeighborCells==0 && this->location[i]->countFreeExtendedNeighborCells==0; i++);
			if( i<this->countLocations){
				//agent has free neighbors
				//fprintf( stderr, "set agent %i to ACTIVE\n", this->index);
				this->state = ACTIVE;
				//if()
#ifdef REFINE
				GetAgent( this->location[0]->coarseParent)->countActive++;
				GetAgent( this->location[0]->coarseParent)->countNonactive--;
#endif
			}
			break;

		case ACTIVE:
			int countFreeNeighbors = 0;
			for( i=0; i<this->countLocations; i++)
				 countFreeNeighbors += this->location[i]->countFreeNeighborCells + this->location[i]->countFreeExtendedNeighborCells;
			if( !VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD && !VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD && !countFreeNeighbors){
				//agent has no free neighbors
				this->state = NONACTIVE;
				//fprintf( stderr, "set agent %i (on cell %i) to NONACTIVE (free:%i, Xtended free:%i)\n", this->index, this->location[0]->index, this->location[0]->countFreeNeighborCells, this->location[0]->countFreeExtendedNeighborCells);
#ifdef REFINE
				GetAgent( this->location[0]->coarseParent)->countActive--;
				GetAgent( this->location[0]->coarseParent)->countNonactive++;
#endif
			}
			/*for( i=0; i<this->countLocations && (this->location[i]->countFreeNeighborCells!=0 || this->location[i]->countFreeExtendedNeighborCells!=0); i++);
			if( i<this->countLocations){
				//agent has no free neighbors
				this->state = NONACTIVE;
			}*/
			break;
	}

	//fprintf( stderr, "agent is %s\n", cellTypeToString( this->state));

	//============ ACTUALIZE ACTIONS ================//

	if( !this->actionsInitialized)
		initCellActions( this);

	switch( this->state){

		// COMPARTMENT //

		case COMPARTMENT:
			{
#if USE_DIVISION
			if( !this->isDividing()){
				// NO DIVISION
				if( this->actions[INDEX_DIVISION]->top!=NULL){
					actionTree->deleteAction( this->actions[INDEX_DIVISION]);
					//fprintf( stderr, "Delete Division\n");
					//actionTree->getDepth( actionTree->root);	
				}
			}else{
				// -> DIVISION
				if( this->actions[INDEX_DIVISION]->top==NULL){
					actionTree->addAction( this->actions[INDEX_DIVISION]);
					//fprintf( stderr, "Add Division\n");
					//actionTree->getDepth( actionTree->root);	
				}
			}
#endif
#if USE_GROWTH
			if( !this->isGrowing()){
				// NO GROWTH
				if( this->actions[INDEX_GROWTH]->top!=NULL){
					actionTree->deleteAction( this->actions[INDEX_GROWTH]);
					//fprintf( stderr, "Delete Growth\n");
					//actionTree->getDepth( actionTree->root);	
				}
			}else{
				// -> GROWTH
				if( this->actions[INDEX_GROWTH]->top==NULL){
					actionTree->addAction( this->actions[INDEX_GROWTH]);
					//fprintf( stderr, "Add Growth\n");
					//actionTree->getDepth( actionTree->root);	
				}
			}
#endif
#if USE_MIGRATION
			/*int countFreeNeighbors = 0;
			for( int i=0; i<this->location[0]->countNeighborCells; i++)
				if(this->location[0]->neighborCells[i]->isFree())
					countFreeNeighbors++;

			if( !this->isGrowing() || countFreeNeighbors==0){*/
			if( !this->isGrowing() || this->countFreeNeighbors()==0){
				// NO MIGRATION
				if( this->actions[INDEX_MIGRATION]->top!=NULL){
					actionTree->deleteAction( this->actions[INDEX_MIGRATION]);
				}
			}else{
				// -> MIGRATION
				if( this->actions[INDEX_MIGRATION]->top==NULL){
					actionTree->addAction( this->actions[INDEX_MIGRATION]);
				}
			}
#endif
#if USE_NECROSIS
			if( !this->isDying()){
				// NO NECROSIS
				if( this->actions[INDEX_NECROSIS]->top!=NULL)
					actionTree->deleteAction( this->actions[INDEX_NECROSIS]);
			}else{
				// -> NECROSIS
				if( this->actions[INDEX_NECROSIS]->top==NULL){
					//fprintf( stderr, "Add Necrosis\n");
					actionTree->addAction( this->actions[INDEX_NECROSIS]);
				}
			}
#endif
#if USE_LYSIS
			if( !this->isLysing()){
				// NO LYSIS
				if( this->actions[INDEX_LYSIS]->top!=NULL)
					actionTree->deleteAction( this->actions[INDEX_LYSIS]);
			}else{
				// -> LYSIS
				if( this->actions[INDEX_LYSIS]->top==NULL)
					actionTree->addAction( this->actions[INDEX_LYSIS]);
			}
#endif
			}
			break;


		// FREE //

		case FREE:
			if( this->actionsInitialized)
			{
#if USE_DIVISION
			// NO DIVISION
			if( this->actions[INDEX_DIVISION]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_DIVISION]);
#endif
#if USE_GROWTH
			// NO GROWTH
			if( this->actions[INDEX_GROWTH]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_GROWTH]);
#endif
#if USE_MIGRATION
			// NO MIGRATION
			if( this->actions[INDEX_MIGRATION]->top!=NULL){
				actionTree->deleteAction( this->actions[INDEX_MIGRATION]);
			}
#endif
#if USE_NECROSIS
			// NO NECROSIS
			if( this->actions[INDEX_NECROSIS]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_NECROSIS]);
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_LYSIS]);
#endif
			}
			break;


		// ACTIVE //
//actionList->getDepth( actionList->root);
//							fprintf( stderr, "after selected_action->originalCell->actualize( actionList);\n ");
							
		case ACTIVE:

#if USE_GROWTH
			if( this->countLocations < MAX_SUBCELLULAR_COMPONENTS){
				//fprintf( stderr, "Agent %i is ACTIVE and has %i/%i components!\n", this->index, this->countLocations, 2*SUBCELLULAR_COMPONENTS);
				// -> GROWTH
				if( this->actions[INDEX_GROWTH]->top==NULL /*&&  this->divide == 1.*/){
	//fprintf( stderr, "DIVISION (rate=%lf)\n", this->actions[INDEX_DIVISION]->rate);
					actionTree->addAction( this->actions[INDEX_GROWTH]);
					//fprintf( stderr, "Add GROWTH (rate=%lf)\n", this->actions[INDEX_GROWTH]->rate);
	//fprintf( stderr, "DIVISION (rate=%lf)\n", this->actions[INDEX_DIVISION]->rate);
					//actionTree->getDepth(actionTree->root);
				}
#if USE_DIVISION
				if( this->countLocations >= MIN_SUBCELLULAR_COMPONENTS*2){
					// -> DIVISION
					if( this->actions[INDEX_DIVISION]->top==NULL /*&&  this->divide == 1.*/){
						actionTree->addAction( this->actions[INDEX_DIVISION]);
						//fprintf( stderr, "Add DIVISION (rate=%lf)\n", this->actions[INDEX_DIVISION]->rate);
					}
				}else{
					// NO DIVISION
					if( this->actions[INDEX_DIVISION]->top!=NULL){
						//fprintf( stderr, "Delete DIVISION (rate=%lf)\n", this->actions[INDEX_DIVISION]->rate);
						actionTree->deleteAction( this->actions[INDEX_DIVISION]);
						//actionTree->getDepth(actionTree->root);
					}
				}
#endif
			}else{
				//fprintf( stderr, "Agent %i is ACTIVE and has %i/%i components!\n", this->index, this->countLocations, 2*SUBCELLULAR_COMPONENTS);
#endif
#if USE_DIVISION
				// -> DIVISION
				if( this->actions[INDEX_DIVISION]->top==NULL /*&&  this->divide == 1.*/){
					actionTree->addAction( this->actions[INDEX_DIVISION]);
					//fprintf( stderr, "Add DIVISION (rate=%lf)\n", this->actions[INDEX_DIVISION]->rate);
				}

#endif
#if USE_GROWTH
				// NO GROWTH
				if( this->actions[INDEX_GROWTH]->top!=NULL){
					//fprintf( stderr, "Delete GROWTH (rate=%lf)\n", this->actions[INDEX_GROWTH]->rate);
					actionTree->deleteAction( this->actions[INDEX_GROWTH]);
				}

			}
#endif

#if USE_MIGRATION
			// -> MIGRATION
			if( this->actions[INDEX_MIGRATION]->top==NULL){
				actionTree->addAction( this->actions[INDEX_MIGRATION]);
				//fprintf( stderr, "Add Necrosis (rate=%lf)\n", this->actions[INDEX_NECROSIS]->rate);
			}
#endif

#if USE_NECROSIS
			// -> NECROSIS
			if( this->actions[INDEX_NECROSIS]->top==NULL){
				actionTree->addAction( this->actions[INDEX_NECROSIS]);
				//fprintf( stderr, "Add Necrosis (rate=%lf)\n", this->actions[INDEX_NECROSIS]->rate);
			}
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_LYSIS]);
#endif
			break;


		// NONACTIVE //

		case QUIESCENT:
		case NONACTIVE:
#if USE_GROWTH
			// NO GROWTH
			if( this->actions[INDEX_GROWTH]->top!=NULL){
				//fprintf( stderr, "Delete GROWTH (rate=%lf) address: %p\n", this->actions[INDEX_GROWTH]->rate, this->actions[INDEX_GROWTH]);
				actionTree->deleteAction( this->actions[INDEX_GROWTH]);
			}else{
				//fprintf( stderr, "GROWTH (rate=%lf) alreadty deleted?!\n", this->actions[INDEX_GROWTH]->rate);
			}
#endif
#if USE_DIVISION
			// -> DIVISION
			if( this->countLocations >= 2*MIN_SUBCELLULAR_COMPONENTS){
				if( this->actions[INDEX_DIVISION]->top==NULL){
					//fprintf( stderr, "Add DIVISION (rate=%lf)\n", this->actions[INDEX_DIVISION]->rate);
					actionTree->addAction( this->actions[INDEX_DIVISION]);
				}
			}
			// NO DIVISION
			else{
				if( this->actions[INDEX_DIVISION]->top!=NULL){
					//fprintf( stderr, "Delete DIVISION (rate=%lf)\n", this->actions[INDEX_DIVISION]->rate);
					actionTree->deleteAction( this->actions[INDEX_DIVISION]);
				}
			}
#endif
#if USE_MIGRATION
			// NO MIGRATION
			if( this->actions[INDEX_MIGRATION]->top!=NULL){
				actionTree->deleteAction( this->actions[INDEX_MIGRATION]);
			}
#endif
#if USE_NECROSIS
			// -> NECROSIS
			if( this->actions[INDEX_NECROSIS]->top==NULL){
				actionTree->addAction( this->actions[INDEX_NECROSIS]);
				//fprintf( stderr, "Add Necrosis (%s)\n", "NONACTIVE	");
			}
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_LYSIS]);
#endif
			break;

		case NECROTIC:
#if USE_DIVISION
			// NO DIVISION
			if( this->actions[INDEX_DIVISION]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_DIVISION]);
#endif
#if USE_GROWTH
			// NO GROWTH
			if( this->actions[INDEX_GROWTH]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_GROWTH]);
#endif
#if USE_MIGRATION
			// NO MIGRATION
			if( this->actions[INDEX_MIGRATION]->top!=NULL){
				actionTree->deleteAction( this->actions[INDEX_MIGRATION]);
			}
#endif
#if USE_NECROSIS
			// NO NECROSIS
			if( this->actions[INDEX_NECROSIS]->top!=NULL)
				actionTree->deleteAction( this->actions[INDEX_NECROSIS]);
#endif
#if USE_LYSIS
			// -> LYSIS
			if( this->actions[INDEX_LYSIS]->top==NULL)
				actionTree->addAction( this->actions[INDEX_LYSIS]);
#endif
			break;
	}	
}
#else
void Agent::actualize( ActionList *actionList)
{
	int i;

	//fprintf( stderr, "Agent::actualize( ActionList *actionList)\n");

	//============ ACTUALIZE STATE ================//

	switch( this->state){
		case FREE:
			break;

		case NONACTIVE:
			for( i=0; i<this->countLocations && this->location[i]->countFreeNeighborCells==0 && this->location[i]->countFreeExtendedNeighborCells==0; i++);
			if( i<this->countLocations){
				//agent has free neighbors
				//fprintf( stderr, "set agent %i to ACTIVE\n", this->index);
				this->state = ACTIVE;
			}
			break;

		case ACTIVE:
			for( i=0; i<this->countLocations && (this->location[i]->countFreeNeighborCells!=0 || this->location[i]->countFreeExtendedNeighborCells!=0); i++);
			if( i<this->countLocations){
				//agent has no free neighbors
				this->state = NONACTIVE;
			}
			break;
	}



	//============ ACTUALIZE ACTIONS ================//

	switch( this->state){

		// FREE //

		case FREE:
			if( this->actionsInitialized)
			{
#if USE_DIVISION
			// NO DIVISION
			if( this->actions[INDEX_DIVISION]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_DIVISION]);
#endif
#if USE_GROWTH
			// NO GROWTH
			if( this->actions[INDEX_GROWTH]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_GROWTH]);
#endif
#if USE_NECROSIS
			// NO NECROSIS
			if( this->actions[INDEX_NECROSIS]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_NECROSIS]);
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_LYSIS]);
#endif
			}
			break;


		// ACTIVE //

		case ACTIVE:

#if USE_GROWTH
			if( this->countLocations < MAX_SUBCELLULAR_COMPONENTS){
				//fprintf( stderr, "Agent %i is ACTIVE and has %i/%i components!\n", this->index, this->countLocations, 2*SUBCELLULAR_COMPONENTS);
				// -> GROWTH
				if( this->actions[INDEX_GROWTH]->next==NULL)
					addAction( actionList, this->actions[INDEX_GROWTH]);
#if USE_DIVISION
				// NO DIVISION
				if( this->actions[INDEX_DIVISION]->next!=NULL)
					deleteAction( actionList, this->actions[INDEX_DIVISION]);
#endif
			}else{
#endif
#if USE_DIVISION
				// -> DIVISION
				if( this->actions[INDEX_DIVISION]->next==NULL)
					addAction( actionList, this->actions[INDEX_DIVISION]);
#endif
#if USE_GROWTH
				// NO GROWTH
				if( this->actions[INDEX_GROWTH]->next!=NULL)
					deleteAction( actionList, this->actions[INDEX_GROWTH]);
			}
#endif

#if USE_NECROSIS
			// -> NECROSIS
			if( this->actions[INDEX_NECROSIS]->next==NULL)
				addAction( actionList, this->actions[INDEX_NECROSIS]);
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_LYSIS]);
#endif
			break;


		// NONACTIVE //

		case NONACTIVE:
#if USE_GROWTH
			// NO GROWTH
			if( this->actions[INDEX_GROWTH]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_GROWTH]);
#endif
#if USE_DIVISION
			// -> DIVISION
			if( this->countLocations >= 2*MIN_SUBCELLULAR_COMPONENTS){
				if( this->actions[INDEX_DIVISION]->next==NULL)
					addAction( actionList, this->actions[INDEX_DIVISION]);
			}
			// NO DIVISION
			else{
				if( this->actions[INDEX_DIVISION]->next!=NULL)
					deleteAction( actionList, this->actions[INDEX_DIVISION]);
			}
#endif
#if USE_NECROSIS
			// -> NECROSIS
			if( this->actions[INDEX_NECROSIS]->next==NULL)
				addAction( actionList, this->actions[INDEX_NECROSIS]);
#endif
#if USE_LYSIS
			// NO LYSIS
			if( this->actions[INDEX_LYSIS]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_LYSIS]);
#endif
			break;

		case NECROTIC:
#if USE_DIVISION
			// NO DIVISION
			if( this->actions[INDEX_DIVISION]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_DIVISION]);
#endif
#if USE_GROWTH
			// NO GROWTH
			if( this->actions[INDEX_GROWTH]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_GROWTH]);
#endif
#if USE_NECROSIS
			// NO NECROSIS
			if( this->actions[INDEX_NECROSIS]->next!=NULL)
				deleteAction( actionList, this->actions[INDEX_NECROSIS]);
#endif
#if USE_LYSIS
			// -> LYSIS
			if( this->actions[INDEX_LYSIS]->next==NULL)
				addAction( actionList, this->actions[INDEX_LYSIS]);
#endif
			break;
	}	
}
#endif

void Agent::print()
{
	fprintf( stderr, "index=%i, countLocations=%i\n", this->index, countLocations);
}


void AgentList::print()
{
	int i;
	
	for( i=0; i<this->countAgents; i++){
		this->agents[i]->print();
	}
}


void AgentList::printActiveAgents()
{
	int i;
	
	for( i=0; i<this->countActiveAgents; i++){
		this->agents[i]->print();
	}
}


AgentList* AgentList::newAgentList( int countNewAgents)
{
	int i;
	AgentList *agentList = (AgentList *) malloc( sizeof(AgentList));
	
	agentList->agents = (Agent **) calloc( countNewAgents, sizeof(Agent *));
	for( i=0; i<countNewAgents; i++){
		agentList->agents[i] = new Agent();
		agentList->agents[i]->index = i;
	}
	
	agentList->countActiveAgents = 0;
	agentList->countAgents = countNewAgents;
	
	return agentList;
}


Agent* AgentList::activateAgent()
{
	if( countActiveAgents==countAgents){
		this->agents = (Agent **) realloc( this->agents, (countAgents+1000) * sizeof(Agent *));
		for( int i=countAgents; i<countAgents+1000; i++){
			this->agents[i] = new Agent();
			this->agents[i]->index = i;
		}
		countAgents += 1000;
	}

	if( countActiveAgents<countAgents){
		// remove old action structure
		if( this->agents[ countActiveAgents]->actionsInitialized){
			for( int k=0; k<INDEX_MIGRATION + USE_MIGRATION; k++)
				free( this->agents[ countActiveAgents]->actions[k]);
#if USE_MIGRATION
			//for( int k=0; k < GetVoronoiCell( this->agents[ countActiveAgents])->countNeighborCells; k++)
			//	free( this->agents[ countActiveAgents]->actions[k+INDEX_MIGRATION]);
#endif
			free( this->agents[ countActiveAgents]->actions);
		
			this->agents[ countActiveAgents]->actionsInitialized = FALSE;
		}

		this->agents[ countActiveAgents]->initAgent();
		this->agents[ countActiveAgents]->cellCount = 0;
		this->agents[ countActiveAgents]->maxCellCount = CountCellsPerVoronoiCell;
		this->agents[ countActiveAgents]->growingTumorCellCount = 0;
		this->agents[ countActiveAgents]->dividingTumorCellCount = 0;
		this->agents[ countActiveAgents]->divide = 1;

		return this->agents[ countActiveAgents++];
	}else{
		fprintf( stderr, "ERROR in Agent* AgentList::activateAgent()\nAgentList is already full!\n");
		exit( 0);
	}
	
}

void AgentList::deactivateAgent( Agent* agent)
{
	int tempIndex    = agent->index;
	//Agent* tempAgent = agent;
	
	// decrease number of active agents
	--(this->countActiveAgents);
	
	// replace deactivated agent by another active one
	this->agents[ tempIndex]        = this->agents[ this->countActiveAgents];
	this->agents[ tempIndex]->index = tempIndex;
	
	// move deactivated agent 
	this->agents[ this->countActiveAgents] = agent;
	agent->index = this->countActiveAgents;
}


void AgentList::deleteAgentList()
{}



/*AgentArray * newAgentArray( VoronoiDiagram * voronoiDiagram)
{
	int i;
	AgentArray * cellDataArray;

	cellDataArray = (AgentArray *) malloc( sizeof(AgentArray));

	cellDataArray->countAgent = voronoiDiagram->countVoronoiCells;
	cellDataArray->maxAgent   = voronoiDiagram->maxVoronoiCells;
	cellDataArray->cellData      = (Agent **) calloc( sizeof(Agent*) , cellDataArray->maxAgent);
	for( i=0; i<cellDataArray->maxAgent; i++)
		cellDataArray->cellData[i] = NULL;

	cellDataArray->countSetAgent = 0;
	cellDataArray->maxSetAgent = 0;
	cellDataArray->setAgent = NULL;

	return cellDataArray;
}*/
void Agent::attach( VoronoiCell * correspondingVoronoiCell)
{
	AddLocation( this, correspondingVoronoiCell);
	correspondingVoronoiCell->agent = this;

	//this->countDirectFreeNeighbors    = correspondingVoronoiCell->countNeighborCells;	
	//this->countReachableFreeNeighbors = 0;
	
	
}
/****************************************************************************/


void Agent::detach()
{
	int i;
	
	for( i=0; i<this->countLocations; i++){
		this->location[i]->agent = NULL;
		//this->location[i] = NULL;
	}

	this->countLocations = 0;
}
/****************************************************************************/


void Agent::detach( VoronoiCell * correspondingVoronoiCell)
{
	int i;
	
	for( i=0; i<this->countLocations; i++){
		if( this->location[i] == correspondingVoronoiCell){
			this->location[i] = this->location[ --(this->countLocations)];
			correspondingVoronoiCell->agent = NULL;
			return;
		}
	}

	fprintf( stderr, "ERROR: in function void Agent::detach( VoronoiCell * correspondingVoronoiCell)\n");
	exit( 0);
	//this->countDirectFreeNeighbors    = correspondingVoronoiCell->countNeighborCells;	
	//this->countReachableFreeNeighbors = 0;
	
	
}

/****************************************************************************/


Agent::Agent( VoronoiCell * correspondingVoronoiCell)
{
	//Agent *newAgent = new Agent();
	this->attach( correspondingVoronoiCell);

	//return newAgent;
}
/****************************************************************************/


Agent::Agent()
{
	//Agent *newAgent = (Agent *) malloc( sizeof( Agent));

	this->location = NULL;
	
	this->actions = NULL;
	this->actionsInitialized = FALSE;
	
	this->initAgent();
	//newAgent->actions = NULL;
	this->divide = 1.;
	this->waste = 0;
	intoxication=0;

//	return newAgent;
}
/****************************************************************************/


void Agent::initAgent()
{
//	newAgent->location = correspondingVoronoiCell;
	this->countLocations = 0;
	
	// cell state
	this->state = FREE;				// is there a cell 0 - FREE, 1 - ACTIVE, 2 - NONACTIVE, 3 - NECROTIC, 4 - VESSEL

	// multiscale
	this->maxCellCount   = CountCellsPerVoronoiCell;
	this->cellCount      = 0;
	this->growingTumorCellCount = 0;
	this->dividingTumorCellCount = 0;
	this->necroticCellCount = 0;

	this->countFree = 0;
	this->countActive = 0;
	this->countNonactive = 0;


	if( this->actionsInitialized)
	for(int i=0; i<=M_gro; i++)
		this->actions[INDEX_GROWTH]->internalStateM[i]=0;

	this->waste = 0;

	// cell kinetics
	/*newAgent->mx;
	newAgent->my;
	newAgent->mz;*/

	//newAgent->pTissueData;

	// neighborships
	//newAgent->countDirectNeighbors        = correspondingVoronoiCell->countNeighborCells;	
	//newAgent->countDirectFreeNeighbors    = correspondingVoronoiCell->countNeighborCells;	
	//newAgent->countReachableNeighbors     = 0;	
	//newAgent->countReachableFreeNeighbors = 0;
	//newAgent->extendedNeighborhood          = NULL;	

	// cell dynamics
}
/****************************************************************************/


Agent ** initAgents( VoronoiDiagram* voronoiDiagram)
{
	int i;
	Agent ** agentArray = (Agent **)calloc( sizeof( Agent*), voronoiDiagram->countVoronoiCells);

	for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
		agentArray[i] = new Agent( );
		//agentArray[i] = newAgent( voronoiDiagram->voronoiCells[i]);
		//voronoiDiagram->voronoiCells[i]->agent = agentArray[i];
		//fprintf( stderr, "%i: ");
	}

	return agentArray;
}

/*Agent * GetAgent( VoronoiCell * voronoiCell)
{
	return ( Agent *) voronoiCell->agent;
}


void AddAgent( VoronoiCell * voronoiCell)
{
	voronoiCell->agent = (Agent *) calloc( sizeof(Agent));
}


void RemoveAgent( VoronoiCell * voronoiCell)
{
	voronoiCell->agent = NULL;
}


void DestroyAgent( VoronoiCell * voronoiCell)
{
	free( voronoiCell->agent);
}*/

/****************************************************************************/


const char *cellTypeToString( int type){
	switch( type){
		case 0:
		return "FREE";
		
		case 1:
		return "ACTIVE";
		
		case 2:
		return "NONACTIVE";
		
		case 3:
		return "NECROTIC";
		
		case 4:
		return "VESSEL";
		
		case COMPARTMENT:
		return "COMPARTMENT";
		
	}
	return "";
}

int Agent::isFree()
{
	return (this->state == FREE || (this->state == COMPARTMENT && this->cellCount<this->maxCellCount));
}

int Agent::countFreeNeighbors()
{
	int countFreeNeighbors = 0;

	for( int l=0; l<this->countLocations; l++)
		for( int n=0; n<this->location[l]->countNeighborCells; n++)
			if( this->location[l]->neighborCells[n]->isFree())
				countFreeNeighbors++;

	return countFreeNeighbors;
}

int Agent::isGrowing()
{							
	//fprintf( stderr, "cell %i isGrowing( %i/%i)? %s\n", this->index, this->cellCount, this->maxCellCount, (this->state == COMPARTMENT && this->growingTumorCellCount>0 ? "YES" : "NO"));
	return	( 	
				this->state == ACTIVE || 
				(
	       			this->state == COMPARTMENT && 
	        		this->growingTumorCellCount>0 &&
	        		(
	        			this->cellCount < this->maxCellCount ||
	        			GetVoronoiCell(this)->countFreeNeighborCells > 0 ||
	        			GetVoronoiCell(this)->countFreeExtendedNeighborCells > 0
	        		)
	        	)
	        ); // is there still an "ungrown" tumor cell?
}

int Agent::isDividing()
{
	//fprintf( stderr, "cell %i isDividing( %i/%i)? %s\n", this->index, this->cellCount, this->maxCellCount, (this->state == COMPARTMENT && (this->cellCount < this->maxCellCount || GetVoronoiCell(this)->countFreeNeighborCells>0 || GetVoronoiCell(this)->countFreeExtendedNeighborCells>0) &&  this->dividingTumorCellCount>0 ? "YES" : "NO"));
	return (this->state == ACTIVE || 
	       (this->state == COMPARTMENT && 
	        //this->tumorVolume > 0 &&                    // tumor cell to divide exists?
	        /*(this->cellCount < this->maxCellCount || GetVoronoiCell(this)->countFreeNeighborCells>0 || GetVoronoiCell(this)->countFreeExtendedNeighborCells>0) &&*/ // is there still empty space?
	        this->dividingTumorCellCount>0)); // is there still an "undivided" tumor cell?
//	        this->tumorVolume > this->tumorCellCount)); // is there still an "undivided" tumor cell?
}

int Agent::isDying()
{
	return (this->state == ACTIVE || this->state == NONACTIVE ||
	       (this->state == COMPARTMENT && 
	        this->growingTumorCellCount+this->dividingTumorCellCount > 0));                    // tumor cell to die exists?
}

int Agent::isLysing()
{
	return (this->state == NECROTIC ||
	       (this->state == COMPARTMENT && 
	        this->necroticCellCount > 0));                    // tumor cell to die exists?
}


