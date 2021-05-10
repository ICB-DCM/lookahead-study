#include "CellProcesses.h"

#include <stdio.h>
#include <math.h>

#include "Agent.h"
#include "Substrate.h"
#include "Mathematix.h"
#include "Interpolation.h"

#if USE_VESSEL_GROWTH
#include "VesselNetwork.h"
#endif



// GLOBAL VARIABLES

double MaxCellDivisionRate = (0.023/0.717734625);
double CellDivisionDepth   = 1;
int    M_div = 0;
int    M_gro = 1;
int    M_nec = 0;//10;

double CellMigrationRate 	= 0.;
double CellApoptosisRate	= 0.; //0.01;
double CellNecrosisRate 	= 0.0072;
double CellLysisRate 		= 0.35; // h^-1 => 3h
double VesselGrowthRate 	= 0.;


double alpha_ref = 0.1;
char ReentranceProbabilityFunction = HEAVYSIDE;
double ReentranceProbabilityReferenceLength = 130;

extern int Case;

// FUNCTION DEFINITION

double Action::getActualRate()
{
		float newRate = 0.;
	
		switch( this->type){


#if USE_VESSEL_GROWTH

			case VESSEL_GROWTH:

			if( this->originalCell->growthfactors > GROWTHFACTOR_THRESHOLD){
				newRate = VesselGrowthRate;
			}else{
				newRate = 0.;
			}
			break;
#endif


#if USE_DIVISION

			case DIVISION:		
			newRate = GetDivisionRate( this->originalCell) * (M_div + M_gro + 2.);
			

			break;

#endif //USE_DIVISION


#if USE_GROWTH

			case GROWTH:		
			newRate = GetGrowthRate( this->originalCell) * (M_div + M_gro + 2.);

			break;

#endif //USE_DIVISION


#if USE_MIGRATION

			case MIGRATION:
			
	#if TYPE_OF_MIGRATION == FREEMIGRATION
			//if( this->originalCell->getGlucose() * this->originalCell->getOxygen() < THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)
			newRate = CellMigrationRate;// / GetVoronoiCell(this->originalCell)->countFreeNeighborCells;
			// NOTOX: GRAVITY
			if( Agent::USE_GRAVITY){
			double sumAll=0;
			double sumFree=0;
			for( int n=0; n<this->originalCell->location[0]->countNeighborCells; n++){
				double dx = this->originalCell->location[0]->neighborCells[n]->position[0] - this->originalCell->location[0]->position[0];
				double dy = this->originalCell->location[0]->neighborCells[n]->position[1] - this->originalCell->location[0]->position[1];
#if DIMENSIONS == 3
				double dz = this->originalCell->location[0]->neighborCells[n]->position[2] - this->originalCell->location[0]->position[2];
#else
				double dz=0;
#endif
				double vx=0, vy=-1, vz=0;
				double gravitation = exp( - acos(vy*dy / sqrt(vx*vx+vy*vy+vz*vz) / sqrt(dx*dx+dy*dy+dz*dz)) / alpha_ref);

				sumAll += gravitation;
				if( this->originalCell->location[0]->neighborCells[n]->isFree())
					sumFree += gravitation;
			}
			newRate *= sumFree/sumAll;
			}
			// END NOTOX: GRAVITY
			if( this->originalCell->state == COMPARTMENT)
				newRate *= (double)this->originalCell->growingTumorCellCount / pow(this->originalCell->maxCellCount, 2./3.);

	#else
			newRate = 0.;
	#endif
	

			break;

#endif //USE_MIGRATION


#if USE_NECROSIS

			case NECROSIS:
			newRate = GetDeathRate( this->originalCell) * (M_nec + 1.);
			//printf("INFO: Actualize NECROSIS rate: p = %lf\n", newRate);
			break;

#endif


#if USE_LYSIS

			case LYSIS:
			newRate = CellLysisRate * (double)this->originalCell->necroticCellCount;
			//fprintf( stderr, "Lysis prob = %lf (= %lf * %i)\n", newRate, CellLysisRate, this->originalCell->necroticCellCount);
			break;

#endif


			case GAP:
			newRate = this->rate;
			break;

			default:
			fprintf( stderr, "ERROR: Not Identified Action (NIA) discovert!!! :) C'est grave ca!\n");
			exit( 0);
		}

		return newRate;	
}

int Action::active()
{
	return (this->top!=NULL ? TRUE : FALSE);
}

double Action::actualizeRate()
{
	this->rate = this->getActualRate();

	return this->rate;
}



Action* ActionTree::selectAction( double *time)
{
	//int i;
	double random, sum;

	// CHOOSE ACTION
	Action* p_elem = this->root;
  
	// randomly choose part of rate sum
	double temp_rand = myRand();
	random = temp_rand * this->rateSum;



	// sum probabilities until randomly choosen part is reached
	sum = 0.;
	do{
		//fprintf( stderr, "sum = %lf, random = %lf -> p_elem->rateSumPrev = %lf, p_elem->rateSumNext = %lf, p_elem->rate = %lf => %lf\n", sum, random, p_elem->rateSumPrev, p_elem->rateSumNext, p_elem->rate, p_elem->rateSumPrev+p_elem->rateSumNext+p_elem->rate);
		if( p_elem->sizePrev && random <= p_elem->rateSumPrev + sum){
			// take left branch

			p_elem = p_elem->prev;
		}else{
			sum += p_elem->rateSumPrev;
			if( p_elem->sizeNext && random <= p_elem->rateSumNext + sum){
				// take right branch

				p_elem = p_elem->next;
			}else{
				sum += p_elem->rateSumNext;
				if( random <= p_elem->rate + sum){
					// take this action
					//fprintf( stderr, "-> take this action\n");
					// calculate passed time
					random = myRand();
					*time += -log(1 - random) / this->rateSum;

				

					return p_elem;
					
				}else{
					fprintf( stderr, "ERROR: Something is going wrong in selectAction()\n");
					fprintf(stderr, "INFO: actionList->sum_of_prob=%lf, random=%lf (%lf), sum=%lf\n", this->rateSum, random, temp_rand, sum);

					exit( 0);	
				}
			}
		}
	}while( p_elem!=NULL);
	printf("Something is going wrong in selectAction()\n");
	printf("INFO: actionList->sum_of_prob=%lf, random=%lf (%lf), sum=%lf\n", this->rateSum, random, temp_rand, sum);
	exit(0);

}

void    ActionTree::actualizeAllRates( VoronoiDiagram* voronoiDiagram)
{
	Action* action = this->root;
	if( this->size > 0){

		Action *actionStack[ this->size], *actualAction;
		int stackSize = 0, lastStackSize = 0;
		int i;
		
		for( i=0; i<this->size; i++){
			actionStack[ i] = NULL;
		}
		
		actionStack[stackSize++] = action;
		do{
			actualAction = actionStack[stackSize-1];
			// TEST
			if( actualAction->next == actualAction || actualAction->next == actualAction){
				fprintf( stderr, "ERROR: in ActionTree::actualizeAllRates()\nactualAction->next == actualAction || actualAction->next == actualAction\n");
				exit( 0);
			}
			// root condition
			if( actualAction->top == actualAction && actualAction!=this->root){
				fprintf( stderr, "ERROR: in ActionTree::actualizeAllRates()\nactualAction->top == actualAction && actualAction!=this->root\n");
				exit( 0);
			}
			// TEST
			
			if( lastStackSize < stackSize){
				lastStackSize = stackSize;
				if( actualAction->sizePrev > 0){
					actionStack[stackSize++] = actualAction->prev;
				}

			}else{
				lastStackSize = stackSize;
				if( actualAction->sizeNext > 0 &&  actualAction->next != actionStack[stackSize]){
					//rateDifferenceStack[stackSize] = 0.;
					actionStack[stackSize++] = actualAction->next;
				}
				// decend tree
				else{
					// adapt own rate
					actualAction->actualizeRate();
					
					stackSize--;

					// adapt parent node
					if( stackSize>0){
						if( actionStack[stackSize-1]->prev == actualAction)
							actionStack[stackSize-1]->rateSumPrev = actionStack[stackSize]->rateSumPrev + actionStack[stackSize]->rateSumNext + actionStack[stackSize]->rate;
						else
							actionStack[stackSize-1]->rateSumNext = actionStack[stackSize]->rateSumPrev + actionStack[stackSize]->rateSumNext + actionStack[stackSize]->rate;

					}else{
						this->rateSum = actionStack[stackSize]->rateSumPrev + actionStack[stackSize]->rateSumNext + actionStack[stackSize]->rate;
					}
					
				}
			}
		}while( stackSize>0);
	}
}

void ActionTree::destroyActionTree()
{
	if( this->size > 0){

		Action **actionStack = (Action**) malloc( this->size*sizeof(Action*)), *actualAction;
		int stackSize = 0, lastStackSize = 0;
		int i;
		
		for( i=0; i<this->size; i++){
			actionStack[ i] = NULL;
		}
		
		actionStack[stackSize++] = this->root;
		do{
			actualAction = actionStack[stackSize-1];
		
			if( lastStackSize < stackSize){
				// climbing tree -> prev
				lastStackSize = stackSize;
				if( actualAction->sizePrev > 0){
					actionStack[stackSize++] = actualAction->prev;
				}
			}else{
				lastStackSize = stackSize;
				// climb tree -> next
				if( actualAction->sizeNext > 0 &&  actualAction->next != actionStack[stackSize]){
					actionStack[stackSize++] = actualAction->next;
				}
				// decend tree
				else{
					// reinit action
					actualAction->top = NULL;
					actualAction->prev = NULL;
					actualAction->next = NULL;
					actualAction->sizePrev = 0;
					actualAction->sizeNext = 0;
					actualAction->rateSumPrev = 0.;
					actualAction->rateSumNext = 0.;
					
					stackSize--;
					//fprintf( stderr, "DECENDING: stackSize = %i, from node %i\n", stackSize, actualAction->type);
				}
			}
		}while( stackSize>0);
		free( actionStack);
	}
	free( this);
}

void Action::print()
{
	fprintf( stderr, "[%s: rate=%lf] -> prevSize=%i | nextSize=%i\n", actionTypeToString( this->type), this->rate, this->sizePrev, this->sizeNext);
	if( this->sizePrev>0){
		fprintf( stderr, "prev [size=%i, rateSum=%lf]:\n", this->sizePrev, this->rateSumPrev);
		this->prev->print();
	}
	if( this->sizeNext>0){
		fprintf( stderr, "next [size=%i, rateSum=%lf]:\n", this->sizeNext, this->rateSumNext);
		this->next->print();
	}
}

void ActionTree::print()
{
	fprintf( stderr, "root [size=%i, rateSum=%lf]:\n", this->size, this->rateSum);
	if( this->size>0){
		this->root->print();
	}
}
/*****************************************************************************/



ActionTree *ActionTree::newActionTree()
{
	// memory allocation
	ActionTree* actionTree = ( ActionTree*) malloc( sizeof(ActionTree));

	// attribute initialization
	actionTree->size	= 0;
	//actionTree->depth	= 0;
	actionTree->root	= NULL;

	return actionTree;
}
/*****************************************************************************/



ActionList* newActionList(){

  ActionList* actionList = ( ActionList*) malloc( sizeof(ActionList));

  actionList->length      = 0;
  actionList->sum_of_prob = 0.;
  actionList->head        = NULL;

  return actionList;
}
/*****************************************************************************/



Action* newAction( int type, double rate){

	Action* actions = ( Action*) malloc( sizeof( Action));

    if (actions == NULL) {
        printf("Got NULL pointer in newAction.\n");
    }
	// initialize action attributes	   
	actions->type	= type;
	actions->rate	= rate;
	actions->internalState	= 0;
	actions->originalCell	= NULL;
	actions->destinationCell	= NULL;

	// action list
	actions->next	= NULL;
	actions->prev	= NULL;


#if USE_ACTION_TREE
	actions->top	= NULL;    
	
	actions->rateSumNext	= 0.;
	actions->sizeNext	= 0;

	actions->rateSumPrev	= 0.;
	actions->sizePrev	= 0;
#endif
  
	return actions;  
}
/*****************************************************************************/



void destroyActionList( ActionList* actionList){


        free( actionList);
}
/*****************************************************************************/




void initCellActions( Agent* cell){

	int countInitializedActions = 0;

	if( cell->actionsInitialized == FALSE){
		cell->actions = (Action **) calloc( sizeof( Action*), USE_DIVISION + USE_NECROSIS + USE_LYSIS + USE_GROWTH + USE_MIGRATION

		);

#if USE_DIVISION 		
		// division
		cell->actions[INDEX_DIVISION] = newAction( DIVISION, 0.);
		cell->actions[INDEX_DIVISION]->originalCell = cell;
		countInitializedActions++;
#endif // USE_DIVISION 		


#if USE_MIGRATION 		
		// diffusion
		cell->actions[INDEX_MIGRATION] = newAction( MIGRATION, 0.);
		cell->actions[INDEX_MIGRATION]->originalCell    = cell;
		countInitializedActions++;

#endif // USE_MIGRATION 		

#if USE_LYSIS 		
		// lysis
		cell->actions[INDEX_LYSIS] = newAction( LYSIS, 0.);
		cell->actions[INDEX_LYSIS]->originalCell = cell;
		countInitializedActions++;
#endif // USE_LYSIS 		
		
#if USE_NECROSIS 		
		// necrosis
		cell->actions[INDEX_NECROSIS] = newAction( NECROSIS, 0.);//new_prob_necrosis_element( cell);
		cell->actions[INDEX_NECROSIS]->originalCell = cell;
		countInitializedActions++;
#endif // USE_NECROSIS 		

#if USE_GROWTH 		
		// growth 
		cell->actions[INDEX_GROWTH] = newAction( GROWTH, 0.);//new_prob_necrosis_element( cell);
		cell->actions[INDEX_GROWTH]->originalCell = cell;
		countInitializedActions++;
		for(int i=0; i<=M_gro; i++)
			cell->actions[INDEX_GROWTH]->internalStateM[i]=0;

#endif // USE_GROWTH 		

		cell->actionsInitialized = TRUE;
		
	}
}

/*****************************************************************************/



#if USE_ACTION_TREE


void ActionTree::addAction( Action* action)
{

	if( action->top != NULL){
		return;
	}
		
	/********        add action to prob_list    **********/

	// init action
	action->actualizeRate();


	action->next = NULL;     // next action in probability list
	action->rateSumNext = 0.;
	action->sizeNext = 0;

	action->prev = NULL;     // previous action in probability list
	action->rateSumPrev = 0.;
	action->sizePrev = 0;

	
	if( this->size == 0){
		this->root = action;
		this->root->top = this->root;
		this->rateSum = action->rate;
		this->size++;
		return;
	}else{
		Action *actualAction = this->root;
		do{
			if( actualAction->sizePrev <= actualAction->sizeNext){
				actualAction->sizePrev++;

				/* FAST */ actualAction->rateSumPrev += action->rate;
				if( actualAction->sizePrev > 1)				
					actualAction = actualAction->prev;
				else{
					actualAction->prev = action;
					action->top = actualAction;
				}
				
			}else{
				actualAction->sizeNext++;
				/* FAST */ actualAction->rateSumNext += action->rate;

				if( actualAction->sizeNext > 1)				
					actualAction = actualAction->next;
				else{
					actualAction->next = action;
					action->top = actualAction;
				}
			}
		}while( action->top == NULL);
		/* FAST */ this->rateSum += action->rate;
		this->size++;
	}

}


int ActionTree::getDepth( Action* action)
{
	fprintf( stderr, "ActionTree::getDepth()\n");
				
	if( action == NULL){
		return 0;
	}//if( this->size == 1)
	//	return 1;
	else{
		Action *actionStack[ this->size], *actualAction;
		int stackSize = 0, lastStackSize = 0;
		int maxDepth = 0;
		//stackSize++;
		int i;
		for( i=0; i<this->size; i++)
			actionStack[ i] = NULL;
			
		actionStack[stackSize++] = action;
		do{
			actualAction = actionStack[stackSize-1];
			// TEST
			if( (actualAction->sizePrev == 0 && actualAction->rateSumPrev!=0.) || (actualAction->sizeNext == 0 && actualAction->rateSumNext!=0.)){
				fprintf( stderr, "ERROR: in ActionTree::getDepth()\nactualAction->sizePrev == 0 && actualAction->rateSumPrev!=0. || actualAction->sizeNext == 0 && actualAction->rateSumNext!=0.\n");
				exit( 0);
			}
			if( actualAction->next == actualAction || actualAction->next == actualAction){
				fprintf( stderr, "ERROR: in ActionTree::getDepth()\nactualAction->next == actualAction || actualAction->next == actualAction\n");
				exit( 0);
			}
			// root condition
			if( actualAction->top == actualAction && actualAction!=this->root){
				fprintf( stderr, "ERROR: in ActionTree::getDepth()\nactualAction->top == actualAction && actualAction!=this->root\n");
				exit( 0);
			}
			// rateSums 
			double temp_rateSum = 0.;
			if( actualAction->sizePrev > 0)
				temp_rateSum = actualAction->prev->rateSumPrev + actualAction->prev->rateSumNext + actualAction->prev->rate;
			if( actualAction->rateSumPrev!=temp_rateSum) {
				fprintf( stderr, "ERROR: in ActionTree::getDepth()\nactualAction->rateSumPrev(%.10lf)!=temp_rateSum(%.10lf) =>error: %e\n", actualAction->rateSumPrev, temp_rateSum, actualAction->rateSumPrev - temp_rateSum);
				exit( 0);
			}
				
			temp_rateSum = 0.;
			if( actualAction->sizeNext > 0)
				temp_rateSum = actualAction->next->rateSumPrev + actualAction->next->rateSumNext + actualAction->next->rate;
			if( actualAction->rateSumNext!=temp_rateSum) {
				fprintf( stderr, "ERROR: in ActionTree::getDepth()\nactualAction->rateSumNext(%.10lf)!=temp_rateSum(%.10lf) =>error: %e\n", actualAction->rateSumNext, temp_rateSum, actualAction->rateSumNext - temp_rateSum);
				exit( 0);
			}
				
			// TEST
			
			if( lastStackSize < stackSize){

				lastStackSize = stackSize;

				if( actualAction->sizePrev > 0){

					actionStack[stackSize++] = actualAction->prev;

				}
				
			}else{
				lastStackSize = stackSize;

				if( actualAction->sizeNext > 0 &&  actualAction->next != actionStack[stackSize]){
					actionStack[stackSize++] = actualAction->next;

				}
				// decend tree
				else{
					stackSize--;

				}
			}
			if( stackSize > maxDepth)
				maxDepth = stackSize;
		}while( stackSize>0);
		

		return maxDepth;
	}
}

void ActionTree::actualizeRate( Action* action, double newRate)
{

	/* FAST */ double rateDifference = newRate - action->rate;
	
	if(action->top == NULL)
	return;


	// actualize action
	action->rate = newRate;
	
	// actualize tree
	Action *actualAction = action,
	       *lastAction;
	
	
	while( actualAction != this->root){
		lastAction = actualAction;
		actualAction = actualAction->top;
		if( actualAction->prev == lastAction){
			/* FAST */ actualAction->rateSumPrev += rateDifference;
		}else{ 
			if(actualAction->next == lastAction)
				/* FAST */ actualAction->rateSumNext += rateDifference;
			else{
				fprintf(stderr, "ERROR in ActionTree::actualizeRate()\n");
				exit( 0);
			}
		}
	}

	/* FAST */ this->rateSum += rateDifference;


}

void ActionTree::deleteAction( Action* action)
{

	/********        add action to prob_list    **********/
	if( action->top != NULL){
		// tree until root
		Action *actualAction = action;
		while( actualAction != this->root){
			if( actualAction == actualAction->top->prev){
				// actualize prev attributes
				actualAction->top->sizePrev--;
				/* FAST */ actualAction->top->rateSumPrev -= action->rate;
			}
			else{
				// actualize next attributes
				actualAction->top->sizeNext--;
				/* FAST */ actualAction->top->rateSumNext -= action->rate;
			}
			actualAction = actualAction->top;
		}
		
		// rearange subtrees
		Action* shortSubtree, *longSubtree;

		if( action->sizePrev <= action->sizeNext){
			shortSubtree = action->prev;
			longSubtree  = action->next;

		}else{
			shortSubtree = action->next;
			longSubtree  = action->prev;

		}

		/* FAST */ this->rateSum -= action->rate;
		if( --(this->size) == 0){
			// tree is empty
			this->root = NULL;

		}else{

		
			if( shortSubtree == NULL){
				// only one subtree has to be added
				if( this->root == action)
					this->root = longSubtree;
				else{	
					if( action == action->top->prev){
						// actualize prev attributes
						action->top->prev = longSubtree;
					}else{
						// actualize next attributes
						action->top->next = longSubtree;
					}
				}
				if( longSubtree != NULL)		
					longSubtree->top = action->top;
					
				
			}else{
				// attach first subtree
				if( this->root == action)
					this->root = shortSubtree;
				else{	
					if( action == action->top->prev){
						// actualize prev attributes
						action->top->prev = shortSubtree;
					}else{
						// actualize next attributes
						action->top->next = shortSubtree;
					}
				}		
				shortSubtree->top = action->top;
				
				// attach seconde subtree
				Action* lastAction = actualAction = shortSubtree;
				
				while( actualAction!=NULL){
					lastAction = actualAction;
					if( actualAction->sizePrev <= actualAction->sizeNext){
						actualAction->sizePrev    += longSubtree->sizePrev    + longSubtree->sizeNext    + 1;
						actualAction->rateSumPrev += longSubtree->rateSumPrev + longSubtree->rateSumNext + longSubtree->rate;
						actualAction = actualAction->prev;
						if( actualAction==NULL)
							lastAction->prev = longSubtree;
					}else{
						actualAction->sizeNext    += longSubtree->sizePrev    + longSubtree->sizeNext    + 1;
						actualAction->rateSumNext += longSubtree->rateSumPrev + longSubtree->rateSumNext + longSubtree->rate;
						actualAction = actualAction->next;
						if( actualAction==NULL)
							lastAction->next = longSubtree;
					}					
				}

				longSubtree->top = lastAction;

			}
				
		}
		
		// detach from tree
		action->top = NULL;
		action->prev = NULL;
		action->next = NULL;
		action->sizePrev = 0;
		action->sizeNext = 0;
		action->rateSumPrev = 0.;
		action->rateSumNext = 0.;
	}

}
#endif
void addNewAction( ActionList *actionList, int type, double rate, Agent *cell){
  
	Action* temp_prev;
	Action* temp_next;
	Action* new_elem;

	// add actions to prob_list
	//printf("addAction(): type:%i, prob:%lf", type, probability_of_actions);
	new_elem = (Action*) newAction( type, rate);
	new_elem->originalCell = cell;

	//new_elem = cell->actions[0];

	if( actionList->length == 0){
		actionList->head = new_elem;
		actionList->head->next = actionList->head;
		actionList->head->prev = actionList->head;
		actionList->length++;  

		return;
	}

	temp_next = actionList->head;
	temp_prev = actionList->head->prev;
      
	/* include new element as head of prob_list */
	new_elem->prev = temp_prev;
	temp_prev->next = new_elem;
	new_elem->next = temp_next;
	temp_next->prev = new_elem;
	/* add division actions to cell actions list */
	actionList->head = new_elem;
      
	/* add probability to sum of probabilities */
	actionList->sum_of_prob += new_elem->rate;
	/* increase actionList_length */
	actionList->length++;  
}
/*****************************************************************************/




void addAction( ActionList *probList, Action* action){

	

	Action* temp_prev;
	Action* temp_next;

	/********        add action to prob_list    **********/
	if( probList->length == 0){
		probList->head = action;
		probList->head->next = probList->head;
		probList->head->prev = probList->head;
		probList->length++;  
		return;
	}

	temp_next = probList->head;
	temp_prev = probList->head->prev;
      
	/* include new element as head of prob_list */
	action->prev = temp_prev;
	temp_prev->next = action;
	action->next = temp_next;
	temp_next->prev = action;

	/* add division action to cell action list */
	probList->head = action;
      
	/* add probability to sum of probabilities */
	probList->sum_of_prob += action->rate;
	probList->length++;  
}
/*****************************************************************************/




/*void addActions( ActionList *actionList, Agent* originalCell, Agent* destinationCell){

	int i;
	Action* temp_prev;
	Action* temp_next;
	Action* new_elem;

	initCellActions( originalCell);

//	for(i = 0; i < originalCell->locationcountDirectNeighbors + 1; i++){
		for(i = 0; i < originalCell->countNeighborCells + 1; i++){
			if( originalCell->actions[i]->destinationCell == destinationCell){
				// ********        add action to prob_list    **********
				new_elem = originalCell->actions[i];
				temp_next = actionList->head;
				temp_prev = actionList->head->prev;
      
				// include new element as head of prob_list 
				new_elem->prev = temp_prev;
				temp_prev->next = new_elem;
				new_elem->next = temp_next;
				temp_next->prev = new_elem;
				// add division action to cell action list 
				actionList->head = new_elem;
      
				// add probability to sum of probabilities 
				actionList->sum_of_prob += new_elem->rate;
				// increase actionList_length 
				//printf("actionList->length++\n");
				actionList->length++;  
			}
		}
//	}
}*/
/*****************************************************************************/




void deleteAction( ActionList *probList, Action* cellAction){
	//printf("DELETE ACTIONS - begin (cell->actions_initialized = %i)\n", cellAction->type);
  
	Action* tempPrev;
	Action* tempNext;

	// changing rateList    
	probList->sum_of_prob -= cellAction->rate;
	probList->length--;

	tempPrev = cellAction->prev;
	tempNext = cellAction->next;
	tempPrev->next = tempNext;
	tempNext->prev = tempPrev;
	cellAction->prev = NULL;
	cellAction->next = NULL;
	//cellAction->actions_active = FALSE;

	if(cellAction == probList->head){
		probList->head = tempNext;
	}
}
/*****************************************************************************/




void ActionTree::deleteAllActions( Agent* cell){

	//if( cell != NULL && cell->actionsInitialized && cell->actions!= NULL)
	if( cell->actionsInitialized)
	for( int a=0; a<INDEX_MIGRATION; a++){
		//fprintf( stdout, "Delete %i\n", a);
		//if( cell->actions[a]!= NULL)
		this->deleteAction( cell->actions[a]);
	}
}
/*****************************************************************************/



void addedVesselAndUpdateSurrounding( ActionTree *actionList, AgentList* agentArray, VoronoiCell* addedVessel, VoronoiDiagram* voronoiDiagram ){
	int i;
	int oldState = addedVessel->getState();

	fprintf( stdout, "INFO: addedVesselAndUpdateSurrounding(...) -> cell %i\n", addedVessel->index);
	//exit( 0);


	// ******     UPDATE ADDED CELL      ******

	// attach agent if not already done
	if( GetAgent( addedVessel)==NULL){
		Agent* daughterCell = agentArray->activateAgent();
		daughterCell->attach( addedVessel);
		initCellActions( daughterCell);
	}

	// set status
	GetAgent( addedVessel)->state = VESSEL;
	GetAgent( addedVessel)->actualize( actionList);


	// ******      UPDATE NEIGHBORS      ******

	if( oldState == FREE){
	//printf("UPDATE NEIGHBORS\n");
	for (i = 0; i < addedVessel->countNeighborCells; i++){
		//printf("UPDATE NEIGHBOR: %i (%i/%i), state: %s\n", added_cell->neighborCells[i]->index, i, added_cell->countNeighborCells, cellTypeToString( added_cell->neighborCells[i]->getState()));

		addedVessel->neighborCells[i]->countFreeNeighborCells--;

		switch( addedVessel->neighborCells[i]->getState()){
			case ACTIVE:

			// reduce number of free neighborCells
			//addedVessel->neighborCells[i]->countFreeNeighborCells--;
			if( addedVessel->neighborCells[i]->countFreeNeighborCells == 0 && 
			    addedVessel->neighborCells[i]->countFreeExtendedNeighborCells == 0){

				GetAgent( addedVessel->neighborCells[i])->actualize( actionList);
			}
			break;

			default:
			break;
		}
	}


	/******     UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k      ******/

	//printf("UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k\n");
	for (i = 0; i < addedVessel->countExtendedNeighborCells; i++){
		//printf("UPDATE EXTENDED NEIGHBOR: %i (%i/%i), state: %s, countFreeExtendedNeighborCells: %i\n", added_cell->extendedNeighborhood[i]->index, i, added_cell->countExtendedNeighborCells, cellTypeToString( added_cell->extendedNeighborhood[i]->getState()), added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells);

		// reduce number of free neighborCells
		addedVessel->extendedNeighborhood[i]->countFreeExtendedNeighborCells--;

		switch( addedVessel->extendedNeighborhood[i]->getState()){
			case ACTIVE:
			if( addedVessel->extendedNeighborhood[i]->countFreeNeighborCells == 0 && 
			    addedVessel->extendedNeighborhood[i]->countFreeExtendedNeighborCells == 0){

				GetAgent( addedVessel->extendedNeighborhood[i])->actualize( actionList);
				
			}
			break;

			default:
			break;
		}
	}
	}else{
		fprintf( stderr, "WARNING: addedVesselAndUpdateSurrounding(...)\n");
		fprintf( stderr, "An cell occupied by a %s agent replaced by VESSEL\n", cellTypeToString(oldState));
	}


}
/*****************************************************************************/



//void update_surrounding_added_cell(ActionList *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram ){
void update_surrounding_added_cell(ActionTree *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram ){
	int i;



	/******      UPDATE AGENT          ******/

	//GetAgent( added_cell)->state = ACTIVE;
	initCellActions( GetAgent( added_cell));
	GetAgent( added_cell)->actualize( actionList);



	/******      UPDATE NEIGHBORS      ******/

	if( !added_cell->isFree()){
	//printf("UPDATE NEIGHBORS\n");
	for (i = 0; i < added_cell->countNeighborCells; i++){
		//printf("UPDATE NEIGHBOR: %i (%i/%i), state: %s\n", added_cell->neighborCells[i]->index, i, added_cell->countNeighborCells, cellTypeToString( added_cell->neighborCells[i]->getState()));

		added_cell->neighborCells[i]->countFreeNeighborCells--;

		switch( added_cell->neighborCells[i]->getState()){
			case ACTIVE:

			// reduce number of free neighborCells
			//added_cell->neighborCells[i]->countFreeNeighborCells--;
			if( added_cell->neighborCells[i]->countFreeNeighborCells == 0 && 
			    added_cell->neighborCells[i]->countFreeExtendedNeighborCells == 0){

				GetAgent( added_cell->neighborCells[i])->actualize( actionList);
			}
			break;

			default:
			break;
		}
	}


	/******     UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k      ******/

	//printf("UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k\n");
	for (i = 0; i < added_cell->countExtendedNeighborCells; i++){
		//printf("UPDATE EXTENDED NEIGHBOR: %i (%i/%i), state: %s, countFreeExtendedNeighborCells: %i\n", added_cell->extendedNeighborhood[i]->index, i, added_cell->countExtendedNeighborCells, cellTypeToString( added_cell->extendedNeighborhood[i]->getState()), added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells);

		// reduce number of free neighborCells
		added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells--;

		switch( added_cell->extendedNeighborhood[i]->getState()){
			case ACTIVE:
			if( added_cell->extendedNeighborhood[i]->countFreeNeighborCells == 0 && 
			    added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells == 0){

				GetAgent( added_cell->extendedNeighborhood[i])->actualize( actionList);
				
			}
			break;

			default:
			break;
		}
	}
	}else{
		fprintf( stderr, "WARNING: update_surrounding_added_cell(...)\n");
		fprintf( stderr, "The cell %i is declared as %s!\n", added_cell->index, cellTypeToString(added_cell->getState()));
	}

}
/*****************************************************************************/



void update_expanded_cell(ActionList *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram ){



	/******     UPDATE ADDED CELL      ******/



	if( added_cell->countFreeNeighborCells == 0 && added_cell->countFreeExtendedNeighborCells == 0)
	{
		// set status
		GetAgent( added_cell)->state = NONACTIVE;
		//exit(0);
	}

}


void update_surrounding_expanded_cell(ActionTree *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram ){
	int i;



	/******      UPDATE NEIGHBORS      ******/

	//printf("update_surrounding_expanded_cell( %i)\n", added_cell->index);
	//printf("UPDATE NEIGHBORS\n");
	if( !added_cell->isFree()){
	//printf("UPDATE NEIGHBORS\n");
	for (i = 0; i < added_cell->countNeighborCells; i++){
		//printf("UPDATE NEIGHBOR: %i->%i (%i/%i), state: %s\n", added_cell->index, added_cell->neighborCells[i]->index, i, added_cell->countNeighborCells, cellTypeToString( added_cell->neighborCells[i]->getState()));

			// reduce number of free neighborCells
			added_cell->neighborCells[i]->countFreeNeighborCells--;

		switch( added_cell->neighborCells[i]->getState()){
			case COMPARTMENT:
			case ACTIVE:

			if( added_cell->neighborCells[i]->countFreeNeighborCells == 0 && 
			    added_cell->neighborCells[i]->countFreeExtendedNeighborCells == 0){

				GetAgent( added_cell->neighborCells[i])->actualize( actionList);

			}
			break;

			default:
			break;
		}
	}


	/******     UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k      ******/

	//printf("UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k\n");
	for (i = 0; i < added_cell->countExtendedNeighborCells; i++){
		//printf("UPDATE EXTENDED NEIGHBOR: %i (%i/%i), state: %s, countFreeExtendedNeighborCells: %i\n", added_cell->extendedNeighborhood[i]->index, i, added_cell->countExtendedNeighborCells, cellTypeToString( added_cell->extendedNeighborhood[i]->getState()), added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells);

		if( added_cell->extendedNeighborhood[i]->countExtendedNeighborCells != 0){
			added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells--;
		
			switch( added_cell->extendedNeighborhood[i]->getState()){
				case COMPARTMENT:
				case ACTIVE:
				// reduce number of free neighborCells
				if( added_cell->extendedNeighborhood[i]->countFreeNeighborCells == 0 && 
				    added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells == 0){
	
					GetAgent( added_cell->extendedNeighborhood[i])->actualize( actionList);
					
				}
				break;

				default:
				break;
			}
		}
	}
	}
	
}
/*****************************************************************************/



void update_surrounding_expanded_compartment(ActionTree *actionList, VoronoiCell* added_cell, VoronoiDiagram* voronoiDiagram ){
	int i;

	//printf("update_surrounding_expanded_cell( %i)\n", GetAgent(added_cell)->index);

	if( !added_cell->isFree()){

		/******      UPDATE NEIGHBORS      ******/

		//printf("UPDATE NEIGHBORS\n");
		for (i = 0; i < added_cell->countNeighborCells; i++){
			// reduce number of free neighborCells
			added_cell->neighborCells[i]->countFreeNeighborCells--;

			//if( !added_cell->neighborCells[i]->isFree()){
			if( GetAgent(added_cell->neighborCells[i])!=NULL){
				//printf("update %i\n", GetAgent(added_cell->neighborCells[i])->index);
				GetAgent( added_cell->neighborCells[i])->actualize( actionList);
			}
		}


		/******     UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k      ******/

		//printf("UPDATE NEIGHBORS WITHIN DIVIDING DEPTH k\n");
		for (i = 0; i < added_cell->countExtendedNeighborCells; i++){
			added_cell->extendedNeighborhood[i]->countFreeExtendedNeighborCells--;
		
			if( !added_cell->extendedNeighborhood[i]->isFree()){
				GetAgent( added_cell->extendedNeighborhood[i])->actualize( actionList);
			}
		}
	}
	
}
/*****************************************************************************/



/*void update_deleted_cell(ActionList *actionList, Agent* deleted_cell, VoronoiDiagram* voronoiDiagram ){


	// ******    update deleted cell    ****** //

	if( deleted_cell->state != NECROTIC){
		if( GetVoronoiCell( deleted_cell)->countFreeNeighborCells != 0 || GetVoronoiCell( deleted_cell)->countFreeExtendedNeighborCells != 0 ){
#if USE_DIVISION
			// delete division
			deleteAction( actionList, deleted_cell->actions[0]);
#endif // USE_DIVISION

#if USE_MIGRATION
			// delete diffusions
			int i;
			for (i = 0; i < deleted_cell->countNeighborCells; i++){
				if( deleted_cell->neighborCells[i]->state == FREE){
					deleteAction( actionList, deleted_cell->actions[INDEX_MIGRATION+i]);
				}
			}
#endif // USE_MIGRATION
		}

#if USE_LYSIS
		// delete lysis
		addAction( actionList, deleted_cell->actions[ GetVoronoiCell( deleted_cell)->countNeighborCells + 1]);
#endif //USE_LYSIS

#if USE_NECROSIS
		// delete necrosis
		deleteAction( actionList, deleted_cell->actions[ GetVoronoiCell( deleted_cell)->countNeighborCells + 2]);
#endif //USE_NECROSIS
	}


	// set status
	deleted_cell->state = NECROTIC;
}*/
/*****************************************************************************/




/*void update_surrounding_deleted_cell( ActionList *actionList, Agent* deleted_cell, VoronoiDiagram* voronoiDiagram ){
	int i;
#if USE_MIGRATION
	int j;
#endif


	// ******    update neighborCells    ****** /

	for (i = 0; i < GetVoronoiCell( deleted_cell)->countNeighborCells; i++){

		switch( GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i])->state){
			case ACTIVE:
			// add diffusion actions from neighbor to deleted cell
#if USE_MIGRATION
			for(j = 0; j < deleted_cell->neighborCells[i]->countNeighborCells; j++){
				if( deleted_cell->neighborCells[i]->actions[j]->destinationCell->index == deleted_cell->index){
					addActions( actionList, deleted_cell->neighborCells[i], deleted_cell->neighborCells[i]->actions[j]->destinationCell);
				}
			}
#endif
			// increase number of free neighborCells
			GetVoronoiCell( deleted_cell)->neighborCells[i]->countFreeNeighborCells++;
			break;

			case NONACTIVE:
			// add diffusion actions from neighbor to deleted cell
#if USE_MIGRATION
			for(j = 1; j <= GetAgent( deleted_cell->neighborCells[i])->countNeighborCells; j++){
				if deleted_cell->neighborCells[i])->actions[j]->destinationCell->nr == GetAgent( deleted_cell)->nr){
					addActions( actionList, GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i])), GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i])->actions[j]->destinationCell);
				}
			}
#endif
			// add division actions
#if USE_DIVISION
			//addActions( actionList, GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i]), GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i])->actions[0]->destinationCell);
			//addDivisionAction( actionList, GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i]), voronoiDiagram);
			addAction( actionList, GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i])->actions[0]);
#endif // USE_MIGRATION
			// increase number of free neighborCells
			GetVoronoiCell( deleted_cell)->neighborCells[i]->countFreeNeighborCells++;
			// change state
			GetAgent( GetVoronoiCell( deleted_cell)->neighborCells[i])->state = ACTIVE;
			break;
    
			default:
			// reduce number of free neighborCells
			GetVoronoiCell( deleted_cell)->neighborCells[i]->countFreeNeighborCells++;
		}
	}


  // ******    update neighborCells within kr   ****** //
  for (i = 0; i < deleted_cell->location[0]->countExtendedNeighborCells; i++){

    switch( GetExtendedNeighborAgent( deleted_cell, i)->state){
      case ACTIVE:
        // increase number of free kr neighborCells
        GetExtendedNeighborCell( deleted_cell, i)->countFreeExtendedNeighborCells++;
        break;

      case NONACTIVE:
        // add division actions
        //addActions( actionList, GetExtendedNeighborCell( deleted_cell, i), GetExtendedNeighborCell( deleted_cell, i)->actions[0]->destinationCell);
		//addDivisionAction( actionList, GetExtendedNeighborCell( deleted_cell, i), voronoiDiagram);
		addAction( actionList, GetExtendedNeighborAgent( deleted_cell, i)->actions[0]);
        // increase number of free kr neighborCells
        GetExtendedNeighborCell( deleted_cell, i)->countFreeExtendedNeighborCells++;
        // change state
        GetExtendedNeighborAgent( deleted_cell, i)->state = ACTIVE;
        break;
    
      default:
       // reduce number of free kr neighborCells
        GetExtendedNeighborCell( deleted_cell, i)->countFreeExtendedNeighborCells++;
    }
  }
}*/
/*****************************************************************************/




void update_surrounding_removed_cell( ActionTree *actionList, Agent* deleted_cell, VoronoiDiagram* voronoiDiagram ){
	int i, ii;


	/******    update neighborCells    ******/

	for(ii=0; ii<deleted_cell->countLocations; ii++)
	for (i = 0; i < deleted_cell->location[ii]->countNeighborCells; i++){

		VoronoiCell *currentNeighbor = deleted_cell->location[ii]->neighborCells[i];
		
		// increase number of free neighborCells
		currentNeighbor->countFreeNeighborCells++;
			
		switch( currentNeighbor->getState()){

			case NONACTIVE:
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}
	}


	/******    update neighborCells    ******/

	for(ii=0; ii<deleted_cell->countLocations; ii++)
	for (i = 0; i < deleted_cell->location[ii]->countExtendedNeighborCells; i++){

		VoronoiCell *currentNeighbor = deleted_cell->location[ii]->extendedNeighborhood[i];
		
		// increase number of free neighborCells
		currentNeighbor->countFreeExtendedNeighborCells++;
			
		switch( currentNeighbor->getState()){

			case NONACTIVE:
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}
	}
}
/*****************************************************************************/




void update_surrounding_reduced_compartment( ActionTree *actionList, Agent* deleted_cell, VoronoiDiagram* voronoiDiagram ){
	int i, ii = 0;


	/******    update neighborCells    ******/

	if( deleted_cell->cellCount + 1 == deleted_cell->maxCellCount){
	for (i = 0; i < deleted_cell->location[ii]->countNeighborCells; i++){

		VoronoiCell *currentNeighbor = deleted_cell->location[ii]->neighborCells[i];
		
		// increase number of free neighborCells
		currentNeighbor->countFreeNeighborCells++;

		if( currentNeighbor->agent != NULL)
		GetAgent( currentNeighbor)->actualize( actionList);
			
		/*switch( currentNeighbor->getState()){

			case NONACTIVE:
			//fprintf( stderr, "INFO: NONACTIVE\n");
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}*/
	}


	/******    update neighborCells    ******/

	//for(ii=0; ii<deleted_cell->countLocations; ii++)
	for (i = 0; i < deleted_cell->location[ii]->countExtendedNeighborCells; i++){

		VoronoiCell *currentNeighbor = deleted_cell->location[ii]->extendedNeighborhood[i];
		
		// increase number of free neighborCells
		//currentNeighbor->countFreeExtendedNeighborCells++;
			
		if( currentNeighbor->agent != NULL)
		GetAgent( currentNeighbor)->actualize( actionList);
			
		/*switch( currentNeighbor->getState()){

			case NONACTIVE:
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}*/
	}
	}
}
/*****************************************************************************/




void update_surrounding_reduced_cell( ActionTree *actionList, VoronoiCell* deleted_cell, VoronoiDiagram* voronoiDiagram ){
	int i;//, ii;


	/******    update neighborCells    ******/

	for (i = 0; i < deleted_cell->countNeighborCells; i++){

		//fprintf( stderr, "location: %i, neighbor: %i\n", ii, i);
		VoronoiCell *currentNeighbor = deleted_cell->neighborCells[i];
		
		// increase number of free neighborCells
		currentNeighbor->countFreeNeighborCells++;
			
		switch( currentNeighbor->getState()){

			case NONACTIVE:
			//fprintf( stderr, "INFO: NONACTIVE\n");
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}
	}


	/******    update neighborCells    ******/

	for (i = 0; i < deleted_cell->countExtendedNeighborCells; i++){

		VoronoiCell *currentNeighbor = deleted_cell->extendedNeighborhood[i];
		
		// increase number of free neighborCells
		currentNeighbor->countFreeExtendedNeighborCells++;
			
		switch( currentNeighbor->getState()){

			case NONACTIVE:
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}
	}
}
/*****************************************************************************/



void removeVesselAndUpdateSurrounding( ActionTree *actionList, AgentList* agentArray, VoronoiCell* deletedVessel, VoronoiDiagram* voronoiDiagram ){
//void addedVesselAndUpdateSurrounding(  ){
	int i;

	
	// UPDATE REMOVED VESSEL SITE
	fprintf(stderr, "\nUPDATE REMOVED VESSEL SITE\n");
	if( GetAgent( deletedVessel)->countLocations == 1)
		agentArray->deactivateAgent( GetAgent( deletedVessel));
	GetAgent( deletedVessel)->detach( deletedVessel);
								

	/******    update neighborCells    ******/

	//for(ii=0; ii<deleted_cell->countLocations; ii++)
	for (i = 0; i < deletedVessel->countNeighborCells; i++){

		//fprintf( stderr, "location: %i, neighbor: %i\n", ii, i);
		VoronoiCell *currentNeighbor = deletedVessel->neighborCells[i];
		
		// increase number of free neighborCells
		currentNeighbor->countFreeNeighborCells++;
			
		switch( currentNeighbor->getState()){

			case NONACTIVE:
			//fprintf( stderr, "INFO: NONACTIVE\n");
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}
	}


	/******    update neighborCells    ******/

	//for(ii=0; ii<deleted_cell->countLocations; ii++)
	for (i = 0; i < deletedVessel->countExtendedNeighborCells; i++){

		VoronoiCell *currentNeighbor = deletedVessel->extendedNeighborhood[i];
		
		// increase number of free neighborCells
		currentNeighbor->countFreeExtendedNeighborCells++;
			
		switch( currentNeighbor->getState()){

			case NONACTIVE:
			GetAgent( currentNeighbor)->actualize( actionList);
			break;
    
			default:
			break;
		}
	}
}
/*****************************************************************************/



void addDivisionAction( ActionList *actionList, Agent* cell, VoronoiDiagram* voronoiDiagram){

	if(cell->actionsInitialized == FALSE){
		initCellActions( cell);
	}

	addAction( actionList, cell->actions[0]);
}
/*****************************************************************************/



void addNecrosisAction( ActionList *actionList, Agent* cell, VoronoiDiagram* voronoiDiagram){

	if(cell->actionsInitialized == FALSE){
		initCellActions( cell);
	}

	addAction( actionList, cell->actions[ GetVoronoiCell( cell)->countNeighborCells + 2]);
}
/*****************************************************************************/

double heightCircleSegment( double radius, double height, 
                            double angle)
{
	return radius * (1. - cos( angle / 2.)) - height;
}


double angleCircleSegment( double radius, double height)
{
	//return radius * (1. - cos( angle / 2.)) = height;
	return 2. * acos( 1. - height / radius); 
}


double areaCircleSegment( double radius, double area, 
                          double angle)
{
	return pow( radius, 2) * (angle - sin( angle)) - area;
}
/*double angleCircleSegment( double radius, double height){
	z * z;
}*/


/*void actualizeProcessRates( VoronoiDiagram* voronoiDiagram, ActionList *actionList, double ***Oxygen_Concentration, double ***GrosseFactors_Concentration, double ***Glucose_Concentration){

	int i;
	Action* p_elem;


	//printf("INFO: select_actions()\n");


	// RESET rate LIST
	actionList->sum_of_prob = 0.;

	// SUM ACTION PROBABILITIES
	p_elem = actionList->head;
  
	for(i = 0; i < actionList->length; i++){

		// refresh runtime dependend actions probabilities
		switch( p_elem->type){


#if USE_VESSEL_GROWTH

			case VESSEL_GROWTH:
			//if( GrosseFactors_Concentration[p_elem->originalCell->mx][p_elem->originalCell->my][p_elem->originalCell->mz] > GROWTHFACTOR_THRESHOLD){
			//	p_elem->rate = VesselGrowthRate;
			//}else{
			//	p_elem->rate = 0.;
			//}
			if( p_elem->originalCell->growthfactors > GROWTHFACTOR_THRESHOLD){
				p_elem->rate = VesselGrowthRate;
			}else{
				p_elem->rate = 0.;
			}
			break;
#endif


#if USE_DIVISION

			case DIVISION:		
			p_elem->rate = GetGrowthRate( p_elem->originalCell) * (M_div + M_gro + 3.);
			//fprintf( stderr, "rate: %lf\n", p_elem->rate);
			//p_elem->rate = GetGrowthRate( p_elem->originalCell) * double(M_div + 1) * (SUBCELLULAR_COMPONENTS + 1.);
			//p_elem->rate *= double(M_div + 1);
			break;

#endif //USE_DIVISION


#if USE_GROWTH

			case GROWTH:		
			p_elem->rate = GetGrowthRate( p_elem->originalCell) * (M_div + M_gro + 3.);
			//p_elem->rate *= double(M_div + 1);
			break;

#endif //USE_DIVISION


#if USE_MIGRATION

			case MIGRATION:
			p_elem->rate = CellMigrationRate;
			break;

#endif //USE_MIGRATION


#if USE_NECROSIS

			case NECROSIS:
			p_elem->rate = GetDeathRate( p_elem->originalCell);
			//printf("INFO: Actualize NECROSIS rate: p = %lf\n", p_elem->rate);
			break;

#endif


#if USE_LYSIS

			case LYSIS:
			p_elem->rate = CellLysisRate * (double)p_elem->originalCell->necroticCellCount;
			break;

#endif

			default:
			fprintf( stderr, "ERROR: Not Identified Action (NIA) discovert!!! :) C'est grave ca!\n");
			exit( 0);

		}

		if(p_elem->rate < 0){
			fprintf(stderr, "\nERROR: p_elem->rate = %lf < 0!\n"\
			                "       p_elem->type = %i\n", p_elem->rate, p_elem->type);
			exit(0);                          
		}

		if( isnan(p_elem->rate) || isinf(p_elem->rate)){
			fprintf( stderr, "INFO: rate fucked up: %lf (cell state: %s, actions: %s)\n", p_elem->rate, cellTypeToString( p_elem->originalCell->state), actionTypeToString( p_elem->type));
			exit( 0);
		}

		// add actions rate to rate sum
		actionList->sum_of_prob += p_elem->rate;
		//printf("i:%i  type:%i, prob:%lf   -> actionList->sum_of_prob:%lf\n", i, p_elem->type, p_elem->rate, actionList->sum_of_prob);

		// set pointer to next actions element
		p_elem = p_elem->next;
	}
//	printf("INFO: actionList->sum_of_prob=%lf\n", actionList->sum_of_prob);


//	printActionList( actionList);
}*/
/*****************************************************************************/



Action* selectAction( ActionList *actionList, double* time){

	int i;
	double random, sum;

	// CHOOSE ACTION
	Action* p_elem = actionList->head;
  
	// randomly choose part of rate sum
	double temp_rand = myRand();
	random = temp_rand * actionList->sum_of_prob;

	// sum probabilities until randomly choosen part is reached
	sum = 0.;
	for(i = 0; i < actionList->length; i++){
		sum += p_elem->rate;
		if (random <= sum) {
			// calculate passed time
			random = myRand();
			*time += -log(1 - random) / actionList->sum_of_prob;

			return p_elem;
		}
		p_elem = p_elem->next;
	}
	printf("Something is going wrong in selectAction()\n");
	printf("INFO: actionList->sum_of_prob=%lf, random=%lf (%lf), sum=%lf\n", actionList->sum_of_prob, random, temp_rand, sum);
	printActionList( actionList);
	exit(0);
}
/*****************************************************************************/

VoronoiCell* growCell( VoronoiDiagram *voronoiDiagram, Agent* cell, int &shift, VoronoiCell **sourceLocation){
	
	if( VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD){
		sourceLocation[0] = cell->location[0];

		VoronoiCell *destination;
		if( VoronoiCell::SHIFT_TO_UNOCCUPIED)
			destination = voronoiDiagram->searchClosestUnoccupiedVoronoiCell( cell->location[0], VoronoiCell::extendedNeighborDepth, VoronoiCell::symbolicExtendedNeighborhoodSize, VoronoiCell::symbolicExtendedNeighborhood);
		else
			destination = voronoiDiagram->searchClosestFreeVoronoiCell( cell->location[0], VoronoiCell::extendedNeighborDepth, VoronoiCell::symbolicExtendedNeighborhoodSize, VoronoiCell::symbolicExtendedNeighborhood);

		if( destination!=0 && cell->location[0]->countFreeNeighborCells == 0)
			shift = TRUE;
		return destination;
	}else
	if( VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD){
		sourceLocation[0] = cell->location[0];
		VoronoiCell *destination;
		if( VoronoiCell::SHIFT_TO_UNOCCUPIED)
			destination = voronoiDiagram->searchClosestUnoccupiedVoronoiCell( cell->location[0], VoronoiCell::extendedNeighborDepth);
		else
			destination = voronoiDiagram->searchClosestFreeVoronoiCell( cell->location[0], VoronoiCell::extendedNeighborDepth);
		if( destination!=0 && cell->location[0]->countFreeNeighborCells == 0)
			shift = TRUE;
		return destination;
	}

	int i;

	// SAME CMPARTMENT
	if( cell->state == COMPARTMENT && cell->cellCount<cell->maxCellCount)
		return GetVoronoiCell( cell);

	// DIRECT NEIGHBOR
	int countFreeNeighbors = 0;
	for( i=0; i<cell->countLocations; i++)
		countFreeNeighbors += cell->location[i]->countFreeNeighborCells;

	int which_free_neigbor = (int) (myRand() * (double) countFreeNeighbors) + 1;
	for( i=0; i<cell->countLocations; i++){
		if( cell->location[i]->countFreeNeighborCells>0){
			for (int ii = 0; ii < cell->location[i]->countNeighborCells; ii++){
				if( cell->location[i]->neighborCells[ii]->isFree()){
					if(--which_free_neigbor == 0){

						return cell->location[i]->neighborCells[ii];
					}
				}
			}
		}
	}
	
	// EXTENDED NEIGHBORHOOD

	if( cell->state != NONACTIVE){
		int ii;
		double distance,
		       minDistance = 0;
		VoronoiCell* closestExtendedNeighborCell = NULL;
		VoronoiCell* closestLocation = NULL;

		for(i=0; i<cell->countLocations; i++){
			for(ii=0; ii<cell->location[i]->countExtendedNeighborCells; ii++)
				if( cell->location[i]->extendedNeighborhood[ii]->isFree()){
				distance = pow( cell->location[i]->extendedNeighborhood[ii]->position[0] - cell->location[i]->position[0], 2.) 
				         + pow( cell->location[i]->extendedNeighborhood[ii]->position[1] - cell->location[i]->position[1], 2.)
#if DIMENSIONS == 3
				         + pow( cell->location[i]->extendedNeighborhood[ii]->position[2] - cell->location[i]->position[2], 2.)
#endif
				         ;
				if( distance < minDistance || closestExtendedNeighborCell==NULL){
					minDistance = distance;
					closestExtendedNeighborCell = cell->location[i]->extendedNeighborhood[ii];
					closestLocation = cell->location[i];
				}
			}
		}
		
#if INITIAL_GROWTH_DECISION
		// SEARCH FREE VORONOI CELL
		if( closestExtendedNeighborCell==NULL){
			closestLocation = cell->location[0];
			closestExtendedNeighborCell = voronoiDiagram->searchClosestFreeVoronoiCell( cell->location[0]);
		}		
#endif		
		if( closestLocation==NULL)
			fprintf( stderr, "closestLocation==NULL\n");
		if( closestExtendedNeighborCell==NULL){

			for(i=0; i<cell->countLocations; i++){
				fprintf( stderr, "%ith EXTENDED NEIGHBOR COUNT (free: %i/%i)?\n", i+1, cell->location[i]->countFreeExtendedNeighborCells, cell->location[i]->countExtendedNeighborCells);
				for(ii=0; ii<cell->location[i]->countExtendedNeighborCells; ii++){
					if( cell->location[i]->extendedNeighborhood[ii]->isFree()){
						distance = pow( cell->location[i]->extendedNeighborhood[ii]->position[0], 2.) 
					                 + pow( cell->location[i]->extendedNeighborhood[ii]->position[1], 2.)
					                 #if(DIMENSIONS == 3)
					      	         + pow( cell->location[i]->extendedNeighborhood[ii]->position[2], 2.)
					      	         #endif
					      	         ;
						if( distance < minDistance || closestExtendedNeighborCell==NULL){
							minDistance = distance;
							closestExtendedNeighborCell = cell->location[i]->extendedNeighborhood[ii];
						}
					}
				}
			}

		}
		//fprintf( stderr, "closest EXTENDED NEIGHBOR to location %i is %i\n", closestLocation->index, closestExtendedNeighborCell->index);
			
		shift = TRUE;
		*sourceLocation = closestLocation;
		//fprintf(stderr,"sourceLocation: %i!\n",   sourceLocation->index);
		return closestExtendedNeighborCell;
		
	}else{
		fprintf( stderr, "ERROR: NONACTIVE cell (%i) tries to divide!!\n", cell->index);
		fprintf( stderr, "corresponding Voronoi cell: %i\n", cell->location[0]->index);
		fprintf( stderr, "corresponding action: %p\n", cell->actions[INDEX_GROWTH]);
		exit( 0);
	}
	
	
	fprintf( stderr, "Could not attach cell to any neighboring Voronoi cell: %i!!!\n", cell->countLocations);
	fprintf( stderr, "cell state: %s, free neighbors: %i/%i, free extended neighbors: %i/%i !!!\n", cellTypeToString( cell->state), GetVoronoiCell( cell)->countFreeNeighborCells, GetVoronoiCell( cell)->countNeighborCells, GetVoronoiCell( cell)->countFreeExtendedNeighborCells, GetVoronoiCell( cell)->countExtendedNeighborCells);
	return NULL;
}


/*****************************************************************************/
Agent* divideCell( ActionList *actionList, Action* action, VoronoiDiagram* voronoiDiagram){

	int shift_flag = 0;
	//Agent** path_for_shift = NULL;
	//int path_length = 0;
	Agent* child_cell;
	Agent* parent_cell;
	//double division_time = 0.0;
	//double mutation_variance = 0.0,
	//       random_div_rate = 0.0;

	parent_cell = action->originalCell;
	action->internalState = 0;

	// get child_cell
	child_cell = GetAgent( getChildCell(parent_cell, voronoiDiagram, &shift_flag, actionList));
	//fprintf( stderr, "VoronoiCell %i(%i->%i) is %s(%s->%s)\n", child_cell->index, GetVoronoiCell( child_cell)->index, GetAgent( GetVoronoiCell( child_cell))->index, cellTypeToString( child_cell->state), cellTypeToString( GetVoronoiCell( child_cell)->getState()), cellTypeToString( GetAgent( GetVoronoiCell( child_cell))->state));
	//fprintf(stderr,"DIVISION -> (cell nr %i state %i cells %i)\n", child_cell->nr, child_cell->state, child_cell->cellCount);
	//printf("div_cell()\n");

	// shift
	/*if(shift_flag){
		//printf("NEW_getShiftPath()\n");
		//exit( 0);
		path_for_shift = getShiftPath(parent_cell,child_cell, &path_length);
		//printf("...done\n");
	}*/

#if MULTISCALE
	if( parent_cell->cellCount != 0 && child_cell->state == FREE){
#else
	if( parent_cell->state == ACTIVE && child_cell->state == FREE){
#endif


	}
	else{
		fprintf(stderr, "\nWARNING: Cell division of a(n) %s cell to a(n) %s cell not possible!\n", cellTypeToString( parent_cell->state), cellTypeToString( child_cell->state));
		exit(0);
	}
  
	return child_cell;
}
/****************************************************************************/


const char *actionTypeToString( int type){
	switch( type){
		case GROWTH:
		return "GROWTH";
		
		case DIVISION:
		return "DIVISION";
		
		case MIGRATION:
		return "MIGRATION";
		
		case NECROSIS:
		return "NECROSIS";
		
		case LYSIS:
		return "LYSIS";

		case VESSEL_GROWTH:
		return "VESSEL_GROWTH";

		default:
		return "UNKNOWN";
	}
	return "";
}
/****************************************************************************/



double getDistance( VoronoiCell *a, VoronoiCell *b)
{
	int i;
	double dist = 0.;
	
	for( i=0; i<DIMENSIONS; i++)
		dist += pow( a->position[i] - b->position[i], 2);
		
	return sqrt( dist);
}
/****************************************************************************/



VoronoiCell* getChildCell( Agent *parent_cell, VoronoiDiagram *voronoiDiagram, int *shift_flag, ActionList *actionList){
	int i, which_free_neigbor;
	//Agent* min_radius_cell = parent_cell;
	//int free_kr_neighbor_found = FALSE;
  
	//double min_radius = 0, radius; 
  
  
  
  /* search in neighbor-list for free voronoipoints */
#if IMPORT_COARSE_GRAINED_BEHAVIOR == 1    
	int temp_which_free_neigbor;

	double randomDivisionRateSum = myRand() * (
		parent_cell->action[0]->prob_of_action
		- (double)parent_cell->cellCount 
		/ (double)parent_cell->nrn 
		* IntracompartmentDivisionRateAverages[ 
			getDiscreteIndex( Resolution, 0., 1., (double)parent_cell->cellCount/(double)parent_cell->maxCellCount) 
		]
	);
  					
	for( i=0; i<parent_cell->nrn; i++){
		randomDivisionRateSum -= (double)parent_cell->cellCount / parent_cell->nrn
		* IntercompartmentDivisionRateAverages
        	[ getDiscreteIndex( Resolution, 0., 1., (double)parent_cell->cellCount/(double)parent_cell->maxCellCount) ]
        	[ getDiscreteIndex( Resolution, 0., 1., (double)parent_cell->neighbors[i]->cellCount/(double)parent_cell->neighbors[i]->maxCellCount) ];
		
		if(  randomDivisionRateSum<0.)
			return parent_cell->neighbors[i];
	}

#else

	int temp_which_free_neigbor = (int) (myRand() * (double) GetVoronoiCell(parent_cell)->countFreeNeighborCells) + 1;
	which_free_neigbor = temp_which_free_neigbor;
	//printf("free neighbors = %d, which_free_neigbor = %d\n", GetVoronoiCell( parent_cell)->countFreeNeighborCells, which_free_neigbor);
	for (i = 0; i < GetVoronoiCell( parent_cell)->countNeighborCells; i++){
		if ( GetNeighborCell( parent_cell, i)->isFree()){
			if(--which_free_neigbor == 0){
				//fprintf( stderr, "VoronoiCell %i is %s\n", GetNeighborCell( parent_cell, i)->index, cellTypeToString( GetNeighborCell( parent_cell, i)->getState()));
				return GetNeighborCell( parent_cell, i);
			}
		}
	}

#endif

	printf("WARNING: Necessary to shift!\n");
	exit(0);

	/* search in radius list for free voronoipoints, take the nearest */
/*	shift_flag[0] = 1;
	for (i = 0; i < parent_cell->location[0]->countExtendedNeighborCells; i++){
		if ( GetExtendedNeighborCell( parent_cell, i)->getState() == FREE){
			radius = getDistance( parent_cell, GetExtendedNeighborAgent( parent_cell, i));
			if(free_kr_neighbor_found == FALSE){
				min_radius = radius;
				min_radius_cell = GetExtendedNeighborAgent( parent_cell, i);
				free_kr_neighbor_found = TRUE;
			}
			if(radius < min_radius){
				min_radius = radius;
				min_radius_cell = GetExtendedNeighborAgent( parent_cell, i);
			}
		}
	}

	if(free_kr_neighbor_found){
		return min_radius_cell;
	}
	else{
		printf("Something is going wrong in getChildCell-function!\n");
		printf("Division of cell %i, state %s, cells %i (tumor:%i, necrotic:%i), free neighbors %i, free k neighbors %i\n", 
			GetVoronoiCell(parent_cell)->index, 
			cellTypeToString( parent_cell->state), 
			parent_cell->cellCount, 
			parent_cell->tumorCellCount, 
			parent_cell->necroticCellCount, 
			GetVoronoiCell( parent_cell)->countFreeNeighborCells, 
			GetVoronoiCell( parent_cell)->countFreeExtendedNeighborCells);
		printf("free neighbors = %d, temp_which_free_neigbor = %d, which_free_neigbor = %d\n", GetVoronoiCell( parent_cell)->countFreeNeighborCells, temp_which_free_neigbor, which_free_neigbor);
		int count_free_neighbors = 0;
		for (i = 0; i < GetVoronoiCell( parent_cell)->countNeighborCells; i++){
			printf("%i. neighbor: cell %i, state %s\n", i+1, GetVoronoiCell( parent_cell)->neighborCells[i]->index, cellTypeToString( GetAgent( GetVoronoiCell( parent_cell)->neighborCells[i])->state));
			
			if( GetAgent( GetVoronoiCell( parent_cell)->neighborCells[i])->state == FREE)
				count_free_neighbors ++;
		}
		for (i = 0; i < parent_cell->location[0]->countExtendedNeighborCells; i++)
			printf("%i. reachable neighbor: cell %i, state %s\n", i+1, GetExtendedNeighborCell( parent_cell, i)->index, cellTypeToString( GetExtendedNeighborCell( parent_cell, i)->getState()));
		printf("free neighbors = %d ? = %d = count_free_neighbors\n", GetVoronoiCell( parent_cell)->countFreeNeighborCells, count_free_neighbors);

		printActionList( actionList);

		exit(0);
		//return 0;
	}*/
}
/*****************************************************************************/



VoronoiCell** getShiftPath( Agent *expandingAgent, VoronoiCell **end_point, VoronoiCell *sourceLocation, int &pathLength){
   
	double path_vector_abs;
	double neighbor_vector_abs;
	VoronoiCell* temp_point;
	VoronoiCell* candidate = NULL;
	int i, length = 0;
	VoronoiCell** path = NULL;
	//int    min_neighbor_index = 0;
	double min_neighbor_vector_abs;

	//printf( "Shifting!!\n");
      

	// find start segment of agent
	/*VoronoiCell *start_point = expandingAgent->location[0];
	double min_path_vector_abs = getDistance( end_point, expandingAgent->location[0]);
	for( i=1; i<expandingAgent->countLocations; i++){
		//fprintf( stderr, "Shift starts from Voronoi cell %i?\n", expandingAgent->location[i]->index);
		path_vector_abs = getDistance( end_point, expandingAgent->location[i]);
		if( path_vector_abs < min_path_vector_abs){
			min_path_vector_abs = path_vector_abs;
			start_point = expandingAgent->location[i];
		}
	}
	path_vector_abs = min_path_vector_abs;
	//fprintf( stderr, "Shift from Voronoi cell %i to %i\n", start_point->index, end_point->index);


	temp_point = start_point;
	*/
	temp_point = sourceLocation;
	path_vector_abs = getDistance( *end_point, temp_point);
	// initialize path_vector data
	//path_vector_abs = getDistance( end_point, start_point);

	//fprintf( stderr, "[ %i (%i)", sourceLocation->index, GetAgent(sourceLocation)->index);

	do{
		min_neighbor_vector_abs = path_vector_abs;
		for(i = 0; i < temp_point->countNeighborCells; i++){
		
			// take vector from neighbor of actual point to end point
			neighbor_vector_abs = getDistance( temp_point->neighborCells[i], *end_point);
	//		fprintf( stderr, "%i. neighbor of %i is %i, distance = %lf", i+1, temp_point->index, temp_point->neighborCells[i]->index, neighbor_vector_abs);

			// search nearest neighbor to end point  
			if(min_neighbor_vector_abs > neighbor_vector_abs){
				min_neighbor_vector_abs = neighbor_vector_abs;
				//min_neighbor_index = i;
				candidate = temp_point->neighborCells[i];
	//			fprintf( stderr, " => candidate!");
			}
	//		fprintf( stderr, "\n");
			
	//		for( int ii = 0; ii < temp_point->neighborCells[i]->countNeighborCells; ii++)
	//			fprintf( stderr, "%i%s ", temp_point->neighborCells[i]->neighborCells[ii]->index, (temp_point->neighborCells[i]->neighborCells[ii]==end_point?"!!!!":""));
	//		fprintf( stderr, "\b\b\n");
		}
    

		// add next path element
		length++;
		path = ( VoronoiCell**) realloc( path, sizeof( VoronoiCell*) * length);
		//path[ length-1] = temp_point = GetAgent( GetVoronoiCell( temp_point)->neighborCells[min_neighbor_index]);   
		path[ length-1] = temp_point = candidate;  
		//printf( "Next Point: %d\n", path[length-1]->index);

		if( length>2 && path[length-1] == path[length-3]){
			fprintf( stderr, "new path element %i is already contained in the path\n", temp_point->index);
			exit( 0);
		}
		if( length>1 && path[length-1] == path[length-2]){
			for(i = 0; i < temp_point->countNeighborCells; i++)
				fprintf( stderr, "%i. neighbor of %i is %i\n", i+1, temp_point->index, temp_point->neighborCells[i]->index);
			exit( 0);
		}
	//	fprintf( stderr, ", %i (%i)", temp_point->index, (!temp_point->isFree()?GetAgent(temp_point)->index:-1));
	}while( temp_point != *end_point && !temp_point->isFree());
	//fprintf( stderr, "\n");
	
	//if(temp_point->isFree()){
	//	fprintf( stderr, "test: path length = %i\n", length);
		*end_point = temp_point;
	//}

	//fprintf( stderr, "path length = %i\n", length);

	//length++;

  // copy path to allocated memory
  /*if((path =( int *)calloc(length,sizeof(int)))==NULL){
    fprintf(stderr,"Cannot callocate memory for new shift path! Exiting.");
  }*/
   	
  /*for(i = 0; i < length; i++){
    path[i] = temp_path[i];
  }*/

  // returning results
  pathLength = length;
  return path;
}
/*****************************************************************************/



/*void shiftPath( Agent** path, int length){

	int i = length;

	VoronoiCell* tempLocation         = path[0]->location[0];
	int tempCountDirectNeighbors        = path[0]->countDirectNeighbors,	
	    tempCountDirectFreeNeighbors    = path[0]->countDirectFreeNeighbors,	
	    tempCountReachableNeighbors     = path[0]->location[0]->countExtendedNeighborCells,	
	    tempCountReachableFreeNeighbors = path[0]->countReachableFreeNeighbors;
	Agent** tempReachableNeighbors   = path[0]->location[0]->extendedNeighborhood;

	for(i = 0; i < length - 1; i++){
		//GetAgent( GetVoronoiCell( path[i-1])) = path[i];

		//path[i+1]->location->data = path[i];

		// exchange of data
		path[i]->location[0] = path[i+1]->location[0];
		path[i]->state = path[i+1]->state;	
		path[i]->countDirectNeighbors = path[i+1]->countDirectNeighbors;	
		path[i]->countDirectFreeNeighbors = path[i+1]->countDirectFreeNeighbors;
		path[i]->location[0]->countExtendedNeighborCells = path[i+1]->location[0]->countExtendedNeighborCells;
		path[i]->countReachableFreeNeighbors = path[i+1]->countReachableFreeNeighbors;
		path[i]->location[0]->extendedNeighborhood = path[i+1]->location[0]->extendedNeighborhood;

		path[i]->location[0]->data = path[i];
	}

	path[i]->location[0] = tempLocation;
	path[i]->state = FREE;	
	path[i]->countDirectNeighbors = tempCountDirectNeighbors;	
	path[i]->countDirectFreeNeighbors = tempCountDirectFreeNeighbors;
	path[i]->location[0]->countExtendedNeighborCells = tempCountReachableNeighbors;
	path[i]->countReachableFreeNeighbors = tempCountReachableFreeNeighbors;
	path[i]->location[0]->extendedNeighborhood = tempReachableNeighbors;

	path[i]->location[0]->data = path[i];

}*/

void shiftPath( VoronoiCell** path, int length){

	int i;
	
	//if( length>2){
	//	fprintf(stderr,"path length %i!\n", length);
		//exit( 0);		
	//}
	
	/*if( path[length - 1]->agent!=NULL && path[length - 1]->getState() == FREE){
		fprintf(stderr,"Last Element of shift path (%i) containes free agent (%i)!\n", path[length - 1]->index, GetAgent( path[length - 1])->index);
		exit( 0);
	}*/
	
	for(i = 0; i < length - 1; i++){
		//GetAgent( GetVoronoiCell( path[i-1])) = path[i];

		//path[i+1]->location->data = path[i];

		//fprintf(stderr,"Have to shift: (%i/%i) %i -> %i!\n", i+1, length-1, path[length - i - 2]->index, path[length - i - 1]->index);
		if( GetAgent( path[length - i - 2])==NULL){
			fprintf(stderr,"Voronoi cell %i containes no agent!\n", path[length - i - 2]->index);
			fprintf(stderr,"path length %i!\n", length);
			exit( 0);	
		}
		
		// exchange of data
		//path[length - i - 1]->agent = path[length - i - 2];
		Agent *tempAgent = GetAgent( path[length - i - 2]);
		tempAgent->detach( path[length - i - 2]);
		tempAgent->attach( path[length - i - 1]);
	}
}
/*****************************************************************************/



void printActionList( ActionList *p_list){
	int countActions = 0;

	Action* test_prob;
	//Action* test_head;

	printf("\n\nProb_liste:\n[ ");
	test_prob = p_list->head;
	//test_head = p_list->head->prev;

	//for( i = 0; i < p_list->length; i++){
	if(p_list->length!=0){
		do{
			countActions++;

			printf("[ type: %s, prob: %lf, cell: nr %i type %s]\n", 
				actionTypeToString( test_prob->type),
				test_prob->rate,
				GetVoronoiCell( test_prob->originalCell)->index,
				cellTypeToString( test_prob->originalCell->state));
			if( test_prob->originalCell->state == 2 && test_prob->rate == 0){
				fprintf( stderr, "ERROR: NONACTIVE cell (%i) should not be able to divide!!!\n", GetVoronoiCell( test_prob->originalCell)->index);
				//exit( 0);
			}

		test_prob = test_prob->next;

		}while( test_prob != p_list->head);

			printf("[ (head) type: %s, prob: %lf, cell: nr %i type %s]\n", 
				actionTypeToString( test_prob->type),
				test_prob->rate,
				GetVoronoiCell( test_prob->originalCell)->index,
				cellTypeToString( test_prob->originalCell->state));
	}

	printf("Summe: %lf\n",p_list->sum_of_prob);

	if(countActions != p_list->length){
		printf("ERROR: number of actions (p_list->length=%i) is incorrect: %i counted!!!\n(PRESS ENTER TO EXIT)\n", p_list->length, countActions);
		getchar();
		exit(0);
	}
	fprintf( stderr, "length: %i elements\n", p_list->length);
}
/*****************************************************************************/


//#define Heavyside( a)    (a < 0 ? 0 : 1)
//#define Heavyside( a, c=0) ( a==0 ? c : (a < 0 ? 0 : 1))
inline float Heavyside( float a, float c=0){
	return ( a==0 ? c : (a < 0 ? 0 : 1));
}

double GetDeathRate( Agent * agent){

	// APOPTOSIS
	double rate = CellApoptosisRate;

	//apoptosisRate = 0.01 + (agent->location[0]->ecm > 0.07 ? 0.:0.);//0.002/4.; //REF:0.002;

	// GLUCOSE & OXYGEN DEPENDENT NECROSIS
/*#if (MONOD_KINETICS == 5)
	if( 	//true ||
			GiveMeTheATP( agent) < VoronoiCell::ATP_THRESHOLD_DEATH ||
			agent->getGlucose()  < VoronoiCell::GLUCOSE_THRESHOLD_DEATH ||
		    agent->getOxygen()   < VoronoiCell::OXYGEN_THRESHOLD_DEATH ||
		    agent->getGlucose() * agent->getOxygen() < VoronoiCell::GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_DEATH
	){
			rate += CellNecrosisRate; // cells/s

	}else if( agent->location[0]->lactate >= VoronoiCell::LACTATE_THRESHOLD_DEATH){
			//rate += CellNecrosisRate * MIN( 1, agent->location[0]->lactate / VoronoiCell::LACTATE_THRESHOLD_DEATH); // cells/s
			rate += CellNecrosisRate * pow(agent->location[0]->lactate,3.) / (pow(agent->location[0]->lactate,3.) + pow(VoronoiCell::LACTATE_THRESHOLD_DEATH,3)); // cells/s
	}else{
		agent->actions[INDEX_NECROSIS]->internalState = 0;
	}
#endif //NO MONOD_KINETICS*/

	float necrotic = 0;
	necrotic = MAX( necrotic, Heavyside( VoronoiCell::ATP_THRESHOLD_DEATH     - GiveMeTheATP( agent), 0));
	necrotic = MAX( necrotic, Heavyside( VoronoiCell::GLUCOSE_THRESHOLD_DEATH - agent->getGlucose(), 0));
	necrotic = MAX( necrotic, Heavyside( VoronoiCell::OXYGEN_THRESHOLD_DEATH - agent->getOxygen(), 0));
	necrotic = MAX( necrotic, Heavyside( VoronoiCell::GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_DEATH - agent->getGlucose() * agent->getOxygen(), 0));

	// HEAVYSIDE
	//necrotic = MAX( necrotic, Heavyside( agent->location[0]->lactate - VoronoiCell::LACTATE_THRESHOLD_DEATH, 1));
	// LINEAR
	necrotic = MAX( necrotic, agent->location[0]->lactate/100.);
	// HILL
	float n=2;
	necrotic = MAX( necrotic, pow(agent->location[0]->lactate,n)/(pow(agent->location[0]->lactate,n)+pow(VoronoiCell::LACTATE_THRESHOLD_DEATH,n)));

	necrotic = MAX( necrotic, Heavyside( agent->location[0]->waste - VoronoiCell::WASTE_THRESHOLD_DEATH, 1));
	rate += CellNecrosisRate * necrotic; // cells/s
	/*		*
			agent->getGlucose()  < VoronoiCell::GLUCOSE_THRESHOLD_DEATH ||
		    agent->getOxygen()   < VoronoiCell::OXYGEN_THRESHOLD_DEATH ||
		    agent->getGlucose() * agent->getOxygen() < VoronoiCell::GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_DEATH
	){


	}else if( agent->location[0]->lactate >= VoronoiCell::LACTATE_THRESHOLD_DEATH){
			//rate += CellNecrosisRate * MIN( 1, agent->location[0]->lactate / VoronoiCell::LACTATE_THRESHOLD_DEATH); // cells/s
			rate += CellNecrosisRate * pow(agent->location[0]->lactate,3.) / (pow(agent->location[0]->lactate,3.) + pow(VoronoiCell::LACTATE_THRESHOLD_DEATH,3)); // cells/s
	}else{
		agent->actions[INDEX_NECROSIS]->internalState = 0;
	}*/



	//
	if( agent->state == COMPARTMENT){
		rate *= //(double)( agent->isGrowing() ? agent->growingTumorCellCount //: 0);
#if DIMENSIONS == 1
				MIN( agent->growingTumorCellCount,	agent->maxCellCount * (1.-Interpolation::getDiscritizedVolumeAboveThreshold1D( GetVoronoiCell(agent)->index, 0, THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)));
#elif DIMENSIONS == 2
				MIN( agent->growingTumorCellCount,	agent->maxCellCount * (1.-Interpolation::getDiscritizedVolumeAboveThreshold2D( GetVoronoiCell(agent)->index, 0, THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)));
#elif DIMENSIONS == 3
				MIN( agent->growingTumorCellCount,	agent->maxCellCount * (1.-Interpolation::getDiscritizedVolumeAboveThreshold3D( GetVoronoiCell(agent)->index, 0, THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)));
#endif
		//: 0);
	}

	return rate;
}
/*****************************************************************************/

bool IsQuiescent( Agent * agent)
{
	if( agent){
		if( agent->actions[INDEX_GROWTH]->internalState == 0 &&
			( 	GiveMeTheATP( agent) < VoronoiCell::ATP_THRESHOLD_QUIESCENCE ||
				agent->getGlucose()  < VoronoiCell::GLUCOSE_THRESHOLD_QUIESCENCE ||
				agent->getOxygen()   < VoronoiCell::OXYGEN_THRESHOLD_QUIESCENCE ||
				agent->getGlucose() * agent->getOxygen() < VoronoiCell::GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_QUIESCENCE ||
				agent->location[0]->lactate >= VoronoiCell::LACTATE_THRESHOLD_QUIESCENCE ||
				agent->location[0]->waste >= VoronoiCell::WASTE_THRESHOLD_QUIESCENCE ||
				agent->waste >= Agent::WASTE_THRESHOLD_QUIESCENCE
			)
		)
			return true;
	}
	return false;
}

double GetGrowthRate( Agent * agent){
	
	double rate = 0.;

	if( agent== NULL){
		fprintf( stderr, "ERROR: specified agent doesn't exist!!!\n");
		exit( 0);
	}

	if( agent->actions== NULL){
		fprintf( stderr, "ERROR: division of specified agent is not allocated!!!\nactions initialized: %i\n", agent->actionsInitialized);
		fprintf( stderr, "INFO:  agent %i: cellCount %i, type %s\n", GetIndex( agent), agent->cellCount, cellTypeToString( agent->state));

		exit( 0);
	}


	// MAXIMAL GROWTH RATE OF A SINGLE CELL
	#if (MONOD_KINETICS==5)
	if( agent->state != COMPARTMENT){
		if( IsQuiescent( agent)
			/*agent->actions[INDEX_GROWTH]->internalState == 0 &&
			( 	GiveMeTheATP( agent) < VoronoiCell::ATP_THRESHOLD_QUIESCENCE ||
				agent->getGlucose()  < VoronoiCell::GLUCOSE_THRESHOLD_QUIESCENCE ||
				agent->getOxygen()   < VoronoiCell::OXYGEN_THRESHOLD_QUIESCENCE ||
				agent->getGlucose() * agent->getOxygen() < VoronoiCell::GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_QUIESCENCE ||
				agent->location[0]->lactate >= VoronoiCell::LACTATE_THRESHOLD_QUIESCENCE ||
				agent->location[0]->waste >= VoronoiCell::WASTE_THRESHOLD_QUIESCENCE
			)*/
		){
			rate = 0.;
			//fprintf(stderr, "TEST (agent: %i, oxygen: %lf)\n", agent->index, GetVoronoiCell(agent)->oxygen);
			//agent->actions[INDEX_GROWTH]->internalState = 0;
		}else
	#endif //MONOD_KINETICS

	//rate = agent->divide * MaxCellDivisionRate;//MAX_SPECIFIC_GROWTH_RATE;
	if( agent->divide)
		rate = MaxCellDivisionRate;//MAX_SPECIFIC_GROWTH_RATE;
	else
		rate = Agent::ReentranceRate;

	if(	agent->location[0]->waste >= Agent::WASTE_THRESHOLD_SLOWED_GROWTH)
		//rate/=2.;
		rate*=0.45;
	if( agent->getOxygen() < 0.07)
		rate/=3;



	// NUMBER OF GROWING CELLS
	}else{
	//if( agent->state == COMPARTMENT){
		rate = MaxCellDivisionRate//MAX_SPECIFIC_GROWTH_RATE
				* (double)( agent->isGrowing() ?
		//		agent->growingTumorCellCount
				MIN( agent->growingTumorCellCount,
//					agent->maxCellCount * Interpolation::getRectangleAreaAboveThreshold3D( GetVoronoiCell(agent)->index, 0, THRESHOLD_QUIESCENCE_GLUCOSE_OXYGEN))
#if DIMENSIONS == 1
					floor( agent->maxCellCount * Interpolation::getDiscritizedVolumeAboveThreshold1D( GetVoronoiCell(agent)->index, 0, THRESHOLD_QUIESCENCE_GLUCOSE_OXYGEN)))
#endif
#if DIMENSIONS == 2
					( agent->maxCellCount * Interpolation::getDiscritizedVolumeAboveThreshold2D( GetVoronoiCell(agent)->index, 0, THRESHOLD_QUIESCENCE_GLUCOSE_OXYGEN)))
#endif
#if DIMENSIONS == 3
					/*ceil*/( agent->maxCellCount * Interpolation::getDiscritizedVolumeAboveThreshold3D( GetVoronoiCell(agent)->index, 0, THRESHOLD_QUIESCENCE_GLUCOSE_OXYGEN)))
#endif
					: 0);
	}

	// NUTRIENT LIMITATION
	#if (MONOD_KINETICS == 1)

		rate *= GetVoronoiCell(agent)->oxygen  / (/*Ko*/ MONOD_PARAMETER_GROWTH + GetVoronoiCell(agent)->oxygen);
		rate *= GetVoronoiCell(agent)->glucose / (/*Kg*/ MONOD_PARAMETER_GROWTH + GetVoronoiCell(agent)->glucose);
		if( isnan(rate) || isinf(rate)){
			fprintf( stderr, "INFO: rate fucked up: %lf\n", rate);
			exit( 0);
		}
		//fprintf( stderr, "INFO: rate: %lf\n", rate);

	#endif //MONOD_KINETICS

	#if MONOD_KINETICS==2
		//rate = ( agent->oxygen < R_CELL_GROWTH_THRESHOLD ? 0. : R_MAX_SPECIFIC_GROWTH_RATE);
		rate = ( agent->oxygen < R_CELL_GROWTH_THRESHOLD ? 0. : MaxCellDivisionRate);
	#endif

	#if MONOD_KINETICS==3
		/*rate *= (1. - 
		        (1. - pow( GetVoronoiCell(agent)->oxygen, MONOD_PARAMETER_GROWTH_OXYGEN_N)/(pow( GetVoronoiCell(agent)->oxygen, MONOD_PARAMETER_GROWTH_OXYGEN_N) + MONOD_PARAMETER_GROWTH_OXYGEN))*
		        (1. - pow( GetVoronoiCell(agent)->glucose, MONOD_PARAMETER_GROWTH_GLUCOSE_N)/(pow( GetVoronoiCell(agent)->glucose, MONOD_PARAMETER_GROWTH_GLUCOSE_N) + MONOD_PARAMETER_GROWTH_GLUCOSE)) );*/
		//printf( "%lf %lf %lf\n", GetVoronoiCell(agent)->oxygen, GetVoronoiCell(agent)->glucose, rate);
	#endif
		//fprintf( stderr, "INFO: rate: %lf\n", rate);

	return rate;
}
/*****************************************************************************/



double GetDivisionRate( Agent * agent)
{
	double rate = 0.;

	if( agent== NULL){
		fprintf( stderr, "ERROR: specified agent doesn't exist!!!\n");
		exit( 0);
	}

	if( agent->actions== NULL){
		fprintf( stderr, "ERROR: division of specified agent is not allocated!!!\nactions initialized: %i\n", agent->actionsInitialized);
		fprintf( stderr, "INFO:  agent %i: cellCount %i, type %s\n", GetIndex( agent), agent->cellCount, cellTypeToString( agent->state));

		exit( 0);
	}


	// MAXIMAL GROWTH RATE OF A SINGLE CELL
	/*#if (MONOD_KINETICS==5)
	if( agent->actions[INDEX_DIVISION]->internalState == 0 && 
	    ( GiveMeTheATP( agent) < THRESHOLD_QUIESCENCE_ATP || 
	      agent->getGlucose()  < THRESHOLD_QUIESCENCE_GLUCOSE || 
	      agent->getOxygen()   < THRESHOLD_QUIESCENCE_OXYGEN || 
	      agent->getGlucose() * agent->getOxygen() < THRESHOLD_QUIESCENCE_GLUCOSE_OXYGEN
	    )
	  )
		rate = 0.;
	else
	#endif*/ //MONOD_KINETICS
		rate = MaxCellDivisionRate;//MAX_SPECIFIC_GROWTH_RATE;

	if(	agent->location[0]->waste >= Agent::WASTE_THRESHOLD_SLOWED_GROWTH)
		//rate/=2.;
		rate*=0.45;
	if( agent->getOxygen() < 0.07)
		rate/=3; //  3


	// NUTRIENT LIMITATION
	#if (MONOD_KINETICS==1)

		rate *= GetVoronoiCell(agent)->oxygen  / (/*Ko*/ MONOD_PARAMETER_GROWTH + GetVoronoiCell(agent)->oxygen);
		rate *= GetVoronoiCell(agent)->glucose / (/*Kg*/ MONOD_PARAMETER_GROWTH + GetVoronoiCell(agent)->glucose);
		if( isnan(rate) || isinf(rate)){
			fprintf( stderr, "INFO: rate fucked up: %lf\n", rate);
			exit( 0);
		}
		//fprintf( stderr, "INFO: rate: %lf\n", rate);

	#endif //MONOD_KINETICS

	#if MONOD_KINETICS==2
		//rate = ( agent->oxygen < R_CELL_GROWTH_THRESHOLD ? 0. : R_MAX_SPECIFIC_GROWTH_RATE);
		rate = ( agent->oxygen < R_CELL_GROWTH_THRESHOLD ? 0. : MaxCellDivisionRate);
	#endif

	#if MONOD_KINETICS==3
		//rate *= (1. - (1. - GetVoronoiCell(agent)->oxygen/(GetVoronoiCell(agent)->oxygen + MONOD_PARAMETER_GROWTH_OXYGEN))*(1. - GetVoronoiCell(agent)->glucose/(GetVoronoiCell(agent)->glucose + MONOD_PARAMETER_GROWTH_GLUCOSE)));
	#endif
	
	// NUMBER OF DIVIDING CELLS
	//fprintf( stderr, "%i\n", ( agent->isDividing() ? agent->tumorVolume - agent->tumorCellCount : 0));
	if( agent->state == COMPARTMENT)
	rate *= (double)( agent->isDividing() ? (int)(agent->dividingTumorCellCount) : 0);
	//fprintf( stderr, "isDividing()? : %i\n", agent->isDividing());


	return rate;
	
}
/*****************************************************************************/



double GetDecayRateGrowthFactors( Agent * agent){
	//if( Case == 1)
		return 0.;
	//else
	//	return GetGrowthRate( agent) / CELL_GROWTH_YIELD_OXYGEN + MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN;
}
/*****************************************************************************/



double GetConsumptionRateGrowthFactors( Agent * agent){
	//if( Case == 1)
		return 0.;
	//else
	//	return GetGrowthRate( agent) / CELL_GROWTH_YIELD_OXYGEN + MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN;
}
/*****************************************************************************/



double GetConsumptionRateOxygen( VoronoiCell * cell){
	//if( Case == 1)
	//	return 0.;
	//else{
		if( cell->getState() == ACTIVE)
			return GetGrowthRate( GetAgent(cell)) / CELL_GROWTH_YIELD_OXYGEN +  GetAgent(cell)->cellCount * MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN;
		else
			return 0.;
	//}
}
/*****************************************************************************/



double GetConsumptionRateGlucose( VoronoiCell * cell){
	//if( Case == 1)
	//	return 0.;
	//else
		if( cell->getState() == ACTIVE)
			return GetGrowthRate( GetAgent(cell)) / CELL_GROWTH_YIELD_OXYGEN + GetAgent(cell)->cellCount * MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE;
		else
			return 0.;
}
/*****************************************************************************/

#if SCHALLER == true
double GiveMeTheOxygenRate(VoronoiCell *thecell)
	{
	if(thecell->getState() != NONACTIVE && thecell->getState() != ACTIVE)
		{return 0.;}

	return Schaller_OxygenUptake;
	}
	
double GiveMeTheGlucoseRate(VoronoiCell *thecell)
	{
	if(thecell->getState() != NONACTIVE && thecell->getState() != ACTIVE)
		return 0.;

	return Schaller_OxygenUptake;
	}
#endif
/*****************************************************************************/

#if JAGIELLA == true
double GiveMeTheOxygenRate(VoronoiCell *thecell)
	{
	//if(thecell->getState() != NONACTIVE && thecell->getState() != ACTIVE)
	//	{return 0.;}

	double oxygenUptake = MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR
        * 1./(thecell->oxygen + Agent::NICK_O_CRITICAL_OXY)*(Agent::NICK_O_MAX - (Agent::NICK_O_MAX - Agent::NICK_O_MIN)*(thecell->glucose/(thecell->glucose + Agent::NICK_O_CRITICAL_GLU)));

	// rescale from median to minimal cell volume
	oxygenUptake *= 2/3.;

	if( thecell->getState() == COMPARTMENT)
		oxygenUptake *= (double)GetAgent(thecell)->growingTumorCellCount/(double)GetAgent(thecell)->maxCellCount;

	return oxygenUptake;
	}

double GiveMeTheOxygenRate(Agent *thecell)
	{
	//if( cell->getState() != COMPARTMENT)
	//if(thecell->state != NONACTIVE && thecell->state != ACTIVE)
	//	{return 0.;}

	double oxygenUptake = MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR
        * thecell->getOxygen()/(thecell->getOxygen() + Agent::NICK_O_CRITICAL_OXY)*(Agent::NICK_O_MAX - (Agent::NICK_O_MAX - Agent::NICK_O_MIN)*(thecell->getGlucose()/(thecell->getGlucose() + Agent::NICK_O_CRITICAL_GLU)));

	// rescale from median to minimal cell volume
	oxygenUptake *= 2/3.;

	if( thecell->state == COMPARTMENT)
		oxygenUptake *= (double)thecell->growingTumorCellCount/(double)thecell->maxCellCount;

	return oxygenUptake;
	}
double GiveMeTheOxygenRate(VoronoiCell *thecell, float glu, float oxy)
	{
	//if( cell->getState() != COMPARTMENT)
	//if(thecell->getState() != NONACTIVE && thecell->getState() != ACTIVE)
	//	{return 0.;}

	double oxygenUptake = MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR
        * oxy/(oxy + Agent::NICK_O_CRITICAL_OXY)*(Agent::NICK_O_MAX - (Agent::NICK_O_MAX - Agent::NICK_O_MIN)*(glu/(glu + Agent::NICK_O_CRITICAL_GLU)));

	// rescale from median to minimal cell volume
	oxygenUptake *= 2/3.;

	if( thecell->getState() == COMPARTMENT)
		oxygenUptake *= (double)GetAgent(thecell)->growingTumorCellCount/(double)GetAgent(thecell)->maxCellCount;

	return oxygenUptake;
	}
	 
double GiveMeTheGlucoseRate(VoronoiCell *thecell)
	{
		//if(thecell->getState() != NONACTIVE && thecell->getState() != ACTIVE)
		//return 0.;

		double glucoseUptake = MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR/(thecell->glucose + Agent::NICK_G_CRITICAL_GLU)*(Agent::NICK_G_MAX - (Agent::NICK_G_MAX - Agent::NICK_G_MIN)*(thecell->oxygen/(thecell->oxygen + Agent::NICK_G_CRITICAL_OXY)));

	// rescale from median to minimal cell volume
	glucoseUptake *= 2/3.;
	
	if( thecell->getState() == COMPARTMENT)
		glucoseUptake *= (double)GetAgent(thecell)->growingTumorCellCount/(double)GetAgent(thecell)->maxCellCount;

	return glucoseUptake;
	}
double GiveMeTheGlucoseRate(Agent *thecell)
	{
	//	if(thecell->state != NONACTIVE && thecell->state != ACTIVE)
	//	return 0.;

		double glucoseUptake = MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR*thecell->getGlucose()/(thecell->getGlucose() + Agent::NICK_G_CRITICAL_GLU)*(Agent::NICK_G_MAX - (Agent::NICK_G_MAX - Agent::NICK_G_MIN)*(thecell->getOxygen()/(thecell->getOxygen() + Agent::NICK_G_CRITICAL_OXY)));

	// rescale from median to minimal cell volume
	glucoseUptake *= 2/3.;

	if( thecell->state == COMPARTMENT)
		glucoseUptake *= (double)thecell->growingTumorCellCount/(double)thecell->maxCellCount;

	return glucoseUptake;
	}
double GiveMeTheGlucoseRate( VoronoiCell *thecell, float glu, float oxy)
	{
	//	if(thecell->getState() != NONACTIVE && thecell->getState() != ACTIVE)
	//	return 0.;

		double glucoseUptake = MOL_PER_CELL_PER_SEC_TO_MILLIMOLAR_PER_HOUR*glu/(glu + Agent::NICK_G_CRITICAL_GLU)*(Agent::NICK_G_MAX - (Agent::NICK_G_MAX - Agent::NICK_G_MIN)*(oxy/(oxy + Agent::NICK_G_CRITICAL_OXY)));

	// rescale from median to minimal cell volume
	glucoseUptake *= 2/3.;

	if( thecell->getState() == COMPARTMENT)
		glucoseUptake *= (double)GetAgent(thecell)->growingTumorCellCount/(double)GetAgent(thecell)->maxCellCount;

	return glucoseUptake;
	}
#endif
/*****************************************************************************/

double GiveMeTheATP(Agent *thecell)
{
	double consumptionO = GiveMeTheOxygenRate ( thecell );
	double consumptionG = GiveMeTheGlucoseRate( thecell );
	double cATP = 2. * consumptionG + 34. * ( thecell->getGlucose() > 1e-2 ? consumptionO/6. : 0. );
	
	return cATP;
}

double GiveMeTheATP(VoronoiCell *thecell)
	{

	double consumptionO = thecell->oxygen  * GiveMeTheOxygenRate ( thecell );
	double consumptionG = thecell->glucose * GiveMeTheGlucoseRate( thecell );
	double cATP = 2. * consumptionG + 34. * ( thecell->glucose > 1e-2 ? consumptionO/6. : 0. );
	
	return cATP;
	}



