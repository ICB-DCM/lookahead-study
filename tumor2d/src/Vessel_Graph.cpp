#include "Vessel_Graph.h"

#include <stdio.h>

Vessel_Graph::Vessel_Graph()
		{
		x = -1;
		y = -1;
		z = -1;
		index = -1;
		pressure = 0;
		viscosity = 0;
		radius = 0;
		outflow = 0;
		circulation = true;
		InDomain = true;
		CountBranches = 0;

		d_max = (xmaxDOMAIN - xminDOMAIN)*xmaxDOMAIN;
		d_max += (ymaxDOMAIN - yminDOMAIN)*ymaxDOMAIN;
		d_max += (zmaxDOMAIN - zminDOMAIN)*zmaxDOMAIN;

		d_min = (xmaxDOMAIN - xminDOMAIN)*xminDOMAIN;
		d_min += (ymaxDOMAIN - yminDOMAIN)*yminDOMAIN;
		d_min += (zmaxDOMAIN - zminDOMAIN)*zminDOMAIN;
		}


void Vessel_Graph::PressureComputation()
	{
	//Pourquoi utiliser la fonction InDomain? Pour ne pas pr�calculer, ou utliser mx!=NONE?
	if(InDomain)
		{
		double Sum1 = 0;
		double Sum2 = 0;
		int j;

		for(j = 0; j < CountBranches; j++)
			{
			double lenght = (branch[j]->x - x)*(branch[j]->x - x);
			lenght += (branch[j]->y - y)*(branch[j]->y - y);
			lenght += (branch[j]->z - z)*(branch[j]->z - z);
			lenght = sqrt(lenght);
			double meanviscosity = 0.5*(branch[j]->viscosity + viscosity);

			Sum1 += pow((branch[j]->radius + radius)/2.0,4)*(branch[j]->pressure)/(meanviscosity*lenght);
			Sum2 += pow((branch[j]->radius + radius)/2.0,4)/(meanviscosity*lenght);
				
			}
		pressure = Sum1/Sum2;
		}
	else
		{
		//if(z < zminDOMAIN + BoundaryLayer|| z > zmaxDOMAIN - BoundaryLayer){pressure = WallPRESSURE;}
		/*
		if(x < xminDOMAIN + BoundaryLayer){pressure = InletPRESSURE;}
		if(y < yminDOMAIN + BoundaryLayer){pressure = InletPRESSURE;}
		if(z < zminDOMAIN + BoundaryLayer){pressure = InletPRESSURE;}

		if(x > xmaxDOMAIN - BoundaryLayer){pressure = OutletPRESSURE;}
		if(y > ymaxDOMAIN - BoundaryLayer){pressure = OutletPRESSURE;}
		if(z > zmaxDOMAIN - BoundaryLayer){pressure = OutletPRESSURE;}
		*/

		d = (xmaxDOMAIN - xminDOMAIN)*x;
		d += (ymaxDOMAIN - yminDOMAIN)*y;
		d += (zmaxDOMAIN - zminDOMAIN)*z;
		
		pressure = ( OutletPRESSURE - InletPRESSURE ) * ( d - d_min) / ( d_max - d_min ) + InletPRESSURE;
		}

	}

int Vessel_Graph::FindBranchNumber(Vessel_Graph *b)
	{
	int j;
	for(j = 0; j < CountBranches; j++)
		{
		if(branch[j] == b){return j;}
		}
	printf("Erreur in FindBranchNumber: I didn't find the branch.\n");
	return -7;
	}

bool Vessel_Graph::DefineInDomain()
	{
	if(x < xminDOMAIN + BoundaryLayer){return false;}
	if(y < yminDOMAIN + BoundaryLayer){return false;}
	if(z < zminDOMAIN + BoundaryLayer){return false;}
	if(x > xmaxDOMAIN - BoundaryLayer){return false;}
	if(y > ymaxDOMAIN - BoundaryLayer){return false;}
	if(z > zmaxDOMAIN - BoundaryLayer){return false;}
	
	
	return true;
	}



/*
//Quicksort
int partitionner(Vessel_Graph** tableau, int p, int r)
	{
        double pivot = tableau[p]->x;
	int i = p-1, j = r+1;
        Vessel_Graph* temp;
        while(1)
		{
                do
                        j--;
                	while(tableau[j]->x > pivot);
                do
                        i++;
                	while(tableau[i]->x < pivot);
                	if(i<j)
				{
                        	temp = tableau[i];
                        	tableau[i] = tableau[j];
                        	tableau[j] = temp;
                		}
                	else return j;
        	}
        return j;
	}

void quickSort(Vessel_Graph** tab, int p, int r)
	{
        int q;
        if(p<r) {
                q = partitionner(tab, p, r);
                quickSort(tab, p, q);
                quickSort(tab, q+1, r);
        	}
	}
//end of quicksort
*/
