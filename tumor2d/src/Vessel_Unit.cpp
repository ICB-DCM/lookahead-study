#include "Vessel_Unit.h"


	
void Vessel_Unit::Compute_radius()
	{

	radius = 0.5*(node1->radius + node2->radius);
	}

void Vessel_Unit::Compute_viscosity()
	{
	viscosity = 0.5*(node1->viscosity + node2->viscosity);
	}
void Vessel_Unit::Compute_flow()
     {
     flow = fabs(node1->pressure - node2->pressure);
     flow *= Pi*pow(radius,4);
     flow /= 8*viscosity*length;
     }

void Vessel_Unit::Compute_velocity()
     {
     velocity = flow;
     velocity /= Pi*radius*radius;
     }

void Vessel_Unit::Compute_shearstress()
     {
     shear = 0.5*length;
     shear *= radius;
     shear *= fabs(node1->pressure - node2->pressure);
     }
void Vessel_Unit::Compute_theta()
	{
	double l, X, Y, Z;
	X = (node2->x - node1->x);
	Y = (node2->y - node1->y);
	Z = (node2->z - node1->z);
	l = X*X + Y*Y;
	if( l <= 0 ){ l = Z*Z; }
	
	l = sqrt(l);
	
	if(X >= 0 && Y>= 0)
		{
		theta = acos(X/l) - Pi/2;

		}
	if(X < 0 && Y>= 0)
		{
		theta =  - Pi/2 + acos(X/l);
		}
	
	if(X >= 0 && Y < 0)
		{
		theta = -acos(X/l) - Pi/2;
		}

	if(X < 0 && Y < 0)
		{
		theta = -acos(X/l) - Pi/2;
		}	

	}
void Vessel_Unit::Compute_phi()
	{
	double X = (node2->x - node1->x);
	double Y = (node2->y - node1->y);
	double Z = (node2->z - node1->z);
	double l = X*X + Y*Y + Z*Z;
	l = sqrt(l);
	phi = -acos(Z/l);
	}

void Vessel_Unit::Compute_length()
	{
	length = (node1->x - node2->x)*(node1->x - node2->x);
	length += (node1->y - node2->y)*(node1->y - node2->y);
	length += (node1->z - node2->z)*(node1->z - node2->z);
	length = sqrt(length);
	}

