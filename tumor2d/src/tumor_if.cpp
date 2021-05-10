#include "Montecarlo.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <vector>

void tumor2d_interface(double InitialRadius, double InitialQuiescentFraction, double MaxCellDivisionRate, double DivisionDepth,
						double ECMThresholdQuiescence, double ECMProductionRate, double ECMDegradationRate,
						double EndTime, double OutputRate, double profileTime, int profileDepth, int rand_seed,
						std::vector<double> &gc_out, std::vector<double> &ecm_out, std::vector<double> &prolif_out)
{
	gc_out.reserve((int)(EndTime / OutputRate));  
	ecm_out.reserve(profileDepth);
	prolif_out.reserve(profileDepth);	
	double epsilon = montecarlo(InitialRadius, InitialQuiescentFraction, MaxCellDivisionRate, DivisionDepth,
							ECMThresholdQuiescence, ECMProductionRate, ECMDegradationRate,
							EndTime, OutputRate, profileTime, profileDepth, rand_seed,
							gc_out, ecm_out, prolif_out);
}

void tumor2d_default(std::vector<double> &gc_out, std::vector<double> &ecm_out, std::vector<double> &prolif_out)
{
	std:srand(time(0));
	int rand_seed = rand();
	double InitialRadius = 12.0;
	double InitialQuiescentFraction = 0.75;
	double MaxCellDivisionRate = 0.0417;
	double ECMThresholdQuiescence = 0.010;
	double ECMProductionRate = 0.005;
	double ECMDegradationRate = 0.0008;
	double EndTime     = 1000;
	double OutputRate = 24;
	double profileTime  = 800/24;
	double DivisionDepth = 100;
	int profileDepth = 1000;
	gc_out.reserve((int)(EndTime / OutputRate));  
	ecm_out.reserve(profileDepth);
	prolif_out.reserve(profileDepth);	
	double epsilon = montecarlo(InitialRadius, InitialQuiescentFraction, MaxCellDivisionRate, DivisionDepth,
							ECMThresholdQuiescence, ECMProductionRate, ECMDegradationRate,
							EndTime, OutputRate, profileTime, profileDepth, rand_seed,
							gc_out, ecm_out, prolif_out);
}
	
	
