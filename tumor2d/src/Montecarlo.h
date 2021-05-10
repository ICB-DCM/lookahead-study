#ifndef __MONTECARLO_H
#define __MONTECARLO_H
#include <vector>
/****************************************************************************
 * Defines                                                                  *
 ****************************************************************************/
// UNIX/WINDOWS PORTABILITY
#if defined( __unix__) || defined(__MACH__)
   #include <sys/types.h>
   #include <sys/stat.h>
   #define MODUS     ,S_IRWXU|S_IRWXG|S_IRWXO
   #define SEPERATOR '/'
#elif defined(__WIN32__) || defined(_MS_DOS_)
    #include <dir.h>
    #define MODUS  
    #define SEPERATOR '\\'
#else
    #include <direct.h>  /* Visual C++ */
    #define MODUS  
    #define SEPERATOR '\\'
#endif

// OUTPUT SETTINGS
#define SINGLE_OUTPUT		FALSE
#define GLOBAL_OUTPUT		TRUE
#define RADIAL_OUTPUT		FALSE
#define PROB_OUTPUT			FALSE
#define SLICE_OUTPUT		TRUE
#define ON_THE_FLY_OUTPUT	TRUE

#define LOGARITHMIC_OUTPUT 0
#define CORRECTION_FACTOR_EXP 1.0 //1.389 // default: 1.0 

// CONSTANTS
#define CUBICMICROMETER_TO_MILLILITER (CUBICMICROMETER_TO_CUBICMETER * CUBICMETER_TO_LITER * LITER_TO_MILLILITER)    //1 000 000 000 000.
#define CUBICMETER_TO_LITER 1000.
#define CUBICMICROMETER_TO_CUBICMETER 0.000000000000000001
#define LITER_TO_MILLILITER 1000.

// VALUES
#define CELL_DIAMETER 15.
#define CELL_RADIUS (CELL_DIAMETER/2.)
#define CELL_HEIGTH_ON_SURFACE 3.5
#define CELL_MC_CONTACT_AREA 500.
#define CELL_DISTRIBUTION 1	// 0 - randomly among all microcarriers
				// 1 - uniformly among all microcarriers

// DETACHMENT SETTINGS
#define NECROTIC_DETACHMENT_ONLY 0


// ENVIRONMENT
#define REFILL_MEDIUM	FALSE

// GLOBAL VARIABLES
extern int Case;	

double montecarlo( double InitialRadius, double InitialQuiescentFraction, double MaxCellDivisionRate,double DivisionDepth,
				double ECMThresholdQuiescence, double ECMProductionRate, double ECMDegradationRate,
				double EndTime, double OutputRate, double profileTime, int profileDepth, int rand_seed,
				std::vector<double> &gc_out, std::vector<double> &ecm_out, std::vector<double> &prolif_out);

double DESubstrate( double substrate, double cells);
double rungeKutta( double (*equation)(double,double), double xStart, double xEnd, double yStart, double z, int steps);
double rungeKuttaStep( double (*equation)(double,double), double xStart, double xEnd, double yStart, double z);

#endif

