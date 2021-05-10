#ifndef __MATHEMATIX_H
#define __MATHEMATIX_H



// DEFINES

#ifndef TRUE
#define TRUE             1
#endif

#ifndef FALSE
#define FALSE            0
#endif
//#endif

// Constantes
#define PI 3.141592653589793238462643383279502

// Macros
#define absOfVector( x, y, z) \
        ( sqrt( x*x + y*y + z*z))

//#define min( a, b)	(a<b?a:b)
#define MIN( a, b)	(a<b?a:b)
//#define max( a, b)	(a>b?a:b)
#define MAX( a, b)	(a>b?a:b)

#if LOGARITHMIC_OUTPUT
    #define nrOfSteps( BeginningTime, EndTime, OutputRate) \
        ( (int) ( (log10(EndTime)-log10(BeginningTime))/OutputRate) + 3)
       // ( (int) ( (log10(EndTime)-log10(BeginningTime))/OutputRate) + 3 )
    #define indexOfTime( Time, BeginningTime, Step) \
        ( Time  > BeginningTime ? (int)(1. + (log10(Time)-log10(BeginningTime)) / Step + 0.5)   : 0)
         //( Time  > BeginningTime ? (int)((log10(Time)-log10(BeginningTime)) / Step + 1.5)   : 0)
   #define timeOfIndex( Index, BeginningTime, Step) \
        ( Index != 0  ? pow(10, log10(BeginningTime) + (double)(Index-1)*Step) : 0 )
#else
    #define nrOfSteps( BeginningTime, EndTime, Step) \
        ( (int) (EndTime/OutputRate) + 1 )
    #define indexOfTime( Time, BeginningTime, Step) \
        ( (int) (Time/Step) )
    #define timeOfIndex( Index, BeginningTime, Step) \
        ( (double) Index * Step )
#endif



// FUNCTION PROTOTYPES

// Computation
double myPow( double base, int exponent);

// Algebra
float ** newMatrix( int dimi, int dimj);
double ** newDoubleMatrix( int dimi, int dimj);
long double ** newLongDoubleMatrix( int dimi, int dimj);

void deleteMatrix( float ** matrix, int dimi);
void deleteDoubleMatrix( double ** matrix, int dimi);
void deleteLongDoubleMatrix( long double ** matrix, int dimi);

void solveLinearSystem( float **A, float *b, float *x, int dim);
void solveLinearSystem( double **A, double *b, double *x, int dim);

void solveLinearSystemB( float **A, float *b, float *x, int dim, float **B);
void solveLinearSystemB( double **A, double *b, double *x, int dim, double **B);
void solveLinearSystemB( long double **A, long double *b, long double *x, int dim, long double **B);

double solveDeterminant( double **A, int dim);
double getDeterministe( double **A, int dim);
void matrixVectorProduct( float **A, float *b, float *x, int dim);
void matrixVectorProduct( double **A, double *b, double *x, int dim);
void matrixTranspose( double **Ai, double **Ao, int n, int m);
void matrixInversion( double **Ai, double **Ao, int n);
void vectorSum( float *vectorA, float *vectorB, float *vectorA_B, int dim);
void vectorSum( double *vectorA, double *vectorB, double *vectorA_B, int dim);
void vectorDifference( float *vectorA, float *vectorB, float *vectorA_B, int dim);
void vectorDifference( double *vectorA, double *vectorB, double *vectorA_B, int dim);
void vectorCopy( float *vectorA, float *vectorB, int dim);
float dotProduct( float *vectorA, float *vectorB, int dim);
double dotProduct( double *vectorA, double *vectorB, int dim);
void vectorScale( float *vector, float scalar, float *vectorxscalar, int dim);
void vectorScale( double *vector, double scalar, double *vectorxscalar, int dim);
void matrixProduct( float **A, float **B, float **C, int n, int m, int p);
void matrixProduct( double **A, double **B, double **C, int n, int m, int p);
void matrixScale( double **Ai, double **Ao, double c, int n, int m);




// Numerics
void solveLinearSystem( float **A, float *b, float *x, int dim);
void solveLinearSystemTridiagonalMatrix( float **A, float *d, float *x, int dim);
double BiSection( double ( *function)( double), double a, double b, double error);
double BiSection2_1( double ( *function)( double, double, double), double c0, double c1, double xa, double xb, double error);

//double BiSectionFunction( double);
double NumericalIntegrationSimpson( double ( *function)( double), double xa, double xb, int n);
double NumericalIntegrationSimpson1_1( double ( *function)( double, double), double c, double a, double b, int n);
double NumericalIntegrationSimpson2_1( double ( *function)( double, double, double), double c0, double c1, double a, double b, int n);
double NumericalIntegrationQuadrature( double ( *function)( double), double xa, double xb, int n);

// Geometrie
void rotate( double* x, double* y, double* z, double x_r, double y_r, double z_r, double angle_r);
void rotateX( double& xi, double& yi, double& zi, double& xo, double& yo, double& zo, double& angle);
void rotateY( double& xi, double& yi, double& zi, double& xo, double& yo, double& zo, double& angle);
void rotateZ( double& xi, double& yi, double& zi, double& xo, double& yo, double& zo, double& angle);
void rotateXYZ( double xi, double yi, double zi, double& xo, double& yo, double& zo, double a, double b, double c);
void rotateZYX( double xi, double yi, double zi, double& xo, double& yo, double& zo, double a, double b, double c);

double dotProduct2D( double *vectorA, double *vectorB);
double dotProduct3D( double *vectorA, double *vectorB);
void crossProduct3D( double *vectorA, double *vectorB, double *vectorAxB);
void vectorDifference3D( float *vectorA, float *vectorB, float *vectorA_B);
void vectorDifference3D( double *vectorA, double *vectorB, double *vectorA_B);
void vectorScale3D( double *vector, double scalar, double *vectorxscalar);

// Stochastics
double myRand();
double myRandIE( double A, double B);
double myRandE( double B);
double myRandII( double A, double B);
double myRandI( double B);

double ProbabilityDesnityFunctionExponential( double x, double mean);
double CumulativeDistributionFunctionExponential( double x, double mean);

double ProbabilityDesnityFunctionErlang( double x, double mean, int k);
double CumulativeDistributionFunctionErlang( double x, double mean, int k);

double ProbabilityFunctionPoisson0( int k, double mean); // k>=0
double ProbabilityFunctionPoisson1( int k, double mean); // k>0
double CumulativeDistributionFunctionPoisson0( int k, double mean); // k>=0
double CumulativeDistributionFunctionPoisson1( int k, double mean); // k>0


void memoryAllocationError( char* const variableName);

#include "Mathematix.ipp"

#endif
