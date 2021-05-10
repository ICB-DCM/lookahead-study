/*
 * statistics.hpp
 *
 *  Created on: Oct 7, 2014
 *      Author: jagiella
 */

#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

template <class T> T unifrnd()                    { return (T)rand()      /(T)(RAND_MAX); }
template <class T> T unifrnd( unsigned int *seed) { return (T)rand_r(seed)/(T)(RAND_MAX); }
void    normrnd( double *rnd, int n, unsigned int *p_seed);
double  normrnd( unsigned int *p_seed);
double  normrnd();
void mvnrnd( double *x, double *mean, double **s2, int n, unsigned int *p_seed);
double mvnpdf( const double *x, double *mean, double **inv_s2, double det_s2, int n);
double kdepdf( const double *x, double **samples, double **inv_s2, double det_s2, int n, int nsamples);



typedef struct{
	bool   H;
	double pValue;
	double KSstatistic;
} KSTestResult;

KSTestResult KolmogorovSmirnoffTest( int n1, double *x1, int n2, double *x2);

double logLikelihood( double *v1, double *v2, double *s21, double *s22, int n);
double logLikelihood( double v1, double v2, double s21);




enum VectorOrientation {rowwise = 1, columnwise = 2};
double *mean( double **x, int n, int d, char orientation = 1);
void    mean( double **x, int d, int n, double *mean_x, char orientation = 1);
double  mean( double  *x, int n);
double *std2(   double **x, double *m, int n, int d, char orientation = 1);
void    stdmean( double *x, double *dx, double d, int n, double &SEM_x, double &mean_x);
double median( double  *x, int n);
void   median( double **x, int d, int n, double *&median_x);
double mad(    double  *x, int n);
void findmin( double *x, int n, double &o_min, int &o_idx);
void findmax( double *x, int n, double &max, int &idx);
double  min( double  *x, int n);
void    min( double *x, int n, double &xmin, int &imin);
template <class T>
T min( T x1, T x2){
	if( x1 < x2)
		return x1;
	else
		return x2;
}
double  max( double  *x, int n);
double quantile( double *x, int n, double f);

void cov( double **A, int n, int m, double *mean, double **cov);


#endif /* STATISTICS_HPP_ */
