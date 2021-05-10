#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "statistics.hpp"

int cmpfunc (const void * a, const void * b){
   return ( *(double*)a < *(double*)b ? -1 : 1 );
}

KSTestResult KolmogorovSmirnoffTest( int n1, double *x1, int n2, double *x2){
	int i1=0, i2=0;
	double y1 = 0, y2 = 0;

	// sort x1 & x2
	qsort( x1, n1, sizeof(double), cmpfunc);
	qsort( x2, n2, sizeof(double), cmpfunc);

	KSTestResult result;

	// Statistics
	result.KSstatistic = 0;
	do{
		// which is next?
		if( i1 < n1 && (i2==n2 || x1[i1] < x2[i2])){
			y1 += 1./n1;
			i1++;
		}else{
			y2 += 1./n2;

			i2++;
		}

		result.KSstatistic = fmax( fabs(y1-y2), result.KSstatistic);

	}while(i1<n1 || i2<n2);

	// p-Value
	result.pValue = 0;
	int    n      =  n1 * n2 /(n1 + n2);
	double lambda =  fmax((sqrt(n) + 0.12 + 0.11/sqrt(n)) * result.KSstatistic , 0.);
	for( int j=1; j<=101; j++)
		result.pValue  +=  2 * pow(-1, j-1)*exp(-2*lambda*lambda*j*j);
	result.pValue  =  fmin( fmax(result.pValue, 0.), 1.);

	// Significance Test
	double alpha = 0.05;
	result.H = ( alpha >= result.pValue ? true : false);

	return result;
}


double logLikelihood( double v1, double v2, double s21){
	if(s21==0)
		return 0;
	else
		return - 0.5 * pow( v1 - v2, 2) / s21;// - 0.5 * log(2*M_PI*s21);
}

double logLikelihood( double *v1, double *v2, double *s21, double *s22, int n)
{
	double logL = 0;

	/*if( s22){
		for( int i=0; i<n; i++) if(s21[i])
			logL += - pow( v1[i] - v2[i], 2) / s21[i] - pow(s21[i] - s22[i], 2) / s21[i];
	}else*/
	if( s21){
		for( int i=0; i<n; i++)
			if( s21[i]>0){
				double last = logL;
				logL += logLikelihood( v1[i], v2[i], s21[i]);
				if( isinf(logL) || isnan(logL) || isinf(-logL)){
					fprintf( stderr, "%i: v1=%e, v2=%e, s21=%e -> %e + %e = %e -> %e\n", i, v1[i], v2[i], s21[i], - 0.5 * log(2*M_PI*s21[i]), - 0.5 * pow( v1[i] - v2[i], 2) / s21[i], - 0.5 * log(2*M_PI*s21[i]) - 0.5 * pow( v1[i] - v2[i], 2) / s21[i], last);
										exit(0);
				}
			}
	}else
		for( int i=0; i<n; i++)
			logL += - pow( (v1[i] - v2[i]), 2);

	return logL;
}


double *mean( double **x, int n1, int n2, char orientation){
	double *m;

	switch( orientation){
	case rowwise: // row wise
		m = (double*) malloc( n2 * sizeof(double));
		for( int j=0; j<n2; j++){
			m[j] = 0;
			for( int i=0; i<n1; i++)
				m[j] += x[i][j];
			m[j] /= (double) n1;
		}
		break;

	case columnwise: // column wise
		m = (double*) malloc( n1 * sizeof(double));
		for( int j=0; j<n1; j++){
			m[j] = 0;
			for( int i=0; i<n2; i++)
				m[j] += x[j][i];
			m[j] /= (double) n2;
		}
		break;

	}

	return m;
}

void stdmean( double *x, double *dx, double d, int n, double &SEM_x, double &mean_x){
	mean_x = 0;
	double s2_x  = 0;
	SEM_x = 0;

	if( dx)
		for( int i=0; i<n; i++){
			mean_x +=    (x[i] + dx[i]*d);
			s2_x   += pow(x[i] + dx[i]*d, 2.);
		}
	else
		for( int i=0; i<n; i++){
			mean_x += x[i];
			s2_x   += x[i]*x[i];
		}

	mean_x = mean_x / (double) n;
	s2_x   = s2_x   / (double) n - mean_x*mean_x;
	SEM_x  = sqrt(s2_x / (double) n);
}

void mean( double **x, int d, int n, double *mean_x, char orientation){

	if( orientation == 1)
		for( int id=0; id<d; id++){
			mean_x[id] = 0;
			for( int i=0; i<n; i++)
				mean_x[id] += x[id][i];
			mean_x[id] = mean_x[id] / (double) n;
		}
	else
    	for( int i=0; i<n; i++){
	    	mean_x[i] = 0;
	    	for( int id=0; id<d; id++)
	    		mean_x[i] += x[id][i];
	    	mean_x[i] = mean_x[i] / (double) d;
	    }

}
double mean( double *x, int n){
	double mean_x = 0;
	for( int i=0; i<n; i++)
		mean_x += x[i];
	return mean_x / (double) n;
}

void cov( double **A, int n, int m, double *mean, double **cov){
	// covariance
	for( int j=0;  j<m;  j++)
		for( int jj=0; jj<m; jj++){
			cov[j][jj] = 0;
			for( int i=0;  i<n;  i++)
				cov[j][jj] += (A[i][j] - mean[j])*(A[i][jj] - mean[jj]);
			cov[j][jj] /= (double)(n-1);
		}
}



double min( double *x, int n){
	double min = x[0];
	for( int i=1; i<n; i++)
		min = fmin( min, x[i]);
	return min;
}
void min( double *x, int n, double &xmin, int &imin){
	imin = 0;
	xmin = x[0];
	for( int i=1; i<n; i++)
	if( x[i] < xmin){
		xmin = x[i];
		imin =   i ;
	}
}

double max( double *x, int n){
	double max = x[0];
	for( int i=1; i<n; i++)
		max = fmax( max, x[i]);
	return max;
}

double *std2( double **x, double *m, int n, int d, char orientation){
	double *s2 = (double*) malloc( d * sizeof(double));

	if( orientation == 1)
		for( int j=0; j<d; j++){
			s2[j] = 0;
			for( int i=0; i<n; i++){
				s2[j] += pow( x[i][j] - m[j], 2);

			}
			s2[j] /= (double) n;
		}
	else
		for( int i=0; i<n; i++){
			s2[i] = 0;
			for( int j=0; j<d; j++){
				s2[i] += pow( x[i][j] - m[i], 2);

			}
			s2[i] /= (double) d;
		}


	return s2;
}

int compare_dbl (const void * a, const void * b)
{
   return ( *(double*)a - *(double*)b );
}
void median( double **x, int d, int n, double *&median_x){

    for( int id=0; id<d; id++){
    	double x_sorted[n];
    	for( int i=0; i<n; i++)	x_sorted[i] = x[id][i];
    	qsort( x_sorted, n, sizeof(double), compare_dbl);
    	median_x[id] = x_sorted[n/2];
    }
}

double median( double *x, int n){

	double x_sorted[n];
	for( int i=0; i<n; i++)	x_sorted[i] = x[i];
	qsort( x_sorted, n, sizeof(double), compare_dbl);

    return x_sorted[n/2];
}

double mad( double *x, int n){
	double r[n];
	double median_x = median( x, n);
	for( int i=0; i<n; i++)
		r[i] = fabs( x[i] - median_x);
	return median( r, n);
}

void findmin( double *x, int n, double &o_min, int &o_idx){
	o_min = x[0];
	o_idx = 0;
	for( int i=1; i<n; i++)
		if( o_min > x[i]){
			o_min = x[i];
			o_idx = i;
		}
}

void findmax( double *x, int n, double &max, int &idx){
	max = x[0];
	idx = 0;
	for( int i=1; i<n; i++)
		if( max < x[i]){
			max = x[i];
			idx = i;
		}
}

void normrnd( double *rnd, int n, unsigned int *p_seed)
{
	// Box-Muller

	for( int i=0; i<n; i+=2){
		double u1=unifrnd<double>( p_seed), u2=unifrnd<double>( p_seed);
		rnd[i]   = sqrt(-2*log(u1))*cos(2*M_PI*u2);
		if(i+1<n)
		rnd[i+1] = sqrt(-2*log(u1))*sin(2*M_PI*u2);
	}
}
double normrnd( unsigned int *p_seed)
{
	// Box-Muller
	double u1=unifrnd<double>( p_seed), u2=unifrnd<double>( p_seed);
	return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}
double normrnd()
{
	// Box-Muller
	double u1=unifrnd<double>(), u2=unifrnd<double>();
	return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

double mvnpdf_helper( double *a, double **B, double *c, int n){

	double rtn=0, tmpi;
	for( int i=0; i<n; i++){
		tmpi = 0;
		for( int j=0; j<n; j++)
			tmpi += a[j]*B[j][i];
		rtn += tmpi*c[i];
	}

	return rtn;
}

double mvnpdf( const double *x, double *mean, double **inv_s2, double det_s2, int n)
{
	double x0[n];
	for( int i=0; i<n; i++)
		x0[i] = x[i] - mean[i];
	return 1. / sqrt( pow(2*M_PI, n) * det_s2) * exp(-0.5*mvnpdf_helper( x0, inv_s2, x0,n));
}

double kdepdf( const double *x, double **samples, double **inv_s2, double det_s2, int n, int nsamples)
{
	double pdf = 0;
	for( int i=0; i<nsamples; i++)
		pdf += mvnpdf( x, samples[i], inv_s2, det_s2, n);
	return pdf / (double)nsamples;
}

void mvnrnd( double *x, double *mean, double **s2, int n, unsigned int *p_seed)
{
	// sample from MVnorm dist around chosen point
	double z[n];
	normrnd( z, n, p_seed);


	for( int i=0; i<n; i++){
		x[i] = mean[i];
		double tmp = 0;
		for( int j=0; j<n; j++){
			x[i] += s2[i][j] * z[j];
			tmp += s2[i][j] * z[j];
		}

	}



}
