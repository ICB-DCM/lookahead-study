#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#ifdef OMP
  #include <omp.h>
/*#else
  #warning "no openmp support"
#endif*/

#include "Mathematix.hpp"


/****************************************************************************
 * ALGEBRA                                                                  *
 ****************************************************************************/

double myPow( double base, int exponent)
{
	if( exponent == 0)
		return 1.;
	else
		return base * myPow( base, exponent-1 );
}
/*****************************************************************************/


double getDeterministe( double **A, int dim)
{
	if(dim==1)
		return A[0][0];
	if(dim==2)
		return A[0][0]*A[1][1] - A[0][1]*A[1][0];
	if(dim==3)
		return A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[1][2]*A[2][1];

	double **B = newDoubleMatrix( dim-1, dim-1);
	if( B == NULL) exit( 0);

	double det = 0.;
	for( int i=0; i<dim; i++)
		for( int j=0; j<dim; j++){
			for( int ii=0; ii<dim; ii++)
				for( int jj=0; jj<dim; jj++){
					B[(ii<=i?ii:ii-1)][(jj<=j?jj:jj-1)] = A[i][j];
				}
			det += pow(-1, i+j) * A[i][j] * getDeterministe( B, dim-1);
		}
	deleteDoubleMatrix( B, dim-1);

	return det;
}
/*****************************************************************************/


float getDeterministe( float **A, int dim)
{
	if(dim==3)
		return A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[1][2]*A[2][1];

	float **B = newMatrix( dim-1, dim-1);
	if( B == NULL) exit( 0);

	float det = 0.;
	for( int i=0; i<dim; i++)
		for( int j=0; j<dim; j++){
			for( int ii=0; ii<dim; ii++)
				for( int jj=0; jj<dim; jj++){
					B[(ii<=i?ii:ii-1)][(jj<=j?jj:jj-1)] = A[i][j];
				}
			det += pow(-1, i+j) * A[i][j] * getDeterministe( B, dim-1);
		}
	deleteMatrix( B, dim-1);

	return det;
}
/*****************************************************************************/


double solveDeterminant( double **A, int dim)
{
	double ret = 0.;
	double term;

	// positive diagonals
	for( int i=0; i<dim; i++){
		term = 0.;
		for( int j=0; j<dim; j++){
			term *= A[(j+i)%dim][j%dim];
		}
		ret += term;
	}

	// negative diagonals
	for( int i=0; i<dim; i++){
		term = 0.;
		for( int j=0; j<dim; j++){
			term *= A[(dim+i-j)%dim][j%dim];
		}
		ret -= term;
	}

	return ret;
}
/*****************************************************************************/


void solveLinearSystem( float **A, float *b, float *x, int dim)
{
	float **B = newMatrix( dim, dim);
	if( B == NULL){
		exit( 0);	
	}
	solveLinearSystemB( A, b, x, dim, B);
	deleteMatrix( B, dim);
}
/*****************************************************************************/


void solveLinearSystem( double **A, double *b, double *x, int dim)
{
	double **B = newDoubleMatrix( dim, dim);
	if( B == NULL){
		exit( 0);
	}
	solveLinearSystemB( A, b, x, dim, B);
	deleteDoubleMatrix( B, dim);
}
/*****************************************************************************/


// PIVOT GAUSSE //
void solveLinearSystemB( float **A, float *b, float *x, int dim, float **B)
{
	int i, // index of equation
	    j; // index of column
	int k;
	//double B[dim][dim];
	//float **B = newMatrix( dim, dim);

	// copy matrix
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			B[i][j] = A[i][j];
	//		printf("%lf  ", B[i][j]);
		}
	//	printf("| %lf\n", b[i]);
	}

	// solving the linear system

	// forward reduction
	for( k=0; k<dim; k++){
		//printf("forward reduction: line %i/%i!\n", k+1, dim);
		if(B[k][k]==0){
			// find better line
//			printf("looking for better line!\n");
			for( i=k+1; i<dim && B[i][k]==0; i++);

			// resort lines
			if(i<dim && B[i][k]!=0){
//				printf("resort!\n");
				double temp;
				for( j=k; j<dim; j++){
					temp = B[i][j];
					B[i][j] = B[k][j];
					B[k][j] = temp;
				}
				temp = b[i];
				b[i] = b[k];
				b[k] = temp;
			}
		}
		if(B[k][k]!=0){
			// normalize first row
			for( j=k+1; j<dim; j++)
				B[k][j]=B[k][j]/B[k][k];
			b[k]=b[k]/B[k][k];
			B[k][k]=1.;

			// reduce following rows
			for( i=k+1; i<dim; i++){
				for( j=k+1; j<dim; j++)
					B[i][j]=B[i][j]-B[i][k]*B[k][j];
				b[i]=b[i]+b[k]*-B[i][k];
				B[i][k]=0;
			}
		}
	}

	/*printf("----------------------------\n");
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			printf("%lf  ", B[i][j]);
		}
		printf("| %lf\n", b[i]);
	}//*/

	// backward reduction
	for( k=dim-1; k>=0; k--){
		if( B[k][k]!=0)
		for( i=0; i<k; i++){
			b[i]=b[i]+b[k]*-B[i][k];
			B[i][k]=0;
		}
	}
	
	/*printf("----------------------------\n");
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			printf("%lf  ", B[i][j]);
		}
		printf("| %lf\n", b[i]);
	}//*/

	// copy solution 
	for( i=0; i<dim; i++)
		x[i] = b[i];
		//x[i] = B[i][i];
}
/*****************************************************************************/


void solveLinearSystemB( double **A, double *b, double *x, int dim, double **B)
{
	int i, // index of equation
	    j; // index of column
	int k;
	//double B[dim][dim];
	//float **B = newMatrix( dim, dim);

	// copy matrix
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			B[i][j] = A[i][j];
	//		printf("%lf  ", B[i][j]);
		}
	//	printf("| %lf\n", b[i]);
	}

	// solving the linear system

	// forward reduction
	for( k=0; k<dim; k++){
		//printf("forward reduction: line %i/%i!\n", k+1, dim);
		if(B[k][k]==0){
			// find better line
//			printf("looking for better line!\n");
			for( i=k+1; i<dim && B[i][k]==0; i++);

			// resort lines
			if(i<dim && B[i][k]!=0){
//				printf("resort!\n");
				double temp;
				for( j=k; j<dim; j++){
					temp = B[i][j];
					B[i][j] = B[k][j];
					B[k][j] = temp;
				}
				temp = b[i];
				b[i] = b[k];
				b[k] = temp;
			}
		}
		if(B[k][k]!=0){
			// normalize first row
			for( j=k+1; j<dim; j++)
				B[k][j]=B[k][j]/B[k][k];
			b[k]=b[k]/B[k][k];
			B[k][k]=1.;

			// reduce following rows
			for( i=k+1; i<dim; i++){
				for( j=k+1; j<dim; j++)
					B[i][j]=B[i][j]-B[i][k]*B[k][j];
				b[i]=b[i]+b[k]*-B[i][k];
				B[i][k]=0;
			}
		}
	}

	/*printf("----------------------------\n");
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			printf("%lf  ", B[i][j]);
		}
		printf("| %lf\n", b[i]);
	}//*/

	// backward reduction
	for( k=dim-1; k>=0; k--){
		if( B[k][k]!=0)
		for( i=0; i<k; i++){
			b[i]=b[i]+b[k]*-B[i][k];
			B[i][k]=0;
		}
	}

	/*printf("----------------------------\n");
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			printf("%lf  ", B[i][j]);
		}
		printf("| %lf\n", b[i]);
	}//*/

	// copy solution
	for( i=0; i<dim; i++)
		x[i] = b[i];
		//x[i] = B[i][i];
}
/*****************************************************************************/


void solveLinearSystemB( long double **A, long double *b, long double *x, int dim, long double **B)
{
	int i, // index of equation
	    j; // index of column
	int k;
	//double B[dim][dim];
	//float **B = newMatrix( dim, dim);

	// copy matrix
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			B[i][j] = A[i][j];
	//		printf("%lf  ", B[i][j]);
		}
	//	printf("| %lf\n", b[i]);
	}

	// solving the linear system

	// forward reduction
	for( k=0; k<dim; k++){
		//printf("forward reduction: line %i/%i!\n", k+1, dim);
		if(B[k][k]==0){
			// find better line
//			printf("looking for better line!\n");
			for( i=k+1; i<dim && B[i][k]==0; i++);

			// resort lines
			if(i<dim && B[i][k]!=0){
//				printf("resort!\n");
				long double temp;
				for( j=k; j<dim; j++){
					temp = B[i][j];
					B[i][j] = B[k][j];
					B[k][j] = temp;
				}
				temp = b[i];
				b[i] = b[k];
				b[k] = temp;
			}
		}
		if(B[k][k]!=0){
			// normalize first row
			for( j=k+1; j<dim; j++)
				B[k][j]=B[k][j]/B[k][k];
			b[k]=b[k]/B[k][k];
			B[k][k]=1.;

			// reduce following rows
			for( i=k+1; i<dim; i++){
				for( j=k+1; j<dim; j++)
					B[i][j]=B[i][j]-B[i][k]*B[k][j];
				b[i]=b[i]+b[k]*-B[i][k];
				B[i][k]=0;
			}
		}
	}

	/*printf("----------------------------\n");
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			printf("%lf  ", B[i][j]);
		}
		printf("| %lf\n", b[i]);
	}//*/

	// backward reduction
	for( k=dim-1; k>=0; k--){
		if( B[k][k]!=0)
		for( i=0; i<k; i++){
			b[i]=b[i]+b[k]*-B[i][k];
			B[i][k]=0;
		}
	}

	/*printf("----------------------------\n");
	for( i=0; i<dim; i++){
		for( j=0; j<dim; j++){
			printf("%lf  ", B[i][j]);
		}
		printf("| %lf\n", b[i]);
	}//*/

	// copy solution
	for( i=0; i<dim; i++)
		x[i] = b[i];
		//x[i] = B[i][i];
}
/*****************************************************************************/


void TridiagonalSolve(float *a, float *b, float *c, float *d, float *x, int n){
	int i;
 
	/* Modify the coefficients. */
	c[0] /= b[0];				/* Division by zero risk. */
	d[0] /= b[0];				/* Division by zero would imply a singular matrix. */
	for(i = 1; i < n; i++){
		double id = (b[i] - c[i-1] * a[i]);	/* Division by zero risk. */
		c[i] /= id;				/* Last value calculated is redundant. */
		d[i] = (d[i] - d[i-1] * a[i])/id;
	}
 
	/* Now back substitute. */
	x[n - 1] = d[n - 1];
	for(i = n - 2; i >= 0; i--)
		x[i] = d[i] - c[i] * x[i + 1];
}
/*****************************************************************************/


void solveLinearSystemTridiagonalMatrix( float **A, float *d, float *x, int dim){
 	float a[dim];
	float b[dim];
	float c[dim];
	
	for(int m=0; m<dim; m++){
		b[m] = A[m][m];
		if(m>0)
			a[m-1] = A[m][m-1];
		if(m<dim-1)
			c[m+1] = A[m][m+1];			
	}

	TridiagonalSolve( a, b, c, d, x, dim);
}
/*****************************************************************************/


double ** newDoubleMatrix( int dimi, int dimj)
{
	int i;
	double **A;

	A = (double**) calloc( sizeof(double*), dimi );
	if( A == NULL) memoryAllocationError( (char*)"A");

	for( i=0; i<dimi; i++){
		A[i] = (double*) calloc( sizeof(double), dimj);
		if( A[i] == NULL) memoryAllocationError( (char*)"A[i]");
	}

	return A;
}
/*****************************************************************************/


long double ** newLongDoubleMatrix( int dimi, int dimj)
{
	int i;
	long double **A;

	A = (long double**) calloc( sizeof(long double*), dimi );
	if( A == NULL) memoryAllocationError( (char*)"A");

	for( i=0; i<dimi; i++){
		A[i] = (long double*) calloc( sizeof(long double), dimj);
		if( A[i] == NULL) memoryAllocationError( (char*)"A[i]");
	}

	return A;
}
/*****************************************************************************/


float ** newMatrix( int dimi, int dimj)
{
	int i;
	float **A;

	A = (float**) calloc( sizeof(float*), dimi );
	if( A == NULL) memoryAllocationError( (char*)"A");

	for( i=0; i<dimi; i++){
		A[i] = (float*) calloc( sizeof(float), dimj);
		if( A[i] == NULL) memoryAllocationError( (char*)"A[i]");
	}

	return A;
}
/*****************************************************************************/


void deleteMatrix( float ** matrix, int dimi)
{
	int i;
	
	for( i=0; i<dimi; i++)
		free( matrix[i]);
	free( matrix);
}
/*****************************************************************************/


void deleteDoubleMatrix( double ** matrix, int dimi)
{
	int i;

	for( i=0; i<dimi; i++)
		free( matrix[i]);
	free( matrix);
}
/*****************************************************************************/


void deleteLongDoubleMatrix( long double ** matrix, int dimi)
{
	int i;

	for( i=0; i<dimi; i++)
		free( matrix[i]);
	free( matrix);
}
/*****************************************************************************/




/****************************************************************************
 * GEOMETRIE                                                                *
 ****************************************************************************/

/*****************************************************************************/


void rotate( double* x, double* y, double* z, double x_r, double y_r, double z_r, double angle_r){
  double value_r;
  double a[3][3];
  double p[3] = {*x, *y, *z};

  // Normalize rotation vector
  value_r = sqrt(x_r*x_r + y_r*y_r + z_r*z_r);
  x_r = x_r/value_r;
  y_r = y_r/value_r;
  z_r = z_r/value_r;

  // Rotation matrix
  a[0][0] = cos(angle_r) + x_r*x_r*(1-cos(angle_r));
  a[1][0] = x_r*y_r*(1-cos(angle_r)) - z_r*sin(angle_r);
  a[2][0] = x_r*z_r*(1-cos(angle_r)) + y_r*sin(angle_r);

  a[0][1] = x_r*y_r*(1-cos(angle_r)) + z_r*sin(angle_r);
  a[1][1] = cos(angle_r) + y_r*y_r*(1-cos(angle_r));
  a[2][1] = y_r*z_r*(1-cos(angle_r)) - x_r*sin(angle_r);

  a[0][2] = x_r*z_r*(1-cos(angle_r)) - y_r*sin(angle_r);
  a[1][2] = y_r*z_r*(1-cos(angle_r)) + x_r*sin(angle_r);
  a[2][2] = cos(angle_r) + z_r*z_r*(1-cos(angle_r));

  // matrix-vector-product
  *x = p[0]*a[0][0] + p[1]*a[1][0] + p[2]*a[2][0];
  *y = p[0]*a[0][1] + p[1]*a[1][1] + p[2]*a[2][1];
  *z = p[0]*a[0][2] + p[1]*a[1][2] + p[2]*a[2][2];
}
/*****************************************************************************/


void rotateX( double& xi, double& yi, double& zi, double& xo, double& yo, double& zo, double& angle){
	//temp
	double temp = yi*sin(angle) + zi*cos(angle);
	
	yo = yi*cos(angle) - zi*sin(angle);
	zo = temp; //yi*sin(angle) + zi*cos(angle)
	
	xo = xi;
}
/*****************************************************************************/


void rotateY( double& xi, double& yi, double& zi, double& xo, double& yo, double& zo, double& angle){
	// temp
	double temp = zi*cos(angle) - xi*sin(angle);
	
	xo = zi*sin(angle) + xi*cos(angle);
	zo = temp;// zi*cos(angle) - xi*sin(angle)

	yo = yi;
}
/*****************************************************************************/


void rotateZ( double& xi, double& yi, double& zi, double& xo, double& yo, double& zo, double& angle){
	//temp
	double temp = xi*sin(angle) + yi*cos(angle);

	xo = xi*cos(angle) - yi*sin(angle);
	yo = temp; //xi*sin(angle) + yi*cos(angle)
	
	zo = zi;
}
/*****************************************************************************/


void rotateXYZ( double xi, double yi, double zi, double& xo, double& yo, double& zo, double a, double b, double c){
	xo = xi * cos(b)*cos(c)
	   + yi * cos(b)*sin(c)
	   - zi * sin(b);
	yo = xi * (sin(a)*sin(b)*cos(c) - cos(a)*sin(c))
	   + yi * (sin(a)*sin(b)*sin(c) + cos(a)*cos(c))
	   + zi * (sin(a)*cos(b));
	zo = xi * (cos(a)*sin(b)*cos(c) + sin(a)*sin(c))
	   + yi * (cos(a)*sin(b)*sin(c) - sin(a)*cos(c))
	   + zi * (cos(a)*cos(b));
}
/*****************************************************************************/


void rotateZYX( double xi, double yi, double zi, double& xo, double& yo, double& zo, double a, double b, double c){
	xo = xi *  cos(b)*cos(c)
	   + yi * (sin(a)*sin(b)*cos(c) + cos(a)*sin(c))
	   + zi *(-cos(a)*sin(b)*cos(c) + sin(a)*sin(c));

	yo = xi * -cos(b)*sin(c)
	   + yi *(-sin(a)*sin(b)*sin(c) + cos(a)*cos(c))
	   + zi * (cos(a)*sin(b)*sin(c) + sin(a)*cos(c));

	zo = xi * (sin(b))
	   + yi *(-sin(a)*cos(b))
	   + zi * (cos(a)*cos(b));
}
/*****************************************************************************/


double dotProduct2D( double *vectorA, double *vectorB)
{
	return vectorA[0]*vectorB[0] + vectorA[1]*vectorB[1];
}
/*****************************************************************************/


double dotProduct3D( double *vectorA, double *vectorB)
{
	return vectorA[0]*vectorB[0] + vectorA[1]*vectorB[1] + vectorA[2]*vectorB[2];
}
/*****************************************************************************/


float dotProduct( float *vectorA, float *vectorB, int dim)
{
	float ret = 0.;

	for( int i=0; i<dim; i++){
		ret += vectorA[i]*vectorB[i];
#if DEBUG > 0
		if( isnan(ret)){
			fprintf( stderr, "nan occures in dotProduct: %lf * %lf = %lf\n", vectorA[i], vectorB[i], vectorA[i]*vectorB[i]);
			//exit( 0);
			return 1.;
		}
#endif
	}
		
	return ret;
}
/*****************************************************************************/


double dotProduct( double *vectorA, double *vectorB, int dim)
{
	double ret = 0.;

#pragma omp parallel for reduction(+:ret)
	for( int i=0; i<dim; i++){
		ret += vectorA[i]*vectorB[i];
#if DEBUG > 0
		if( isnan(ret)){
			fprintf( stderr, "nan occures in dotProduct: %lf * %lf = %lf\n", vectorA[i], vectorB[i], vectorA[i]*vectorB[i]);
			//exit( 0);
			return 1.;
		}
#endif
	}

	return ret;
}
/*****************************************************************************/


void crossProduct3D( double *vectorA, double *vectorB, double *vectorAxB)
{
	//(a2b3 − a3b2, a3b1 − a1b3, a1b2 − a2b1)
	vectorAxB[0] = vectorA[1]*vectorB[2] - vectorA[2]*vectorB[1];
	vectorAxB[1] = vectorA[2]*vectorB[0] - vectorA[0]*vectorB[2];
	vectorAxB[2] = vectorA[0]*vectorB[1] - vectorA[1]*vectorB[0];
}
/*****************************************************************************/


void vectorDifference3D( double *vectorA, double *vectorB, double *vectorA_B)
{
	vectorA_B[0] = vectorA[0]-vectorB[0];
	vectorA_B[1] = vectorA[1]-vectorB[1];
	vectorA_B[2] = vectorA[2]-vectorB[2];
}
/*****************************************************************************/


void vectorDifference3D( float *vectorA, float *vectorB, float *vectorA_B)
{
	vectorA_B[0] = vectorA[0]-vectorB[0];
	vectorA_B[1] = vectorA[1]-vectorB[1];
	vectorA_B[2] = vectorA[2]-vectorB[2];
}
/*****************************************************************************/


void vectorDifference( float *vectorA, float *vectorB, float *vectorA_B, int dim)
{
	for( int d=0; d<dim; d++)
		vectorA_B[d] = vectorA[d] - vectorB[d];
}
/*****************************************************************************/


void vectorDifference( double *vectorA, double *vectorB, double *vectorA_B, int dim)
{
#pragma omp parallel for
	for( int d=0; d<dim; d++)
		vectorA_B[d] = vectorA[d] - vectorB[d];
}
/*****************************************************************************/


void vectorSum( float *vectorA, float *vectorB, float *vectorA_B, int dim)
{
	for( int d=0; d<dim; d++)
		vectorA_B[d] = vectorA[d] + vectorB[d];
}
/*****************************************************************************/


void vectorSum( double *vectorA, double *vectorB, double *vectorA_B, int dim)
{
#pragma omp parallel for
	for( int d=0; d<dim; d++)
		vectorA_B[d] = vectorA[d] + vectorB[d];
}
/*****************************************************************************/


void vectorCopy( float *vectorA, float *vectorB, int dim)
{
	for( int d=0; d<dim; d++)
		vectorB[d] = vectorA[d];
}
/*****************************************************************************/


void vectorScale3D( double *vector, double scalar, double *vectorxscalar)
{
	vectorxscalar[0] = vector[0]*scalar;
	vectorxscalar[1] = vector[1]*scalar;
	vectorxscalar[2] = vector[2]*scalar;
}
/*****************************************************************************/


void vectorScale( double *vector, double scalar, double *vectorxscalar, int dim)
{
#pragma omp parallel for
	for( int i=0; i<dim; i++){
		vectorxscalar[i] = vector[i]*scalar;

	}
}
/*****************************************************************************/


void vectorScale( float *vector, float scalar, float *vectorxscalar, int dim)
{
	for( int i=0; i<dim; i++){
		vectorxscalar[i] = vector[i]*scalar;
	}
}
/*****************************************************************************/


void matrixVectorProduct( float **A, float *b, float *x, int dim)
{
	for( int m=0; m<dim; m++){
		x[m] = 0.;
		for( int n=0; n<dim; n++){
			x[m] += A[m][n] * b[n];
		}	
	}
}
/*****************************************************************************/


void matrixVectorProduct( double **A, double *b, double *x, int dim)
{
//#pragma omp parallel for
	for( int m=0; m<dim; m++){
		x[m] = 0.;
		for( int n=0; n<dim; n++){
			x[m] += A[m][n] * b[n];
		}
	}
}
/*****************************************************************************/


void matrixProduct( float **A, float **B, float **C, int n, int m, int p)
{
	// A*B = C
	// (n x m) * (m x p) = (n x p)

	for( int i=0; i<n; i++)
	for( int j=0; j<p; j++)
	for( int k=0; k<m; k++){
		C[i][j] += A[n][m] * B[m][p];	
	}
}
/*****************************************************************************/


void matrixProduct( double **A, double **B, double **C, int n, int m, int p)
{
	// A*B = C
	// (n x m) * (m x p) = (n x p)

//#pragma omp parallel for
	for( int i=0; i<n; i++)
		for( int k=0; k<p; k++)
			for( int j=0; j<m; j++)
				if( j==0){
					C[i][k] = A[i][j] * B[j][k];
				}
				else{
					C[i][k] += A[i][j] * B[j][k];
				}
}
/*****************************************************************************/


void matrixScale( double **Ai, double **Ao, double c, int n, int m)
{
	// A*B = C
	// (n x m) * (m x p) = (n x p)

//#pragma omp parallel for
	for( int i=0; i<n; i++)
		for( int j=0; j<m; j++)
			Ao[i][j] = Ai[i][j] * c;

}
/*****************************************************************************/


void matrixTranspose( double **Ai, double **Ao, int n, int m)
{
//#pragma omp parallel for
	for( int i=0; i<n; i++)
		for( int j=0; j<m; j++)
			Ao[j][i] = Ai[i][j];

}
/*****************************************************************************/


void matrixInversion( double **Ai, double **Ao, int n)
{
	//double **A = newDoubleMatrix( n, n);//[n][n];
	double A[n][n];

	// Init Ao & A
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++){
			if(i==j)
				Ao[i][j] = 1;
			else
				Ao[i][j] = 0;
			A[i][j] = Ai[i][j];
		}

	// Gauss-Jordan-Algorithm
	//forward
	//actual diagonal element
	for( int d=0; d<n; d++){
		// rescale first line
		for( int j=d+1; j<n; j++)
			A[d][j] /= A[d][d];
		for( int j=0; j<n; j++)
			Ao[d][j] /= A[d][d];
		A[d][d] = 1.;

		//
		for( int i=d+1; i<n; i++)
			if(A[i][d]!=0){
				// adapt folowing lines
				//double c = (A[i][d]/A[d][d]);
				for( int j=d+1; j<n; j++)
					//B[i][j]=B[i][j]-B[i][k]*B[k][j];
					A[i][j] = A[i][j] -A[d][j] *A[i][d];
				for( int j=0; j<n; j++)
					Ao[i][j]= Ao[i][j]-Ao[d][j]*A[i][d];
				A[i][d] = 0.; // NOT NECESSARY!
			}
	}


	// backward
	for( int d=n-1; d>0; d--){
		for( int i=0; i<d; i++)
			if(A[i][d]!=0){
				for( int j=0; j<n; j++)
					Ao[i][j]= Ao[i][j]-Ao[d][j]*A[i][d];
				//fprintf( stderr , "test\n");
				A[i][d] = 0.; // NOT NECESSARY!
			}
			//else
				//fprintf( stderr , "test2 \n");
	}

}
/*****************************************************************************/





/****************************************************************************
 * NUMERIC                                                                  *
 ****************************************************************************/

double BiSection( double ( *function)( double), double a, double b, double error)
{
	if( b - a < error)
		return (a + b) / 2.;

	if( function(a) * function(b) < 0.){
		return BiSection( function, a, (a + b) / 2., error);
	}else{
		return BiSection( function, (a + b) / 2., b, error);
	}
}
double BiSection2_1( double ( *function)( double, double, double), double c0, double c1, double xa, double xb, double error)
{
	//fprintf( stderr, "a:%lf (f(a):%lf), b:%lf (f(b):%lf)\n", xa, function(c0, c1, xa), xb, function(c0, c1, xb));
	if( xb - xa < error)
		return (xa + xb) / 2.;

	double xab = (xa + xb) / 2.;
	if( function(c0, c1, xa) * function(c0, c1, xab) < 0.){
		return BiSection2_1( function, c0, c1, xa, xab, error);
	}else{
		return BiSection2_1( function, c0, c1, xab, xb, error);
	}
}

double NumericalIntegrationQuadrature( double ( *function)( double), double a, double b, int n)
{
	double ans = (function( a) + function( b))/2.;
	
	int k;
	
	for( k=1; k<n; k++)
		ans += function( a + k*(b-a)/n);
		
	return ans * (b-a)/n;
}		


/*double NumericalIntegrationStepSimpson( double ( *function)( double), double a, double b)
{
	return (b-a)/6. * ( function(a) + 4.*function( (a+b)/2) + function(b));
}*/

double NumericalIntegrationSimpson( double ( *function)( double), double a, double b, int n)
{
	int k;
	double ans = 0.;
	double h = (b-a)/(double)n;
	double x0 = a;
	double x1; //= a + h;
	
	ans += 0.5 * (function(a) + function(b));

	for( k=1; k<n; k++){
		x1 = x0 + h;
		ans += ( 2.*function( (x0+x1)/2.) + function(x1));
		x0 += h;
	}
	ans += ( 2.*function( (x0+b)/2.));

	return h/3. * ans;
}

double NumericalIntegrationSimpson2_1( double ( *function)( double, double, double), double c0, double c1, double a, double b, int n)
{
	int k;
	double ans = 0.;
	double h = (b-a)/(double)n;
	double x0 = a;
	double x1; //= a + h;
	
	ans += 0.5 * (function(c0, c1, a) + function(c0, c1, b));

	for( k=1; k<n; k++){
		x1 = x0 + h;
		ans += ( 2.*function( c0, c1, (x0+x1)/2.) + function(c0, c1, x1));
		x0 += h;
	}
	ans += ( 2.*function( c0, c1, (x0+b)/2.));

	return h/3. * ans;
}


/****************************************************************************
 * STOCHASTIC                                                               *
 ****************************************************************************/

double myRand()
{
	return (double)rand()/(double)RAND_MAX;	
}

double myRandIE( double A, double B)
{
	return A + (B-A)*(double)rand()/(double)(RAND_MAX + 1.);	
}

double myRandE( double B)
{
	return B*(double)rand()/(double)(RAND_MAX + 1.);	
}

double myRandII( double A, double B)
{
	return A + (B-A)*(double)rand()/(double)(RAND_MAX);	
}

double myRandI( double B)
{
	return B*(double)rand()/(double)(RAND_MAX);	
}

double ProbabilityDesnityFunctionExponential( double x, double mean)
{
	return exp( -x/mean)/mean;	
}

double CumulativeDistributionFunctionExponential( double x, double mean)
{
	return 1. - exp( -x/mean);	
}

double ProbabilityDesnityFunctionErlang( double x, double mean, int k)
{
	double lamda = k/mean;
	double ans = 1.;
	//double ans = exp( -lamda * x) * lamda;
	
	int i;
	for( i=1; i<k; i++)
		ans *= lamda * x / i;
	
	return ans * exp( -lamda * x) * lamda;
}

double CumulativeDistributionFunctionErlang( double x, double mean, int k)
{
	double ans = 0.;
	double term;
	int n, i;
	//double lamba = k/mean;
	double lamda_x = x * k/mean;
	//double exp_lamda_x = exp( -lamda_x);
	double exp_lamda_x_times_lamda_x = lamda_x * exp( -lamda_x);
	
	ans = 0.;
	
	for( n=1; n<k; n++){
		term = exp_lamda_x_times_lamda_x;
		
		for( i=2; i<=n; i++)
			term *= lamda_x / i;
		fprintf( stderr, "term: %lf\n", term);
		ans += term;
	}
	
	return 1. - (ans);	
}

double ProbabilityFunctionPoisson1( int k, double mean)
{
	double prob = exp(-mean);
	double i;
	for( i=1; i<=k-1; i++ )
		prob *= mean/i;

	return prob;
}

double ProbabilityFunctionPoisson0( int k, double mean)
{
	double prob = exp(-mean);
	double i;
	for( i=1; i<=k; i++ )
		prob *= mean/i;

	return prob;
}

double CumulativeDistributionFunctionPoisson1( int k, double mean)
{
	int i;
	double prob = 0.;

	for( i=1; i<=k; i++)
		prob += ProbabilityFunctionPoisson1( i, mean);
	return prob;
}

double CumulativeDistributionFunctionPoisson0( int k, double mean)
{
	int i;
	double prob = 0.;

	for( i=0; i<=k; i++)
		prob += ProbabilityFunctionPoisson0( i, mean);
	return prob;
}
/*****************************************************************************/


void memoryAllocationError( char* const variableName)
{
	fprintf( stderr, "Can't allocate memory for %s\n", variableName);
	exit( 0);
}

