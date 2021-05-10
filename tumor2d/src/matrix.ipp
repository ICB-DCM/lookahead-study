#ifndef MATRIX_IPP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// Definition of a singular matrix
#define SINGULAR 1e-16

template <class T>
T **allocMatrix( int n, int m)
{
	T **A = (T**) malloc( n*sizeof(T*));
	for( int i=0; i<n; i++){
		A[i] = (T*) malloc( m*sizeof(T));
		for( int j=0; j<m; j++)
			A[i][j] = 0.;
	}
	return A;
}

template <class T>
void freeMatrix( T** &A, int n)
{

	for( int i=0; i<n; i++)
		free(A[i]);
	free(A);
}


template <class T>
void vectorMatrixProduct( T* &v, T** &A, T* &u, int m, int n)
{
	for( int j=0; j<n; j++){
		u[j] = 0;
		for( int i=0; i<m; i++)
			u[j] += v[i] * A[i][j];
	}
}
template <class T>
void matrixVectorProduct( T** &A, T* &v, T* &u, int m, int n)
{
	for( int i=0; i<m; i++){
		u[i] = 0;
		for( int j=0; j<n; j++)
			u[i] += A[i][j] * v[j];
	}
}
template <class T>
T dotProduct( T* &v1, T* &v2, int n)
{
	T dotProd = 0;
	for( int i=0; i<n; i++){
		dotProd += v1[i] * v2[i];
	}
	return dotProd;
}





template <class T>
void decomposeGaussJordan( T** &A, T** &LR, int n) {
	// Gauss Elimination
	for( int i=0; i<n; n++){
	       // Bestimmen von R
	       for( int j=i; j<n; j++){
	           for( int k=0; k<i-1; k++){
	               A[i][j] -= A[i][k] * A[k][j];
	           }
	       }
	       // Bestimmen von L
	       for( int j=i+1; j<n; j++){
	           for( int k=0; k<i-1; k++){
	               A[j][i] -= A[j][k] * A[k][i];
	           }
	           A[j][i] /= A[i][i];
	       }
	}
}


template <class T>
void decompCholeskyOwn( T** &A, T** &L, int n) {

	for( int i=0; i<n; i++)
		for( int j=0; j<=i; j++){
			double Summe = A[i][j];
			for( int k=0; k<=j-1; k++)
				Summe = Summe - L[i][k] * L[j][k];
			if( i > j){
				L[i][j] = Summe / L[j][j];   // Untere Dreiecksmatrix
				L[j][i] = 0.;                // Obere Dreiecksmatrix
				//L[j][i] = Summe / A[j][j]; // Obere Dreiecksmatrix
			}
			else if( Summe > 0)          // Diagonalelement
				L[i][i] = sqrt( Summe);       // ... ist immer groesser Null
			else
	            fprintf( stderr, "Matrix not positive definite\n");   // ERROR
		}
}

template <class T>
void invertGaussJordan( T** &a, T** &ainv, int n) {
	int i, j;                    // Zeile, Spalte
	int s;                       // Elimininationsschritt
	int pzeile;                  // Pivotzeile
	int fehler = 0;              // Fehlerflag
	double f;                      // Multiplikationsfaktor
	const double Epsilon = 0.01;   // Genauigkeit
	double Maximum;                // Zeilenpivotisierung
	FILE *fout = stdout;
	int pivot = 1;

	// erg�nze die Matrix a um eine Einheitsmatrix (rechts anh�ngen)
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][n + j] = 0.0;
			if (i == j)
				a[i][n + j] = 1.0;
		}
	}*/
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i][j] = 0.0;
			if (i == j)
				ainv[i][j] = 1.0;
		}
	}
#if DEBUG
	MatOut (stdout, a, n, 2*n);
#endif

	// die einzelnen Eliminationsschritte
	s = 0;
	do {
		// Pivotisierung vermeidet unn�tigen Abbruch bei einer Null in der Diagnonalen und
		// erh�ht die Rechengenauigkeit
		Maximum = fabs(a[s][s]);
		if (pivot) {
			pzeile = s;
			for (i = s + 1; i < n; i++)
				if (fabs(a[i][s]) > Maximum) {
					Maximum = fabs(a[i][s]);
					pzeile = i;
				}
		}
		fehler = (Maximum < Epsilon);

		if (fehler)
			break;           // nicht l�sbar

		if (pivot) {
			if (pzeile != s)  // falls erforderlich, Zeilen tauschen
					{
				double h;
				/*for (j = s; j < 2 * n; j++) {
					h = a[s][j];
					a[s][j] = a[pzeile][j];
					a[pzeile][j] = h;
				}*/
				for (j = s; j < n; j++) {
					h = a[s][j];
					a[s][j] = a[pzeile][j];
					a[pzeile][j] = h;
				}
				for (j = 0; j < n; j++) {
					h = ainv[s][j];
					ainv[s][j] = ainv[pzeile][j];
					ainv[pzeile][j] = h;
				}
			}
		}

		// Eliminationszeile durch Pivot-Koeffizienten f = a[s][s] dividieren
		f = a[s][s];
		/*for (j = s; j < 2 * n; j++)
			a[s][j] = a[s][j] / f;*/
		for (j = s; j < n; j++)
			a[s][j] = a[s][j] / f;
		for (j = 0; j < n; j++)
			ainv[s][j] = ainv[s][j] / f;

		// Elimination --> Nullen in Spalte s oberhalb und unterhalb der Diagonalen
		// durch Addition der mit f multiplizierten Zeile s zur jeweiligen Zeile i
		for (i = 0; i < n; i++) {
			if (i != s) {
				f = -a[i][s];                 // Multiplikationsfaktor
				/*for (j = s; j < 2 * n; j++)    // die einzelnen Spalten
					a[i][j] += f * a[s][j];*/       // Addition der Zeilen i, s
				for (j = s; j < n; j++)    // die einzelnen Spalten
					a[i][j] += f * a[s][j];       // Addition der Zeilen i, s
				for (j = 0; j < n; j++)    // die einzelnen Spalten
					ainv[i][j] += f * ainv[s][j];       // Addition der Zeilen i, s
			}
		}
#if DEBUG
		fprintf(stdout, "Nach %1i-tem Eliminationschritt:\n", s+1);
		MatOut (stdout, a, n, 2*n);
#endif
		s++;
	} while (s < n);

	if (fehler) {
		fprintf(fout, "Inverse: Matrix ist singulaer\n");
		//return 0;
	}
	// Die angeh�ngte Einheitsmatrix Matrix hat sich jetzt in die inverse Matrix umgewandelt
	// Umkopieren auf die Zielmatrix
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i][j] = a[i][n + j];
		}
	}*/
}


template <class T>
void invertGaussJordan( T* &a, T* &ainv, int n) {
	int i, j;                    // Zeile, Spalte
	int s;                       // Elimininationsschritt
	int pzeile;                  // Pivotzeile
	int fehler = 0;              // Fehlerflag
	double f;                      // Multiplikationsfaktor
	const double Epsilon = 0.01;   // Genauigkeit
	double Maximum;                // Zeilenpivotisierung
	FILE *fout = stdout;
	int pivot = 1;

	// erg�nze die Matrix a um eine Einheitsmatrix (rechts anh�ngen)
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][n + j] = 0.0;
			if (i == j)
				a[i][n + j] = 1.0;
		}
	}*/
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i*n+j] = 0.0;
			if (i == j)
				ainv[i*n+j] = 1.0;
		}
	}
#if DEBUG
	MatOut (stdout, a, n, 2*n);
#endif

	// die einzelnen Eliminationsschritte
	s = 0;
	do {
		// Pivotisierung vermeidet unn�tigen Abbruch bei einer Null in der Diagnonalen und
		// erh�ht die Rechengenauigkeit
		Maximum = fabs(a[s*n+s]);
		if (pivot) {
			pzeile = s;
			for (i = s + 1; i < n; i++)
				if (fabs(a[i*n+s]) > Maximum) {
					Maximum = fabs(a[i*n+s]);
					pzeile = i;
				}
		}
		fehler = (Maximum < Epsilon);

		if (fehler)
			break;           // nicht l�sbar

		if (pivot) {
			if (pzeile != s)  // falls erforderlich, Zeilen tauschen
					{
				double h;
				/*for (j = s; j < 2 * n; j++) {
					h = a[s][j];
					a[s][j] = a[pzeile][j];
					a[pzeile][j] = h;
				}*/
				for (j = s; j < n; j++) {
					h = a[s*n+j];
					a[s*n+j] = a[pzeile*n+j];
					a[pzeile*n+j] = h;
				}
				for (j = 0; j < n; j++) {
					h = ainv[s*n+j];
					ainv[s*n+j] = ainv[pzeile*n+j];
					ainv[pzeile*n+j] = h;
				}
			}
		}

		// Eliminationszeile durch Pivot-Koeffizienten f = a[s][s] dividieren
		f = a[s*n+s];
		/*for (j = s; j < 2 * n; j++)
			a[s][j] = a[s][j] / f;*/
		for (j = s; j < n; j++)
			a[s*n+j] = a[s*n+j] / f;
		for (j = 0; j < n; j++)
			ainv[s*n+j] = ainv[s*n+j] / f;

		// Elimination --> Nullen in Spalte s oberhalb und unterhalb der Diagonalen
		// durch Addition der mit f multiplizierten Zeile s zur jeweiligen Zeile i
		for (i = 0; i < n; i++) {
			if (i != s) {
				f = -a[i*n+s];                 // Multiplikationsfaktor
				/*for (j = s; j < 2 * n; j++)    // die einzelnen Spalten
					a[i][j] += f * a[s][j];*/       // Addition der Zeilen i, s
				for (j = s; j < n; j++)    // die einzelnen Spalten
					a[i*n+j] += f * a[s*n+j];       // Addition der Zeilen i, s
				for (j = 0; j < n; j++)    // die einzelnen Spalten
					ainv[i*n+j] += f * ainv[s*n+j];       // Addition der Zeilen i, s
			}
		}
#if DEBUG
		fprintf(stdout, "Nach %1i-tem Eliminationschritt:\n", s+1);
		MatOut (stdout, a, n, 2*n);
#endif
		s++;
	} while (s < n);

	if (fehler) {
		fprintf(fout, "Inverse: Matrix ist singulaer\n");
		//return 0;
	}
	// Die angeh�ngte Einheitsmatrix Matrix hat sich jetzt in die inverse Matrix umgewandelt
	// Umkopieren auf die Zielmatrix
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i][j] = a[i][n + j];
		}
	}*/
}





template <class T>
void solveLinearSystemB( T **A, T *b, T *x, int dim, T **B)
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
				T temp;
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

template <class T>
void solveLinearSystem( T **A, T *b, T *x, int dim)
{
	T **B = allocMatrix<T>( dim, dim);
	if( B == NULL){
		exit( 0);
	}
	solveLinearSystemB<T>( A, b, x, dim, B);
	freeMatrix( B, dim);
}

template <class T>
void normalizeVector( T *x, int m){

	T sum = 0;
	for( int j=0; j<m; j++)
		sum += x[j];
	for( int j=0; j<m; j++)
		x[j] /= sum;
}

template <class T>
T norm( T *x, int dim, int order = 0)
{
	if( order == 0){
		// MAX NORM
		T max = 0;
		for( int i=0; i<dim; i++)
			max = fmax( fabs(x[i]), max);
		return max;
	}else{
		// MAX NORM
		T sum = 0;
		for( int i=0; i<dim; i++)
			sum += pow( x[i], order);
		return pow( sum, 1./order);
	}
}

template <class T>
T norm_diff( T *x1, T *x2, int dim, int order = 0)
{
	if( order == 0){
		// MAX NORM
		T max = 0;
		for( int i=0; i<dim; i++)
			max = fmax( fabs(x1[i]-x2[i]), max);
		return max;
	}else{
		// MAX NORM
		T sum = 0;
		for( int i=0; i<dim; i++)
			sum += pow( x1[i]-x2[i], order);
		return pow( sum, 1./order);
	}
}

template <class T>
bool inbound( T *x, T *lb, T *ub, int d)
{
	for( int i=0; i<d; i++)
		if( x[i] < lb[i] || x[i] > ub[i])
			return false;
	return true;
}
template <class T>
bool inbound( T *x, T lb, T ub, int d)
{
	for( int i=0; i<d; i++)
		if( x[i] < lb || x[i] > ub)
			return false;
	return true;
}
template <class T>
bool inbound( T *x, T lb, T ub, int d, char mode)
{
	if( mode == 0){ // AND (all in bound)
		for( int i=0; i<d; i++)
			if( x[i] < lb || x[i] > ub)
				return false;
		return true;
	}

	if( mode == 1){ // OR (at least one in bound)
		for( int i=0; i<d; i++)
			if( x[i] >= lb && x[i] <= ub)
				return true;
		return false;
	}
}


#define MATRIX_IPP
#endif
