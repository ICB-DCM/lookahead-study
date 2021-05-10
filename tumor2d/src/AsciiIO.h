/*
 * AsciiIO.h
 *
 *  Created on: Dec 5, 2014
 *      Author: jagiella
 */

#ifndef SRC_TUMOR3D_ASCIIIO_H_
#define SRC_TUMOR3D_ASCIIIO_H_



typedef struct {
	int 	dim, size;
	double *x; // x
	double **y;// y's
	double *m; // mean of y
	double *s; // standard deviation of y
} comparison_t;

comparison_t create_comparison();
int readFileColumn( const char* filename, double *&col, int icol);
int readFileColumn( char** filename, int n, double **&col, int icol);

enum   compare_mode{ mean_vs_mean, mean_vs_single};
//double compare( comparison_t d1, comparison_t d2, char mode );

#endif /* SRC_TUMOR3D_ASCIIIO_H_ */
