/*
 * AsciiIO.cpp
 *
 *  Created on: Dec 5, 2014
 *      Author: jagiella
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h> // strcpy
#include <math.h> // pow



#include "statistics.hpp"
#include "matrix.hpp"
#include "AsciiIO.h"

comparison_t create_comparison(){
	comparison_t c;
	c.size = c.dim = 0;
	c.x = c.m = c.s = 0;
	c.y = 0;
	return c;
}

int readFile( const char* filename, double **&cols, int ncol)
{
	//printf( "[ %i ]\n", sizeof(cols));
	if( !cols){
		cols = (double**)realloc( cols, sizeof(double*)*ncol);
		for( int i=0; i<ncol; i++)
			cols[i] = 0;
	}
	int nrow = 0;

	FILE *fp = fopen( filename, "r");
	char buffer[1024], *ptr;

	while( fgets( buffer, 1024, fp ))if(buffer[0] != '#'){
		ptr = buffer;
		nrow++;

		for( int i=0; i<ncol; i++){
			cols[i] = (double*)realloc( cols[i], sizeof(double)*nrow);
			cols[i][nrow-1]  = strtof( ptr, &ptr) ;
		}
	}
	fclose(fp);

	return nrow;
}
int readFileColumn( const char* filename, double *&col, int icol)
{
	int nrow = 0;

	FILE *fp = fopen( filename, "r");
	char buffer[1024], *ptr;

	while( fgets( buffer, 1024, fp ))
	if(buffer[0] != '#'){
		ptr = buffer;
		nrow++;

		for( int i=0; i<icol; i++)
			strtof( ptr, &ptr) ;

		col = (double*)realloc( col, sizeof(double)*nrow);
		col[nrow-1]  = strtof( ptr, &ptr) ;
	}
	fclose(fp);

	return nrow;
}

int readFileColumn( char** filename, int n, double **&col, int icol)
{
	int nrow = 0;
	double **old=col;
	col = (double**)realloc( col, sizeof(double*)*n);
	for( int i=0; i<n; i++){
		if( old != col) col[i] = 0;
		nrow = readFileColumn( filename[i], col[i], icol);
	}

	return nrow;
}

#define inrange( a, x, b) ( a<=x && x<=b )


