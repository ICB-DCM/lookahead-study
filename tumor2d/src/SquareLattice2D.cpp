#include <stdlib.h>
#include <stdio.h>

#include "SquareLattice2D.h"


/*****************************************************************************/

VoronoiDiagram* newSquareLattice( int countPoints[DIMENSIONS], int periodic[DIMENSIONS]){
//VoronoiDiagram3D* newSquareLattice2D( int countPointsX, int countPointsY, int periodicX, int periodicY){
/****************************************************************************
 * Definieren der Variablen                                                 *
 ****************************************************************************/

	VoronoiDiagram	*newVoronoiDiagram = VoronoiDiagram::newVoronoiDiagram();
	int		i, index, countNeighborCells;

	//int countPointsX = countPoints[0];
	//int countPointsY = countPoints[1];
	//int periodicX = periodic[0];
	//int periodicY = periodic[1];


/****************************************************************************
 * Dreieck- in Punktbeziehung wandeln                                       *
 ****************************************************************************/

	/* Anzahl der Punkte speichern */
	newVoronoiDiagram->countVoronoiCells = 1;            // Anzahl der Punkte speichern
	for( i=0; i<DIMENSIONS; i++)
		newVoronoiDiagram->countVoronoiCells *= countPoints[i];
	newVoronoiDiagram->maxVoronoiCells = newVoronoiDiagram->countVoronoiCells;



	newVoronoiDiagram->voronoiCells = ( VoronoiCell**) malloc ( newVoronoiDiagram->countVoronoiCells * sizeof( VoronoiCell*));
	for( i=0; i<newVoronoiDiagram->countVoronoiCells; i++)
		newVoronoiDiagram->voronoiCells[i] = ( VoronoiCell*) malloc ( sizeof( VoronoiCell));




/****************************************************************************
 *  alle Punkte in Datei schreiben                                          *
 ****************************************************************************/

	/* Gitter initialisieren*/

	
	
	/* Punkte mit Nachbarn in Datei sichern */

#if DIMENSIONS == 3 || DIMENSIONS == 2 || DIMENSIONS == 1

	int index_i[DIMENSIONS];
	for( i=0; i<DIMENSIONS; i++)
		index_i[i] = 0;


	
	for( index = 0; index < newVoronoiDiagram->countVoronoiCells; index++ ){

		// set number of voronoi cell
		newVoronoiDiagram->voronoiCells[index]->index = index;

		newVoronoiDiagram->voronoiCells[index]->agent = NULL;

		// set coordinates

		for( i=0; i<DIMENSIONS; i++){
			newVoronoiDiagram->voronoiCells[index]->position[i] = (double) index_i[i];

		}

			
		// set neighbors of voronoi cell
		newVoronoiDiagram->voronoiCells[index]->neighborCells = ( VoronoiCell**) calloc ( 2 * DIMENSIONS, sizeof( VoronoiCell*)); // allocate memory for list of neighbors

		// number of neighbours
		countNeighborCells = 0;
		newVoronoiDiagram->voronoiCells[index]->countNeighborCells = 2 * DIMENSIONS;

		// set neighbors for each dimension
		int ni, 
		    nindex[2],
		    factor;
		for( i=0; i<DIMENSIONS; i++){

			// determine indexes of 2 neighbors in actual dimension
			nindex[0] = nindex[1] = 0;
			factor = 1;
			for( ni=0; ni<DIMENSIONS; ni++){
				if( i!=ni){
					nindex[0] += index_i[ni]*factor;
					nindex[1] += index_i[ni]*factor;
				}else{
					nindex[0] += ((index_i[ni]+countPoints[ni]-1)%countPoints[ni] )*factor;
					nindex[1] += ((index_i[ni]                +1)%countPoints[ni] )*factor;
				}
				factor *= countPoints[ni];
			}

			// set neighbors for actual dimension
			if( index_i[i]>0 || periodic[i]){

				newVoronoiDiagram->voronoiCells[ index]->neighborCells[ countNeighborCells++] = newVoronoiDiagram->voronoiCells[ nindex[0]];
			}else
				newVoronoiDiagram->voronoiCells[ index]->countNeighborCells--;

			if( index_i[i]<countPoints[i]-1 || periodic[i]){

				newVoronoiDiagram->voronoiCells[ index]->neighborCells[ countNeighborCells++] = newVoronoiDiagram->voronoiCells[ nindex[1]];
			}else
				newVoronoiDiagram->voronoiCells[ index]->countNeighborCells--;

		}
		newVoronoiDiagram->voronoiCells[ index]->countFreeNeighborCells = newVoronoiDiagram->voronoiCells[ index]->countNeighborCells;

		// extended neighbors
		newVoronoiDiagram->voronoiCells[index]->countExtendedNeighborCells = 0;
		newVoronoiDiagram->voronoiCells[index]->countFreeExtendedNeighborCells = 0;
		newVoronoiDiagram->voronoiCells[index]->extendedNeighborhood = NULL;

		// increase indexes
		i = 0;
		do{
			index_i[i] = (index_i[i]+1)%countPoints[i];
			i++;
		}while( i<DIMENSIONS && index_i[i-1]==0);
			
	}
#else
	for( j = 0; j<countPointsY; j++ ){
		for( i = 0; i<countPointsX; i++ ){
			
			index = i + j*countPointsX;
			
			// set number of voronoi cell
			newVoronoiDiagram->voronoiCells[index]->index = index;

			// set coordinates
			newVoronoiDiagram->voronoiCells[index]->position[0] = (double) i;
			newVoronoiDiagram->voronoiCells[index]->position[1] = (double) j;
			//newVoronoiDiagram->voronoiCells[index]->position[2] = strtod( ptr, &ptr );
		
			// set neighbors of voronoi cell
			newVoronoiDiagram->voronoiCells[index]->neighborCells = ( VoronoiCell3D**) calloc ( 4, sizeof( VoronoiCell3D*)); // allocate memory for list of neighbors
			
			// number of neighbours
			countNeighborCells = 0;
			newVoronoiDiagram->voronoiCells[index]->countNeighborCells = 4;

			if( i>0 ){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[(i-1) + j*countPointsX];
			}else if( periodicX){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[countPointsX-1 + j*countPointsX];
			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;

			if( i<countPointsX-1 ){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[(i+1) + j*countPointsX];
			}else if( periodicX){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[0 + j*countPointsX];
			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;

			if( j>0 ){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[i + (j-1)*countPointsX];
			}else if( periodicY){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[i + (countPointsY-1)*countPointsX];
			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;

			if( j<countPointsY-1 ){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[i + (j+1)*countPointsX];
			}else if( periodicY){
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[i + (0)*countPointsX];
			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;

		}
	}
#endif
	
	/*for( i = 0; i<lastpoint; i++ ){
		printf( "point %i: ", newVoronoiDiagram->voronoiCells[i].index);
		for( j = 0; j < newVoronoiDiagram->voronoiCells[i].countNeighborCells; j++)
			printf( " %i", newVoronoiDiagram->voronoiCells[i].neighborCells[j]->index);
		printf( "\n");
	}*/
	

	return newVoronoiDiagram;
}
/*****************************************************************************/
