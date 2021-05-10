#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "HexagonalLattice2D.h"


/*****************************************************************************/


VoronoiDiagram* newHexagonalLattice2D( int countPointsX, int countPointsY, int periodicX, int periodicY){
/****************************************************************************
 * Definieren der Variablen                                                 *
 ****************************************************************************/

	VoronoiDiagram	*newVoronoiDiagram = VoronoiDiagram::newVoronoiDiagram();
	int		i, j, index, countNeighborCells;

	// make even
	if( countPointsX%2 == 1) countPointsX++; 
	if( countPointsY%2 == 1) countPointsY++; 


/****************************************************************************
 * Dreieck- in Punktbeziehung wandeln                                       *
 ****************************************************************************/

	/* Anzahl der Punkte speichern */
	newVoronoiDiagram->countVoronoiCells = countPointsX * countPointsY;            // Anzahl der Punkte speichern
	newVoronoiDiagram->maxVoronoiCells   = countPointsX * countPointsY;            // Anzahl der Punkte speichern

	newVoronoiDiagram->voronoiCells = ( VoronoiCell**) malloc ( newVoronoiDiagram->countVoronoiCells * sizeof( VoronoiCell*));
	for( i=0; i<newVoronoiDiagram->countVoronoiCells; i++)
		newVoronoiDiagram->voronoiCells[i] = ( VoronoiCell*) malloc ( sizeof( VoronoiCell));




/****************************************************************************
 *  alle Punkte in Datei schreiben                                          *
 ****************************************************************************/

	/* Gitter initialisieren*/
	/*double a = 2. / pow( 3., 0.25);
	double stepX = a;
	double stepY = sqrt( 3.)/2. * a;*/

	double stepX = 2. / pow( 3., 0.25);
	double stepY = sqrt( 3.) / pow( 3., 0.25);
	int ii, jj;
	
	/* Punkte mit Nachbarn in Datei sichern */
	for( j = 0; j<countPointsY; j++ ){
		for( i = 0; i<countPointsX; i++ ){
			
			index = j*countPointsX + i;

			
			// set number of voronoi cell
			newVoronoiDiagram->voronoiCells[index]->index = index;

			// set coordinates
			newVoronoiDiagram->voronoiCells[index]->position[0] = ( i%2 == 0 ? (double)i * stepX : ((double)i+0.5) * stepX);
			newVoronoiDiagram->voronoiCells[index]->position[1] = (double) j * stepY;
			//newVoronoiDiagram->voronoiCells[index]->position[2] = strtod( ptr, &ptr );
		
			// set neighbors of voronoi cell
			newVoronoiDiagram->voronoiCells[index]->neighborCells = ( VoronoiCell**) calloc ( 6, sizeof( VoronoiCell*)); // allocate memory for list of neighbors
			
			// number of neighbours
			countNeighborCells = 0;
			newVoronoiDiagram->voronoiCells[index]->countNeighborCells = 6;



			// x-direction left neighbor
			ii = i-1;
			jj = j;
			if( i>0 || periodicX){
				ii = (ii+countPointsX)%countPointsX;
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[ii + jj*countPointsX];

			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;


			// x-direction right neighbor
			ii = i+1;
			jj = j;
			if( ii<countPointsX || periodicX ){
				ii = ii%countPointsX;
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[ii + jj*countPointsX];

			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;


			// y-direction lower-left neighbor
			ii = ( j%2 ? i : i-1);
			jj = j-1;
			if( (ii>=0 || periodicX) && (jj>=0 || periodicY)){
				ii = (ii+countPointsX)%countPointsX;
				jj = (jj+countPointsY)%countPointsY;
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[ii + jj*countPointsX];

			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;


			// y-direction lower-right neighbor
			ii = ( j%2 ? i+1 : i);
			jj = j-1;
			if( (ii<countPointsX || periodicX) && (jj>=0 || periodicY)){
				ii = ii%countPointsX;
				jj = (jj+countPointsY)%countPointsY;
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[ii + jj*countPointsX];

			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;


			// y-direction upper-left neighbor
			ii = ( j%2 ? i : i-1);
			jj = j+1;
			if( (ii>=0 || periodicX) && (jj<countPointsY || periodicY)){
				ii = (ii+countPointsX)%countPointsX;
				jj = jj%countPointsY;
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[ii + jj*countPointsX];

			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;

			// y-direction upper-right neighbor
			ii = ( j%2 ? i+1 : i);
			jj = j+1;
			if( (ii<countPointsX || periodicX) && (jj<countPointsY || periodicY)){
				ii = ii%countPointsX;
				jj = jj%countPointsY;
				newVoronoiDiagram->voronoiCells[index]->neighborCells[countNeighborCells++] = newVoronoiDiagram->voronoiCells[ii + jj*countPointsX];

			}else
				newVoronoiDiagram->voronoiCells[index]->countNeighborCells--;

		}
	}
	

	

	return newVoronoiDiagram;
}
/*****************************************************************************/
