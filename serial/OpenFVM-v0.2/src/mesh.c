/***************************************************************************
 *   Copyright (C) 2004-2006 by OpenFVM team                               *
 *   http://sourceforge.net/projects/openfvm/                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "globals.h"
#include "mesh.h"
#include "bcond.h"
#include "geocalc.h"
#include "ioutils.h"
#include "octree.h"

void MshGetElementTypes()
{

	int i;

	int element, face, pair;

	nbtris = 0;
	nbquads = 0;
	nbtetras = 0;
	nbhexas = 0;
	nbprisms = 0;

	for (i = 0; i < nbfaces; i++)
	{
		
		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue;

		if (faces[face].type == TRIANGLE)
		{
			nbtris++;	
		}

		if (faces[face].type == QUADRANGLE)
		{
			nbquads++;	
		}

	}

	for (i = 0; i < nbelements; i++)
	{	

		element = i;

		if (elements[element].type == TETRAHEDRON)
		{
			nbtetras++;	
		}

		if (elements[element].type == HEXAHEDRON)
		{
			nbhexas++;	
		}

		if (elements[element].type == PRISM)
		{
			nbprisms++;	
		}
	}

}

void MshCreateFaces()
{

	int i, j;

	int element;

	msh_vector d;

	// Create faces
	faces = realloc(faces, nbelements * 8 * sizeof(msh_face));

	nbfaces = 0;

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type == TRIANGLE)
		{

			elements[element].nbfaces = 1;

			elements[element].face[0] = nbfaces;

			faces[nbfaces].type = TRIANGLE;
			faces[nbfaces].nbnodes = 3;

			faces[nbfaces].node[0] = elements[element].node[0];
			faces[nbfaces].node[1] = elements[element].node[1];
			faces[nbfaces].node[2] = elements[element].node[2];

			// Calculate area
			faces[nbfaces].Aj = GeoCalcTriArea(nodes[faces[nbfaces].node[0]],
						        nodes[faces[nbfaces].node[1]],
						        nodes[faces[nbfaces].node[2]]);

			// Calculate centroid
			faces[nbfaces].cface = GeoCalcCentroid3(nodes[faces[nbfaces].node[0]],
													nodes[faces[nbfaces].node[1]],
													nodes[faces[nbfaces].node[2]]);

			// Calculate normal of the face				
			faces[nbfaces].n = GeoCalcNormal(nodes[faces[nbfaces].node[0]], 
											 nodes[faces[nbfaces].node[1]], 
											 nodes[faces[nbfaces].node[2]]); 

			faces[nbfaces].element = element;
			faces[nbfaces].physreg = -1;
			faces[nbfaces].elemreg = -1;
			faces[nbfaces].pair = -1;
			faces[nbfaces].bc = NONE;
		
			nbfaces++;
		
		}

		if (elements[element].type == QUADRANGLE)
		{

			elements[element].nbfaces = 1;

			elements[element].face[0] = nbfaces;

			faces[nbfaces].type = QUADRANGLE;
			faces[nbfaces].nbnodes = 4;

			faces[nbfaces].node[0] = elements[element].node[0];
			faces[nbfaces].node[1] = elements[element].node[1];
			faces[nbfaces].node[2] = elements[element].node[2];
			faces[nbfaces].node[3] = elements[element].node[3];

			// Calculate area
			faces[nbfaces].Aj = GeoCalcQuadArea(nodes[faces[nbfaces].node[0]],
					                  nodes[faces[nbfaces].node[1]],
					                  nodes[faces[nbfaces].node[2]],
					                  nodes[faces[nbfaces].node[3]]);

			// Calculate centroid			
			faces[nbfaces].cface = GeoCalcCentroid4(nodes[faces[nbfaces].node[0]],
													nodes[faces[nbfaces].node[1]],
													nodes[faces[nbfaces].node[2]],
													nodes[faces[nbfaces].node[3]]);

			// Calculate normal of the face				
			faces[nbfaces].n = GeoCalcNormal(nodes[faces[nbfaces].node[0]], 
											 nodes[faces[nbfaces].node[1]], 
											 nodes[faces[nbfaces].node[2]]); 

			faces[nbfaces].element = element;
			faces[nbfaces].physreg = -1;
			faces[nbfaces].elemreg = -1;
			faces[nbfaces].pair = -1;
			faces[nbfaces].bc = NONE;
	
			nbfaces++;
		
		}

		if (elements[element].type == TETRAHEDRON)
		{

			elements[element].nbfaces = 4;

			for (j = 0; j < 4; j++)
			{

				elements[element].face[j] = nbfaces;

				faces[nbfaces].type = TRIANGLE;
				faces[nbfaces].nbnodes = 3;

				switch (j)
				{
				case 0:

					faces[nbfaces].node[0] = elements[element].node[0];
					faces[nbfaces].node[1] = elements[element].node[1];
					faces[nbfaces].node[2] = elements[element].node[2];

					break;

				case 1:

					faces[nbfaces].node[0] = elements[element].node[1];
					faces[nbfaces].node[1] = elements[element].node[3];
					faces[nbfaces].node[2] = elements[element].node[2];

					break;

				case 2:

					faces[nbfaces].node[0] = elements[element].node[2];
					faces[nbfaces].node[1] = elements[element].node[3];
					faces[nbfaces].node[2] = elements[element].node[0];

					break;

				case 3:

					faces[nbfaces].node[0] = elements[element].node[3];
					faces[nbfaces].node[1] = elements[element].node[1];
					faces[nbfaces].node[2] = elements[element].node[0];

					break;
				}

				// Calculate area
				faces[nbfaces].Aj = GeoCalcTriArea(nodes[faces[nbfaces].node[0]],
								  nodes[faces[nbfaces].node[1]],
								  nodes[faces[nbfaces].node[2]]);

				// Calculate centroid
				faces[nbfaces].cface = GeoCalcCentroid3(nodes[faces[nbfaces].node[0]],
														nodes[faces[nbfaces].node[1]],
														nodes[faces[nbfaces].node[2]]);

				// Calculate normal of the face				
				faces[nbfaces].n = GeoCalcNormal(nodes[faces[nbfaces].node[0]], 
												 nodes[faces[nbfaces].node[1]], 
												 nodes[faces[nbfaces].node[2]]); 

				d = GeoSubVectorVector(faces[nbfaces].cface, elements[element].celement);

				// Normal should point to neighbour
				// Flip normal if necessary
				if (GeoDotVectorVector(faces[nbfaces].n, d) < 0.0)
				{
					faces[nbfaces].n.x *= -1;
					faces[nbfaces].n.y *= -1;
					faces[nbfaces].n.z *= -1;
				}

				faces[nbfaces].element = element;
				faces[nbfaces].physreg = -1;
				faces[nbfaces].elemreg = -1;
				faces[nbfaces].pair = -1;
				faces[nbfaces].bc = NONE;

				nbfaces++;

			}	

		}

		if (elements[element].type == HEXAHEDRON)
		{

			elements[element].nbfaces = 6;

			for (j = 0; j < 6; j++)
			{

				elements[element].face[j] = nbfaces;

				faces[nbfaces].type = QUADRANGLE;
				faces[nbfaces].nbnodes = 4;

				switch (j)
				{
				case 0:

					faces[nbfaces].node[0] = elements[element].node[0];
					faces[nbfaces].node[1] = elements[element].node[1];
					faces[nbfaces].node[2] = elements[element].node[2];
					faces[nbfaces].node[3] = elements[element].node[3];

					break;

				case 1:

					faces[nbfaces].node[0] = elements[element].node[7];
					faces[nbfaces].node[1] = elements[element].node[6];
					faces[nbfaces].node[2] = elements[element].node[5];
					faces[nbfaces].node[3] = elements[element].node[4];

					break;

				case 2:

					faces[nbfaces].node[0] = elements[element].node[5];
					faces[nbfaces].node[1] = elements[element].node[6];
					faces[nbfaces].node[2] = elements[element].node[2];
					faces[nbfaces].node[3] = elements[element].node[1];

					break;

				case 3:

					faces[nbfaces].node[0] = elements[element].node[3];
					faces[nbfaces].node[1] = elements[element].node[7];
					faces[nbfaces].node[2] = elements[element].node[4];
					faces[nbfaces].node[3] = elements[element].node[0];

					break;

				case 4:

					faces[nbfaces].node[0] = elements[element].node[0];
					faces[nbfaces].node[1] = elements[element].node[4];
					faces[nbfaces].node[2] = elements[element].node[5];
					faces[nbfaces].node[3] = elements[element].node[1];

					break;

				case 5:

					faces[nbfaces].node[0] = elements[element].node[7];
					faces[nbfaces].node[1] = elements[element].node[3];
					faces[nbfaces].node[2] = elements[element].node[2];
					faces[nbfaces].node[3] = elements[element].node[6];

					break;

				}

				// Calculate area
				faces[nbfaces].Aj = GeoCalcQuadArea(nodes[faces[nbfaces].node[0]],
							nodes[faces[nbfaces].node[1]],
							nodes[faces[nbfaces].node[2]],
							nodes[faces[nbfaces].node[3]]);

				// Calculate centroid	
				faces[nbfaces].cface = GeoCalcCentroid4(nodes[faces[nbfaces].node[0]],
														nodes[faces[nbfaces].node[1]],
														nodes[faces[nbfaces].node[2]],
														nodes[faces[nbfaces].node[3]]);

				// Calculate normal of the face				
				faces[nbfaces].n = GeoCalcNormal(nodes[faces[nbfaces].node[0]], 
												 nodes[faces[nbfaces].node[1]], 
												 nodes[faces[nbfaces].node[2]]); 

				d = GeoSubVectorVector(faces[nbfaces].cface, elements[element].celement);

				// Normal should point to neighbour
				// Flip normal if necessary
				if (GeoDotVectorVector(faces[nbfaces].n, d) < 0.0)
				{
					faces[nbfaces].n.x *= -1;
					faces[nbfaces].n.y *= -1;
					faces[nbfaces].n.z *= -1;
				}

				faces[nbfaces].element = element;
				faces[nbfaces].physreg = -1;
				faces[nbfaces].elemreg = -1;
				faces[nbfaces].pair = -1;
				faces[nbfaces].bc = NONE;

				nbfaces++;

			}	
		}
		
		if (elements[element].type == PRISM)
		{

			elements[element].nbfaces = 5;

			for (j = 0; j < 2; j++)
			{

				elements[element].face[j] = nbfaces;

				faces[nbfaces].type = TRIANGLE;
				faces[nbfaces].nbnodes = 3;

				switch (j)
				{
				case 0:

					faces[nbfaces].node[0] = elements[element].node[0];
					faces[nbfaces].node[1] = elements[element].node[1];
					faces[nbfaces].node[2] = elements[element].node[2];

					break;

				case 1:

					faces[nbfaces].node[0] = elements[element].node[3];
					faces[nbfaces].node[1] = elements[element].node[5];
					faces[nbfaces].node[2] = elements[element].node[4];

					break;

				}

				// Calculate area
				faces[nbfaces].Aj = GeoCalcTriArea(nodes[faces[nbfaces].node[0]], 
								nodes[faces[nbfaces].node[1]],
							          nodes[faces[nbfaces].node[2]]);

				// Calculate centroid	
				faces[nbfaces].cface = GeoCalcCentroid3(nodes[faces[nbfaces].node[0]],
														nodes[faces[nbfaces].node[1]],
														nodes[faces[nbfaces].node[2]]);

				// Calculate normal of the face				
				faces[nbfaces].n = GeoCalcNormal(nodes[faces[nbfaces].node[0]], 
												 nodes[faces[nbfaces].node[1]], 
												 nodes[faces[nbfaces].node[2]]); 

				d = GeoSubVectorVector(faces[nbfaces].cface, elements[element].celement);

				// Normal should point to neighbour
				// Flip normal if necessary
				if (GeoDotVectorVector(faces[nbfaces].n, d) < 0.0)
				{
					faces[nbfaces].n.x *= -1;
					faces[nbfaces].n.y *= -1;
					faces[nbfaces].n.z *= -1;
				}

				faces[nbfaces].element = element;
				faces[nbfaces].physreg = -1;
				faces[nbfaces].elemreg = -1;
				faces[nbfaces].pair = -1;

				nbfaces++;

			}	

			for (j = 2; j < 5; j++)
			{

				elements[element].face[j] = nbfaces;

				faces[nbfaces].type = QUADRANGLE;
				faces[nbfaces].nbnodes = 4;

				switch (j)
				{
				case 2:

					faces[nbfaces].node[0] = elements[element].node[0];
					faces[nbfaces].node[1] = elements[element].node[2];
					faces[nbfaces].node[2] = elements[element].node[5];
					faces[nbfaces].node[3] = elements[element].node[3];

					break;

				case 3:

					faces[nbfaces].node[0] = elements[element].node[1];
					faces[nbfaces].node[1] = elements[element].node[4];
					faces[nbfaces].node[2] = elements[element].node[5];
					faces[nbfaces].node[3] = elements[element].node[2];

					break;

				case 4:

					faces[nbfaces].node[0] = elements[element].node[0];
					faces[nbfaces].node[1] = elements[element].node[3];
					faces[nbfaces].node[2] = elements[element].node[4];
					faces[nbfaces].node[3] = elements[element].node[1];

					break;

				}

				// Calculate area 
				faces[nbfaces].Aj = GeoCalcQuadArea(nodes[faces[nbfaces].node[0]], nodes[faces[nbfaces].node[1]],
							          nodes[faces[nbfaces].node[2]], nodes[faces[nbfaces].node[3]]);

				// Calculate centroid	
				faces[nbfaces].cface = GeoCalcCentroid4(nodes[faces[nbfaces].node[0]],
								 nodes[faces[nbfaces].node[1]],
								 nodes[faces[nbfaces].node[2]],
								 nodes[faces[nbfaces].node[3]]);

				// Calculate normal of the face				
				faces[nbfaces].n = GeoCalcNormal(nodes[faces[nbfaces].node[0]], 
												 nodes[faces[nbfaces].node[1]], 
												 nodes[faces[nbfaces].node[2]]); 

				d = GeoSubVectorVector(faces[nbfaces].cface, elements[element].celement);

				// Normal should point to neighbour
				// Flip normal if necessary
				if (GeoDotVectorVector(faces[nbfaces].n, d) < 0.0)
				{
					faces[nbfaces].n.x *= -1;
					faces[nbfaces].n.y *= -1;
					faces[nbfaces].n.z *= -1;
				}

				faces[nbfaces].element = element;
				faces[nbfaces].physreg = -1;
				faces[nbfaces].elemreg = -1;
				faces[nbfaces].pair = -1;
				faces[nbfaces].bc = NONE;

				nbfaces++;

			}	

		}

	}

	faces = realloc(faces, nbfaces * sizeof(msh_face));

}

void MshAddBoundaryFaces()
{

	int i;

	nbfaces += nbpatches;
	
	faces = realloc(faces, nbfaces * sizeof(msh_face));

	for (i = nbfaces - nbpatches; i < nbfaces; i++)
	{
  
		faces[i] = patches[i - nbfaces + nbpatches];

	}

}

void MshConnectFaces()
{

	int i, j, k;

	double min[3], max[3];

	int nb;
	int nbv1, nbv2;
	int face, another_face;
	double dist;
	msh_vector cent[2];

	oct_data *tab = NULL;

	// Find pairs of all faces 
	// Boundary faces have no pair: -1
	// Create an octree with face cfaces

	if (nbfaces == 0)
		return;

	tab = malloc(nbfaces * sizeof(oct_data));

	min[0] = faces[0].cface.x;
	min[1] = faces[0].cface.y;
	min[2] = faces[0].cface.z;

	max[0] = faces[0].cface.x;
	max[1] = faces[0].cface.y;
	max[2] = faces[0].cface.z;
		
	for (i = 0; i < nbfaces; i++)
	{

		tab[i].x = faces[i].cface.x;
		tab[i].y = faces[i].cface.y;
		tab[i].z = faces[i].cface.z;

		min[0] = MIN(min[0], faces[i].cface.x);
		min[1] = MIN(min[1], faces[i].cface.y);
		min[2] = MIN(min[2], faces[i].cface.z);

		max[0] = MAX(max[0], faces[i].cface.x);
		max[1] = MAX(max[1], faces[i].cface.y);
		max[2] = MAX(max[2], faces[i].cface.z);

    	}
  
	OctCreateOctree(min, max, tab, nbfaces);

	for(i = 0; i < nbleafs; i++)
	{
	
		nb = leafs[i].nbentities;

		// For each node, we find the two corresponding faces
		for(j = 0; j < nb; j++)
		{

			face = leafs[i].entities[j];
			nbv1 = faces[face].nbnodes;
			cent[0] = faces[face].cface;

			if (faces[face].pair != -1) continue;

			for(k = j + 1; k < nb; k++)
			{

				another_face = leafs[i].entities[k];
				nbv2 = faces[another_face].nbnodes;
				cent[1] = faces[another_face].cface;
				
				dist = GeoMagVector(GeoSubVectorVector(cent[0], cent[1]));

				if (dist < EPSILON && nbv1 == nbv2 && face != another_face)
				{

					faces[face].pair = another_face;
					faces[another_face].pair = face;
     
		        	break;
				}
			

			}
		
			leafs[i].entities[k] = leafs[i].entities[nb-j-1];

		}
  	}

	OctDestroyOctree();

	free(tab);
}

void MshRemoveBoundaryFaces()
{

	int i;

	int face, pair;

	nbfaces -= nbpatches;

	for (i = 0; i < nbfaces; i++)
	{
		face = i;
	
		pair = faces[face].pair;

		if (pair != -1)
		{
			faces[face].physreg = faces[pair].physreg;  
			faces[face].elemreg = faces[pair].elemreg; 
			faces[face].bc = faces[pair].bc; 	
		}
	}

	for (i = 0; i < nbfaces; i++)
	{
		face = i;
		
		faces[face].pair = -1;
	}

	faces = realloc(faces, nbfaces * sizeof(msh_face));

}

void MshCalcPropMesh()
{

	int i, j;

	int element, face, pair;

	double sum_area;
	double sum_volume;

	double total_area;
	double total_volume;

	total_area = 0.0;
	total_volume = 0.0;

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		sum_area = 0.0;
		sum_volume = 0.0;

		// Volume calculation using Gauss theorem
		
		for (j = 0; j < elements[element].nbfaces; j++)
		{

			face = elements[element].face[j];

			sum_volume += faces[face].cface.x * faces[face].A.x;

			pair = faces[face].pair;

			if (pair == -1)
				sum_area += faces[face].Aj;

		}

		if (sum_volume <= 0.0)
		{
       		printf("\nError: Element with zero or negative volume\n");
			exit(LOGICAL_ERROR);
		}

		elements[element].Vp = sum_volume;

		total_volume += elements[i].Vp;
		total_area += sum_area; 

	}

	printf("Total surface area: \t\t\t\t%+.3E m^2\n", total_area);
	printf("Total volume: \t\t\t\t\t%+.3E m^3\n", total_volume);

}

void MshCalcPropFaces()
{

	int i, j;	

	int element, face, pair, neighbor;

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		faces[face].rpl = GeoSubVectorVector(faces[face].cface, GeoMultScalarVector(GeoDotVectorVector(GeoSubVectorVector(faces[face].cface, elements[element].celement), faces[face].n), faces[face].n));

		if (pair != -1)
		{

			neighbor = faces[pair].element;

			faces[face].rnl = GeoSubVectorVector(faces[face].cface, GeoMultScalarVector(GeoDotVectorVector(GeoSubVectorVector(faces[face].cface, elements[neighbor].celement), faces[face].n), faces[face].n));
	
			faces[face].A = GeoMultScalarVector(faces[face].Aj, faces[face].n);

			faces[face].d = GeoSubVectorVector(faces[face].rnl, faces[face].rpl);

			faces[face].dj = GeoMagVector(faces[face].d);
				
		}
		else
		{

			faces[face].A = GeoMultScalarVector(faces[face].Aj, faces[face].n);

			faces[face].d = GeoSubVectorVector(faces[face].cface, faces[face].rpl);

			faces[face].dj = GeoMagVector(faces[face].d);
							
		}

	}			

}

int MshImportMSH(char *file)
{
	int i, j;

	int nnod;
	int nele;

	int index;
	int maxindex;
	int etype;

	int physreg, elemreg;

	int ival;
	
	float fval;

	FILE *fp;
	char descr[512];

	fp = fopen(file, "r");

	if (fp == NULL)
	{
		printf("\nError: Mesh file not found!\n");
		printf("%s\n\n", file);
		exit(LOGICAL_ERROR);
	}

	printf("\nReading mesh...\n", file);

	maxindex = 0;
	
	do
	{

		do
		{
		
			strcpy(descr, "");
			
			fscanf(fp, "%s", descr);

			if (strcmp(descr, "$NOD") == 0)
				break;

			if (strcmp(descr, "$ELM") == 0)
				break;

		} while (!feof(fp));

		if (strcmp(descr, "$NOD") == 0)
		{
			fscanf(fp, "%d", &nnod);

			nodes = realloc(nodes, nnod * sizeof(msh_vector));
			
			nod_correlation_malloced = nnod;
			nod_correlation = realloc(nod_correlation, nod_correlation_malloced * sizeof(int));
	
			nbnodes = 0;
			
			for (i = 0; i < nnod; i++)
			{

				fscanf(fp, "%d", &ival);
				
				if (ival >= nod_correlation_malloced)
				{ 
					nod_correlation_malloced += ival + 1;
					nod_correlation = realloc(nod_correlation, nod_correlation_malloced * sizeof(int));
				}

				nod_correlation[ival] = i;
							
				fscanf(fp, "%f", &fval);
				nodes[nbnodes].x = fval;

				fscanf(fp, "%f", &fval);
				nodes[nbnodes].y = fval;

				fscanf(fp, "%f", &fval);
				nodes[nbnodes].z = fval;

				nbnodes++;

			}
		}

		if (strcmp(descr, "$ELM") == 0)
		{
		
			fscanf(fp, "%d", &nele);
			
			patches = realloc(patches, nele * sizeof(msh_face));
			elements = realloc(elements, nele * sizeof(msh_element));

			nbpatches = 0;
			nbelements = 0;
			
			for (i = 0; i < nele; i++)
			{
				
				fscanf(fp, "%d %d", &index, &etype);	

				if (etype < 1)
				{
					// Invalid element

         			printf("\nError: Invalid element\n");
					exit(LOGICAL_ERROR);

				}

				if (etype == 1)
				{

					// Beam element

         			printf("\nError: Invalid element\n");
					exit(LOGICAL_ERROR);

				}
				
				if (etype == 2 || etype == 3)
				{

					// Shell element or surface element

					fscanf(fp, "%d %d %d", &physreg, &elemreg, &ival);
										
					patches[nbpatches].nbnodes = ival;

					patches[nbpatches].cface.x = 0.0;
					patches[nbpatches].cface.y = 0.0;
					patches[nbpatches].cface.z = 0.0;

					for (j = 0; j < patches[nbpatches].nbnodes; j++)
					{ 
						fscanf(fp, "%d", &ival);

						patches[nbpatches].node[j] = nod_correlation[ival];

						patches[nbpatches].cface.x += nodes[patches[nbpatches].node[j]].x;
						patches[nbpatches].cface.y += nodes[patches[nbpatches].node[j]].y;
						patches[nbpatches].cface.z += nodes[patches[nbpatches].node[j]].z;
						
					}

					patches[nbpatches].cface.x /= patches[nbpatches].nbnodes;
					patches[nbpatches].cface.y /= patches[nbpatches].nbnodes;
					patches[nbpatches].cface.z /= patches[nbpatches].nbnodes;
					
					patches[nbpatches].physreg = physreg;
					patches[nbpatches].elemreg = elemreg;
					patches[nbpatches].element = -1;
					patches[nbpatches].pair = -1;
					patches[nbpatches].bc = NONE;
						
					if (etype == 2)
						patches[nbpatches].type = TRIANGLE;

					if (etype == 3)
						patches[nbpatches].type = QUADRANGLE;

					patches[nbpatches].index = index;

					nbpatches++;			

				}

				if (etype >= 4 && etype <= 6)
				{

					// Solid elements or volume element
					
					fscanf(fp, "%d %d %d", &physreg, &elemreg, &ival);
										
					elements[nbelements].nbnodes = ival;
					
					elements[nbelements].celement.x = 0.0;
					elements[nbelements].celement.y = 0.0;
					elements[nbelements].celement.z = 0.0;

					for (j = 0; j < elements[nbelements].nbnodes; j++)
					{ 
						fscanf(fp, "%d", &ival);

						elements[nbelements].node[j] = nod_correlation[ival];

						elements[nbelements].celement.x += nodes[elements[nbelements].node[j]].x;
						elements[nbelements].celement.y += nodes[elements[nbelements].node[j]].y;
						elements[nbelements].celement.z += nodes[elements[nbelements].node[j]].z;

					}

					elements[nbelements].celement.x /= elements[nbelements].nbnodes;
					elements[nbelements].celement.y /= elements[nbelements].nbnodes;
					elements[nbelements].celement.z /= elements[nbelements].nbnodes;

					elements[nbelements].physreg = physreg;
					elements[nbelements].elemreg = elemreg;

					if (etype == 4)
						elements[nbelements].type = TETRAHEDRON;

					if (etype == 5)
						elements[nbelements].type = HEXAHEDRON;

					if (etype == 6)
						elements[nbelements].type = PRISM;

					elements[nbelements].index = index;
					
					nbelements++;

					maxindex = MAX(index, maxindex);
					
				}

				if (etype >= 7)
				{

					// Unkown element

         			printf("\nError: Unkown element\n");
					exit(LOGICAL_ERROR);

				}
				
			}
				
		}

	}
	while (!feof(fp));

	fclose(fp);

	nodes = realloc(nodes, nbnodes * sizeof(msh_vector));
	patches = realloc(patches, nbpatches * sizeof(msh_face));
	elements = realloc(elements, nbelements * sizeof(msh_element));

	free(nod_correlation);

	// Assign patches for parallel processing
				
	for (i = 0; i < nbpatches; i++)
	{
		if (patches[i].index > maxindex) 
		{
			patches[i].bc = PROCESSOR;
		}
	}	
	
	printf("Done.\n");

	MshCreateFaces();
	MshAddBoundaryFaces();
	MshConnectFaces();
	MshRemoveBoundaryFaces();
	MshConnectFaces();
	MshCalcPropFaces();
	
	printf("\n");
	printf("Number of nodes: \t\t\t\t%d\n", nbnodes);
	printf("Number of faces: \t\t\t\t%d\n", nbfaces);
	printf("Number of boundary faces: \t\t\t%d\n", nbpatches);
	printf("Number of elements: \t\t\t\t%d\n", nbelements);

	if (nbelements == 0)
	{
		printf("\nError: No solid elements.\n");
		exit(LOGICAL_ERROR);
	}

	MshCalcPropMesh();
	MshGetElementTypes();
		
	return LOGICAL_TRUE;
	
}

int MshExportMSH(char *file)
{

	int i, j;

	int node, patch, element;
	
	FILE *fp;
			
	fp = fopen(file, "w");

	if (fp == NULL) return LOGICAL_FALSE;

	fprintf(fp, "$NOD\n");
	fprintf(fp, "%d\n", nbnodes);

	for (i = 0; i < nbnodes; i++)
	{
		node = i;

		fprintf(fp, "%d %f %f %f", node + 1, nodes[node].x, nodes[node].y, nodes[node].z);
		fprintf(fp, "\n");
	}

	fprintf(fp, "$ENDNOD\n");

	fprintf(fp, "$ELM\n");	

	fprintf(fp, "%d\n", nbpatches + nbelements);
		
	for (i = 0; i < nbpatches; i++)
	{

		patch = i;

		fprintf(fp, "%d", patches[patch].index);

		if (patches[patch].type == TRIANGLE)
			fprintf(fp, " %d", 2);

		if (patches[patch].type == QUADRANGLE)
			fprintf(fp, " %d", 3);

		fprintf(fp, " %d", patches[patch].physreg);
		fprintf(fp, " %d", patches[patch].elemreg);

		fprintf(fp, " %d", patches[patch].nbnodes);

		for (j = 0; j < patches[patch].nbnodes; j++)
		{
			fprintf(fp, " %d", patches[patch].node[j] + 1);
		}

		fprintf(fp, "\n");
	}

	for (i = 0; i < nbelements; i++)
	{

		element = i;
		
		fprintf(fp, "%d", elements[element].index);

		if (elements[element].type == TETRAHEDRON)
			fprintf(fp, " %d", 4);
	
		if (elements[element].type == HEXAHEDRON)
			fprintf(fp, " %d", 5);
		
		if (elements[element].type == PRISM)
			fprintf(fp, " %d", 6);

		fprintf(fp, " %d", elements[element].physreg);
		fprintf(fp, " %d", elements[element].elemreg);

		fprintf(fp, " %d", elements[element].nbnodes);

		for (j = 0; j < elements[element].nbnodes; j++)
		{
			fprintf(fp, " %d", elements[element].node[j] + 1);
		}

		fprintf(fp, "\n");
	}

	fprintf(fp, "$ENDELM\n");

	fclose(fp);

	return LOGICAL_TRUE;

}

int MshExportDecomposedMSH(char *file, int region)
{

	int i, j, n;

	int node, face, pair, patch, element, neighbor;
	
	FILE *fp;

	int *nflag;

	int nbnewnodes;
	int nbnewpatches;
	int nbnewelements;
	
	msh_vector *newnodes;
	msh_face *newpatches;
	msh_element *newelements;
	
	int *newindex;
						
	// Identify elements in region
	
	newelements = calloc(nbelements, sizeof(msh_element));
		
	nbnewelements = 0;
	
	for (i = 0; i < nbelements; i++)
	{
		
		element = i;
		
		if (elements[element].elemreg != region) continue;
		
		newelements[nbnewelements] = elements[element];
		nbnewelements++;
		 
	}

	if (nbnewelements == 0)
	{
		printf("\nError: Number of elements is zero\n"); 
		return LOGICAL_ERROR;
	}

	realloc(newelements, nbnewelements * sizeof(msh_element));
	
	// Remove unconnected and merged nodes
		
	nflag = calloc(nbnodes, sizeof(int));

	for (i = 0; i < nbnodes; i++)
	{
		node = i;
	
		nflag[node] = LOGICAL_FALSE;
	}

	for (i = 0; i < nbnewelements; i++)
	{
	
		element = i;

		for (j = 0; j < newelements[element].nbnodes; j++)
		{
			node = newelements[element].node[j];
			
			nflag[node] = LOGICAL_TRUE;
		}

	}
	
	newnodes = calloc(nbnodes, sizeof(msh_vector));
	
	newindex = calloc(nbnodes, sizeof(int));

	nbnewnodes = 0;

	for (i = 0; i < nbnodes; i++)
	{
	
		node = i;
	
		if (nflag[node] == LOGICAL_TRUE)
		{

			newnodes[nbnewnodes] = nodes[node];

			newindex[node] = nbnewnodes;

			nbnewnodes++;

		}
	}
	
	realloc(newnodes, nbnewnodes * sizeof(msh_vector));
	
	// Identify patches in region
	
	newpatches = calloc(nbpatches + nbfaces, sizeof(msh_face));
	
	nbnewpatches = 0;
	
	for (i = 0; i < nbpatches; i++)
	{
		
		patch = i;
	
		n = LOGICAL_TRUE;

		for (j = 0; j < patches[patch].nbnodes; j++)
		{
			node = patches[patch].node[j];
			
			if (nflag[node] == LOGICAL_FALSE) n = LOGICAL_FALSE;
		}
		
		if (n == LOGICAL_FALSE) continue;

		newpatches[nbnewpatches] = patches[patch];
		nbnewpatches++;		
	}
	
	// Add patches which belong to cut 
	
	n = 0;
	
	for (i = 0; i < nbfaces; i++)
	{
	
		face = i;
		
		element = faces[face].element;
		
		pair = faces[face].pair;
		
		if (pair != -1)
		{
		
			neighbor = faces[pair].element;
			
			if (elements[element].elemreg != region) continue;
			
			if (elements[element].elemreg == elements[neighbor].elemreg) continue;
			
			newpatches[nbnewpatches] = faces[face];

			// Set physical region pointing to neighbor
			newpatches[nbnewpatches].physreg = elements[neighbor].index;
				
			// Set element region pointing to neighbor region
			newpatches[nbnewpatches].elemreg = elements[neighbor].elemreg; 
				
			// Set index to the end of the list of elements
			newpatches[nbnewpatches].index = nbelements + nbpatches + n + 1;
				
			n++;
				
			nbnewpatches++;
				
		}
			
	}
	
	realloc(newpatches, nbnewpatches * sizeof(msh_face));
		
	// Renumber patches and elements 	

	for (i = 0; i < nbnewelements; i++)
	{
	
		element = i;
		
		for (j = 0; j < newelements[element].nbnodes; j++)
		{
			node = newelements[element].node[j];

			newelements[element].node[j] = newindex[node]; 

		}
	}

	for (i = 0; i < nbnewpatches; i++)
	{
	
		patch = i;
	
		for (j = 0; j < newpatches[patch].nbnodes; j++)
		{
			node = newpatches[patch].node[j];

			newpatches[patch].node[j] = newindex[node]; 

		}
	}

	fp = fopen(file, "w");

	if (fp == NULL) return LOGICAL_FALSE;
	
	fprintf(fp, "$NOD\n");
	fprintf(fp, "%d\n", nbnewnodes);
	
	for (i = 0; i < nbnewnodes; i++)
	{
		node = i;

		fprintf(fp, "%d %f %f %f", node + 1, newnodes[node].x, newnodes[node].y, newnodes[node].z);
				
		fprintf(fp, "\n");
	}

	fprintf(fp, "$ENDNOD\n"); 

	fprintf(fp, "$ELM\n");	

	fprintf(fp, "%d\n", nbnewpatches + nbnewelements);
		
	for (i = 0; i < nbnewpatches; i++)
	{

		patch = i;

		if (newpatches[patch].index >= nbelements + nbpatches) continue;
		
		fprintf(fp, "%d", newpatches[patch].index);

		if (newpatches[patch].type == TRIANGLE)
			fprintf(fp, " %d", 2);

		if (newpatches[patch].type == QUADRANGLE)
			fprintf(fp, " %d", 3);

		fprintf(fp, " %d", newpatches[patch].physreg);
		fprintf(fp, " %d", newpatches[patch].elemreg);

		fprintf(fp, " %d", newpatches[patch].nbnodes);

		for (j = 0; j < newpatches[patch].nbnodes; j++)
		{
			fprintf(fp, " %d", newpatches[patch].node[j] + 1);
		}

		fprintf(fp, "\n");
	}

	for (i = 0; i < nbnewelements; i++)
	{

		element = i;
		
		fprintf(fp, "%d", newelements[element].index);

		if (newelements[element].type == TETRAHEDRON)
			fprintf(fp, " %d", 4);
	
		if (newelements[element].type == HEXAHEDRON)
			fprintf(fp, " %d", 5);
		
		if (newelements[element].type == PRISM)
			fprintf(fp, " %d", 6);

		fprintf(fp, " %d", newelements[element].physreg);
		fprintf(fp, " %d", newelements[element].elemreg);

		fprintf(fp, " %d", newelements[element].nbnodes);

		for (j = 0; j < newelements[element].nbnodes; j++)
		{
			fprintf(fp, " %d", newelements[element].node[j] + 1);
		}

		fprintf(fp, "\n");
	}

	for (i = 0; i < nbnewpatches; i++)
	{

		patch = i;

		if (newpatches[patch].index < nbelements + nbpatches) continue;
		
		fprintf(fp, "%d", newpatches[patch].index);

		if (newpatches[patch].type == TRIANGLE)
			fprintf(fp, " %d", 2);

		if (newpatches[patch].type == QUADRANGLE)
			fprintf(fp, " %d", 3);

		fprintf(fp, " %d", newpatches[patch].physreg);
		fprintf(fp, " %d", newpatches[patch].elemreg);

		fprintf(fp, " %d", newpatches[patch].nbnodes);

		for (j = 0; j < newpatches[patch].nbnodes; j++)
		{
			fprintf(fp, " %d", newpatches[patch].node[j] + 1);
		}

		fprintf(fp, "\n");
	}
	
	fprintf(fp, "$ENDELM\n");

	fclose(fp);

	free(nflag);
	free(newnodes);
	free(newpatches);
	free(newelements);
	free(newindex);	
	
	return LOGICAL_TRUE;

}


