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
#include <math.h>

#include "ioutils.h"
#include "globals.h"
#include "mesh.h"
#include "param.h"
#include "bcond.h"
#include "geocalc.h"
#include "gradient.h"

#include "laspack/itersolv.h"
#include "laspack/rtc.h"
#include "laspack/errhandl.h"

// Original streamFunction.C developed by OpenFOAM
// translated to C by the OpenFVM team in 08/02/2006

void CalculateStreamFunction(Vector *uf, double *streamFunction)
{

	int i, j;
	
	int face, pair, element;

	int found, finished;
	
	int *visitedPoint;
	int nVisited, nVisitedOld;
	
	int bPointFound, pointFound;
	
	double currentBStream, currentStream;
	
	double length;
		
	msh_vector currentBStreamPoint, currentStreamPoint;
	msh_vector edgeHat;

	double vmin, vmax;
	
	visitedPoint = calloc(nbnodes, sizeof(int));
	
	nVisited = 0;
	nVisitedOld = 0;
	
	finished = LOGICAL_TRUE;
		
	do
	{

		// Find the boundary face with zero flux. set the stream function
		// to zero on that face
			
		found = LOGICAL_FALSE;
		
		// Boundary faces
		
		for (i = 0; i < nbfaces; i++)
		{
		
			face = i;
			
			pair = faces[face].pair;
			
			if (pair != -1) continue;

			if (faces[face].bc == EMPTY) continue;
			
			if (ABS(V_GetCmp(uf, face + 1)) < SMALL)
			{

				// Zero flux face found
				found = LOGICAL_TRUE;
				
				for (j = 0; j < faces[face].nbnodes; j++)
				{
					
					if (visitedPoint[faces[face].node[j]] == 1)
					{
						found = LOGICAL_FALSE;
						break;
					}
				
				}
				
				if (found == LOGICAL_TRUE)
				{
				
					//printf("Zero face: %d\n", face);		
					
					for (j = 0; j < faces[face].nbnodes; j++)
					{
						streamFunction[faces[face].node[j]] = 0.0;
						visitedPoint[faces[face].node[j]] = 1;
						nVisited++;
					}	
					
					break;
				}
				
				
			}
			
			if (found == LOGICAL_TRUE) break;
			
		}
		
		if (found == LOGICAL_FALSE)
		{
			
			for (i = 0; i < nbelements; i++)
			{
			
				element = i;
				
				found = LOGICAL_TRUE;
			
				for (j = 0; j < elements[element].nbnodes; j++)
				{
					
					if (visitedPoint[elements[element].node[j]] == 1)
					{
						found = LOGICAL_FALSE;
						break;
					}
				
				}
				
				if (found == LOGICAL_TRUE)
				{
				
					for (j = 0; j < elements[element].nbnodes; j++)
					{
						streamFunction[elements[element].node[j]] = 0.0;
						visitedPoint[elements[element].node[j]] = 1;
						nVisited++;
					}
					
					break;	
				}
				
				
			}
		
		}
		
		// Loop through all faces. If one of the points on
		// the face has the streamFunction value different
		// from -1, all points with -1 on that face have the
		// streamFunction value equal to the face flux in
		// that point plus the value in the visited point
		
		do
		{
			
			finished = LOGICAL_TRUE;
		     
			// Boundary faces
		     
			for (i = 0; i < nbfaces; i++)
			{
		
				face = i;
			
				pair = faces[face].pair;
			
				if (pair != -1) continue;
			
				bPointFound = LOGICAL_FALSE;
				
				currentBStream = 0.0;
				
				currentBStreamPoint.x = 0.0;
				currentBStreamPoint.y = 0.0;
				currentBStreamPoint.z = 0.0;

				for (j = 0; j < faces[face].nbnodes; j++)
				{
					
					// Check if the point has been visited
					
					if (visitedPoint[faces[face].node[j]] == 1)
					{
						// The point has been visited
						currentBStream = streamFunction[faces[face].node[j]];
						currentBStreamPoint = nodes[faces[face].node[j]];

						bPointFound = LOGICAL_TRUE;

						break;
					}
				
				    
				}
				
				if (bPointFound == LOGICAL_TRUE)
				{
					
					// Sort out other points on the face				
					for (j = 0; j < faces[face].nbnodes; j++)
					{
					
						// Check if the point has been visited
						if (visitedPoint[faces[face].node[j]] == 0)
						{

							if (faces[face].bc != EMPTY)
							{
												
								edgeHat = GeoSubVectorVector(nodes[faces[face].node[j]], currentBStreamPoint);
							
								edgeHat.z = 0.0;
					
								edgeHat = GeoNormalizeVector(edgeHat);
														
								if (edgeHat.y > VSMALL)
								{
									visitedPoint[faces[face].node[j]] = 1;
									nVisited++;
							
									streamFunction[faces[face].node[j]] = currentBStream + V_GetCmp(uf, face + 1) * faces[face].Aj * SGN(faces[face].n.x);
								}
								else if (edgeHat.y < -VSMALL)
								{
									visitedPoint[faces[face].node[j]] = 1;
									nVisited++;

									streamFunction[faces[face].node[j]] = currentBStream - V_GetCmp(uf, face + 1) * faces[face].Aj * SGN(faces[face].n.x);
								}
								else
								{
									if (edgeHat.x > VSMALL)
									{
										visitedPoint[faces[face].node[j]] = 1;
										nVisited++;
									
										streamFunction[faces[face].node[j]] = currentBStream + V_GetCmp(uf, face + 1) * faces[face].Aj * SGN(faces[face].n.y);
									}
									else if (edgeHat.x < -VSMALL)
									{
										visitedPoint[faces[face].node[j]] = 1;
										nVisited++;

										streamFunction[faces[face].node[j]] = currentBStream - V_GetCmp(uf, face + 1) * faces[face].Aj * SGN(faces[face].n.y);
									}
								}
							}
						}
						
					}
					
					
				}
				else
				{
					finished = LOGICAL_FALSE;
				}
			
			}
				
			// Internal faces
			
			for (i = 0; i < nbfaces; i++)
			{
		
				face = i;
								
				pair = faces[face].pair;
			
				if (pair == -1) continue;
					
				pointFound = LOGICAL_FALSE;
				
				currentStream = 0.0;
				
				currentStreamPoint.x = 0.0;
				currentStreamPoint.y = 0.0;
				currentStreamPoint.z = 0.0;
				
				for (j = 0; j < faces[face].nbnodes; j++)
				{
				
					// Check if the point has been visited
					if (visitedPoint[faces[face].node[j]] == 1)
					{
						// The point has been visited
						currentStream = streamFunction[faces[face].node[j]];
						currentStreamPoint = nodes[faces[face].node[j]];

						pointFound = LOGICAL_TRUE;

						break;
					}
										
				}
				
				if (pointFound == LOGICAL_TRUE)
				{
					
					// Sort out other points on the face				
					for (j = 0; j < faces[face].nbnodes; j++)
					{
					
						// Check if the point has been visited
						if (visitedPoint[faces[face].node[j]] == 0)
						{
						
							edgeHat = GeoSubVectorVector(nodes[faces[face].node[j]], currentStreamPoint);
							
							edgeHat.z = 0.0;
					
							edgeHat = GeoNormalizeVector(edgeHat);
																		
							if (edgeHat.y > VSMALL)
							{
								visitedPoint[faces[face].node[j]] = 1;
								nVisited++;
							
								streamFunction[faces[face].node[j]] = currentStream + V_GetCmp(uf, face + 1) * faces[face].Aj * SGN(faces[face].n.x);
							}
							else if (edgeHat.y < -VSMALL)
							{
								visitedPoint[faces[face].node[j]] = 1;
								nVisited++;

								streamFunction[faces[face].node[j]] = currentStream - V_GetCmp(uf, face + 1) * faces[face].Aj * SGN(faces[face].n.x);
							}
						}
						
					}
					
					
				}
				else
				{
					finished = LOGICAL_FALSE;
				}
						
				
			}
			
			if (nVisited == nVisitedOld)
			{
			
				//printf("Exhausted a seed. Looking for new seed.\n");
			     
				break;
			}
			else
			{
				nVisitedOld = nVisited;
			}			
				
		} while (finished == LOGICAL_FALSE);
				
		
	} while (finished == LOGICAL_FALSE);	

	free(visitedPoint);

	// Get maximum and minimum values	
	vmin = +VGREAT;
	vmax = -VGREAT;

	for (i = 0; i < nbnodes; i++)
	{
		
		vmin = MIN(vmin, streamFunction[i]);		
		vmax = MAX(vmax, streamFunction[i]);
	
	}
	
	// Normalize stream function
	
	if (vmax != 0.0)
	{
		for (i = 0; i < nbnodes; i++)
		{
		
			streamFunction[i] = (streamFunction[i] - vmin) / (vmax - vmin); 
			//streamFunction[i] /= vmax; 
		
		}
	}
	
}

void GetNodeScalarFace(double *sn, int *si, double *sf)
{

	int i, j;

	int face, pair;
	
	for (i = 0; i < nbfaces; i++)
	{	

		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue;
		
		for (j = 0; j < faces[face].nbnodes; j++)
	    {	  
			sn[faces[face].node[j]] += sf[face];
			si[faces[face].node[j]]++;
	    }
	}

}

void GetNodeScalarElement(double *sn, int *si, double *se)
{

	int i, j;

	int element;
	
	for (i = 0; i < nbelements; i++)
	{	

		element = i;

		for (j = 0; j < elements[element].nbnodes; j++)
		{	  
			sn[elements[element].node[j]] += se[element];
			si[elements[element].node[j]]++;
	    	}
	}

}

void GetNodeScalarElementFace(double *sn, int *si, double *se, double *sf, int k)
{

	int i, j;

	int face, pair, element;
	
	for (i = 0; i < nbelements; i++)
	{	

		element = i;

		for (j = 0; j < elements[element].nbnodes; j++)
		{	  
			sn[elements[element].node[j]] += se[element];
			si[elements[element].node[j]]++;
	    	}
	}

	for (i = 0; i < nbfaces; i++)
	{	

		face = i;

		// Do not consider EMPTY bc and ADIABATIC WALL bc when probing temperature
		if (faces[face].bc == INLET || 
		    faces[face].bc == MOVINGWALL || 
		    faces[face].bc == WALL ||
		    faces[face].bc == SURFACE ||
		    (faces[face].bc == ADIABATICWALL && k != 4))
		{
			pair = faces[face].pair;
		
			if (pair != -1) continue;
		
			for (j = 0; j < faces[face].nbnodes; j++)
	 		{	  
				sn[faces[face].node[j]] = 0.0;
				si[faces[face].node[j]] = 0;
			}
	    	}
				
	}	
	for (i = 0; i < nbfaces; i++)
	{	

		face = i;

		if (faces[face].bc == INLET || 
		    faces[face].bc == MOVINGWALL || 
		    faces[face].bc == WALL ||
		    faces[face].bc == SURFACE ||
		    (faces[face].bc == ADIABATICWALL && k != 4))
		{
							
			pair = faces[face].pair;
		
			if (pair != -1) continue;
		
			for (j = 0; j < faces[face].nbnodes; j++)
	 		{	  
				sn[faces[face].node[j]] += sf[face];
				si[faces[face].node[j]]++;
			}
	    	}
		
	}
	
}

void GetNodeVectorElement(Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
						  double *sn, int *si, int iphix, int iphiy, int iphiz)
{

	int i, j;

	int element;
	
	double vx, vy, vz;

	for (i = 0; i < nbelements; i++)
	{	

		element = i;

		if (iphix == 0) vx = V_GetCmp(xu, element + 1);
		if (iphix == 1) vx = V_GetCmp(xv, element + 1);
		if (iphix == 2) vx = V_GetCmp(xw, element + 1);
		if (iphix == 3) vx = V_GetCmp(xp, element + 1);
		if (iphix == 4) vx = V_GetCmp(xT, element + 1);
		if (iphix == 5) vx = V_GetCmp(xs, element + 1);

		if (iphiy == 0) vy = V_GetCmp(xu, element + 1);
		if (iphiy == 1) vy = V_GetCmp(xv, element + 1);
		if (iphiy == 2) vy = V_GetCmp(xw, element + 1);
		if (iphiy == 3) vy = V_GetCmp(xp, element + 1);
		if (iphiy == 4) vy = V_GetCmp(xT, element + 1);
		if (iphiy == 5) vy = V_GetCmp(xs, element + 1);

		if (iphiz == 0) vz = V_GetCmp(xu, element + 1);
		if (iphiz == 1) vz = V_GetCmp(xv, element + 1);
		if (iphiz == 2) vz = V_GetCmp(xw, element + 1);
		if (iphiz == 3) vz = V_GetCmp(xp, element + 1);
		if (iphiz == 4) vz = V_GetCmp(xT, element + 1);
		if (iphiz == 5) vz = V_GetCmp(xs, element + 1);

		for (j = 0; j < elements[element].nbnodes; j++)
	    {	  
			sn[elements[element].node[j]] += sqrt(vx * vx + vy * vy + vz * vz);
			si[elements[element].node[j]]++;
	    }
	}

}

void WriteProbeViews(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
			       Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf,
			       char *var, double time)
{

	int i, j, k;

	int node, face, element;
	
	double vs;
	double sd;

	double *se, *sf;
	double *sn;
	int *si;

	se = calloc(nbelements, sizeof(double));
	sf = calloc(nbfaces, sizeof(double));
	sn = calloc(nbnodes, sizeof(double));
	si = calloc(nbnodes, sizeof(int));

	// Create probe views
	for (k = 0; k < nphi; k++)
	{
	
		if (parameter.probe[k] == LOGICAL_FALSE) continue; 

		for (i = 0; i < nbelements; i++)
		{
	
			element = i;
			
			if (k == 0) vs = V_GetCmp(xu, element + 1);
			if (k == 1) vs = V_GetCmp(xv, element + 1);
			if (k == 2) vs = V_GetCmp(xw, element + 1);
			if (k == 3) vs = V_GetCmp(xp, element + 1);
			if (k == 4) vs = V_GetCmp(xT, element + 1);
			if (k == 5) vs = V_GetCmp(xs, element + 1);

			se[element] = vs;
		
		}

		for (i = 0; i < nbfaces; i++)
		{
	
			face = i;
			
			if (k == 0) vs = V_GetCmp(xuf, face + 1);
			if (k == 1) vs = V_GetCmp(xvf, face + 1);
			if (k == 2) vs = V_GetCmp(xwf, face + 1);
			if (k == 3) vs = V_GetCmp(xpf, face + 1);
			if (k == 4) vs = V_GetCmp(xTf, face + 1);
			if (k == 5) vs = V_GetCmp(xsf, face + 1);

			sf[face] = vs;
		
		}
				
		GetNodeScalarElementFace(sn, si, se, sf, k);
				
		fprintf(fp, "View \"Variable: %c\"{\n", var[k]);
		
		fprintf(fp, "TIME { %f };\n", time);
	
		for (i = 0; i < nbelements; i++)
		{

			element = i;
			
			if (elements[element].type == TETRAHEDRON)
			{
				fprintf(fp, "SS(");
				
				for (j = 0; j < elements[element].nbnodes; j++)
				{
					node = elements[element].node[j];
					
					fprintf(fp, "%f, %f, %f", nodes[node].x, nodes[node].y, nodes[node].z);
					
					if (j != elements[element].nbnodes - 1)
						fprintf(fp, ",");

				}
				
				fprintf(fp, ")");
				
				fprintf(fp, "{");
				
				for (j = 0; j < elements[element].nbnodes; j++)
	    			{
		
					sd = sn[elements[element].node[j]] / si[elements[element].node[j]];

					fprintf(fp, " %f", sd);
					
					if (j != elements[element].nbnodes - 1)
						fprintf(fp, ",");
					
				}
					
				fprintf(fp, "};\n");				
			}
			
			if (elements[element].type == HEXAHEDRON)	
			{
				fprintf(fp, "SH(");
				
				for (j = 0; j < elements[element].nbnodes; j++)
				{
					node = elements[element].node[j];
					
					fprintf(fp, " %f, %f, %f", nodes[node].x, nodes[node].y, nodes[node].z);
					
					if (j != elements[element].nbnodes - 1)
						fprintf(fp, ",");					
				}
				
				fprintf(fp, ")");
				
				fprintf(fp, "{");
				
				for (j = 0; j < elements[element].nbnodes; j++)
	    			{
		
					sd = sn[elements[element].node[j]] / si[elements[element].node[j]];

					fprintf(fp, " %f", sd);
					
					if (j != elements[element].nbnodes - 1)
						fprintf(fp, ",");
					
				}
					
				fprintf(fp, "};\n");
				
			}
					
			if (elements[element].type == PRISM)
			{
			
				fprintf(fp, "SI(");
				
				for (j = 0; j < elements[element].nbnodes; j++)
				{
					node = elements[element].node[j];
					
					fprintf(fp, " %f, %f, %f", nodes[node].x, nodes[node].y, nodes[node].z);
					
					if (j != elements[element].nbnodes - 1)
						fprintf(fp, ",");
				}
				
				fprintf(fp, ")");
				
				fprintf(fp, "{");
				
				for (j = 0; j < elements[element].nbnodes; j++)
	    			{
		
					sd = sn[elements[element].node[j]] / si[elements[element].node[j]];

					fprintf(fp, " %f", sd);
					
					if (j != elements[element].nbnodes - 1)
						fprintf(fp, ",");
					
				}
					
				fprintf(fp, "};\n");
			}
		}
				
		fprintf(fp, "};\n");

	}

	free(se);
	free(sf);
	free(sn);
	free(si);
	
}

void PrintAsciiHeaderScalarFace(FILE *fp, char *label)
{

	fprintf(fp, "%s %d\n", label, 1);
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-points nb-vector-points nb-tensor-points
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-lines nb-vector-lines nb-tensor-lines
	fprintf(fp, "%d %d %d\n", nbtris, 0, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
	fprintf(fp, "%d %d %d\n", nbquads, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
	fprintf(fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void PrintAsciiHeaderScalarElement(FILE *fp, char *label)
{


	fprintf(fp, "%s %d\n", label, 1);
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-points nb-vector-points nb-tensor-points
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-lines nb-vector-lines nb-tensor-lines
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
	fprintf(fp, "%d %d %d\n", nbtetras, 0, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
	fprintf(fp, "%d %d %d\n", nbhexas, 0, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
	fprintf(fp, "%d %d %d\n", nbprisms, 0, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
	fprintf(fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void PrintAsciiHeaderVectorElement(FILE *fp, char *label)
{


	fprintf(fp, "%s %d\n", label, 1);
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-points nb-vector-points nb-tensor-points
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-lines nb-vector-lines nb-tensor-lines
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
	fprintf(fp, "%d %d %d\n", 0, nbtetras, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
	fprintf(fp, "%d %d %d\n", 0, nbhexas, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
	fprintf(fp, "%d %d %d\n", 0, nbprisms, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
	fprintf(fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void PrintAsciiHeaderVectorFace(FILE *fp, char *label)
{

	fprintf(fp, "%s %d\n", label, 1);
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-points nb-vector-points nb-tensor-points
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-lines nb-vector-lines nb-tensor-lines
	fprintf(fp, "%d %d %d\n", 0, nbtris, 0);	//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
	fprintf(fp, "%d %d %d\n", 0, nbquads, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
	fprintf(fp, "%d %d %d\n", 0, 0, 0);		//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
	fprintf(fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

}

void PrintBinaryHeaderScalarFace(FILE *fp, char *label)
{

	int one = 1;

	fprintf(fp, "%s %d ", label, 1);		//view-name nb-time-steps,
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-points nb-vector-points nb-tensor-points
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-lines nb-vector-lines nb-tensor-lines
	fprintf(fp, "%d %d %d ", nbtris, 0, 0);		//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
	fprintf(fp, "%d %d %d ", nbquads, 0, 0);	//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
	fprintf(fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

	fwrite(&one, sizeof(int), 1, fp);

}

void PrintBinaryHeaderScalarElement(FILE *fp, char *label)
{

	int one = 1;

	fprintf(fp, "%s %d ", label, 1);		//view-name nb-time-steps,
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-points nb-vector-points nb-tensor-points
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-lines nb-vector-lines nb-tensor-lines
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
	fprintf(fp, "%d %d %d ", nbtetras, 0, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
	fprintf(fp, "%d %d %d ", nbhexas, 0, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
	fprintf(fp, "%d %d %d ", nbprisms, 0, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
	fprintf(fp, "%d %d %d ", 0, 0, 0);		//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
	fprintf(fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

	fwrite(&one, sizeof(int), 1, fp);

}

void PrintBinaryHeaderVectorElement(FILE *fp, char *label)
{

	int one = 1;

	fprintf(fp, "%s %d ", label, 1);			//view-name nb-time-steps,
	fprintf(fp, "%d %d %d ", 0, 0, 0);			//nb-scalar-points nb-vector-points nb-tensor-points
	fprintf(fp, "%d %d %d ", 0, 0, 0);			//nb-scalar-lines nb-vector-lines nb-tensor-lines
	fprintf(fp, "%d %d %d ", 0, 0, 0);			//nb-scalar-triangles nb-vector-triangles nb-tensor-triangles
	fprintf(fp, "%d %d %d ", 0, 0, 0);			//nb-scalar-quadrangles nb-vector-quadrangles nb-tensor-quadrangles
	fprintf(fp, "%d %d %d ", 0, nbtetras, 0);	//nb-scalar-tetrahedra nb-vector-tetrahedra nb-tensor-tetrahedra
	fprintf(fp, "%d %d %d ", 0, nbhexas, 0);	//nb-scalar-hexahedra nb-vector-hexahedra nb-tensor-hexahedra
	fprintf(fp, "%d %d %d ", 0, nbprisms, 0);	//nb-scalar-prisms nb-vector-prisms nb-tensor-prisms
	fprintf(fp, "%d %d %d ", 0, 0, 0);			//nb-scalar-pyramids nb-vector-pyramids nb-tensor-pyramids
	fprintf(fp, "%d %d %d %d\n", 0, 0, 0, 0);	//nb-text2d nb-text2d-chars nb-text3d nb-text3d-chars

	fwrite(&one, sizeof(int), 1, fp);

}
	
void PrintAsciiTime(FILE *fp, double time)
{

	fprintf(fp, "%E\n", time);

}

void PrintBinaryTime(FILE *fp, double time)
{

	double *times;

	times = calloc(1, sizeof(double));

	times[0] = time;

	fwrite(times, sizeof(double), 1, fp);

	free(times);

}

void PrintAsciiScalarNode(FILE *fp, double *sn)
{

	int i, j;

	int element;
	
	double sd;

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != TETRAHEDRON) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
	    	{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{
		
			sd = sn[elements[element].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != HEXAHEDRON) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{
		
			sd = sn[elements[element].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != PRISM) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{
		
			sd = sn[elements[element].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

}

void PrintAsciiScalarElement(FILE *fp, double *se, int smooth)
{

	int i, j;

	int element;
	
	double sd;

	double *sn;
	int *si;

	sn = calloc(nbnodes, sizeof(double));
	si = calloc(nbnodes, sizeof(int));

	GetNodeScalarElement(sn, si, se);

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != TETRAHEDRON) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{
		
			if (smooth == LOGICAL_FALSE)
				sd = se[element];
			else
				sd = sn[elements[element].node[j]] / si[elements[element].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != HEXAHEDRON) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{
		
			if (smooth == LOGICAL_FALSE)
				sd = se[element];
			else
				sd = sn[elements[element].node[j]] / si[elements[element].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != PRISM) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{
		
			if (smooth == LOGICAL_FALSE)
				sd = se[element];
			else
				sd = sn[elements[element].node[j]] / si[elements[element].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

	free(sn);
	free(si);

}

void PrintAsciiVectorElement(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
							 int iphix, int iphiy, int iphiz, int smooth)
{

	int i, j;

	int element;
	
	double sd;

	double *sn;
	int *si;

	double vx, vy, vz;

	sn = calloc(nbnodes, sizeof(double));
	si = calloc(nbnodes, sizeof(int));
		
	GetNodeVectorElement(xu, xv, xw, xp, xT, xs, sn, si, iphix, iphiy, iphiz);

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != TETRAHEDRON) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
	    }

		for (j = 0; j < elements[element].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		if (iphix == 0) vx = V_GetCmp(xu, element + 1);
		if (iphix == 1) vx = V_GetCmp(xv, element + 1);
		if (iphix == 2) vx = V_GetCmp(xw, element + 1);
		if (iphix == 3) vx = V_GetCmp(xp, element + 1);
		if (iphix == 4) vx = V_GetCmp(xT, element + 1);
		if (iphix == 5) vx = V_GetCmp(xs, element + 1);

		if (iphiy == 0) vy = V_GetCmp(xu, element + 1);
		if (iphiy == 1) vy = V_GetCmp(xv, element + 1);
		if (iphiy == 2) vy = V_GetCmp(xw, element + 1);
		if (iphiy == 3) vy = V_GetCmp(xp, element + 1);
		if (iphiy == 4) vy = V_GetCmp(xT, element + 1);
		if (iphiy == 5) vy = V_GetCmp(xs, element + 1);

		if (iphiz == 0) vz = V_GetCmp(xu, element + 1);
		if (iphiz == 1) vz = V_GetCmp(xv, element + 1);
		if (iphiz == 2) vz = V_GetCmp(xw, element + 1);
		if (iphiz == 3) vz = V_GetCmp(xp, element + 1);
		if (iphiz == 4) vz = V_GetCmp(xT, element + 1);
		if (iphiz == 5) vz = V_GetCmp(xs, element + 1);

		for (j = 0; j < elements[element].nbnodes; j++)
	    {
			fprintf(fp, " %E %E %E", vx, vy, vz);
		}

		fprintf(fp, "\n");  
	}

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != HEXAHEDRON) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
	    }

		for (j = 0; j < elements[element].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		if (iphix == 0) vx = V_GetCmp(xu, element + 1);
		if (iphix == 1) vx = V_GetCmp(xv, element + 1);
		if (iphix == 2) vx = V_GetCmp(xw, element + 1);
		if (iphix == 3) vx = V_GetCmp(xp, element + 1);
		if (iphix == 4) vx = V_GetCmp(xT, element + 1);
		if (iphix == 5) vx = V_GetCmp(xs, element + 1);

		if (iphiy == 0) vy = V_GetCmp(xu, element + 1);
		if (iphiy == 1) vy = V_GetCmp(xv, element + 1);
		if (iphiy == 2) vy = V_GetCmp(xw, element + 1);
		if (iphiy == 3) vy = V_GetCmp(xp, element + 1);
		if (iphiy == 4) vy = V_GetCmp(xT, element + 1);
		if (iphiy == 5) vy = V_GetCmp(xs, element + 1);

		if (iphiz == 0) vz = V_GetCmp(xu, element + 1);
		if (iphiz == 1) vz = V_GetCmp(xv, element + 1);
		if (iphiz == 2) vz = V_GetCmp(xw, element + 1);
		if (iphiz == 3) vz = V_GetCmp(xp, element + 1);
		if (iphiz == 4) vz = V_GetCmp(xT, element + 1);
		if (iphiz == 5) vz = V_GetCmp(xs, element + 1);

		for (j = 0; j < elements[element].nbnodes; j++)
	    {
			fprintf(fp, " %E %E %E", vx, vy, vz);
		}

		fprintf(fp, "\n");  
	}

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type != PRISM) continue;

		for (j = 0; j < elements[element].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].x);
	    }

		for (j = 0; j < elements[element].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].y);
		}

		for (j = 0; j < elements[element].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[elements[element].node[j]].z);
		}

		if (iphix == 0) vx = V_GetCmp(xu, element + 1);
		if (iphix == 1) vx = V_GetCmp(xv, element + 1);
		if (iphix == 2) vx = V_GetCmp(xw, element + 1);
		if (iphix == 3) vx = V_GetCmp(xp, element + 1);
		if (iphix == 4) vx = V_GetCmp(xT, element + 1);
		if (iphix == 5) vx = V_GetCmp(xs, element + 1);

		if (iphiy == 0) vy = V_GetCmp(xu, element + 1);
		if (iphiy == 1) vy = V_GetCmp(xv, element + 1);
		if (iphiy == 2) vy = V_GetCmp(xw, element + 1);
		if (iphiy == 3) vy = V_GetCmp(xp, element + 1);
		if (iphiy == 4) vy = V_GetCmp(xT, element + 1);
		if (iphiy == 5) vy = V_GetCmp(xs, element + 1);

		if (iphiz == 0) vz = V_GetCmp(xu, element + 1);
		if (iphiz == 1) vz = V_GetCmp(xv, element + 1);
		if (iphiz == 2) vz = V_GetCmp(xw, element + 1);
		if (iphiz == 3) vz = V_GetCmp(xp, element + 1);
		if (iphiz == 4) vz = V_GetCmp(xT, element + 1);
		if (iphiz == 5) vz = V_GetCmp(xs, element + 1);

		for (j = 0; j < elements[element].nbnodes; j++)
	    {
			fprintf(fp, " %E %E %E", vx, vy, vz);
		}

		fprintf(fp, "\n");  
	}

	free(sn);
	free(si);

}

void PrintAsciiScalarFace(FILE *fp, double *sf, int smooth)
{

	int i, j;

	int face, pair;
	
	double sd;

	double *sn;
	int *si;

	sn = calloc(nbnodes, sizeof(double));
	si = calloc(nbnodes, sizeof(int));
		
	GetNodeScalarFace(sn, si, sf);

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue; 
		
		if (faces[face].type != TRIANGLE) continue;

		for (j = 0; j < faces[face].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[faces[face].node[j]].x);
	    }

		for (j = 0; j < faces[face].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[faces[face].node[j]].y);
		}

		for (j = 0; j < faces[face].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[faces[face].node[j]].z);
		}

		for (j = 0; j < faces[face].nbnodes; j++)
	    {

			if (smooth == LOGICAL_FALSE)
				sd = sf[face];
			else
				sd = sn[faces[face].node[j]] / si[faces[face].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue; 

		if (faces[face].type != QUADRANGLE) continue;

		for (j = 0; j < faces[face].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[faces[face].node[j]].x);
	    }

		for (j = 0; j < faces[face].nbnodes; j++)
	    {	 
			fprintf(fp, " %E", nodes[faces[face].node[j]].y);
		}

		for (j = 0; j < faces[face].nbnodes; j++)
		{	 
			fprintf(fp, " %E", nodes[faces[face].node[j]].z);
		}

		for (j = 0; j < faces[face].nbnodes; j++)
	    {
			if (smooth == LOGICAL_FALSE)
				sd = sf[face];
			else
				sd = sn[faces[face].node[j]] / si[faces[face].node[j]];

			fprintf(fp, " %E", sd);
		}

		fprintf(fp, "\n");  
	}

	free(sn);
	free(si);

}

void PrintBinaryScalarNode(FILE *fp, double *sn)
{

	int i, j;
	
	int element;
	
	double *tetras;
	double *hexas;
	double *prisms;

	int nt, nh, np;

	tetras = calloc(4 * nbtetras * 4, sizeof(double));
	hexas = calloc(4 * nbhexas * 8, sizeof(double));
	prisms = calloc(4 * nbprisms * 6, sizeof(double));

	nt = 0;
	nh = 0;
	np = 0;

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type == TETRAHEDRON)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].x;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].y;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].z;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				tetras[nt] = sn[elements[element].node[j]];
				nt++;

			}
		}

		if (elements[element].type == HEXAHEDRON)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].x;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].y;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].z;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				hexas[nh] = sn[elements[element].node[j]];
				nh++;

			}
		}

		if (elements[element].type == PRISM)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].x;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].y;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].z;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				prisms[np] = sn[elements[element].node[j]];
				np++;

			}
		}
	}
	
	fwrite(tetras, sizeof(double), nt, fp);
	fwrite(hexas, sizeof(double), nh, fp);
	fwrite(prisms, sizeof(double), np, fp);

	free(tetras);
	free(hexas);
	free(prisms);

}

void PrintBinaryScalarFace(FILE *fp, double *sf, int smooth)
{

	int i, j;

	int face, pair;

	double *tris;
	double *quads;

	int nt, nq;

	double *sn;
	int *si;

	sn = calloc(nbnodes, sizeof(double));
	si = calloc(nbnodes, sizeof(int));
		
	GetNodeScalarFace(sn, si, sf);

	tris = calloc(4 * nbtris * 3, sizeof(double));
	quads = calloc(4 * nbquads * 4, sizeof(double));

	nt = 0;
	nq = 0;

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue;
		
		if (faces[face].type == TRIANGLE)
		{
		
			for (j = 0; j < faces[face].nbnodes; j++)
			{	 
				tris[nt] = nodes[faces[face].node[j]].x;
				nt++;
			}

			for (j = 0; j < faces[face].nbnodes; j++)
			{	 
				tris[nt] = nodes[faces[face].node[j]].y;
				nt++;
			}

			for (j = 0; j < faces[face].nbnodes; j++)
			{	 
				tris[nt] = nodes[faces[face].node[j]].z;
				nt++;
			}

			for (j = 0; j < faces[face].nbnodes; j++)
			{
		
				if (smooth == LOGICAL_FALSE)
					tris[nt] = sf[face];
				else
					tris[nt] = sn[faces[face].node[j]] / si[faces[face].node[j]];

				nt++;

			}
		}

		if (faces[face].type == QUADRANGLE)
		{
		
			for (j = 0; j < faces[face].nbnodes; j++)
			{	 
				quads[nq] = nodes[faces[face].node[j]].x;
				nq++;
			}

			for (j = 0; j < faces[face].nbnodes; j++)
			{	 
				quads[nq] = nodes[faces[face].node[j]].y;
				nq++;
			}

			for (j = 0; j < faces[face].nbnodes; j++)
			{	 
				quads[nq] = nodes[faces[face].node[j]].z;
				nq++;
			}

			for (j = 0; j < faces[face].nbnodes; j++)
			{
		
				if (smooth == LOGICAL_FALSE)
					quads[nq] = sf[face];
				else
					quads[nq] = sn[faces[face].node[j]] / si[faces[face].node[j]];

				nq++;

			}
		}

	}

	fwrite(tris, sizeof(double), nt, fp);
	fwrite(quads, sizeof(double), nq, fp);

	free(tris);
	free(quads);

	free(sn);
	free(si);

}

void PrintBinaryScalarElement(FILE *fp, double *se, int smooth)
{

	int i, j;

	int element;

	double *tetras;
	double *hexas;
	double *prisms;

	int nt, nh, np;

	double *sn;
	int *si;

	sn = calloc(nbnodes, sizeof(double));
	si = calloc(nbnodes, sizeof(int));
		
	GetNodeScalarElement(sn, si, se);

	tetras = calloc(4 * nbtetras * 4, sizeof(double));
	hexas = calloc(4 * nbhexas * 8, sizeof(double));
	prisms = calloc(4 * nbprisms * 6, sizeof(double));

	nt = 0;
	nh = 0;
	np = 0;

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type == TETRAHEDRON)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].x;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].y;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].z;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				if (smooth == LOGICAL_FALSE)
					tetras[nt] = se[element];
				else
					tetras[nt] = sn[elements[element].node[j]] / si[elements[element].node[j]];

				nt++;

			}
		}

		if (elements[element].type == HEXAHEDRON)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].x;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].y;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].z;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				if (smooth == LOGICAL_FALSE)
					hexas[nh] = se[element];
				else
					hexas[nh] = sn[elements[element].node[j]] / si[elements[element].node[j]];

				nh++;

			}
		}

		if (elements[element].type == PRISM)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].x;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].y;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].z;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				if (smooth == LOGICAL_FALSE)
					prisms[np] = se[element];
				else
					prisms[np] = sn[elements[element].node[j]] / si[elements[element].node[j]];

				np++;

			}
		}
	}

	fwrite(tetras, sizeof(double), nt, fp);
	fwrite(hexas, sizeof(double), nh, fp);
	fwrite(prisms, sizeof(double), np, fp);

	free(tetras);
	free(hexas);
	free(prisms);

	free(sn);
	free(si);

}

void PrintBinaryVectorElement(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
							  int iphix, int iphiy, int iphiz, int smooth)
{

	int i, j;

	int element;

	double vx, vy, vz;

	double *tetras;
	double *hexas;
	double *prisms;

	int nt, nh, np;

	tetras = calloc(6 * nbtetras * 4, sizeof(double));
	hexas = calloc(6 * nbhexas * 8, sizeof(double));
	prisms = calloc(6 * nbprisms * 6, sizeof(double));

	nt = 0;
	nh = 0;
	np = 0;

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (elements[element].type == TETRAHEDRON)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].x;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].y;
				nt++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				tetras[nt] = nodes[elements[element].node[j]].z;
				nt++;
			}

			if (iphix == 0) vx = V_GetCmp(xu, element + 1);
			if (iphix == 1) vx = V_GetCmp(xv, element + 1);
			if (iphix == 2) vx = V_GetCmp(xw, element + 1);
			if (iphix == 3) vx = V_GetCmp(xp, element + 1);
			if (iphix == 4) vx = V_GetCmp(xT, element + 1);
			if (iphix == 5) vx = V_GetCmp(xs, element + 1);

			if (iphiy == 0) vy = V_GetCmp(xu, element + 1);
			if (iphiy == 1) vy = V_GetCmp(xv, element + 1);
			if (iphiy == 2) vy = V_GetCmp(xw, element + 1);
			if (iphiy == 3) vy = V_GetCmp(xp, element + 1);
			if (iphiy == 4) vy = V_GetCmp(xT, element + 1);
			if (iphiy == 5) vy = V_GetCmp(xs, element + 1);

			if (iphiz == 0) vz = V_GetCmp(xu, element + 1);
			if (iphiz == 1) vz = V_GetCmp(xv, element + 1);
			if (iphiz == 2) vz = V_GetCmp(xw, element + 1);
			if (iphiz == 3) vz = V_GetCmp(xp, element + 1);
			if (iphiz == 4) vz = V_GetCmp(xT, element + 1);
			if (iphiz == 5) vz = V_GetCmp(xs, element + 1);

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				tetras[nt] = vx;
				nt++;

				tetras[nt] = vy;
				nt++;

				tetras[nt] = vz;
				nt++;

			}
		}

		if (elements[element].type == HEXAHEDRON)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].x;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].y;
				nh++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				hexas[nh] = nodes[elements[element].node[j]].z;
				nh++;
			}

			if (iphix == 0) vx = V_GetCmp(xu, element + 1);
			if (iphix == 1) vx = V_GetCmp(xv, element + 1);
			if (iphix == 2) vx = V_GetCmp(xw, element + 1);
			if (iphix == 3) vx = V_GetCmp(xp, element + 1);
			if (iphix == 4) vx = V_GetCmp(xT, element + 1);
			if (iphix == 5) vx = V_GetCmp(xs, element + 1);

			if (iphiy == 0) vy = V_GetCmp(xu, element + 1);
			if (iphiy == 1) vy = V_GetCmp(xv, element + 1);
			if (iphiy == 2) vy = V_GetCmp(xw, element + 1);
			if (iphiy == 3) vy = V_GetCmp(xp, element + 1);
			if (iphiy == 4) vy = V_GetCmp(xT, element + 1);
			if (iphiy == 5) vy = V_GetCmp(xs, element + 1);

			if (iphiz == 0) vz = V_GetCmp(xu, element + 1);
			if (iphiz == 1) vz = V_GetCmp(xv, element + 1);
			if (iphiz == 2) vz = V_GetCmp(xw, element + 1);
			if (iphiz == 3) vz = V_GetCmp(xp, element + 1);
			if (iphiz == 4) vz = V_GetCmp(xT, element + 1);
			if (iphiz == 5) vz = V_GetCmp(xs, element + 1);

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				hexas[nh] = vx;
				nh++;

				hexas[nh] = vy;
				nh++;

				hexas[nh] = vz;
				nh++;

			}
		}

		if (elements[element].type == PRISM)
		{
		
			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].x;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].y;
				np++;
			}

			for (j = 0; j < elements[element].nbnodes; j++)
			{	 
				prisms[np] = nodes[elements[element].node[j]].z;
				np++;
			}

			if (iphix == 0) vx = V_GetCmp(xu, element + 1);
			if (iphix == 1) vx = V_GetCmp(xv, element + 1);
			if (iphix == 2) vx = V_GetCmp(xw, element + 1);
			if (iphix == 3) vx = V_GetCmp(xp, element + 1);
			if (iphix == 4) vx = V_GetCmp(xT, element + 1);
			if (iphix == 5) vx = V_GetCmp(xs, element + 1);

			if (iphiy == 0) vy = V_GetCmp(xu, element + 1);
			if (iphiy == 1) vy = V_GetCmp(xv, element + 1);
			if (iphiy == 2) vy = V_GetCmp(xw, element + 1);
			if (iphiy == 3) vy = V_GetCmp(xp, element + 1);
			if (iphiy == 4) vy = V_GetCmp(xT, element + 1);
			if (iphiy == 5) vy = V_GetCmp(xs, element + 1);

			if (iphiz == 0) vz = V_GetCmp(xu, element + 1);
			if (iphiz == 1) vz = V_GetCmp(xv, element + 1);
			if (iphiz == 2) vz = V_GetCmp(xw, element + 1);
			if (iphiz == 3) vz = V_GetCmp(xp, element + 1);
			if (iphiz == 4) vz = V_GetCmp(xT, element + 1);
			if (iphiz == 5) vz = V_GetCmp(xs, element + 1);

			for (j = 0; j < elements[element].nbnodes; j++)
			{
		
				prisms[np] = vx;
				np++;

				prisms[np] = vy;
				np++;

				prisms[np] = vz;
				np++;

			}
		}
	}

	fwrite(tetras, sizeof(double), nt, fp);
	fwrite(hexas, sizeof(double), nh, fp);
	fwrite(prisms, sizeof(double), np, fp);

	free(tetras);
	free(hexas);
	free(prisms);

}

void WriteAsciiScalarElement(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
							 char *label, int iphi, double time, int smooth)
{

	int i;

	int element;

	double vs;

	double *se;

	fprintf(fp, "$View\n");

	se = calloc(nbelements, sizeof(double));

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (iphi == 0) vs = V_GetCmp(xu, element + 1);
		if (iphi == 1) vs = V_GetCmp(xv, element + 1);
		if (iphi == 2) vs = V_GetCmp(xw, element + 1);
		if (iphi == 3) vs = V_GetCmp(xp, element + 1);
		if (iphi == 4) vs = V_GetCmp(xT, element + 1);
		if (iphi == 5) vs = V_GetCmp(xs, element + 1);

		se[element] = vs;

	}

	PrintAsciiHeaderScalarElement(fp, label);
	PrintAsciiTime(fp, time);
	PrintAsciiScalarElement(fp, se, smooth);

	fprintf(fp, "$EndView\n");

	free(se);

}

void WriteAsciiVectorElement(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
							 char *label, int iphix, int iphiy, int iphiz, double time, int smooth)
{

	fprintf(fp, "$View\n");
	
	PrintAsciiHeaderVectorElement(fp, label);
	PrintAsciiTime(fp, time);
	PrintAsciiVectorElement(fp, xu, xv, xw, xp, xT, xs, iphix, iphiy, iphiz, smooth);

	fprintf(fp, "$EndView\n");
	
}

void WriteAsciiStreamFunction(FILE *fp, Vector *uf, char *label, double time)
{

	double *sn;
			
	fprintf(fp, "$View\n");
	
	sn = calloc(nbnodes, sizeof(double));
		
	CalculateStreamFunction(uf, sn);
	
	PrintAsciiHeaderScalarElement(fp, label);
	PrintAsciiTime(fp, time);
	PrintAsciiScalarNode(fp, sn);
	
	fprintf(fp, "$EndView\n");
	
	free(sn);

}


void WriteAsciiVorticity(FILE *fp, Vector *xu, Vector *xv, Vector *xw, 
				    Vector *xuf, Vector *xvf, Vector *xwf,
			            char *label, int axis, double time, int smooth)
{

	int i;

	int element;

	double vs;

	double *se;

	msh_vector gradu, gradv, gradw;
	
	fprintf(fp, "$View\n");

	se = calloc(nbelements, sizeof(double));
			
	for (i = 0; i < nbelements; i++)
	{

		element = i;
		
		gradu = Gradient(xu, xuf, LOGICAL_FALSE, element);
		gradv = Gradient(xv, xvf, LOGICAL_FALSE, element);
		gradw = Gradient(xw, xwf, LOGICAL_FALSE, element);
		
		if (axis == 0)
		{
			vs = gradw.y - gradv.z;
		}

		if (axis == 1)
		{
			vs = gradu.z - gradw.x;
		}

		if (axis == 2)
		{
			vs = gradv.x - gradu.y;
		}
				
		se[element] = vs;

	}

	PrintAsciiHeaderScalarElement(fp, label);
	PrintAsciiTime(fp, time);
	PrintAsciiScalarElement(fp, se, smooth);

	fprintf(fp, "\n$EndView\n");

	free(se);

}

void WriteAsciiScalarFace(FILE *fp, Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf, 
			  char *label, int iphi, double time, int smooth)
{

	int i;

	int face;

	double vs;

	double *sf;

	fprintf(fp, "$View\n");

	sf = calloc(nbfaces, sizeof(double));

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		if (iphi == 0) vs = V_GetCmp(xuf, face + 1);
		if (iphi == 1) vs = V_GetCmp(xvf, face + 1);
		if (iphi == 2) vs = V_GetCmp(xwf, face + 1);
		if (iphi == 3) vs = V_GetCmp(xpf, face + 1);
		if (iphi == 4) vs = V_GetCmp(xTf, face + 1);
		if (iphi == 5) vs = V_GetCmp(xsf, face + 1);

		sf[face] = vs;

	}

	PrintAsciiHeaderScalarFace(fp, label);
	PrintAsciiTime(fp, time);
	PrintAsciiScalarFace(fp, sf, smooth);

	fprintf(fp, "\n$EndView\n");

	free(sf);

}

void WriteAsciiVectorMagnitudeFace(FILE *fp, Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf, 
								   char *label, int iphix, int iphiy, int iphiz, double time, int smooth)
{

	int i;

	int face, pair;

	double vx, vy, vz;

	double *sf;

	sf = calloc(nbfaces, sizeof(double));

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue;

		if (iphix == 0) vx = V_GetCmp(xuf, face + 1);
		if (iphix == 1) vx = V_GetCmp(xvf, face + 1);
		if (iphix == 2) vx = V_GetCmp(xwf, face + 1);
		if (iphix == 3) vx = V_GetCmp(xpf, face + 1);
		if (iphix == 4) vx = V_GetCmp(xTf, face + 1);
		if (iphix == 5) vx = V_GetCmp(xsf, face + 1);

		if (iphiy == 0) vy = V_GetCmp(xuf, face + 1);
		if (iphiy == 1) vy = V_GetCmp(xvf, face + 1);
		if (iphiy == 2) vy = V_GetCmp(xwf, face + 1);
		if (iphiy == 3) vy = V_GetCmp(xpf, face + 1);
		if (iphiy == 4) vy = V_GetCmp(xTf, face + 1);
		if (iphiy == 5) vy = V_GetCmp(xsf, face + 1);

		if (iphiz == 0) vz = V_GetCmp(xuf, face + 1);
		if (iphiz == 1) vz = V_GetCmp(xvf, face + 1);
		if (iphiz == 2) vz = V_GetCmp(xwf, face + 1);
		if (iphiz == 3) vz = V_GetCmp(xpf, face + 1);
		if (iphiz == 4) vz = V_GetCmp(xTf, face + 1);
		if (iphiz == 5) vz = V_GetCmp(xsf, face + 1);

		sf[face] = sqrt(vx * vx + vy * vy + vz * vz);

	}

	fprintf(fp, "$View\n");

	PrintAsciiHeaderScalarFace(fp, label);
	PrintAsciiTime(fp, time);
	PrintAsciiScalarFace(fp, sf, smooth);

	fprintf(fp, "\n$EndView\n");

	free(sf);

}
void WriteBinaryScalarFace(FILE *fp, Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf, 
						   char *label, int iphi, double time, int smooth)
{

	int i;

	int face, pair;

	double vs;

	double *sf;

	fprintf(fp, "$View\n");

	sf = calloc(nbfaces, sizeof(double));

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue;

		if (iphi == 0) vs = V_GetCmp(xuf, face + 1);
		if (iphi == 1) vs = V_GetCmp(xvf, face + 1);
		if (iphi == 2) vs = V_GetCmp(xwf, face + 1);
		if (iphi == 3) vs = V_GetCmp(xpf, face + 1);
		if (iphi == 4) vs = V_GetCmp(xTf, face + 1);
		if (iphi == 5) vs = V_GetCmp(xsf, face + 1);

		sf[face] = vs;

	}

	PrintBinaryHeaderScalarFace(fp, label);
	PrintBinaryTime(fp, time);
	PrintBinaryScalarFace(fp, sf, smooth);

	fprintf(fp, "\n$EndView\n");

	free(sf);

}

void WriteBinaryVectorMagnitudeFace(FILE *fp, Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf, 
									char *label, int iphix, int iphiy, int iphiz, double time, int smooth)
{

	int i;

	int face, pair;

	double vx, vy, vz;

	double *sf;

	sf = calloc(nbfaces, sizeof(double));

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		pair = faces[face].pair;

		if (pair != -1) continue;

		if (iphix == 0) vx = V_GetCmp(xuf, face + 1);
		if (iphix == 1) vx = V_GetCmp(xvf, face + 1);
		if (iphix == 2) vx = V_GetCmp(xwf, face + 1);
		if (iphix == 3) vx = V_GetCmp(xpf, face + 1);
		if (iphix == 4) vx = V_GetCmp(xTf, face + 1);
		if (iphix == 5) vx = V_GetCmp(xsf, face + 1);

		if (iphiy == 0) vy = V_GetCmp(xuf, face + 1);
		if (iphiy == 1) vy = V_GetCmp(xvf, face + 1);
		if (iphiy == 2) vy = V_GetCmp(xwf, face + 1);
		if (iphiy == 3) vy = V_GetCmp(xpf, face + 1);
		if (iphiy == 4) vy = V_GetCmp(xTf, face + 1);
		if (iphiy == 5) vy = V_GetCmp(xsf, face + 1);

		if (iphiz == 0) vz = V_GetCmp(xuf, face + 1);
		if (iphiz == 1) vz = V_GetCmp(xvf, face + 1);
		if (iphiz == 2) vz = V_GetCmp(xwf, face + 1);
		if (iphiz == 3) vz = V_GetCmp(xpf, face + 1);
		if (iphiz == 4) vz = V_GetCmp(xTf, face + 1);
		if (iphiz == 5) vz = V_GetCmp(xsf, face + 1);

		sf[face] = sqrt(vx * vx + vy * vy + vz * vz); 

	}

	fprintf(fp, "$View\n");

	PrintBinaryHeaderScalarFace(fp, label);
	PrintBinaryTime(fp, time);
	PrintBinaryScalarFace(fp, sf, smooth);

	fprintf(fp, "\n$EndView\n");

	free(sf);

}

void WriteBinaryScalarElement(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
							  char *label, int iphi, double time, int smooth)
{

	int i;

	int element;

	double vs;

	double *se;

	fprintf(fp, "$View\n");

	se = calloc(nbelements, sizeof(double));

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		if (iphi == 0) vs = V_GetCmp(xu, element + 1);
		if (iphi == 1) vs = V_GetCmp(xv, element + 1);
		if (iphi == 2) vs = V_GetCmp(xw, element + 1);
		if (iphi == 3) vs = V_GetCmp(xp, element + 1);
		if (iphi == 4) vs = V_GetCmp(xT, element + 1);
		if (iphi == 5) vs = V_GetCmp(xs, element + 1);

		se[element] = vs;

	}

	PrintBinaryHeaderScalarElement(fp, label);
	PrintBinaryTime(fp, time);
	PrintBinaryScalarElement(fp, se, smooth);

	fprintf(fp, "\n$EndView\n");

	free(se);

}

void WriteBinaryVectorElement(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
							  char *label, int iphix, int iphiy, int iphiz, double time, int smooth)
{

	fprintf(fp, "$View\n");

	PrintBinaryHeaderVectorElement(fp, label);
	PrintBinaryTime(fp, time);
	PrintBinaryVectorElement(fp, xu, xv, xw, xp, xT, xs, iphix, iphiy, iphiz, smooth);

	fprintf(fp, "\n$EndView\n");

}

void WriteBinaryStreamFunction(FILE *fp, Vector *uf, char *label, double time)
{

	double *sn;
			
	fprintf(fp, "$View\n");
	
	sn = calloc(nbnodes, sizeof(double));
		
	CalculateStreamFunction(uf, sn);
	
	PrintBinaryHeaderScalarElement(fp, label);
	PrintBinaryTime(fp, time);
	PrintBinaryScalarNode(fp, sn);
	
	fprintf(fp, "\n$EndView\n");
	
	free(sn);
}

void WriteBinaryVorticity(FILE *fp, Vector *xu, Vector *xv, Vector *xw, 
			  Vector *xuf, Vector *xvf, Vector *xwf,
			  char *label, int axis, double time, int smooth)
{

	int i;

	int element;

	double omega;

	double *se;

	msh_vector gradu, gradv, gradw;
	
	fprintf(fp, "$View\n");

	se = calloc(nbelements, sizeof(double));
			
	for (i = 0; i < nbelements; i++)
	{

		element = i;
		
		gradu = Gradient(xu, xuf, LOGICAL_FALSE, element);
		gradv = Gradient(xv, xvf, LOGICAL_FALSE, element);
		gradw = Gradient(xw, xwf, LOGICAL_FALSE, element);
		
		if (axis == 0)
		{
			omega = gradw.y - gradv.z;
		}

		if (axis == 1)
		{
			omega = gradu.z - gradw.x;
		}

		if (axis == 2)
		{
			omega = gradv.x - gradu.y;
		}
						
		se[element] = omega;

	}
	
	PrintBinaryHeaderScalarElement(fp, label);
	PrintBinaryTime(fp, time);
	PrintBinaryScalarElement(fp, se, smooth);

	fprintf(fp, "\n$EndView\n");

	free(se);

}

void WriteResults(FILE *fp, Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs,
			    Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf,
			    Vector *uf, 
			    int onFace, int inElement, double time)
{

	if (onFace == LOGICAL_TRUE)
	{
	
		if (parameter.wbinary == LOGICAL_TRUE)
		{
			fprintf(fp, "$PostFormat\n");
			fprintf(fp, "%g %d %d\n", 1.3, 1, sizeof(double));
			fprintf(fp, "$EndPostFormat\n");

			// Write velocity-U to file 
			if (parameter.fsav[iu] == LOGICAL_TRUE)
				WriteBinaryScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "VelocityU-Face_[m/s]", iu, 0.0, LOGICAL_FALSE);

			// Write velocity-V to file 
			if (parameter.fsav[iv] == LOGICAL_TRUE)
				WriteBinaryScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "VelocityV-Face_[m/s]", iv, 0.0, LOGICAL_FALSE);

			// Write velocity-W to file 
			if (parameter.fsav[iw] == LOGICAL_TRUE)
				WriteBinaryScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "VelocityW-Face_[m/s]", iw, 0.0, LOGICAL_FALSE);

			// Write pressure to file
			if (parameter.fsav[ip] == LOGICAL_TRUE)
				WriteBinaryScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Pressure-Face_[kg/(m*s*s)]", ip, 0.0, LOGICAL_FALSE);

			// Write indicator function to file
			if (parameter.fsav[iT] == LOGICAL_TRUE)
				WriteBinaryScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Temperature-Face_[K]", iT, 0.0, LOGICAL_FALSE);

			// Write indicator function to file
			if (parameter.fsav[is] == LOGICAL_TRUE)
				WriteBinaryScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Gamma-Face_[]", is, 0.0, LOGICAL_FALSE);

			// Write velocity to file 
			if (parameter.fvec == LOGICAL_TRUE)
				WriteBinaryVectorMagnitudeFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Velocity-Face_[m/s]", iu, iv, iw, 0.0, LOGICAL_FALSE);

		}
		else
		{

			fprintf(fp, "$PostFormat\n");
			fprintf(fp, "%g %d %d\n", 1.2, 0, sizeof(double));
			fprintf(fp, "$EndPostFormat\n");

			// Write velocity-U to file 
			if (parameter.fsav[iu] == LOGICAL_TRUE)
				WriteAsciiScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "VelocityU-Face_[m/s]", iu, 0.0, LOGICAL_FALSE);

			// Write velocity-V to file 
			if (parameter.fsav[iv] == LOGICAL_TRUE)
				WriteAsciiScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "VelocityV-Face_[m/s]", iv, 0.0, LOGICAL_FALSE);

			// Write velocity-W to file 
			if (parameter.fsav[iw] == LOGICAL_TRUE)
				WriteAsciiScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "VelocityW-Face_[m/s]", iw, 0.0, LOGICAL_FALSE);

			// Write pressure to file
			if (parameter.fsav[ip] == LOGICAL_TRUE)
				WriteAsciiScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Pressure-Face_[kg/(m*s*s)]", ip, 0.0, LOGICAL_FALSE);

			// Write indicator function to file
			if (parameter.fsav[iT] == LOGICAL_TRUE)
				WriteAsciiScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Temperature-Face_[K]", iT, 0.0, LOGICAL_FALSE);

			// Write indicator function to file
			if (parameter.fsav[is] == LOGICAL_TRUE)
				WriteAsciiScalarFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Gamma-Face_[]", is, 0.0, LOGICAL_FALSE);

			// Write velocity to file 
			if (parameter.fvec == LOGICAL_TRUE)
				WriteBinaryVectorMagnitudeFace(fp, xuf, xvf, xwf, xpf, xTf, xsf, "Velocity-Face_[m/s]", iu, iv, iw, 0.0, LOGICAL_FALSE);

		}
	}

	if (inElement == LOGICAL_TRUE)
	{

		if (parameter.wbinary == LOGICAL_TRUE)
		{
			
			// Write velocity-U to file 
			if (parameter.csav[iu] == LOGICAL_TRUE)
				WriteBinaryScalarElement(fp, xu, xv, xw, xp, xT, xs, "VelocityU-Cell_[m/s]", iu, time, LOGICAL_TRUE);

			// Write velocity-V to file 
			if (parameter.csav[iv] == LOGICAL_TRUE)
				WriteBinaryScalarElement(fp, xu, xv, xw, xp, xT, xs, "VelocityV-Cell_[m/s]", iv, time, LOGICAL_TRUE);

			// Write velocity-W to file 
			if (parameter.csav[iw] == LOGICAL_TRUE)
				WriteBinaryScalarElement(fp, xu, xv, xw, xp, xT, xs, "VelocityW-Cell_[m/s]", iw, time, LOGICAL_TRUE);
			
			// Write pressure to file
			if (parameter.csav[ip] == LOGICAL_TRUE)
				WriteBinaryScalarElement(fp, xu, xv, xw, xp, xT, xs, "Pressure-Cell_[kg/(m*s*s)]", ip, time, LOGICAL_TRUE);		

			// Write indicator function to file
			if (parameter.csav[iT] == LOGICAL_TRUE)
				WriteBinaryScalarElement(fp, xu, xv, xw, xp, xT, xs, "Temperature-Cell_[K]", iT, time, LOGICAL_TRUE);		

			// Write indicator function to file
			if (parameter.csav[is] == LOGICAL_TRUE)
				WriteBinaryScalarElement(fp, xu, xv, xw, xp, xT, xs, "Gamma-Cell_[]", is, time, LOGICAL_TRUE);
				
			// Write velocity to file 
			if (parameter.cvec == LOGICAL_TRUE)
				WriteBinaryVectorElement(fp, xu, xv, xw, xp, xT, xs, "Velocity-Cell_[m/s]", iu, iv, iw, time, LOGICAL_TRUE);
			
			// Write vorticity to file			
			if (parameter.vortex[0] == LOGICAL_TRUE)
				WriteBinaryVorticity(fp, xu, xv, xw, xuf, xvf, xwf, "VorticityU-Cell_[1/s]", 0, time, LOGICAL_TRUE);

			if (parameter.vortex[1] == LOGICAL_TRUE)
				WriteBinaryVorticity(fp, xu, xv, xw, xuf, xvf, xwf, "VorticityV-Cell_[1/s]", 1, time, LOGICAL_TRUE);
								
			if (parameter.vortex[2] == LOGICAL_TRUE)
				WriteBinaryVorticity(fp, xu, xv, xw, xuf, xvf, xwf, "VorticityW-Cell_[1/s]", 2, time, LOGICAL_TRUE);

			if (parameter.streamf == LOGICAL_TRUE)
				WriteBinaryStreamFunction(fp, uf, "StreamFunction-Cell[]", time);
								
		}
		else
		{

			// Write velocity-U to file 
			if (parameter.csav[iu] == LOGICAL_TRUE)
				WriteAsciiScalarElement(fp, xu, xv, xw, xp, xT, xs, "VelocityU-Cell_[m/s]", iu, time, LOGICAL_TRUE);

			// Write velocity-V to file 
			if (parameter.csav[iv] == LOGICAL_TRUE)
				WriteAsciiScalarElement(fp, xu, xv, xw, xp, xT, xs, "VelocityV-Cell_[m/s]", iv, time, LOGICAL_TRUE);

			// Write velocity-W to file 			
			if (parameter.csav[iw] == LOGICAL_TRUE)
				WriteAsciiScalarElement(fp, xu, xv, xw, xp, xT, xs, "VelocityW-Cell_[m/s]", iw, time, LOGICAL_TRUE);
						
			// Write pressure to file
			if (parameter.csav[ip] == LOGICAL_TRUE)
				WriteAsciiScalarElement(fp, xu, xv, xw, xp, xT, xs, "Pressure-Cell_[kg/(m*s*s)]", ip, time, LOGICAL_TRUE);		

			// Write indicator function to file
			if (parameter.csav[iT] == LOGICAL_TRUE)
				WriteAsciiScalarElement(fp, xu, xv, xw, xp, xT, xs, "Temperature-Cell_[K]", iT, time, LOGICAL_TRUE);		

			// Write indicator function to file
			if (parameter.csav[is] == LOGICAL_TRUE)
				WriteAsciiScalarElement(fp, xu, xv, xw, xp, xT, xs, "Gamma-Cell_[]", is, time, LOGICAL_TRUE);		
				
			// Write velocity to file 
			if (parameter.cvec == LOGICAL_TRUE)
				WriteAsciiVectorElement(fp, xu, xv, xw, xp, xT, xs, "Velocity-Cell_[m/s]", iu, iv, iw, time, LOGICAL_TRUE);
				
			// Write vorticity to file			
			if (parameter.vortex[0] == LOGICAL_TRUE)
				WriteAsciiVorticity(fp, xu, xv, xw, xuf, xvf, xwf, "VorticityU-Cell_[1/s]", 0, time, LOGICAL_TRUE);
								
			if (parameter.vortex[1] == LOGICAL_TRUE)
				WriteAsciiVorticity(fp, xu, xv, xw, xuf, xvf, xwf, "VorticityV-Cell_[1/s]", 1, time, LOGICAL_TRUE);

			if (parameter.vortex[2] == LOGICAL_TRUE)
				WriteAsciiVorticity(fp, xu, xv, xw, xuf, xvf, xwf, "VorticityW-Cell_[1/s]", 2, time, LOGICAL_TRUE);

			if (parameter.streamf == LOGICAL_TRUE)
				WriteAsciiStreamFunction(fp, uf, "StreamFunction-Cell[]", time);
			
		}

	}

}

void WriteResidual(FILE *fp, int iter, double *res)
{

	int i;

	fprintf(fp, "%d", iter);

	for (i = 0; i < nphi; i++)
	{
		fprintf(fp, " \t %f", res[i]);
	}

	fprintf(fp, "\n");

}
