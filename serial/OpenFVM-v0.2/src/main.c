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
#include <malloc.h>
#include <math.h>
#include <time.h>

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"
 
#include "gradient.h"
#include "globals.h"
#include "ioutils.h"
#include "geocalc.h"
#include "parser.h"

#include "laspack/operats.h"
#include "laspack/itersolv.h"
#include "laspack/rtc.h"
#include "laspack/errhandl.h"

int verbose = LOGICAL_FALSE;

void WriteMatrix(QMatrix *A, Vector *b, int ni)
{

	FILE *fp;

	int i, j, nj;

	fp = fopen("matrix.out", "w");
	
	fprintf(fp, "%d %d\n", ni, ni);

	for (i = 0; i < ni; i++)
	{

		nj = Q_GetLen(A, i + 1);

		for (j = 0; j < nj; j++)
		{
			fprintf(fp, "%d %d %+E\n", i + 1, Q_GetPos(A, i + 1, j), Q_GetVal(A, i + 1, j));
		}
	}

	for (i = 0; i < ni; i++)
	{
		fprintf(fp, "%d %+E\n", i + 1, V_GetCmp(b, i + 1));
	}

	fclose(fp);

}

int CheckIfDiagonalMatrix(QMatrix *A, int n)
{

	int i, j, nj;

	int cond1, cond2;

	double sum_ap_parcial, sum_an_parcial;
	double sum_ap_total, sum_an_total;

	sum_ap_total = 0.0;
	sum_an_total = 0.0;

	cond1 = LOGICAL_FALSE;
	cond2 = LOGICAL_FALSE;

	for (i = 0; i < n; i++)
	{

		nj = Q_GetLen(A, i + 1);

		if (nj > 0) sum_ap_parcial = ABS(Q_GetVal(A, i + 1, 0));

		sum_an_parcial = 0.0;

		for (j = 1; j < nj; j++)
		{
			sum_an_parcial += ABS(Q_GetVal(A, i + 1, j));
		}

		if (sum_ap_parcial >= sum_an_parcial) cond1 = LOGICAL_TRUE;

		sum_ap_total += sum_ap_parcial;

		sum_an_total += sum_an_parcial;

	}

	if (sum_ap_total >= sum_an_total) 
		cond2 = LOGICAL_TRUE;

	if (!cond1 && !cond2)
		return LOGICAL_FALSE;
	else
		return LOGICAL_TRUE;

}

void CorrectFaceP(Vector *xp, Vector *xpf)
{

	int i;

	int face, pair;
	int element, neighbor;

	double apj;

	double dNf, dPf;
	double lambda;

	double ppl;
	double pnl;
		
	msh_vector gradpp;
	msh_vector gradpn;
			
	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		gradpp = Gradient(xp, xpf, LOGICAL_TRUE, element);

		if (pair != -1)
		{

			neighbor = faces[pair].element;

			dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[face].cface)); 
			dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 
	
			lambda = dPf / (dPf + dNf);
										
			V_SetCmp(xpf, face + 1, V_GetCmp(xp, neighbor + 1) * lambda + V_GetCmp(xp, element + 1) * (1.0 - lambda));

		}
		else
		{
			
			ppl = V_GetCmp(xp, element + 1) + GeoDotVectorVector(gradpp, GeoSubVectorVector(faces[face].rpl, elements[element].celement)); 
			
			if (faces[face].bc == INLET || 
				faces[face].bc == MOVINGWALL || 
				faces[face].bc == WALL ||
				faces[face].bc == ADIABATICWALL ||
				faces[face].bc == SURFACE)
			{

				V_SetCmp(xpf, face + 1, ppl);
								
			}

		}

	}

}

void CorrectFaceUVW(Vector *ap, 
				    Vector *xu, Vector *xv, Vector *xw, Vector *xp, 
				    Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf,
				    Vector *uf)
{

	int i;

	int face, pair;
	int element, neighbor;

	double apj;

	double dNf, dPf;
	double lambda;

	double ppl;
	double pnl;
		
	msh_vector gradpp;
	msh_vector gradpn;
			
	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		gradpp = Gradient(xp, xpf, LOGICAL_TRUE, element);

		if (pair != -1)
		{

			neighbor = faces[pair].element;

			dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[face].cface)); 
			dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 
	
			lambda = dPf / (dPf + dNf);

			apj = V_GetCmp(ap, neighbor + 1) * lambda + V_GetCmp(ap, element + 1) * (1.0 - lambda);
				
			gradpn = Gradient(xp, xpf, LOGICAL_TRUE, neighbor);
							
			ppl = V_GetCmp(xp, element + 1) + GeoDotVectorVector(gradpp, GeoSubVectorVector(faces[face].rpl, elements[element].celement)); 
			pnl = V_GetCmp(xp, neighbor + 1) + GeoDotVectorVector(gradpn, GeoSubVectorVector(faces[face].rnl, elements[neighbor].celement)); 
			
			V_SetCmp(uf, face + 1, V_GetCmp(uf, face + 1) - 1 / (apj * faces[face].dj) * (pnl - ppl));

			V_SetCmp(xuf, face + 1, V_GetCmp(xu, neighbor + 1) * lambda + V_GetCmp(xu, element + 1) * (1.0 - lambda));
			V_SetCmp(xvf, face + 1, V_GetCmp(xv, neighbor + 1) * lambda + V_GetCmp(xv, element + 1) * (1.0 - lambda));
			V_SetCmp(xwf, face + 1, V_GetCmp(xw, neighbor + 1) * lambda + V_GetCmp(xw, element + 1) * (1.0 - lambda));

			V_SetCmp(xpf, face + 1, V_GetCmp(xp, neighbor + 1) * lambda + V_GetCmp(xp, element + 1) * (1.0 - lambda));

		}
		else
		{

			apj = V_GetCmp(ap, element + 1);
	
			ppl = V_GetCmp(xp, element + 1) + GeoDotVectorVector(gradpp, GeoSubVectorVector(faces[face].rpl, elements[element].celement)); 

			if (faces[face].bc == OUTLET)
			{

				// velocity gradient = 0
				// specified pressure

				V_SetCmp(uf, face + 1, V_GetCmp(uf, face + 1) - 1 / (apj * faces[face].dj) * (V_GetCmp(xpf, face + 1) - ppl));

				V_SetCmp(xuf, face + 1, V_GetCmp(xu, element + 1)); 
				V_SetCmp(xvf, face + 1, V_GetCmp(xv, element + 1)); 
				V_SetCmp(xwf, face + 1, V_GetCmp(xw, element + 1)); 
				
				
			}

			if (faces[face].bc == PRESSUREINLET)
			{

				// velocity gradient = 0
				// specified pressure

  				V_SetCmp(uf, face + 1, V_GetCmp(uf, face + 1) - 1 / (apj * faces[face].dj) * (V_GetCmp(xpf, face + 1) - ppl));

				V_SetCmp(xuf, face + 1, V_GetCmp(xu, element + 1)); 
				V_SetCmp(xvf, face + 1, V_GetCmp(xv, element + 1)); 
				V_SetCmp(xwf, face + 1, V_GetCmp(xw, element + 1)); 
			
			}
			
			if (faces[face].bc == INLET || 
				faces[face].bc == MOVINGWALL || 
				faces[face].bc == WALL ||
				faces[face].bc == ADIABATICWALL ||
				faces[face].bc == SURFACE)
			{

				// pressure gradient = 0
				// specified velocity

				V_SetCmp(uf, face + 1, V_GetCmp(xuf, face + 1) * faces[face].n.x + 
									   V_GetCmp(xvf, face + 1) * faces[face].n.y + 
									   V_GetCmp(xwf, face + 1) * faces[face].n.z);
				
				V_SetCmp(xpf, face + 1, ppl);
								
			}

			if (faces[face].bc == SLIP)
			{

				// pressure gradient = 0
				// velocity gradient = 0

				V_SetCmp(xuf, face + 1, V_GetCmp(xu, element + 1)); 
				V_SetCmp(xvf, face + 1, V_GetCmp(xv, element + 1)); 
				V_SetCmp(xwf, face + 1, V_GetCmp(xw, element + 1)); 

				V_SetCmp(uf, face + 1, 0.0);
				
				V_SetCmp(xpf, face + 1, ppl);
								
			}

		}

	}

}

void CorrectFaceT(Vector *xT, Vector *xTf)
{

	int i;

	int face, pair;
	int element, neighbor;

	double apj;

	double dNf, dPf;
	double lambda;

	double Tpl;
	double Tnl;
		
	msh_vector gradTp;
	msh_vector gradTn;
			
	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		gradTp = Gradient(xT, xTf, LOGICAL_TRUE, element);

		if (pair != -1)
		{

			neighbor = faces[pair].element;

			dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[face].cface)); 
			dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 
	
			lambda = dPf / (dPf + dNf);
										
			V_SetCmp(xTf, face + 1, V_GetCmp(xT, neighbor + 1) * lambda + V_GetCmp(xT, element + 1) * (1.0 - lambda));

		}
		else
		{
			
			Tpl = V_GetCmp(xT, element + 1) + GeoDotVectorVector(gradTp, GeoSubVectorVector(faces[face].rpl, elements[element].celement)); 
			
			if (faces[face].bc == ADIABATICWALL)
			{

				V_SetCmp(xTf, face + 1, Tpl);
								
			}

		}

	}

}

void CorrectFaceS(Vector *xs, Vector *xsf, Vector *xsm, Vector *xsmf, Vector *betaf, Vector *uf)
{

	int i;

	int face, pair;
	int element, neighbor;

	double betaj;

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		if (pair != -1)
		{

			neighbor = faces[pair].element;

			if (V_GetCmp(uf, face + 1) > 0.0)
				betaj = V_GetCmp(betaf, face + 1);
			else
				betaj = 1.0 - V_GetCmp(betaf, face + 1);
			 
			V_SetCmp(xsf, face + 1, MAX(MIN((1.0 - betaj) * V_GetCmp(xs, element + 1) + betaj * V_GetCmp(xs, neighbor + 1), 1.0), 0.0));
			
			V_SetCmp(xsmf, face + 1, MAX(MIN((1.0 - betaj) * V_GetCmp(xsm, element + 1) + betaj * V_GetCmp(xsm, neighbor + 1), 1.0), 0.0));
			
		}
		else
		{
		
			if (faces[face].bc == OUTLET)
			{
			
				// zero gradient
				
				V_SetCmp(xsf, face + 1, V_GetCmp(xs, element + 1));
			
			}
		
		}

	}

}

void CalculateVelocityField(Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xpf,
					        Vector *hu, Vector *hv, Vector *hw, Vector *ap, Vector *xs)
{

	int i, j;

	int element;
	
	msh_vector gradp;
						
	for (i = 0; i < nbelements; i++)
	{

		element = i;

		gradp = Gradient(xp, xpf, LOGICAL_TRUE, element);

		V_SetCmp(xu, element + 1, V_GetCmp(hu, element + 1) / V_GetCmp(ap, element + 1) - gradp.x / V_GetCmp(ap, element + 1));
		V_SetCmp(xv, element + 1, V_GetCmp(hv, element + 1) / V_GetCmp(ap, element + 1) - gradp.y / V_GetCmp(ap, element + 1));
		V_SetCmp(xw, element + 1, V_GetCmp(hw, element + 1) / V_GetCmp(ap, element + 1) - gradp.z / V_GetCmp(ap, element + 1));
		
	}

}

void SetMaterialProperties(Vector *dens, Vector *visc, Vector *spheat, Vector *thcond, Vector *xs)
{

	int i;

	int element;

	double fr[2];

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		fr[0] = (1.0 - V_GetCmp(xs, element + 1));
		fr[1] = V_GetCmp(xs, element + 1);

		V_SetCmp(dens, element + 1, material.dens[0].constant * fr[0] + material.dens[1].constant* fr[1]);
		V_SetCmp(visc, element + 1, material.visc[0].constant * fr[0] + material.visc[1].constant* fr[1]);
		V_SetCmp(spheat, element + 1, material.therm[0].spheat * fr[0] + material.therm[1].spheat* fr[1]);
		V_SetCmp(thcond, element + 1, material.therm[0].thcond * fr[0] + material.therm[1].thcond* fr[1]);

	}

}

void ParseErrorVolume(int phi, int volume)
{

	printf("\nError: Problem evaluating variable %d function f(x,y,z) in volume: %d\n", phi, volume); 
	exit(LOGICAL_ERROR);

}

void SetInitialConditions(Vector *xu0, Vector *xv0, Vector *xw0, Vector *xp0, Vector *xT0, Vector *xs0, 
						  Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs, 
						  Vector *dens, Vector *visc, Vector *shear, Vector *thcond, Vector *spheat)
{

	int i, j;
	
	int element, volume;

	int rv;
	double value;
	
	for (i = 0; i < nbelements; i++)
	{

		element = i;

		elements[element].bc = NONE;

		V_SetCmp(xu0, element + 1, 0.0);
		V_SetCmp(xv0, element + 1, 0.0);
		V_SetCmp(xw0, element + 1, 0.0);
		V_SetCmp(xp0, element + 1, 0.0);
		V_SetCmp(xT0, element + 1, 0.0);
		V_SetCmp(xs0, element + 1, 0.0);

		V_SetCmp(xu, element + 1, 0.0);
		V_SetCmp(xv, element + 1, 0.0);
		V_SetCmp(xw, element + 1, 0.0);
		V_SetCmp(xp, element + 1, 0.0);
		V_SetCmp(xT, element + 1, 0.0);
		V_SetCmp(xs, element + 1, 0.0);

		V_SetCmp(dens, element + 1, 0.0);
		V_SetCmp(visc, element + 1, 0.0);
		V_SetCmp(shear, element + 1, 0.0);
		V_SetCmp(thcond, element + 1, 0.0);
		V_SetCmp(spheat, element + 1, 0.0);

	}

	for (j = 0; j < nbbcvolumes; j++)
	{
		
		volume = j;

		for (i = 0; i < nbelements; i++)
		{

			element = i;

			if (elements[element].physreg == bcvolumes[volume].physreg)
			{

				elements[element].bc = bcvolumes[volume].bc;

				x = elements[element].celement.x;
				y = elements[element].celement.y;
				z = elements[element].celement.z;
				
				strcpy(gs, bcvolumes[volume].fu);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xu, element + 1, value);

				strcpy(gs, bcvolumes[volume].fv);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xv, element + 1, value);

				strcpy(gs, bcvolumes[volume].fw);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xw, element + 1, value);
												
				strcpy(gs, bcvolumes[volume].fp);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xp, element + 1, value);

				strcpy(gs, bcvolumes[volume].fT);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xT, element + 1, value);
												
				strcpy(gs, bcvolumes[volume].fs);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xs, element + 1, value);
							
			}  
			
		}
		
	}
	

}

void SetInitialFlux(Vector *xu, Vector *xv, Vector *xw, Vector *uf)
{

	int i;

	int face, pair;
	int element, neighbor;

	double dNf, dPf;
	double lambda;
	
	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		if (pair != -1)
		{

			neighbor = faces[pair].element;

			dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[pair].cface)); 
			dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 
	
			lambda = dPf / (dPf + dNf);

			// Element face velocity 
			V_SetCmp(uf, face + 1, (V_GetCmp(xu, neighbor + 1) * lambda + V_GetCmp(xu, element + 1) * (1.0 - lambda)) * faces[face].n.x + 
								   (V_GetCmp(xv, neighbor + 1) * lambda + V_GetCmp(xv, element + 1) * (1.0 - lambda)) * faces[face].n.y + 
								   (V_GetCmp(xw, neighbor + 1) * lambda + V_GetCmp(xw, element + 1) * (1.0 - lambda)) * faces[face].n.z);

		}
		else
		{

			// Element face velocity 
			V_SetCmp(uf, face + 1, V_GetCmp(xu, element + 1) * faces[face].n.x + 
								   V_GetCmp(xv, element + 1) * faces[face].n.y + 
								   V_GetCmp(xw, element + 1) * faces[face].n.z);
		
		}
	}

}

void ParseErrorSurface(int phi, int surface)
{

	printf("\nError: Problem evaluating variable %d function f(x,y,z) on surface: %d\n", phi, surface); 
	exit(LOGICAL_ERROR);

}

void SetBoundary(Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf)
{

	int i, j, n;

	int face, pair, surface;

	int rv;	
	double value;
	
	int cyclic[2];
	
	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		pair = faces[face].pair;

		if (pair != -1)
		{
			faces[face].bc = NONE;
		}
			
		V_SetCmp(xuf, face + 1, 0.0);
		V_SetCmp(xvf, face + 1, 0.0);
		V_SetCmp(xwf, face + 1, 0.0);
		V_SetCmp(xpf, face + 1, 0.0);
		V_SetCmp(xTf, face + 1, 0.0);
		V_SetCmp(xsf, face + 1, 0.0);

	}

	for (j = 0; j < nbbcsurfaces; j++)
	{
		
		surface = j;

		for (i = 0; i < nbfaces; i++)
		{

			face = i;
			
			pair = faces[face].pair;
					
			if (pair != -1) continue;  

			if (faces[face].physreg == bcsurfaces[surface].physreg)
			{
			
				faces[face].bc = bcsurfaces[surface].bc;
				
				x = faces[face].cface.x;
				y = faces[face].cface.y;				
				z = faces[face].cface.z;
								
				strcpy(gs, bcsurfaces[surface].fu);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xuf, face + 1, value);
								
				strcpy(gs, bcsurfaces[surface].fv);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xvf, face + 1, value);

				strcpy(gs, bcsurfaces[surface].fw);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xwf, face + 1, value);
				
				strcpy(gs, bcsurfaces[surface].fp);
				rv = evaluate(gs, &value);
				strcat(gs, "\n");
				V_SetCmp(xpf, face + 1, value);
									
				strcpy(gs, bcsurfaces[surface].fT);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xTf, face + 1, value);

				strcpy(gs, bcsurfaces[surface].fs);
				strcat(gs, "\n");
				rv = evaluate(gs, &value);
				V_SetCmp(xsf, face + 1, value);
				
			}  
			
		}
		
	}

	n = 0;
	
	for (i = 0; i < nbfaces; i++)
	{
		face = i;
		
		if (faces[face].bc == CYCLIC)
		{
			if (n < 2)
			{
				cyclic[n] = face;
				n++;
			}
		}
	}
		
	if (n == 2)
	{
		faces[cyclic[0]].pair = cyclic[1];
		faces[cyclic[1]].pair = cyclic[0];
				
		//printf("element: %d, neighbor: %d\n", faces[cyclic[0]].element, faces[cyclic[1]].element);
	}

}

void CalculateCorrectionFactors(QMatrix *Am,
				Vector *xu, Vector *xv, Vector *xw, 
				Vector *hu, Vector *hv, Vector *hw)
{

	AddAsgn_VV(hu, Mul_QV(Sub_QQ(Diag_Q(Am), Am), xu));
	AddAsgn_VV(hv, Mul_QV(Sub_QQ(Diag_Q(Am), Am), xv));
	AddAsgn_VV(hw, Mul_QV(Sub_QQ(Diag_Q(Am), Am), xw));
	
}

void CheckMassConservation(Vector *mc, Vector *uf, double dt)
{

	int i, j;

	int face, pair;
	
	int element, neighbor;

	double mcp;

	double sum;

	sum = 0.0;

	for (i = 0; i < nbelements; i++)
	{

		element = i;
	
		mcp = 0.0;

		for (j = 0; j < elements[element].nbfaces; j++)
		{
			
			face = elements[element].face[j];

			pair = faces[face].pair;

			mcp += V_GetCmp(uf, face + 1) * faces[face].Aj;		
	
		}

		sum += ABS(mcp);

		V_SetCmp(mc, element + 1, ABS(mcp));

	}

	if (verbose == LOGICAL_TRUE)	
		printf("\nMass conservation error: %+E kg\n", sum);

}

void BoundScalar(Vector *xs, double min, double max)
{

	int i;

	int element;

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		V_SetCmp(xs, element + 1, MAX(MIN(V_GetCmp(xs, element + 1), max), min));

	}

}

void SmoothScalar(Vector *xsm, Vector *xs, int n)
{

	int i, j, k;

	int face, pair;
	
	int element, neighbor;

	double sj;

	double dNf, dPf;
	double lambda;

	double sum1, sum2;

	double *sa;

	sa = calloc(nbelements, sizeof(double));
	
	for (i = 0; i < nbelements; i++)
	{

		element = i;

		sa[i] = V_GetCmp(xs, element + 1);

	}

	for (k = 0; k < n; k++)
	{

		for (i = 0; i < nbelements; i++)
		{

			element = i;

			sum1 = 0.0;

			sum2 = 0.0;

			for (j = 0; j < elements[element].nbfaces; j++)
			{
			
				face = elements[element].face[j];

				pair = faces[face].pair;

				if (pair != -1)
				{
					neighbor = faces[pair].element;

					dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[pair].cface)); 
					dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 
	
					lambda = dPf / (dPf + dNf);
					
					sj = sa[neighbor] * lambda + sa[element] * (1.0 - lambda);

				}
				else
				{

					sj = sa[element];
				
				}

				sum1 += sj * faces[face].Aj;
	
				sum2 += faces[face].Aj;

			}

			V_SetCmp(xsm, element + 1, sum1 / sum2);

		}

		for (i = 0; i < nbelements; i++)
		{

			element = i;

			sa[i] = V_GetCmp(xsm, element + 1);
		}

	}

	free(sa);

}

void CalculateMassFraction(Vector *xs, Vector *dens)
{

	int i;

	int element;

	double f[2];
	double vol[2];

	vol[0] = 0.0;
	vol[1] = 0.0;

	for (i = 0; i < nbelements; i++)
	{
		element = i;

		f[0] = (1.0 - V_GetCmp(xs, element + 1));
		f[1] = V_GetCmp(xs, element + 1);

		vol[0] += f[0] * V_GetCmp(dens, element + 1) * elements[element].Vp; 
		vol[1] += f[1] * V_GetCmp(dens, element + 1) * elements[element].Vp; 
	}

	printf ("\nMass of fluid %d: %+E kg\n", 0, vol[0]);
	printf ("\nMass of fluid %d: %+E kg\n", 1, vol[1]);

}

double CalculateMaxCourantNumber(Vector *Co, Vector *uf, double dt)
{

	int i, j;

	int element, face;
	
	double Cpp;
	double Cj;

	double maxCp;

	for (i = 0; i < nbelements; i++)
	{

		element = i;
		
		Cpp = 0.0;
	
		for (j = 0; j < elements[element].nbfaces; j++)
		{
			
			face = elements[element].face[j];
	
			Cj = MAX(-V_GetCmp(uf, face + 1) * faces[face].Aj * dt / elements[element].Vp, 0.0);
	
			Cpp += Cj;

		}

		maxCp = MAX(maxCp, Cpp);
			
		V_SetCmp(Co, element + 1, Cpp);

	}

	if (verbose == LOGICAL_TRUE)	
		printf("\nMaximum Courant number: %.3f\n", maxCp);

	return maxCp;

}

void PredictBeta(Vector *betaf, Vector *xsf, Vector *xs, Vector *Co, Vector *uf)
{
	
	int i;

	int face, pair;
	int element, neighbor;
	int donor, acceptor;

	double dot, l1, l2;

	double su, sdn, sjnCBC, sjnUQ, sjn;

	double qj, tetaj;

	double betaj;

	msh_vector grads;

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		betaj = 0.5;

		if (pair != -1)
		{

			neighbor = faces[pair].element;

			if (V_GetCmp(uf, face + 1) != 0.0)
			{

				if (V_GetCmp(uf, face + 1) > 0.0)
				{

					acceptor = neighbor;
					donor = element;

					grads = Gradient(xs, xsf, LOGICAL_FALSE, donor);

					dot = grads.x * faces[face].d.x +
						  grads.y * faces[face].d.y +
						  grads.z * faces[face].d.z;

					l1 = GeoMagVector(grads);

					l2 = GeoMagVector(faces[face].d); 

				}
				else
				{

					acceptor = element;
					donor = neighbor;

					grads = Gradient(xs, xsf, LOGICAL_FALSE, donor);

					dot = grads.x * faces[pair].d.x +
						  grads.y * faces[pair].d.y +
		  				  grads.z * faces[pair].d.z;
	
					l1 = GeoMagVector(grads);

					l2 = GeoMagVector(faces[pair].d); 

				}	

				su = MIN(MAX(V_GetCmp(xs, acceptor + 1) - 2 * dot, 0.0), 1.0);  

				if ((V_GetCmp(xs, acceptor + 1) - su) != 0.0)
				{

					sdn = (V_GetCmp(xs, donor + 1) - su) / (V_GetCmp(xs, acceptor + 1) - su);

					if (sdn >= 0.0 && sdn <= 1.0)
						sjnCBC = MIN(1, sdn / V_GetCmp(Co, donor + 1));		
					else
						sjnCBC = sdn;		

					if (sdn >= 0.0 && sdn <= 1.0)
						sjnUQ = MIN((8.0 * V_GetCmp(Co, donor + 1) * sdn + (1.0 - V_GetCmp(Co, donor + 1)) * (6.0 * sdn + 3.0)) / 8.0, sjnCBC);		
					else
						sjnUQ = sdn;		
			
					tetaj = acos(ABS(dot / (l1 * l2)));

					qj = MIN(parameter.kq * 0.5 * (cos(2 * tetaj) + 1.0), 1.0);
	
					sjn = qj * sjnCBC + (1 - qj) * sjnUQ;

					if ((1.0 - sdn) != 0.0)
					{
						betaj = MIN(MAX((sjn - sdn) / (1.0 - sdn), 0.0), 1.0);
					}
					else
					{
						if ((sjn - sdn) == 0.0)
							betaj = 0.0;
						else
							betaj = 1.0;
					}

				}
			
			}

		}
		
		V_SetCmp(betaf, face + 1, betaj);
		
	}

}

void CorrectBeta(Vector *betaf, Vector *xsf, Vector *xs, Vector *xs0, Vector *Co, Vector *uf, double dt)
{

	int i;

	int face, pair;	
	int element, neighbor;
	int donor, acceptor;

	double Cj;

	double cbetaj, betaj;

	double ds, Ep, Em, Ec;

	for (i = 0; i < nbfaces; i++)
	{

		face = i;

		element = faces[face].element;

		pair = faces[face].pair;

		cbetaj = 0.0;

		betaj = V_GetCmp(betaf, face + 1);
	
		if (pair != -1)
		{

			neighbor = faces[pair].element;

			if (V_GetCmp(uf, face + 1) != 0.0)
			{
		
				if (V_GetCmp(uf, face + 1) > 0.0)
				{

					acceptor = neighbor;
					donor = element;

				}
				else
				{
					acceptor = element;
					donor = neighbor;

				}

				Cj = MAX(-V_GetCmp(uf, face + 1) * faces[face].Aj * dt / elements[element].Vp, 0.0);

				ds = 0.5 * (V_GetCmp(xs0, acceptor + 1) + V_GetCmp(xs, acceptor + 1)) - 
			 	     0.5 * (V_GetCmp(xs0, donor + 1) + V_GetCmp(xs, donor + 1));

				if (V_GetCmp(xs, donor + 1) < 0.0)
				{ 
					Em = MAX(-V_GetCmp(xs, donor + 1), 0.0);

					// Donor value < 0.0 Ex: sd = -0.1 -> Em = +0.1 
					if (Em > 0.0)
					{
						if (ds > Em)
						{
							cbetaj = Em * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (ds - Em));
							
							cbetaj = MIN(cbetaj, betaj);
						}
					}

				}
				
				
				if (V_GetCmp(xs, donor + 1) > 1.0)
				{ 
				
					Ep = MAX(V_GetCmp(xs, donor + 1) - 1.0, 0.0);

					// Donor value > 1.0 Ex: sd = 1.1 -> Ep = +0.1 
					if (Ep > 0.0)
					{
						if (ds < -Ep)
						{
							cbetaj = Ep * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (-ds - Ep));
							
							cbetaj = MIN(cbetaj, betaj);
						}
					}
				}
				
			}
						
		}

		betaj -= cbetaj;
		
		betaj = MAX(betaj, 0.0);
		
		V_SetCmp(betaf, face + 1, betaj);

	}

}

/**************************************** START: SOLVE MATRIX ****************************************/

void SolveMatrix(QMatrix *A, Vector *x, Vector *b, int *iter, double *res, int msolver, int mprecond)
{

	
	PrecondProcType PrecondProc;
	
	// (0-Null, 1-Jacobi, 2-SOR, 3-ILU)
	
	switch (mprecond)
	{
	case 0:
	
		PrecondProc = NULL;
		break;
	
	case 1:

		PrecondProc = JacobiPrecond;
		break;

	case 2:
	
		PrecondProc = SSORPrecond;
		break;

	case 3:
	
		PrecondProc = ILUPrecond;
		break;
					
	}
	
	// (0-Jacobi, 1-SOR, 2-CGN, 3-GMRES, 4-BiCG, 5-QMR, 6-CGS, 7-BiCGStab, 8-BiCGStabM) 
	
	switch (msolver)
	{
	case 0:
	
		// Jacboi
		
		JacobiIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;
	
	case 1:
	
		// SOR
		
		SSORIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;

	case 2:
	
		// CGN
		
		CGNIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;

	case 3:
	
		// GMRES
		
		GMRESIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;

	case 4:
	
		// BiCG
		
		BiCGIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;

	case 5:
	
		// QMR
		
		QMRIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;

	case 6:
	
		// CGS
		
		CGSIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;

	case 7:
	
		// BiCGStab
		
		BiCGSTABIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
		break;
													
	case 8:

		// BiCGStabM
		
		BiCGSTABIter(A, x, b, parameter.miter, PrecondProc, 1.0);

		*iter = GetLastNoIter();
		*res  = GetLastAccuracy();
	
		if (*iter == parameter.miter)
		{

			V_SetAllCmp(x, 0.0);

			BiCGSTABIter(A, x, b, parameter.miter, NULL, 1.0);

			*iter = GetLastNoIter();
			*res  = GetLastAccuracy();
		
		}
		break;
	
	}

}

/**************************************** END: SOLVE MATRIX ****************************************/

/**************************************** START: BUILD MATRIX ****************************************/

void BuildMomentumMatrix(QMatrix Am, Vector bu, Vector bv, Vector bw, Vector hu, Vector hv, Vector hw, Vector ap, 
						 Vector *xu0, Vector *xv0, Vector *xw0, Vector *xp0, Vector *xT0, Vector *xs0,
						 Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs,
						 Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf,
						 Vector *xsm, Vector *xsmf,
						 Vector *dens, Vector *visc,
						 Vector *uf, double dt)
{
	
	int i, j, k, n;

	int face, pair, element, neighbor;

	double app;

	double apn[MAXFACES];
	int    ani[MAXFACES];

	double bpu, bpv, bpw;

	double densj;
	double viscj;

	msh_vector g;

	msh_vector gradup, gradvp, gradwp;
	msh_vector gradun, gradvn, gradwn;

	msh_vector gradvisc;
	msh_vector grads;
	msh_vector gradp;
	msh_vector graddens;

	msh_vector gradsmp;
	msh_vector gradsmn;
	msh_vector gradsj;
			
	double na;

	double Kp;
 
	double dNf, dPf;
	double lambda;
	double xsi;
	
	g.x = parameter.g[0];
	g.y = parameter.g[1];
	g.z = parameter.g[2];

	// Equation: dU/dt + div(rho*U*U) - div(mi*grad(U)) = qU

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		bpu = 0.0;
		bpv = 0.0;
		bpw = 0.0;
		
		Kp = 0.0;

		app = 0.0;

		n = 0;

		gradup = Gradient(xu, xuf, LOGICAL_TRUE, element);
		gradvp = Gradient(xv, xvf, LOGICAL_TRUE, element);
		gradwp = Gradient(xw, xwf, LOGICAL_TRUE, element);

		for (j = 0; j < elements[element].nbfaces; j++)
		{
			
			face = elements[element].face[j];
			
			pair = faces[face].pair;

			if (pair != -1)
			{

				neighbor = faces[pair].element;

				dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[pair].cface)); 
				dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 
	
				lambda = dPf / (dPf + dNf);
				
				densj = material.dens[0].constant * (1.0 - V_GetCmp(xsf, face + 1)) + material.dens[1].constant * V_GetCmp(xsf, face + 1);
				viscj = material.visc[0].constant * (1.0 - V_GetCmp(xsf, face + 1)) + material.visc[1].constant * V_GetCmp(xsf, face + 1);

				if (parameter.scheme == UPWIND)
				{
					// Upwind scheme
					if (V_GetCmp(uf, face + 1) > 0.0)   
						xsi = 0.0;
					else
						xsi = 1.0;
				
				}

				if (parameter.scheme == CDS)
				{
					// Central-differencing scheme
					xsi = lambda;
				}

				app += (1.0 - xsi) * densj * V_GetCmp(uf, face + 1) * faces[face].Aj / elements[element].Vp + viscj * faces[face].Aj / faces[face].dj / elements[element].Vp;

				apn[n] = xsi * densj * V_GetCmp(uf, face + 1) * faces[face].Aj / elements[element].Vp - viscj * faces[face].Aj / faces[face].dj / elements[element].Vp;
		
				ani[n] = neighbor;
				n++;

				gradun = Gradient(xu, xuf, LOGICAL_TRUE, neighbor);
				gradvn = Gradient(xv, xvf, LOGICAL_TRUE, neighbor);
				gradwn = Gradient(xw, xwf, LOGICAL_TRUE, neighbor);
				
				// Non-orthogonal correction terms
				bpu += -viscj  * faces[face].Aj / faces[face].dj / elements[element].Vp * 
					   (GeoDotVectorVector(gradun, GeoSubVectorVector(faces[face].rnl, elements[neighbor].celement)) -
				        GeoDotVectorVector(gradup, GeoSubVectorVector(faces[face].rpl, elements[element].celement)));

				bpv += -viscj  * faces[face].Aj / faces[face].dj / elements[element].Vp * 
					   (GeoDotVectorVector(gradvn, GeoSubVectorVector(faces[face].rnl, elements[neighbor].celement)) -
				        GeoDotVectorVector(gradvp, GeoSubVectorVector(faces[face].rpl, elements[element].celement)));

				bpw += -viscj  * faces[face].Aj / faces[face].dj / elements[element].Vp * 
					   (GeoDotVectorVector(gradwn, GeoSubVectorVector(faces[face].rnl, elements[neighbor].celement)) -
				        GeoDotVectorVector(gradwp, GeoSubVectorVector(faces[face].rpl, elements[element].celement)));
			
				// Calculate only if surface tension is non-zero
				if (material.tens != 0.0)
				{

					// Curvature divergence (smoothed indicator function)
					
					gradsmp = Gradient(xs, xsmf, LOGICAL_TRUE, element);
					gradsmn = Gradient(xs, xsmf, LOGICAL_TRUE, neighbor);

					gradsj.x = gradsmn.x * lambda + gradsmp.x * (1.0 - lambda);
					gradsj.y = gradsmn.y * lambda + gradsmp.y * (1.0 - lambda);
					gradsj.z = gradsmn.z * lambda + gradsmp.z * (1.0 - lambda);
					
					na = (gradsj.x * faces[face].A.x + gradsj.y * faces[face].A.y + gradsj.z * faces[face].A.z) / (GeoMagVector(gradsj) + 1E-8);
									
					Kp += -1.0 / elements[element].Vp * na;
				}

			
			}
			else
			{

				densj = V_GetCmp(dens, element + 1); 
				viscj = V_GetCmp(visc, element + 1);

				app += viscj * faces[face].Aj / faces[face].dj / elements[element].Vp;

				// Source - face velocity
				bpu += -densj * V_GetCmp(uf, face + 1) * faces[face].Aj * V_GetCmp(xuf, face + 1) / elements[element].Vp + viscj * faces[face].Aj / faces[face].dj * V_GetCmp(xuf, face + 1) / elements[element].Vp;
				bpv += -densj * V_GetCmp(uf, face + 1) * faces[face].Aj * V_GetCmp(xvf, face + 1) / elements[element].Vp + viscj * faces[face].Aj / faces[face].dj * V_GetCmp(xvf, face + 1) / elements[element].Vp;
				bpw += -densj * V_GetCmp(uf, face + 1) * faces[face].Aj * V_GetCmp(xwf, face + 1) / elements[element].Vp + viscj * faces[face].Aj / faces[face].dj * V_GetCmp(xwf, face + 1) / elements[element].Vp;

				// Non-orthogonal correction terms
				bpu += viscj  * faces[face].Aj / faces[face].dj / elements[element].Vp * 
					    GeoDotVectorVector(gradup, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

				bpv += viscj  * faces[face].Aj / faces[face].dj / elements[element].Vp * 
					    GeoDotVectorVector(gradvp, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

				bpw += viscj  * faces[face].Aj / faces[face].dj / elements[element].Vp * 
					    GeoDotVectorVector(gradwp, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

			}

		}

		// Unsteady term - Euler		
		if (dt > 0) 
		{
			app += V_GetCmp(dens, element + 1) / dt;

			bpu += V_GetCmp(dens, element + 1) / dt * V_GetCmp(xu0, element + 1);
			bpv += V_GetCmp(dens, element + 1) / dt * V_GetCmp(xv0, element + 1);
			bpw += V_GetCmp(dens, element + 1) / dt * V_GetCmp(xw0, element + 1);
		}

		// Source - viscous term
		gradvisc = Gradient(visc, NULL, LOGICAL_FALSE, element);
	
		bpu += GeoDotVectorVector(gradup, gradvisc);
		bpv += GeoDotVectorVector(gradvp, gradvisc);
		bpw += GeoDotVectorVector(gradwp, gradvisc);
	
		// Source - gravity	
						
		graddens = Gradient(dens, NULL, LOGICAL_FALSE, element);
				
		bpu += -(g.x * elements[element].celement.x + g.y * elements[element].celement.y + g.z * elements[element].celement.z) * graddens.x;
		bpv += -(g.x * elements[element].celement.x + g.y * elements[element].celement.y + g.z * elements[element].celement.z) * graddens.y;
		bpw += -(g.x * elements[element].celement.x + g.y * elements[element].celement.y + g.z * elements[element].celement.z) * graddens.z;
			
		/*				
		bpu += g.x * V_GetCmp(dens, element + 1);
		bpv += g.y * V_GetCmp(dens, element + 1);
		bpw += g.z * V_GetCmp(dens, element + 1);
		*/
					
		// Source - surface tension
		grads = Gradient(xs, xsf, LOGICAL_TRUE, element);

		bpu += material.tens * Kp * grads.x;
		bpv += material.tens * Kp * grads.y;
		bpw += material.tens * Kp * grads.z;
		
		// Initialize H with source contribution without pressure 

		V_SetCmp(&hu, element + 1, bpu);
		V_SetCmp(&hv, element + 1, bpv);
		V_SetCmp(&hw, element + 1, bpw);

		// Source - pressure 
		gradp = Gradient(xp, xpf, LOGICAL_TRUE, element);

		bpu += -gradp.x;
		bpv += -gradp.y;
		bpw += -gradp.z;

		if (app == 0.0)
		{
			printf("\nError: Problem setting up momentum matrix\n");
			exit(LOGICAL_ERROR);
		}

		V_SetCmp(&ap, element + 1, app);

		Q_SetLen(&Am, element + 1, n + 1);

		Q_SetEntry(&Am, element + 1, 0, element + 1, app);	

		for (j = 0; j < n; j++)
		{
			Q_SetEntry(&Am, element + 1, j + 1, ani[j] + 1, apn[j]);	
		}

		V_SetCmp(&bu, element + 1, bpu);
		V_SetCmp(&bv, element + 1, bpv);
		V_SetCmp(&bw, element + 1, bpw);

	}

}

void BuildContinuityMatrix(QMatrix Ac, Vector bp, Vector *hu, Vector *hv, Vector *hw, Vector *ap, 
				Vector *xu0, Vector *xv0, Vector *xw0, Vector *xp0, Vector *xT0, Vector *xs0,  
				Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs,  
				Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf,  
				Vector *uf,
				double dt)
{

	int i, j, n;

	int element, neighbor, face, pair;

	double acp;
	double acn[MAXFACES];
	int    ani[MAXFACES];
	double bcp;

	double apj;
	
	double Huj, Hvj, Hwj;
	double Hf;
				
	double dNf, dPf;
	double lambda;
	
	msh_vector gradpp;
	msh_vector gradpn;
				
	// Equation: div(U) = 0

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		acp = 0.0;

		bcp = 0.0;

		n = 0;

		gradpp = Gradient(xp, xpf, LOGICAL_TRUE, element);

		for (j = 0; j < elements[element].nbfaces; j++)
		{

			face = elements[element].face[j];

			pair = faces[face].pair;

			if (pair != -1)
			{

				neighbor = faces[pair].element;

				dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[face].cface)); 
				dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 

				lambda = dPf / (dPf + dNf);
					
				apj = V_GetCmp(ap, neighbor + 1) * lambda + V_GetCmp(ap, element + 1) * (1.0 - lambda);
												
				Huj = V_GetCmp(hu, neighbor + 1) * lambda + V_GetCmp(hu, element + 1) * (1.0 - lambda);
				Hvj = V_GetCmp(hv, neighbor + 1) * lambda + V_GetCmp(hv, element + 1) * (1.0 - lambda);
				Hwj = V_GetCmp(hw, neighbor + 1) * lambda + V_GetCmp(hw, element + 1) * (1.0 - lambda);
				
				Hf = 1 / apj * (Huj * faces[face].n.x + Hvj * faces[face].n.y + Hwj * faces[face].n.z); 

				gradpn = Gradient(xp, xpf, LOGICAL_TRUE, neighbor);

				acp += -1 / (apj * faces[face].dj) * faces[face].Aj;

				acn[n] = 1 / (apj * faces[face].dj) * faces[face].Aj;

				ani[n] = neighbor;
				n++;

				V_SetCmp(uf, face + 1, Hf);

				bcp += V_GetCmp(uf, face + 1) * faces[face].Aj;
				
				// Non-orthogonal correction term
				bcp += -1 / (apj * faces[face].dj) * faces[face].Aj * 
					  (GeoDotVectorVector(gradpn, GeoSubVectorVector(faces[face].rnl, elements[neighbor].celement)) -
				       GeoDotVectorVector(gradpp, GeoSubVectorVector(faces[face].rpl, elements[element].celement))); 
								
			}
			else
			{

				apj = V_GetCmp(ap, element + 1);

				Huj = V_GetCmp(hu, element + 1);
				Hvj = V_GetCmp(hv, element + 1);
				Hwj = V_GetCmp(hw, element + 1);
	
				Hf = Huj / apj * faces[face].n.x + Hvj / apj * faces[face].n.y + Hwj / apj * faces[face].n.z; 
				
				if (faces[face].bc == OUTLET)
				{

					// velocity gradient = 0
					// specified pressure

					acp += -1 / (apj * faces[face].dj) * faces[face].Aj;
					bcp += -1 / (apj * faces[face].dj) * faces[face].Aj * V_GetCmp(xpf, face + 1);

					V_SetCmp(uf, face + 1, Hf);
	
					bcp += V_GetCmp(uf, face + 1) * faces[face].Aj;

					// Non-orthogonal correction term
					bcp += 1 / (apj * faces[face].dj) * (GeoDotVectorVector(gradpp, GeoSubVectorVector(faces[face].rpl, elements[element].celement))) * faces[face].Aj;
					
				}

				if (faces[face].bc == PRESSUREINLET)
				{

					// specified pressure 
					// velocity gradient = 0

					acp += -1 / (apj * faces[face].dj) * faces[face].Aj;
					bcp += -1 / (apj * faces[face].dj) * faces[face].Aj * V_GetCmp(xpf, face + 1);

					V_SetCmp(uf, face + 1, Hf);

					bcp += V_GetCmp(uf, face + 1) * faces[face].Aj;

					// Non-orthogonal correction term
					bcp += 1 / (apj * faces[face].dj) * (GeoDotVectorVector(gradpp, GeoSubVectorVector(faces[face].rpl, elements[element].celement))) * faces[face].Aj;

				}

				if (faces[face].bc == INLET || 
					faces[face].bc == MOVINGWALL || 
					faces[face].bc == WALL ||
					faces[face].bc == ADIABATICWALL ||
					faces[face].bc == SURFACE)
				{

					// pressure gradient = 0
					// specified velocity

					V_SetCmp(uf, face + 1, V_GetCmp(xuf, face + 1) * faces[face].n.x + 
										   V_GetCmp(xvf, face + 1) * faces[face].n.y + 
										   V_GetCmp(xwf, face + 1) * faces[face].n.z);

					bcp += V_GetCmp(uf, face + 1) * faces[face].Aj;

				}
				
			}
			
		}

		if (acp == 0.0)
		{
			printf("\nError: Problem setting up continuity matrix\n");
			exit(LOGICAL_ERROR);
		}

		Q_SetLen(&Ac, element + 1, n + 1);

		Q_SetEntry(&Ac, element + 1, 0, element + 1, acp);

		for (j = 0; j < n; j++)
		{
			Q_SetEntry(&Ac, element + 1, j + 1, ani[j] + 1, acn[j]);	
		}

		V_SetCmp(&bp, element + 1, bcp);

	}

}

void BuildEnergyMatrix(QMatrix Ae, Vector bT, 
					   Vector *xu0, Vector *xv0, Vector *xw0, Vector *xp0, Vector *xT0, Vector *xs0,  
					   Vector *xu, Vector *xv, Vector *xw, Vector *xp, Vector *xT, Vector *xs,  
					   Vector *xuf, Vector *xvf, Vector *xwf, Vector *xpf, Vector *xTf, Vector *xsf,  
					   Vector *dens, Vector *visc, Vector *shear, Vector *spheat, Vector *thcond,
					   Vector *uf,
					   double dt)
{
	
	int i, j, k, n;

	int face, pair, element, neighbor;

	double aep;
	double aen[MAXFACES];
	int    ani[MAXFACES];
 
	double bep;

	double dNf, dPf;
	double lambda;

	msh_vector gradTp;
	msh_vector gradTn;

	msh_vector gradup, gradvp, gradwp;

	double thcondj;
	double spheatj;
	double densj;

	double vterm[3];

	for (i = 0; i < nbelements; i++)
	{

		element = i;

		aep =0.0;
		bep = 0.0;

		n = 0;

		gradTp = Gradient(xT, xTf, LOGICAL_TRUE, element);

		for (j = 0; j < elements[element].nbfaces; j++)
		{
			
			face = elements[element].face[j];
			
			pair = faces[face].pair;

			if (pair != -1)
			{

				neighbor = faces[pair].element;

				dNf = GeoMagVector(GeoSubVectorVector(elements[neighbor].celement, faces[pair].cface)); 
				dPf = GeoMagVector(GeoSubVectorVector(elements[element].celement, faces[face].cface)); 
	
				lambda = dPf / (dPf + dNf);

				densj = material.dens[0].constant * (1.0 - V_GetCmp(xsf, face + 1)) + material.dens[1].constant * V_GetCmp(xsf, face + 1);
				thcondj = material.therm[0].thcond * (1.0 - V_GetCmp(xsf, face + 1)) + material.therm[1].thcond * V_GetCmp(xsf, face + 1);
				spheatj = material.therm[0].spheat * (1.0 - V_GetCmp(xsf, face + 1)) + material.therm[1].spheat * V_GetCmp(xsf, face + 1);

				// Convection 
				aep += (1.0 - lambda) * densj * spheatj * V_GetCmp(xu, element + 1) * faces[face].A.x / elements[element].Vp;
				aep += (1.0 - lambda) * densj * spheatj * V_GetCmp(xv, element + 1) * faces[face].A.y / elements[element].Vp;
				aep += (1.0 - lambda) * densj * spheatj * V_GetCmp(xw, element + 1) * faces[face].A.z / elements[element].Vp;

				// Conduction 
				aep += thcondj * faces[face].Aj / faces[face].dj / elements[element].Vp;

				aen[n] = 0.0;

				// Convection 
				aen[n] += lambda * densj * spheatj * V_GetCmp(xu, element + 1) * faces[face].A.x / elements[element].Vp;
				aen[n] += lambda * densj * spheatj * V_GetCmp(xv, element + 1) * faces[face].A.y / elements[element].Vp;
				aen[n] += lambda * densj * spheatj * V_GetCmp(xw, element + 1) * faces[face].A.z / elements[element].Vp;

				// Conduction
				aen[n] += -thcondj * faces[face].Aj / faces[face].dj / elements[element].Vp;
				
				ani[n] = neighbor;

				n++;

				gradTn = Gradient(xT, xTf, LOGICAL_TRUE, neighbor);

				// Non-orthogonal correction term
				bep += thcondj * faces[face].Aj / faces[face].dj / elements[element].Vp * 
					   (GeoDotVectorVector(gradTn, GeoSubVectorVector(faces[face].rnl, elements[neighbor].celement)) -
				        GeoDotVectorVector(gradTp, GeoSubVectorVector(faces[face].rpl, elements[element].celement)));
			}
			else
			{
	
				// Boundary conductivity				
				thcondj = material.bthcond;

				densj = V_GetCmp(dens, element + 1);
				spheatj = V_GetCmp(spheat, element + 1);

				if (faces[face].bc != EMPTY && faces[face].bc != ADIABATICWALL)
				{
					
					// Conduction
					aep += thcondj * faces[face].Aj / faces[face].dj / elements[element].Vp;

					// Convection
					bep += -densj * spheatj * V_GetCmp(xu, element + 1) * faces[face].A.x / elements[element].Vp * V_GetCmp(xTf, face + 1);
					bep += -densj * spheatj * V_GetCmp(xv, element + 1) * faces[face].A.y / elements[element].Vp * V_GetCmp(xTf, face + 1);
					bep += -densj * spheatj * V_GetCmp(xw, element + 1) * faces[face].A.z / elements[element].Vp * V_GetCmp(xTf, face + 1);

					// Conduction
					bep += thcondj * faces[face].Aj / faces[face].dj / elements[element].Vp * V_GetCmp(xTf, face + 1);
				
					// Non-orthogonal correction term
					bep += thcondj * faces[face].Aj / faces[face].dj / elements[element].Vp * 
						    GeoDotVectorVector(gradTp, GeoSubVectorVector(faces[face].rpl, elements[element].celement));

				}

			}	

		}


		// Unsteady term
		if (dt > 0) 
		{
			aep += V_GetCmp(spheat, element + 1) * V_GetCmp(dens, element + 1) / dt;

			bep += V_GetCmp(spheat, element + 1) * V_GetCmp(dens, element + 1) / dt * V_GetCmp(xT0, element + 1);
		}

		gradup = Gradient(xu, xuf, LOGICAL_TRUE, element);
		gradvp = Gradient(xv, xvf, LOGICAL_TRUE, element);
		gradwp = Gradient(xw, xwf, LOGICAL_TRUE, element);

		// Source - viscous heat dissipation
		vterm[0] = 2.0 * ((gradup.x * gradup.x) + 
			              (gradvp.y * gradvp.y) + 
				          (gradwp.z * gradwp.z));

		vterm[1] = (gradvp.x + gradup.y) * (gradvp.x + gradup.y) + 
				   (gradwp.y + gradvp.z) * (gradwp.y + gradvp.z) +
			       (gradup.z + gradwp.x) * (gradup.z + gradwp.x); 

		vterm[2] = 2.0 / 3.0 * (gradup.x + gradvp.y + gradwp.z) * (gradup.x + gradvp.y + gradwp.z); 

		bep += V_GetCmp(visc, element + 1) * (vterm[0] + vterm[1] + vterm[2]);

		// Shear rate 
		V_SetCmp(shear, element + 1, sqrt(vterm[0] + vterm[1] + vterm[2]));

		if (aep == 0.0)
		{
			printf("\nError: Problem setting up energy matrix\n");
			exit(LOGICAL_ERROR);
		}

		Q_SetLen(&Ae, element + 1, n + 1);

		Q_SetEntry(&Ae, element + 1, 0, element + 1, aep);	

		for (j = 0; j < n; j++)
		{
			Q_SetEntry(&Ae, element + 1, j + 1, ani[j] + 1, aen[j]);	
		}

		V_SetCmp(&bT, element + 1, bep);

	}
	
}

void BuildVolumeOfFluidMatrix(QMatrix As, Vector bs, 
							  Vector *xs0,
							  Vector *xs,
							  Vector *xsf,
							  Vector *betaf, 
							  Vector *uf,
							  double dt)
{

	int i, j, n;

	int element, neighbor, face, pair;

	double aip;
	double ain[MAXFACES];
	int    ani[MAXFACES];
	double bip;

	double dNf, dPf;
	double lambda;

	double betaj;

	// Equation: ds/dt + div(s*U) = 0
	
	for (i = 0; i < nbelements; i++)
	{

		element = i;

		aip = 0;

		bip = 0;

		n = 0;
				
		for (j = 0; j < elements[element].nbfaces; j++)
		{

			face = elements[element].face[j];

			pair = faces[face].pair;

			if (V_GetCmp(uf, face + 1) > 0.0)
				betaj = V_GetCmp(betaf, face + 1);
			else
				betaj = 1.0 - V_GetCmp(betaf, face + 1);

			if (pair != -1)
			{

				neighbor = faces[pair].element;

				aip += 0.5 * (1.0 - betaj) * V_GetCmp(uf, face + 1) * faces[face].Aj;
			
				ain[n] = 0.5 * betaj * V_GetCmp(uf, face + 1) * faces[face].Aj;

				ani[n] = neighbor;
				n++;

				bip += -0.5 * (1.0 - betaj) * V_GetCmp(uf, face + 1) * faces[face].Aj * V_GetCmp(xs, element + 1);
				bip += -0.5 * betaj * V_GetCmp(uf, face + 1) * faces[face].Aj * V_GetCmp(xs, neighbor + 1); 
				
			}
			else
			{
			
				bip += -1.0 * V_GetCmp(uf, face + 1) * faces[face].Aj * V_GetCmp(xsf, face + 1);
		
			}
						 
		}

		if (dt > 0.0)
		{

			aip += elements[element].Vp / dt;
		
			bip += elements[element].Vp / dt * V_GetCmp(xs0, element + 1);   
	
		}

		if (aip == 0.0)
		{
			printf("\nError: Problem setting up volume-of-fluid matrix\n");
			exit(LOGICAL_ERROR);
		}

		Q_SetLen(&As, element + 1, n + 1);

		Q_SetEntry(&As, element + 1, 0, element + 1, aip);

		for (j = 0; j < n; j++)
		{
			Q_SetEntry(&As, element + 1, j + 1, ani[j] + 1, ain[j]);	
		}

		V_SetCmp(&bs, element + 1, bip);

	}
	
}

/**************************************** END: BUILD MATRIX ****************************************/

/**************************************** START: SIMULATION ****************************************/

int Simulation(char *path)
{

	int i, j, n;
	
	char var[6];
	
	double mres[6];
	int    miter[6];

	double fres[6];
	int    fiter[6];

	int iter;
	double pres, temp;

	double time, dt;
	double wtime, wdt;

	double vmin, vmax;

	double maxCp;

	char *file;

	FILE *fpresults;		// Output to gmsh post-processing file (results)
	FILE *fpprobe;			// Output to gnuplot file (probe)
	FILE *fpresiduals;		// Output to gnuplot file (residuals)

	Vector  uf;
	
	QMatrix Am, Ac, Ae, As;
	Vector  bu, bv, bw, bp, bT, bs;
	
	Vector  xu0, xv0, xw0, xp0, xT0, xs0;
	Vector  xu, xv, xw, xp, xT, xs;
	
	Vector  hu, hv, hw;
	Vector  ap;

	Vector  xuf, xvf, xwf, xpf, xTf, xsf;

	Vector  betaf;

	Vector  dens, visc, shear, thcond, spheat;

	Vector xpp, xTp;
	
	Vector xsm;

	Vector xsmf;

	Vector Co;

	Vector mc;

	clock_t start, finish;

	var[iu] = 'u';
	var[iv] = 'v';
	var[iw] = 'w';
	var[ip] = 'p';
	var[iT] = 'T';
	var[is] = 's';

	// Allocate memory
	printf("\n");
	printf("Allocating memory...\n");

	V_Constr(&dens, "Density", nbelements, Normal, True);
	V_Constr(&visc, "Dynamic viscosity", nbelements, Normal, True);
	V_Constr(&shear, "Shear rate", nbelements, Normal, True);
	V_Constr(&thcond, "Thermal conductivity", nbelements, Normal, True);
	V_Constr(&spheat, "Specific heat", nbelements, Normal, True);

	V_Constr(&mc, "Mass conservation", nbelements, Normal, True);

	V_Constr(&Co, "Courant number", nbelements, Normal, True);
	
	V_Constr(&uf, "Face flux velocity", nbfaces, Normal, True);
	
	V_Constr(&bu, "Momentum source x-component", nbelements, Normal, True);
	V_Constr(&bv, "Momentum source y-component", nbelements, Normal, True);
	V_Constr(&bw, "Momentum source z-component", nbelements, Normal, True);
	V_Constr(&bp, "Continuity source", nbelements, Normal, True);

	V_Constr(&xu0, "Velocity x-component at cell center (previous time step)", nbelements, Normal, True);
	V_Constr(&xu, "Velocity x-component at cell center", nbelements, Normal, True);
	V_Constr(&xuf, "Velocity x-component at face center", nbfaces, Normal, True);
		
	V_Constr(&xv0, "Velocity y-component at cell center (previous time step)", nbelements, Normal, True);
	V_Constr(&xv, "Velocity y-component at cell center", nbelements, Normal, True);
	V_Constr(&xvf, "Velocity y-component at face center", nbfaces, Normal, True);
		
	V_Constr(&xw0, "Velocity z-component at cell center (previous time step)", nbelements, Normal, True);
	V_Constr(&xw, "Velocity z-component at cell center ", nbelements, Normal, True);
	V_Constr(&xwf, "Velocity z-component at face center ", nbfaces, Normal, True);
		
	V_Constr(&xp0, "Pressure at cell center (previous time step)", nbelements, Normal, True);
	V_Constr(&xp, "Pressure at cell center", nbelements, Normal, True);
	V_Constr(&xpf, "Pressure at face center", nbfaces, Normal, True);

	V_Constr(&xpp, "Pressure at cell center - previous iteration", nbelements, Normal, True);

	V_Constr(&hu, "Momentum matrix source x-component without pressure", nbelements, Normal, True);
	V_Constr(&hv, "Momentum matrix source y-component without pressure", nbelements, Normal, True);
	V_Constr(&hw, "Momentum matrix source z-component without pressure", nbelements, Normal, True);

	V_Constr(&ap, "Momentum matrix diagonal", nbelements, Normal, True);
	
	V_Constr(&bT, "Energy source", nbelements, Normal, True);

	V_Constr(&xT0, "Temperature at cell center (previous time step)", nbelements, Normal, True);
	V_Constr(&xT, "Temperature at cell center", nbelements, Normal, True);
	V_Constr(&xTf, "Temperature at face center", nbfaces, Normal, True);

	V_Constr(&xTp, "Temperature at cell center - previous iteration", nbelements, Normal, True);
	
	V_Constr(&bs, "Gamma source", nbelements, Normal, True);

	V_Constr(&xs0, "Gamma at cell center (previous time step)", nbelements, Normal, True);
	V_Constr(&xs, "Gamma at cell center", nbelements, Normal, True);
	V_Constr(&xsf, "Gamma at face center", nbfaces, Normal, True);
		
	V_Constr(&xsm, "Smoothed gamma at cell center", nbelements, Normal, True);
	V_Constr(&xsmf, "Smoothed gamma at face center", nbfaces, Normal, True);

	V_Constr(&betaf, "CICSAM interpolation factor", nbfaces, Normal, True);
	
	printf("Memory allocated.\n");
		
	// Set initial conditions
	SetInitialConditions(&xu0, &xv0, &xw0, &xp0, &xT0, &xs0, &xu, &xv, &xw, &xp, &xT, &xs, &dens, &visc, &shear, &thcond, &spheat);
	
	// Smoothing...
	SmoothScalar(&xsm, &xs, 2);
	
	// Set initial flux
	SetInitialFlux(&xu, &xv, &xw, &uf);
		
	// Set boundary velocity and pressure
	SetBoundary(&xuf, &xvf, &xwf, &xpf, &xTf, &xsf);
	
	// Set matrix solution accuracy
	SetRTCAccuracy(parameter.mtol);

	// Set time intervals
	time = parameter.t0;
	dt = parameter.dt;
	
	// Set write time intervals
	wtime = parameter.t0;

	if (parameter.nsav != 0)
		wdt = (parameter.t1 - parameter.t0) / parameter.nsav;
	else
	{
		wdt = 2 * (parameter.t1 - parameter.t0);
	}
	
	// Allocate memory for file string
	file = calloc(strlen(path) + 9, sizeof(char));
	
	if (parameter.steady == LOGICAL_TRUE)
	{
	
		sprintf(file, "%s.res", path);
	
		// Open output files for residuals
		fpresiduals = fopen(file, "w");
	}

	sprintf(file, "%s.pos", path);
	
	// Open output file for results
	fpresults = fopen(file, "w");
	
	n = 0;
	
	WriteResults(fpresults, &xu, &xv, &xw, &xp, &xT, &xs, &xuf, &xvf, &xwf, &xpf, &xTf, &xsf, &uf, LOGICAL_TRUE, LOGICAL_TRUE, time);
		
	wtime = wdt;
	
	iter = 0;

	fiter[iu] = 0;
	fiter[iv] = 0;
	fiter[iw] = 0;
	fiter[ip] = 0;
	fiter[iT] = 0;
	fiter[is] = 0;

	// Start clock to measure simulation time
	start = clock();
			
	do
	{

		time += dt;

		iter++;
	
		printf("\n**** TIME = %f ****\n", time);

		// Set material properties 
		SetMaterialProperties(&dens, &visc, &spheat, &thcond, &xs);

		// Store previous time step values
		Asgn_VV(&xu0, &xu);
		Asgn_VV(&xv0, &xv);
		Asgn_VV(&xw0, &xw);
		Asgn_VV(&xp0, &xp);
		Asgn_VV(&xT0, &xT);
		Asgn_VV(&xs0, &xs);
		
		Q_Constr(&Am, "Am", nbelements, False, Rowws, Normal, True);

		V_SetAllCmp(&hu, 0.0);
		V_SetAllCmp(&hv, 0.0);
		V_SetAllCmp(&hw, 0.0);
		
		V_SetAllCmp(&ap, 1.0);
				
		if (parameter.calc[iu] == LOGICAL_TRUE || parameter.calc[iv] == LOGICAL_TRUE || parameter.calc[iw] == LOGICAL_TRUE)
		{
			
			// Build three momentum matrices for u, v, w velocity components
			BuildMomentumMatrix(Am, bu, bv, bw, hu, hv, hw, ap, 
								&xu0, &xv0, &xw0, &xp0, &xT0, &xs0, 
								&xu, &xv, &xw, &xp, &xT, &xs, 
								&xuf, &xvf, &xwf, &xpf, &xTf, &xsf, 
								&xsm, &xsmf,
								&dens, &visc,
								&uf,
								dt);
		
			if (!CheckIfDiagonalMatrix(&Am, nbelements))
			{
				printf("\nWarning: Momentum matrix is not diagonal dominant\n");
				WriteMatrix(&Am, &bu, nbelements);
				exit(LOGICAL_ERROR);
			}
			
			if (parameter.calc[iu] == LOGICAL_TRUE)
			{

				fiter[iu]++;

				// Solve matrix for u velocity component
				SolveMatrix(&Am, &xu, &bu, &miter[iu], &mres[iu], parameter.msolver[iu], parameter.mprecond[iu]);

				if (verbose == LOGICAL_TRUE)
					printf("\nMatrix %c Number of iterations: %d Residual: %e\n", var[iu], miter[iu], mres[iu]);

				if ((mres[iu] > parameter.mtol && miter[iu] == parameter.miter) || LASResult() != LASOK)
				{
					printf("\nError: Problem solving matrix %c\n", var[iu]);
					WriteMatrix(&Am, &bu, nbelements);
					exit(LOGICAL_ERROR);
				}
			
			}

			if (parameter.calc[iv] == LOGICAL_TRUE)
			{

				fiter[iv]++;

				// Solve matrix for v velocity component
				SolveMatrix(&Am, &xv, &bv, &miter[iv], &mres[iv], parameter.msolver[iv], parameter.mprecond[iv]);

				if (verbose == LOGICAL_TRUE)
					printf("\nMatrix %c Number of iterations: %d Residual: %e\n", var[iv], miter[iv], mres[iv]);

				if (mres[iv] > parameter.mtol && miter[iv] == parameter.miter)
				{
					printf("\nError: Problem solving matrix %c\n", var[iv]);
					WriteMatrix(&Am, &bv, nbelements);
					exit(LOGICAL_ERROR);
				}
				
			}

			if (parameter.calc[iw] == LOGICAL_TRUE)
			{

				fiter[iw]++;

				// Solve matrix for w velocity component
				SolveMatrix(&Am, &xw, &bw, &miter[iw], &mres[iw], parameter.msolver[iw], parameter.mprecond[iw]);

				if (verbose == LOGICAL_TRUE)
					printf("\nMatrix %c Number of iterations: %d Residual: %e\n", var[iw], miter[iw], mres[iw]);

				if ((mres[iw] > parameter.mtol && miter[iw] == parameter.miter) || LASResult() != LASOK)
				{
					printf("\nProblem solving matrix %c\n", var[iw]);
					WriteMatrix(&Am, &bw, nbelements);
					exit(LOGICAL_ERROR);
				}
						
			}

			// Calculate correction factors
			CalculateCorrectionFactors(&Am, &xu, &xv, &xw, &hu, &hv, &hw);

		}
		
		Q_Constr(&Ac, "Ac", nbelements, False, Rowws, Normal, True);
		
		if (parameter.calc[ip] == LOGICAL_TRUE)
		{
		
			fiter[ip]++;

			for (i = 0; i <= parameter.northocor; i++)
			{

				// Store previous iteration values
				Asgn_VV(&xpp, &xp);

				// Build the continuity matrix (mass conservation)
				BuildContinuityMatrix(Ac, bp, &hu, &hv, &hw, &ap, 
									  &xu0, &xv0, &xw0, &xp0, &xT0, &xs0,  
									  &xu, &xv, &xw, &xp, &xT, &xs,  
									  &xuf, &xvf, &xwf, &xpf, &xTf, &xsf,  
									  &uf, dt);

				if (!CheckIfDiagonalMatrix(&Ac, nbelements))
				{
					printf("\nWarning: Continuity matrix is not diagonal dominant\n");
					WriteMatrix(&Ac, &bp, nbelements);
					exit(LOGICAL_ERROR);
				}

				// Solve matrix to get pressure p
				SolveMatrix(&Ac, &xp, &bp, &miter[ip], &mres[ip], parameter.msolver[ip], parameter.mprecond[ip]);

				if (verbose == LOGICAL_TRUE)
					printf("\nMatrix %c Number of iterations: %d Residual: %+E\n", var[ip], miter[ip], mres[ip]);

				if ((mres[ip] > parameter.mtol && miter[ip] == parameter.miter) || LASResult() != LASOK)
				{
					printf("\nError: Problem solving matrix %c\n", var[ip]);
					WriteMatrix(&Ac, &bp, nbelements);
					exit(LOGICAL_ERROR);
				}	

				// Calculate pressure convergence
				pres = l2Norm_V(Sub_VV(&xpp, &xp));
				
				if (verbose == LOGICAL_TRUE)
					printf("\nNon-orthogonality error (continuity): %+E\n", pres);

				if (pres < parameter.mtol)
					break;

				CorrectFaceP(&xp, &xpf);

			}

			// Correct face values
			CorrectFaceUVW(&ap, &xu, &xv, &xw, &xp, &xuf, &xvf, &xwf, &xpf, &uf);   
			
			// Correct cell center
			CalculateVelocityField(&xu, &xv, &xw, &xp, &xpf, &hu, &hv, &hw, &ap, &xs);
					
		}

		if (verbose == LOGICAL_TRUE)
		{
			// Check mass conservation
			CheckMassConservation(&mc, &uf, dt);
		}
		
		Q_Destr(&Am);
		Q_Destr(&Ac);

		if (parameter.calc[iT] == LOGICAL_TRUE)
		{

			Q_Constr(&Ae, "Ae", nbelements, False, Rowws, Normal, True);

			fiter[iT]++;

			for (i = 0; i <= parameter.northocor; i++)
			{

				// Store previous iteration values
				Asgn_VV(&xTp, &xT);
				
				// Build energy matrix
				BuildEnergyMatrix(Ae, bT, 
								  &xu0, &xv0, &xw0, &xp0, &xT0, &xs0,
								  &xu, &xv, &xw, &xp, &xT, &xs,
								  &xuf, &xvf, &xwf, &xpf, &xTf, &xsf,  
								  &dens, &visc, &shear, &spheat, &thcond,
								  &uf,
								  dt);

				if (!CheckIfDiagonalMatrix(&Ae, nbelements))
				{
					printf("\nWarning: Energy matrix is not diagonal dominant\n");
					WriteMatrix(&Ae, &bT, nbelements);
					exit(LOGICAL_ERROR);
				}

				// Solve matrix to get temperature T
				SolveMatrix(&Ae, &xT, &bT, &miter[iT], &mres[iT], parameter.msolver[iT], parameter.mprecond[iT]);

				if (verbose == LOGICAL_TRUE)			
					printf("\nMatrix %c Number of iterations: %d Residual: %+E\n", var[iT], miter[iT], mres[iT]);

				if ((mres[iT] > parameter.mtol && miter[iT] == parameter.miter) || LASResult() != LASOK)
				{
					printf("\nError: Problem solving matrix %c\n", var[iT]);
					WriteMatrix(&Ae, &bT, nbelements);
					exit(LOGICAL_ERROR);
				}

				// Calculate temperature convergence
				temp = l2Norm_V(Sub_VV(&xTp, &xT));
				
				if (verbose == LOGICAL_TRUE)
					printf("\nNon-orthogonality error (energy): %+E\n", temp);

				if (temp < parameter.mtol)
					break;

				CorrectFaceT(&xT, &xTf);

			}
			
			Q_Destr(&Ae);

		}

		// Calculate maximum Courant number 
		maxCp = CalculateMaxCourantNumber(&Co, &uf, dt);

		if (maxCp >= 1.0)
		{
			
			if (parameter.calc[is] == LOGICAL_TRUE) 
			{
				printf("\nError: Courant number >= 1.0\n");
				exit(1);
			}
			else
			{
				printf("\nWarning: Courant number >= 1.0\n");
			}
		}
		
		if (parameter.calc[is] == LOGICAL_TRUE)
		{

			Q_Constr(&As, "As", nbelements, False, Rowws, Normal, True);
		
			fiter[is]++;

			// Predict beta - CICSAM 
			PredictBeta(&betaf, &xsf, &xs, &Co, &uf);

			for (j = 0; j <= parameter.ncicsamcor; j++)
			{
	
				// Build VOF matrix
				BuildVolumeOfFluidMatrix(As, bs, &xs0, &xs, &xsf, &betaf, &uf, dt);

				if (!CheckIfDiagonalMatrix(&As, nbelements))
				{
					printf("\nWarning: Volume-of-fluid matrix is not diagonal dominant\n");
					WriteMatrix(&As, &bs, nbelements);
					exit(LOGICAL_ERROR);
				}
			
				// Solve matrix to get indicator function s	
				SolveMatrix(&As, &xs, &bs, &miter[is], &mres[is], parameter.msolver[is], parameter.mprecond[is]);

				if (verbose == LOGICAL_TRUE)	
					printf("\nMatrix %c Number of iterations: %d Residual: %+E\n", var[is], miter[is], mres[is]);
		
				if ((mres[is] > parameter.mtol && miter[is] == parameter.miter) || LASResult() != LASOK)
				{
					printf("\nError: Problem solving matrix %c\n", var[is]);
					printf("\nExiting...\n");
					exit(LOGICAL_ERROR);
				}
		
				// Correct beta 
				CorrectBeta(&betaf, &xsf, &xs, &xs0, &Co, &uf, dt);
				
				// Correct face and boundary 
				CorrectFaceS(&xs, &xsf, &xsm, &xsmf, &betaf, &uf);
					
			}
			
			Q_Destr(&As);

			// Bound indicator function
			BoundScalar(&xs, 0.0, 1.0);

			// Calculate mass fractions 
			CalculateMassFraction(&xs, &dens);
						
			// Smooth indicator function (two times)
			if (material.tens != 0.0)
			{
				SmoothScalar(&xsm, &xs, 2);
			}
			
		}

		if (parameter.adjdt == LOGICAL_TRUE)
		{
			if (maxCp > parameter.maxCp) dt *= 0.85 * parameter.maxCp / maxCp;
		
			if (1.25 * maxCp < parameter.maxCp) dt *= 1.05;
		}

		if ((time + dt) >= wtime)
		{

			WriteResults(fpresults, &xu, &xv, &xw, &xp, &xT, &xs, &xuf, &xvf, &xwf, &xpf, &xTf, &xsf, &uf, LOGICAL_TRUE, LOGICAL_TRUE, time);

			fflush(fpresults);

			sprintf(file, "%s.prb", path);
							
			fpprobe = fopen(file, "w");

			WriteProbeViews(fpprobe, &xu, &xv, &xw, &xp, &xT, &xs, &xuf, &xvf, &xwf, &xpf, &xTf, &xsf, var, time);

			fclose(fpprobe);

			wtime += wdt;

		}

		if (parameter.steady == LOGICAL_TRUE)
		{
			
		  // Get residual
			
		  if (parameter.calc[iu] == LOGICAL_TRUE) 
		    fres[iu] = l2Norm_V(Sub_VV(&xu, &xu0));
		  else
		    fres[iu] = 0.0;

		  if (parameter.calc[iv] == LOGICAL_TRUE) 
		    fres[iv] = l2Norm_V(Sub_VV(&xv, &xv0));
		  else
		    fres[iv] = 0.0;
		  
		  if (parameter.calc[iw] == LOGICAL_TRUE) 
		    fres[iw] = l2Norm_V(Sub_VV(&xw, &xw0));
		  else
		    fres[iw] = 0.0;
		  
		  if (parameter.calc[ip] == LOGICAL_TRUE) 
		    fres[ip] = l2Norm_V(Sub_VV(&xp, &xp0));
		  else
		    fres[ip] = 0.0;
		  
		  if (parameter.calc[iT] == LOGICAL_TRUE) 
		    fres[iT] = l2Norm_V(Sub_VV(&xT, &xT0));
		  else
		    fres[iT] = 0.0;
		  
		  if (parameter.calc[is] == LOGICAL_TRUE) 
		    fres[iv] = l2Norm_V(Sub_VV(&xs, &xs0));
		  else	
		    fres[is] = 0.0;
		  
				  
		  n = 0;

		  for (i = 0; i < nphi; i++)
		    {
		      
		      if (verbose == LOGICAL_TRUE)
			printf("\nVariable: %c Iteration: %d Final residual: %+E\n", var[i], fiter[i], fres[i]);
		      
		      if (fres[i] > parameter.ftol)
			{	
			  n++;
			}
		    }
		  
		  if (n == 0)
		    {
		      printf("\nSteady state reached.\n");
		      break;
		    }
		  
		  WriteResidual(fpresiduals, iter, fres);
		  
		  fflush(fpresiduals);
		  
		}
		else
		{

			if (time + 0.5 * dt > parameter.t1)
				break;
		
		}

	} while (dt > 0.0);

	finish = clock();

	printf("\nSimulation time: %f seconds.\n", (float) (finish - start) / CLOCKS_PER_SEC);
	
	WriteResults(fpresults, &xu, &xv, &xw, &xp, &xT, &xs, &xuf, &xvf, &xwf, &xpf, &xTf, &xsf, &uf, LOGICAL_TRUE, LOGICAL_TRUE, time);
	
	// Close output files

	fclose(fpresults);

	if (parameter.steady == LOGICAL_TRUE)
	{
	
		WriteResidual(fpresiduals, iter, fres);

		fflush(fpresiduals);
	
		fclose (fpresiduals);
	}

	// Open output file for probes
	sprintf(file, "%s.prb", path);
	
	fpprobe = fopen(file, "w");
	
	WriteProbeViews(fpprobe, &xu, &xv, &xw, &xp, &xT, &xs, &xuf, &xvf, &xwf, &xpf, &xTf, &xsf, var, time);

	fclose(fpprobe);
			
	// Release memory 

	free(file);
		
	V_Destr(&uf);

	V_Destr(&bu);
	V_Destr(&bv);
	V_Destr(&bw);
	V_Destr(&bp);
	V_Destr(&bT);
	V_Destr(&bs);

	V_Destr(&xu0);
	V_Destr(&xv0);
	V_Destr(&xw0);
	V_Destr(&xp0);
	V_Destr(&xT0);
	V_Destr(&xs0);
	
	V_Destr(&xu);
	V_Destr(&xuf);

	V_Destr(&xv);
	V_Destr(&xvf);

	V_Destr(&xw);
	V_Destr(&xwf);

	V_Destr(&xp);
	V_Destr(&xpf);

	V_Destr(&xT);
	V_Destr(&xTf);
	V_Destr(&xs);
	V_Destr(&xsf);

	V_Destr(&xpp);
	V_Destr(&xTp);

	V_Destr(&xsm);
	V_Destr(&xsmf);

	V_Destr(&hu);
	V_Destr(&hv);
	V_Destr(&hw);

	V_Destr(&ap);

	V_Destr(&betaf);

	V_Destr(&dens);
	V_Destr(&visc);
	V_Destr(&shear);
	V_Destr(&thcond);
	V_Destr(&spheat);

	V_Destr(&mc);

	V_Destr(&Co);
	
	return LOGICAL_TRUE;

}

/**************************************** END: SIMULATION ****************************************/

void usage()
{

	printf("\n");
	printf("** Linux **\n");
	printf("Usage:       ./OpenFVM [options] [n] [file]\n");
	printf("Example 1:   ./OpenFVM -f  1 lid/lid\n");
	printf("Example 3:   ./OpenFVM -fv 1 lid/lid\n");
	
	printf("\n");
	printf("** Windows **\n");
	printf("Usage:       OpenFVM [option] [n] [file]\n");
	printf("Example 1:   OpenFVM -f  1 lid/lid\n");
	printf("Example 3:   OpenFVM -fv 1 lid/lid\n");
	
	printf("\n");
	printf("OPTIONS:\n");
	printf("  f   Start simulation\n");
	printf("  v   Verbose mode\n");
	printf("\n");

}

int main(int argc, char *argv[])
{

	char *path;
	char *file;

	char *ptr;
			
	printf("\n");
	printf("*****************************************\n");
	printf("*                                       *\n");
	printf("*    OpenFVM v0.2                       *\n");
	printf("*                                       *\n");
	printf("*****************************************\n");
	printf("\n");
		
	if (argc != 4)
	{
		usage();
		printf("\nError: Wrong number of arguments, expected %d found %d.\n\n", 4 - 1, argc - 1);
		return LOGICAL_ERROR;
	}

	if (strchr(argv[1], 'v') != NULL)
		verbose = LOGICAL_TRUE;
	else
		verbose = LOGICAL_FALSE;
		
	// Simulate
	if (strchr(argv[1], 'f') != NULL)
	{
						
		path = calloc(strlen(argv[3]), sizeof(char));
		file = calloc(strlen(argv[3]) + 9, sizeof(char));

		strcpy(path, argv[3]);
							
		// Read mesh file
		sprintf(file, "%s.msh", path);
		MshImportMSH(file);
				
		if (nbpatches > 0) 
		{
		  nbpatches = 0;
		  free(patches);
		}
	
		// Read boundary conditions file
		sprintf(file, "%s.bcd", path);
		BcdImportBCD(file);

		// Read material file
		sprintf(file, "%s.mtl", path);		
		MtlImportMTL(file);

		// Read processing conditions file
		sprintf(file, "%s.par", path);
		ParImportPAR(file);
		
		// Start simulation
		Simulation(path);
	
	}
	
	printf("Done.\n\n");
		
	return 0;
}

