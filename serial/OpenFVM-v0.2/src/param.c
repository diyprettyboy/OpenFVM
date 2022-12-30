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

#include "ioutils.h"
#include "globals.h"
#include "param.h"

void SetDefaults()
{

	parameter.scheme = 1;
	
	parameter.steady = 0;
	parameter.ftol = 1E-6;

	parameter.wbinary = 0;
	parameter.nsav = 1;
		
	parameter.calc[0] = 1;
	parameter.calc[1] = 1;
	parameter.calc[2] = 1;
	parameter.calc[3] = 1;
	parameter.calc[4] = 0;
	parameter.calc[5] = 0;
	
	parameter.fsav[0] = 0;
	parameter.fsav[1] = 0;
	parameter.fsav[2] = 0;
	parameter.fsav[3] = 0;
	parameter.fsav[4] = 0;
	parameter.fsav[5] = 0;

	parameter.csav[0] = 0;
	parameter.csav[1] = 0;
	parameter.csav[2] = 0;
	parameter.csav[3] = 0;
	parameter.csav[4] = 0;
	parameter.csav[5] = 0;
	
	parameter.probe[0] = 0;
	parameter.probe[1] = 0;
	parameter.probe[2] = 0;
	parameter.probe[3] = 0;
	parameter.probe[4] = 0;
	parameter.probe[5] = 0;
	
	parameter.vortex[0] = 0;
	parameter.vortex[1] = 0;
	parameter.vortex[2] = 0;

	parameter.streamf = 0;
		
	parameter.fvec = 0;
	parameter.cvec = 0;

	parameter.kq = 2.0;
	parameter.ncicsamcor = 3;

	parameter.g[0] = 0.0;
	parameter.g[1] = 0.0;
	parameter.g[2] = 0.0;

	parameter.msolver[0] = 8;
	parameter.msolver[1] = 8;
	parameter.msolver[2] = 8;
	parameter.msolver[3] = 8;
	parameter.msolver[4] = 8;
	parameter.msolver[5] = 3;

	parameter.mprecond[0] = 3;
	parameter.mprecond[1] = 3;
	parameter.mprecond[2] = 3;
	parameter.mprecond[3] = 3;
	parameter.mprecond[4] = 3;
	parameter.mprecond[5] = 3;

	
	parameter.northocor = 10;
	parameter.mtol = 1E-8;
	parameter.miter = 5000;

	parameter.adjdt = 0;
	parameter.maxCp = 0.25;

	parameter.t0 = 0.0;
	parameter.t1 = 0.001;
	parameter.dt = 0.001;

}

int ParImportPAR(char *file)
{

	int i, j, n;

	int inull;

	int tcode;
	int nbpar;

	FILE *fp;
	char descr[512];

	//Set default values
	SetDefaults();
	
	fp = fopen(file, "r");

	if (fp == NULL)
	{
		printf("\nError: Parameter file not found!\n");
		printf("%s\n\n", file);
		exit(LOGICAL_ERROR);
	}

	printf("\nReading parameters...\n");

	printf("\n");

	do
	{

		do
		{
			fscanf(fp, "%s", descr);

			if (strcmp(descr, "$PRMT") == 0)
				break;

			if (strcmp(descr, "$ENDF") == 0)
				break;

		} while (!feof(fp));

		if (strcmp(descr, "$ENDF") == 0)
			break;

		if (strcmp(descr, "$PRMT") == 0)
		{

			fscanf(fp, "%d %d ", &inull, &nbpar);

			GetLine(fp); 
	
			for (i = 0; i < nbpar; i++)
			{

				fscanf(fp, "%d %d", &tcode, &n);

				switch (tcode)
				{
				case 30005:

					// Convection interpolation scheme
					GetLine(fp);
					fscanf(fp, "%d", &parameter.scheme);
				
					break;

				case 30020:

					// Binary output
					GetLine(fp);
					fscanf(fp, "%d", &parameter.wbinary);
		
					break;

				case 30040:

					// Calculate nth variable 
					GetLine(fp);
					for (j = 0; j < nphi; j++) 
						fscanf(fp, "%d", &parameter.calc[j]);
					
					break;

				case 30100:

					// Steady state
					GetLine(fp);
					fscanf(fp, "%d", &parameter.steady);

					break;

				case 30105:

					// Convergence criterion for steady state 			
					GetLine(fp);
					fscanf(fp, "%f", &parameter.ftol);

					break;

				case 30200:

					// Adjust time interval

					GetLine(fp);
					fscanf(fp, "%d", &parameter.adjdt);

					break;

				case 30201:

					// Maximum Courant number
					GetLine(fp);
					fscanf(fp, "%f", &parameter.maxCp);

					break;

				case 30400:

					// Number of saves
					GetLine(fp);
					fscanf(fp, "%d", &parameter.nsav);
	
					break;

				case 30450:

					// Write face scalars (u v w p T s)
					GetLine(fp);
					for (j = 0; j < nphi; j++) 
						fscanf(fp, "%d", &parameter.fsav[j]);

					break;

				case 30455:

					//Write face vectors (uvw)
					GetLine(fp);
					fscanf(fp, "%d", &parameter.fvec);

					break;

				case 30460:

					// Write element scalars (u v w p T s)
					GetLine(fp);					
					for (j = 0; j < nphi; j++) 
						fscanf(fp, "%d", &parameter.csav[j]);

					break;

				case 30465:
	
					// Write element vectors (uvw)
					GetLine(fp);
					fscanf(fp, "%d", &parameter.cvec);
		
					break;

				case 30470:

					// Write vorticity (x y z)
					GetLine(fp);
					for (j = 0; j < 3; j++) 
						fscanf(fp, "%d", &parameter.vortex[j]);

					break;

				case 30475:

					// Write stream function (xy)
					GetLine(fp);
					fscanf(fp, "%d", &parameter.streamf);

					break;
									
				case 30485:
					
					// Probe (u v w p T s)
					GetLine(fp);
					for (j = 0; j < nphi; j++) 
						fscanf(fp, "%d", &parameter.probe[j]);

					break;
				
				case 30550:

					// Maximum number of non-othorgonal corrections
					GetLine(fp);
					fscanf(fp, "%d", &parameter.northocor);

					break;

				case 30600:
				
					// Convergence criterion (matrix solution)					
					GetLine(fp);
					fscanf(fp, "%f", &parameter.mtol);

					break;

				case 30601:

					// Maximum number of iterations (matrix solution)
					GetLine(fp);
					fscanf(fp, "%d", &parameter.miter);
	
					break;

				case 30650:

					// Matrix solver (u v w p T s) (0-Jacobi, 1-SOR, 2-CGN, 3-GMRES, 4-BiCG, 5-QMR, 6-CGS, 7-BiCGStab, 8-BiCGStabM) 
					GetLine(fp);
					for (j = 0; j < nphi; j++) 
						fscanf(fp, "%d", &parameter.msolver[j]);

					break;

				case 30651:

					// Matrix preconditioner (0-Null, 1-Jacobi, 2-SOR, 3-ILU)
					GetLine(fp);			
					for (j = 0; j < nphi; j++) 
						fscanf(fp, "%d", &parameter.mprecond[j]);

					break;

				case 30800:

					// Interface scheme factor - CICSAM
					GetLine(fp);
					fscanf(fp, "%f", &parameter.kq);

					break;

				case 30900:

					// Maximum number of CICSAM corrections

					GetLine(fp);
					fscanf(fp, "%d", &parameter.ncicsamcor);

					break;
				
				case 32000:

					// Start time
					GetLine(fp);
					fscanf(fp, "%f", &parameter.t0);

					break;

				case 32001:

					// End time
					GetLine(fp);
					fscanf(fp, "%f", &parameter.t1);

					break;

				case 32002:

					// Time interval
					GetLine(fp);
					fscanf(fp, "%f", &parameter.dt);

					break;

				case 34000:

					// Gravity vector
					GetLine(fp);
					fscanf(fp, "%f %f %f", &parameter.g[0], &parameter.g[1], &parameter.g[2]);

					break;

				default:

					printf("\nError: Unknown parameter code.\n");
					exit(LOGICAL_ERROR);
					break;

				}
			
			}

		}

	} while (!feof(fp));

	// Convection interpolation scheme
	printf("\n");	
	if (parameter.scheme == 0)
		printf("Use upwind differencing scheme (UDS) for interpolation of convection term.\n");
	if (parameter.scheme == 1)
		printf("Use central differencing scheme (CDS) for interpolation of convection term.\n");
	
	// Steady-state
	printf("\n");	
	if (parameter.steady == LOGICAL_TRUE)
		printf("Calculate until steady state is reached.\n");	
	else
		printf("Stop simulation at end time.\n");	
	printf("Convergence criterion for steady state: \t%+.3E\n", parameter.ftol);

	// Calculate variables & save & probe parameters 

	printf("\n");	
	if (parameter.wbinary == LOGICAL_TRUE)		
		printf("Save results in binary format.\n");
	else
		printf("Save results in ascii format.\n");
	
	printf("Number of saves: \t\t\t\t%d\n", parameter.nsav);
		
	printf("\n");	
	printf("Variable: \t\t\t\t\t[ u v w p T s]\n");
	
	printf("Calculate: \t\t\t\t\t[");
	for (j = 0; j < nphi; j++) printf(" %d", parameter.calc[j]);
	printf("]\n");
	
	printf("Save scalars on face: \t\t\t\t[");
	for (j = 0; j < nphi; j++) printf(" %d", parameter.fsav[j]);
	printf("]\n");
	
	printf("Save scalars in cell: \t\t\t\t[");
	for (j = 0; j < nphi; j++) printf(" %d", parameter.csav[j]);
	printf("]\n");
	
	printf("Probe options: \t\t\t\t\t[");
	for (j = 0; j < nphi; j++) printf(" %d", parameter.probe[j]);
	printf("]\n");

	printf("\n");	
	printf("Axis: \t\t\t\t\t\t[ x y z]\n");
	
	printf("Save vorticity: \t\t\t\t[");
	for (j = 0; j < 3; j++) printf(" %d", parameter.vortex[j]);
	printf("]\n");
			
	printf("\n");	
	if (parameter.fvec == LOGICAL_TRUE)
		printf("Save face - vector magnitude: \t\t\t[yes]\n");
	else
		printf("Save face - vector magnitude: \t\t\t[no]\n");
	
	if (parameter.cvec == LOGICAL_TRUE)
		printf("Save cell center - vector: \t\t\t[yes]\n");
	else
		printf("Save cell center - vector: \t\t\t[no]\n");
	
	// VOF		
	printf("\n");	
	printf("Interface scheme factor: \t\t\t%.3f\n", parameter.kq);	
	printf("Number of CICSAM corrections: \t\t\t%d\n", parameter.ncicsamcor);	
	printf("Gravity vector (x): \t\t\t\t%+.3E m/s^2\n", parameter.g[0]);
	printf("Gravity vector (y): \t\t\t\t%+.3E m/s^2\n", parameter.g[1]);
	printf("Gravity vector (z): \t\t\t\t%+.3E m/s^2\n", parameter.g[2]);
	
	// Matrix solver		
	printf("\n");	
	printf("Variable: \t\t\t\t\t[ u v w p T s]\n");
	printf("Solver: a) \t\t\t\t\t[");
	for (j = 0; j < nphi; j++) printf(" %d", parameter.msolver[j]);
	printf("]\n");
	printf("Pre-conditioner: b) \t\t\t\t[");
	for (j = 0; j < nphi; j++) printf(" %d", parameter.mprecond[j]);
	printf("]\n");
	
	printf("Convergence criterion of matrix solution: \t%+.3E\n", parameter.mtol);	
	printf("Maximum number of matrix iterations: \t\t%d\n", parameter.miter);	
	printf("Maximum number of non-orthogonal corrections: \t%d\n", parameter.northocor);

	// Adjust time interval	
	printf("\n");	
	if (parameter.adjdt == LOGICAL_TRUE)
		printf("Adjust time intervals.\n");
	else
		printf("Do not adjust time intervals.\n");
	
	printf("Maximum Courant number: \t\t\t%.3f\n", parameter.maxCp);
			
	// Time data							
	printf("\n");	
	printf("Time interpolation method: \t\t\tIMPLICIT EULER\n");
	printf("Start time: \t\t\t\t\t%+.3E s\n", parameter.t0);
	printf("End time: \t\t\t\t\t%+.3E s\n", parameter.t1);
	printf("Time interval: \t\t\t\t\t%+.3E s\n", parameter.dt);

	printf("\n");	
	printf("a) Solvers list\n");
	printf("0-Jacobi\n");
	printf("1-SOR\n");
	printf("2-CGN\n");
	printf("3-GMRES\n");
	printf("4-BiCG\n");
	printf("5-QMR\n");
	printf("6-CGS\n");
	printf("7-BiCGStab\n");
	printf("8-BiCGStabM\n");
	printf("\n");	
	printf("b) Pre-conditioners list\n");
	printf("0-Null\n");
	printf("1-Jacobi\n");
	printf("2-SOR\n");
	printf("3-ILU\n");
	printf("\n");	

	printf("Done.\n");
		
	return LOGICAL_TRUE;

}

