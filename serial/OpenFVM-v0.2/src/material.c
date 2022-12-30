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
#include "material.h"

int MtlImportMTL(char *file)
{

	int i, n;

	int inull;

	int tcode;
	int nbmat;

	FILE *fp;
	char descr[512];

	fp = fopen(file, "r");

	if (fp == NULL)
	{
		printf("\nError: Material file not found!\n");
		printf("%s\n\n", file);
		exit(LOGICAL_ERROR);
	}

	SetDefaults();

	printf("\nReading materials...\n", file);

	do
	{

		do
		{
			fscanf(fp, "%s", descr);

			if (strcmp(descr, "$MTRL") == 0)
				break;

			if (strcmp(descr, "$ENDF") == 0)
				break;

		} while (!feof(fp));

		if (strcmp(descr, "$ENDF") == 0)
			break;

		if (strcmp(descr, "$MTRL") == 0)
		{

			fscanf(fp, "%d %d", &inull, &nbmat);

			GetLine(fp); 
	
			for (i = 0; i < nbmat; i++)
			{

				fscanf(fp, "%d %d", &tcode, &n);

				switch (tcode)
				{
				case 20012:

					// Compressibility of fluid 0
					GetLine(fp);
					fscanf(fp, "%f", &material.psi[0].constant);

					break;

				case 20015:

					// Density of fluid 0
					GetLine(fp);
					fscanf(fp, "%f", &material.dens[0].constant);

					break;

				case 20021:

					// Viscosity of fluid 0
					GetLine(fp);
					fscanf(fp, "%f", &material.visc[0].constant);
					
					break;

				case 20022:

					// Specific heat of fluid 0
					GetLine(fp);
					fscanf(fp, "%f", &material.therm[0].spheat);

					break;

				case 20023:

					// Thermal conductivity of fluid 0
					GetLine(fp);
					fscanf(fp, "%f", &material.therm[0].thcond);

					break;

				case 21012:

					// Compressibility of fluid 1
					GetLine(fp);
					fscanf(fp, "%f", &material.psi[1].constant);

					break;

				case 21015:

					// Density of fluid 1
					GetLine(fp);
					fscanf(fp, "%f", &material.dens[1].constant);

					break;

				case 21021:

					// Viscosity of fluid 1
					GetLine(fp);
					fscanf(fp, "%f", &material.visc[1].constant);
					
					break;

				case 21022:

					// Specific heat of fluid 1
					GetLine(fp);
					fscanf(fp, "%f", &material.therm[1].spheat);

					break;

				case 21023:

					// Thermal conductivity of fluid 1
					GetLine(fp);
					fscanf(fp, "%f", &material.therm[1].thcond);

					break;

				case 26050:

					// Constant surface tension (fluid 0 / fluid 1)
					GetLine(fp);
					fscanf(fp, "%f", &material.tens);

					break;

				case 26060:

					// Thermal conductivity (boundary)
					GetLine(fp);
					fscanf(fp, "%f", &material.bthcond);

					break;

				default:

					printf("\nError: Unknown material code.\n");
					exit(LOGICAL_ERROR);
					break;

				}

			}

		}

	} while (!feof(fp));

	printf("\n");
	printf("- Fluid 0\n");
	printf("Density: \t\t\t\t\t%+.3E kg/m^3\n", material.dens[0].constant);
	printf("Viscosity: \t\t\t\t\t%+.3E kg/(m.s)\n", material.visc[0].constant);
	printf("Specific heat: \t\t\t\t\t%+.3E J/(kg.K)\n", material.therm[0].spheat);
	printf("Thermal conductivity: \t\t\t\t%+.3E W/(m.K)\n", material.therm[0].thcond);

	printf("\n");
	printf("- Fluid 1\n");
	printf("Density: \t\t\t\t\t%+.3E kg/m^3\n", material.dens[1].constant);
	printf("Viscosity: \t\t\t\t\t%+.3E kg(m.s)\n", material.visc[1].constant);
	printf("Specific heat: \t\t\t\t\t%+.3E J/(kg.K)\n", material.therm[1].spheat);
	printf("Thermal conductivity: \t\t\t\t%+.3E W/(m.K)\n", material.therm[1].thcond);

	printf("\n");
	printf("- Fluid 0 / Fluid 1\n");
	printf("Surace tension: \t\t\t\t%+.3E kg(m.s^2)\n", material.tens);

	printf("\n");
	printf("- Boundary\n");
	printf("Thermal conductivity: \t\t\t\t%+.3E W/(m.K)\n", material.bthcond);

	printf("Done.\n");
	
	return LOGICAL_TRUE;

}

