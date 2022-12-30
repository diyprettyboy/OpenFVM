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

#define CONSTANT  0

#define IDEALGAS  1
#define TAIT      2

#define CROSSWLF  1

#define GAS       0
#define LIQUID    1

typedef struct
{

	float constant;

}mat_psi;

typedef struct
{

	float constant;

}mat_dens;

typedef struct
{
	float constant;

}mat_visc;

typedef struct
{
	float spheat;
	float thcond;
	
}mat_therm;

typedef struct
{

	mat_psi psi[2];
	mat_dens dens[2];
	mat_visc visc[2];
	mat_therm therm[2];

	float tens;
    float bthcond;

}mat_material;

mat_material material;

int MtlImportMTL(char *file);

