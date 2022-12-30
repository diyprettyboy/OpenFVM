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

#include "mesh.h"
#include "geocalc.h"
#include "globals.h"

#include "laspack/itersolv.h"
#include "laspack/rtc.h"
#include "laspack/errhandl.h"

msh_vector Gradient(Vector *phi, Vector *phif, int bound, int element)
{

	int j;

	int neighbor, face, pair;
	
	double phij;

	double dNf, dPf;
	double lambda;

	msh_vector rv;

	rv.x = 0.0;
	rv.y = 0.0;
	rv.z = 0.0;

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

			// Element face variable
			
			phij = V_GetCmp(phi, neighbor + 1) * lambda + V_GetCmp(phi, element + 1) * (1.0 - lambda);

			// Element center gradient

			rv.x += 1 / elements[element].Vp * phij * faces[face].A.x;
			rv.y += 1 / elements[element].Vp * phij * faces[face].A.y;
			rv.z += 1 / elements[element].Vp * phij * faces[face].A.z;

		}
		else
		{

			// Element face variable

			if (bound == LOGICAL_TRUE) 
				phij = V_GetCmp(phif, face + 1);
			else 
				phij = V_GetCmp(phi, element + 1);

			// Element center gradient

			rv.x += 1 / elements[element].Vp * phij * faces[face].A.x;
			rv.y += 1 / elements[element].Vp * phij * faces[face].A.y;
			rv.z += 1 / elements[element].Vp * phij * faces[face].A.z;

		}

	}

	return rv;

}
