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

enum
{
	UPWIND = 0,
	CDS

}par_scheme;

typedef struct
{
	
	int phi;
	
	int fixed;
	int along;

	float value;

}par_probe;

typedef struct
{

	int wbinary;

	int steady;
	int scheme;
	int adjdt;

	float maxCp;

	float mtol;
	int miter;

	int   northocor;
	float ftol;

	float kq;
	int ncicsamcor;

	int nsav;
	
	int smooth;
	
	int calc[6];

	int fsav[6];
	int csav[6];

	int fvec;
	int cvec;

	int vortex[3];
	int streamf;
	
	float t0, t1, dt;
	float g[3];

	int probe[6];
	
	int msolver[6];
	int mprecond[6];
	
}par_parameter;

par_parameter parameter;

int ParImportPAR(char *file);
