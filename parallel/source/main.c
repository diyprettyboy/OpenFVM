/***************************************************************************
 *   Copyright (C) 2004-2006 by OpenCAE team                               *
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

#include "variables.h"
#include "vector.h"

#include "mesh.h"
#include "material.h"
#include "bcond.h"
#include "param.h"

#include "globals.h"
#include "geocalc.h"
#include "parser.h"
#include "post.h"
#include "decomp.h"
#include "reorder.h"
#include "parallel.h"
#include "restart.h"
#include "setup.h"
#include "solve.h"
#include "fill.h"

int
Simulation (char *path)
{

  int i, n;

  int d;

  char var[6];

  double fres[6];
  int fiter[6];

  int iter;

  int irestart;

  double starttime, endtime, curtime, dt;
  double wtime, wdt;

  double maxCp;

  double Vc, Vt;
  double Vcp, Vtp;

  double pf;

  char *file;

  FILE *fpresults;		// Output to gmsh post-processing file (results)
  FILE *fpprobe;		// Output to gnuplot file (probe)
  FILE *fpresiduals;		// Output to gnuplot file (residuals)
  FILE *fprestart;		// Input/Output for restart

  var[iu] = 'u';
  var[iv] = 'v';
  var[iw] = 'w';
  var[ip] = 'p';
  var[iT] = 'T';
  var[is] = 's';

  // Set ghosts
  SetGhosts ();

  // Allocate memory
  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Allocating memory...\n");

  AllocateMemory ();

  PetscPrintf (PETSC_COMM_WORLD, "Memory allocated.\n");

  // Set elements center
  SetCenters ();

  VecGhostGetLocalForm (cex, &cexl);
  VecGhostGetLocalForm (cey, &ceyl);
  VecGhostGetLocalForm (cez, &cezl);

  // Set initial conditions
  SetInitialConditions ();

  // Set initial flux
  SetInitialFlux ();

  // Set boundary velocity and pressure
  SetBoundary ();

  // Set material properties       
  SetMaterialProperties ();

  // Set time intervals
  starttime = parameter.t0;
  endtime = parameter.t1;
  curtime = parameter.t0;
  dt = parameter.dt;

  // Set restart intervals
  irestart = parameter.restart;

  // Allocate memory for file string
  file = calloc (strlen (path) + 9, sizeof (char));

  // Read restart file
  sprintf (file, "%s.%03d.ini", path, processor);

  fprestart = fopen (file, "rb");

  if (fprestart != NULL)
    {

      ReadRestart (fprestart, &curtime);

      endtime += curtime;

      fclose (fprestart);

    }

  if (parameter.steady == LOGICAL_TRUE)
    {

      sprintf (file, "%s.res", path);

      // Open output files for residuals
      PetscFOpen (PETSC_COMM_WORLD, file, "w", &fpresiduals);
    }

  sprintf (file, "%s.%03d.pos", path, processor);

  // Open output file for results
  fpresults = fopen (file, "w");

  d = sizeof (double);

  if (parameter.wbinary == LOGICAL_TRUE)
    {
      fprintf (fpresults, "$PostFormat\n");
      fprintf (fpresults, "%g %d %d\n", 1.3, 1, d);
      fprintf (fpresults, "$EndPostFormat\n");
    }
  else
    {
      fprintf (fpresults, "$PostFormat\n");
      fprintf (fpresults, "%g %d %d\n", 1.2, 0, d);
      fprintf (fpresults, "$EndPostFormat\n");
    }

  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

  // Set write time intervals
  if (parameter.nsav != 0)
    wdt = (endtime - starttime) / parameter.nsav;
  else
    wdt = 2 * (endtime - starttime);

  wtime = wdt;

  iter = 0;

  fiter[iu] = 0;
  fiter[iv] = 0;
  fiter[iw] = 0;
  fiter[ip] = 0;
  fiter[iT] = 0;
  fiter[is] = 0;

  if (parameter.fill == LOGICAL_TRUE)
    {
      Vtp = VolumeTotal ();

      Vc = 0.0;

      Vt = Vtp;

      if (processor == 0)
	{
	  for (i = 1; i < nbprocessors; i++)
	    {
	      MPI_Recv (&Vtp, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, 0);

	      Vt += Vtp;
	    }
	}
      else
	{
	  MPI_Send (&Vtp, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

      MPI_Barrier (MPI_COMM_WORLD);
    }

  do
    {

      curtime += dt;

      iter++;

      PetscPrintf (PETSC_COMM_WORLD, "\nTime = %.3E\n", curtime);

      Solve (var, fiter, dt, &maxCp, verbose, pchecks);

      if (parameter.fill == LOGICAL_TRUE)
	{

	  Vcp = VolumeFilled ();

	  Vc = Vcp;

	  if (processor == 0)
	    {

	      for (i = 1; i < nbprocessors; i++)
		{
		  MPI_Recv (&Vcp, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, 0);

		  Vc += Vcp;
		}

	    }
	  else
	    {
	      MPI_Send (&Vcp, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	    }

	  MPI_Barrier (MPI_COMM_WORLD);

	  if (processor == 0)
	    {
	      pf = Vc / Vt * 100;

	      for (i = 1; i < nbprocessors; i++)
		{
		  MPI_Send (&pf, 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
		}
	    }
	  else
	    {
	      MPI_Recv (&pf, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, 0);
	    }

	  MPI_Barrier (MPI_COMM_WORLD);

	  PetscPrintf (PETSC_COMM_WORLD, "\nPercentage filled: %.2f%%\n", pf);

	  if (pf > parameter.pf)
	    break;

	}

      if (parameter.adjdt == LOGICAL_TRUE)
	{

	  if (maxCp > parameter.maxCp)
	    dt *= 0.85 * parameter.maxCp / maxCp;

	  if (1.25 * maxCp < parameter.maxCp)
	    dt *= 1.05;

	}

      if (iter >= irestart)
	{

	  // Write restart file
	  sprintf (file, "%s.%03d.ini", path, processor);

	  fprestart = fopen (file, "wb");

	  if (fprestart != NULL)
	    {

	      WriteRestart (fprestart, curtime);

	    }

	  fclose (fprestart);

	  irestart += parameter.restart;

	}

      if (curtime > wtime + dt)
	{

	  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

	  //fflush (fpresults);

	  // Open output file for probes
	  sprintf (file, "%s.%03d.prb", path, processor);

	  fpprobe = fopen (file, "w");

	  WriteProbeViews (fpprobe, var, curtime);

	  fclose (fpprobe);

	  wtime += wdt;

	}

      if (parameter.steady == LOGICAL_TRUE)
	{

	  // Get residual

	  if (parameter.calc[iu] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xu, xu0);
	      VecNorm (temp1, NORM_2, &fres[iu]);

	    }
	  else
	    fres[iu] = 0.0;

	  if (parameter.calc[iv] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xv, xv0);
	      VecNorm (temp1, NORM_2, &fres[iv]);

	    }
	  else
	    fres[iv] = 0.0;

	  if (parameter.calc[iw] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xw, xw0);
	      VecNorm (temp1, NORM_2, &fres[iw]);

	    }
	  else
	    fres[iw] = 0.0;

	  if (parameter.calc[ip] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xp, xp0);
	      VecNorm (temp1, NORM_2, &fres[ip]);

	    }
	  else
	    fres[ip] = 0.0;

	  if (parameter.calc[iT] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xT, xT0);
	      VecNorm (temp1, NORM_2, &fres[iT]);

	    }
	  else
	    fres[iT] = 0.0;

	  if (parameter.calc[is] == LOGICAL_TRUE)
	    {

	      VecWAXPY (temp1, -1.0, xs, xs0);
	      VecNorm (temp1, NORM_2, &fres[is]);

	    }
	  else
	    fres[is] = 0.0;

	  n = 0;

	  for (i = 0; i < nphi; i++)
	    {

	      if (verbose == LOGICAL_TRUE)
		PetscPrintf (PETSC_COMM_WORLD,
			     "\nVariable: %c Iteration: %d Final residual: %+E\n",
			     var[i], fiter[i], fres[i]);

	      if (fres[i] > parameter.ftol[i])
		n++;

	    }

	  if (n == 0)
	    {
	      PetscPrintf (PETSC_COMM_WORLD, "\nSteady state reached.\n");
	      break;
	    }

	  WriteResidual (fpresiduals, iter, fres);

	  fflush (fpresiduals);

	}
      else
	{

	  if (curtime + 0.5 * dt > endtime)
	    break;

	}

    }
  while (dt > 0.0);

  WriteResults (fpresults, LOGICAL_TRUE, LOGICAL_TRUE, curtime);

  // Close output files

  fclose (fpresults);

  if (parameter.steady == LOGICAL_TRUE)
    {

      WriteResidual (fpresiduals, iter, fres);

      fflush (fpresiduals);

      PetscFClose (PETSC_COMM_WORLD, fpresiduals);

    }

  // Open output file for probes
  sprintf (file, "%s.%03d.prb", path, processor);

  fpprobe = fopen (file, "w");

  WriteProbeViews (fpprobe, var, curtime);

  fclose (fpprobe);

  VecGhostRestoreLocalForm (cex, &cexl);
  VecGhostRestoreLocalForm (cey, &ceyl);
  VecGhostRestoreLocalForm (cez, &cezl);

  // Release memory 
  free (file);

  DeallocateMemory ();

  return LOGICAL_TRUE;

}

void
usage ()
{

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Usage:     mpirun -np [n] ../OpenFVM [problem] [options] [n] [petsc-options]\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Example 1: mpirun -np [n] ../OpenFVM lid f 1\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Example 2: mpirun -np [n] ../OpenFVM lid f 1 > lid.log\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "Example 3: mpirun -np [n] ../OpenFVM lid d 4\n");

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "Main options:\n");
  PetscPrintf (PETSC_COMM_WORLD, "  r - Reorder the mesh\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "  d - Decompose the mesh into n partitions/processors\n");
  PetscPrintf (PETSC_COMM_WORLD,
	       "  f - Start simulation with n partitions/processors\n");
  PetscPrintf (PETSC_COMM_WORLD, "\n");

  PetscPrintf (PETSC_COMM_WORLD, "Additional options:\n");
  PetscPrintf (PETSC_COMM_WORLD, "  v - Verbose mode\n");
  PetscPrintf (PETSC_COMM_WORLD, "  c - Perform checks\n");
  PetscPrintf (PETSC_COMM_WORLD, "\n");

}

int
main (int argc, char **argv)
{

  char *path;
  char *file;

  double start, end;
  double elapsed;

  static char help[] =
    "Three-dimensional unstructured finite-volume implicit flow solver.\n";

  // Start MPI
  PetscInitialize (&argc, &argv, (char *) 0, help);

  MPI_Comm_size (PETSC_COMM_WORLD, &nbprocessors);
  MPI_Comm_rank (PETSC_COMM_WORLD, &processor);

  //PetscPrintf (PETSC_COMM_WORLD, "Processor: %d\n", processor);

  PetscPrintf (PETSC_COMM_WORLD, "\n");
  PetscPrintf (PETSC_COMM_WORLD, "*****************************************\n");
  PetscPrintf (PETSC_COMM_WORLD, "*                                       *\n");
  PetscPrintf (PETSC_COMM_WORLD, "*    OpenFVM-Flow v1.0 - Parallel       *\n");
  PetscPrintf (PETSC_COMM_WORLD, "*                                       *\n");
  PetscPrintf (PETSC_COMM_WORLD, "*****************************************\n");
  PetscPrintf (PETSC_COMM_WORLD, "\n");

  if (argc < 4)
    {
      usage ();
      PetscPrintf (PETSC_COMM_WORLD,
		   "\nError: Wrong number of arguments, expected %d found %d.\n\n",
		   4 - 1, argc - 1);

      // end MPI
      PetscFinalize ();

      return LOGICAL_ERROR;
    }

  if (strchr (argv[2], 'c') != NULL)
    pchecks = LOGICAL_TRUE;
  else
    pchecks = LOGICAL_FALSE;

  if (strchr (argv[2], 'v') != NULL)
    verbose = LOGICAL_TRUE;
  else
    verbose = LOGICAL_FALSE;

  path = calloc (strlen (argv[1]) + 1, sizeof (char));
  file = calloc (strlen (argv[1]) + 12, sizeof (char));

  // Simulate
  if (strchr (argv[2], 'f') != NULL)
    {

      if (nbprocessors < 1)
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Number of processors cannot be lower than 1.\n\n");

	  // end MPI
	  PetscFinalize ();
	  
	  return LOGICAL_ERROR;
	}

      if (atoi (argv[3]) < 1)
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Number of partitions cannot be lower than 1.\n\n");

	  // end MPI
	  PetscFinalize ();

	  return LOGICAL_ERROR;
	}

      if (nbprocessors != atoi (argv[3]))
	{
	  PetscPrintf (PETSC_COMM_WORLD,
		       "\nError: Number of processors and partitions are different.\n\n");

	  // end MPI
	  PetscFinalize ();

	  return LOGICAL_ERROR;
	}

      strcpy (path, argv[1]);

      // Read parameter file
      sprintf (file, "%s.par", path);
      ParImportPAR (file);

      // Read mesh file
      strcpy (file, path);
      sprintf (file, "%s.%03d.msh", path, processor);
      MshImportMSH (file);

      // Read boundary conditions file
      sprintf (file, "%s.bcd", path);
      BcdImportBCD (file);

      // Read material file
      sprintf (file, "%s.mtl", path);
      MtlImportMTL (file);

      MPI_Barrier (MPI_COMM_WORLD);

      // Start clock to measure simulation time
      start = MPI_Wtime ();

      // Start simulation
      Simulation (path);

      end = MPI_Wtime ();
      elapsed = end - start;
      printf ("\nProcessor: %d, Elapsed time: %.2f seconds.\n", processor,
	      elapsed);

    }

  // Decompose
  if (strchr (argv[2], 'd') != NULL)
    {

      strcpy (path, argv[1]);

      // Read mesh file
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Mesh file not found!\n");
	  PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      // Read boundary conditions file
      sprintf (file, "%s.bcd", path);
      BcdImportBCD (file);

      // Allocate memory
      PetscPrintf (PETSC_COMM_WORLD, "\n");
      PetscPrintf (PETSC_COMM_WORLD, "Allocating memory...\n");

      AllocateMemory ();

      PetscPrintf (PETSC_COMM_WORLD, "Memory allocated.\n");

      SetBoundary ();

      DeallocateMemory ();

      DecomposeMesh (path, atoi (argv[3]));

    }

  // Reorder
  if (strchr (argv[2], 'r') != NULL)
    {

      strcpy (path, argv[1]);

      // Read mesh file
      sprintf (file, "%s.msh", path);

      if (MshImportMSH (file) == LOGICAL_ERROR)
	{
	  PetscPrintf (PETSC_COMM_WORLD, "\nError: Mesh file not found!\n");
	  PetscPrintf (PETSC_COMM_WORLD, "%s\n\n", file);
	  return LOGICAL_ERROR;
	}

      ReorderMesh (path);

    }

  free (path);
  free (file);

  // Free memory          
  MshFreeMemory ();

  MPI_Barrier (MPI_COMM_WORLD);

  PetscPrintf (PETSC_COMM_WORLD, "Done.\n\n");

  // end MPI
  PetscFinalize ();

  return 0;
}
