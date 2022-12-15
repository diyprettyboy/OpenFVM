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

#include "petscksp.h"

#include "itersolv.h"

void
GMRESIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	   double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPSetType (ksp, KSPFGMRES);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}

void
CRIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	 double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPCR);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}
void
CGIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPCG);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);
    
}

void
CGSIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter, double maxtol,
	 double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPCGS);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}

void
BiCGIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter,
	      double maxtol, double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPBICG);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}

void
BiCGSTABIter (Mat * A, Vec * x, Vec * b, int maxiter, int *iter,
	      double maxtol, double *res, PCType PrecondProc)
{

  KSP ksp;
  PC pc;

  KSPCreate (PETSC_COMM_WORLD, &ksp);

  KSPSetOperators (ksp, *A, *A, DIFFERENT_NONZERO_PATTERN);

  KSPGetPC (ksp, &pc);

  PCSetType (pc, PrecondProc);

  KSPSetType (ksp, KSPBCGS);

  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  KSPSetTolerances (ksp, 1E-50, maxtol, PETSC_DEFAULT, maxiter);

  KSPSetFromOptions (ksp);

  KSPSolve (ksp, *b, *x);

  KSPGetIterationNumber (ksp, iter);

  KSPGetResidualNorm (ksp, res);

  KSPDestroy (ksp);

}
