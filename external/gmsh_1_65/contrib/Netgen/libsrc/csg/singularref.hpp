#ifndef FILE_SINGULARREF
#define FILE_SINGULARREF

/**************************************************************************/
/* File:   singularref.hh                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Sep. 99                                                    */
/**************************************************************************/

/**
   Control for local refinement
*/



/**
   Singular Face.
   Causes a bounday layer mesh refinement.
   All elements in subdomain domnr will get a boundary layer
   on faces sharing the solid sol
 */
class SingularFace 
{
public:
  int domnr;
  const Solid *sol;
  // ARRAY<Point<3> > points;
  // ARRAY<INDEX_2> segms;
public:
  SingularFace (int adomnr, const Solid * asol)
    : domnr(adomnr), sol(asol) { ; }
  const Solid * GetSolid() const { return sol; }
  int GetDomainNr () const { return domnr; }
};


///
class SingularEdge 
{
public:
  double beta;
  const Solid *sol1, *sol2;
  ARRAY<Point<3> > points;
  ARRAY<INDEX_2> segms;
public:
  SingularEdge (double abeta, const Solid * asol1, const Solid * asol2);
  void FindPointsOnEdge (class Mesh & mesh);
  void SetMeshSize (class Mesh & mesh, double globalh);
};


///
class SingularPoint
{
public:
  double beta;
  const Solid *sol1, *sol2, *sol3;
  ARRAY<Point<3> > points;
 
public:
  SingularPoint (double abeta, const Solid * asol1, const Solid * asol2,
		 const Solid * asol3);
  void FindPoints (class Mesh & mesh);
  void SetMeshSize (class Mesh & mesh, double globalh);
};


#endif
