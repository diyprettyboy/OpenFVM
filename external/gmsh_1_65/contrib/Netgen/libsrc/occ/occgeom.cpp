#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <occgeom.hpp>  
#include "ShapeAnalysis_ShapeTolerance.hxx"
#include "ShapeAnalysis_ShapeContents.hxx"
#include "ShapeAnalysis_CheckSmallFace.hxx"
#include "ShapeAnalysis_DataMapOfShapeListOfReal.hxx"
#include "BRepAlgoAPI_Fuse.hxx"
#include "BRepCheck_Analyzer.hxx"
#include "BRepLib.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "ShapeFix.hxx"
#include "ShapeFix_FixSmallFace.hxx"


namespace netgen
{

void OCCGeometry :: PrintNrShapes ()
{
  TopExp_Explorer e;
  int count = 0;
  for (e.Init(shape, TopAbs_COMPSOLID); e.More(); e.Next()) count++;
  cout << "CompSolids: " << count << endl;

  cout << "Solids    : " << somap.Extent() << endl;
  cout << "Shells    : " << shmap.Extent() << endl;
  cout << "Faces     : " << fmap.Extent() << endl;
  cout << "Edges     : " << emap.Extent() << endl;
  cout << "Vertices  : " << vmap.Extent() << endl;
}


void PrintContents (OCCGeometry * geom)
{
  ShapeAnalysis_ShapeContents cont;
  cont.Clear();
  cont.Perform(geom->shape);

  (*testout) << "OCC CONTENTS" << endl;
  (*testout) << "============" << endl;
  (*testout) << "SOLIDS   : " << cont.NbSolids() << endl;
  (*testout) << "SHELLS   : " << cont.NbShells() << endl;
  (*testout) << "FACES    : " << cont.NbFaces() << endl;
  (*testout) << "WIRES    : " << cont.NbWires() << endl;
  (*testout) << "EDGES    : " << cont.NbEdges() << endl;
  (*testout) << "VERTICES : " << cont.NbVertices() << endl;

  TopExp_Explorer e;
  int count = 0;
  for (e.Init(geom->shape, TopAbs_COMPOUND); e.More(); e.Next())
    count++;
  (*testout) << "Compounds: " << count << endl;

  count = 0;
  for (e.Init(geom->shape, TopAbs_COMPSOLID); e.More(); e.Next())
    count++;
  (*testout) << "CompSolids: " << count << endl;

  (*testout) << endl;

  cout << "Highest entry in topology hierarchy: " << endl;
  if (count)
    cout << count << " composite solid(s)" << endl;
  else
    if (geom->somap.Extent())
      cout << geom->somap.Extent() << " solid(s)" << endl;
    else
      if (geom->shmap.Extent())
	cout << geom->shmap.Extent() << " shells(s)" << endl;
      else
	if (geom->fmap.Extent())
	  cout << geom->fmap.Extent() << " face(s)" << endl;
	else
	  if (geom->wmap.Extent())
	    cout << geom->wmap.Extent() << " wire(s)" << endl;
	  else
	    if (geom->emap.Extent())
	      cout << geom->emap.Extent() << " edge(s)" << endl;
	    else
	      if (geom->vmap.Extent())
		cout << geom->vmap.Extent() << " vertices(s)" << endl;
	      else
		cout << "no entities" << endl;

}



void OCCGeometry :: HealGeometry ()
{
  int nrc = 0, nrcs = 0,
    nrso = somap.Extent(),
    nrsh = shmap.Extent(),
    nrf = fmap.Extent(),
    nrw = wmap.Extent(),
    nre = emap.Extent(),
    nrv = vmap.Extent();

  TopExp_Explorer e;
  for (e.Init(shape, TopAbs_COMPOUND); e.More(); e.Next()) nrc++;
  for (e.Init(shape, TopAbs_COMPSOLID); e.More(); e.Next()) nrcs++;

  double surfacecont = 0;

  for (int i = 1; i <= fmap.Extent(); i++)
    {
      GProp_GProps system;
      BRepGProp::LinearProperties(fmap(i), system);
      surfacecont += system.Mass();
    }

  cout << "Starting geometry healing procedure (tolerance: " << tolerance << ")" << endl
       << "-----------------------------------" << endl;

  if (fixsmalledges)
    {
      cout << endl << "- fixing small edges" << endl;

      Handle(ShapeFix_Wire) sfw;
      Handle_ShapeBuild_ReShape rebuild = new ShapeBuild_ReShape;
      rebuild->Apply(shape);

      for (int i = 1; i <= fmap.Extent(); i++)
	{
	  TopExp_Explorer exp1;
	  for (exp1.Init (fmap(i), TopAbs_WIRE); exp1.More(); exp1.Next())
	    {
	      TopoDS_Wire oldwire = TopoDS::Wire(exp1.Current());
	      sfw = new ShapeFix_Wire (oldwire, TopoDS::Face(fmap(i)),tolerance);
	      sfw->ModifyTopologyMode() = Standard_True;

	      if (sfw->FixSmall (false, tolerance))
		{
		  cout << "Fixed small edge in wire " << wmap.FindIndex (oldwire) << endl;
		  TopoDS_Wire newwire = sfw->Wire();
		  rebuild->Replace(oldwire, newwire, Standard_False);
		}
	      if ((sfw->StatusSmall(ShapeExtend_FAIL1)) ||
		  (sfw->StatusSmall(ShapeExtend_FAIL2)) ||
		  (sfw->StatusSmall(ShapeExtend_FAIL3)))
		cout << "Failed to fix small edge in wire " << wmap.FindIndex (oldwire) << endl;

	      
	    }
	}

      shape = rebuild->Apply(shape);



      {
      Handle_ShapeBuild_ReShape rebuild = new ShapeBuild_ReShape;
      rebuild->Apply(shape);
      TopExp_Explorer exp1;
      for (exp1.Init (shape, TopAbs_EDGE); exp1.More(); exp1.Next())
	{
	  TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
	  if (vmap.FindIndex(TopExp::FirstVertex (edge)) == 
	      vmap.FindIndex(TopExp::LastVertex (edge)))
	    {
	      GProp_GProps system;
	      BRepGProp::LinearProperties(edge, system);
	      if (system.Mass() < tolerance)
		{
		  cout << "removing degenerated edge " << emap.FindIndex(edge) << endl;
		  rebuild->Remove(edge, false);
		}
	    }
	}
      shape = rebuild->Apply(shape);
      }


      Handle(ShapeFix_Wireframe) sfwf = new ShapeFix_Wireframe;
      sfwf->SetPrecision(tolerance);
      sfwf->Load (shape);

      if (sfwf->FixSmallEdges())
	{
	  cout << endl << "- fixing wire frames" << endl;  
	  if (sfwf->StatusSmallEdges(ShapeExtend_OK)) cout << "no small edges found" << endl;
	  if (sfwf->StatusSmallEdges(ShapeExtend_DONE1)) cout << "some small edges fixed" << endl;
	  if (sfwf->StatusSmallEdges(ShapeExtend_FAIL1)) cout << "failed to fix some small edges" << endl;
	}
  

      if (sfwf->FixWireGaps())
	{
	  cout << endl << "- fixing wire gaps" << endl;
	  if (sfwf->StatusWireGaps(ShapeExtend_OK)) cout << "no gaps found" << endl;
	  if (sfwf->StatusWireGaps(ShapeExtend_DONE1)) cout << "some 2D gaps fixed" << endl;
	  if (sfwf->StatusWireGaps(ShapeExtend_DONE2)) cout << "some 3D gaps fixed" << endl;
	  if (sfwf->StatusWireGaps(ShapeExtend_FAIL1)) cout << "failed to fix some 2D gaps" << endl;
	  if (sfwf->StatusWireGaps(ShapeExtend_FAIL2)) cout << "failed to fix some 3D gaps" << endl;
	}
      

      shape = sfwf->Shape();
    }





  if (fixspotstripfaces)
    {
  
      cout << endl << "- fixing spot and strip faces" << endl;
      Handle(ShapeFix_FixSmallFace) sffsm = new ShapeFix_FixSmallFace();
      sffsm -> Init (shape);
      sffsm -> SetPrecision (tolerance);
      sffsm -> Perform();
      
      shape = sffsm -> FixShape();
    }

  if (sewfaces)
    {
      cout << endl << "- sewing faces" << endl;

      TopExp_Explorer exp0;

      BRepOffsetAPI_Sewing sewedObj(tolerance);

      for (exp0.Init (shape, TopAbs_FACE); exp0.More(); exp0.Next())
	{
	  TopoDS_Face face = TopoDS::Face (exp0.Current());
	  sewedObj.Add (face);
	}
      
      sewedObj.Perform();
  
      if (!sewedObj.SewedShape().IsNull())
	shape = sewedObj.SewedShape();
      else
	cout << " not possible";
    }

  if (makesolids)
    {  
      cout << endl << "- making solids" << endl;
      
      TopExp_Explorer exp0;

      BRepBuilderAPI_MakeSolid ms;
      int count = 0;
      for (exp0.Init(shape, TopAbs_SHELL); exp0.More(); exp0.Next())
	{
	  count++;
	  ms.Add (TopoDS::Shell(exp0.Current()));
	}
      
      if (!count)
	{
	  cout << " not possible (no shells)" << endl;
	}
      else
	{
	  BRepCheck_Analyzer ba(ms);
	  if (ba.IsValid ())
	    {
	      Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
	      sfs->Init (ms);
	      sfs->SetPrecision(tolerance);
	      sfs->SetMaxTolerance(tolerance);
	      sfs->Perform();
	      shape = sfs->Shape();
	      
	      for (exp0.Init(shape, TopAbs_SOLID); exp0.More(); exp0.Next())
		{
		  TopoDS_Solid solid = TopoDS::Solid(exp0.Current());
		  TopoDS_Solid newsolid = solid;
		  BRepLib::OrientClosedSolid (newsolid);
		  Handle_ShapeBuild_ReShape rebuild = new ShapeBuild_ReShape;
		  //		  rebuild->Apply(shape);
		  rebuild->Replace(solid, newsolid, Standard_False);
		  TopoDS_Shape newshape = rebuild->Apply(shape, TopAbs_COMPSOLID, 1);
		  //		  TopoDS_Shape newshape = rebuild->Apply(shape);
		  shape = newshape;
		}
	    }
	  else
	    cout << " not possible" << endl;
	}
    }

  BuildFMap();

  double newsurfacecont = 0;

  for (int i = 1; i <= fmap.Extent(); i++)
    {
      GProp_GProps system;
      BRepGProp::LinearProperties(fmap(i), system);
      newsurfacecont += system.Mass();
    }

  int nnrc = 0, nnrcs = 0,
    nnrso = somap.Extent(),
    nnrsh = shmap.Extent(),
    nnrf = fmap.Extent(),
    nnrw = wmap.Extent(),
    nnre = emap.Extent(),
    nnrv = vmap.Extent();

  for (e.Init(shape, TopAbs_COMPOUND); e.More(); e.Next()) nnrc++;
  for (e.Init(shape, TopAbs_COMPSOLID); e.More(); e.Next()) nnrcs++;

  cout << "-----------------------------------" << endl;
  cout << "Compounds       : " << nnrc << " (" << nrc << ")" << endl;
  cout << "Composite solids: " << nnrcs << " (" << nrcs << ")" << endl;
  cout << "Solids          : " << nnrso << " (" << nrso << ")" << endl;
  cout << "Shells          : " << nnrsh << " (" << nrsh << ")" << endl;
  cout << "Wires           : " << nnrw << " (" << nrw << ")" << endl;
  cout << "Faces           : " << nnrf << " (" << nrf << ")" << endl;
  cout << "Edges           : " << nnre << " (" << nre << ")" << endl;
  cout << "Vertices        : " << nnrv << " (" << nrv << ")" << endl;
  cout << endl;
  cout << "Totol surface area : " << newsurfacecont << " (" << surfacecont << ")" << endl;
  cout << endl;

}
 



void OCCGeometry :: BuildFMap()
{
  somap.Clear();
  shmap.Clear();
  fmap.Clear();
  wmap.Clear();
  emap.Clear();
  vmap.Clear();
  
  TopExp_Explorer exp0, exp1, exp2, exp3, exp4, exp5;
  
  for (exp0.Init(shape, TopAbs_SOLID);
       exp0.More(); exp0.Next())
    {
      TopoDS_Solid solid = TopoDS::Solid (exp0.Current());
      
      if (somap.FindIndex(TopoDS::Solid (exp0.Current())) < 1)
	{
	  somap.Add (TopoDS::Solid (exp0.Current()));
	  
	  for (exp1.Init(exp0.Current(), TopAbs_SHELL);
	       exp1.More(); exp1.Next())
	    {
	      TopoDS_Shell shell = TopoDS::Shell (exp1.Current().Composed (exp0.Current().Orientation()));
	      if (shmap.FindIndex(shell) < 1)
		{
		  shmap.Add (shell);
		  
		  for (exp2.Init(shell, TopAbs_FACE);
		       exp2.More(); exp2.Next())
		    {
		      TopoDS_Face face = TopoDS::Face(exp2.Current().Composed(shell.Orientation()));
		      if (fmap.FindIndex(face) < 1)
			{
			  fmap.Add (face);
			  
			  for (exp3.Init(exp2.Current(), TopAbs_WIRE);
			       exp3.More(); exp3.Next())
			    {
			      TopoDS_Wire wire = TopoDS::Wire (exp3.Current().Composed(face.Orientation()));
			      if (wmap.FindIndex(wire) < 1)
				{
				  wmap.Add (wire);
				  
				  for (exp4.Init(exp3.Current(), TopAbs_EDGE);
				       exp4.More(); exp4.Next())
				    {
				      TopoDS_Edge edge = TopoDS::Edge(exp4.Current().Composed(wire.Orientation()));
				      if (emap.FindIndex(edge) < 1)
					{
					  emap.Add (edge);
					  for (exp5.Init(exp4.Current(), TopAbs_VERTEX);
					       exp5.More(); exp5.Next())
					    {
					      TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
					      if (vmap.FindIndex(vertex) < 1)
						vmap.Add (vertex);
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  
  // Free Shells
  for (exp1.Init(exp0.Current(), TopAbs_SHELL, TopAbs_SOLID);
       exp1.More(); exp1.Next())
    {
      TopoDS_Shape shell = exp1.Current().Composed (exp0.Current().Orientation());
      if (shmap.FindIndex(shell) < 1)
	{
	  shmap.Add (shell);
	  
	  for (exp2.Init(shell, TopAbs_FACE);
	       exp2.More(); exp2.Next())
	    {
	      TopoDS_Face face = TopoDS::Face(exp2.Current().Composed(shell.Orientation()));
	      if (fmap.FindIndex(face) < 1)
		{
		  fmap.Add (face);
		  
		  for (exp3.Init(exp2.Current(), TopAbs_WIRE);
		       exp3.More(); exp3.Next())
		    {
		      TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
		      if (wmap.FindIndex(wire) < 1)
			{
			  wmap.Add (wire);
			  
			  for (exp4.Init(exp3.Current(), TopAbs_EDGE);
			       exp4.More(); exp4.Next())
			    {
			      TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
			      if (emap.FindIndex(edge) < 1)
				{
				  emap.Add (edge);
				  for (exp5.Init(exp4.Current(), TopAbs_VERTEX);
				       exp5.More(); exp5.Next())
				    {
				      TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
				      if (vmap.FindIndex(vertex) < 1)
					vmap.Add (vertex);
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  
  
  // Free Faces
  
  for (exp2.Init(shape, TopAbs_FACE, TopAbs_SHELL);
       exp2.More(); exp2.Next())
    {
      TopoDS_Face face = TopoDS::Face(exp2.Current());
      if (fmap.FindIndex(face) < 1)
	{
	  fmap.Add (face);
	  
	  for (exp3.Init(exp2.Current(), TopAbs_WIRE);
	       exp3.More(); exp3.Next())
	    {
	      TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
	      if (wmap.FindIndex(wire) < 1)
		{
		  wmap.Add (wire);
		  
		  for (exp4.Init(exp3.Current(), TopAbs_EDGE);
		       exp4.More(); exp4.Next())
		    {
		      TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
		      if (emap.FindIndex(edge) < 1)
			{
			  emap.Add (edge);
			  for (exp5.Init(exp4.Current(), TopAbs_VERTEX);
			       exp5.More(); exp5.Next())
			    {
			      TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
			      if (vmap.FindIndex(vertex) < 1)
				vmap.Add (vertex);
			    }
			}
		    }
		}
	    }
	}
    }


  // Free Wires
  
  for (exp3.Init(shape, TopAbs_WIRE, TopAbs_FACE);
       exp3.More(); exp3.Next())
    {
      TopoDS_Wire wire = TopoDS::Wire (exp3.Current());
      if (wmap.FindIndex(wire) < 1)
	{
	  wmap.Add (wire);
	  
	  for (exp4.Init(exp3.Current(), TopAbs_EDGE);
	       exp4.More(); exp4.Next())
	    {
	      TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
	      if (emap.FindIndex(edge) < 1)
		{
		  emap.Add (edge);
		  for (exp5.Init(exp4.Current(), TopAbs_VERTEX);
		       exp5.More(); exp5.Next())
		    {
		      TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
		      if (vmap.FindIndex(vertex) < 1)
			vmap.Add (vertex);
		    }
		}
	    }
	}
    }


  // Free Edges
  
  for (exp4.Init(shape, TopAbs_EDGE, TopAbs_WIRE);
       exp4.More(); exp4.Next())
    {
      TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
      if (emap.FindIndex(edge) < 1)
	{
	  emap.Add (edge);
	  for (exp5.Init(exp4.Current(), TopAbs_VERTEX);
	       exp5.More(); exp5.Next())
	    {
	      TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
	      if (vmap.FindIndex(vertex) < 1)
		vmap.Add (vertex);
	    }
	}
    }


  // Free Vertices
  
  for (exp5.Init(shape, TopAbs_VERTEX, TopAbs_EDGE);
       exp5.More(); exp5.Next())
    {
      TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
      if (vmap.FindIndex(vertex) < 1)
	vmap.Add (vertex);
    }



  
  facemeshstatus.SetSize (fmap.Extent());
  facemeshstatus = 0;

  fvispar.SetSize (fmap.Extent());
  evispar.SetSize (emap.Extent());
  vvispar.SetSize (vmap.Extent());
}



void OCCGeometry :: SewFaces ()
{
  (*testout) << "Trying to sew faces ..." << endl;
  cout << "Trying to sew faces ..." << flush;

  BRepOffsetAPI_Sewing sewedObj(1);
  //  BRepOffsetAPI_Sewing sewedObj(healingtolerance);

  for (int i = 1; i <= fmap.Extent(); i++)
    {
      TopoDS_Face face = TopoDS::Face (fmap(i));
      sewedObj.Add (face);
    }
  
  sewedObj.Perform();
  
  if (!sewedObj.SewedShape().IsNull())
    {
      shape = sewedObj.SewedShape();
      cout << " done" << endl;
    }
  else
    cout << " not possible";
  
  /*
  ShapeUpgrade_ShellSewing sewing;
  TopoDS_Shape sh = sewing.ApplySewing (shape);
  shape = sh;
  */
}





void OCCGeometry :: MakeSolid ()
{
  TopExp_Explorer exp0;

  (*testout) << "Trying to build solids ..." << endl;
  cout << "Trying to build solids ..." << flush;

  BRepBuilderAPI_MakeSolid ms;
  int count = 0;
  for (exp0.Init(shape, TopAbs_SHELL); exp0.More(); exp0.Next())
    {
      count++;
      ms.Add (TopoDS::Shell(exp0.Current()));
    }

  if (!count)
    {
     cout << " not possible (no shells)" << endl;
     return;
    }

  BRepCheck_Analyzer ba(ms);
  if (ba.IsValid ())
    {
      Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
      sfs->Init (ms);
  
      sfs->SetPrecision(1e-5);
      sfs->SetMaxTolerance(1e-5);
      
      sfs->Perform();

      shape = sfs->Shape();
      
      for (exp0.Init(shape, TopAbs_SOLID); exp0.More(); exp0.Next())
	{
	  TopoDS_Solid solid = TopoDS::Solid(exp0.Current());
	  TopoDS_Solid newsolid = solid;
	  BRepLib::OrientClosedSolid (newsolid);
	  Handle_ShapeBuild_ReShape rebuild = new ShapeBuild_ReShape;
	  //	  rebuild->Apply(shape);
	  rebuild->Replace(solid, newsolid, Standard_False);
	  //	  TopoDS_Shape newshape = rebuild->Apply(shape);
	  
	  TopoDS_Shape newshape = rebuild->Apply(shape, TopAbs_SHAPE, 1);
	  shape = newshape;
	}
      
      cout << " done" << endl;
    }
  else
    cout << " not possible" << endl;
}


void OCCGeometry :: BuildVisualizationMesh ()
{
  cout << "Preparing visualization (deflection = " << vispar.occdeflection << ") ... " << flush;
  BRepTools::Clean (shape);
  BRepMesh_IncrementalMesh::BRepMesh_IncrementalMesh (shape, vispar.occdeflection, true);
  cout << "done" << endl;
  
  Bnd_Box bb;
  BRepBndLib::Add (shape, bb);
  
  double x1,y1,z1,x2,y2,z2;
  bb.Get (x1,y1,z1,x2,y2,z2);
  Point<3> p1 = Point<3> (x1,y1,z1);
  Point<3> p2 = Point<3> (x2,y2,z2);
  
  (*testout) << "Bounding Box = [" << p1 << " - " << p2 << "]" << endl;
  boundingbox = Box<3> (p1,p2);
  SetCenter();
}



  bool OCCGeometry :: FastProject (int surfi, Point<3> & ap, double& u, double& v) const
  {
    gp_Pnt p(ap(0), ap(1), ap(2));
  
    Handle(Geom_Surface) surface = BRep_Tool::Surface(TopoDS::Face(fmap(surfi)));
  
    gp_Pnt x = surface->Value (u,v);
  
    if (p.SquareDistance(x) <= sqr(PROJECTION_TOLERANCE)) return true;
  
    gp_Vec du, dv;
  
    surface->D1(u,v,x,du,dv);
  
    int count = 0;
  
    gp_Pnt xold;
    gp_Vec n;
    double det, lambda, mu;
  
    do {
      count++;
  
      n = du^dv;
  
      det = Det3 (n.X(), du.X(), dv.X(),
		  n.Y(), du.Y(), dv.Y(),
		  n.Z(), du.Z(), dv.Z());
  
      if (det < 1e-15) return false; 
  
      lambda = Det3 (n.X(), p.X()-x.X(), dv.X(),
		     n.Y(), p.Y()-x.Y(), dv.Y(),
		     n.Z(), p.Z()-x.Z(), dv.Z())/det;
  
      mu     = Det3 (n.X(), du.X(), p.X()-x.X(),
		     n.Y(), du.Y(), p.Y()-x.Y(),
		     n.Z(), du.Z(), p.Z()-x.Z())/det;
    
      u += lambda;
      v += mu;
  
      xold = x;
      surface->D1(u,v,x,du,dv);
  
    } while (xold.SquareDistance(x) > sqr(PROJECTION_TOLERANCE) && count < 50);

//    (*testout) << "FastProject count: " << count << endl;
  
    if (count == 50) return false;
  
    ap = Point<3> (x.X(), x.Y(), x.Z());
  
    return true;
  }


OCCGeometry * LoadOCC_IGES (const char * filename)
{
  OCCGeometry * occgeo;
  occgeo = new OCCGeometry;
  
  IGESControl_Reader reader;
  
#ifdef OCC52
  Standard_Integer stat = reader.ReadFile((char*)filename);
#else
  Standard_Integer stat = reader.LoadFile((char*)filename);
  reader.Clear();
#endif
  
#ifdef OCC52
  reader.TransferRoots(); // Tranlate IGES -> OCC
#else
  reader.TransferRoots(Standard_False); // Tranlate IGES -> OCC
#endif

  occgeo->shape = reader.OneShape();
  occgeo->changed = 1;
  occgeo->BuildFMap();
  occgeo->BuildVisualizationMesh();
  PrintContents (occgeo);

  return occgeo;
}

OCCGeometry * LoadOCC_STEP (const char * filename)
{
  OCCGeometry * occgeo;
  occgeo = new OCCGeometry;
  
  STEPControl_Reader reader;
  Standard_Integer stat = reader.ReadFile((char*)filename);
  Standard_Integer nb = reader.NbRootsForTransfer();
  reader.TransferRoots (); // Tranlate STEP -> OCC
  occgeo->shape = reader.OneShape();
  occgeo->changed = 1;
  occgeo->BuildFMap();
  occgeo->BuildVisualizationMesh();
  PrintContents (occgeo);

  return occgeo;
}

char * shapesname[] =
  {" ", "CompSolids", "Solids", "Shells",
   "Faces", "Wires", "Edges", "Vertices"};

char * shapename[] =
  {" ", "CompSolid", "Solid", "Shell",
   "Face", "Wire", "Edge", "Vertex"};

char * orientationstring[] =
  {"+", "-"};

void OCCGeometry :: RecursiveTopologyTree (const TopoDS_Shape & sh,
					   stringstream & str,
					   TopAbs_ShapeEnum l,
					   bool isfree,
					   const char * lname)
{
  if (l > TopAbs_VERTEX) return;

  TopExp_Explorer e;
  int count = 0;
  int count2;

  if (isfree)
    e.Init(sh, l, TopAbs_ShapeEnum(l-1));
  else
    e.Init(sh, l);

  for (; e.More(); e.Next())
    {
      count++;

      stringstream lname2;
      lname2 << lname << "/" << shapename[l] << count;
      str << lname2.str() << " ";

      switch (e.Current().ShapeType())
	{
	case TopAbs_SOLID:
	  count2 = somap.FindIndex(TopoDS::Solid(e.Current())); break;
	case TopAbs_SHELL:
	  count2 = shmap.FindIndex(TopoDS::Shell(e.Current())); break;
	case TopAbs_FACE:
	  count2 = fmap.FindIndex(TopoDS::Face(e.Current())); break;
	case TopAbs_WIRE:
	  count2 = wmap.FindIndex(TopoDS::Wire(e.Current())); break;
	case TopAbs_EDGE:
	  count2 = emap.FindIndex(TopoDS::Edge(e.Current())); break;
	case TopAbs_VERTEX:
	  count2 = vmap.FindIndex(TopoDS::Vertex(e.Current())); break;
	}
      
      int nrsubshapes = 0;
      
      if (l <= TopAbs_WIRE)
	{
	  TopExp_Explorer e2;
	  for (e2.Init (e.Current(), TopAbs_ShapeEnum (l+1));
	       e2.More(); e2.Next())
	    nrsubshapes++;
	}
      
      str << "{" << shapename[l] << " " << count2;
      
      if (l <= TopAbs_EDGE)
	{
	  str << " (" << orientationstring[e.Current().Orientation()];
	  if (nrsubshapes != 0) str << ", " << nrsubshapes;
	  str << ") } ";
	}
      else
	str << " } ";

      RecursiveTopologyTree (e.Current(), str, TopAbs_ShapeEnum (l+1),
			     false, (char*)lname2.str().c_str());
      
    }
}

void OCCGeometry :: GetTopologyTree (stringstream & str)
{
  cout << "Building topology tree ... " << flush;
  RecursiveTopologyTree (shape, str, TopAbs_COMPSOLID, false, "CompSolids");
  RecursiveTopologyTree (shape, str, TopAbs_SOLID, true, "FreeSolids");
  RecursiveTopologyTree (shape, str, TopAbs_SHELL, true, "FreeShells");
  RecursiveTopologyTree (shape, str, TopAbs_FACE, true, "FreeFaces");
  RecursiveTopologyTree (shape, str, TopAbs_WIRE, true, "FreeWires");
  RecursiveTopologyTree (shape, str, TopAbs_EDGE, true, "FreeEdges");
  RecursiveTopologyTree (shape, str, TopAbs_VERTEX, true, "FreeVertices");
  str << flush;
  //  cout << "done" << endl;
}

void OCCGeometry :: CheckIrregularEntities(stringstream & str)
{
  ShapeAnalysis_CheckSmallFace csm;

  csm.SetTolerance (1e-6);

  TopTools_DataMapOfShapeListOfShape mapEdges;
  ShapeAnalysis_DataMapOfShapeListOfReal mapParam;
  TopoDS_Compound theAllVert;

  int spotfaces = 0;
  int stripsupportfaces = 0;
  int singlestripfaces = 0;
  int stripfaces = 0;
  int facessplitbyvertices = 0;
  int stretchedpinfaces = 0;
  int smoothpinfaces = 0;
  int twistedfaces = 0;
  int edgessamebutnotidentified = 0;

  cout << "checking faces ... " << flush;

  int i;
  for (i = 1; i <= fmap.Extent(); i++)
    {
      TopoDS_Face face = TopoDS::Face (fmap(i));
      TopoDS_Edge e1, e2;

      if (csm.CheckSpotFace (face))
	{
	  if (!spotfaces++)
	    str << "SpotFace {Spot face} ";

	  (*testout) << "Face " << i << " is a spot face" << endl;
	  str << "SpotFace/Face" << i << " ";
	  str << "{Face " << i << " } ";
	}

      if (csm.IsStripSupport (face))
	{
	  if (!stripsupportfaces++)
	    str << "StripSupportFace {Strip support face} ";

	  (*testout) << "Face " << i << " has strip support" << endl;
	  str << "StripSupportFace/Face" << i << " ";
	  str << "{Face " << i << " } ";
	}

      if (csm.CheckSingleStrip(face, e1, e2))
	{
	  if (!singlestripfaces++)
	    str << "SingleStripFace {Single strip face} ";

	  (*testout) << "Face " << i << " is a single strip (edge " << emap.FindIndex(e1)
	       << " and edge " << emap.FindIndex(e2) << " are identical)" << endl;
	  str << "SingleStripFace/Face" << i << " ";
	  str << "{Face " << i << " (edge " << emap.FindIndex(e1)
	       << " and edge " << emap.FindIndex(e2) << " are identical)} ";
	}

      if (csm.CheckStripFace(face, e1, e2))
	{
	  if (!stripfaces++)
	    str << "StripFace {Strip face} ";

	  (*testout) << "Face " << i << " is a strip (edge " << emap.FindIndex(e1)
		     << " and edge " << emap.FindIndex(e2)
		     << " are identical)" << endl;
	  str << "StripFace/Face" << i << " ";
	  str << "{Face " << i << " (edge " << emap.FindIndex(e1)
	      << " and edge " << emap.FindIndex(e2) << " are identical)} ";
	}

      if (int count = csm.CheckSplittingVertices(face, mapEdges, mapParam, theAllVert))
	{
	  if (!facessplitbyvertices++)
	    str << "FaceSplitByVertices {Face split by vertices} ";

	  (*testout) << "Face " << i << " is split by " << count
		     << " vertex/vertices " << endl;
	  str << "FaceSplitByVertices/Face" << i << " ";
	  str << "{Face " << i << " (split by " << count << "vertex/vertices)} ";
	}
      
      int whatrow, sens;
      if (int type = csm.CheckPin (face, whatrow, sens))
	{
	  if (type == 1)
	    {
	      if (!smoothpinfaces++)
		str << "SmoothPinFace {Smooth pin face} ";

	      (*testout) << "Face " << i << " is a smooth pin" << endl;
	      str << "SmoothPinFace/Face" << i << " ";
	      str << "{Face " << i << " } ";
	    }
	  else
	    {
	      if (!stretchedpinfaces++)
		str << "StretchedPinFace {Stretched pin face} ";

	      (*testout) << "Face " << i << " is a streched pin" << endl;
	      str << "StretchedPinFace/Face" << i << " ";
	      str << "{Face " << i << " } ";
	    }
	}

      double paramu, paramv;
      if (csm.CheckTwisted (face, paramu, paramv))
	{
	  if (!twistedfaces++)
	    str << "TwistedFace {Twisted face} ";

	  (*testout) << "Face " << i << " is twisted" << endl;
	  str << "TwistedFace/Face" << i << " ";
	  str << "{Face " << i << " } ";
	}
    }

  cout << "done" << endl;
  cout << "checking edges ... " << flush;

  double dmax;
  int cnt = 0;
  ARRAY <double> edgeLengths;
  ARRAY <int> order;
  edgeLengths.SetSize (emap.Extent());
  order.SetSize (emap.Extent());

  for (i = 1; i <= emap.Extent(); i++)
    {
      TopoDS_Edge edge1 = TopoDS::Edge (emap(i));
      GProp_GProps system;
      BRepGProp::LinearProperties(edge1, system);
      edgeLengths[i-1] = system.Mass();
      /*
      int j;
      for (j = i+1; j <= emap.Extent(); j++)
	{
	  TopoDS_Edge edge2 = TopoDS::Edge (emap(j));

	  if (csm.CheckStripEdges(edge1, edge2, csm.Tolerance(), dmax))
	    {
	      if (!edgessamebutnotidentified++)
		str << "EdgesSameButNotIdentified {Edges same but not identified} ";

	      cnt++;
	      (*testout) << "Edge " << i << " and edge " << j
			 << " are on one strip (same but not identified)" << endl;
	      str << "EdgesSameButNotIdentified/Edge" << cnt << " ";
	      str << "{Edge " << i << " and Edge " << j << "} ";
	    }
	}
      */
    }

  Sort (edgeLengths, order);

  str << "ShortestEdges {Shortest edges} ";
  for (i = 1; i <= min(20, emap.Extent()); i++)
    {
      str << "ShortestEdges/Edge" << i;
      str << " {Edge " << order[i-1] << " (L=" << edgeLengths[order[i-1]-1] << ")} ";
    }

  str << flush;

  cout << "done" << endl;

  /*
  for (i = 1; i <= shmap.Extent(); i++)
    {
      TopoDS_Shell shell = TopoDS::Shell (shmap(i));
      if (!shell.Closed()) 
	cout << "Shell " << i << " is not closed" << endl;
      if (shell.Infinite()) 
	cout << "Shell " << i << " is infinite" << endl;

      BRepCheck_Analyzer ba(shell);
      if (!ba.IsValid ())
	cout << "Shell " << i << " is not valid" << endl;
    }

  for (i = 1; i <= somap.Extent(); i++)
    {
      TopoDS_Solid solid = TopoDS::Solid (somap(i));
      if (!solid.Closed()) 
	cout << "Solid " << i << " is not closed" << endl;
      if (solid.Infinite()) 
	cout << "Solid " << i << " is infinite" << endl;

      BRepCheck_Analyzer ba(solid);
      if (!ba.IsValid ())
	cout << "Solid " << i << " is not valid" << endl;
    }
  */


}


void OCCGeometry :: GetUnmeshedFaceInfo (stringstream & str)
{
  for (int i = 1; i <= fmap.Extent(); i++)
    {
      if (facemeshstatus[i-1] == -1)
	str << "Face" << i << " {Face " << i << " } ";
    }
  str << flush;
}

void OCCGeometry :: GetNotDrawableFaces (stringstream & str)
{
  for (int i = 1; i <= fmap.Extent(); i++)
    {
      if (!fvispar[i-1].IsDrawable())
	str << "Face" << i << " {Face " << i << " } ";
    }
  str << flush;
}

bool OCCGeometry :: ErrorInSurfaceMeshing ()
{
  for (int i = 1; i <= fmap.Extent(); i++)
    if (facemeshstatus[i-1] == -1)
      return true;

  return false;
}

}


#endif
