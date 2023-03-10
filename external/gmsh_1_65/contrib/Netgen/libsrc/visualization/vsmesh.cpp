#include <mystdlib.h>
#include "incvis.hpp"


#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

namespace netgen
{

#include "mvdraw.hpp"


  // #define FAST3DELEMENTS



  extern AutoPtr<Mesh> mesh;
  extern STLGeometry * stlgeometry;
  VisualSceneMesh vsmesh;



  VisualSceneMesh :: VisualSceneMesh ()
    : VisualScene()
  {
    filledlist = 0;
    linelist = 0;
    badellist = 0;
    tetlist = 0;
    prismlist = 0;
    hexlist = 0;
    pyramidlist = 0;
    identifiedlist = 0;
    pointnumberlist = 0;
    domainsurflist = 0;

    vstimestamp = GetTimeStamp();
    selecttimestamp = GetTimeStamp();
    filledtimestamp = GetTimeStamp();
    linetimestamp = GetTimeStamp();
    pointnumbertimestamp = GetTimeStamp();
  
    tettimestamp = GetTimeStamp();
    prismtimestamp = GetTimeStamp();
    hextimestamp = GetTimeStamp();
    pyramidtimestamp = GetTimeStamp();
  
    badeltimestamp = GetTimeStamp();
    identifiedtimestamp = GetTimeStamp();
    domainsurftimestamp = GetTimeStamp();


    selface = -1;
    selelement = -1;
    locpi = 1;
    selpoint = -1;
    selpoint2 = -1;
    seledge = -1;
  }

  VisualSceneMesh :: ~VisualSceneMesh ()
  {
    ;
  }

  
  // ARRAY<Point3d> drawel;

  void VisualSceneMesh :: DrawScene ()
  {
    if (!mesh) 
      {
	VisualScene::DrawScene();
	return;
      }
    
    lock = NULL;
    
    BuildScene();
    
    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glEnable (GL_COLOR_MATERIAL);
    glColor3f (1.0f, 1.0f, 1.0f);
    glLineWidth (1.0f);

    SetLight();

    glPushMatrix();
    glMultMatrixf (transformationmat);

    GLdouble projmat[16];
    glGetDoublev (GL_PROJECTION_MATRIX, projmat);
  
  
    glInitNames ();
    glPushName (0);

    //    glEnable (GL_LINE_SMOOTH);
    //    glEnable (GL_BLEND);
    //    glEnable (GL_POLYGON_SMOOTH);
    //    glDisable (GL_DEPTH_TEST);
    //    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //    glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  
    glDisable (GL_COLOR_MATERIAL);
  
    GLfloat matcol0[] = { 0, 0, 0, 1 };
    GLfloat matcol1[] = { 1, 1, 1, 1 };
    GLfloat matcolf[] = { 0, 1, 0, 1 };
    GLfloat matcolb[] = { 0.5, 0, 0, 1 };
    GLfloat matcolblue[] = { 0, 0, 1, 1 };
  
    glMatrixMode (GL_MODELVIEW); 
  
    glMaterialfv(GL_FRONT, GL_EMISSION, matcol0);
    glMaterialfv(GL_BACK, GL_EMISSION, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matcol1);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcolf);
    glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, matcolb);
  
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    // glPolygonOffset (1,10);
    glPolygonOffset (2,2);
    glEnable (GL_POLYGON_OFFSET_FILL);

    SetClippingPlane ();

    if (vispar.drawfilledtrigs)
      {
	if (filledtimestamp < mesh->GetTimeStamp () ||
	    filledtimestamp < selecttimestamp)
	  {
	    BuildFilledList ();
	  }
	glCallList (filledlist);
      }

    if (vispar.drawbadels)
      glCallList (badellist);

    if (vispar.drawprisms)
      {
	BuildPrismList ();
	static float prismcol[] = { 1.0f, 1.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, prismcol);
	glLineWidth (1.0f);
	glCallList (prismlist);
      }

    if (vispar.drawpyramids)
      {
	BuildPyramidList ();
	static float pyramidcol[] = { 1.0f, 1.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, pyramidcol);
	glLineWidth (1.0f);
	glCallList (pyramidlist);
      }

    if (vispar.drawhexes)
      {
	BuildHexList ();
	static float hexcol[] = { 1.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, hexcol);
	glLineWidth (1.0f);
	glCallList (hexlist);
      }

    if (vispar.drawtets)
      {
	BuildTetList ();
	static float tetcol[] = { 1.0f, 1.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, tetcol);
	glLineWidth (1.0f);
	glCallList (tetlist);
      }

    if (vispar.drawdomainsurf)
      {
	BuildDomainSurfList();
	glCallList (domainsurflist);
      }

    glDisable (GL_POLYGON_OFFSET_FILL);
  
    // draw lines

    glMatrixMode (GL_MODELVIEW); 
  
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matcol0);
  
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth (1.0f);
    glColor3f (0.0f, 0.0f, 0.0f);
    glDisable (GL_LINE_SMOOTH);
  
    if (vispar.drawoutline)
      {
	glPolygonOffset (1, 1);
	glEnable (GL_POLYGON_OFFSET_LINE);

	if (linetimestamp < mesh->GetTimeStamp ())
	  BuildLineList ();

	glCallList (linelist);
	glDisable (GL_POLYGON_OFFSET_LINE);
      }

    if (vispar.drawidentified)
      {
	glPolygonOffset (1, -1);
	glEnable (GL_POLYGON_OFFSET_LINE);
	glCallList (identifiedlist);
	glDisable (GL_POLYGON_OFFSET_LINE);
      }
  
    if (vispar.drawpointnumbers ||
	vispar.drawedgenumbers ||
	vispar.drawfacenumbers ||
	vispar.drawelementnumbers)
      glCallList (pointnumberlist);
  
  
    glPopName();

    if (vispar.drawedges)
      {
	GLfloat matcoledge[] = { 0, 0, 1, 1 };
	GLfloat matcolsingedge[] = { 1, 0, 1, 1 };

	glEnable (GL_POLYGON_OFFSET_LINE);
	glPolygonOffset (1, -1);
	glLineWidth (2);
      

	for (int i = 1; i <= mesh->GetNSeg(); i++)
	  {
	    const Segment & seg = mesh->LineSegment(i);
	    const Point3d & p1 = (*mesh)[seg.p1];
	    const Point3d & p2 = (*mesh)[seg.p2];

	    if (seg.singedge_left || seg.singedge_right)
	      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, 
			    matcolsingedge);
	    else
	      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, 
			    matcoledge);

	    if (seg.edgenr == seledge)
	      glLineWidth(5);
	    else
	      glLineWidth(2);

	    if (mesh->GetCurvedElements().IsHighOrder()) {

	      int j;
	      int hoplotn = 1 << vispar.subdivisions; 
	      // mesh->GetCurvedElements().GetNVisualSubsecs();

	      Point<3> x;
	      glBegin (GL_LINE_STRIP);

	      for (int j = 0; j <= hoplotn; j++) 
		{
		  mesh->GetCurvedElements().CalcSegmentTransformation ((double) j/hoplotn, i-1, x);
		  glVertex3d (x(0), x(1), x(2));
		}
	      
	      glEnd();
	      
	    } else {
	      
	      glBegin (GL_LINES);
	      glVertex3f (p1.X(), p1.Y(), p1.Z());
	      glVertex3f (p2.X(), p2.Y(), p2.Z());
	      glEnd();

	    }
	  }

	glLineWidth (2);
	glDisable (GL_POLYGON_OFFSET_LINE);
      }


    if (selpoint > 0 && selpoint <= mesh->GetNP())
      {
	glPointSize (10);
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcolblue);
	glBegin (GL_POINTS);
      
	const Point3d p = mesh->Point(selpoint);
	glVertex3f (p.X(), p.Y(), p.Z());
	glEnd();
      }


    glDisable(GL_CLIP_PLANE0);

    glPopMatrix();

    if (vispar.colormeshsize)
      DrawColorBar (minh, maxh, 1);

    DrawCoordinateCross ();
    DrawNetgenLogo ();

    if (lock)
      {
	lock -> UnLock();
	delete lock;
      }

    glFinish();  
  }


  void VisualSceneMesh :: BuildScene (int zoomall)
  {
    if (!mesh)
      {
	VisualScene::BuildScene (zoomall);
	return;
      }
      
    int i, j;
	
	
    Point3d pmin, pmax;
    static double oldrad = 0;
	
    ARRAY<Element2d> faces;

    int meshtimestamp = mesh->GetTimeStamp();
    if (meshtimestamp > vstimestamp || zoomall)
      {
	mesh->GetBox (pmin, pmax, SURFACEPOINT);
	

	if (selpoint >= 1 && zoomall == 2)
	  center = mesh->Point (selpoint);
	else if (vispar.use_center_coords && zoomall == 2)
	  {
	    center.X() = vispar.centerx; center.Y() = vispar.centery; center.Z() = vispar.centerz;
	  }
	else if (vispar.centerpoint >= 1 && zoomall == 2)
	  center = mesh->Point (vispar.centerpoint);
	else
	  center = Center (pmin, pmax);
	        
	rad = 0.5 * Dist (pmin, pmax);
      
      
	if (rad > 1.5 * oldrad || 
	    mesh->GetMajorTimeStamp() > vstimestamp || 
	    zoomall)
	  {
	    CalcTransformationMatrices();
	    oldrad = rad;
	  }
      }

    glEnable (GL_NORMALIZE);

    if (pointnumberlist)
      {
	glDeleteLists (pointnumberlist, 1);      
	pointnumberlist = 0;
      }

    if (badellist)
      {
	glDeleteLists (badellist, 1);
	badellist = 0;
      }
    /*
      if (prismlist)
      {
      glDeleteLists (prismlist, 1);
      prismlist = 0;
      }

      if (pyramidlist)
      {
      glDeleteLists (pyramidlist, 1);
      pyramidlist = 0;
      }

      if (hexlist)
      {
      glDeleteLists (hexlist, 1);
      hexlist = 0;
      }
    */
    if (identifiedlist)
      {
	glDeleteLists (identifiedlist, 1);
	identifiedlist = 0;
      }


    pointnumberlist = glGenLists (1);
    glNewList (pointnumberlist, GL_COMPILE);

    if (vispar.drawpointnumbers ||
	vispar.drawedgenumbers ||
	vispar.drawfacenumbers ||
	vispar.drawelementnumbers)
      {
	glEnable (GL_COLOR_MATERIAL);
	GLfloat textcol[3] = { 1 - backcolor,
			       1 - backcolor,
			       1 - backcolor };
	glColor3fv (textcol);
	glNormal3d (0, 0, 1);
	glPushAttrib (GL_LIST_BIT);
	glListBase (fontbase);
      
	char buf[30];

	if (vispar.drawpointnumbers)
	  for (i = 1; i <= mesh->GetNP(); i++)
	    {
	      const Point3d & p = mesh->Point(i);
	      glRasterPos3d (p.X(), p.Y(), p.Z());

	      sprintf (buf, "%d", i);

	      glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
	    }

	if (vispar.drawedgenumbers)
	  {
	    const MeshTopology & top = mesh->GetTopology();
	    for (i = 1; i <= top.GetNEdges(); i++)
	      {
		int v1, v2;
		top.GetEdgeVertices (i, v1, v2);
		const Point3d & p1 = mesh->Point(v1);
		const Point3d & p2 = mesh->Point(v2);
		const Point3d p = Center (p1, p2);
		glRasterPos3d (p.X(), p.Y(), p.Z());

		sprintf (buf, "%d", i);

		glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
	      }
	  }      


	if (vispar.drawfacenumbers)
	  {
	    const MeshTopology & top = mesh->GetTopology();
	    ARRAY<int> v;
	    for (i = 1; i <= top.GetNFaces(); i++)
	      {
		top.GetFaceVertices (i, v);
		const Point3d & p1 = mesh->Point(v.Elem(1));
		const Point3d & p2 = mesh->Point(v.Elem(2));
		const Point3d & p3 = mesh->Point(v.Elem(3));
		Point3d p;
		if (v.Elem(4) == 0)
		  {
		    p = Center (p1, p2, p3);
		  }
		else
		  {
		    const Point3d & p4 = mesh->Point(v.Elem(4));
		    Point3d hp1 = Center (p1, p2);
		    Point3d hp2 = Center (p3, p4);
		    p = Center (hp1, hp2);
		  }

		glRasterPos3d (p.X(), p.Y(), p.Z());
		sprintf (buf, "%d", i);
		glCallLists (strlen (buf), GL_UNSIGNED_BYTE, buf);
	      }
	  }      


	glPopAttrib ();
	glDisable (GL_COLOR_MATERIAL);
      }
    glEndList ();













    badellist = glGenLists (1);
    glNewList (badellist, GL_COMPILE);

    if (vispar.drawbadels)
      {
	//  SetClippingPlane ();

	static float badelcol[] = { 1.0f, 0.0f, 1.0f, 1.0f };
	glLineWidth (1.0f);

	for (i = 1; i <= mesh->GetNE(); i++)
	  {
	    if (mesh->VolumeElement(i).flags.badel || 
		mesh->VolumeElement(i).flags.illegal || 
		(i == vispar.drawelement))
	      {
		// copy to be thread-safe
		Element el = mesh->VolumeElement (i);
		el.GetSurfaceTriangles (faces);

		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, badelcol);


		//	  if ( (el.GetNP() == 4) || (el.GetNP() == 10))
		if (el.PNum(1))
		  {
		    glBegin (GL_TRIANGLES);
	      
		    for (j = 1; j <= faces.Size(); j++)
		      {
			Element2d & face = faces.Elem(j);
			const Point3d & lp1 = mesh->Point (el.PNum(face.PNum(1)));
			const Point3d & lp2 = mesh->Point (el.PNum(face.PNum(2)));
			const Point3d & lp3 = mesh->Point (el.PNum(face.PNum(3)));
			Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
			n /= (n.Length()+1e-12);
			glNormal3d (n.X(), n.Y(), n.Z());
			glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
			glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
			glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
		      }
	      
		    glEnd();
		  }
	      }
	  }



	for (i = 1; i <= mesh->GetNE(); i++)
	  {
	    if (mesh->VolumeElement(i).flags.badel)
	      {
		// copy to be thread-safe
		Element el = mesh->VolumeElement (i);
		if ( (el.GetNP() == 4) || (el.GetNP() == 10))
		  {
		    glBegin (GL_LINES);
		    glVertex3d (0,0,0);
		    const Point3d & p = mesh->Point(el.PNum(1));
		    glVertex3d (p.X(), p.Y(), p.Z());
		    glEnd();
		  }
	      }
	  }
  

	for (i = 1; i <= mesh->GetNE(); i++)
	  {
	    Element el = mesh->VolumeElement (i);
	    int hascp = 0;
	    for (j = 1; j <= el.GetNP(); j++)
	      if (el.PNum(j) == vispar.centerpoint)
		hascp = 1;

	    if (hascp)
	      {
		(*testout) << "draw el " << i << " : ";
		for (j = 1; j <= el.GetNP(); j++)
		  (*testout) << el.PNum(j) << " ";
		(*testout) << endl;

		if (el.GetNP() == 4)
		  {
		    int et[6][2] = 
		      { { 1, 2 },
			{ 1, 3 },
			{ 1, 4 },
			{ 2, 3 },
			{ 2, 4 },
			{ 3, 4 } } ;

		    for (j = 0; j < 6; j++)
		      {
			glBegin (GL_LINES);
			const Point3d & p1 = mesh->Point (el.PNum(et[j][0]));
			const Point3d & p2 = mesh->Point (el.PNum(et[j][1]));
			glVertex3d (p1.X(), p1.Y(), p1.Z());
			glVertex3d (p2.X(), p2.Y(), p2.Z());
			glEnd ();
		      }
		  }


		if (el.GetNP() == 10)
		  {
		    int et[12][2] = 
		      { { 1, 5 },
			{ 2, 5 },
			{ 1, 6 },
			{ 3, 6 },
			{ 1, 7 },
			{ 4, 7 },
			{ 2, 8 },
			{ 3, 8 },
			{ 2, 9 },
			{ 4, 9 },
			{ 3, 10 },
			{ 4, 10 } };

		    for (j = 0; j < 12; j++)
		      {
			glBegin (GL_LINES);
			const Point3d & p1 = mesh->Point (el.PNum(et[j][0]));
			const Point3d & p2 = mesh->Point (el.PNum(et[j][1]));
			glVertex3d (p1.X(), p1.Y(), p1.Z());
			glVertex3d (p2.X(), p2.Y(), p2.Z());
			glEnd ();
		      }
		  }
	      }
	  }


	for (i = 1; i <= mesh->GetNSE(); i++)
	  {
	    Element2d el = mesh->SurfaceElement(i);
	    if (!el.BadElement())
	      continue;

	    int drawel = 1;
	    for (j = 1; j <= el.GetNP(); j++)
	      if (!el.PNum(j))
		drawel = 0;

	    if (!drawel)
	      continue;

	    cout << int (el.GetType()) << " " << flush;
	    switch (el.GetType())
	      {
	      case TRIG:
		{
		  glBegin (GL_TRIANGLES);
	    
		  Point3d & lp1 = mesh->Point (el.PNum(1));
		  Point3d & lp2 = mesh->Point (el.PNum(2));
		  Point3d & lp3 = mesh->Point (el.PNum(3));
		  Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
		  n /= (n.Length() + 1e-12);
		  glNormal3dv (&n.X());
		  glVertex3dv (&lp1.X());
		  glVertex3dv (&lp2.X());
		  glVertex3dv (&lp3.X());
		  glEnd();
		  break;
		}
	      case QUAD:
		{
		  glBegin (GL_QUADS);
	    
		  const Point3d & lp1 = mesh->Point (el.PNum(1));
		  const Point3d & lp2 = mesh->Point (el.PNum(2));
		  const Point3d & lp3 = mesh->Point (el.PNum(4));
		  const Point3d & lp4 = mesh->Point (el.PNum(3));
		  Vec3d n = Cross (Vec3d (lp1, lp2), 
				   Vec3d (lp1, Center (lp3, lp4)));
		  n /= (n.Length() + 1e-12);
		  glNormal3d (n.X(), n.Y(), n.Z());
		  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		  glVertex3d (lp4.X(), lp4.Y(), lp4.Z());
		  glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
		  glEnd();
		  break;
		}
	      case TRIG6:
		{
		  int lines[6][2] = {
		    { 1, 6 }, { 2, 6 },
		    { 1, 5 }, { 3, 5 },
		    { 2, 4 }, { 3, 4 } };
	      
		  glBegin (GL_LINES);
		  for (j = 0; j < 6; j++)
		    {
		      glVertex3dv (&mesh->Point (el.PNum(lines[j][0])).X());
		      glVertex3dv (&mesh->Point (el.PNum(lines[j][0])).X());
		    }
		  glEnd();
		  break;
		}

	      case QUAD6:
		{
		  int lines[6][2] = {
		    { 1, 5 }, { 2, 5 },
		    { 3, 6 }, { 4, 6 },
		    { 1, 4 }, { 2, 3 } };
	      
		  glBegin (GL_LINES);
	    
		  for (j = 0; j < 6; j++)
		    {
		      const Point3d & lp1 = mesh->Point (el.PNum(lines[j][0]));
		      const Point3d & lp2 = mesh->Point (el.PNum(lines[j][1]));
		
		      glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		      glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		    }
		  glEnd ();
		  break;
		}
	      default:
		PrintSysError ("Cannot draw surface element of type ", 
			       int(el.GetType()));
	      }
	  }
	glLoadName (0);  

      }
    glEndList ();
  


  

    if (1)
      {
      
	identifiedlist = glGenLists (1);
	glNewList (identifiedlist, GL_COMPILE);
  
	GLfloat identifiedcol[] = { 1, 0, 1, 1 };
  
	glLineWidth (3);
      

	//  for (i = 1; i <= mesh->GetNSeg(); i++)
	INDEX_2_HASHTABLE<int> & idpts = 
	  mesh->GetIdentifications().GetIdentifiedPoints();
	if (&idpts)
	  for (i = 1; i <= idpts.GetNBags(); i++)
	    for (j = 1; j <= idpts.GetBagSize(i); j++)
	      {
		INDEX_2 pts;
		int val;

		idpts.GetData (i, j, pts, val);
		const Point3d & p1 = mesh->Point(pts.I1());
		const Point3d & p2 = mesh->Point(pts.I2());
      
		glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, 
			      identifiedcol);

		glBegin (GL_LINES);
		glVertex3f (p1.X(), p1.Y(), p1.Z());
		glVertex3f (p2.X(), p2.Y(), p2.Z());
		glEnd(); 
	      }

	glEndList ();
      }

    vstimestamp = meshtimestamp;
  }




  void VisualSceneMesh :: BuildFilledList()
  {
    // clock_t starttime, endtime;
    // starttime = clock();

    if (!lock)
      {
	lock = new NgLock (mesh->Mutex());
	lock -> Lock();
      }

    filledtimestamp = NextTimeStamp();

    if (filledlist)
      glDeleteLists (filledlist, 1);

    filledlist = glGenLists (1);
    glNewList (filledlist, GL_COMPILE);


    bool checkvicinity = (stlgeometry != NULL) && stldoctor.showvicinity;

    glEnable (GL_NORMALIZE);
      
    glLineWidth (1.0f);
  
    Vector locms;

    if (vispar.colormeshsize)
      {
	glEnable (GL_COLOR_MATERIAL);
	locms.SetSize (mesh->GetNP());
	double maxh = -1;
	double minh = 1e99;
	for (int i = 1; i <= locms.Size(); i++)
	  {
	    Point3d p = mesh->Point(i);
	    locms.Elem(i) = mesh->GetH (p);
	    if (locms.Elem(i) > maxh) maxh = locms.Elem(i);
	    if (locms.Elem(i) < minh) minh = locms.Elem(i);
	  }
	if (!locms.Size())
	  { minh = 1; maxh = 10; }
      }
    else
      glDisable (GL_COLOR_MATERIAL);


    GLfloat matcol[] = { 0, 1, 0, 1 };
    GLfloat matcolsel[] = { 1, 0, 0, 1 };

    CurvedElements & curv = mesh->GetCurvedElements();
    int hoplotn = 1 << vispar.subdivisions; 

	
    for (int col = 1; col <= 2; col++)
      {
	if (col == 2)
	  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcolsel);
	else
	  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcol);


	for (SurfaceElementIndex sei = 0; sei < mesh->GetNSE(); sei++)
	  {
	    const Element2d & el = (*mesh)[sei];
	    
	    bool drawel = !el.IsDeleted();

	    if (checkvicinity)
	      for (int j = 0; j < el.GetNP(); j++)
		if (!stlgeometry->Vicinity(el.GeomInfoPi(j+1).trignum))
		  drawel = 0;

	    if (!drawel)
	      continue;

	    if (vispar.colormeshsize && col == 2)
	      continue;
	    if (!vispar.colormeshsize &&
		(col == 2) != (el.GetIndex() == selface))
	      continue;

	    glLoadName (sei+1);

	    switch (el.GetType())
	      {
	      case TRIG:
		{
		  if (curv.IsHighOrder() && curv.IsSurfaceElementCurved(sei))
		    {
		      Point<2> xr[3];
		      Point<3> xg;
		      Vec<3> dx, dy, n;

		      glBegin (GL_TRIANGLES);

		      for (int i = 0; i < hoplotn; i++)
			for (int j = 0; j < hoplotn-i; j++)
			  for (int k = 0; k < 2; k++)
			    {
			      if (k == 0)
				{
				  xr[0](0) = (double)    i/hoplotn; xr[0](1) = (double)    j/hoplotn;
				  xr[1](0) = (double)(i+1)/hoplotn; xr[1](1) = (double)    j/hoplotn;
				  xr[2](0) = (double)    i/hoplotn; xr[2](1) = (double)(j+1)/hoplotn;
				} 
			      else
				{
				  if (j == hoplotn-i-1) continue;
				  xr[0](0) = (double)(i+1)/hoplotn; xr[0](1) = (double)    j/hoplotn;
				  xr[1](0) = (double)(i+1)/hoplotn; xr[1](1) = (double)(j+1)/hoplotn;
				  xr[2](0) = (double)    i/hoplotn; xr[2](1) = (double)(j+1)/hoplotn;
				};
			      
			      for (int l=0; l<3; l++)
				{
				  Mat<3,2> dxdxi;

				  curv.CalcSurfaceTransformation (xr[l], sei, xg, dxdxi);
				  for (int i = 0; i < 3; i++)
				    {
				      dx(i) = dxdxi(i,0);
				      dy(i) = dxdxi(i,1);
				    }
				  n = Cross (dx, dy);
				  n.Normalize();
				  glNormal3d (n(0), n(1), n(2));
				  glVertex3d (xg(0), xg(1), xg(2));
				}
			    }
		      
		      glEnd();
		    } 
		  else // not high order
		    {
		      glBegin (GL_TRIANGLES);
		      
		      const Point<3> & lp0 = (*mesh) [el[0]];
		      const Point<3> & lp1 = (*mesh) [el[1]];
		      const Point<3> & lp2 = (*mesh) [el[2]];
		      
		      Vec<3> n = Cross (lp1-lp0, lp2-lp0);
		      glNormal3dv (n);

		      if (vispar.colormeshsize)
			{
			  SetOpenGlColor  (locms.Get(el[0]), minh, maxh, 1);
			  glVertex3dv (lp0);
			  SetOpenGlColor  (locms.Get(el[1]), minh, maxh, 1);
			  glVertex3dv (lp1);
			  SetOpenGlColor  (locms.Get(el[2]), minh, maxh, 1);
			  glVertex3dv (lp2);
			}
		      else
			{
			  glVertex3dv (lp0);
			  glVertex3dv (lp1);
			  glVertex3dv (lp2);
			}

		      glEnd();
		    }
		  
		  break;
		}
	      case QUAD:
		{
		  // cout << "BuildFilledList: QUAD" << endl;
		  // CurvedElements & curv = mesh->GetCurvedElements();
		  if (curv.IsHighOrder() && curv.IsSurfaceElementCurved(sei))
		    {
		      Point<2> xr[4];
		      Point<3> xg;
		      Vec<3> dx, dy, n;

		      glBegin (GL_QUADS);

		      for (int i = 0; i < hoplotn; i++)
			for (int j = 0; j < hoplotn; j++)
			  {
			    xr[0](0) = (double)    i/hoplotn; xr[0](1) = (double)    j/hoplotn;
			    xr[1](0) = (double)(i+1)/hoplotn; xr[1](1) = (double)    j/hoplotn;
			    xr[2](0) = (double)(i+1)/hoplotn; xr[2](1) = (double)(j+1)/hoplotn;
			    xr[3](0) = (double)    i/hoplotn; xr[3](1) = (double)(j+1)/hoplotn;

			    for (int l=0; l<4; l++)
			      {
				Mat<3,2> dxdxi;

				curv.CalcSurfaceTransformation (xr[l], sei, xg, dxdxi);
				for (int i = 0; i < 3; i++)
				  {
				    dx(i) = dxdxi(i,0);
				    dy(i) = dxdxi(i,1);
				  }

				n = Cross (dx, dy);
				n.Normalize();
				glNormal3d (n(0), n(1), n(2));
				glVertex3d (xg(0), xg(1), xg(2));
			      }

			  }
                    
		      glEnd();
		    } 

		  else // not high order

		    {
		      glBegin (GL_QUADS);
		      
		      const Point<3> & lp1 = mesh->Point (el.PNum(1));
		      const Point<3> & lp2 = mesh->Point (el.PNum(2));
		      const Point<3> & lp3 = mesh->Point (el.PNum(4));
		      const Point<3> & lp4 = mesh->Point (el.PNum(3));

		      Vec<3> n = Cross (lp2-lp1,  Center (lp3, lp4)-lp1);
		      glNormal3dv (n);

		      glVertex3dv (lp1);
		      glVertex3dv (lp2);
		      glVertex3dv (lp4);
		      glVertex3dv (lp3);
		      
		      glEnd ();
		    }
		  break;
		}

	      case TRIG6:
		{
		  glBegin (GL_TRIANGLES);

		  static int trigs[4][3] = {
		    { 1, 6, 5 },
		    { 2, 4, 6 },
		    { 3, 5, 4 },
		    { 4, 5, 6 } };
		
		  for (int j = 0; j < 4; j++)
		    {
		      Point3d & lp1 = mesh->Point (el.PNum(trigs[j][0]));
		      Point3d & lp2 = mesh->Point (el.PNum(trigs[j][1]));
		      Point3d & lp3 = mesh->Point (el.PNum(trigs[j][2]));
		      Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
		      glNormal3dv (&n.X());

		      glVertex3dv (&lp1.X());
		      glVertex3dv (&lp2.X());
		      glVertex3dv (&lp3.X());
		    }
		  glEnd();
		  break;
		}

	      case QUAD6:
		{
		  glBegin (GL_QUADS);
		  static int quads[2][4] = {
		    { 1, 5, 6, 4 },
		    { 5, 2, 3, 6 } };
		
		  for (int j = 0; j < 2; j++)
		    {
		      Point3d & lp1 = mesh->Point (el.PNum(quads[j][0]));
		      Point3d & lp2 = mesh->Point (el.PNum(quads[j][1]));
		      Point3d & lp3 = mesh->Point (el.PNum(quads[j][2]));
		      Point3d & lp4 = mesh->Point (el.PNum(quads[j][3]));
		      Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
		      n /= (n.Length() + 1e-12);
		      glNormal3dv (&n.X());
		      glVertex3dv (&lp1.X());
		      glVertex3dv (&lp2.X());
		      glVertex3dv (&lp3.X());
		      glVertex3dv (&lp4.X());
		    }
		  glEnd();
		  break;
		}

	      case QUAD8:
		{
		  glBegin (GL_TRIANGLES);
		  static int boundary[] = 
		    { 1, 5, 2, 8, 3, 6, 4, 7, 1 };
		
		  Point3d c(0,0,0);
		  for (int j = 0; j < 4; j++)
		    {
		      Point3d & hp = mesh->Point (el[j]);
		      c.X() -= 0.25 * hp.X();
		      c.Y() -= 0.25 * hp.Y();
		      c.Z() -= 0.25 * hp.Z();
		    }
		  for (int j = 4; j < 8; j++)
		    {
		      Point3d & hp = mesh->Point (el[j]);
		      c.X() += 0.5 * hp.X();
		      c.Y() += 0.5 * hp.Y();
		      c.Z() += 0.5 * hp.Z();
		    }

		  for (int j = 0; j < 8; j++)
		    {
		      Point3d & lp1 = mesh->Point (el.PNum(boundary[j]));
		      Point3d & lp2 = mesh->Point (el.PNum(boundary[j+1]));

		      Vec3d n = Cross (Vec3d (c, lp1), Vec3d (c, lp2));
		      n /= (n.Length() + 1e-12);
		      glNormal3dv (&n.X());
		      glVertex3dv (&lp1.X());
		      glVertex3dv (&lp2.X());
		      glVertex3dv (&c.X());
		    }
		  glEnd();
		  break;
		}


	      default:
		PrintSysError ("Cannot draw (2) surface element of type ", 
			       int(el.GetType()));
	      }
	  }
      }
    glLoadName (0);
    glEndList ();

    // endtime = clock();
    // cout << "BuildFillList time = " << double(endtime - starttime)/CLOCKS_PER_SEC << endl;
  }


  void VisualSceneMesh :: BuildLineList()
  {
    SurfaceElementIndex sei;

    if (!lock)
      {
	lock = new NgLock (mesh->Mutex());
	lock -> Lock();
      }

    linetimestamp = NextTimeStamp();


    bool checkvicinity = (stlgeometry != NULL) && stldoctor.showvicinity;

    if (linelist)
      glDeleteLists (linelist, 1);
  
    linelist = glGenLists (1);
    glNewList (linelist, GL_COMPILE);


    glLineWidth (1.0f);

    int hoplotn = 1 << vispar.subdivisions; 
  
    for (sei = 0; sei < mesh->GetNSE(); sei++)
      {
	const Element2d & el = (*mesh)[sei];
      
	bool drawel = !el.IsDeleted();
	if (checkvicinity)
	  for (int j = 0; j < el.GetNP(); j++)
	    if (!stlgeometry->Vicinity(el.GeomInfoPi(j+1).trignum))
	      drawel = 0;

	if (!drawel)
	  continue;

	switch (el.GetType())
	  {
	  case TRIG:
	    {
	      CurvedElements & curv = mesh->GetCurvedElements();
	      if (curv.IsHighOrder() && curv.IsSurfaceElementCurved(sei))
		{
		  Point<3> xg;
		  glBegin (GL_LINE_LOOP);

		  for (int i = 0; i < hoplotn; i++)
		    {
		      Point<2> xr (double(i) / hoplotn, 0);
		      curv.CalcSurfaceTransformation (xr, sei, xg);
		      glVertex3dv (xg);
		    }
		  for (int i = 0; i < hoplotn; i++)
		    {
		      Point<2> xr (double(hoplotn-i) / hoplotn, double(i)/hoplotn);
		      curv.CalcSurfaceTransformation (xr, sei, xg);
		      glVertex3dv (xg);
		    }
		  for (int i = 0; i < hoplotn; i++)
		    {
		      Point<2> xr (0, double(hoplotn-i) / hoplotn);
		      curv.CalcSurfaceTransformation (xr, sei, xg);
		      glVertex3dv (xg);
		    }

		  glEnd();
		} 
	      else 
		{
		  glBegin (GL_TRIANGLES);

		  const Point<3> & lp0 = (*mesh) [el[0]];
		  const Point<3> & lp1 = (*mesh) [el[1]];
		  const Point<3> & lp2 = (*mesh) [el[2]];

		  glVertex3dv (lp0);
		  glVertex3dv (lp1);
		  glVertex3dv (lp2);

		  glEnd();
		}

	      break;

	    }     
 
	  case QUAD:
	    {
	      CurvedElements & curv = mesh->GetCurvedElements();
	      if (curv.IsHighOrder() && curv.IsSurfaceElementCurved(sei))
		{
		  Point<2> xr;
		  Point<3> xg;

		  glBegin (GL_LINE_STRIP);

		  for (int side = 0; side < 4; side++)
		    {
		      for (int i = 0; i <= hoplotn; i++)
			{
			  switch (side)
			    {
			    case 0:
			      xr(0) = (double) i/hoplotn;
			      xr(1) = 0.;
			      break;
			    case 1:
			      xr(0) = 1.;
			      xr(1) = (double) i/hoplotn;
			      break;
			    case 2:
			      xr(0) = (double) (hoplotn-i)/hoplotn;
			      xr(1) = 1.;
			      break;
			    case 3:
			      xr(0) = 0.;
			      xr(1) = (double) (hoplotn-i)/hoplotn;
			      break;
			    }

			  curv.CalcSurfaceTransformation (xr, sei, xg);
			  glVertex3d (xg(0), xg(1), xg(2));

			}

		    }
		  glEnd();

		} else {

		  glBegin (GL_QUADS);

		  const Point3d & lp1 = mesh->Point (el.PNum(1));
		  const Point3d & lp2 = mesh->Point (el.PNum(2));
		  const Point3d & lp3 = mesh->Point (el.PNum(4));
		  const Point3d & lp4 = mesh->Point (el.PNum(3));
		  Vec3d n = Cross (Vec3d (lp1, lp2),
				   Vec3d (lp1, Center (lp3, lp4)));
		  glNormal3d (n.X(), n.Y(), n.Z());
		  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		  glVertex3d (lp4.X(), lp4.Y(), lp4.Z());
		  glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
		  glEnd();

		}
	    
	      break;
	    
	    }
	    
	  case TRIG6:
	    {
	      int lines[6][2] = {
		{ 1, 6 }, { 2, 6 },
		{ 1, 5 }, { 3, 5 },
		{ 2, 4 }, { 3, 4 } };
	    
	      glBegin (GL_LINES);
	      for (int j = 0; j < 6; j++)
		{
		  const Point3d & lp1 = mesh->Point (el.PNum(lines[j][0]));
		  const Point3d & lp2 = mesh->Point (el.PNum(lines[j][1]));
		
		  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		}
	    
	      glEnd();
	      break;
	    }
	  
	  case QUAD6:
	    {
	      int lines[6][2] = {
		{ 1, 5 }, { 2, 5 },
		{ 3, 6 }, { 4, 6 },
		{ 1, 4 }, { 2, 3 } };
	    
	      glBegin (GL_LINES);
	    
	      for (int j = 0; j < 6; j++)
		{
		  const Point3d & lp1 = mesh->Point (el.PNum(lines[j][0]));
		  const Point3d & lp2 = mesh->Point (el.PNum(lines[j][1]));
		
		  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		}
	      glEnd ();
	      break;
	    }

	  case QUAD8:
	    {
	      int lines[8][2] = {
		{ 1, 5 }, { 2, 5 }, { 3, 6 }, { 4, 6 },
		{ 1, 7 }, { 4, 7 }, { 2, 8 }, { 3, 8 }
	      };
	    
	      glBegin (GL_LINES);
	    
	      for (int j = 0; j < 8; j++)
		{
		  const Point3d & lp1 = mesh->Point (el.PNum(lines[j][0]));
		  const Point3d & lp2 = mesh->Point (el.PNum(lines[j][1]));
		
		  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		}
	      glEnd ();
	      break;
	    }



	  default:
	    PrintSysError ("Cannot draw (4) surface element of type ", 
			   int(el.GetType()));
	  }
      }
    
    glEndList ();
  }

  void VisualSceneMesh :: BuildPointNumberList()
  {
    ;
  }






  inline long int Fact (int n)
  {
    long int res = 1;
    for (int i = 2; i <= n; i++)
      res *= i;
    return res;
  }
  inline int Binom (int n, int i)
  {
    return Fact(n) / Fact(i) / Fact(n-i);
  }
  
  void ToBernstein (int order, Point<3> * pts, int stride)
  {
    static DenseMatrix mat, inv;
    static Vector vec1, vec2;

    if (mat.Height () != order+1)
      {
	mat.SetSize (order+1);
	inv.SetSize (order+1);
	vec1.SetSize (order+1);
	vec2.SetSize (order+1);
	for (int i = 0; i <= order; i++)
	  {
	    double x = double(i) / order;
	    for (int j = 0; j <= order; j++)
	      mat(i,j) = Binom (order, j) * pow (x, j) * pow (1-x, order-j);
	  }

	CalcInverse (mat, inv);
      }

    for (int i = 0; i < 3; i++)
      {
	for (int j = 0; j <= order; j++)
	  vec1(j) = pts[j*stride](i);

	inv.Mult (vec1, vec2);

	for (int j = 0; j <= order; j++)
	  pts[j*stride](i) = vec2(j);
      }
  }














  void VisualSceneMesh :: BuildTetList()
  {


#ifdef FAST3DELEMENTS
    

    cout << "start fast test" << endl;

    int i, j;
    ARRAY<Element2d> faces;

    if (tettimestamp > mesh->GetTimeStamp () &&
	tettimestamp > vispar.clipplanetimestamp )
      return;

    if (!lock)
      {
	lock = new NgLock (mesh->Mutex());
	lock -> Lock();
      }

    tettimestamp = NextTimeStamp();

    if (tetlist)
      glDeleteLists (tetlist, 1);


    tetlist = glGenLists (1);
    glNewList (tetlist, GL_COMPILE);

  
    BitArray shownode(mesh->GetNP());
    if (vispar.clipenable)
      {
	shownode.Clear();
	for (i = 1; i <= shownode.Size(); i++)
	  {
	    Point3d p = mesh->Point(i);
	  
	    double val =
	      p.X() * clipplane[0] +
	      p.Y() * clipplane[1] +
	      p.Z() * clipplane[2] +
	      clipplane[3];
	  
	    if (val > 0)
	      shownode.Set (i);
	  }
      }
    else
      shownode.Set();


    static float tetcols[][4] = 
      {
	{ 1.0f, 1.0f, 0.0f, 1.0f },
	{ 1.0f, 0.0f, 0.0f, 1.0f },
	{ 0.0f, 1.0f, 0.0f, 1.0f },
	{ 0.0f, 0.0f, 1.0f, 1.0f }
      };


    ARRAY<int> elfaces;

    const MeshTopology & top = mesh->GetTopology();
    CurvedElements & curv = mesh->GetCurvedElements();

    ARRAY<int> displayface(top.GetNFaces());
    for (i = 0; i < top.GetNFaces(); i++)
      displayface[i] = -1;
    

    for (i = 1; i <= mesh->GetNE(); i++)
      {
	if (vispar.drawtetsdomain > 0 &&
	    vispar.drawtetsdomain != mesh->VolumeElement(i).GetIndex())
	  continue;

	Element el = (*mesh)[(ElementIndex) (i-1)];

	if (el.GetType() == TET)
	  {
	    if (el.PNum(1))
	      {
		bool drawtet = 1;
		for (j = 1; j <= el.GetNP(); j++)
		  if (!shownode.Test(el.PNum(j)))
		    drawtet = 0;
		if (!drawtet) continue;
		
		{
		  top.GetElementFaces (i, elfaces);
		  
		  for (j = 0; j < 4; j++)
		    displayface[elfaces[j]-1] *= -1;
		}
	      }
	  }
      }

    
    for (i = 1; i <= mesh->GetNE(); i++)
      {
	if (vispar.drawtetsdomain > 0 &&
	    vispar.drawtetsdomain != mesh->VolumeElement(i).GetIndex())
	  continue;

	if (mesh->VolumeElement(i).GetType() == TET)
	  {
	    // copy to be thread-safe
	    Element el = mesh->VolumeElement (i);
	    el.GetSurfaceTriangles (faces);

	    if (el.PNum(1))
	      {
		bool drawtet = 1;
		for (j = 1; j <= el.GetNP(); j++)
		  if (!shownode.Test(el.PNum(j)))
		    drawtet = 0;
		if (!drawtet) continue;

		int ind = el.GetIndex() % 4;
		if (vispar.drawmetispartition && (el.GetPartition()!=-1))
		  ind = el.GetPartition() % 4;
		(*testout) << "ind = " << ind << endl;
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, tetcols[ind]);


		top.GetElementFaces (i, elfaces);

		glBegin (GL_TRIANGLES);

		for (j = 0; j < faces.Size(); j++)
		  {
		    if (displayface[elfaces[j]-1] == -1) continue;

		    Element2d & face = faces.Elem(j+1);

		    if (curv.IsHighOrder() && curv.IsElementCurved(i-1))
		      {
			int hoplotn = 1 << vispar.subdivisions; 
			// int hoplotn = curv.GetNVisualSubsecs();
			    
			const Point3d * facepoint = MeshTopology :: GetVertices (TET);
			const ELEMENT_FACE * elface = MeshTopology :: GetFaces(TET);
                    
			Vec<3> x0,x1,d0,d1;
			Point<3> xg;
			x0 = facepoint[face.PNum(3)-1] - facepoint[face.PNum(1)-1];
			x1 = facepoint[face.PNum(2)-1] - facepoint[face.PNum(1)-1];
			x0.Normalize();
			x1.Normalize();

			for (int m0 = 0; m0 < hoplotn; m0++)
			  for (int m1 = 0; m1 < hoplotn-m0; m1++)
			    for (int k = 0; k < 2; k++)
			      {
				Vec<3> dx, dy, dz, n;
				Point<4> la[3];
				int l;
				for (l = 0; l<3; l++) la[l] = Point<4>(0.,0.,0.,0.);

				if (k == 0)
				  {
				    la[0](face.PNum(1)-1) = (m0  )/(double)hoplotn;
				    la[0](face.PNum(2)-1) = (m1  )/(double)hoplotn;
				    la[0](face.PNum(3)-1) = 1-la[0](face.PNum(1)-1)-la[0](face.PNum(2)-1);

				    la[1](face.PNum(1)-1) = (m0+1)/(double)hoplotn;
				    la[1](face.PNum(2)-1) = (m1  )/(double)hoplotn;
				    la[1](face.PNum(3)-1) = 1-la[1](face.PNum(1)-1)-la[1](face.PNum(2)-1);

				    la[2](face.PNum(1)-1) = (m0  )/(double)hoplotn;
				    la[2](face.PNum(2)-1) = (m1+1)/(double)hoplotn;
				    la[2](face.PNum(3)-1) = 1-la[2](face.PNum(1)-1)-la[2](face.PNum(2)-1);
				  } else
				    {
				      if (m1 == hoplotn-m0-1) continue;
				      la[0](face.PNum(1)-1) = (m0+1)/(double)hoplotn;
				      la[0](face.PNum(2)-1) = (m1+1)/(double)hoplotn;
				      la[0](face.PNum(3)-1) = 1-la[0](face.PNum(1)-1)-la[0](face.PNum(2)-1);

				      la[1](face.PNum(1)-1) = (m0  )/(double)hoplotn;
				      la[1](face.PNum(2)-1) = (m1+1)/(double)hoplotn;
				      la[1](face.PNum(3)-1) = 1-la[1](face.PNum(1)-1)-la[1](face.PNum(2)-1);

				      la[2](face.PNum(1)-1) = (m0+1)/(double)hoplotn;
				      la[2](face.PNum(2)-1) = (m1  )/(double)hoplotn;
				      la[2](face.PNum(3)-1) = 1-la[2](face.PNum(1)-1)-la[2](face.PNum(2)-1);
				    }

				for (l = 0; l<3; l++)
				  {
				    Mat<3,3> dxdxi;
				    Point<3> xr( la[l](0), la[l](1), la[l](2) );
				    curv.CalcElementTransformation (xr, i-1, xg, dxdxi);
				    for (int i = 0; i < 3; i++)
				      {
					dx(i) = dxdxi(i,0);
					dy(i) = dxdxi(i,1);
					dz(i) = dxdxi(i,2);
				      }

				    d0 = x0(0)*dx + x0(1)*dy + x0(2)*dz;
				    d1 = x1(0)*dx + x1(1)*dy + x1(2)*dz;
				    n = Cross (d0, d1);
				    glNormal3d ( n(0),  n(1),  n(2));
				    glVertex3d (xg(0), xg(1), xg(2));
				  }
			      }

		      } else {
			const Point3d & lp1 = mesh->Point (el.PNum(face.PNum(1)));
			const Point3d & lp2 = mesh->Point (el.PNum(face.PNum(2)));
			const Point3d & lp3 = mesh->Point (el.PNum(face.PNum(3)));
			Vec3d n = Cross (Vec3d (lp1, lp3), Vec3d (lp1, lp2));
			n /= (n.Length()+1e-12);
			glNormal3d (n.X(), n.Y(), n.Z());
			glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
			glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
			glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
		      }
		  }
	     
		glEnd();
	      }

	  }
      }
    glEndList ();




#else

    if (tettimestamp > mesh->GetTimeStamp () &&
	tettimestamp > vispar.clipplanetimestamp )
      return;

    if (!lock)
      {
	lock = new NgLock (mesh->Mutex());
	lock -> Lock();
      }

    tettimestamp = NextTimeStamp();

    if (tetlist)
      glDeleteLists (tetlist, 1);


    tetlist = glGenLists (1);
    glNewList (tetlist, GL_COMPILE);



    int i, j, k, l;
    ARRAY<Element2d> faces;

  
    BitArray shownode(mesh->GetNP());
    if (vispar.clipenable)
      {
	shownode.Clear();
	for (i = 1; i <= shownode.Size(); i++)
	  {
	    Point3d p = mesh->Point(i);
	  
	    double val =
	      p.X() * clipplane[0] +
	      p.Y() * clipplane[1] +
	      p.Z() * clipplane[2] +
	      clipplane[3];
	  
	    if (val > 0)
	      shownode.Set (i);
	  }
      }
    else
      shownode.Set();


    static float tetcols[][4] = 
      {
	{ 1.0f, 1.0f, 0.0f, 1.0f },
	{ 1.0f, 0.0f, 0.0f, 1.0f },
	{ 0.0f, 1.0f, 0.0f, 1.0f },
	{ 0.0f, 0.0f, 1.0f, 1.0f }
      };


    CurvedElements & curv = mesh->GetCurvedElements();

    for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
      {
	i = ei + 1;

	if (vispar.drawtetsdomain > 0 &&
	    vispar.drawtetsdomain != mesh->VolumeElement(i).GetIndex())
	  continue;

	const Element & el = (*mesh)[ei];

	if (el.GetType() == TET && !el.IsDeleted())
	  {

	    bool drawtet = 1;
	    for (j = 0; j < 4; j++)
	      if (!shownode.Test(el[j]))
		drawtet = 0;
	    if (!drawtet) continue;
	    
	    
	    int ind = el.GetIndex() % 4;
	    if (vispar.drawmetispartition && (el.GetPartition()!=-1))
	      ind = el.GetPartition() % 4;
	    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, tetcols[ind]);
	    
	    
	    if (curv.IsHighOrder() && curv.IsElementCurved(ei))
	      {
		const ELEMENT_FACE * faces = MeshTopology :: GetFaces (TET);
		const Point3d * vertices = MeshTopology :: GetVertices (TET);
		
		Point<3> grid[11][11];
		Point<3> fpts[3];
		int order = vispar.subdivisions+1;
		
		for (int trig = 0; trig < 4; trig++)
		  {   
		    for (int j = 0; j < 3; j++)
		      fpts[j] = vertices[faces[trig][j]-1];
		    
		    static Point<3> c(0.25, 0.25, 0.25);
		    if (vispar.shrink < 1)
		      for (int j = 0; j < 3; j++)
			fpts[j] += (1-vispar.shrink) * (c-fpts[j]);
			
		    for (int ix = 0; ix <= order; ix++)
		      for (int iy = 0; iy <= order; iy++)
			{
			  double lami[3] = 
			    { (1-double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (1-double(iy)/order),
			      double(iy)/order };
			      
			  Point<3> xl;
			  for (int l = 0; l < 3; l++)
			    xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
			      lami[2] * fpts[2](l);
			      
			  curv.CalcElementTransformation (xl, i-1, grid[ix][iy]);
			}
			
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);
			
		    glMap2d(GL_MAP2_VERTEX_3, 
			    0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
			    0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
			    &grid[0][0](0));
		    glEnable(GL_MAP2_VERTEX_3);
		    glEnable(GL_AUTO_NORMAL); 
			
		    glMapGrid2f(8, 0.0, 0.999, 8, 0.0, 1.0);
		    glEvalMesh2(GL_FILL, 0, 8, 0, 8);
			
		    glDisable (GL_AUTO_NORMAL);
		    glDisable (GL_MAP2_VERTEX_3);
		  }		    
	      }


	    else

	      {
		el.GetSurfaceTriangles (faces);

		Point3d c;
		if (vispar.shrink < 1)
		  c = Center (Center (mesh->Point (el.PNum(1)),
				      mesh->Point (el.PNum(2))),
			      Center (mesh->Point (el.PNum(3)),
				      mesh->Point (el.PNum(4))));


		glBegin (GL_TRIANGLES);

		for (j = 0; j < faces.Size(); j++)
		  {
		    const Element2d & face = faces[j];

		    /*
		      if (curv.IsHighOrder() && curv.IsElementCurved(i-1))
		      {
		      int hoplotn = 1 << vispar.subdivisions; 
		      // int hoplotn = curv.GetNVisualSubsecs();
			    
		      const Point3d * facepoint = MeshTopology :: GetVertices (TET);
		      const ELEMENT_FACE * elface = MeshTopology :: GetFaces(TET);
                    
		      Vec<3> x0,x1,d0,d1;
		      Point<3> xg;
		      x0 = facepoint[face.PNum(3)-1] - facepoint[face.PNum(1)-1];
		      x1 = facepoint[face.PNum(2)-1] - facepoint[face.PNum(1)-1];
		      x0.Normalize();
		      x1.Normalize();

		      for (int m0 = 0; m0 < hoplotn; m0++)
		      for (int m1 = 0; m1 < hoplotn-m0; m1++)
		      for (k = 0; k < 2; k++)
		      {
		      Vec<3> dx, dy, dz, n;
		      Point<4> la[3];
					
		      for (l = 0; l<3; l++) la[l] = Point<4>(0.,0.,0.,0.);

		      if (k == 0)
		      {
		      la[0](face.PNum(1)-1) = (m0  )/(double)hoplotn;
		      la[0](face.PNum(2)-1) = (m1  )/(double)hoplotn;
		      la[0](face.PNum(3)-1) = 1-la[0](face.PNum(1)-1)-la[0](face.PNum(2)-1);

		      la[1](face.PNum(1)-1) = (m0+1)/(double)hoplotn;
		      la[1](face.PNum(2)-1) = (m1  )/(double)hoplotn;
		      la[1](face.PNum(3)-1) = 1-la[1](face.PNum(1)-1)-la[1](face.PNum(2)-1);

		      la[2](face.PNum(1)-1) = (m0  )/(double)hoplotn;
		      la[2](face.PNum(2)-1) = (m1+1)/(double)hoplotn;
		      la[2](face.PNum(3)-1) = 1-la[2](face.PNum(1)-1)-la[2](face.PNum(2)-1);
		      } else
		      {
		      if (m1 == hoplotn-m0-1) continue;
		      la[0](face.PNum(1)-1) = (m0+1)/(double)hoplotn;
		      la[0](face.PNum(2)-1) = (m1+1)/(double)hoplotn;
		      la[0](face.PNum(3)-1) = 1-la[0](face.PNum(1)-1)-la[0](face.PNum(2)-1);

		      la[1](face.PNum(1)-1) = (m0  )/(double)hoplotn;
		      la[1](face.PNum(2)-1) = (m1+1)/(double)hoplotn;
		      la[1](face.PNum(3)-1) = 1-la[1](face.PNum(1)-1)-la[1](face.PNum(2)-1);

		      la[2](face.PNum(1)-1) = (m0+1)/(double)hoplotn;
		      la[2](face.PNum(2)-1) = (m1  )/(double)hoplotn;
		      la[2](face.PNum(3)-1) = 1-la[2](face.PNum(1)-1)-la[2](face.PNum(2)-1);
		      }

		      for (l = 0; l<3; l++)
		      {
		      Mat<3,3> dxdxi;
		      Point<3> xr( la[l](0), la[l](1), la[l](2) );
		      curv.CalcElementTransformation (xr, i-1, xg, dxdxi);
		      for (int i = 0; i < 3; i++)
		      {
		      dx(i) = dxdxi(i,0);
		      dy(i) = dxdxi(i,1);
		      dz(i) = dxdxi(i,2);
		      }

		      d0 = x0(0)*dx + x0(1)*dy + x0(2)*dz;
		      d1 = x1(0)*dx + x1(1)*dy + x1(2)*dz;
		      n = Cross (d0, d1);
		      glNormal3d ( n(0),  n(1),  n(2));
		      glVertex3d (xg(0), xg(1), xg(2));
		      }
		      }

		      } else 
		    */
		    {
		      Point3d lp1 = mesh->Point (el.PNum(face.PNum(1)));
		      Point3d lp2 = mesh->Point (el.PNum(face.PNum(2)));
		      Point3d lp3 = mesh->Point (el.PNum(face.PNum(3)));
		      Vec3d n = Cross (Vec3d (lp1, lp3), Vec3d (lp1, lp2));
		      n /= (n.Length()+1e-12);
		      glNormal3d (n.X(), n.Y(), n.Z());

		      if (vispar.shrink < 1)
			{
			  lp1 = c + vispar.shrink * (lp1 - c);
			  lp2 = c + vispar.shrink * (lp2 - c);
			  lp3 = c + vispar.shrink * (lp3 - c);
			}

		      glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		      glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		      glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
		    }
		  }
	     
		glEnd();
	      }
	  }
      }
    
    glEndList ();
    
#endif
  }




  void VisualSceneMesh :: BuildPrismList()
  {
    if (prismtimestamp > mesh->GetTimeStamp () &&
	prismtimestamp > vispar.clipplanetimestamp )
      return;

    if (!lock)
      {
	lock = new NgLock (mesh->Mutex());
	lock -> Lock();
      }

    prismtimestamp = NextTimeStamp();



    if (prismlist)
      glDeleteLists (prismlist, 1);

    prismlist = glGenLists (1);
    glNewList (prismlist, GL_COMPILE);

    static float prismcol[] = { 0.0f, 1.0f, 1.0f, 1.0f };
    glLineWidth (1.0f);  

    ARRAY<Element2d> faces;

    
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, prismcol);

    for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
      {
	const Element & el = (*mesh)[ei];
	if (el.GetType() == PRISM && !el.IsDeleted())
	  {
	    int j;
	    int i = ei + 1;

	    CurvedElements & curv = mesh->GetCurvedElements();
	    if (curv.IsHighOrder() && curv.IsElementCurved(ei))
	      {
		const ELEMENT_FACE * faces = MeshTopology :: GetFaces (PRISM);
		const Point3d * vertices = MeshTopology :: GetVertices (PRISM);
		
		Point<3> grid[11][11];
		Point<3> fpts[4];
		int order = vispar.subdivisions+1;
		
		for (int trig = 0; trig < 2; trig++)
		  {   
		    for (int j = 0; j < 3; j++)
		      fpts[j] = vertices[faces[trig][j]-1];
		    
		    static Point<3> c(1.0/3.0, 1.0/3.0, 0.5);
		    if (vispar.shrink < 1)
		      for (int j = 0; j < 3; j++)
			fpts[j] += (1-vispar.shrink) * (c-fpts[j]);
			
		    for (int ix = 0; ix <= order; ix++)
		      for (int iy = 0; iy <= order; iy++)
			{
			  double lami[3] = 
			    { (1-double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (1-double(iy)/order),
			      double(iy)/order };
			      
			  Point<3> xl;
			  for (int l = 0; l < 3; l++)
			    xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
			      lami[2] * fpts[2](l);
			      
			  curv.CalcElementTransformation (xl, i-1, grid[ix][iy]);
			}
			
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);
			
		    glMap2d(GL_MAP2_VERTEX_3, 
			    0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
			    0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
			    &grid[0][0](0));
		    glEnable(GL_MAP2_VERTEX_3);
		    glEnable(GL_AUTO_NORMAL); 
			
		    glMapGrid2f(8, 0.0, 0.999, 8, 0.0, 1.0);
		    glEvalMesh2(GL_FILL, 0, 8, 0, 8);
			
		    glDisable (GL_AUTO_NORMAL);
		    glDisable (GL_MAP2_VERTEX_3);
		  }		    

		for (int quad = 2; quad < 5; quad++)
		  {   
		    for (int j = 0; j < 4; j++)
		      fpts[j] = vertices[faces[quad][j]-1];
		    
		    static Point<3> c(1.0/3.0, 1.0/3.0, 0.5);
		    if (vispar.shrink < 1)
		      for (int j = 0; j < 4; j++)
			fpts[j] += (1-vispar.shrink) * (c-fpts[j]);
			
		    for (int ix = 0; ix <= order; ix++)
		      for (int iy = 0; iy <= order; iy++)
			{
			  double lami[4] = 
			    { (1-double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (  double(iy)/order),
			      (1-double(ix)/order) * (  double(iy)/order) };
			      
			  Point<3> xl;
			  for (int l = 0; l < 3; l++)
			    xl(l) = 
			      lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
			      lami[2] * fpts[2](l) + lami[3] * fpts[3](l);
			      
			  curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
			}
			
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);
			
		    glMap2d(GL_MAP2_VERTEX_3, 
			    0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
			    0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
			    &grid[0][0](0));
		    glEnable(GL_MAP2_VERTEX_3);
		    glEnable(GL_AUTO_NORMAL); 
			
		    glMapGrid2f(8, 0.0, 1.0, 8, 0.0, 1.0);
		    glEvalMesh2(GL_FILL, 0, 8, 0, 8);
			
		    glDisable (GL_AUTO_NORMAL);
		    glDisable (GL_MAP2_VERTEX_3);
		  }		    





		/*
		int hoplotn = 1 << vispar.subdivisions; 
		// int hoplotn = curv.GetNVisualSubsecs();

		const Point3d * facepoint = MeshTopology :: GetVertices (TRIG);
		const ELEMENT_FACE * elface = MeshTopology :: GetFaces(TRIG);

		glBegin (GL_TRIANGLES);

		for (int trig = 0; trig<2; trig++)
		  {
                    
		    Vec<3> x0,x1,d0,d1;
		    x0 = facepoint[1] - facepoint[2];
		    x1 = facepoint[0] - facepoint[2];
		    x0.Normalize();
		    x1.Normalize();
		    if (trig == 1) swap (x0,x1);

		    Point<3> xr[3];
		    Point<3> xg;
		    Vec<3> dx, dy, dz, n;

		    for (int i1 = 0; i1 < hoplotn; i1++)
		      for (int j1 = 0; j1 < hoplotn-i1; j1++)
			for (int k = 0; k < 2; k++)
			  {
			    if (k == 0)
			      {
				xr[0](0) = (double)    i1/hoplotn; xr[0](1) = (double)    j1/hoplotn;
				xr[1](0) = (double)(i1+1)/hoplotn; xr[1](1) = (double)    j1/hoplotn;
				xr[2](0) = (double)    i1/hoplotn; xr[2](1) = (double)(j1+1)/hoplotn;
			      } else
				{
				  if (j1 == hoplotn-i1-1) continue;
				  xr[0](0) = (double)(i1+1)/hoplotn; xr[0](1) = (double)    j1/hoplotn;
				  xr[1](0) = (double)(i1+1)/hoplotn; xr[1](1) = (double)(j1+1)/hoplotn;
				  xr[2](0) = (double)    i1/hoplotn; xr[2](1) = (double)(j1+1)/hoplotn;
				};
				    
			    for (int l=0; l<3; l++)
			      {
				Mat<3,3> dxdxi;
				xr[l](2) = (double) trig;
				curv.CalcElementTransformation (xr[l], i-1, xg, dxdxi);
				for (int i = 0; i < 3; i++)
				  {
				    dx(i) = dxdxi(i,0);
				    dy(i) = dxdxi(i,1);
				    dz(i) = dxdxi(i,2);
				  }

				Vec<3> d0 = x0(0)*dx + x0(1)*dy + x0(2)*dz;
				Vec<3> d1 = x1(0)*dx + x1(1)*dy + x1(2)*dz;
				n = Cross (d1, d0);
				glNormal3d (n(0), n(1), n(2));
				glVertex3d (xg(0), xg(1), xg(2));
			      }
			  }
			
		  }

		glEnd ();

		glBegin (GL_QUADS);

		for (int quad = 0; quad<3; quad++)
		  {   
		    const Point3d * facepoint = MeshTopology :: GetVertices (PRISM);

		    Vec<3> x0,x1;
		    int xyz;

		    switch (quad)
		      {
		      case 0:
			x0 = facepoint[5] - facepoint[2];
			x1 = facepoint[0] - facepoint[2];
			xyz = 0;
			break;
		      case 1:
			x0 = facepoint[4] - facepoint[0];
			x1 = facepoint[1] - facepoint[0];
			xyz = 0;
			break;
		      case 2:
			x0 = facepoint[1] - facepoint[2];
			x1 = facepoint[5] - facepoint[2];
			xyz = 1;
			break;
		      }

		    x0.Normalize();
		    x1.Normalize();

		    swap (x0,x1);

		    Point<3> xr[4];
		    Point<3> xg;
		    Vec<3> dx, dy, dz, n;

		    for (int i1 = 0; i1 < hoplotn; i1++)
		      for (int j1 = 0; j1 < hoplotn; j1++)
			{
			  xr[0](xyz) = (double)    i1/hoplotn; xr[0](2) = (double)    j1/hoplotn;
			  xr[1](xyz) = (double)(i1+1)/hoplotn; xr[1](2) = (double)    j1/hoplotn;
			  xr[2](xyz) = (double)(i1+1)/hoplotn; xr[2](2) = (double)(j1+1)/hoplotn;
			  xr[3](xyz) = (double)    i1/hoplotn; xr[3](2) = (double)(j1+1)/hoplotn;
				    
			  for (int l=0; l<4; l++)
			    {
			      switch (quad)
				{
				case 0: xr[l](1) = 0; break;
				case 1: xr[l](1) = 1-xr[l](0); break;
				case 2: xr[l](0) = 0; break;
				}

			      Mat<3,3> dxdxi;
			      curv.CalcElementTransformation (xr[l], i-1, xg, dxdxi);
			      for (int i = 0; i < 3; i++)
				{
				  dx(i) = dxdxi(i,0);
				  dy(i) = dxdxi(i,1);
				  dz(i) = dxdxi(i,2);
				}

			      Vec<3> d0 = x0(0)*dx + x0(1)*dy + x0(2)*dz;
			      Vec<3> d1 = x1(0)*dx + x1(1)*dy + x1(2)*dz;
			      n = Cross (d1, d0);
			      glNormal3d (n(0), n(1), n(2));
			      glVertex3d (xg(0), xg(1), xg(2));
			    }
			}
		  }
		glEnd ();
		*/
	      } 
	    else
	      { 
		Point3d c(0,0,0);
		if (vispar.shrink < 1)
		  {
		    for (j = 1; j <= 6; j++)
		      {
			Point3d p = mesh->Point(el.PNum(j));
			c.X() += p.X() / 6;
			c.Y() += p.Y() / 6;
			c.Z() += p.Z() / 6;
		      }
		  }

		el.GetSurfaceTriangles (faces);
		glBegin (GL_TRIANGLES);
		for (j = 1; j <= faces.Size(); j++)
		  {
		    Element2d & face = faces.Elem(j);
		    Point3d lp1 = mesh->Point (el.PNum(face.PNum(1)));
		    Point3d lp2 = mesh->Point (el.PNum(face.PNum(2)));
		    Point3d lp3 = mesh->Point (el.PNum(face.PNum(3)));
		    Vec3d n = Cross (Vec3d (lp1, lp3), Vec3d (lp1, lp2));
		    n /= (n.Length()+1e-12);
		    glNormal3d (n.X(), n.Y(), n.Z());
		    if (vispar.shrink < 1)
		      {
			lp1 = c + vispar.shrink * (lp1 - c);
			lp2 = c + vispar.shrink * (lp2 - c);
			lp3 = c + vispar.shrink * (lp3 - c);
		      }
		    glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		    glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		    glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
		  }
	      
		glEnd();
	      }
	  }
      }
    glEndList ();
  }
 



  void VisualSceneMesh :: BuildHexList()
  {
    if (hextimestamp > mesh->GetTimeStamp () &&
	hextimestamp > vispar.clipplanetimestamp )
      return;

    if (!lock)
      {
	lock = new NgLock (mesh->Mutex());
	lock -> Lock();
      }

    hextimestamp = NextTimeStamp();

    if (hexlist) glDeleteLists (hexlist, 1);

    hexlist = glGenLists (1);
    glNewList (hexlist, GL_COMPILE);


    static float hexcol[] = { 1.0f, 1.0f, 0.0f, 1.0f };
    glLineWidth (1.0f);  
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, hexcol);

    ARRAY<Element2d> faces;
    int hoplotn = 1 << vispar.subdivisions; 

    for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
      {
	const Element & el = (*mesh)[ei];
	if (el.GetType() == HEX && !el.IsDeleted())
	  {
	    CurvedElements & curv = mesh->GetCurvedElements();
	    if (curv.IsHighOrder() && curv.IsElementCurved(ei))
	      {
		/* // classical 
		   glBegin (GL_QUADS);
		
		   const ELEMENT_FACE * faces = MeshTopology :: GetFaces (HEX);
		   const Point3d * vertices = MeshTopology :: GetVertices (HEX);

		   Point<3> grid[33][33];
		   Vec<3> gridn[33][33];
		   Point<3> fpts[4];
		   for (int quad = 0; quad<6; quad++)
		   {   
		   for (int j = 0; j < 4; j++)
		   fpts[j] = vertices[faces[quad][j]-1];

		   static Point<3> c(0.5, 0.5, 0.5);
		   if (vispar.shrink < 1)
		   for (int j = 0; j < 4; j++)
		   fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

		   Vec<3> taux = fpts[1]-fpts[0];
		   Vec<3> tauy = fpts[3]-fpts[0];

		   for (int ix = 0; ix <= hoplotn; ix++)
		   for (int iy = 0; iy <= hoplotn; iy++)
		   {
		   Point<3> xl;
		   Mat<3,3> dxdxi;
		   double lami[4] = 
		   { (1-double(ix)/hoplotn) * (1-double(iy)/hoplotn),
		   (  double(ix)/hoplotn) * (1-double(iy)/hoplotn),
		   (  double(ix)/hoplotn) * (  double(iy)/hoplotn),
		   (1-double(ix)/hoplotn) * (  double(iy)/hoplotn) };
		   for (int l = 0; l < 3; l++)
		   xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
		   lami[2] * fpts[2](l) + lami[3] * fpts[3](l);

		   curv.CalcElementTransformation (xl, ei, grid[ix][iy], dxdxi);
			  
		   Vec<3> gtaux = dxdxi * taux;
		   Vec<3> gtauy = dxdxi * tauy;
		   gridn[ix][iy] = Cross (gtauy, gtaux).Normalize();
		   }
		    
		   for (int ix = 0; ix < hoplotn; ix++)
		   for (int iy = 0; iy < hoplotn; iy++)
		   {
		   glNormal3dv (gridn[ix][iy]);
		   glVertex3dv (grid[ix][iy]);

		   glNormal3dv (gridn[ix+1][iy]);
		   glVertex3dv (grid[ix+1][iy]);

		   glNormal3dv (gridn[ix+1][iy+1]);
		   glVertex3dv (grid[ix+1][iy+1]);

		   glNormal3dv (gridn[ix][iy+1]);
		   glVertex3dv (grid[ix][iy+1]);
		   }
		   }
		
		   glEnd ();
		*/

		const ELEMENT_FACE * faces = MeshTopology :: GetFaces (HEX);
		const Point3d * vertices = MeshTopology :: GetVertices (HEX);

		Point<3> grid[11][11];
		Point<3> fpts[4];
		int order = vispar.subdivisions+1;

		for (int quad = 0; quad<6; quad++)
		  {   
		    for (int j = 0; j < 4; j++)
		      fpts[j] = vertices[faces[quad][j]-1];

		    static Point<3> c(0.5, 0.5, 0.5);
		    if (vispar.shrink < 1)
		      for (int j = 0; j < 4; j++)
			fpts[j] += (1-vispar.shrink) * (c-fpts[j]);

		    for (int ix = 0; ix <= order; ix++)
		      for (int iy = 0; iy <= order; iy++)
			{
			  double lami[4] = 
			    { (1-double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (  double(iy)/order),
			      (1-double(ix)/order) * (  double(iy)/order) };

			  Point<3> xl;
			  for (int l = 0; l < 3; l++)
			    xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
			      lami[2] * fpts[2](l) + lami[3] * fpts[3](l);

			  curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
			}
		    
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);

		    glMap2d(GL_MAP2_VERTEX_3, 
			    0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
			    0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
			    &grid[0][0](0));
		    glEnable(GL_MAP2_VERTEX_3);
		    glEnable(GL_AUTO_NORMAL); 

		    glMapGrid2f(8, 0.0, 1.0, 8, 0.0, 1.0);
		    glEvalMesh2(GL_FILL, 0, 8, 0, 8);

		    glDisable (GL_AUTO_NORMAL);
		    glDisable (GL_MAP2_VERTEX_3);
		  }
	      } 
	    else
	      { 
		Point3d c(0,0,0);
		if (vispar.shrink < 1)
		  {
		    for (int j = 1; j <= 8; j++)
		      {
			Point3d p = mesh->Point(el.PNum(j));
			c.X() += p.X();
			c.Y() += p.Y();
			c.Z() += p.Z();
		      }
		    c.X() /= 8;
		    c.Y() /= 8;
		    c.Z() /= 8;
		  }

		glBegin (GL_TRIANGLES);

		el.GetSurfaceTriangles (faces);
		for (int j = 1; j <= faces.Size(); j++)
		  {
		    Element2d & face = faces.Elem(j);
		    Point<3> lp1 = mesh->Point (el.PNum(face.PNum(1)));
		    Point<3> lp2 = mesh->Point (el.PNum(face.PNum(2)));
		    Point<3> lp3 = mesh->Point (el.PNum(face.PNum(3)));
		    Vec<3> n = Cross (lp3-lp1, lp2-lp1);  
		    n.Normalize();
		    glNormal3dv (n);
		  
		    if (vispar.shrink < 1)
		      {
			lp1 = c + vispar.shrink * (lp1 - c);
			lp2 = c + vispar.shrink * (lp2 - c);
			lp3 = c + vispar.shrink * (lp3 - c);
		      }

		    glVertex3dv (lp1);
		    glVertex3dv (lp2);
		    glVertex3dv (lp3);
		  }

		glEnd();
	      }
	  }
      }
    glEndList ();
  }








 
  void VisualSceneMesh :: BuildPyramidList()
  {
    if (pyramidtimestamp > mesh->GetTimeStamp () &&
	pyramidtimestamp > vispar.clipplanetimestamp )
      return;

    if (!lock)
      {
	lock = new NgLock (mesh->Mutex());
	lock -> Lock();
      }

    pyramidtimestamp = NextTimeStamp();


    if (pyramidlist)
      glDeleteLists (pyramidlist, 1);




    pyramidlist = glGenLists (1);
    glNewList (pyramidlist, GL_COMPILE);
  
    static float pyramidcol[] = { 1.0f, 0.0f, 1.0f, 1.0f };
    glLineWidth (1.0f);
    ARRAY<Element2d> faces;

    for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
      {
	const Element & el = (*mesh)[ei];
	if (el.GetType() == PYRAMID && !el.IsDeleted())
	  {
	    int j;
	    int i = ei + 1;

	    CurvedElements & curv = mesh->GetCurvedElements();
	    if (curv.IsHighOrder() && curv.IsElementCurved(ei))
	      {

		const ELEMENT_FACE * faces = MeshTopology :: GetFaces (PYRAMID);
		const Point3d * vertices = MeshTopology :: GetVertices (PYRAMID);
		
		Point<3> grid[11][11];
		Point<3> fpts[4];
		int order = vispar.subdivisions+1;
		
		for (int trig = 0; trig < 4; trig++)
		  {   
		    for (int j = 0; j < 3; j++)
		      fpts[j] = vertices[faces[trig][j]-1];
		    
		    static Point<3> c(0.375, 0.375, 0.25);
		    if (vispar.shrink < 1)
		      for (int j = 0; j < 3; j++)
			fpts[j] += (1-vispar.shrink) * (c-fpts[j]);
			
		    for (int ix = 0; ix <= order; ix++)
		      for (int iy = 0; iy <= order; iy++)
			{
			  double lami[3] = 
			    { (1-double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (1-double(iy)/order),
			      double(iy)/order };
			      
			  Point<3> xl;
			  for (int l = 0; l < 3; l++)
			    xl(l) = lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
			      lami[2] * fpts[2](l);
			      
			  curv.CalcElementTransformation (xl, i-1, grid[ix][iy]);
			}
			
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);
			
		    glMap2d(GL_MAP2_VERTEX_3, 
			    0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
			    0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
			    &grid[0][0](0));
		    glEnable(GL_MAP2_VERTEX_3);
		    glEnable(GL_AUTO_NORMAL); 
			
		    glMapGrid2f(8, 0.0, 0.999, 8, 0.0, 1.0);
		    glEvalMesh2(GL_FILL, 0, 8, 0, 8);
			
		    glDisable (GL_AUTO_NORMAL);
		    glDisable (GL_MAP2_VERTEX_3);
		  }		    

		for (int quad = 4; quad < 5; quad++)
		  {   
		    for (int j = 0; j < 4; j++)
		      fpts[j] = vertices[faces[quad][j]-1];
		    
		    static Point<3> c(0.375, 0.375, 0.25);
		    if (vispar.shrink < 1)
		      for (int j = 0; j < 4; j++)
			fpts[j] += (1-vispar.shrink) * (c-fpts[j]);
			
		    for (int ix = 0; ix <= order; ix++)
		      for (int iy = 0; iy <= order; iy++)
			{
			  double lami[4] = 
			    { (1-double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (1-double(iy)/order),
			      (  double(ix)/order) * (  double(iy)/order),
			      (1-double(ix)/order) * (  double(iy)/order) };
			      
			  Point<3> xl;
			  for (int l = 0; l < 3; l++)
			    xl(l) = 
			      lami[0] * fpts[0](l) + lami[1] * fpts[1](l) +
			      lami[2] * fpts[2](l) + lami[3] * fpts[3](l);
			      
			  curv.CalcElementTransformation (xl, ei, grid[ix][iy]);
			}
			
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[j][0], &grid[0][1]-&grid[0][0]);
		    for (int j = 0; j <= order; j++)
		      ToBernstein (order, &grid[0][j], &grid[1][0]-&grid[0][0]);
			
		    glMap2d(GL_MAP2_VERTEX_3, 
			    0.0, 1.0, &grid[0][1](0)-&grid[0][0](0), order+1,
			    0.0, 1.0, &grid[1][0](0)-&grid[0][0](0), order+1,
			    &grid[0][0](0));
		    glEnable(GL_MAP2_VERTEX_3);
		    glEnable(GL_AUTO_NORMAL); 
			
		    glMapGrid2f(8, 0.0, 1.0, 8, 0.0, 1.0);
		    glEvalMesh2(GL_FILL, 0, 8, 0, 8);
			
		    glDisable (GL_AUTO_NORMAL);
		    glDisable (GL_MAP2_VERTEX_3);
		  }		    






		/*
		int hoplotn = 1 << vispar.subdivisions; 

		const ELEMENT_FACE * faces = MeshTopology :: GetFaces (PYRAMID);
		const Point3d * vertices = MeshTopology :: GetVertices (PYRAMID);

		Point<3> grid[33][33];
		Vec<3> gridn[33][33];


		glBegin (GL_TRIANGLES);
		
		for (int trig = 0; trig < 4; trig++)
		  {   
		    Point<3> p0 = vertices[faces[trig][0]-1];
		    Point<3> p1 = vertices[faces[trig][1]-1];
		    Point<3> p2 = vertices[faces[trig][2]-1];

		    if (vispar.shrink < 1)
		      {
			static Point<3> c(0.375, 0.375, 0.25);
			p0 = c + vispar.shrink * (p0 - c);
			p1 = c + vispar.shrink * (p1 - c);
			p2 = c + vispar.shrink * (p2 - c);
		      }
		   

		    Vec<3> taux = p0-p2;
		    Vec<3> tauy = p1-p2;
		    Vec<3> gtaux, gtauy;

		    Point<3> xl;
		    Mat<3,3> dxdxi;

		    for (int ix = 0; ix <= hoplotn; ix++)
		      for (int iy = 0; iy <= hoplotn-ix; iy++)
			{
			  for (int l = 0; l < 3; l++)
			    xl(l) = 
			      (1-double(ix+iy)/hoplotn) * p2(l) +
			      (double(ix)/hoplotn) * p0(l) +
			      (double(iy)/hoplotn) * p1(l);
			  
			  curv.CalcElementTransformation (xl, i-1, grid[ix][iy], dxdxi);
			  
			  gtaux = dxdxi * taux;
			  gtauy = dxdxi * tauy;
			  gridn[ix][iy] = Cross (gtauy, gtaux).Normalize();
			}

		    for (int ix = 0; ix < hoplotn; ix++)
		      for (int iy = 0; iy < hoplotn-ix; iy++)
			{
			  glNormal3dv (gridn[ix][iy]);
			  glVertex3dv (grid[ix][iy]);

			  glNormal3dv (gridn[ix+1][iy]);
			  glVertex3dv (grid[ix+1][iy]);

			  glNormal3dv (gridn[ix][iy+1]);
			  glVertex3dv (grid[ix][iy+1]);

			  if (iy < hoplotn-ix-1)
			    {
			      glNormal3dv (gridn[ix][iy+1]);
			      glVertex3dv (grid[ix][iy+1]);
			      
			      glNormal3dv (gridn[ix+1][iy]);
			      glVertex3dv (grid[ix+1][iy]);

			      glNormal3dv (gridn[ix+1][iy+1]);
			      glVertex3dv (grid[ix+1][iy+1]);
			    }
			}
		  }

		glEnd ();




		glBegin (GL_QUADS);
		
		for (int quad = 4; quad < 5; quad++)
		  {   
		    Point<3> p0 = vertices[faces[quad][0]-1];
		    Point<3> p1 = vertices[faces[quad][1]-1];
		    Point<3> p2 = vertices[faces[quad][2]-1];
		    Point<3> p3 = vertices[faces[quad][3]-1];

		    if (vispar.shrink < 1)
		      {
			static Point<3> c(0.375, 0.375, 0.25);
			p0 = c + vispar.shrink * (p0 - c);
			p1 = c + vispar.shrink * (p1 - c);
			p2 = c + vispar.shrink * (p2 - c);
			p3 = c + vispar.shrink * (p3 - c);
		      }

		    Vec<3> taux = p1-p0;
		    Vec<3> tauy = p3-p0;
		    Vec<3> gtaux, gtauy;

		    Point<3> xl, xg;
		    Mat<3,3> dxdxi;

		    for (int ix = 0; ix <= hoplotn; ix++)
		      for (int iy = 0; iy <= hoplotn; iy++)
			{
			  Point<3> xl;
			  for (int l = 0; l < 3; l++)
			    xl(l) = 
			      (1-double(ix)/hoplotn)*(1-double(iy)/hoplotn) * p0(l) +
			      (  double(ix)/hoplotn)*(1-double(iy)/hoplotn) * p1(l) +
			      (  double(ix)/hoplotn)*(  double(iy)/hoplotn) * p2(l) +
			      (1-double(ix)/hoplotn)*(  double(iy)/hoplotn) * p3(l);
			  
			  curv.CalcElementTransformation (xl, i-1, grid[ix][iy], dxdxi);
			  
			  gtaux = dxdxi * taux;
			  gtauy = dxdxi * tauy;
			  gridn[ix][iy] = Cross (gtauy, gtaux).Normalize();
			}

		    for (int ix = 0; ix < hoplotn; ix++)
		      for (int iy = 0; iy < hoplotn; iy++)
			{
			  glNormal3dv (gridn[ix][iy]);
			  glVertex3dv (grid[ix][iy]);

			  glNormal3dv (gridn[ix+1][iy]);
			  glVertex3dv (grid[ix+1][iy]);

			  glNormal3dv (gridn[ix+1][iy+1]);
			  glVertex3dv (grid[ix+1][iy+1]);

			  glNormal3dv (gridn[ix][iy+1]);
			  glVertex3dv (grid[ix][iy+1]);
			}
		  }

		glEnd ();
		*/


	      } 
	    else
	      { 



		Point3d c(0,0,0);
		if (vispar.shrink < 1)
		  {
		    for (int j = 1; j <= 5; j++)
		      {
			Point3d p = mesh->Point(el.PNum(j));
			c.X() += p.X() / 5;
			c.Y() += p.Y() / 5;
			c.Z() += p.Z() / 5;
		      }
		  }


		el.GetSurfaceTriangles (faces);
	  
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, pyramidcol);
		if (el.PNum(1))
		  {
		    glBegin (GL_TRIANGLES);
	      
		    for (int j = 1; j <= faces.Size(); j++)
		      {
			Element2d & face = faces.Elem(j);
			Point3d lp1 = mesh->Point (el.PNum(face.PNum(1)));
			Point3d lp2 = mesh->Point (el.PNum(face.PNum(2)));
			Point3d lp3 = mesh->Point (el.PNum(face.PNum(3)));
			Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
			n /= (n.Length()+1e-12);
			n *= -1;
			glNormal3d (n.X(), n.Y(), n.Z());

			if (vispar.shrink < 1)
			  {
			    lp1 = c + vispar.shrink * (lp1 - c);
			    lp2 = c + vispar.shrink * (lp2 - c);
			    lp3 = c + vispar.shrink * (lp3 - c);
			  }

			glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
			glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
			glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
		      }
	      
		    glEnd();
		  }
	      }
	  }
      }
    glEndList ();
  }

  void VisualSceneMesh :: BuildBadelList()
  {
    ;
  }

  void VisualSceneMesh :: BuildIdentifiedList()
  {
    ;
  }
  
  void VisualSceneMesh :: BuildDomainSurfList()
  {
    if (domainsurflist)
      glDeleteLists (domainsurflist, 1);

    domainsurflist = glGenLists (1);
    glNewList (domainsurflist, GL_COMPILE);

    int i, j;
    glLineWidth (1.0f);
  
    glDisable (GL_COLOR_MATERIAL);
  
    for (i = 1; i <= mesh->GetNSE(); i++)
      {
	Element2d el = mesh->SurfaceElement (i);
      
	int drawel = 1;
	for (j = 1; j <= el.GetNP(); j++)
	  {
	    if (!el.PNum(j))
	      drawel = 0;
	  }
      
	if (!drawel)
	  continue;
      
	if (el.GetIndex() < 1 || el.GetIndex() > mesh->GetNFD())
	  continue;
	int domin = mesh->GetFaceDescriptor(el.GetIndex()).DomainIn();
	int domout = mesh->GetFaceDescriptor(el.GetIndex()).DomainOut();

	int fac;
	if (domin == vispar.drawdomainsurf)
	  fac = 1;
	else if (domout == vispar.drawdomainsurf)
	  fac = -1;
	else
	  continue;
      
      
	GLfloat matcol[] = { 1, 0, 0, 1 };
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcol);
      
      
	if (el.GetNP() == 3)
	  {
	    glBegin (GL_TRIANGLES);
	      
	    const Point3d & lp1 = mesh->Point (el.PNum(1));
	    const Point3d & lp2 = mesh->Point (el.PNum(2));
	    const Point3d & lp3 = mesh->Point (el.PNum(3));
	    Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
	    n /= ( fac * (n.Length()+1e-12));
	    glNormal3d (n.X(), n.Y(), n.Z());
	  
	    if (!vispar.colormeshsize)
	      {
		glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	      }
	    glEnd();
	  }
	else if (el.GetNP() == 4)
	  {
	    glBegin (GL_QUADS);
	  
	    const Point3d & lp1 = mesh->Point (el.PNum(1));
	    const Point3d & lp2 = mesh->Point (el.PNum(2));
	    const Point3d & lp3 = mesh->Point (el.PNum(4));
	    const Point3d & lp4 = mesh->Point (el.PNum(3));
	    Vec3d n = Cross (Vec3d (lp1, lp2), 
			     Vec3d (lp1, Center (lp3, lp4)));
	    n /= (fac * (n.Length()+1e-12));
	    glNormal3d (n.X(), n.Y(), n.Z()); 
	    glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	    glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	    glVertex3d (lp4.X(), lp4.Y(), lp4.Z());
	    glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	    glEnd();
	  }
	else if (el.GetNP() == 6)
	  {
	    glBegin (GL_TRIANGLES);
	    static int trigs[4][3] = {
	      { 1, 6, 5 },
	      { 2, 4, 6 },
	      { 3, 5, 4 },
	      { 4, 5, 6 } };
	  
	    for (j = 0; j < 4; j++)
	      {
		const Point3d & lp1 = mesh->Point (el.PNum(trigs[j][0]));
		const Point3d & lp2 = mesh->Point (el.PNum(trigs[j][1]));
		const Point3d & lp3 = mesh->Point (el.PNum(trigs[j][2]));
		Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
		n /= (fac * (n.Length() + 1e-12));
		glNormal3d (n.X(), n.Y(), n.Z());
		glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
		glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
		glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	      }
	    glEnd();
	  }
      }
    glEndList ();
  }












  void VisualSceneMesh :: MouseDblClick (int px, int py)
  {
    int i, hits;

    // select surface triangle by mouse click

    GLuint selbuf[10000];
    glSelectBuffer (10000, selbuf);


    glRenderMode (GL_SELECT);

    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);


    glMatrixMode (GL_PROJECTION); 
    glPushMatrix();

    GLdouble projmat[16];
    glGetDoublev (GL_PROJECTION_MATRIX, projmat);

    glLoadIdentity(); 
    gluPickMatrix (px, viewport[3] - py, 1, 1, viewport); 
    glMultMatrixd (projmat);
  


    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode (GL_MODELVIEW); 

    glPushMatrix();
    glMultMatrixf (transformationmat);


    //  SetClippingPlane();

    glInitNames();
    glPushName (1);

    glPolygonOffset (1, 1);
    glEnable (GL_POLYGON_OFFSET_FILL);

    glDisable(GL_CLIP_PLANE0);
  
    if (vispar.clipenable)
      {
	Vec<3> n(clipplane[0], clipplane[1], clipplane[2]);
	double len = Abs(n);
	double mu = -clipplane[3] / (len*len);
	Point<3> p (mu * n);
	n /= len;
	Vec<3> t1 = n.GetNormal ();
	Vec<3> t2 = Cross (n, t1);
      
	double xi1mid = (center - p) * t1;
	double xi2mid = (center - p) * t2;
      
	glLoadName (0);  
	glBegin (GL_QUADS);
	glVertex3dv (p + (xi1mid-rad) * t1 + (xi2mid-rad) * t2);
	glVertex3dv (p + (xi1mid+rad) * t1 + (xi2mid-rad) * t2);
	glVertex3dv (p + (xi1mid+rad) * t1 + (xi2mid+rad) * t2);
	glVertex3dv (p + (xi1mid-rad) * t1 + (xi2mid+rad) * t2);
	glEnd ();
      }

    //  SetClippingPlane();

    glCallList (filledlist);

    glDisable (GL_POLYGON_OFFSET_FILL);
  
    glPopName();

    glMatrixMode (GL_PROJECTION); 
    glPopMatrix();

    glMatrixMode (GL_MODELVIEW); 
    glPopMatrix();

    glFlush();  

	
    hits = glRenderMode (GL_RENDER);

    //  cout << "hits = " << hits << endl;

    int minname = 0;
    GLuint mindepth = 0;

    // find clippingplane
    GLuint clipdepth = 0; // GLuint(-1);

    for (i = 0; i < hits; i++)
      {
	int curname = selbuf[4*i+3];
	if (!curname) clipdepth = selbuf[4*i+1];
      }

    for (i = 0; i < hits; i++)
      {
	int curname = selbuf[4*i+3];
	GLuint curdepth = selbuf[4*i+1];
	/*
	  cout << selbuf[4*i] << " " << selbuf[4*i+1] << " " 
	  << selbuf[4*i+2] << " " << selbuf[4*i+3] << endl;
	*/
	if (curname && (curdepth > clipdepth) &&
	    (curdepth < mindepth || !minname))
	  {
	    mindepth = curdepth;
	    minname = curname;
	  }
      }

    seledge = -1;
    if (minname)
      {
	const Element2d & sel = mesh->SurfaceElement(minname);


	cout << "select element " << minname
	     << " on face " << sel.GetIndex() << endl;
	cout << "Nodes: ";
	for (i = 1; i <= sel.GetNP(); i++)
	  cout << sel.PNum(i) << " ";
	cout << endl;

	selelement = minname;
	selface = mesh->SurfaceElement(minname).GetIndex();

	locpi = (locpi % sel.GetNP()) + 1;
	selpoint2 = selpoint;
	selpoint = sel.PNum(locpi);
	cout << "selected point " << selpoint 
	     << ", pos = " << mesh->Point (selpoint) 
	     << endl;

	for (i = 1; i <= mesh->GetNSeg(); i++)
	  {
	    const Segment & seg = mesh->LineSegment(i);
	    if (seg.p1 == selpoint && seg.p2 == selpoint2 || 
		seg.p2 == selpoint && seg.p1 == selpoint2)
	      {
		seledge = seg.edgenr;
		cout << "seledge = " << seledge << endl;
	      }
	  }
      
      }
    else
      {
	selface = -1;
	selelement = -1;
	selpoint = -1;
	selpoint2 = -1;
      }

    glDisable(GL_CLIP_PLANE0);

    selecttimestamp = NextTimeStamp();
  }





}









  
