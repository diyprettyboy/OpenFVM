/// Function Cylindric: 
// create a cylindre - of diameter dCyl and height hCyl
//                   - of center xCyl, yCyl, zCyl
//                   - characteristic length lcCyl
//Function cylindric

  p1 = newp; Point(p1) = {xCyl       ,    yCyl    , zCyl,  lcCyl} ;
  p2 = newp; Point(p2) = {xCyl+dCyl/2,    yCyl    , zCyl,  lcCyl} ;
  p3 = newp; Point(p3) = {xCyl       , yCyl+dCyl/2, zCyl,  lcCyl} ;
 
  p4 = newp; Point(p4) = {xCyl       ,    yCyl    , zCyl+hCyl,  lcCyl} ;
  p5 = newp; Point(p5) = {xCyl+dCyl/2,    yCyl    , zCyl+hCyl,  lcCyl} ;
  p6 = newp; Point(p6) = {xCyl       , yCyl+dCyl/2, zCyl+hCyl,  lcCyl} ;

  Cyles_Point1[iCyl] = p1 ;
  Cyles_Point2[iCyl] = p2 ;
  Cyles_Point3[iCyl] = p3 ;
  Cyles_Point4[iCyl] = p4 ;
  Cyles_Point5[iCyl] = p5 ;
  Cyles_Point6[iCyl] = p6 ;

  c1 = newreg; Line(c1) = {p1,p4};
  c2 = newreg; Line(c2) = {p2,p5};
  c3 = newreg; Line(c3) = {p3,p6};
  c4 = newreg; Line(c4) = {p1,p2};
  c5 = newreg; Line(c5) = {p1,p3};
  c6 = newreg; Line(c6) = {p4,p5};
  c7 = newreg; Line(c7) = {p4,p6};

  c8 = newreg; Circle(c8) = {p2,p1,p3};
  c9 = newreg; Circle(c9) = {p5,p4,p6};
  
  Cyles_Linep1p2[iCyl] = c4 ;
  Cyles_Linep1p3[iCyl] = c5 ;
  Cyles_Linep4p5[iCyl] = c6 ;
  Cyles_Linep5p2[iCyl] = -c2 ;
  Cyles_Linep4p6[iCyl] = c7 ;
  Cyles_Linep6p3[iCyl] = -c3 ;
  Cyles_Curvep3p2[iCyl] = -c8 ; 

  l1 = newreg; Line Loop(l1) = {c4,c8,-c5}; Plane Surface(l1+1) = {l1};
  l2 = newreg; Line Loop(l2) = {c6,c9,-c7}; Plane Surface(l2+1) = {l2};
  l3 = newreg; Line Loop(l3) = {c4,c2,-c6,-c1}; Plane Surface(l3+1) = {l3};
  l4 = newreg; Line Loop(l4) = {c5,c3,-c7,-c1}; Plane Surface(l4+1) = {l4};
  l5 = newreg; Line Loop(l5) = {c8,c3,-c9,-c2}; Ruled Surface(l5+1) = {l5};

  Cyles_BaseSurf[iCyl] = l1+1 ;
  Cyles_HautSurf[iCyl] = l2+1 ;
  Cyles_RuleSurf1[iCyl] = l5+1 ;
  Cyles_Surf3[iCyl] = l3+1 ;
  Cyles_Surf4[iCyl] = l4+1 ;
  
  s = newreg; Surface Loop(s) = {l1+1,l2+1,l3+1,l4+1,l5+1}; Volume(s+1) = s ;
  
  Cyles_SurfaceLoop[iCyl] = s ;
  Cyles_Volume[iCyl] = s+1 ; 

//Return


