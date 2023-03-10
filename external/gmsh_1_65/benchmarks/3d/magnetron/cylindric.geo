/// Function Cylindric: 
// create a cylindre - of diameter dCyl and height hCyl
//                   - of center xCyl, yCyl, zCyl
//                   - characteristic length lcCyl
//Function cylindric

  p1 = newp; Point(p1) = {xCyl       ,    yCyl    , zCyl,  lcCyl} ;
  p2 = newp; Point(p2) = {xCyl+dCyl/2,    yCyl    , zCyl,  lcCyl} ;
  p3 = newp; Point(p3) = {xCyl       , yCyl+dCyl/2, zCyl,  lcCyl} ;
  p4 = newp; Point(p4) = {xCyl-dCyl/2,    yCyl    , zCyl,  lcCyl} ;
  p5 = newp; Point(p5) = {xCyl       , yCyl-dCyl/2, zCyl,  lcCyl} ;

  p6 = newp; Point(p6) = {xCyl       ,    yCyl    , zCyl+hCyl,  lcCyl} ;
  p7 = newp; Point(p7) = {xCyl+dCyl/2,    yCyl    , zCyl+hCyl,  lcCyl} ;
  p8 = newp; Point(p8) = {xCyl       , yCyl+dCyl/2, zCyl+hCyl,  lcCyl} ;
  p9 = newp; Point(p9) = {xCyl-dCyl/2,    yCyl    , zCyl+hCyl,  lcCyl} ;
  p10 = newp; Point(p10) = {xCyl       ,yCyl-dCyl/2, zCyl+hCyl,  lcCyl} ;

  c1 = newreg; Line(c1) = {p2,p7};
  c2 = newreg; Line(c2) = {p3,p8};
  c3 = newreg; Line(c3) = {p4,p9};
  c4 = newreg; Line(c4) = {p5,p10};

  c5 = newreg; Circle(c5) = {p2,p1,p3};
  c6 = newreg; Circle(c6) = {p3,p1,p4};
  c7 = newreg; Circle(c7) = {p4,p1,p5};
  c8 = newreg; Circle(c8) = {p5,p1,p2};

  c9 = newreg; Circle(c9) = {p7,p6,p8};
  c10 = newreg; Circle(c10) = {p8,p6,p9};
  c11 = newreg; Circle(c11) = {p9,p6,p10};
  c12 = newreg; Circle(c12) = {p10,p6,p7};

  l1 = newreg; Line Loop(l1) = {c1,c9,-c2,-c5}; Ruled Surface(l1+1) = {l1};
  l2 = newreg; Line Loop(l2) = {c2,c10,-c3,-c6};Ruled Surface(l2+1) = {l2};
  l3 = newreg; Line Loop(l3) = {c3,c11,-c4,-c7};Ruled Surface(l3+1) = {l3};
  l4 = newreg; Line Loop(l4) = {c4,c12,-c1,-c8};Ruled Surface(l4+1) = {l4};
  l5 = newreg; Line Loop(l5) = {c5,c6,c7,c8}; Plane Surface(l5+1) = {l5};
  l6 = newreg; Line Loop(l6) = {c9,c10,c11,c12}; Plane Surface(l6+1) = {l6};

  Cyles_BaseLoop[iCyl] = l5 ;
  Cyles_BaseSurf[iCyl] = l5+1;
  Cyles_HautSurf[iCyl] = l6+1;
  Cyles_RuleSurf1[iCyl] = l1+1 ;
  Cyles_RuleSurf2[iCyl] = l2+1 ;
  Cyles_RuleSurf3[iCyl] = l3+1 ;
  Cyles_RuleSurf4[iCyl] = l4+1 ;


  s = newreg; Surface Loop(s) = {l2+1,l4+1,l5+1,l6+1,l3+1,l1+1}; Volume(s+1) = s ;
/*
  Cyles_SurfaceLoop[iCyl] = s ;
  Cyles_Volume[iCyl] = s+1 ; 
*/
//Return


