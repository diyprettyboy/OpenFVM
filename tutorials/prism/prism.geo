Mesh.MshFileVersion=1;
nT=6;
nC=15;
nt=8;
T=1;
lc=0.2;
lc2=0.01;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {0.4,0.0,0.0,lc};
Point(3) = {0.4,0.1,0.0,lc2};
Point(4) = {0.0,0.1,0.0,lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Transfinite Line {3,1} = 40 Using Progression 1.0;
Transfinite Line {2,4} = 10 Using Progression 1.0;
Transfinite Surface {6} = {3,4,1,2};

Extrude Surface {6, {0.0,0.0,0.01}}
  { Recombine; Layers { {nT}, {9000}, {T} }; 
};

// Inlet
Physical Surface(73) = {6};
// Outlet
Physical Surface(74) = {28};
// Walls
Physical Surface(75) = {15,27,23,19};
Physical Volume(76) = {9000};


