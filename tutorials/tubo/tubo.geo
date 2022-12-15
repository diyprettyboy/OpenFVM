Mesh.MshFileVersion=1;

dz=1.0;
nT=50;
lc=0.015;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {0.0,0.05,0.0,lc};
Point(3) = {0.05,0.0,0.0,lc};
Point(4) = {-0.05,0,0.0,lc};
Point(5) = {0,-0.05,0.0,lc};

Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,4};
Circle(4) = {4,1,2};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude Surface {6, {0.0,0.0,dz}}
  { Recombine; Layers { {nT}, {9000}, {1.0} }; };

// Left
Physical Surface(73) = {28};
// Right
Physical Surface(74) = {6};
// Walls
Physical Surface(75) = {15,27,23,19};

Physical Volume(76) = {9000};

