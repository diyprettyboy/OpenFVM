Mesh.MshFileVersion=1;

lc = 0.5;

dx = 0.05;
dy = 0.05;
dz = 1.0;

nz1 = 100;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {dx,0.0,0.0,lc};
Point(3) = {dx,dy,0.0,lc};
Point(4) = {0.0,dy,0.0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {1,2,3,4};

Plane Surface(6) = {5};

Extrude Surface {6, {0,0,dz} } {
  Layers { {nz1}, {9000}, {1.0} };
  Recombine;
};

Recombine Surface {6};

// Hot surface
Physical Surface(30) = {6};

// Cold surface
Physical Surface(31) = {28};

// Side walls
Physical Surface(32) = {15,19,23,27};

// Box volume
Physical Volume (35) = {9000};

