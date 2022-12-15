Mesh.MshFileVersion=1;
dx = 1.0;
dy = 1.0;
dz = 5.0;

nx = 20;
ny = 20;
nz = 1;

cx = dx/nx;
cy = dy/ny;
cz = (cx + cy)/2;

Point(1) = {0.0,0.0,0.0,1.0};
Point(2) = {dx,0.0,0.0,1.0};
Point(3) = {dx,dy,0.0,1.0};
Point(4) = {0.0,dy,0.0,1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude Surface {6, {0,0,dz} } {
  Layers { {1}, {32}, {1.0} };
  Recombine;
};

Transfinite Line {4,10,2,8} = nx + 1 Using Progression 1.0;
Transfinite Line {3,9,1,11} = ny + 1 Using Progression 1.0;
Transfinite Line {14,18,13,22} = nz + 1 Using Progression 1.0;

Transfinite Surface {6} = {3,2,1,4};
Transfinite Surface {27} = {5,14,2,3};
Transfinite Surface {15} = {5,3,4,6};
Transfinite Surface {28} = {6,10,14,5};
Transfinite Surface {23} = {14,2,1,10};
Transfinite Surface {19} = {6,10,1,4};

Recombine Surface {27,23,6,19,15,28};

// Top surface
Physical Surface(33) = {23};
// Bottom surface
Physical Surface(34) = {15};
// Left surface
Physical Surface(35) = {27};
// Right surface
Physical Surface(36) = {19};
// Front surface
Physical Surface(37) = {28};
// Back surface
Physical Surface(38) = {6};
// Box volume
Physical Volume (39) = {32};
