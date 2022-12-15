Mesh.MshFileVersion=1;

lc = 0.05;

dx = 1.0;
dy = 1.0;
dz = 5.0;

nx = 25;
ny = 25;
nz = 1;

cx = dx/nx;
cy = dy/ny;
cz = (cx + cy)/2;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {dx,0.0,0.0,lc};
Point(3) = {dx,dy,0.0,lc};
Point(4) = {0.0,dy,0.0,lc};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,1};
Line(4) = {1,2};

Line Loop(5) = {2,3,4,1};

Plane Surface(6) = {5};

Extrude Surface {6, {0,0,dz} } {
  Layers { {nz}, {32}, {1.0} };
  //Recombine;
};

/*
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
*/

// Box volume
Physical Volume (39) = {32};

// Top surface
Physical Surface(33) = {15};

// Bottom surface
Physical Surface(34) = {23};

// Left surface
Physical Surface(35) = {19};

// Right surface
Physical Surface(36) = {27};

// Front surface
Physical Surface(37) = {28};

// Back surface
Physical Surface(38) = {6};
