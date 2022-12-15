Mesh.MshFileVersion=1;

lc = 0.5;

dx = 1.0;
dy = 1.0;
dz = 1.0;

nx = 20;
ny = 20;
nz = 20;

cx = dx/nx;
cy = dy/ny;
cz = (cx + cy)/2;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {dx,0.0,0.0,lc};
Point(3) = {dx,dy,0.0,lc};
Point(4) = {0.0,dy,0.0,lc};
Point(5) = {0.0,dy*0.25,0.0,lc};
Point(6) = {dx,dy*0.25,0.0,lc};

Line(1) = {4,5};
Line(2) = {5,6};
Line(3) = {6,3};
Line(4) = {3,4};
Line(5) = {5,1};
Line(6) = {1,2};
Line(7) = {2,6};
Line Loop(8) = {1,2,3,4};
Plane Surface(9) = {8};
Line Loop(10) = {2,-7,-6,-5};
Plane Surface(11) = {10};

Transfinite Line {3,1} = ny*0.75+1 Using Progression 1.0;
Transfinite Line {6,2,4} = nx+1 Using Progression 1.0;
Transfinite Line {7,5} = ny*0.25+1 Using Progression 1.0;

Transfinite Surface {9} = {6,3,4,5};
Transfinite Surface {11} = {2,6,5,1};
Recombine Surface {9,11};


Extrude Surface {9, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {11, {0,0,dz} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

// Top surface
Physical Surface(33) = {32};

// Bottom surface
Physical Surface(34) = {50};

// Left surface
Physical Surface(35) = {20,54};

// Right surface
Physical Surface(36) = {46,28};

// Front surface
Physical Surface(37) = {33,55};

// Back surface
Physical Surface(38) = {9,11};

// Box volume
Physical Volume (39) = {9000};
Physical Volume (40) = {9001};

