Mesh.MshFileVersion=1;

lc = 0.5;

dx = 0.05;
dy = 0.05;
dz = 1.0;

nz = 80;

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

Recombine Surface {6};

Extrude Surface {6, {0,0,dz/3} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {28, {0,0,dz/3} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {50, {0,0,dz/3} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

/*
Extrude Surface {72, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {94, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {116, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {138, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {160, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {182, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {204, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {226, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {248, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {270, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {292, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {314, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {336, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {358, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {380, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {402, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {424, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};
*/

// Inlet
Physical Surface(30) = {6};

// Outlet
Physical Surface(32) = {72};

// Box volume
Physical Volume (35) = {9000};
Physical Volume (36) = {9001};


