Mesh.MshFileVersion=1;

// Dimensions in mm

n = 4;

l = 2;

dz = 0.0025;

Point(1) = {0.000000,-0.0150000,0.000000,0.1};
Point(2) = {0.000000,0.0150000,0.000000,0.1};
Point(3) = {0.02500000,0.0150000,0.000000,0.1};
Point(4) = {0.02500000,-0.0150000,0.000000,0.1};

Point(5) = {0.025000000,-0.0050000,0.000000,0.1};
Point(6) = {0.025000000,+0.0050000,0.000000,0.1};
Point(7) = {0.0000000,-0.0050000,0.000000,0.1};
Point(8) = {0.0000000,+0.0050000,0.000000,0.1};

Line(1) = {5,4};
Line(2) = {4,1};
Line(3) = {1,7};
Line(4) = {7,5};
Line(5) = {5,6};
Line(6) = {6,8};
Line(7) = {8,7};
Line(8) = {6,3};
Line(9) = {3,2};
Line(10) = {2,8};

Line Loop(11) = {1,2,3,4};
Plane Surface(12) = {11};

Line Loop(13) = {6,7,4,5};
Plane Surface(14) = {13};

Line Loop(15) = {8,9,10,-6};
Plane Surface(16) = {15};

Transfinite Line {2,4,6,9} = 5*n+1 Using Progression 1.0;
Transfinite Line {3,1,8,10} = 2*n + 1 Using Progression 1.0;
Transfinite Line {7,5} = 2*n + 1 Using Progression 1.0;

Transfinite Surface {16} = {3,2,8,6};
Transfinite Surface {14} = {6,8,7,5};
Transfinite Surface {12} = {5,7,1,4};

Recombine Surface {12,14,16};

Extrude Surface {12, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {14, {0,0,dz} } {
  Layers { {l}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {16, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {38, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {60, {0,0,dz} } {
  Layers { {l}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {82, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {117, {-0.005,0,0} } {
  Layers { {2*l}, {9002}, {1.0} };
  Recombine;
};

Extrude Surface {170, {-0.01,0,0} } {
  Layers { {4*l}, {9002}, {1.0} };
  Recombine;
};

Extrude Surface {179, {0,0,-0.0025} } {
  Layers { {l}, {9002}, {1.0} };
  Recombine;
};

// Inlet
Physical Surface(301) = {192,209};

// Outlet
Physical Surface(302) = {59,125};

// Walls
Physical Surface(303) = {205,183,161,213,191,169,47,113,37,103,214,201,157,51,14,126,165,187};

Physical Volume(400) = {9001};
Physical Volume(401) = {9002};
