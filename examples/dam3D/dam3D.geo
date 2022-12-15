Mesh.MshFileVersion=1;
a = 0.148;
d = 0.024;

n = 100;
c = 4*a/n;

nz = 2*a/c;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {4*a, 0, 0, 0.1};
Point(3) = {4*a, 3*a, 0, 0.1};
Point(4) = {0, 3*a, 0, 0.1};

Point(5) = {2*a, 0, 0, 0.1};
Point(6) = {2*a+d, 0, 0, 0.1};

Point(7) = {2*a, 2*d, 0, 0.1};
Point(8) = {2*a+d, 2*d, 0, 0.1};

Point(9) = {0, 2*d, 0, 0.1};
Point(10) = {4*a, 2*d, 0, 0.1};

Point(11) = {2*a, 3*a, 0, 0.1};
Point(12) = {2*a+d, 3*a, 0, 0.1};

Point(13) = {a, 0, 0, 0.1};
Point(14) = {a, 2*d, 0, 0.1};
Point(15) = {a, 3*a, 0, 0.1};

Point(16) = {a, 2*a, 0, 0.1};
Point(17) = {0, 2*a, 0, 0.1};
Point(18) = {2*a, 2*a, 0, 0.1};
Point(19) = {2*a+d, 2*a, 0, 0.1};
Point(20) = {4*a, 2*a, 0, 0.1};

Line(1) = {4,17};
Line(2) = {17,16};
Line(3) = {16,15};
Line(4) = {4,15};
Line(5) = {15,11};
Line(6) = {11,18};
Line(7) = {18,16};
Line(8) = {11,12};
Line(9) = {19,12};
Line(10) = {18,19};
Line(11) = {19,8};
Line(12) = {8,7};
Line(13) = {7,18};
Line(14) = {7,14};
Line(15) = {14,16};
Line(16) = {14,9};
Line(17) = {9,17};
Line(18) = {9,1};
Line(19) = {1,13};
Line(20) = {13,14};
Line(21) = {5,7};
Line(22) = {13,5};
Line(23) = {8,10};
Line(24) = {10,2};
Line(25) = {2,6};
Line(26) = {6,8};
Line(27) = {20,10};
Line(28) = {19,20};
Line(29) = {3,20};
Line(30) = {12,3};

Line Loop(31) = {1,2,3,-4};
Plane Surface(32) = {31};
Line Loop(33) = {7,3,5,6};
Plane Surface(34) = {33};
Line Loop(35) = {9,-8,6,10};
Plane Surface(36) = {35};
Line Loop(37) = {28,-29,-30,-9};
Plane Surface(38) = {37};
Line Loop(39) = {11,23,-27,-28};
Plane Surface(40) = {39};
Line Loop(41) = {25,26,23,24};
Plane Surface(42) = {41};
Line Loop(43) = {11,12,13,10};
Plane Surface(44) = {43};
Line Loop(45) = {7,-15,-14,13};
Plane Surface(46) = {45};
Line Loop(47) = {2,-15,16,17};
Plane Surface(48) = {47};
Line Loop(49) = {16,18,19,20};
Plane Surface(50) = {49};
Line Loop(51) = {22,21,14,-20};
Plane Surface(52) = {51};

Transfinite Line {1,3,6,9,29} = a/c+1 Using Progression 1.0;
Transfinite Line {27,11,13,15,17} = (2*a-2*d)/c+1 Using Progression 1.0;
Transfinite Line {25,23,28,30} = (2*a-d)/c+1 Using Progression 1.0;
Transfinite Line {12,10,8} = d/c+1 Using Progression 1.0;
Transfinite Line {22,14,7,5} = a/c+1 Using Progression 1.0;
Transfinite Line {4,2,16,19} = a/c+1 Using Progression 1.0;
Transfinite Line {18,20,21,26,24} = 2*d/c+1 Using Progression 1.0;

Transfinite Surface {32} = {4,15,16,17};
Transfinite Surface {34} = {11,18,16,15};
Transfinite Surface {36} = {11,12,19,18};
Transfinite Surface {38} = {12,3,20,19};
Transfinite Surface {40} = {19,20,10,8};
Transfinite Surface {44} = {18,19,8,7};
Transfinite Surface {46} = {16,18,7,14};
Transfinite Surface {48} = {17,16,14,9};
Transfinite Surface {50} = {9,14,13,1};
Transfinite Surface {52} = {14,7,5,13};
Transfinite Surface {42} = {8,10,2,6};
Recombine Surface {32,34,36,38,40,44,46,48,50,52,42};

Extrude Surface {32, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {34, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {36, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {38, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {40, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {42, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {44, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {46, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {48, {0,0,2*a} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {50, {0,0,2*a} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {52, {0,0,2*a} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

// Top
Physical Surface(137) = {73,91,109,135};

// Bottom 
Physical Surface(138) = {267,281,285,197,175,171};

// Right 
Physical Surface(139) = {183,157,131};

// Left 
Physical Surface(140) = {263,249,61};

// Walls
Physical Surface(141) = {74,96,118,140,250,228,206,162,272,294,184,38,40,44,36,34,32,46,48,153,223,50,42,52};

Physical Volume(150) = {9000};
Physical Volume(151) = {9001};
