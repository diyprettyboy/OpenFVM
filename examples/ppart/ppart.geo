Mesh.MshFileVersion=1;

lc = 0.005;

n = 1;

l = 1;

// Metric conversion (m: c = 1.0)
c = 1;

dz = 0.005*c;

Point(2) = {-0.100000*c,-0.0550000*c,0.000000,lc};
Point(45) = {-0.020000*c,-0.0550000*c,0.000000,lc};
Point(49) = {0.020000*c,-0.0550000*c,0.000000,lc};
Point(4) = {0.100000*c,-0.0550000*c,0.000000,lc};
Point(3) = {0.100000*c,0.0550000*c,0.000000,lc};
Point(19) = {0.020000*c,0.0550000*c,0.000000,lc};
Point(23) = {-0.020000*c,0.0550000*c,0.000000,lc};
Point(1) = {-0.100000*c,0.0550000*c,0.000000,lc};
Point(79) = {-0.020000*c,-0.030000*c,0.000000,lc};
Point(134) = {0.020000*c,-0.030000*c,0.000000,lc};
Point(137) = {0.020000*c,0.030000*c,0.000000,lc};
Point(76) = {-0.020000*c,0.030000*c,0.000000,lc};
Point(31) = {0.100000*c,-0.030000*c,0.000000,lc};
Point(37) = {0.100000*c,0.030000*c,0.000000,lc};
Point(5) = {-0.100000*c,0.030000*c,0.000000,lc};
Point(11) = {-0.100000*c,-0.030000*c,0.000000,lc};

Point(100) = {-0.100000*c,-0.005000*c,0.000000,lc};
Point(101) = {-0.100000*c,+0.0050000*c,0.000000,lc};
Point(102) = {-0.0200000*c,-0.0050000*c,0.000000,lc};
Point(103) = {-0.0200000*c,+0.0050000*c,0.000000,lc};

Point(104) = {-0.100000*c-0.005*c,-0.005000*c,0.000000,lc};
Point(105) = {-0.100000*c-0.005*c,+0.0050000*c,0.000000,lc};

Line(1) = {1,23};
Line(2) = {23,76};
Line(3) = {76,5};
Line(4) = {5,1};
Line(5) = {5,101};
Line(6) = {101,103};
Line(7) = {103,76};
Line(8) = {76,137};
Line(9) = {137,19};
Line(10) = {19,23};
Line(11) = {103,102};
Line(12) = {102,100};
Line(13) = {100,101};
Line(14) = {100,11};
Line(15) = {11,79};
Line(16) = {79,102};
Line(17) = {79,45};
Line(18) = {45,2};
Line(19) = {11,2};
Line(20) = {79,134};
Line(21) = {134,49};
Line(22) = {49,45};
Line(23) = {134,137};
Line(24) = {137,37};
Line(25) = {37,31};
Line(26) = {31,134};
Line(27) = {31,4};
Line(28) = {4,49};
Line(29) = {19,3};
Line(30) = {3,37};


Line Loop(31) = {1,2,3,4};
Plane Surface(32) = {31};
Line Loop(33) = {3,5,6,7};
Plane Surface(34) = {-33};
Line Loop(35) = {6,11,12,13};
Plane Surface(36) = {35};
Line Loop(37) = {14,15,16,12};
Plane Surface(38) = {-37};
Line Loop(39) = {15,17,18,-19};
Plane Surface(40) = {39};

Line Loop(41) = {20,21,22,-17};
Plane Surface(42) = {41};
Line Loop(43) = {8,9,10,2};
Plane Surface(44) = {-43};
Line Loop(45) = {23,24,25,26};
Plane Surface(46) = {45};
Line Loop(47) = {27,28,-21,-26};
Plane Surface(48) = {47};
Line Loop(49) = {29,30,-24,9};
Plane Surface(50) = {49};

/*
Transfinite Line {1,3,6,12,15,18,26,28,24,29} = 16*n + 1 Using Progression 1.0;
Transfinite Line {13,11} = 2*n + 1 Using Progression 1.0;
Transfinite Line {30,9,2,4,27,21,17,19} = 5*n + 1 Using Progression 1.0;
Transfinite Line {7,5,14,16} = 5*n + 1 Using Progression 1.0;
Transfinite Line {25,23} = 12*n + 1 Using Progression 1.0;
Transfinite Line {8,10,20,22} = 7*n + 1 Using Progression 1.0;

Transfinite Surface {50} = {3,37,137,19};
Transfinite Surface {46} = {137,37,31,134};
Transfinite Surface {48} = {31,4,49,134};
Transfinite Surface {42} = {79,134,49,45};
Transfinite Surface {44} = {19,137,76,23};
Transfinite Surface {32} = {23,76,5,1};
Transfinite Surface {34} = {5,76,103,101};
Transfinite Surface {36} = {103,102,100,101};
Transfinite Surface {38} = {102,79,11,100};
Transfinite Surface {40} = {11,79,45,2};

Recombine Surface {44,50,46,48,42,36,34,32,38,40};
*/

Extrude Surface {32, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {34, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {36, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {38, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {40, {0,0,dz} } {
  Layers { {l}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {42, {0,0,dz} } {
  Layers { {l}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {44, {0,0,dz} } {
  Layers { {l}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {46, {0,0,dz} } {
  Layers { {l}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {48, {0,0,dz} } {
  Layers { {l}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {50, {0,0,dz} } {
  Layers { {l}, {9001}, {1.0} };
  Recombine;
};

Line(402) = {101,105};
Line(403) = {105,104};
Line(404) = {104,100};

Line Loop(405) = {402,403,404,13};
Plane Surface(406) = {405};

Extrude Surface {406, {0,0,dz} } {
  Layers { {n}, {9002}, {1.0} };
  Recombine;
};

// Inlet
Physical Surface(301) = {419};

// Walls
Physical Surface(302) = {235,261,221};
Physical Surface(303) = {283,291,89,137,159,155,71,59,195,257,239,177,169,129,107,81,203,213};
Physical Surface(304) = {34,94,32,72,36,116,138,38,160,40,182,42,44,204,50,270,46,226,248,248,48,279,287};
Physical Surface(305) = {415,423,428,406};

Physical Volume(399) = {9000};
Physical Volume(400) = {9001};
Physical Volume(401) = {9002};
