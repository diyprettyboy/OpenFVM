Mesh.MshFileVersion=1;
cl = 0.5;

n = 81;

dx = 1.0;
dy = dx;
dz = dx / n;

Point(1) = {0.0,0.0,0.0,cl};
Point(2) = {dx,0.0,0.0,cl};
Point(3) = {dx,dy,0.0,cl};
Point(4) = {0.0,dy,0.0,cl};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Transfinite Line {4,3,2,1} = n Using Progression 1;
Transfinite Surface {6} = {4,3,2,1};
Recombine Surface {6};

Extrude Surface {6, {0,0,dz} } {
  Layers { {1}, {1}, {1.0} };
  Recombine;
};

Physical Volume(1) = {1};
Physical Volume(16) = {16};

/*
cl = 0.5;

dx = 1.0;
dy = 1.0;
dz = 0.05;

n1 = 20;
n2 = 10;

nx = 41;
ny = 41;
nz = 1;

cx = dx / nx;
cy = dy / ny;
cz = dz / nz;

Point(1) = {0.0,0.0,0.0,cl};
Point(2) = {dx,0.0,0.0,cl};
Point(3) = {dx,dy,0.0,cl};
Point(4) = {0.0,dy,0.0,cl};

Point(5) = {dx*0.5,dy*0.5,0.0,cl};

Point(6) = {dx*0.4,dy*0.4,0.0,cl};
Point(7) = {dx*0.4,dy*0.6,0.0,cl};
Point(8) = {dx*0.6,dy*0.6,0.0,cl};
Point(9) = {dx*0.6,dy*0.4,0.0,cl};

Point(10) = {0.0,dy*0.4,0.0,cl};
Point(11) = {dx,dy*0.4,0.0,cl};

Point(12) = {0.0,dy*0.6,0.0,cl};
Point(13) = {dx,dy*0.6,0.0,cl};

Point(14) = {dx*0.4,0,0.0,cl};
Point(15) = {dx*0.4,dy,0.0,cl};

Point(16) = {dx*0.6,0,0.0,cl};
Point(17) = {dx*0.6,dy,0.0,cl};

Line(1) = {4,15};
Line(2) = {15,7};
Line(3) = {7,12};
Line(4) = {12,4};
Line(5) = {15,17};
Line(6) = {17,8};
Line(7) = {8,7};
Line(8) = {17,3};
Line(9) = {3,13};
Line(10) = {13,8};
Line(11) = {8,9};
Line(12) = {9,11};
Line(13) = {11,13};
Line(14) = {9,16};
Line(15) = {16,2};
Line(16) = {2,11};
Line(17) = {9,6};
Line(18) = {6,14};
Line(19) = {14,16};
Line(20) = {6,7};
Line(21) = {12,10};
Line(22) = {10,6};
Line(23) = {10,1};
Line(24) = {1,14};
Line Loop(25) = {8,9,10,-6};
Plane Surface(26) = {25};
Line Loop(27) = {13,10,11,12};
Plane Surface(28) = {27};
Line Loop(29) = {7,-2,5,6};
Plane Surface(30) = {29};
Line Loop(31) = {1,2,3,4};
Plane Surface(32) = {31};
Line Loop(33) = {11,17,20,-7};
Plane Surface(34) = {33};
Line Loop(35) = {14,15,16,-12};
Plane Surface(36) = {35};
Line Loop(37) = {17,18,19,-14};
Plane Surface(38) = {37};
Line Loop(39) = {24,-18,-22,23};
Plane Surface(40) = {39};
Line Loop(41) = {21,22,20,3};
Plane Surface(42) = {41};

Transfinite Line {6,9,2,4,23,18,14,16,15,15,12,10,8,1,3,22,24} = n1 + 1 Using Progression 1.0;
Transfinite Line {13,11,7,20,17,5,19,21} = n2 + 1 Using Progression 1.0;

Transfinite Surface {26} = {17,3,13,8};
Transfinite Surface {28} = {8,13,11,9};
Transfinite Surface {36} = {9,11,2,16};
Transfinite Surface {38} = {9,16,14,6};
Transfinite Surface {34} = {8,9,6,7};
Transfinite Surface {30} = {17,8,7,15};
Transfinite Surface {32} = {15,7,12,4};
Transfinite Surface {42} = {12,7,6,10};
Transfinite Surface {40} = {10,6,14,1};
Recombine Surface {26,30,32,42,34,28,36,38,40};

Extrude Surface {26, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {28, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {30, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {32, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {34, {0,0,dz} } {
  Layers { {nz}, {16}, {1.0} };
  Recombine;
};

Extrude Surface {36, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {38, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {40, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {42, {0,0,dz} } {
  Layers { {nz}, {1}, {1.0} };
  Recombine;
};

// Top
Physical Surface(241) = {51,103,117};
Physical Surface(242) = {130,32,108,30,64,26,55};
Physical Surface(243) = {64,26,55,108,30,130,32,129,227,240,42,152,34,86,28,73,169,174,36,196,38,218,40,217,205,191,165};

Physical Volume(1) = {1};
Physical Volume(16) = {16};
*/
