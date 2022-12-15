Mesh.MshFileVersion=1;

n = 50;
c = 1.0/n;

nz = 25;
dz = 0.5;

h1  = c*n*0.1;

h2a = c*n*0.1;
h2b = c*n*0.6;

w1 = c*n*0.1;
w2 = c*n*0.8;
ws = c;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {w1, 0, 0, 0.1};
Point(3) = {w1+ws, 0, 0, 0.1};
Point(4) = {w1+ws+w2, 0, 0, 0.1};

Point(5) = {0, h1, 0, 0.1};
Point(6) = {w1, h1, 0, 0.1};
Point(7) = {w1+ws, h1, 0, 0.1};
Point(8) = {w1+ws+w2, h1, 0, 0.1};

Point(9) = {0, h1+h2a, 0, 0.1};
Point(10) = {w1, h1+h2a, 0, 0.1};
Point(11) = {w1+ws, h1+h2b, 0, 0.1};
Point(12) = {w1+ws+w2, h1+h2b, 0, 0.1};

Line(1) = {1,5};
Line(2) = {5,9};
Line(3) = {9,10};
Line(4) = {10,6};
Line(5) = {6,7};
Line(6) = {7,11};
Line(7) = {11,12};
Line(8) = {12,8};
Line(9) = {8,4};
Line(10) = {4,3};
Line(11) = {3,2};
Line(12) = {2,1};
Line(13) = {2,6};
Line(14) = {3,7};
Line(15) = {7,8};
Line(16) = {6,5};

Line Loop(17) = {3,4,16,2};
Plane Surface(18) = {17};
Line Loop(19) = {16,-1,-12,13};
Plane Surface(20) = {19};
Line Loop(21) = {14,-5,-13,-11};
Plane Surface(22) = {21};
Line Loop(23) = {15,9,10,14};
Plane Surface(24) = {23};
Line Loop(25) = {6,7,8,-15};
Plane Surface(26) = {25};

Transfinite Line {7,15,10} = w2/c+1 Using Progression 1.0;
Transfinite Line {4,2} = h2a/c+1 Using Progression 1.0;
Transfinite Line {8,6} = h2b/c+1 Using Progression 1.0;
Transfinite Line {9,14,13,1} = h1/c+1 Using Progression 1.0;
Transfinite Line {3,16,12} = w1/c+1 Using Progression 1.0;
Transfinite Line {5,11} = ws/c+1 Using Progression 1.0;

Transfinite Surface {26} = {12,11,7,8};
Transfinite Surface {18} = {10,9,5,6};
Transfinite Surface {20} = {6,5,1,2};
Transfinite Surface {22} = {7,6,2,3};
Transfinite Surface {24} = {3,4,8,7};

Recombine Surface {18,26,20,22,24};

Extrude Surface {18, {0,0,dz} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {26, {0,0,dz} } {
  Layers { {nz}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {20, {0,0,dz} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {22, {0,0,dz} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {24, {0,0,dz} } {
  Layers { {nz}, {9001}, {1.0} };
  Recombine;
};

// Empty
Physical Surface(137) = {48,18,20,92,136,24,70,26};

// Walls
Physical Surface(138) = {47,39,87,83,113,105,57,65,127,131};

// Inlet 
Physical Surface(139) = {35};

// Outlet
Physical Surface(140) = {61};

Physical Volume(160) = {9000};
Physical Volume(161) = {9001};

