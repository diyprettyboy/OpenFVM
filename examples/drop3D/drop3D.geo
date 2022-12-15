//Mesh.MshFileVersion=1;

lc = 0.5;

n = 5;

d = 1.0/5*n;

w = 1.0;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {w,0.0,0.0,lc};
Point(3) = {w,0.0,w,lc};
Point(4) = {0.0,0.0,w,lc};

Point(5) = {w*0.5,0.0,w*0.5,lc};

Point(6) = {0.4*w,0,0.6*w,0.1};
Point(7) = {0.6*w,0,0.6*w,0.1};
Point(8) = {0.6*w,0,0.4*w,0.1};
Point(9) = {0.4*w,0,0.4*w,0.1};
Point(10) = {0,0,0.4*w,0.1};
Point(11) = {0.4*w,0,0,0.1};
Point(12) = {0.6*w,0,0,0.1};
Point(13) = {1*w,0,0.4*w,0.1};
Point(14) = {1*w,0,0.6,0.1};
Point(15) = {0,0,0.6*w,0.1};
Point(16) = {0.4*w,0,1*w,0.1};
Point(17) = {0.6*w,0,1*w,0.1};

Line(1) = {6,9};
Line(2) = {9,8};
Line(3) = {8,7};
Line(4) = {7,6};
Line(5) = {17,16};
Line(6) = {16,6};
Line(7) = {7,17};
Line(8) = {17,3};
Line(9) = {3,14};
Line(10) = {14,7};
Line(11) = {14,13};
Line(12) = {13,8};
Line(13) = {8,12};
Line(14) = {12,2};
Line(15) = {2,13};
Line(16) = {12,11};
Line(17) = {11,9};
Line(18) = {9,10};
Line(19) = {10,1};
Line(20) = {1,11};
Line(21) = {10,15};
Line(22) = {6,15};
Line(23) = {15,4};
Line(24) = {4,16};
Line Loop(25) = {23,24,6,22};
Plane Surface(26) = {25};
Line Loop(27) = {6,-4,7,5};
Plane Surface(28) = {27};
Line Loop(29) = {4,1,2,3};
Plane Surface(30) = {29};
Line Loop(31) = {1,18,21,-22};

Plane Surface(32) = {31};
Line Loop(33) = {19,20,17,18};
Plane Surface(34) = {33};
Line Loop(35) = {2,13,16,17};
Plane Surface(36) = {35};
Line Loop(37) = {14,15,12,13};
Plane Surface(38) = {37};
Line Loop(39) = {3,-10,11,12};
Plane Surface(40) = {39};
Line Loop(41) = {9,10,7,8};
Plane Surface(42) = {41};

Transfinite Line {21,1,3,11,4,2,16,5} = n+1 Using Progression 1;
Transfinite Line {-23,6,-7,9,-8,24,15,12,10,-13,-14,17,20,-19,-18,-22} = 5*n*0.5+1 Using Progression 1.0;

Transfinite Surface {26} = {16,4,15,6};
Transfinite Surface {28} = {17,16,6,7};
Transfinite Surface {42} = {3,17,7,14};
Transfinite Surface {40} = {14,7,8,13};
Transfinite Surface {30} = {7,6,9,8};
Transfinite Surface {32} = {6,15,10,9};
Transfinite Surface {38} = {13,8,12,2};
Transfinite Surface {36} = {8,9,11,12};
Transfinite Surface {34} = {9,10,1,11};

Recombine Surface {26,28,42,40,30,32,34,36,38};

Extrude Surface {42, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {40, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {30, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {28, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {38, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {36, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {26, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {32, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {34, {0,w*0.05,0} } {
  Layers { {5*n*0.05}, {2}, {1.0} };
  Recombine;
};

Extrude Surface {152, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {174, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {218, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {108, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {86, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {64, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {130, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {196, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {240, {0,w*0.45,0} } {
  Layers { {5*n*0.45}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {262, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {284, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {438, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {306, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {328, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {350, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {372, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {394, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

Extrude Surface {416, {0,w*0.5,0} } {
  Layers { {5*n*0.5}, {1}, {1.0} };
  Recombine;
};

// Surfaces
Physical Surface(1) = {504,460,482,570,614,548,592,636,526};
Physical Surface(2) = {627,613,591,371,63,129,393,407,187};
Physical Surface(3) = {579,565,451,253,143,345,81,359,51};
Physical Surface(4) = {447,249,139,169,279,477,495,429,231};
Physical Surface(5) = {491,521,623,183,403,301,213,425,227};
Physical Surface(6) = {34,32,26,28,42,40,30,36,38};

// Box volume
Physical Volume (1) = {1};
Physical Volume (2) = {2};





