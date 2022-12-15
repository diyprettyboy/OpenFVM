Mesh.MshFileVersion=1;

ri = 0.03;
re = 0.06;
h = 0.03;

nc = 15;
nl = 6;
nh = 6;

Point(1) = {0.0,0.0,0.0,0.1};

Point(2) = {ri,0.,0,0.1};
Point(3) = {0,ri,0,0.1};
Point(4) = {-ri,0.,0,0.1};
Point(5) = {0,-ri,0,0.1};

Point(6) = {re,0.,0,0.1};
Point(7) = {0,re,0,0.1};
Point(8) = {-re,0.,0,0.1};
Point(9) = {0,-re,0,0.1};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};
Line(9) = {6,2};
Line(10) = {5,9};
Line(11) = {4,8};
Line(12) = {3,7};

Line Loop(13) = {5,-12,-1,-9};
Plane Surface(14) = {13};
Line Loop(15) = {6,-11,-2,12};
Plane Surface(16) = {15};
Line Loop(17) = {3,10,-7,-11};
Plane Surface(18) = {17};
Line Loop(19) = {4,-9,-8,-10};
Plane Surface(20) = {19};

Transfinite Line {4,8,5,1,2,6,3,7} = nc Using Progression 1.0;
Transfinite Line {10,9,12,11} = nl Using Progression 1.0;

Transfinite Surface {14} = {2,6,7,3};
Transfinite Surface {16} = {3,7,8,4};
Transfinite Surface {18} = {4,8,9,5};
Transfinite Surface {20} = {5,9,6,2};
Recombine Surface {14,16,18,20};

Extrude Surface {14, {0,0,h} } {
  Layers { {nh}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {16, {0,0,h} } {
  Layers { {nh}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {18, {0,0,h} } {
  Layers { {nh}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {20, {0,0,h} } {
  Layers { {nh}, {9000}, {1.0} };
  Recombine;
};

// Top
Physical Surface(100000) = {42,64,108,86};

// Inner wall
Physical Surface(100001) = {59,37,95,73};

// Outer wall
Physical Surface(100002) = {81,103,29,51};

// Bottom
Physical Surface(100003) = {14,20,18,16};

Physical Volume(9001) = {9000};

