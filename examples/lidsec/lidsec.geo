Mesh.MshFileVersion=1;
n = 1;

Point(1) = {-0.043118,-0.055012,0.000000,0.1};
Point(23) = {0.006882,-0.055012,0.000000,0.1};
Point(71) = {0.006882,-0.020012,0.000000,0.1};
Point(45) = {0.056882,-0.020012,0.000000,0.1};
Point(3) = {0.056882,0.044988,0.000000,0.1};
Point(2) = {-0.043118,0.044988,0.000000,0.1};
Point(10) = {-0.043118,-0.020012,0.000000,0.1};
Point(55) = {0.006882,0.044988,0.000000,0.1};
Line(1) = {55,3};
Line(2) = {3,45};
Line(3) = {45,71};
Line(4) = {71,55};
Line(5) = {55,2};
Line(6) = {2,10};
Line(7) = {10,71};
Line(8) = {71,23};
Line(9) = {23,1};
Line(10) = {1,10};
Line Loop(11) = {5,6,7,4};
Plane Surface(12) = {11};
Line Loop(13) = {8,9,10,7};
Plane Surface(14) = {13};
Line Loop(15) = {4,1,2,3};
Plane Surface(16) = {15};

Transfinite Line {9,7,5,1,3} = 10*n + 1 Using Progression 1.0;
Transfinite Line {8,10} = 7*n + 1 Using Progression 1.0;
Transfinite Line {2,4,6} = 13*n + 1 Using Progression 1.0;

Transfinite Surface {16} = {3,55,71,45};
Transfinite Surface {12} = {55,2,10,71};
Transfinite Surface {14} = {10,1,23,71};
Recombine Surface {16,12,14};

Extrude Surface {12, {0,0,0.5/(10*n)} } {
  Layers { {1}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {14, {0,0,0.5/(10*n)} } {
  Layers { {1}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {16, {0,0,0.5/(10*n)} } {
  Layers { {1}, {9000}, {1.0} };
  Recombine;
};

// Walls
Physical Surface(87) = {38,12,82,16,77,81,47,60,14,51,55,29};

// lid
Physical Surface(88) = {25,73};

Physical Volume(300) = {9000};
