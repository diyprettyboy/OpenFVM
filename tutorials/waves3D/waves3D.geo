Mesh.MshFileVersion=1;

n = 40;

dimx = 1.0;
dimy = 1.0;
dimz = 0.5;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {dimx, 0, 0, 0.1};
Point(3) = {dimx, dimy, 0, 0.1};
Point(4) = {0, dimy, 0, 0.1};

Point(5) = {dimx*0.5, 0, 0, 0.1};
Point(6) = {dimx, dimy*0.75, 0, 0.1};
Point(7) = {dimx*0.5, dimy, 0, 0.1};
Point(8) = {0, dimy*0.75, 0, 0.1};

Point(9) = {dimx*0.5, dimy*0.75, 0, 0.1};

Line(1) = {1,8};
Line(2) = {8,9};
Line(3) = {9,5};
Line(4) = {5,1};
Line(5) = {5,2};
Line(6) = {2,6};
Line(7) = {6,9};
Line(8) = {6,3};
Line(9) = {3,7};
Line(10) = {7,9};
Line(11) = {7,4};
Line(12) = {4,8};

Line Loop(13) = {10,-2,-12,-11};
Plane Surface(14) = {13};
Line Loop(15) = {7,-10,-9,-8};
Plane Surface(16) = {15};
Line Loop(17) = {3,5,6,7};
Plane Surface(18) = {17};
Line Loop(19) = {2,3,4,1};
Plane Surface(20) = {19};

Transfinite Line {12,10,8} = 0.25*n+1 Using Progression 1.0;
Transfinite Line {1,3,6} = 0.75*n+1 Using Progression 1.0;
Transfinite Line {9,7,5} = 0.5*n+1 Using Progression 1.0;
Transfinite Line {4,2,11} = 0.5*n+1 Using Progression 1.0;

Transfinite Surface {14} = {7,4,8,9};
Transfinite Surface {20} = {8,1,5,9};
Transfinite Surface {18} = {6,9,5,2};
Transfinite Surface {16} = {6,3,7,9};

Recombine Surface {14,16,18,20};

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

Extrude Surface {14, {0,0,dimz} } {
  Layers { {0.5*n}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {16, {0,0,dimz} } {
  Layers { {0.5*n}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {18, {0,0,dimz} } {
  Layers { {0.5*n}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {20, {0,0,dimz} } {
  Layers { {0.5*n}, {9001}, {1.0} };
  Recombine;
};

Extrude Surface {42, {0,0,dimz} } {
  Layers { {0.5*n}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {64, {0,0,dimz} } {
  Layers { {0.5*n}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {108, {0,0,dimz} } {
  Layers { {0.5*n}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {86, {0,0,dimz} } {
  Layers { {0.5*n}, {9000}, {1.0} };
  Recombine;
};

// Top 
Physical Surface(140) = {59,147,41,129};

// Bottom
Physical Surface(141) = {169,187,103,77};

// Left
Physical Surface(142) = {173,125,37,107};

// Right
Physical Surface(143) = {151,63,81,191};

// Front
Physical Surface(144) = {174,130,152,196};

// Back
Physical Surface(145) = {18,20,14,16};

Physical Volume(150) = {9000};
Physical Volume(151) = {9001};

