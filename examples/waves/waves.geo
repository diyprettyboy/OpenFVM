//Mesh.MshFileVersion=1;

n = 6;

H = 0.05175;

dimx = 3*H;
dimy = 1.2*H;
dimz = 0.05;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {dimx, 0, 0, 0.1};
Point(3) = {dimx, dimy, 0, 0.1};
Point(4) = {0, dimy, 0, 0.1};

Point(5) = {H, 0, 0, 0.1};
Point(6) = {dimx, H, 0, 0.1};
Point(7) = {H, dimy, 0, 0.1};
Point(8) = {0, H, 0, 0.1};

Point(9) = {H, H, 0, 0.1};

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

Transfinite Line {12,10,8} = n+1 Using Progression 1.0;
Transfinite Line {1,3,6} = 5*n+1 Using Progression 1.0;
Transfinite Line {9,7,5} = 10*n+1 Using Progression 1.0;
Transfinite Line {4,2,11} = 5*n+1 Using Progression 1.0;

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
  Layers { {1}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {16, {0,0,dimz} } {
  Layers { {1}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {18, {0,0,dimz} } {
  Layers { {1}, {9000}, {1.0} };
  Recombine;
};

Extrude Surface {20, {0,0,dimz} } {
  Layers { {1}, {9001}, {1.0} };
  Recombine;
};

// Top 
Physical Surface(140) = {41,59};

// Bottom
Physical Surface(141) = {103,77};

// Left
Physical Surface(142) = {37,107};

// Right
Physical Surface(143) = {81,63};

// Front
Physical Surface(144) = {42,64,108,86};

// Back
Physical Surface(145) = {14,16,18,20};

Physical Volume(150) = {9000};
Physical Volume(151) = {9001};




