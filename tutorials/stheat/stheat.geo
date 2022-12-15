Mesh.MshFileVersion=1;
lc = 0.02;
//lc = 0.01;

n = 14;

dimx = 0.2;
dimy = 0.1;
dimz = 0.04;

Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, dimz, lc};
Point(3) = {0, dimy, 0, lc};
Point(4) = {0, dimy, dimz, lc};
Point(5) = {dimx, 0, 0, lc};
Point(6) = {dimx, 0, dimz, lc};
Point(7) = {dimx, dimy, 0, lc};
Point(8) = {dimx, dimy, dimz, lc};
Line (9) = {2, 4};
Line (10) = {1, 2};
Line (11) = {2, 6};
Line (12) = {5, 6};
Line (13) = {1, 5};
Line (14) = {3, 4};
Line (15) = {4, 8};
Line (16) = {7, 8};
Line (17) = {3, 7};
Line (18) = {6, 8};
Line (19) = {5, 7};
Line (20) = {1, 3};
Line Loop (1000027) = {10, 11, -12, -13};
Plane Surface (27) = {1000027};
Line Loop (1000028) = {12, 18, -16, -19};
Plane Surface (28) = {1000028};
Line Loop (1000029) = {14, -17, -16, 15};
Plane Surface (29) = {1000029};
Line Loop (1000030) = {9, 10, -20, -14};
Plane Surface (30) = {1000030};
Line Loop (1000031) = {9, 15, -18, -11};
Plane Surface (31) = {-1000031};
Line Loop (1000032) = {13, 19, -17, -20};
Plane Surface (32) = {1000032};

Surface Loop (1000034) = {27, 28, 29, 30, 31, 32};
Volume (34) = {1000034};

//Transfinite Line {15,17,19,18,11,13,9,20} = n+1 Using Progression 1.0;
//Transfinite Line {12,16,14,10} = 1+1 Using Progression 1.0;

//Transfinite Surface {32} = {7,3,1,5};
//Transfinite Surface {28} = {7,5,6,8};
//Transfinite Surface {31} = {6,2,4,8};
//Transfinite Surface {30} = {4,2,1,3};
//Transfinite Surface {27} = {1,5,6,2};
//Transfinite Surface {29} = {8,7,3,4};

//Transfinite Volume{34} = {7,8,4,3,5,6,2,1};

Physical Surface (35) = {29};
Physical Surface (36) = {30};
Physical Surface (37) = {28};
Physical Surface (38) = {32};
Physical Surface (39) = {27};
Physical Surface (40) = {31};

Physical Volume (41) = {34};
