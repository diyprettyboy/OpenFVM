Mesh.MshFileVersion=1;
lc1=0.1;
lc=0.05;
lc2=0.02;
nT=15;

Point(1) = {0.0,0.0,0.0,lc2};
Point(2) = {1,0.0,0.0,lc};
Point(3) = {1,0.5,0.0,lc2};
Point(4) = {1,0.7,0.0,lc2};
Point(5) = {0,0.7,0.0,lc1};
Point(6) = {0,0.1,0.0,lc2};
Point(7) = {-0.3,0.1,0.0,lc2};
Point(8) = {-0.3,0.,0.0,lc2};
Point(9) = {1.2,0.7,0.0,lc2};
Point(10) = {1.2,0.5,0.0,lc2};
Line(1) = {5,6};
Line(2) = {6,7};
Line(3) = {7,8};
Line(4) = {8,1};
Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,10};
Line(8) = {10,9};
Line(9) = {9,4};
Line(10) = {4,5};
Line Loop(11) = {10,1,2,3,4,5,6,7,8,9};
Plane Surface(12) = {11};
Extrude Surface {12, {0.0,0.0,0.5}} { Recombine; Layers { {nT}, {9000}, {1} }; };

// Empty : 64
Physical Surface(64) = {64,12};

// Walls : 65
Physical Surface(65) = {27,31,47,51,55,63,35,43};

// Inlet : 67
Physical Surface(67) = {39};

// Outlet : 66
Physical Surface(66) = {59};

Physical Volume(70) = {9000};
