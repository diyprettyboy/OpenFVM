Mesh.MshFileVersion=1;
lc = 0.05;
T = 0.5;
nT = 1;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.6, 0, lc};
Point(4) = {1, 0.8, 0, lc};
Point(5) = {0, 0.8, 0, lc};
Point(6) = {1.2, 0.8, 0, lc};
Point(7) = {1.2, 0.6, 0, lc};
Point(8) = {0, 0.6, 0, lc};
Point(9) = {0.8, 0.6, 0, lc};
Point(10) = {0.8, 0, 0, lc};
Point(11) = {0.8, 0.8, 0, lc};
Point(12) = {0.8, -0.4, 0, lc};
Point(13) = {1.0, -0.4, 0, lc};
Line(1) = {6,4};
Line(2) = {4,3};
Line(3) = {3,7};
Line(4) = {7,6};
Line(5) = {4,11};
Line(6) = {11,9};
Line(7) = {9,3};
Line(8) = {9,10};
Line(9) = {10,2};
Line(10) = {2,3};
Line(11) = {9,8};
Line(12) = {8,5};
Line(13) = {5,11};
Line(14) = {8,1};
Line(15) = {1,10};
Line(16) = {10,12};
Line(17) = {12,13};
Line(18) = {13,2};
Line Loop(19) = {3,4,1,2};
Plane Surface(20) = {19};
Line Loop(21) = {7,-2,5,6};
Plane Surface(22) = {21};
Line Loop(23) = {13,6,11,12};
Plane Surface(24) = {23};
Line Loop(25) = {14,15,-8,11};
Plane Surface(26) = {25};
Line Loop(27) = {8,9,10,-7};
Plane Surface(28) = {27};
Line Loop(29) = {16,17,18,-9};
Plane Surface(30) = {29};
Transfinite Surface {20} = {7,6,4,3};
Transfinite Surface {22} = {3,4,11,9};
Transfinite Surface {28} = {2,3,9,10};
Transfinite Surface {30} = {13,2,10,12};
Transfinite Surface {26} = {1,10,9,8};
Transfinite Surface {24} = {9,11,5,8};
Recombine Surface {20,22,24,26,28,30};
Transfinite Line {13,11,15} = 32 Using Progression 1.0;
Transfinite Line {12,6,2,4} = 8 Using Progression 1.0;
Transfinite Line {5,7,9,17} = 8 Using Progression 1.0;
Transfinite Line {1,3} = 8 Using Progression 1.0;
Transfinite Line {14,8,10} = 24 Using Progression 1.0;
Transfinite Line {16,18} = 16 Using Progression 1.0;

Extrude Surface {20, {0,0,T}}
{ Recombine; Layers { {nT}, {9000}, {nT} }; };
Extrude Surface {22, {0,0,T}}
{ Recombine; Layers { {nT}, {9001}, {nT} }; };
Extrude Surface {24, {0,0,T}}
{ Recombine; Layers { {nT}, {9002}, {nT} }; };
Extrude Surface {26, {0,0,T}}
{ Recombine; Layers { {nT}, {9003}, {nT} }; };
Extrude Surface {28, {0,0,T}}
{ Recombine; Layers { {nT}, {9004}, {nT} }; };
Extrude Surface {30, {0,0,T}}
{ Recombine; Layers { {nT}, {9005}, {nT} }; };

// Walls
Physical Surface(163) = {52,20,47,39,74,22,69,140,28,83,95,96,24,105,118,26,109,135,149,162,30,157};

// Inlet 
Physical Surface(164) = {43};
// Outlet
Physical Surface(165) = {153};

Physical Volume(200) = {9000:9005};
