
size = 0.002;

z = 0.006 ;

Point(1)  = {0.048, 0.008, 0, size};
Point(2)  = {0.048, 0.012, 0, size/2};
Point(3)  = {0.048, 0.014, 0, size/2};
Point(4)  = {0.04,  0.014, 0, size/2};
Point(5)  = {0.04,  0.012, 0, size/2};
Point(6)  = {0.035, 0.012, 0, size};
Point(7)  = {0.021, 0.004, 0, size*2};
Point(8)  = {0,     0.004, 0, size*2};
Point(9)  = {0,     0,     0, size*2};
Point(10) = {0.022, 0,     0, size*2};
Point(11) = {0.036, 0.008, 0, size};
Point(12) = {0.058, 0.008, 0, size};
Point(13) = {0.058, 0.010, 0, size};
Point(14) = {0.048, 0.010, 0, size};

Point(21) = {0.048, 0.008, z, size};
Point(22) = {0.048, 0.012, z, size/2};
Point(23) = {0.048, 0.014, z, size/2};
Point(24) = {0.04,  0.014, z, size/2};
Point(25) = {0.04,  0.012, z, size/2};
Point(26) = {0.035, 0.012, z, size};
Point(27) = {0.021, 0.004, z, size*2};
Point(28) = {0,     0.004, z, size*2};
Point(29) = {0,     0,     z, size*2};
Point(30) = {0.022, 0,     z, size*2};
Point(31) = {0.036, 0.008, z, size};
Point(32) = {0.058, 0.008, z, size};
Point(33) = {0.058, 0.010, z, size};
Point(34) = {0.048, 0.010, z, size};

/* pts suppl dans la pastille */

Point(35) = {0.048, 0.014, z/3, size/2};
Point(36) = {0.04,  0.014, z/3, size/2};
Point(37) = {0.048, 0.014, 2*z/3, size/2};
Point(38) = {0.04,  0.014, 2*z/3, size/2};

/* trou de mesure */

Point(40) = {0.044, 0.008, z/2, size/2};
Point(41) = {0.044, 0.011, z/2, size/2};
Point(42) = {0.044, 0.012, z/2, size/2};

Point(43) = {0.0435, 0.008, z/2, size/4};
Point(44) = {0.0445, 0.008, z/2, size/4};
Point(45) = {0.0435, 0.011, z/2, size/4};
Point(46) = {0.0445, 0.011, z/2, size/4};
Point(47) = {0.0435, 0.012, z/2, size/4};
Point(48) = {0.0445, 0.012, z/2, size/4};

Point(49) = {0.044, 0.008, z/2-0.0005, size/4};
Point(50) = {0.044, 0.008, z/2+0.0005, size/4};
Point(51) = {0.044, 0.011, z/2-0.0005, size/4};
Point(52) = {0.044, 0.011, z/2+0.0005, size/4};
Point(53) = {0.044, 0.012, z/2-0.0005, size/4};
Point(54) = {0.044, 0.012, z/2+0.0005, size/4};


/* trou trou de droite */

Point(60) = {0.053, 0.008, z/2, size/2};
Point(61) = {0.053, 0.01,  z/2, size/2};

Point(62) = {0.052, 0.008, z/2, size/2};
Point(63) = {0.054, 0.008, z/2, size/2};
Point(64) = {0.052, 0.01,  z/2, size/2};
Point(65) = {0.054, 0.01,  z/2, size/2};

Point(66) = {0.053, 0.008, z/2-0.001, size/2};
Point(67) = {0.053, 0.008, z/2+0.001, size/2};
Point(68) = {0.053, 0.01,  z/2-0.001, size/2};
Point(69) = {0.053, 0.01,  z/2+0.001, size/2};


Line(101) = {9,10};
Line(102) = {10,11};
Line(103) = {11,1};
Line(104) = {1,12};
Line(105) = {12,32};
Line(106) = {32,21};
Line(107) = {21,31};
Line(108) = {31,30};
Line(109) = {30,29};
Line(110) = {29,9};
Line(111) = {30,10};
Line(112) = {31,11};
Line(113) = {29,28};
Line(114) = {8,28};
Line(115) = {8,9};
Line(116) = {8,7};
Line(117) = {7,10};
Line(118) = {27,7};
Line(119) = {27,30};
Line(120) = {27,28};
Line(121) = {27,26};
Line(122) = {7,6};
Line(123) = {6,26};
Line(124) = {26,31};
Line(125) = {11,6};
Line(126) = {26,25};
Line(127) = {25,5};
Line(128) = {5,6};
Line(129) = {25,24};
Line(130) = {24,38};
Line(131) = {38,36};
Line(132) = {36,4};
Line(133) = {4,5};
Line(134) = {22,2};
Line(135) = {2,3};
Line(136) = {3,35};
Line(137) = {35,37};
Line(138) = {37,23};
Line(139) = {23,22};
Line(140) = {24,23};
Line(141) = {37,38};
Line(142) = {36,35};
Line(143) = {3,4};
Line(144) = {25,22};
Line(145) = {5,2};
Line(146) = {14,2};
Line(147) = {22,34};
Line(148) = {34,14};
Line(149) = {34,33};
Line(150) = {33,13};
Line(151) = {13,14};
Line(152) = {33,32};
Line(153) = {12,13};
Line(154) = {21,34};
Line(155) = {1,14};
Line(156) = {21,1};
Circle(157) = {48,42,53};
Circle(158) = {53,42,47};
Circle(159) = {47,42,54};
Circle(160) = {54,42,48};
Circle(161) = {51,41,45};
Circle(162) = {45,41,52};
Circle(163) = {52,41,46};
Circle(164) = {46,41,51};
Circle(165) = {44,40,49};
Circle(166) = {49,40,43};
Circle(167) = {43,40,50};
Circle(168) = {50,40,44};
Circle(169) = {65,61,68};
Circle(170) = {68,61,64};
Circle(171) = {64,61,69};
Circle(172) = {69,61,65};
Circle(173) = {66,60,62};
Circle(174) = {62,60,67};
Circle(175) = {67,60,63};
Circle(176) = {63,60,66};
Line(177) = {64,62};
Line(178) = {69,67};
Line(179) = {63,65};
Line(180) = {66,68};
Line(181) = {43,45};
Line(182) = {49,51};
Line(183) = {44,46};
Line(184) = {50,52};
Line(185) = {52,54};
Line(186) = {46,48};
Line(187) = {51,53};
Line(188) = {45,47};

Line Loop(189) = {-114,115,-110,113};
Plane Surface(190) = {189};
Line Loop(191) = {-114,116,-118,120};
Plane Surface(192) = {191};
Line Loop(193) = {-101,-115,116,117};
Plane Surface(194) = {193};
Line Loop(195) = {109,113,-120,119};
Plane Surface(196) = {195};
Line Loop(197) = {-101,-110,-109,111};
Plane Surface(198) = {197};
Line Loop(199) = {-123,-122,-118,121};
Plane Surface(200) = {199};
Line Loop(201) = {-122,117,102,125};
Plane Surface(202) = {201};
Line Loop(203) = {108,-119,121,124};
Plane Surface(204) = {203};
Line Loop(205) = {-112,108,111,102};
Plane Surface(206) = {205};
Line Loop(207) = {-128,145,-146,-155,-103,125};
Plane Surface(208) = {207};
Line Loop(209) = {-103,-112,-107,156};
Line Loop(210) = {168,165,166,167};
Plane Surface(211) = {209,210};
Plane Surface(212) = {210};
Line Loop(213) = {162,163,164,161};
Plane Surface(214) = {213};
Line Loop(215) = {159,160,157,158};
Plane Surface(216) = {215};
Line Loop(217) = {123,126,127,128};
Plane Surface(218) = {217};
Line Loop(219) = {-134,-144,127,145};
Plane Surface(220) = {219,215};
Line Loop(221) = {147,-154,107,-124,126,144};
Plane Surface(222) = {221};
Line Loop(223) = {-148,-154,156,155};
Plane Surface(224) = {223};
Line Loop(225) = {-134,147,148,146};
Plane Surface(226) = {225};
Line Loop(227) = {-127,129,130,131,132,133};
Plane Surface(228) = {227};
Line Loop(229) = {139,-144,129,140};
Plane Surface(230) = {229};
Line Loop(231) = {136,137,138,139,134,135};
Plane Surface(232) = {231};
Line Loop(233) = {145,135,143,133};
Plane Surface(234) = {233};
Line Loop(235) = {-132,142,-136,143};
Plane Surface(236) = {235};
Line Loop(237) = {137,141,131,142};
Plane Surface(238) = {237};
Line Loop(239) = {-130,140,-138,141};
Plane Surface(240) = {239};
Line Loop(241) = {152,106,154,149};
Plane Surface(242) = {241};
Line Loop(243) = {104,105,106,156};
Line Loop(244) = {174,175,176,173};
Plane Surface(245) = {243,244};
Plane Surface(246) = {244};
Line Loop(247) = {-151,-153,-104,155};
Plane Surface(248) = {247};
Line Loop(249) = {-150,152,-105,153};
Plane Surface(250) = {249};
Line Loop(251) = {-151,-150,-149,148};
Line Loop(252) = {171,172,169,170};
Plane Surface(253) = {251,252};
Plane Surface(254) = {252};
Line Loop(255) = {174,-178,-171,177};
Ruled Surface(256) = {255};
Line Loop(257) = {-179,-175,-178,172};
Ruled Surface(258) = {257};
Line Loop(259) = {-169,-179,176,180};
Ruled Surface(260) = {259};
Line Loop(261) = {-173,180,170,177};
Ruled Surface(262) = {261};
Line Loop(263) = {162,-184,-167,181};
Ruled Surface(264) = {263};
Line Loop(265) = {-183,-168,184,163};
Ruled Surface(266) = {265};
Line Loop(267) = {182,-164,-183,165};
Ruled Surface(268) = {267};
Line Loop(269) = {-181,-166,182,161};
Ruled Surface(270) = {269};
Line Loop(271) = {-159,-188,162,185};
Ruled Surface(272) = {271};
Line Loop(273) = {-160,-185,163,186};
Ruled Surface(274) = {273};
Line Loop(275) = {157,-187,-164,186};
Ruled Surface(276) = {275};
Line Loop(277) = {188,-158,-187,161};
Ruled Surface(278) = {277};

Surface Loop(279) = {266,-268,270,264,-272,-220,226,-222,224,211,208,218,200,-202,194,-198,190,-192,-196,-204,-206,-274,276,-278};
Complex Volume(280) = {279};

Surface Loop(281) = {258,-260,-253,248,250,-242,245,-224,-256,262};
Complex Volume(282) = {281};

Surface Loop(283) = {214,-272,-216,-274,276,-278};
Complex Volume(284) = {283};

Surface Loop(285) = {216,-220,-232,-236,-228,230,-240,238,234};
Complex Volume(286) = {285};



/* entree flux */

Physical Surface(1000) = {236,238 /*,240 */ };    

/* pastille */
Physical Volume(2000) = {286}; 

/* cuivre */
Physical Volume(3000) = {280,282 /* ,284 */};

