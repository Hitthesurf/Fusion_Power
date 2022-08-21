lc = 0.1;
//Inner and Outer Circle Points
Point(1) = { 1.5,    0,    0, lc};
Point(2) = {   0,  1.5,    0, lc};
Point(3) = {-1.5,    0,    0, lc};
Point(4) = {   0, -1.5,    0, lc};
Point(5) = { 0.5,    0,    0, lc};
Point(6) = {   0,  0.5,    0, lc};
Point(7) = {-0.5,    0,    0, lc};
Point(8) = {   0, -0.5,    0, lc};

//Circle Above and Below Points
Point(9)  = {  1.,    0,  0.5, lc};
Point(10) = {   0,   1.,  0.5, lc};
Point(11) = { -1.,    0,  0.5, lc};
Point(12) = {   0,  -1.,  0.5, lc};
Point(13) = {  1.,    0, -0.5, lc};
Point(14) = {   0,   1., -0.5, lc};
Point(15) = { -1.,    0, -0.5, lc};
Point(16) = {   0,  -1., -0.5, lc};

//Middle Points
Point(17) = {  0.,   0.,   0., lc};
Point(18) = {  0.,   0.,  0.5, lc};
Point(19) = {  0.,   0., -0.5, lc};
Point(20) = {  1.,   0.,   0., lc};
Point(21) = {  0.,   1.,   0., lc};
Point(22) = { -1.,   0.,   0., lc};
Point(23) = {  0.,  -1.,   0., lc};

//Outer Circle
Circle(1) = {1, 17, 2};
Circle(2) = {2, 17, 3};
Circle(3) = {3, 17, 4};
Circle(4) = {4, 17, 1};

//Inner Circle
Circle(5) = {5, 17, 6};
Circle(6) = {6, 17, 7};
Circle(7) = {7, 17, 8};
Circle(8) = {8, 17, 5};

//Upper Circle
Circle(11) = {9,  18, 10};
Circle(12) = {10, 18, 11};
Circle(13) = {11, 18, 12};
Circle(14) = {12, 18, 9};

//Lower Circle
Circle(15) = {13, 19, 14};
Circle(16) = {14, 19, 15};
Circle(17) = {15, 19, 16};
Circle(18) = {16, 19, 13};

//Circles connecting Outer Upper Inner and Lower
Circle(21) = {1, 20, 9};
Circle(22) = {9, 20, 5};
Circle(23) = {5, 20, 13};
Circle(24) = {13, 20, 1};

Circle(25) = {10 ,21, 2};
Circle(26) = {2, 21, 14};
Circle(27) = {14, 21, 6};
Circle(28) = {6, 21, 10};

Circle(29) = {11, 22, 3};
Circle(30) = {3, 22, 15};
Circle(31) = {15, 22, 7};
Circle(32) = {7, 22, 11};

Circle(33) = {12, 23, 8};
Circle(34) = {8, 23, 16};
Circle(35) = {16, 23, 4};
Circle(36) = {4, 23, 12};


Curve Loop(1) = {35, -3, 30, 17};
Surface(1) = {1};
Curve Loop(2) = {34, -17, 31, 7};
Surface(2) = {2};
Curve Loop(3) = {29, 3, 36, -13};
Surface(3) = {3};
Curve Loop(4) = {7, -33, -13, -32};
Surface(4) = {4};
Curve Loop(5) = {30, -16, -26, 2};
Surface(5) = {5};
Curve Loop(6) = {6, -31, -16, 27};
Surface(6) = {6};
Curve Loop(7) = {5, -27, -15, -23};
Surface(7) = {7};
Curve Loop(8) = {24, 1, 26, -15};
Surface(8) = {8};
Curve Loop(9) = {4, -24, -18, 35};
Surface(9) = {9};
Curve Loop(10) = {34, 18, -23, -8};
Surface(10) = {10};
Curve Loop(11) = {2, -29, -12, 25};
Surface(11) = {11};
Curve Loop(12) = {28, 12, -32, -6};
Surface(12) = {12};
Curve Loop(13) = {22, 5, 28, -11};
Surface(13) = {13};
Curve Loop(14) = {25, -1, 21, 11};
Surface(14) = {14};
Curve Loop(15) = {8, -22, -14, 33};
Surface(15) = {15};
Curve Loop(16) = {36, 14, -21, -4};
Surface(16) = {16};

Curve Loop(17) = {36, 33, 34, 35};
Surface(17) = {17};
Curve Loop(18) = {31, 32, 29, 30};
Surface(18) = {18};
Curve Loop(19) = {28, 25, 26, 27};
Surface(19) = {19};
Curve Loop(20) = {24, 21, 22, 23};
Surface(20) = {20};

Surface Loop(1) = {17, 3, 1, 2, 4, 18};
Volume(1) = {1};
Surface Loop(2) = {18, 11, 5, 6, 12, 19};
Volume(2) = {2};
Surface Loop(3) = {19, 7, 13, 14, 8, 20};
Volume(3) = {3};
Surface Loop(4) = {20, 16, 15, 10, 9, 17};
Volume(4) = {4};

Physical Volume("Inside", 37) = {3, 4, 1, 2};
Physical Surface("MagSplit", 9) = {19};
Physical Surface("OuterSurface", 1) = {11, 12, 5, 14, 13, 15, 16, 3, 4, 9, 10, 8, 7, 6, 2, 1};
