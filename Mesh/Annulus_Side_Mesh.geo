lc = 0.05;
Point(1) = { 2,  0, 0, lc};
Point(2) = {-2,  0, 0, lc};
Point(3) = { 0,  2, 0, lc};
Point(4) = { 0, -2, 0, lc};
Point(5) = { 0,  0, 0, lc};
Point(6) = { 1,  0, 0, lc};
Point(7) = {-1,  0, 0, lc};
Point(8) = { 0,  1, 0, lc};
Point(9) = { 0, -1, 0, lc};
Circle(1) = {3, 5, 1};
Circle(2) = {1, 5, 4};
Circle(3) = {4, 5, 2};
Circle(4) = {2, 5, 3};
Circle(5) = {8, 5, 6};
Circle(6) = {6, 5, 9};
Circle(7) = {9, 5, 7};
Circle(8) = {7, 5, 8};

Line(9) = {6, 1};
Line(10) = {7, 2};
Line(11) = {8, 3};
Line(12) = {9, 4};

Curve Loop(1) = {9, -1, -11, 5};
Curve Loop(2) = {11, -4, -10, 8};
Curve Loop(3) = {10, -3, -12, 7};
Curve Loop(4) = {6, 12, -2, -9};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Physical Curve("Inner", 11) = {8, 7, 6, 5};
Physical Curve("Outer", 12) = {1, 4, 3, 2};
Physical Curve("MagSplit", 1) = {9};

Physical Surface("Annulus", 5) = {1,2,3,4};