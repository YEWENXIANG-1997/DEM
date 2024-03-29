// Gmsh project created on Tue May 19 12:18:13 2020

//radius [m]
r = 0.0025;

//height [m]
h = 0.0027;

//+
Point(1) = {0.0, 0.0, 0.0, 1.0};
//+
Point(2) = {r, 0.0, 0.0, 1.0};
//+
Point(3) = {-r, 0.0, 0.0, 1.0};
//+
Point(4) = {0.0, r, 0.0, 1.0};
//+
Point(5) = {0.0, -r, 0.0, 1.0};
//+
Point(6) = {0.0, -r, 0.0, 1.0};
//+
Point(7) = {0.0, -r, 0.0, 1.0};
//+
Point(8) = {0, 0, -0, 1.0};
//+
Point(9) = {0, 0, -0, 1.0};
//+
Point(10) = {0, 0, -0, 1.0};
//+
Point(11) = {0, 0, -0, 1.0};
//+
Point(12) = {0, 0, -0, 1.0};
//+
Point(13) = {0, 0, -0, 1.0};
//+
Point(14) = {0, 0, -0, 1.0};
//+
Circle(1) = {5, 1, 2};
//+
Circle(2) = {2, 1, 4};
//+
Circle(3) = {4, 1, 3};
//+
Circle(4) = {3, 1, 5};
//+
Extrude {0, 0, h} {
  Curve{2}; Curve{3}; Curve{4}; Curve{1}; 
}
//+
Curve Loop(5) = {6, -7, -5, 2};
//+
Plane Surface(5) = {5};
//+
Curve Loop(7) = {5, -12, -10, 1};
//+
Plane Surface(6) = {7};
//+
Curve Loop(8) = {10, -11, -8, 4};
//+
Plane Surface(7) = {8};
//+
Curve Loop(10) = {8, -9, -6, 3};
//+
Plane Surface(8) = {10};

Mesh.CharacteristicLengthMin = 0.;
Mesh.CharacteristicLengthMax = 0.001;
