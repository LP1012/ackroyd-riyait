SetFactory("OpenCASCADE");

// Outer rectangle
Rectangle(1) = {0, 0, 0, 10, 10, 0}; // Outer material
Rectangle(2) = {0, 0, 0, 5, 5, 0}; // Void
Rectangle(3) = {0, 0, 0, 1.25, 1.25, 0}; // Source

BooleanFragments{Surface{1, 2, 3}; Delete;}{}

// Assign physical regions
Physical Surface("material") = {5};
Physical Surface("void") = {4};
Physical Surface("source") = {3};
//+
Physical Curve("west", 13) = {7, 3, 11};
//+
Physical Curve("north", 14) = {8};
//+
Physical Curve("east", 15) = {9};
//+
Physical Curve("south", 16) = {10, 6, 12};
