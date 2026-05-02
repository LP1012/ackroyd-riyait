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
