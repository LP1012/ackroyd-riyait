SetFactory("OpenCASCADE");

// Outer rectangle
Rectangle(1) = {0, 0, 0, 3.14159, 1.5707963, 0}; // Azimuthal x Polar


// Assign physical regions
Physical Surface("solid_angle") = {1};
