SetFactory("OpenCASCADE");

// Outer rectangle
Rectangle(1) = {0, 0, 0, 6.283185,3.14159, 0}; // Azimuthal x Polar


// Assign physical regions
Physical Surface("solid_angle") = {1};
