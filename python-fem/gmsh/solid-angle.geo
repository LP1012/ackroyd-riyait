SetFactory("OpenCASCADE");

// Outer rectangle
Rectangle(1) = {0, 0, 0, Pi/2, Pi/2, 0}; // Azimuthal x Polar


// Assign physical regions
Physical Surface("solid_angle") = {1};
