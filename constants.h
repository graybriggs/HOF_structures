#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <opencv2/core/core.hpp>

const cv::Point3i O3(0, 0, 0);
const cv::Point3d origin(0, 0, 0);

// THESE VALUES SHOULD BE READ FROM FILE SO THE PROGRAM DOES NOT NEED RECOMPILING IF
// THEY ARE CHANGED
const double distance_error = 1e-2;
const double length_ON = 3.07;
const double max_angle_NHO = 100; // degrees
const double tolerance_length = 1e-2;
//const double tolerance_dist_HA = 1e-4;


inline float deg_to_rad(float angle) {
    const float deg_rad = 3.14159f * 2.0f / 360.0f;
    return angle * deg_rad;
}

inline float rad_to_deg(float radians) {
    const float rad_deg = 3.14159f / (2.0f / 360.0f);
    return radians * rad_deg;
}

/*
#define RAD2DEG(angle) (180 / (M_PI))  * (angle)
#define DEG2RAD(angle) ((M_PI) / 180)  * (angle)
*/

#endif // CONSTANTS_H
