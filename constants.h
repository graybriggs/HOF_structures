#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <opencv2/core/core.hpp>

#include <map>

const cv::Point3i O3(0, 0, 0);
const cv::Point3d origin(0, 0, 0);

// THESE VALUES SHOULD BE READ FROM FILE SO THE PROGRAM DOES NOT NEED RECOMPILING IF
// THEY ARE CHANGED
const double distance_error = 1e-2;
const double length_ON = 3.07;
const double max_angle_NHO = 100; // degrees
const double tolerance_length = 1e-2;
//const double tolerance_dist_HA = 1e-4;


#endif // CONSTANTS_H
