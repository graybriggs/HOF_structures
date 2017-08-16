#ifndef MATH_HELPER_H
#define MATH_HELPER_H

#include <cmath>

#include "opencv2/core/core.hpp"

#include "eigen3/Eigen/Core"

const float PI = 3.14159f;

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

cv::Point3d V(Eigen::VectorXd v) {
    return cv::Point3d(v(0), v(1), v(2));
}

Eigen::VectorXd V (cv::Point3d v) {
    Eigen::VectorXd vector(3);
    vector << v.x, v.y, v.z;
    return vector;
}


inline float deg_to_rad(float angle) {
    const float deg_rad = PI * 2.0f / 360.0f;
    return angle * deg_rad;
}

inline float rad_to_deg(float radians) {
    const float rad_deg = PI / (2.0f / 360.0f);
    return radians * rad_deg;
}

/*
#define RAD2DEG(angle) (180 / (M_PI))  * (angle)
#define DEG2RAD(angle) ((M_PI) / 180)  * (angle)
*/


Matrix cross_product(cv::Point3d v) {

    Matrix m(3,3);
    m <<
    0   ,-v.z,  v.y,
    v.z , 0  , -v.x,
    -v.y, v.x,  0;
    return m;
}


Matrix tensor_product(cv::Point3d v) {

    Matrix m(3,3);
    m <<
    v.x * v.x, v.x * v.y, v.x * v.z,
    v.y * v.x, v.y * v.y, v.y * v.z,
    v.z * v.x, v.z * v.y, v.z * v.z;
    return m;
}


double determinant(cv::Point2d v1, cv::Point2d v2) {
    return v1.x * v2.y - v1.y * v2.x;
}

double angle_positive(cv::Point3d v1, cv::Point3d v2) {
    return std::acos( v1.dot( v2 ) / ( cv::norm(v1) * cv::norm( v2 ))) * 180 / CV_PI;
}

double angle_signed(cv::Point2d v1, cv::Point2d v2) {
    double dot_product = v1.dot(v2) / (cv::norm(v1) * cv::norm(v2));
    double det = determinant(v1, v2);
    double angle = std::acos(dot_product) * 180 / CV_PI;
    if (det >= 0)
        return angle; // range[0,pi)
    else
        return -angle; // range [-pi,0)
}


double dot_product_angle(const cv::Vec3f& v1, const cv::Vec3f& v2) {
    /*
    float len1 = std::sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
    float len2 = std::sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);

    float dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; // v1.dot(v2) / norm(v1)

    float a = dot / (len1 * len2);

    return std::acos(a); // 0..PI
    */
    return v1.ddot(v2);
}

double dot_product_angle(const cv::Point3d &v1, const cv::Point3d &v2) {
    float len1 = std::sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
    float len2 = std::sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z);

    float dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; // v1.dot( v2 ) / norm(v1)

    float a = dot / (len1 * len2);

    return std::acos(a); // 0..PI
}

#endif // MATH_HELPER_H
