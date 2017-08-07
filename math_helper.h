#ifndef MATH_HELPER_H
#define MATH_HELPER_H

#include <cmath>

#include "opencv2/core/core.hpp"

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

Matrix cross_product(cv::Point3d v)
{
    Matrix m(3,3);
    m <<
    0   ,-v.z,  v.y,
    v.z , 0  , -v.x,
    -v.y, v.x,  0;
    return m;
}

Matrix tensor_product(cv::Point3d v)
{
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

#endif // MATH_HELPER_H
