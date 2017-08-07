#ifndef UNIT_CELL_H
#define UNIT_CELL_H

#include <iostream>
#include <cmath>

#include <opencv2/core/core.hpp>

#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseCholesky>
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;


// Point from vector; Vector from point
cv::Point3d V(Eigen::VectorXd v) {
    return cv::Point3d(v(0), v(1), v(2));
}

Eigen::VectorXd V(cv::Point3d v) {
    Eigen::VectorXd vector(3);
    vector << v.x, v.y, v.z;
    return vector;
}



// Box is a representation of the non-rectangular box of the space that the molecule
// under consideration is in.
class UnitCell {
public:
    double a = 0, b = 0, c = 0; // sides
    double alpha = 0, beta = 0, gamma = 0; //angles
    double delta = 0, nu = 0; // delta is the nu is the vertical angle
    Matrix matrix = Eigen::ArrayXXd::Zero(3, 3);
    std::size_t n_molecules = 0;
    UnitCell () {
        matrix << 0, 0, 0, 0, 0, 0, 0, 0, 0;
    }
    void print() {
        std::cout << "\nBox: a=" << a << " b=" << b << " c=" << c
                  << " alpha=" << alpha << " beta=" << beta
                  << " gamma=" << gamma << " de=" << delta << " nu="
                  << nu << " m=" << n_molecules;
    }
    void find_matrix() {
        double al = alpha * M_PI / 180;
        double be = beta * M_PI / 180;
        double ga = gamma * M_PI / 180;
        double cos_delta = (std::cos(al) - std::cos(be) * std::cos(ga)) / (std::sin(be) * std::sin(ga));
        delta = std::acos( cos_delta );
        matrix(0,0) = a;
        matrix(0,1) = b * std::cos(ga);
        matrix(1,1) = b * std::sin(ga);
        matrix(0,2) = c * std::cos(be);
        matrix(1,2) = c * std::sin(be) * cos_delta;
        matrix(2,2) = c * std::sin(be) * std::sin( delta );
        nu = std::acos( sin(be) * std::sin( delta ) ) * 180 / M_PI;
        delta *= 180 / M_PI;
    }
    cv::Point3d abs_position(cv::Point3d point) {
        return V(matrix * V( point) );
    }
};


#endif // UNIT_CELL_H
