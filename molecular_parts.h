#ifndef MOLECULAR_PARTS_H
#define MOLECULAR_PARTS_H

#include <vector>
#include <string>
#include <utility>
#include <iostream>

#include "eigen3/Eigen/Core"

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "constants.h"
#include "math_helper.h"

enum class ElementLabel { Hydrogen, Nitrogen, Carbon, Oxygen, ELEMENT_MAX };

struct ElementIndex {

    ElementIndex() {}
    ElementIndex(ElementLabel lab, std::size_t i)
        : label(lab), index(i) {}

    ElementLabel label;
    std::size_t index;
};

ElementIndex to_element_index(char c, std::size_t index) {
    switch (c) {
    case 'C':
        return ElementIndex(ElementLabel::Carbon, index);
    case 'O':
        return ElementIndex(ElementLabel::Oxygen, index);
    case 'N':
        return ElementIndex(ElementLabel::Nitrogen, index);
    case 'H':
        return ElementIndex(ElementLabel::Hydrogen, index);
    }
}

ElementLabel to_element(char c) {
    switch (c) {
    case 'C':
        return ElementLabel::Carbon;
    case 'O':
        return ElementLabel::Oxygen;
    case 'N':
        return ElementLabel::Nitrogen;
    case 'H':
        return ElementLabel::Hydrogen;
    }
    return ElementLabel::ELEMENT_MAX;
}

std::map<ElementLabel, std::size_t> max_indices;

void init_max_indices() {
    max_indices[ElementLabel::Oxygen] = 3;
    max_indices[ElementLabel::Nitrogen] = 6;
    max_indices[ElementLabel::Carbon] = 23;
    max_indices[ElementLabel::Hydrogen] = 14;
}

std::size_t get_max_indices(ElementLabel label) {

    return max_indices[label];
}

struct Atom {
    ElementIndex element; // C for carbon, O for oxygen
    //std::size_t index; // C1, O2
    cv::Point3d point_abs; // coordinates in he orthogonal system
    cv::Point3d point_box; // fractional coordinates in the box
    Atom() {}
    Atom (ElementIndex ei, double x, double y, double z) {
        //element = ei.label;
        //index = ei.index;
        element = ei;
        point_box = cv::Point3d(x,y,z);\
    }
};

//////////

struct HBond {

    // read from _geom_hbond_publ_flag in .cif file, section #10
    int donor_ind;
    int atom_H_ind;
    int acceptor_ind; // temp

    // desirable
    int donor_site_label;
    int hydrogen_site_label;
    int acceptor_site_label;

    ElementLabel donor_element = ElementLabel::Nitrogen;
    ElementLabel acceptor_element;
    double distance_hydrogen_acceptor;
    double angle_to_carbon_axis;
    cv::Point3i shift = O3;

};



/////////////////////////

struct Edge {
    //cv::Vec3d arrow;
    cv::Point3d arrow;

    bool in_edge = false;
    std::vector<HBond> hbonds; // This Hbond is always size 1 // ??

    cv::Point3d carbon_axis;

    Edge()
        :
       arrow(cv::Point3d(0,0,0))
    {
    }
    explicit Edge(const cv::Point3d a)
        : arrow(a)
    {
    }

    void set_arrow(const cv::Vec3d& arr) {
        arrow = arr;
    }

    void set_arrow(const cv::Point3d& arr) {
        arrow = arr;
    }

    bool hbonds_not_size_one() const {
        return hbonds.size() != 1;
    }

    bool edge_arrow_equal(const Edge& in_arrow) const {
        return arrow == in_arrow.arrow;
    }
};


struct Molecule {

    Molecule () {}
    Molecule (int ind) {
        index = ind;
        internal = true;
    }


    std::size_t index = 0; // the molecule's unique id index in the structure
    cv::Point3i cell = O3; // the cell in the 3D space defined by the given box
    cv::Point3d centre = cv::Point3d(0,0,0);
    cv::Point3d c = cv::Point3d(0,0,0);
    cv::Point3d carbon_axis = cv::Point3d(0,0,0);

    //std::map<char, std::vector<Atom>> atoms;
    std::vector<Atom> associated_atoms;

    std::vector<std::pair<int, int>> indices_NH;

    bool internal = false; // flag of a molecule in the given box
    double carbon_angle_ver = 0;
    double carbon_angle_hor = 0;
    double oxygen_angle = 0;


    std::vector<cv::Point3d> oxygen_rays; // WHAT ARE OXYGEN RAYS?

    std::string id;

    cv::Point3d unit_box; // the unit box that the molecule lies in
    std::vector<Edge> edges;

    std::map<char, std::vector<Atom>> atoms; // stores all the atoms present in the molecule

    void find_angles() {
        carbon_axis = atoms['C'][0].point_abs - atoms['C'][13].point_abs;
        if (carbon_axis.z < 0 or (carbon_axis.z == 0 and carbon_axis.y < 0)
            or (carbon_axis.z == 0 and carbon_axis.y == 0 and carbon_axis.x < 0)) carbon_axis *= -1; // goes upwards
        //std::cout<<"\ncarbon_axis="<<carbon_axis;
        carbon_angle_ver = angle_positive(cv::Point3d(0,0,1), carbon_axis);
        if (carbon_angle_ver > 0) {// only if carbon_axis can be projected to xy-plane
            cv::Point2d carbon_axis_xy(carbon_axis.x, carbon_axis.y);
            carbon_angle_hor = angle_signed(cv::Point2d(1,0), carbon_axis_xy);
        }
        // Find oxygen_rays
        for (std::size_t i = 0; i < atoms['O'].size(); i++ )
            oxygen_rays.push_back( atoms['O'][i].point_abs - centre );

        Matrix rotation = Eigen::MatrixXd::Identity(3,3);
        double c = cos(carbon_angle_ver * CV_PI / 180);
        rotation *= c;

        // Find the rotation matrix to make the carbon axis vertical
        if (carbon_angle_ver > 0) {
            cv::Point3d a( 0, 0, 0 ); // rotation_axis
            a.x = carbon_axis.y;
            a.y = -carbon_axis.x; // a = rotation_axis is orthogonal to carbon_axis
            a *= 1.0 / norm( a );
            double s = std::sin(carbon_angle_ver * CV_PI / 180 );
            rotation += s * cross_product(a) + (1-c) * tensor_product(a); //std::cout<<"\nr="<<rotation;
        }

        ////

        std::vector<Eigen::VectorXd> oxygen_vectors(3);
        for ( int i = 0; i < 3; i++ )
            oxygen_vectors[i] = rotation * V(oxygen_rays[i]);

        // Find 3 angles with oxygen rays
        int min_index = -1;
        double min_abs_angle = 180;
        std::vector<double> oxygen_angles(3);

        for (int i = 0; i < 3; i++) {
            if (fabs( oxygen_vectors[i][2]) > distance_error)
                std::cout << "\nError in Find_Angles: oxygen_vectors[ i ][2]=" << oxygen_vectors[i][2];
            oxygen_angles[i] = angle_signed(cv::Point2d(1,0), cv::Point2d(oxygen_vectors[i][0], oxygen_vectors[i][1]));
            //std::cout<<"\nray"<<i<<"="<<Point2d( oxygen_rays[ i ][0], oxygen_rays[ i ][1] )<<" angle="<<oxygen_angles[ i ];
            if (min_abs_angle > fabs(oxygen_angles[i])) {
                min_abs_angle = fabs(oxygen_angles[i]);
                min_index = i;
            }
        }
        if (min_index < 0) {
            std::cout << "\nError in Find_Angles: min_index="
                      << min_index << " C_axis=" << carbon_axis;
            //<<" rot_axis="<<a<<" mat="<<rotation;
            for (int i = 0; i < 3; i++)
                std::cout << " ray=" << atoms['O'][i].point_abs - centre << "->" << oxygen_rays[i];
        }
        else
            oxygen_angle = oxygen_angles[min_index];
    }
};


void init_molecular_parts() {
    init_max_indices();
}

////// MOLECULE CALCULATIONS

std::vector<double> carbon_axis_hbond_angle(Edge& edge) {

    std::vector<double> angles;
    cv::Point3d c_axis = edge.carbon_axis;

    std::size_t edge_num = 0;
    for (auto& hbond : edge.hbonds) {
        double angle_rad = dot_product_angle(c_axis, edge.arrow);
        double angle_deg = deg_to_rad(angle_rad);

        if (angle_deg > 90) {
            double corrected_angle = 180 - angle_deg;
            angles.push_back(corrected_angle);
            hbond.angle_to_carbon_axis = corrected_angle;
        }
        else {
            hbond.angle_to_carbon_axis = angle_deg;
            angles.push_back(angle_deg);
        }
        ++edge_num;
    }
    return angles;
}


#endif // MOLECULAR_PARTS_H
