#include <fstream>
#include <iostream>
#include <deque>
#include <set>
#include <string>
#include <iterator>
#include <utility>
#include <algorithm>
#include <limits>
#include <cstdlib>

//Boost
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/graph/graphviz.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graph_traits.hpp>
//#include "boost/filesystem/operations.hpp"
//#include "boost/filesystem/path.hpp"
#include <boost/progress.hpp>
#include <boost/config.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>

//Eigen
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseCholesky>
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;

// Constants
const Point3i O3(0, 0, 0);
const Point3d origin(0, 0, 0);
//const double M_PI=3.1472;
//const bool debug = false;
const double distance_error = 1e-2;
//const double length_NH = 1.00577;
//const double length_HO = 2.62; //1.97436;
const double length_ON = 3.07;
const double max_angle_NHO = 100; // degrees
//const double distance_ON = 2.82822; // min distance between N/O atoms in adjacent molecules
//const double big_constant = 1e+8;
const double tolerance_length = 1e-4;
//const double tolerance_dist_HA = 1e-4;

Point3d V (Vector v) { return Point3d( v(0), v(1), v(2) ); }
Vector V (Point3d v) { Vector vector(3); vector << v.x, v.y, v.z; return vector; }

class Index_Value
{
public:
    int index;
    double value;
    Index_Value () {}
    Index_Value (int i, double v) { index=i; value=v; }
};
bool Decreasing_Values (Index_Value const& p1, Index_Value const& p2) { return p1.value > p2.value; }
bool Increasing_Values (Index_Value const& p1, Index_Value const& p2) { return p1.value < p2.value; }

double Angle_Positive (Point3d v1, Point3d v2) { return acos( v1.dot( v2 ) / ( norm( v1 ) * norm( v2 ) ) ) * 180 / M_PI; }
double Det (Point2d v1, Point2d v2) { return v1.x * v2.y - v1.y * v2.x; }
double Angle_Signed (Point2d v1, Point2d v2)
{
    double dot_product = v1.dot( v2 ) / ( norm( v1 ) * norm( v2 ) );
    double det = Det( v1, v2 );
    double angle = acos( dot_product ) * 180 / M_PI;
    if ( det >= 0 ) return angle; // range[0,pi)
    else return -angle; // range [-pi,0)
}

Matrix Cross_Product (Point3d v)
{
    Matrix m(3,3);
    m <<
    0,    -v.z, v.y,
    v.z,  0,    -v.x,
    -v.y, v.x,  0;
    return m;
}

Matrix Tensor_Product (Point3d v)
{
    Matrix m(3,3);
    m <<
    v.x * v.x, v.x * v.y, v.x * v.z,
    v.y * v.x, v.y * v.y, v.y * v.z,
    v.z * v.x, v.z * v.y, v.z * v.z;
    return m;
}

class Box
{
public:
    double a = 0, b = 0, c = 0; // sides
    double alpha = 0, beta = 0, gamma = 0; //angles
    double delta = 0, nu = 0; // delta is the nu is the vertical angle
    Matrix matrix = Eigen::ArrayXXd::Zero( 3, 3 );
    int n_molecules = 0;
    Box () { matrix << 0, 0, 0, 0, 0, 0, 0, 0, 0; }
    void Print()
    {
        std::cout<<"\nBox: a="<<a<<" b="<<b<<" c="<<c<<" al="<<alpha<<" be="<<beta<<" ga="<<gamma<<" de="<<delta<<" nu="<<nu<<" m="<<n_molecules;
    }
    void Find_Matrix ()
    {
        double al = alpha * M_PI / 180;
        double be = beta * M_PI / 180;
        double ga = gamma * M_PI / 180;
        double cos_delta = ( cos(al) - cos(be) * cos(ga) ) / ( sin(be) * sin(ga) );
        delta = acos( cos_delta );
        matrix(0,0) = a;
        matrix(0,1) = b * cos(ga);
        matrix(1,1) = b * sin(ga);
        matrix(0,2) = c * cos(be);
        matrix(1,2) = c * sin(be) * cos_delta;
        matrix(2,2) = c * sin(be) * sin( delta );
        nu = acos( sin(be) * sin( delta ) ) * 180 / M_PI;
        delta *= 180 / M_PI;
    }
    Point3d Abs_Position (Point3d point) { return V( matrix * V( point ) ); }
};

class Atom
{
public:
    char element; // C for carbon, O for oxygen
    int index; // C1, O2
    cv::Point3d point_abs; // coordinates in the orthogonal system
    cv::Point3d point_box; // fractional coordinates in the box
    Atom() {}
    Atom (char e, int i, double x, double y, double z) { element = e; index = i; point_box = Point3d( x, y, z ); }
    void Print() { std::cout<<"\n"<<element<<index<<" p_a="<<point_abs<<" p_b="<<point_box; }
};

class Hbond
{
public:
    char donor_elem = 'N', acceptor_elem;
    int donor_ind, atom_H_ind, acceptor_ind;
    double distance_HA = 0;
    Point3i shift = O3;
    Hbond () {}
    //pb
    void Print ()
    {
        std::cout<<"\n"<<donor_elem<<donor_ind<<"-H"<<atom_H_ind<<"-"<<acceptor_elem<<acceptor_ind;
        std::cout<<" shift="<<shift<<" d_HA="<<distance_HA;
    }
};

class Edge
{
public:
    Point3d arrow;
    bool in = false;
    std::vector<Hbond> bonds;
    Edge () {}
    Edge (Point3d a) { arrow = a; }
};

std::pair<bool,bool> First_edge_smaller (Edge const& e0, Edge const& e1)
{
    if ( e0.bonds.size() < e1.bonds.size() ) return std::make_pair( true, true );
    if ( e0.bonds.size() > e1.bonds.size() ) return std::make_pair( true, false );
    if ( norm( e0.arrow ) + tolerance_length < norm( e1.arrow ) ) return std::make_pair( true, true );
    if ( norm( e0.arrow ) - tolerance_length > norm( e1.arrow ) ) return std::make_pair( true, false );
    for ( int i = 0; i < e0.bonds.size() and i < e1.bonds.size(); i++ )
    {
        if ( e0.bonds[i].distance_HA + tolerance_length < e1.bonds[i].distance_HA )  return std::make_pair( true, true );
        if ( e0.bonds[i].distance_HA - tolerance_length > e1.bonds[i].distance_HA )  return std::make_pair( true, false );
    }
    //std::cout<<" E";
    return std::make_pair( false, false );
}

struct Compare_Edges
{
    bool operator() (Edge const& e0, Edge const& e1) const
    {
        return First_edge_smaller( e0, e1 ).second;
    }
};

class Molecule
{
public:
    int index = 0; //
    Point3i cell = O3; // the cell in the 3G grid defined by the given box
    Point3d centre = Point3d(0,0,0), c = Point3d(0,0,0), carbon_axis = Point3d(0,0,0);
    std::map< char, std::vector<Atom> > atoms;
    std::vector< std::pair <int, int> > indices_NH;
    bool internal = false; // flag of a molecule in the given box
    double carbon_angle_ver = 0, carbon_angle_hor = 0, oxygen_angle = 0;
    std::vector<Point3d> oxygen_rays;
    
    Molecule () {}
    Molecule (int ind) { index = ind; internal = true; }
    Molecule (Molecule& molecule_old, Point3i cell_shift, Point3d abs_shift, int ind)
    {
        index = ind;
        cell = molecule_old.cell + cell_shift;
        // Copy angles
        carbon_angle_ver = molecule_old.carbon_angle_ver;
        carbon_angle_hor = molecule_old.carbon_angle_hor;
        oxygen_angle = molecule_old.oxygen_angle;
        carbon_axis = molecule_old.carbon_axis;
        // Shift atoms
        centre = molecule_old.centre + abs_shift;
        atoms = molecule_old.atoms;
        indices_NH = molecule_old.indices_NH;
        for ( auto it = atoms.begin(); it != atoms.end(); it++ )
            for ( int i = 0; i < it->second.size(); i++ )
                it->second[i].point_abs += abs_shift;
    }
    
    void Print ()
    {
        std::cout<<"\nv"<<index<<cell<<" centre="<<centre<<internal;
        //std::cout<<" C_ver="<<carbon_angle_ver<<" C_hor="<<carbon_angle_hor<<" O_angle="<<oxygen_angle;
        //for ( int i = 0; i < indices_NH.size(); i++ ) std::cout<<"N"<<indices_NH[i].first<<" <->"<<" H"<<indices_NH[i].second<<std::endl;
    }
    
    void Shift (Point3d shift, Molecule& molecule_new)
    {
        // Copy angles
        molecule_new.carbon_angle_ver = carbon_angle_ver;
        molecule_new.carbon_angle_hor = carbon_angle_hor;
        molecule_new.oxygen_angle = oxygen_angle;
        // Shift atoms
        molecule_new.centre = centre + shift;
        molecule_new.atoms = atoms;
        for ( auto it = molecule_new.atoms.begin(); it != molecule_new.atoms.end(); it++ )
            for ( int i = 0; i < it->second.size(); i++ )
                it->second[i].point_abs += shift;
    }
    
    void Find_Angles ()
    {
        carbon_axis = atoms['C'][0].point_abs - atoms['C'][13].point_abs;
        if ( carbon_axis.z < 0 or ( carbon_axis.z == 0 and carbon_axis.y < 0 )
            or ( carbon_axis.z == 0 and carbon_axis.y == 0 and carbon_axis.x < 0 )) carbon_axis *= -1; // goes upwards
        //std::cout<<"\ncarbon_axis="<<carbon_axis;
        carbon_angle_ver = Angle_Positive( cv::Point3d(0,0,1), carbon_axis );
        if ( carbon_angle_ver > 0 ) // only if carbon_axis can be projected to xy-plane
        {
            Point2d carbon_axis_xy( carbon_axis.x, carbon_axis.y );
            carbon_angle_hor = Angle_Signed( Point2d(1,0), carbon_axis_xy );
        }
        // Find oxygen_rays;
        for ( int i = 0; i < atoms['O'].size(); i++ ) oxygen_rays.push_back( atoms['O'][i].point_abs - centre );
        Matrix rotation = Eigen::MatrixXd::Identity(3,3);
        double c = cos( carbon_angle_ver * M_PI / 180 );
        rotation *= c;
        // Find the rotation matrix to make the carbon axis vertical
        if ( carbon_angle_ver > 0 )
        {
            cv::Point3d a( 0, 0, 0 ); // rotation_axis
            a.x = carbon_axis.y;
            a.y = -carbon_axis.x; // a = rotation_axis is orthogonal to carbon_axis
            a *= 1.0 / norm( a );
            double s = sin( carbon_angle_ver * M_PI / 180 );
            rotation += s * Cross_Product( a ) + (1-c) * Tensor_Product( a ); //std::cout<<"\nr="<<rotation;
        }
        std::vector<Vector> oxygen_vectors( 3 );
        for ( int i = 0; i < 3; i++ ) oxygen_vectors[ i ] = rotation * V( oxygen_rays[ i ] );
        // Find 3 angles with oxygen rays
        int min_index = -1;
        double min_abs_angle = 180;
        std::vector<double> oxygen_angles( 3 );
        for ( int i = 0; i < 3; i++ )
        {
            if ( fabs( oxygen_vectors[ i ][2] ) > distance_error )
                std::cout<<"\nError in Find_Angles: oxygen_vectors[ i ][2]="<<oxygen_vectors[ i ][2];
            oxygen_angles[ i ] = Angle_Signed( Point2d(1,0), Point2d( oxygen_vectors[ i ][0], oxygen_vectors[ i ][1] ) );
            //std::cout<<"\nray"<<i<<"="<<Point2d( oxygen_rays[ i ][0], oxygen_rays[ i ][1] )<<" angle="<<oxygen_angles[ i ];
            if ( min_abs_angle > fabs( oxygen_angles[ i ] ) ) { min_abs_angle = fabs( oxygen_angles[ i ] ); min_index = i; }
        }
        if ( min_index < 0 )
        {
            std::cout<<"\nError in Find_Angles: min_index="<<min_index<<" C_axis="<<carbon_axis; //<<" rot_axis="<<a<<" mat="<<rotation;
            for ( int i = 0; i < 3; i++ ) std::cout<<" ray="<<atoms['O'][i].point_abs - centre<<"->"<<oxygen_rays[ i ];
        }
        else oxygen_angle = oxygen_angles[ min_index ];
    }
    
    /*  bool Find_indices_NH()
     {
     indices_NH.assign( atoms['N'].size(), -2 );
     for ( int i = 0; i < atoms['N'].size(); i++ )
     {
     for ( int j = 0; j < atoms['H'].size(); j++ )
     if ( indices_NH[i] < 0 ) // H neighbour not found yet
     {
     double distance_NH = norm( atoms['N'][i].point_abs - atoms['H'][j].point_abs );
     if ( fabs( distance_NH - length_NH ) < distance_error ) indices_NH[i] = j; //atoms['H'][j].index - 1;
     }
     if ( indices_NH[i] < 0 ) { std::cout<<"\nError in Find_indices_NH: no H neighbour of N"<<atoms['N'][i].index; return false; }
     }
     return true;
     }*/
};

Point3d Cross_Product (Point3d v0, Point3d v1)
{
    Point3d v;
    v.x = v0.y * v1.z - v1.y * v0.z;
    v.y = v0.z * v1.y - v1.z * v0.y;
    v.z = v0.x * v1.y - v1.x * v0.y;
    return v;
}

bool Find_Rotation (Point3d v0, Point3d v1, Matrix& rotation)  //What is the purpose of calculation of rotation matrix w.r.t verties
{
    double angle = Angle_Positive( v0, v1 );
    rotation = Eigen::MatrixXd::Identity(3,3);
    if ( angle == 0 ) return true; // identity rotation
    double c = cos( angle * M_PI / 180 );
    double s = sin( angle * M_PI / 180 );
    rotation *= c;
    if ( angle == M_PI ) return true;
    Point3d axis = Cross_Product( v0, v1 );
    axis *= 1.0 / norm( axis );
    rotation += s * Cross_Product( axis ) + (1-c) * Tensor_Product( axis ); //std::cout<<"\nr="<<rotation;
    return true;
}

struct Compare_Points3i
{
    bool operator() (cv::Point3i const& a, cv::Point3i const& b) const
    {
        if ( a.x < b.x ) return true;
        if ( a.x == b.x and a.y < b.y ) return true;
        if ( a.x == b.x and a.y == b.y ) return ( a.z < b.z );
        return false;
    }
};

// tg: type Graph
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS, Molecule, Edge > Graph;
typedef boost::graph_traits<Graph>::edge_iterator Edge_it;
typedef boost::graph_traits<Graph>::vertex_iterator Vertex_it;
typedef boost::property_map< Graph, boost::vertex_index_t >::type Graph_Map;
typedef boost::graph_traits<Graph>::adjacency_iterator Adjacency_it;
typedef std::multimap< Edge, Graph::edge_descriptor, Compare_Edges > Star;

std::pair<bool,bool> First_star_smaller (Star const& s0, Star const& s1)
{
    std::pair<bool,bool> result;
    for ( auto it0 = s0.begin(), it1 = s1.begin(); it0 != s0.end() and it1 != s1.end(); it0++, it1++ )
    {
        result = First_edge_smaller( it0->first, it1->first );
        if ( result.first ) return result;
    }
    return std::make_pair( false, false );
}

class Neighborhood
{
public:
    Star out, in;
    Neighborhood () {}
    Neighborhood (Star out_, Star in_)
    {
        out = out_; in = in_;
    }
};

std::pair<bool,bool> First_neighborhood_smaller (Neighborhood const& n0, Neighborhood const& n1)
{
    if ( n0.out.size() + n0.in.size() < n1.out.size() + n1.in.size() ) return std::make_pair( true, true ); // n0 has fewer out+in edges than n1
    if ( n0.out.size() + n0.in.size() > n1.out.size() + n1.in.size() ) return std::make_pair( true, false );
    if ( n0.out.size() < n1.out.size() ) return std::make_pair( true, true ); // n0 has fewer out edges than n1
    if ( n0.out.size() > n1.out.size() ) return std::make_pair( true, false );
    std::pair<bool,bool> result;
    result = First_star_smaller( n0.out, n1.out );
    if ( result.first ) return result;
    result = First_star_smaller( n0.in, n1.in );
    if ( result.first ) return result;
    std::cout<<"\nEqual stars:";
    return std::make_pair( false, false );
}

struct Compare_Neighborhoods
{
    bool operator() (Neighborhood const& n0, Neighborhood const& n1) const { return First_neighborhood_smaller( n0, n1 ).second; }
};
typedef std::multimap< Neighborhood, Graph::vertex_descriptor, Compare_Neighborhoods > Decoration;

class Vertex_Point
{
public:
    Adjacency_it iter;
    Point3d point;
    Vertex_Point () {}
    Vertex_Point (Adjacency_it it, Point3d p) { iter = it; point = p; }
};
bool Increasing_Points (Vertex_Point const& v1, Vertex_Point const& v2)
{
    if ( v1.point.x < v2.point.x ) return true;
    if ( v1.point.x > v2.point.x ) return false;
    if ( v1.point.y < v2.point.y ) return true;
    if ( v1.point.y > v2.point.y ) return false;
    if ( v1.point.z < v2.point.z ) return true;
    return false;
}

bool Read_Until_Section (std::fstream& file, int section)
{
    std::string line, col1 = "", col2 = "";
    while ( col1 != "#" or col2 != std::to_string( section ) + "." )
    {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> col1 >> col2;
        
    }
    return true;
}

bool Read_Until_String (std::fstream& file, std::string expected_word)
{
    std::string line, current_word = "";
    while ( current_word != expected_word )
    {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> current_word; // the first word in the line
    }
    return true;
}

bool Read_String_Value (std::fstream& file, std::string s, double& v)
{
    std::string line, col1 = "", col2;
    while ( col1 != s )
    {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> col1 >> col2;
    }
    v = atof( col2.c_str() );
    return true;
}

bool Read_String_Value (std::fstream& file, std::string s, int& v)
{
    double value;
    Read_String_Value( file, s, value );
    v = int( value );
    return true;
}

bool Read_2columns (std::fstream& file, std::string& col1, std::string& col2)
{
    std::string line;
    col1 = ""; col2 = "";
    std::getline( file, line );
    std::istringstream stream( line );
    stream >> col1 >> col2;
    return true;
}

bool Read_idxyz (std::fstream& file, std::string& id, double& x, double& y, double& z)
{
    std::string line, col2, col3, col4, col5, col6;
    id = "";
    if ( !getline( file, line ) ) return false;
    std::istringstream stream( line );
    stream >> id >> col2 >> col3 >> col4;
    if ( id == "" ) return false;
    x = atof( col4.c_str() );
    stream >> col5;
    y = atof( col5.c_str() );
    stream >> col6;
    z = atof( col6.c_str() );
    return true;
}

//Read links between atoms within one molecule
bool Read_Link (std::fstream & file, char & element1, char & element2, int & id1, int & id2)
{
    std::string line, col1, col2;
    if ( !getline( file, line ) ) return false;
    std::istringstream stream( line );
    stream >> col1;
    if ( col1 == "" ) return false; // finished reading this section
    stream >> col2;
    element1 = col1[0]; // the first character in the string is the chemical element
    element2 = col2[0];
    id1 = std::stoi( col1.erase(0,1) ); // extract the integer index after the element character
    id2 = std::stoi( col2.erase(0,1) );
    return true;
}

bool Read_Hbond (std::fstream & file, Hbond& bond)
{
    std::string line, col1, col2, col3;
    if ( !getline( file, line ) ) return false;
    std::istringstream stream( line );
    stream >> col1;
    if ( col1 == "" ) return false; // finished reading this section
    stream >> col2 >> col3;
    bond.acceptor_elem = col3[0];
    if ( bond.acceptor_elem != 'O' ) return true; // only N-H-O bonds are needed
    bond.donor_ind = std::stoi( col1.erase(0,1) );
    bond.atom_H_ind = std::stoi( col2.erase(0,1) );
    bond.acceptor_ind = std::stoi( col3.erase(0,1) );
    stream >> col1 >> col2 >> col3; // distances D - H  H...A   D...A
    bond.distance_HA = std::stod( col2 );
    stream >> col1; // angle D - H...A
    stream >> col1; // shift code
    if ( col1 == "." ) { bond.shift = O3; return true; } // both N, O are in the same box
    bond.shift.x = col1.at( 2 ) - '4' - 1; // code 5 means shift 0
    bond.shift.y = col1.at( 3 ) - '4' - 1;
    bond.shift.z = col1.at( 4 ) - '4' - 1;
    return true;
    
}

bool Read_cif (std::string name, std::map<char,int>& max_indices, Box& box, std::vector<Molecule>& molecules, std::vector<Hbond>& Hbonds)
{
    std::fstream file;
    file.open( name.c_str(), std::ios::in );
    if ( ! file.is_open() )
    {
        std::cout<<"\nFile "<<name<<" not found";
        return false;
    }
    std::string line, id, col1, col2, col3, col4;
    Atom atom;
    char element;
    int ind_molecule = 0, ind_atom = 0;
    double x,y,z;
    // initialise molecule parameters
    std::map<char,int> indices;
    for ( auto i : max_indices ) indices.insert( std::make_pair( i.first, 0 ) );
    
    // Read box parameters
    Read_Until_Section( file, 6 );
    Read_String_Value( file, "_cell_length_a", box.a );
    Read_String_Value( file, "_cell_length_b", box.b );
    Read_String_Value( file, "_cell_length_c", box.c );
    Read_String_Value( file, "_cell_angle_alpha", box.alpha );
    Read_String_Value( file, "_cell_angle_beta", box.beta );
    Read_String_Value( file, "_cell_angle_gamma", box.gamma );
    box.Find_Matrix();
    Read_Until_String( file, "_atom_site_refinement_flags" );
    // Read x,y,z coordinates of all atoms
    while ( Read_idxyz( file, id, x, y, z ) )
    {
        element = id[0];
        ind_molecule = int( indices[ element ] / max_indices[ element ] );
        if ( ind_molecule > molecules.size() )
        { std::cout<<"\nError in Read_cif: i_mol="<<ind_molecule<<" mol_size="<<molecules.size(); return false; }
        if ( ind_molecule == molecules.size() ) { molecules.push_back( Molecule( (int)molecules.size() ) ); }
        ind_atom = ++indices[ element ];
        Atom atom( element, ind_atom, x, y ,z );
        
        atom.point_abs = box.Abs_Position( atom.point_box );
        molecules[ ind_molecule ].atoms[ element ].push_back( atom );
        if ( element == 'C' )
        {
            molecules[ ind_molecule ].centre += atom.point_box;
            int ind = ind_atom - ind_molecule * max_indices[ element ];
            if ( ind == 1 or ind == 14 ) molecules[ ind_molecule ].c += atom.point_box;
        }
    }
    
    // Read NH links between atoms
    Read_Until_String( file, "_geom_bond_publ_flag" );
    char element1, element2;
    int id1, id2;
    while( Read_Link( file, element1, element2, id1, id2 ) )
    {
        if ( element1 == 'N' and element2 == 'H' )
        {
            int index_molecule = ( id1 - 1 ) / max_indices['N'];
            id1 -= index_molecule * max_indices[ element1 ];
            id2 -= index_molecule * max_indices[ element2 ];
            molecules[ index_molecule ].indices_NH.push_back( std::make_pair( id1, id2 ) );
        }
    }
    
    // Read hydrogen bonds between molecules
    Read_Until_String( file, "_geom_hbond_publ_flag" ); Read_Until_String( file, "#D" ); Read_Until_String( file, "#" );
    Hbond bond;
    while( Read_Hbond( file, bond ) )
    {
        if ( bond.acceptor_elem != 'O' ) continue;
        Hbonds.push_back( bond );
        //bond.Print();
    }
    
    box.n_molecules = (int)molecules.size();
    box.Print();
    return true;
}

/*
 bool Linked_ON (Molecule& m0, Molecule& m1, double distance_ON)
 {
 bool linked = false;
 for ( auto O : m0.atoms['O'] )
 for ( auto N : m1.atoms['N'] )
 {
 //if ( fabs( norm( O.point_abs - N.point_abs ) - distance_ON ) < 1e-2 ) { return true; }
 if ( norm( O.point_abs - N.point_abs ) < distance_ON ) { std::cout<<" l="<<norm( O.point_abs - N.point_abs ); linked = true; }
 //else std::cout<<" ON="<<norm( O.point_abs - N.point_abs );
 }
 for ( auto O : m1.atoms['O'] )
 for ( auto N : m0.atoms['N'] )
 {
 //std::cout<<" ON="<<norm( O.point_abs - N.point_abs );
 if ( norm( O.point_abs - N.point_abs ) < distance_ON ) { std::cout<<" l="<<norm( O.point_abs - N.point_abs ); linked = true;  }
 //else std::cout<<" ON="<<norm( O.point_abs - N.point_abs );
 }
 return linked;
 }*/

bool Linked_NHO (Molecule& m0, Molecule& m1,double & d)  // I could not understand the operations within this function.
{
    
    bool linked = false;
    for ( int i = 0; i < m0.atoms['N'].size(); i++ )
    {
        int index_H = m0.indices_NH[ i ].second;
        for ( auto j = 0; j < m1.atoms['O'].size(); j++ )
        {
            Point3d vector_HO = m1.atoms['O'][j].point_abs - m0.atoms['H'][ index_H ].point_abs;
            double distance_HO = norm( vector_HO );
            //std::cout<<"\nN"<<m0.atoms['N'][i].index<<"H"<<index_H+1<<"O"<<m1.atoms['O'][j].index<<"="<<distance_HO;
            double distance_ON = norm( m0.atoms['N'][i].point_abs - m1.atoms['O'][j].point_abs );
            if ( distance_ON > length_ON ) continue; //std::cout<<"\nPotential Error: distance_ON="<<distance_ON;
            double angle_NHO = Angle_Positive( m0.atoms['N'][i].point_abs - m0.atoms['H'][ index_H ].point_abs, vector_HO );
            if ( angle_NHO < max_angle_NHO ) continue; //std::cout<<"\nPotential Error: angle_NHO="<<angle_NHO;
            linked = true;
            break;
        }
    }
    return linked;
}

void Order_Neighbours (Graph& s, Vertex_it v, std::vector<Vertex_Point>& c)
{
    c.clear();
    for( auto neighbours = boost::adjacent_vertices( *v, s ); neighbours.first != neighbours.second; ++neighbours.first )
        c.push_back( Vertex_Point( neighbours.first, s[ *neighbours.first ].centre ) );
    sort( c.begin(), c.end(), Increasing_Points );
    std::cout<<"\nSorted:"; for ( auto n : c ) std::cout<<" "<<n.point;
}

bool Neighbours_Equal (std::vector<Vertex_Point>const& neighbours0, std::vector<Vertex_Point>const& neighbours1) // Is it returns equal neighbours within same molecule structure
{
    if ( neighbours0.size() != neighbours1.size() ) { std::cout<<"\nDifferent numbers of neighbours"; return false; }
    for ( int i = 0; i < neighbours0.size(); i++ )
        if ( norm( neighbours0[i].point - neighbours1[i].point ) > distance_error )
        {
            std::cout<<"\nDifferent: "<<neighbours0[i].point<<"!="<<neighbours1[i].point;
            return false;
        }
    std::cout<<"\nEqual neighbours";
    return true;
}

bool Structures_Equal (Graph& s0, Vertex_it v0, Graph& s1, Vertex_it v1) // How it checks equalities of structure through this condition
{
    if ( norm( s0[ *v0 ].centre - s1[ *v1 ].centre ) > distance_error )
    {
        std::cout<<"\nDifferent: "<<s0[ *v0 ].centre<<s1[ *v1 ].centre; // initial vertices
        return false;
    }
    std::cout<<"\nEqual initial vertices: "<<s0[ *v0 ].centre<<s1[ *v1 ].centre;
    std::vector<Vertex_Point> neighbours0, neighbours1;
    Order_Neighbours( s0, v0, neighbours0 );
    Order_Neighbours( s1, v1, neighbours1 );
    if ( ! Neighbours_Equal( neighbours0, neighbours1 ) ) return false;
    
    for ( int i = 0; i < neighbours0.size(); i++ )
    {
        
    }
    return true;
}

//ped
void Print (Graph::edge_descriptor edge, Graph& graph)
{
    auto v0 = boost::source( edge, graph );
    auto v1 = boost::target( edge, graph );
    std::cout<<"\nv"<<graph[ v0 ].index<<graph[ v0 ].cell<<"->v"<<graph[ v1 ].index<<graph[ v1 ].cell<<" arrow="<<graph[ edge ].arrow<<" l="<<norm( graph[ edge ].arrow )<<" bonds="<<graph[ edge ].bonds.size();
}

//pg
void Print (Graph& graph)
{
    for ( auto vertex_pair = vertices( graph ); vertex_pair.first != vertex_pair.second; ++vertex_pair.first)
    {
        auto v = *vertex_pair.first;
        graph[ v ].Print( );
        /*
        std::cout<<" neighbours: ";
        for ( auto neighbors = boost::adjacent_vertices( v, graph ); neighbors.first != neighbors.second; ++neighbors.first )
        {
            neighb++;
            std::cout<<graph[ *neighbors.first ].index<<",";
        }
        */
        int n_out_edges = 0, n_in_edges = 0;
        std::cout<<" out: ";
        for ( auto out_edges = boost::out_edges( v, graph ); out_edges.first != out_edges.second; ++out_edges.first )
        {
            n_out_edges++;
            auto v = graph[ boost::target( *out_edges.first, graph ) ];
            std::cout<<" v"<<v.index<<v.cell; //<<", a="<<graph[ *out_edges.first ].arrow;
        }
        std::cout<<" in: ";
        for ( auto in_edges = boost::in_edges( v, graph ); in_edges.first != in_edges.second; ++in_edges.first )
        {
            n_in_edges++;
            auto v = graph[ boost::source( *in_edges.first, graph ) ];
            std::cout<<" v"<<v.index<<v.cell; //<<", a="<<graph[ *in_edges.first ].arrow;
        }
        std::cout<<" out="<<n_out_edges<<", in="<<n_in_edges<<", total="<<n_out_edges+n_in_edges;
    }
    for ( auto edge_pair = edges( graph ); edge_pair.first != edge_pair.second; ++edge_pair.first)
    {
        Print( *edge_pair.first, graph );
        auto edge = graph[ *edge_pair.first ];
        //std::cout<<" arrow="<<edge.arrow<<" length="<<norm( edge.arrow );
    }
}

typedef std::map< Point3i, std::map< int, Graph::vertex_descriptor >, Compare_Points3i > Vertex_map;

bool Add_Vertex (Graph& graph, Vertex_map& vertices, Point3i cell, int index, Graph::vertex_descriptor& vertex)
{
    auto f = vertices[ cell ].find( index );
    if ( f != vertices[ cell ].end() ) { vertex = f->second; return false; } // vertex found
    vertex = boost::add_vertex( graph );
    graph[ vertex ].index = index;
    graph[ vertex ].cell = cell;
    vertices[ cell ].insert( std::make_pair( index, vertex) );
    return true;
}

int main()
{
    std::string data_folder = "/Users/kurlin/C++/Datasets/MIF/HOF/molGeom/"; // Vitaliy's folder
    Box box;
    std::map< Point3i, std::vector<Molecule>, Compare_Points3i > molecules; // 1st = box position, 2nd = molecules in the box
    std::vector<Hbond> Hbonds;
    int num_structures = 2; //5688
    std::vector<Graph> structures( num_structures );
    std::multimap< Decoration, int > map_structures;
    std::map<char,int> max_indices;
    // T2 molecule info
    max_indices.insert( std::make_pair( 'O', 3 ) );
    max_indices.insert( std::make_pair( 'N', 6 ) );
    max_indices.insert( std::make_pair( 'C', 23 ) );
    max_indices.insert( std::make_pair( 'H', 14 ) );
    
    for ( int ind_structure = 0; ind_structure < num_structures; ind_structure++ )
    {
        Hbonds.clear();
        molecules.clear();
        std::cout<<"\nStructure "<<ind_structure + 1;
        // st1: stage 1 is to extract data from cif, then find centres and angles in the orthogonal system
        Read_cif( data_folder + "T2_" + std::to_string( ind_structure + 1 ) + "_num_molGeom.cif", max_indices, box, molecules[ O3 ], Hbonds );
        for ( auto& m : molecules[ O3 ] )
        {
            //auto m = &molecules[ O3 ][ i ];
            m.internal = true;
            m.centre *= 1.0 / m.atoms['C'].size();
            m.centre = box.Abs_Position( m.centre );
            m.Find_Angles();
        }
        int num_molecules = (int)molecules[ O3 ].size();
    
        Graph graph;
        Vertex_map map_vertices;
        // st2: stage 2 is to build isolated vertices of molecules from 27 boxes in 3x3 neighbourhood
        //for ( auto s : shifts )
        Point3i s = O3; // only internal molecules
        for ( auto& m : molecules[ s ] )
            {
                //molecules[i].Print(); std::cout<<" m"<<i;
                auto vertex = boost::add_vertex( graph ); // new vertex with empty descriptor (Molecule)
                map_vertices[ s ].insert( std::make_pair( m.index, vertex ) );
                graph[ vertex ] = m;
            }
        
        int ind_mol0, ind_mol1;
        Graph::vertex_descriptor v0, v1;
        std::pair<Graph::edge_descriptor, bool> edge;
        // st3: stage 3 is to convert hydrogen bonds into edges
        for ( auto bond : Hbonds )
        {
            ind_mol0 = ( bond.donor_ind - 1 ) / max_indices[ 'N' ];
            ind_mol1 = ( bond.acceptor_ind - 1 ) / max_indices[ 'O' ];
            bond.Print();
            std::cout<<" edge: "<<ind_mol0<<"->"<<ind_mol1; //<<" cells="<<cells;
            std::vector<Point3i> cells{ O3 };
            if ( bond.shift != O3 ) cells.push_back( O3 - bond.shift );
            for ( auto cell : cells )
            {
                //if ( find( shifts.begin(), shifts.end(), cell + bond.shift ) == shifts.end() ) continue; // shifted cell outside 3x3 neighbourhood
                if ( cell == O3 ) v0 = map_vertices[ cell ][ ind_mol0 ];
                else Add_Vertex( graph, map_vertices, cell, ind_mol0, v0 );
                if ( cell + bond.shift == O3 ) { v1 = map_vertices[ cell + bond.shift ][ ind_mol1 ]; }
                else Add_Vertex( graph, map_vertices, cell + bond.shift, ind_mol1, v1 );
                std::cout<<" "<<cell<<"->"<<cell + bond.shift;
                
                edge = boost::edge( v0, v1, graph );
                if ( ! edge.second ) // the edge didn't exist
                {
                    std::cout<<"new";
                    Point3d arrow = molecules[ O3 ][ ind_mol1 ].centre - molecules[ O3 ][ ind_mol0 ].centre;
                    if ( bond.shift != O3 ) arrow += box.Abs_Position( bond.shift );
                    Edge e( arrow );
                    e.bonds.push_back( bond );
                    boost::add_edge( v0, v1, e, graph );
                }
                else
                {
                    std::cout<<"extra";
                    graph[ edge.first ].bonds.push_back( bond );
                }
            }
            
        }
        //boost::write_graphviz( std::cout, graph );
        Print( graph );
    
        Decoration neighborhoods;
        // st4: stage 4 is to build neighborhoods and stars of our/in edges around each vertex
        for ( auto vertex_pair = vertices( graph ); vertex_pair.first != vertex_pair.second; ++vertex_pair.first)
        {
            auto v = *vertex_pair.first;
            if ( ! graph[ v ].internal ) continue; //else graph[ v ].Print();
            Star star_out, star_in;
            for ( auto out_edges = boost::out_edges( v, graph ); out_edges.first != out_edges.second; ++out_edges.first )
            {
                star_out.insert( std::make_pair( graph[ *out_edges.first ], *out_edges.first ) );
            }
            for ( auto in_edges = boost::in_edges( v, graph ); in_edges.first != in_edges.second; ++in_edges.first )
            {
                graph[ *in_edges.first ].in = true;
                star_in.insert( std::make_pair( graph[ *in_edges.first ], *in_edges.first ) );
            }
            neighborhoods.insert( std::make_pair( Neighborhood( star_out, star_in ), v ) );
        }
        std::cout<<"\nNeighborhoods";
        for ( auto n : neighborhoods )
            graph[ n.second ].Print( );
            
        // st5: compute invariants of the graph
        
        structures[ ind_structure ] = graph;
    }
    
    
    std::cout<<"\n";
    
    return 0;
}