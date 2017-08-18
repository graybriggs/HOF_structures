#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "graph_decoration.h"

#include <algorithm>
#include <iterator>
#include <vector>

// neighbourhood analysis
// star analysis


// edge size ordering

struct StarEdgeSize {

    StarEdgeSize(std::size_t sn, std::size_t s)
        : star_num(sn), sz(s)
    {}

    std::size_t star_num;
    std::size_t sz;
};

bool operator<(const StarEdgeSize& se1, const StarEdgeSize& se2) {
    return se1.sz < se2.sz;
}

bool operator>(const StarEdgeSize& se1, const StarEdgeSize& se2) {
    return se1.sz > se2.sz;
}

bool operator<(const Edge& e1, const Edge& e2) {
    return e1.hbonds.size() < e2.hbonds.size();
}

bool operator>(const Edge& e1, const Edge& e2) {
    return e1.hbonds.size() > e2.hbonds.size();
}

///////////////////////////////
///////////////////////////////
///
/// Current analysis implementations
///
///
/// Compare Neighborhoods
/// * For each star within a neighborhood
/// -> Star OUT + IN, LT;GT
/// -> STAR OUT, LT; GT
/// -> STAR IN, LT; GT
///
/// * Compare edges in stars:
/// -> Edge Bonds size, LT;GT
/// -> vector norm of the arrow + a tolerance LT;GT the next without the tolerance
/// -> Check for the bonds size being strictly GT 1 and == 0
/// -> Compare distance of HA: distance_HA + tolerance
/// -> Check bond size > 1; == 0
/// -> Compare the distances of the Hydrogen Acceptor:
///     dist_HA + tolerance GT;LT dist_HA
/// Compare hbonds angle with their associated Carbon atom
///

// Compare (and order) the sum of the in and out edges on each star within a neighborhood
void compare_inout_neighborhood_edge_size(const Decoration& dec) {

    std::vector<StarEdgeSize> size_data;

    std::size_t num = 0;
    for (auto& n : dec) {
        std::size_t out = n.neighborhood.out.star.size();
        std::size_t in = n.neighborhood.in.star.size();
        size_data.push_back(StarEdgeSize(num, out + in));
        ++num;
    }

    std::stable_sort(std::begin(size_data), std::end(size_data));

    std::cout << "Neighborhood IN-OUT order: ";

    for (auto& s : size_data) {
        std::cout << "(" << s.star_num << "," << s.sz << ") ";
    }
    std::cout << std::endl;

}

// Compare (and order) the outgoing edges for each star in the neighborhood.
void compare_out_neighborhood_edge_size(const Decoration& dec) {

    std::vector<StarEdgeSize> size_data;

    std::size_t num = 0;
    for (auto& n : dec) {
        std::size_t out = n.neighborhood.out.star.size();
        size_data.push_back(StarEdgeSize(num, out));
        ++num;
    }

    std::stable_sort(std::begin(size_data), std::end(size_data));

    std::cout << "Edge bonds amount for OUT edges: ";

    for (auto& s : size_data) {
        std::cout << "(" << s.star_num << "," << s.sz << ") ";
    }
    std::cout << std::endl;

}

// Compare (and order) the incoming edges for each star in the neighborhood.
void compare_in_neighborhood_edge_size(const Decoration& dec) {

    std::vector<StarEdgeSize> size_data;

    std::size_t num = 0;
    for (auto& n : dec) {
        std::size_t in = n.neighborhood.in.star.size();
        size_data.push_back(StarEdgeSize(num, in));
        ++num;
    }

    std::stable_sort(std::begin(size_data), std::end(size_data));

    std::cout << "Edge bonds amount for IN edges: ";

    for (auto& s : size_data) {
        std::cout << "(" << s.star_num << "," << s.sz << ") ";
    }
    std::cout << std::endl;
}

/////////////
/////////////

// Edge bond comparisons

//////////////
//////////////
///


// Comapares the number of input and output edges for a star
bool check_neighborhood_inout_edges(Neighborhood& nh) {

    // compare num in edges to out edges
    std::size_t num_in = nh.in.num_edges();
    std::size_t num_out = nh.out.num_edges();

    if (num_in != num_out)
        std::cout << "Edge numbers different, in: " << num_in << ", out: " << num_out << "\n";
    else
        std::cout << "Edge numbers match, in: " << num_in << ", out: " << num_out << "\n";
}

////////////////////////////////
////////////////////////////////

struct EdgeSize {

    EdgeSize(const std::size_t en, const std::size_t bs)
        : edge_num(en), bond_size(bs)
    {}

    std::size_t edge_num;
    std::size_t bond_size;
};

bool operator<(const EdgeSize& es1, const EdgeSize& es2) {
    return es1.bond_size < es2.bond_size;
}

bool operator>(const EdgeSize& es1, const EdgeSize& es2) {
    return es1.bond_size > es2.bond_size;
}

// Edge bonds size LT;GT

void check_edge_out_bonds_size(const Neighborhood& n) {

    std::vector<EdgeSize> edges;

    std::size_t num = 0;
    for (auto& s : n.out.star) {
        edges.push_back(EdgeSize(num, s.edge.hbonds.size()));
        ++num;
    }

    std::stable_sort(std::begin(edges), std::end(edges));

    for (auto& e : edges) {
        std::cout << "(" << e.edge_num << "," << e.bond_size << ") ";
    }
    std::cout << std::endl;
}

void check_edge_in_bonds_size(const Neighborhood& n) {

    std::vector<EdgeSize> edges;

    std::size_t num = 0;
    for (auto& s : n.in.star) {
        edges.push_back(EdgeSize(num, s.edge.hbonds.size()));
        ++num;
    }

    std::stable_sort(std::begin(edges), std::end(edges));

    for (auto& e : edges) {
        std::cout << "(" << e.edge_num << "," << e.bond_size << ") ";
    }
    std::cout << std::endl;
}

struct EdgeNorm {

    std::size_t edge_num;
    Edge edge;

    EdgeNorm(const std::size_t en, const Edge& e)
        : edge_num(en), edge(e)
    {}
};

bool edge_norm(const EdgeNorm& en0, const EdgeNorm& en1) {

    double tol = constants::tolerance_dist_HA;
    if (cv::norm(en0.edge.arrow) + tol < cv::norm(en1.edge.arrow))
        return true;
    else if (cv::norm(en0.edge.arrow) + tol > cv::norm(en1.edge.arrow))
        return false;
}


void out_edge_norm_plus_tolerance(const Neighborhood& n) {

    std::vector<EdgeNorm> edges;

    std::size_t num = 0;
    for (auto& s : n.out.star) {
        edges.push_back(EdgeNorm(num, s.edge));
        ++num;
    }

    std::stable_sort(std::begin(edges), std::end(edges), edge_norm);

    std::cout << "Edge normals comparisons (within tolerance): ";

    for (auto& e : edges) {
        std::cout << e.edge_num << " ";
    }
    std::cout << std::endl;

}


// a Hbonds is a N-O bond with N being the donator and O being the acceptor.
// usually this bond is kept in a list for size of 1. If the size varies then there
// must be more (or less) of these bonds. This function checks for this variance.
void check_num_hbonds(Star& star) {

    for (auto hb : star.star) {
        //Edge e = hb.second;
        Edge e = hb.edge;
        std::size_t hb_sz = e.hbonds.size();
        if (hb_sz != 1) {
            std::cout << "Hbond anomaly" << std::endl;
            if (e.hbonds.size() > hb_sz) {
                std::cout << "Hbond greater than 1. Size: " << hb_sz << "\n";
            }
            else {
                std::cout << "No Hbond exists \n";
            }
        }
    }
    // hbonds normal
}

void check_edge_averages(Star& star) {
    Edge edge;
    if (!star.arrow_average(edge))
        std::cout << "Unable to calculate arrow average - 0 edges\n";

    std::cout << "Star Average Edge: " << edge.arrow << std::endl;
}

void analyse_star(Star& s) {

    check_num_hbonds(s);
    check_edge_averages(s);
}


void analyse(VertexNeighborhood vn) {

    Neighborhood n = vn.neighborhood;

    check_neighborhood_inout_edges(n);

    std::cout << "--- In star ---\n";
    check_edge_in_bonds_size(n);
    analyse_star(n.in);
    std::cout << "--- Out star ---\n";
    check_edge_out_bonds_size(n);
    analyse_star(n.out);
    out_edge_norm_plus_tolerance(n);
}

void do_analysis(Decoration& decoration, std::size_t structure_num) {

    std::cout << "\nAnalysing structure number: " << structure_num << "\n\n";

    std::cout << "Analysing maximum " << decoration.size() << " stars...\n";

    compare_inout_neighborhood_edge_size(decoration);
    compare_out_neighborhood_edge_size(decoration);
    compare_in_neighborhood_edge_size(decoration);



    std::size_t star_num = 0;
    for (auto& d : decoration)  {
        std::cout << "Star num: " << star_num + 1 << "\n";
        analyse(d);
        ++star_num;
    }
}


// --- star sizes
//







// This is interesting dead code found in the original source
/*
void Order_Neighbours (Graph& s, Vertex_it v, std::vector<Vertex_Point>& c) {
    c.clear();
    for( auto neighbours = boost::adjacent_vertices( *v, s ); neighbours.first != neighbours.second; ++neighbours.first )
        c.push_back( Vertex_Point( neighbours.first, s[ *neighbours.first ].centre ) );
    sort( c.begin(), c.end(), Increasing_Points );
    std::cout<<"\nSorted:"; for ( auto n : c ) std::cout<<" "<<n.point;
}

bool Neighbours_Equal (std::vector<Vertex_Point>const& neighbours0, std::vector<Vertex_Point>const& neighbours1) { // Is it returns equal neighbours within same molecule structure
    if ( neighbours0.size() != neighbours1.size() ) { std::cout<<"\nDifferent numbers of neighbours"; return false; }
    for (std::size_t i = 0; i < neighbours0.size(); i++ )
        if ( norm( neighbours0[i].point - neighbours1[i].point ) > distance_error )
        {
            std::cout<<"\nDifferent: " << neighbours0[i].point << "!=" << neighbours1[i].point;
            return false;
        }
    std::cout<<"\nEqual neighbours";
    return true;
}

bool Structures_Equal(Graph& s0, Vertex_it v0, Graph& s1, Vertex_it v1) // How it checks equalities of structure through this condition
{
    if (norm(s0[*v0].centre - s1[*v1].centre) > distance_error) {
        std::cout << "\nDifferent: " << s0[*v0].centre << s1[*v1].centre; // initial vertices
        return false;
    }
    std::cout<<"\nEqual initial vertices: " << s0[*v0].centre << s1[*v1].centre;
    std::vector<Vertex_Point> neighbours0, neighbours1;
    Order_Neighbours(s0, v0, neighbours0);
    Order_Neighbours(s1, v1, neighbours1);
    if (!Neighbours_Equal( neighbours0, neighbours1))
        return false;

    for (std::size_t i = 0; i < neighbours0.size(); i++) {

    }
    return true;
}
*/

#endif // ANALYSIS_H
