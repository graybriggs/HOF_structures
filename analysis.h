#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "graph_decoration.h"

// neighbourhood analysis
// star analysis

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


// a Hbonds is a N-O bond with N being the donator and O being the acceptor.
// usually this bond is kept in a list for size of 1. If the size varies then there
// must be more (or less) of these bonds. This function checks for this variance.
void check_num_hbonds(Star& star) {

    bool is_good = true;

    for (auto hb : star.star) {
        //Edge e = hb.second;
        Edge e = hb.edge;
        std::size_t hb_sz = e.hbonds.size();
        if (hb_sz != 1) {
            std::cout << "Hbond anomaly" << std::endl;
            if (e.hbonds.size() > hb_sz) {
                std::cout << "Hbond greater than 1. Size: " << hb_sz << "\n";
                is_good = false;
            }
            else {
                std::cout << "No Hbond exists \n";
                is_good = false;
            }
        }
    }
    if (is_good)
        std::cout << "HBonds normal." << std::endl;
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

// get all star sizes and sort them into order.
void compare_stars_out_edges(Neighborhood& neighborhood) {

    // = neighborhood.out.num_edges();

}

void compare_stars_in_edges(const Neighborhood& n) {

}

void analyse(VertexNeighborhood vn) {

    Neighborhood n = vn.neighborhood;

    check_neighborhood_inout_edges(n);

    std::cout << "--- In star ---\n";
    analyse_star(n.in);
    std::cout << "--- Out star ---\n";
    analyse_star(n.out);
}

void do_analysis(Decoration& neighbor, std::size_t structure_num) {

    std::cout << "\nAnalysing structure number: " << structure_num << "\n\n";

    std::vector<Neighborhood> n; // hold organised neighborhoods?

    std::cout << "Analysing maximum " << neighbor.size() << " stars...\n";
    std::size_t star_num = 0;
    for (auto& d : neighbor)  {
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
