#ifndef ANALYSIS_H
#define ANALYSIS_H

// neighbourhood analysis
// star analysis


// Neighbourhoods?
void check_num_hbonds(const Star& star) {



}


void do_analysis(const Star& star) {

    check_num_hbonds(star);

}

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
