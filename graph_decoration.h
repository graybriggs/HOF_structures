#ifndef GRAPH_DECORATION_H
#define GRAPH_DECORATION_H

#include "molecular_parts.h"
#include "graph_definition.h"

// link (associate) an actual edge with a vertex descriptor
struct AssociateEdge {
    AssociateEdge(Edge e, Graph::edge_descriptor ed)
        : edge(e), edge_desc(ed) {}
    Edge edge;
    Graph::edge_descriptor edge_desc;
};

class Star {
public:

    Molecule* molecule;
    //Point3d carbon_axis; // initialise this with axis
    //Star_type star;
    std::vector<Edge> edges; // edges are compared as per the predicate in the old code
    Edge average_arrow_vector; // geometrical vector  

    //typedef std::multimap<Edge, Graph::edge_descriptor, Compare_Edges > Star_type;
    //Star_type star;
    std::map<Graph::edge_descriptor, Edge> star;
    //std::vector<AssociateEdge> edge_assoc;

    void set_average_arrow_vector(Edge e) {
        average_arrow_vector = e;
    }

    // std::vector<Edge> edges_out;
    // std::vector<Edge> edges_in;

    Edge arrow_average() {
        cv::Point3d average_vector;

        double x_avg, y_avg, z_avg;
        x_avg = y_avg = z_avg = 0;

        // sum up vector components
        for (auto s = std::begin(edges); s != std::end(edges); ++s) {
            x_avg += s->arrow.x;
            y_avg += s->arrow.y;
            z_avg += s->arrow.z;
        }
        // divide each component by total
        average_vector.x = x_avg / edges.size();
        average_vector.y = y_avg / edges.size();
        average_vector.z = z_avg / edges.size();

        Edge e(average_vector);  // arrow_average
        average_arrow_vector = e;
        //e.print_edge();
        return e;
    }

    // old: First_edge_smaller
    int num_hbonds(const Edge& e0, const Edge& e1) {
        //average_length_edge(e0, e1);
        if (e0.hbonds.size() < e1.hbonds.size()) {
            return 1;
        }
        else if (e0.hbonds.size() > e1.hbonds.size()) {
            return 2;
        }
        else
            return 0;
    }

    int edge_arrow_length(const Edge& e0, const Edge& e1) {
        if (norm(e0.arrow) + tolerance_length < norm(e1.arrow)) {
            return 1;
        }
        if (norm(e0.arrow) - tolerance_length > norm(e1.arrow)) {
            return 2;
        }
        return 0;
    }

    int hbond_per_edge(const Edge& edge) {
        if (edge.hbonds.size() > 1) {
            std::cout << "e0 h-bond size greater than 1" << std::endl;
            //write_to_file("Structure #" + std::to_string(current_mol) + " e0 h-bond size greater than 1\n");
            return 1;
        }
        else if (edge.hbonds.size() == 0) {
            std::cout << "e0 h-bond size == 0" << std::endl;
            //write_to_file("Structure #" + std::to_string(current_mol) + " e0 h-bond size == 0");
            return 0;
        }
    }

    // LOOK AT THIS AGAIN
    int ha_distance(const Edge& e0, const Edge& e1) {
        // compare the distance of the HA
        for (std::size_t i = 0; i < e0.hbonds.size() && i < e1.hbonds.size(); i++ ) {
            if (e0.hbonds[i].distance_hydrogen_acceptor + tolerance_length < e1.hbonds[i].distance_hydrogen_acceptor ) {
                //DBG_PRINT("e0 HA + tolerance < e1 HA")
                return 1;
            }
            if (e0.hbonds[i].distance_hydrogen_acceptor - tolerance_length > e1.hbonds[i].distance_hydrogen_acceptor) {
                //DBG_PRINT("e0 HA + tolerance < e1 HA")
                return 2;
            }
        }
        return 0;
    }

        // compare hbond angles with their associated carbon atom here
    bool angles_within_threshold(const Edge& edge) {
        // function to check if all angles are within some threshhold
        const double angle_threshold = 0.5;

        // this can be converted to an O(log n) complexity routine
        for (auto& a1 : edge.hbonds) {
            for (auto& a2 : edge.hbonds) {
                //FIX THIS!!
                /*
                if (relative_difference(a1.angle_to_c_axis, a2.angle_to_c_axis) < angle_threshold) {
                    return false;
                }
                */
            }
        }
        return true;
    }

};

class Neighborhood {
public:
    Star out, in;
    Neighborhood () {}
    Neighborhood (Star& out_, Star& in_) {
        out = out_;
        in = in_;
    }
};

typedef std::map<Graph::vertex_descriptor, Neighborhood> Decoration;



// st4: stage 4 is to build neighborhoods and stars of out(N->O)/in (O->N) edges around each vertex
void build_neighborhood(Graph& graph, DataStructures& ds, Decoration& decor) {

    std::size_t molecule_number = 0;

    for (auto vertex_pair = vertices(graph); vertex_pair.first != vertex_pair.second; ++vertex_pair.first) {

        auto current_vertex = *vertex_pair.first;
        if (!graph[current_vertex].molecule.internal)
            continue;
        Star star_out, star_in;
        star_out.molecule = &ds.molecules[O3][molecule_number];
        star_in.molecule = &ds.molecules[O3][molecule_number];

        // out_edges(u, g);
        // gets out-going edges from graph g
        typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
        typedef std::pair<EdgeIterator, EdgeIterator> EdgePair;

        for (auto ep = boost::out_edges(current_vertex, graph); ep.first != ep.second; (ep.first)++) {
            Edge e = graph[*ep.first].edge;
            Graph::edge_descriptor edesc = *ep.first;
            //star_out.edge_assoc.push_back(AssociateEdge(e,edesc));
            star_out.star.insert(std::make_pair(edesc, e));
        }

        for (auto ep = boost::in_edges(current_vertex, graph); ep.first != ep.second; (ep.first)++) {
            graph[*ep.first].edge.in_edge = true;
            Graph::edge_descriptor edesc = *ep.first;
            Edge e = graph[*ep.first].edge;
            star_in.star.insert(std::make_pair(edesc, e));
        }
        ++molecule_number;
        decor.insert(std::make_pair(current_vertex, Neighborhood(star_out, star_in)));
    }
}


// output ordered list of vertices + ordered set of egdes for each vertex

#endif // GRAPH_DECORATION_H


