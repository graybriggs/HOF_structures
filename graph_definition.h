#ifndef GRAPH_DEFINITION_H
#define GRAPH_DEFINITION_H

#include <map>
#include <utility>

#include "boost/config.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"

#include "functors.h"
#include "molecular_parts.h"
#include "unit_cell.h"

#include <opencv2/core/core.hpp>

struct VertexProperty {
    Molecule molecule;
    UnitCell cell;
};

struct EdgeProperty {
    //HBond hbonds;

    explicit EdgeProperty(Edge e)
        : edge(e) {}

    Edge edge;
};


struct GraphProperty {

};



typedef boost::adjacency_list<
            boost::vecS,
            boost::vecS,
            boost::bidirectionalS,
            VertexProperty,
            EdgeProperty,
            boost::no_property
        > Graph;

struct EdgePair {
    Graph::vertex_descriptor v0;
    Graph::vertex_descriptor v1;
};

typedef boost::adjacency_list<>::vertex_descriptor VertexDescriptor;
typedef int molecule_index_t;// move this
typedef std::map<cv::Point3i, std::map<molecule_index_t, Graph::vertex_descriptor>, Compare_Points3i> VertexMap;


// if this code fails compilation then check Boost version is > 1.60.
// std::move will need to be used instead of copy assign for, eg., vertex properties assignment

class GraphManager {
public:

    Graph graph;

    bool add_graph_vertex(VertexMap& vert_map, cv::Point3i cell, const std::size_t index, Graph::vertex_descriptor& v) {

        auto found_vertex = vert_map[cell].find(index);

        if (found_vertex != vert_map[cell].end()) {
            v = found_vertex->second;
            return false;
        }
        else {
            v = boost::add_vertex(graph);
            graph[v].molecule.index = index;
            graph[v].molecule.cell = cell;
            vert_map[cell].insert(std::make_pair(index, v));
        }
        return true;
    }

    bool prepare_vertices(DataStructures& ds, VertexMap& map_vertices) {
        // prepare vertices
        for (auto& mol : ds.molecules[O3]) {
            auto vertex = boost::add_vertex(graph); // new vertex with empty descriptor (Molecule)
            map_vertices[O3].insert(std::make_pair(mol.index, vertex));
            graph[vertex].molecule = mol;
            graph[vertex].cell = UnitCell();
        }
    }

    bool add_verts_to_map(HBond& bond, DataStructures& ds, VertexMap& map_vertices, EdgePair& edge_pair, int ind_mol_N, int ind_mol_O) {

        std::vector<cv::Point3i> cells{ O3 };

        if (bond.shift != O3)
            cells.push_back(O3 - bond.shift);

        for (auto& cell : cells) {
            if (cell == O3) {
                edge_pair.v0 = map_vertices[O3][ind_mol_N];
            }
            else {
                add_graph_vertex(map_vertices, cell, ind_mol_N, edge_pair.v0);
            }

            if (cell + bond.shift == O3) {
                edge_pair.v1 = map_vertices[O3 + bond.shift][ind_mol_O];
            }
            else {
                add_graph_vertex(map_vertices, cell + bond.shift, ind_mol_O, edge_pair.v1);
            }
        }
    }

    void add_hbond_edge(HBond& bond, EdgePair& edge_pair, DataStructures& ds, int ind_mol_N, int ind_mol_O) {
        std::pair<Graph::edge_descriptor, bool> edge;
        edge = boost::edge(edge_pair.v0, edge_pair.v1, graph);

        if (!edge.second) { // the edge didn't exist

            cv::Point3d arrow = ds.molecules[O3][ind_mol_O].centre - ds.molecules[O3][ind_mol_N].centre;
            if (bond.shift != O3) {
                arrow += ds.cell.abs_position(bond.shift); // vector computed assuming that the molecules are within the cell
            }
            // compute carbon angle for the arrow
            Edge edge(arrow);
            edge.hbonds.push_back(bond);
            // ADD CARBON AXIS HERE: e.carbon_axis
            //e.carbon_axis = CARBON_AXIS;
            //carbon_axis_hbond_angle(e); // edge is now modified here. Could be setup like: e = carbon_axis_hbond_angle(e);
            //boost::add_edge(v0, v1, e, graph);

            boost::add_edge(edge_pair.v0, edge_pair.v1, EdgeProperty(edge), graph);
        }
        else {
            graph[edge.first].edge.hbonds.push_back(bond);
        }
    }
};

class DataStructures;

void populate_graph_and_vertexmap(GraphManager& gm, DataStructures& ds, VertexMap& vert_map) {
    EdgePair edge_pair;
    gm.prepare_vertices(ds, vert_map);

    for (auto bond : ds.hbonds) {

        int ind_mol_N = (bond.donor_ind - 1) / 6; //get_max_indices(ElementLabel::Nitrogen);
        int ind_mol_O = (bond.acceptor_ind - 1) / 3; //get_max_indices(ElementLabel::Oxygen);

        gm.add_verts_to_map(bond, ds, vert_map, edge_pair, ind_mol_N, ind_mol_O);
        gm.add_hbond_edge(bond, edge_pair, ds, ind_mol_N, ind_mol_O);
    }
}

#endif // GRAPH_DEFINITION_H
