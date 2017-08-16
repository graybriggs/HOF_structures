#ifndef GRAPH_DEFINITION_H
#define GRAPH_DEFINITION_H

#include <algorithm>
#include <map>
#include <utility>

#include "boost/config.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"

#include "datastructures.h"
#include "functors.h"
#include "molecular_parts.h"
#include "unit_cell.h"

#include <opencv2/core/core.hpp>

struct VertexProperty {
    Molecule molecule;
    //UnitCell cell;
    cv::Point3i cell;
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

struct EdgeVertices {
    Graph::vertex_descriptor v0;
    Graph::vertex_descriptor v1;

};

typedef boost::adjacency_list<>::vertex_descriptor VertexDescriptor;
typedef int molecule_index_t;// move this
typedef std::map<cv::Point3i, std::map<int, Graph::vertex_descriptor>, Compare_Points3i> VertexMap;  // maps the vartices to a given cell

struct VerticesInCell {
        int index;
        Graph::vertex_descriptor vert_desc;
};

//typedef std::map<cv::Point3i, VerticesInCell> VertexMap;


// if this code fails compilation then check Boost version is > 1.60.
// std::move will need to be used instead of copy assign for, eg., vertex properties assignment

class GraphManager {
public:

    Graph graph;

    void add_graph_vertex(VertexMap& vert_map, const cv::Point3i cell, const int index, Graph::vertex_descriptor& v) {

        auto found_vertex = vert_map[cell].find(index);

        // found
        if (found_vertex != vert_map[cell].end()) {
            v = found_vertex->second;
            return;
        }

        v = boost::add_vertex(graph);
        graph[v].molecule.index = index;
        graph[v].molecule.cell = cell;
        vert_map[cell].insert(std::make_pair(index, v));

    }
    // update more explicit name
    bool prepare_O3_vertices(DataStructures& ds, VertexMap& map_vertices) {

        for (auto& mol : ds.molecules[O3]) {
            auto vertex = boost::add_vertex(graph); // new vertex with empty descriptor (Molecule)
            map_vertices[O3].insert(std::make_pair(mol.index, vertex));
            graph[vertex].molecule = mol;
            //graph[vertex].cell = UnitCell();
            graph[vertex].cell = O3;
        }
    }

    bool build_graph(DataStructures& ds, VertexMap& map_vertices) {

        Graph::vertex_descriptor v0 = 0;
        Graph::vertex_descriptor v1 = 0;

        for (auto bond : ds.hbonds) {

            int ind_mol_N = (bond.donor_ind - 1) / 6; //get_max_indices(ElementLabel::Nitrogen);
            int ind_mol_O = (bond.acceptor_ind - 1) / 3; //get_max_indices(ElementLabel::Oxygen);

            std::vector<cv::Point3i> cells{ O3 };

            if (bond.shift != O3)
                cells.push_back(O3 - bond.shift);

            for (auto& cell : cells) {
                if (cell == O3) {
                    v0 = map_vertices[O3][ind_mol_N];
                }
                else {
                    add_graph_vertex(map_vertices, cell, ind_mol_N, v0);
                }

                if (cell + bond.shift == O3) {
                    v1 = map_vertices[O3][ind_mol_O];
                }
                else {
                    add_graph_vertex(map_vertices, cell + bond.shift, ind_mol_O, v1);
                }
                // DEBUG
                std::cout << cell << ind_mol_N << ", v0: " << v0 << "->" << cell + bond.shift << ind_mol_O << ", v1: " << v1 << std::endl;

                auto edge = boost::edge(v0, v1, graph);

                if (!edge.second) { // the edge didn't exist

                    cv::Point3d arrow = ds.molecules[O3][ind_mol_O].centre - ds.molecules[O3][ind_mol_N].centre;
                    if (bond.shift != O3) {
                        arrow += ds.cell.abs_position(bond.shift); // vector computed assuming that the molecules are within the cell
                    }
                    Edge e(arrow);
                    e.hbonds.push_back(bond);

                    auto ep = EdgeProperty(e);
                    boost::add_edge(v0, v1, ep, graph);
                }
                else {
                    graph[edge.first].edge.hbonds.push_back(bond);
                }
            }
        }
    }
};

class DataStructures;

void populate_graph_and_vertexmap(GraphManager& gm, DataStructures& ds, VertexMap& vert_map) {
    // IMPORTANT: call prepare_O3_vertices first.
    gm.prepare_O3_vertices(ds, vert_map);
    gm.build_graph(ds, vert_map);

    for (auto& m : vert_map) {
        std::cout << m.first << ": ";
        for (auto& s : m.second) {
            std::cout << s.first << ", " << s.second << std::endl;
        }
    }


}

#endif // GRAPH_DEFINITION_H
