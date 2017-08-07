#include <algorithm>
#include <iostream>
#include <utility>
#include <string>
#include <vector>


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>


#include "unit_cell.h"
#include "molecular_parts.h"
#include "file_reading.h"

#include "datastructures.h"

#include "graph_definition.h"
#include "graph_decoration.h"

int main(int argc, char *argv[]) {

    init_molecular_parts();

    VertexMap map_vertices;

    std::size_t structure_number = 0;
    std::size_t max_files = 1; // 5688 in file

    while (structure_number < max_files) {

        DataStructures ds;
        std::string file_name = "/home/gray/Desktop/molGeom/T2_"
                + std::to_string(structure_number + 1) + "_num_molGeom.cif";

        populate_datastructures(file_name, ds);

        GraphManager graphman;
        populate_graph_and_vertexmap(graphman, ds, map_vertices);

        Decoration neighborhood;
        build_neighborhood(graphman.graph, ds, neighborhood);

        /////////////

        // Analysis starts here

        /////////////

        ++structure_number;
    }

    //box.Find_Matrix();

    // access the edges of the specific vertices?
    //Graph::edge_descriptor oe = *out_edges(v, g).first;

    return 0;
}


// populate all molecules in 0,0,0 box position
// add molecule as graph vertex
