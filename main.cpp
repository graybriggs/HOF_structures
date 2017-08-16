#include <algorithm>
#include <iostream>
#include <utility>
#include <string>
#include <vector>

#include "boost/config.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>


#include "analysis.h"
#include "unit_cell.h"
#include "molecular_parts.h"
#include "file_reading.h"

#include "datastructures.h"

#include "graph_definition.h"
#include "graph_decoration.h"
#include "timer.h"

int main(int argc, char *argv[]) {

    init_molecular_parts();



    std::size_t structure_number = 0;
    std::size_t max_files = 5688; // 5688 in file

    Timer timer;
    timer.start_timer();

    while (structure_number < max_files) {
        VertexMap map_vertices;
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
        std::cout << "Star in edges\n";
        do_analysis(neighborhood, structure_number + 1);
        std::cout << "Star out edges\n";
        //do_analysis(neighborhood[0].out);
        /////////////

        ++structure_number;

        //boost::write_graphviz(std::cout, graphman.graph);
    }

    auto diff = timer.end_timer();
    std::cout << "Time taken to analyse " << max_files << " files: " << std::to_string(diff) << std::endl;

    //box.Find_Matrix();

    // access the edges of the specific vertices?
    //Graph::edge_descriptor oe = *out_edges(v, g).first;

    return 0;
}


// populate all molecules in 0,0,0 box position
// add molecule as graph vertex
