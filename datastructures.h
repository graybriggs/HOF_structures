#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <map>
#include <vector>
#include <string>


#include "opencv2/core/core.hpp"


#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"

#include "functors.h"
#include "constants.h"
#include "unit_cell.h"
#include "molecular_parts.h"
#include "file_reading.h"


struct DataStructures {

    UnitCell cell;
    std::map<char, int> max_indices;
    //std::vector<Molecule> molecules;
    std::map<cv::Point3i, std::vector<Molecule>, Compare_Points3i> molecules;
    std::vector<HBond> hbonds; // holds the vector of molecules within O3
};


void populate_datastructures(const std::string& file_name, DataStructures& ds) {

    UnitCell cell;
    std::map<cv::Point3i, std::vector<Molecule>, Compare_Points3i> molecules; // 1st = box position, 2nd = molecules in the box
    std::vector<HBond> hbonds;
   // std::vector<Graph> structures(num_structures);
    //std::multimap<Decoration, int> map_structures; // what is the int supposed to represent?

    // T2 molecule info
    // moiety
    // this is essentially saying the amount of these types of molecules exist
    std::map<char, int> max_indices;
    max_indices['O'] = 3;
    max_indices['N'] = 6;
    max_indices['C'] = 23;
    max_indices['H'] = 14;

    read_cif(file_name, max_indices, cell, molecules[O3], hbonds); // all molecules

    ds.cell = cell;


    // initialise central molecules
    for (auto& mol : molecules[O3]) {
        mol.internal = true;
        mol.centre *= 1.0 / mol.atoms['C'].size();
        mol.centre = cell.abs_position(mol.centre);
        mol.find_angles();
    }

    ds.molecules[O3] = molecules[O3];
    ds.hbonds = hbonds;
}

#endif // DATASTRUCTURES_H
