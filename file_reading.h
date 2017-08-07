#ifndef FILE_READING_H
#define FILE_READING_H

#include <fstream>
#include <iostream>
#include <utility>
#include <sstream>

#include "molecular_parts.h"

bool Read_Until_Section(std::fstream& file, int section) {
    std::string line, col1 = "", col2 = "";
    while ( col1 != "#" or col2 != std::to_string(section) + "." ) {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> col1 >> col2;

    }
    return true;
}

bool Read_Until_String(std::fstream& file, std::string expected_word) {
    std::string line, current_word = "";
    while ( current_word != expected_word ) {
        if ( ! std::getline( file, line ) ) return false;
        std::istringstream stream( line );
        stream >> current_word; // the first word in the line
    }
    return true;
}

bool Read_String_Value (std::fstream& file, std::string s, double& v) {
    std::string line, col1 = "", col2;
    while (col1 != s) {
        if (!std::getline(file, line)) return false;
        std::istringstream stream(line);
        stream >> col1 >> col2;
    }
    v = atof(col2.c_str());
    return true;
}

bool Read_String_Value (std::fstream& file, std::string s, int& v) {
    double value;
    Read_String_Value( file, s, value );
    v = int(value);
    return true;
}

bool Read_2columns (std::fstream& file, std::string& col1, std::string& col2) {
    std::string line;
    col1 = ""; col2 = "";
    std::getline( file, line );
    std::istringstream stream( line );
    stream >> col1 >> col2;
    return true;
}

bool Read_idxyz (std::fstream& file, std::string& id, double& x, double& y, double& z) {
    std::string line, col2, col3, col4, col5, col6;
    id = "";
    if (!getline(file, line))
        return false;

    std::istringstream stream(line);
    stream >> id >> col2 >> col3 >> col4;

    if (id == "")
        return false;

    x = atof( col4.c_str());
    stream >> col5;
    y = atof( col5.c_str());
    stream >> col6;
    z = atof( col6.c_str());
    return true;
}

std::string get_col_number(std::string& str) {

    std::string data = "";

    for (auto a = std::begin(str) + 1; a != std::end(str); ++a) {
        data += *a;
    }
    return data;
}

//Read links between atoms within one molecule
bool Read_Link(std::fstream& file, char& element1, char& element2, int& id1, int& id2)
{
    std::string line, col1, col2;
    if ( !getline( file, line ) ) return false;
    std::istringstream stream(line);
    stream >> col1;
    if ( col1 == "" ) return false; // finished reading this section
        //std::cout << "--> " << get_col_number(col1) << std::endl;
    stream >> col2;
    element1 = col1[0]; // the first character in the string is the chemical element
    element2 = col2[0];
    //id1 = std::stoi( col1.erase(0,1) ); // extract the integer index after the element character


    std::string str_num1 = get_col_number(col1);
    if (str_num1.size() == 0)
        return false;
    else
        id1 = std::stoi(str_num1);


    std::string str_num2 = get_col_number(col2);
    if (str_num2.size() == 0)
        return false;
    else
        id2 = std::stoi(str_num2);

    return true;
}


bool Read_Hbond(std::fstream& file, HBond& bond)
{
    std::string line, col1, col2, col3;
    if (!getline(file, line))
        return false;

    std::istringstream stream(line);
    stream >> col1;
    if (col1 == "")
        return false; // finished reading this section

    stream >> col2 >> col3;

    //if (col3[0] != 'O')
    //    return true; // only N-H-O bonds are needed

    switch (col3[0]) {
    case 'N':
        bond.acceptor_element = ElementLabel::Nitrogen;
        break;
    case 'H':
        bond.acceptor_element = ElementLabel::Hydrogen;
        break;
    case 'O':
        bond.acceptor_element = ElementLabel::Oxygen;
        break;
    }

    if (bond.acceptor_element != ElementLabel::Oxygen)
        return true;

    // 1
    std::string str_num1 = get_col_number(col1);
    if (str_num1.size() == 0)
        bond.donor_ind = 0;
    else
        bond.donor_ind = std::stoi(str_num1);

    // 2
    std::string str_num2 = get_col_number(col2);
    if (str_num2.size() == 0)
        bond.atom_H_ind = 0;
    else
        bond.atom_H_ind = std::stoi(str_num2);

    // 3
    std::string str_num3 = get_col_number(col3);
    if (str_num3.size() == 0)
        bond.acceptor_ind = 0;
    else
        bond.acceptor_ind = std::stoi(str_num3);


    stream >> col1 >> col2 >> col3; // distances D - H  H...A   D...A
    bond.distance_hydrogen_acceptor = std::stod(col2);
    stream >> col1; // angle D - H...A
    stream >> col1; // shift code
    if (col1 == ".") {
        bond.shift = O3;
        return true;
    } // both N, O are in the same box
    bond.shift.x = col1.at(2) - '4' - 1; // code 5 means shift 0
    bond.shift.y = col1.at(3) - '4' - 1;
    bond.shift.z = col1.at(4) - '4' - 1;
    return true;
}


//bool read_cif(std::string name, Box& box, std::vector<MolecularIndex>& max_indices) {
bool read_cif(std::string name, std::map<char,int>& max_indices, UnitCell& box, std::vector<Molecule>& molecules, std::vector<HBond>& Hbonds) {
    std::fstream file;
    file.open(name.c_str(), std::ios::in);
    if (!file.is_open()) {
        std::cout << "\nFile " << name << " not found";
        return false;
    }
    std::string line, id, col1, col2, col3, col4;
    Atom atom;
    char element;
    std::size_t ind_molecule = 0, ind_atom = 0;
    double x,y,z;
    // initialise molecule parameters


    std::map<char,int> indices;
    for ( auto i : max_indices )
        indices.insert(std::make_pair(i.first, 0));



    // Read box parameters
    Read_Until_Section( file, 6 );
    Read_String_Value( file, "_cell_length_a", box.a );
    Read_String_Value( file, "_cell_length_b", box.b );
    Read_String_Value( file, "_cell_length_c", box.c );
    Read_String_Value( file, "_cell_angle_alpha", box.alpha );
    Read_String_Value( file, "_cell_angle_beta", box.beta );
    Read_String_Value( file, "_cell_angle_gamma", box.gamma );
    box.find_matrix();
    Read_Until_String( file, "_atom_site_refinement_flags" );
    // Read x,y,z coordinates of all atoms


    while (Read_idxyz(file, id, x, y, z )) {
        element = id[0];
        ind_molecule = static_cast<std::size_t>(indices[element] / max_indices[element]); // index of the current element
        if (ind_molecule > molecules.size()) {
            std::cout << "\nError in Read_cif: i_mol=" << ind_molecule << " mol_size=" << molecules.size();
            return false;
        }
        if (ind_molecule == molecules.size()) {
            molecules.push_back(Molecule(molecules.size()));
        }
        ind_atom = ++indices[element];

        Atom atom(to_element_index(element, ind_atom), x, y, z);

        // relative atom positions to box
        atom.point_abs = box.abs_position( atom.point_box );
        molecules[ind_molecule].atoms[element].push_back(atom);
        if (element == 'C') {
            molecules[ind_molecule].centre += atom.point_box;
            int ind = ind_atom - ind_molecule * max_indices[element];
            if (ind == 1 or ind == 14)
                molecules[ind_molecule].c += atom.point_box;
        }
    }

        /*
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
        */



    // Read hydrogen bonds between molecules
    Read_Until_String(file, "_geom_hbond_publ_flag");
    Read_Until_String(file, "#D");
    Read_Until_String(file, "#");
    HBond bond;
    while(Read_Hbond(file, bond)) {
        if (bond.acceptor_element != ElementLabel::Oxygen)
            Hbonds.push_back(bond);
    }

    box.n_molecules = molecules.size();
    //box.Print();


    return true;
}

#endif // FILE_READING_H
