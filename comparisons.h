#ifndef COMPARISONS_H
#define COMPARISONS_H

/*
Comparisons:

Star:

Test for a star being larger or smaller than the next

Neighbourhood:

* Are the point3d values "arrow" of an edge equal (For in/out edges on a star)

* Test for a star being smaller, ie having less edges (this probably won't occur but test anyway)
* Test number of out + number of in edges being LT or GT for two stars in neighbourhood
* Ditto above but without the addition. compare n0 out edges to n1 out edges

* Test the average of the arrows on the in star
* Test the average of the arrows on the out star


* Average arrow vectors on a star

*/


struct ComparisonData {


};


void edge_in_out_equal(ComparisonData& cd) {
    
}



#endif // COMPARISONS_H
