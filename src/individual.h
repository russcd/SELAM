#include <string>
#include <iostream>
#include <time.h> 
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <map>
#include <set>
#include <utility>
#include <list>
#include <ctype.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

/** Now include the SELAM package files */

using namespace std;

/* An individual represented by two vectors of pointers to ancestry blocks CHROMOSOMES. This allows a substantial
 * amount less memory to be passed from generation to generation: instead of passing entire lists of ancestry blocks
 * filled with selected sites, individuals pass on pointers to their ancestry blocks. */
class individual {
    
public:
    /// Vector of pointers to an individuals ancestry_blocks, stored as a chromosome
    vector< vector<ancestry_block*> > chromosomes ;
    /// mtDNA is stored as a single variable
    int mtdna ;
    /// print ancestry information for each individual
    void print( vector<float> &chromosome_lengths, ostream &output_stream, bool male, int index, int gen, int subpop ) ;
} ;

void individual::print ( vector<float> &chromosome_lengths, ostream &output_stream, bool male, int index, int gen, int subpop ) {
    
    for ( int c = 0 ; c < chromosomes.size() ; c ++ ) {
        
        if ( chromosome_lengths.at(c/2) == 0 ) {
            continue ;
        }
    
        float current_start = chromosomes.at(c).at(0)->ancestry_switch.front().first ;
        int current_track = chromosomes.at(c).at(0)->ancestry_switch.front().second ;
        
        // Iterate between ancestral tracks, print tracks separately based on state
        for ( int a = 0 ; a < chromosomes.at(c).size() ; a ++ ) {
            for ( auto t = chromosomes.at(c).at(a)->ancestry_switch.begin() ; t != chromosomes.at(c).at(a)->ancestry_switch.end() ; t++ ) {
                if ( t->second != current_track ) {
                    output_stream << gen << "\t" << subpop << "\t" << male << "\t" << index << "\t" << c/2 << "\t" << c%2 << "\t" << current_track << "\t" << current_start << "\t" << t->first << endl ;
                    current_start = t->first ;
                    current_track = t->second ;
                }
            }
        }
        
        // Print last track of chromosome
        output_stream << gen << "\t" << subpop << "\t" << male << "\t" << index << "\t" << c/2 << "\t" << c%2 << "\t" << current_track << "\t" << current_start << "\t" << chromosome_lengths.at(c/2) << endl ;
        output_stream << "##" << "\t" << gen << "\t" << subpop << "\t" << male << "\t" << index << "\t" << c/2 << "\t" << c%2 ;
        for (int a = 0; a < chromosomes.at(c).size(); a++) {
            for (int m = 0; m < chromosomes.at(c).at(a)->selected_mutations.size(); m++) {
                output_stream << "\t" << chromosomes.at(c).at(a)->selected_mutations.at(m) ;
            }
        }
        output_stream << endl;
    }
    output_stream << "##" << "\t" << gen << "\t" << subpop << "\t" << male << "\t" << index << "\t" << "mtdna" << "\t" << mtdna << endl ;
}
