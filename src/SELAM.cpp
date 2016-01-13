///				 SELAM

/// version 0.31 
/// Copyright (C) 2015 Russell Corbett-Detig and Matt Jones

// This program is free software.
// You can redistribute it and/or modify
// it under the terms of the 
// GNU General Public License.

// This program is distributed without any warranty.

#include <string>
#include <string.h>
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

// Now include the SELAM package files
#include "cmd_line.h"
#include "selection.h"
#include "population.h"

// g++ -std=c++0x -O3 -lgsl -ltcmalloc -o SELAM SELAM_V0.3.cpp

using namespace std;

const gsl_rng *rng;

int main ( int argc, char **argv ) {
    
    clock_t t = clock();

    //// read command line options
    cmd_line options ;
    options.read_cmd_line( argc, argv ) ;
    //initialize rng for gsl lookup table
    
    rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, (long) options.seed);
    // Create population & read in demography.at(0) file
    population pop = population(rng);
    options.hit_error = false;
    
    pop.read_cmd_line(options);
    
    /// check files 
    pop.check_demography_file(options, options.demography_file);
    if ( options.selection_file ) {
        pop.check_selection_file(options, options.selection_file);
    }

    /// read demography and initialize
    pop.read_demography_file(options.demography_file, options);
    pop.initialize_demography(0);
    
    //// store selection
    selection_parameters selection ;
    if (options.selection_file) {
        pop.read_selection_file( options.selection_file) ;
    }
    bool custom_output = false;
    if (options.output_file) {
        pop.read_output_file(options.output_file);
    }
    if (pop.output.size() > 0) {
        custom_output = true;
        pop.curr_output = 0;
    }

    if (options.freq_file) {
        pop.read_freq_file(options.freq_file);
    }
    else {
        pop.default_freq();
    }
    
    //// initialize subpopulations
    pop.initialize_ancestry() ;
        
    if (custom_output) {
        while (pop.output.at(pop.curr_output).gen == 0 && pop.curr_output < pop.output.size())  {
            pop.print_stats_custom();
            pop.curr_output ++ ;
        }
    }
    else {
        pop.print_stats(0);
    }
    
    //// set end point to last output time if output is specified
    if ( pop.output.size() > 0 ) {
        options.generations = pop.output.at(pop.output.size()-2).gen + 1 ;
    }
    
    //// evolve populations
    for ( int g = 1 ; g < options.generations ; g ++ ) {
                
        /// throw out the trash
        if ( g % options.garbage_freq == 0 ) {
            pop.garbage_collect();
        }
        
        pop.compute_fitness();
        pop.update_demography(g);
        pop.migrate();
        pop.select_parents();
        
        // now create new offspring and add them to the appropriate subpopulation
        if (!options.hermaphroditic) {
            pop.create_offspring(pop.total_male_parents, true);
        }
        pop.create_offspring(pop.total_fem_parents, false);
        pop.add_offspring();
        
       if (custom_output) {
            while (pop.output.at(pop.curr_output).gen == g && pop.curr_output < pop.output.size()) {
                pop.print_stats_custom();
                pop.curr_output ++ ;
            }
       }
       else if (g % options.stats_frequency == 0) {
            pop.print_stats(g);
       }
    }

    cout << "SELAM runtime: " << (double) (clock() - t)/1e6 << " seconds" << endl ;
    
    return(0) ;
}
