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

/* A class representing subpopulations constituted by vectors of individuals, MALES & FEMALES, its own selection 
 * parameters SELECTION, and sex-specific cumulative fitnesses, FEMALE_CUMULATIVE_FITNESS & MALE_CUMULATIVE_FITNESS. */
class subpopulation {

public:

    gsl_ran_discrete_t* LT_M;
    gsl_ran_discrete_t* LT_F;

    /// vector of individuals
    vector<individual> males ;          // Vector of males in the subpopulation 
    vector<individual> females ;        // Vector of females in the subpopulation
    int male_size;
    int female_size;
    /// put all selection functions into this; MALE AND FEMALE SPECIFIC?
    selection_parameters selection ;
    /// fitness functions -- need to have male and female specific fitness functions
    vector<individual> new_males;
    vector<individual> new_females;
    map<int, float> num_anc_male;                      // map of new ancestral individuals -- FROM(A*), NUMBER
    map<int, float> num_anc_female;

    void compute_fitness( vector<individual> &pop, vector<vector<float> > &sites_to_track, vector<epistatic_selection> &DMI,
                        vector<single_locus_selection> &sl, float &ancestry_block_length, bool males) ;


    void initialize_gsl(int m, bool male) {
        double M[m]; for(int i = 0; i < m; i++) { M[i] = 1.0; }
        if (male) {
            LT_M = gsl_ran_discrete_preproc(m, M);
        } else {
            LT_F = gsl_ran_discrete_preproc(m, M);
        }
    }
 
} ;

/* Given a particular sex-specific population POP and specific sites SITES_TO_TRACK, computes fitness based on epistatic selection DMI and 
 * single locus selection SL. Ancestry block length ANCESTRY_BLOCK_LENGTH and a list of ancestry blocks RECOMBINANT_ANCESTRY_TRACTS allow
 * one to locate the seleted sites. */
void subpopulation::compute_fitness ( vector<individual> &pop, vector<vector<float> > &sites_to_track, vector<epistatic_selection> &DMI, vector<single_locus_selection> &sl, float &ancestry_block_length, bool males ) {
    
    double fitness[pop.size()] ;
    
    for ( int i = 0 ; i < pop.size() ; i ++ ) {
        
        fitness[i] = 1.0;
        //// begin by counting all selected site genotypes
        map< int, map<float, int> > genotypes ;// map< chrom, map< position, count> >
        
        for ( int s = 0 ; s < sites_to_track.size() ; s ++ ) {      // iterate through sites_to_track
            for ( int c = 0 ; c < 2 ; c ++ ) {                      // to go through both genotypes -- chromosome 1 & 2 
                /// if it's an x chromosome at the very end, push back two of the same genotype
                if ( s * 2 + c == pop.at(i).chromosomes.size() ) {
                    for ( int m = 0 ; m < sites_to_track.at(s).size() ; m ++ ) {
                        genotypes[s][sites_to_track.at(s).at(m)] *= 2 ;
                    }
                    break ;
                }
                /// othwerwise, see if selected sites are in that block and count them
                for ( int m = 0 ; m < sites_to_track.at(s).size() ; m ++ ) {
                    auto it = find( pop.at(i).chromosomes.at(s*2+c).at(floor(sites_to_track.at(s).at(m)/ancestry_block_length))->selected_mutations.begin(), \
                                    pop.at(i).chromosomes.at(s*2+c).at(floor(sites_to_track.at(s).at(m)/ancestry_block_length))->selected_mutations.end(), \
                                    sites_to_track.at(s).at(m) ) ;
                    if ( it != pop.at(i).chromosomes.at(s*2+c).at(floor(sites_to_track.at(s).at(m)/ancestry_block_length))->selected_mutations.end() ) {
                        genotypes[s][sites_to_track.at(s).at(m)] ++ ;
                    }
                }
            }
        }
        
        //// now compute epistatic fitnesses
        for ( int s = 0 ; s < DMI.size() ; s ++ ) {
            fitness[i] *= DMI.at(s).selection_coefficients.at(genotypes[DMI.at(s).chromosomes.at(0)][DMI.at(s).positions.at(0)]).at(genotypes[DMI.at(s).chromosomes.at(1)][DMI.at(s).positions.at(1)]);
        }
        
        /// now compute single locus fitness effects
        for ( int s = 0 ; s < sl.size() ; s ++ ) {
            fitness[i] *= sl.at(s).selection_coefficients.at(genotypes[sl.at(s).chromosome][sl.at(s).position]) ;
        }
    }
    
    /// now update lookup table based on fitnesses
    if ( males ) {
        gsl_ran_discrete_free(LT_M);
        LT_M = gsl_ran_discrete_preproc(pop.size(), fitness);
    } else {
        gsl_ran_discrete_free(LT_F);
        LT_F = gsl_ran_discrete_preproc(pop.size(),fitness);
    }
}
