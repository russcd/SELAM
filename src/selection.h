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

using namespace std;

/* A class representing selection on a single locus specified by the CHROMOSOME, POSITION, TYPE, and 
 * the user-defined SELECTION_COEFFICIENTS. */
class single_locus_selection {
public:
    int chromosome ;                                         /// chromosome it's on
    float position ;                                         /// site on chromosome
    string type ;                                            /// universal selection or population specific
    vector<float> selection_coefficients ;                   /// -> genotype -> fitness
} ;

class mate_choice {
public:
    vector<int> chromosomes;
    vector<float> positions;
    vector< vector<float> > selection_coefficients;
};
/* A class representing epistatic selection specified by a pair of CHROMOSOMES, POSITIONS on the chromosomes, and 
 * the user-defined SELECTION_COEFFICIENTS. */
class epistatic_selection {                                  /// Same organization as single-locus, but vector of 2 b/c comparing two sites (epistasis)
public:
    vector<int> chromosomes ;
    vector<float> positions ;
    vector< vector<float> > selection_coefficients ;        // 3x2 vector of A alleles paired with B alleles (ex. AA & BB, Aa & Bb, etc.)
} ;

class output_parameters {
public:
    int gen;
    int subpop;
    int females;
    int males;
    string file;
};

/* Parameters read in by the external file specified by the user; the class consists of dmi specifications for males
 * and females, DMI_MALE & DMI_FEMALE, single locus specifications for males and females, SL_MALE & SL_FEMALE, and sites
 * to track for males and females, MALE_SITES & FEMALE_SITES. */
class selection_parameters {
    
public:
    
    vector<epistatic_selection> dmi_female ;            // Epistatic and single_locus for females
    vector<single_locus_selection> sl_female ;

    vector<epistatic_selection> dmi_male ;              // Epistatic and single_locus for males
    vector<single_locus_selection> sl_male ;
    
    vector<vector<float> > male_sites ;                 // selected sites for males and females
    vector<vector<float> > female_sites ;

    vector<mate_choice> female_mc;
    
} ;

/* A class representing segments of a chromosome described by a list of ANCESTRY_SWITCHES denoted by a beginning position 
 * and the "state", or ancestral relationship, and a vector of SELECTED_MUTATIONS. */
class ancestry_block {
public:
    list< pair<float,int> > ancestry_switch ;          // Stores beginning sites and the "state," or ancestral population it belongs to
    vector<float> selected_mutations ;                  // Selected mutations on the particular ancestry_block
} ;

/* A lookup table of demography parameters dependent on GENERATION of effectiveness. The class is represented by 
 * sex-specific migration rates MIGRATION_RATES_MALE & MIGRATION_RATES_FEMALE, sex-specific ancestral makeups 
 * ANCESTRAL_RATES_MALE & ANCESTRAL_RATES_FEMALE, and sex-specific population sizes POP_SIZE_MALE &
 * POP_SIZE_FEMALE. */
class demography_parameters {
    
public:
    
    //// UPDATE THIS TO HAVE SEPARATE MALE AND FEMALE RATES
    int generation ;
    int size;
    // map< POPULATION TO,  map< POPULATION FROM, PROPORTION FROM > > -- need for both males and females
    map< int, map<int,float> > migration_rates_male;
    map< int, map<int, float> > migration_rates_female;

    // map< POPLULATION TO, map< ANCESTRAL POPULATION, PROPORTION> > -- need for both males and females
    map<int,map<int,float> > ancestral_rates_male;
    map<int, map<int, float> > ancestral_rates_female;


    /// map<SUBPOPULATION, POP_SIZE> -- need for both males and females 
    map<int,int> pop_size_male ;
    map<int, int> pop_size_female;
    
} ;
