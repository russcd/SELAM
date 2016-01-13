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

class cmd_line {
	
public:
    
    int stats_frequency ;                       /// how often to produce stats estimates
    int number_subpopulations ;                 /// number of island style subpopulations
    float generations ;                         /// generations to run this experiment for
    int seed ;                                  /// random number seed
    float ancestry_block_length ;
    bool track_lengths;
    bool track_freq;
    bool hermaphroditic;

    int garbage_freq ;                          /// frequency of garbage collection in generations
    
    vector<float> chrom_lengths ;          /// total length of chromosomes in morgans
    char* selection_file ;                     /// file that describes sites experiencing selection
    char* demography_file;                     /// file that describes the demography of generations
    char* output_file;
    char* freq_file;
    bool hit_error;
    double male_recomb_scalar;

    void read_cmd_line ( int argc, char *argv[] ) ;
    void error(double type, string line);

} ;

void cmd_line::read_cmd_line ( int argc, char *argv[] ) {
    
    //// DEFAULTS
    stats_frequency = 1 ; 
    generations = 1500 ;                        ////        number of generations to run simulation for
    seed = time(NULL) ;                         ////        if we are running many simultaneous instances of this program THIS NEEDS TO BE SET
    ancestry_block_length = 0.05 ;
    garbage_freq = 1 ;

    selection_file = NULL;
    demography_file = NULL;
    output_file = NULL ; 
    freq_file = NULL;

    track_lengths = false;
    track_freq = false;
    hermaphroditic = false;

    bool hit_error = false;
    male_recomb_scalar = 1.0;
    
    /// DEFAULT CHROMOSOME LENGTHS
    /// THESE ARE APPROXIMATELY MORGANS IN FEMALE MEIOSIS!
    chrom_lengths.push_back(1.0) ;         ////        if specified with -c below, ensure that you say -c #number_chr ## ## ##;

    /// SET BY COMMAND LINE
    for (int i=1; i<argc; i++) {
        if ( strcmp(argv[i],"--garbage") == 0 ) {
            garbage_freq = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--seed") == 0 ) {
            seed = atoi(argv[++i]) ;
        }
        if ( strcmp(argv[i],"-g") == 0 ) {
            generations = atoi( argv[++i] ) ;
        }
        if ( strcmp(argv[i],"--abl") == 0 ) {
            ancestry_block_length = atof(argv[++i]) ;
        }
        if ( strcmp(argv[i],"--sel") == 0 || strcmp(argv[i],"-s") == 0 ) {
            selection_file = argv[++i] ;
        }
        if ( strcmp(argv[i],"--output") == 0 || strcmp(argv[i],"-o") == 0 ) {
            output_file = argv[++i] ;
        }
        if ( strcmp(argv[i], "--dem") == 0 || strcmp(argv[i],"-d") == 0 ) {
            demography_file = argv[++i];
        }
        if (strcmp(argv[i], "--freq") == 0 || strcmp(argv[i], "-f") == 0) {
            freq_file = argv[++i];
        }
        if ( strcmp(argv[i],"-c") == 0 ) {
            chrom_lengths.clear() ;
            int number_chroms = atoi( argv[++i] ) ;
            for ( int l = 0 ; l < number_chroms ; l ++ ) {
                chrom_lengths.push_back(atof(argv[++i])) ;
            }
        }
        if (strcmp(argv[i], "--tf") == 0) {
            track_freq = true;
        }
        if (strcmp(argv[i], "--tl") == 0) {
            track_lengths = true;
        }
        if (strcmp(argv[i], "-h") == 0) {
            hermaphroditic = true;
        }
        if (strcmp(argv[i], "-m") == 0) {
            male_recomb_scalar = atof(argv[++i]);
        }
    }

    if (!demography_file) {
        cerr << endl << "ERROR (demography file): Must specify a demography file." << endl;
    }
}


void cmd_line::error(double type, string line) {
    cerr << endl;

    if (type == 1.0) {
        cerr << "ERROR (demography file): Must specify the 0th generation first";
    }
    else if (type == 1.1) {
        cerr << "ERROR (demography file): Generations specified must be in increasing order:";
    }
    else if (type == 1.2) {
        cerr << "ERROR (demography file): Generations specified in the demography file cannot be repeated:";
    }
    else if (type == 2.0) {
        cerr << "ERROR (demography file): First population specified must be an integer value:";
    }
    else if (type == 3.0) {
        cerr << "ERROR (demography file): Sex specification is incorrect:";
    }
    else if (type == 4.0) {
        cerr << "ERROR (demography file): Must specify population size, not rate of migration:";
    }
    else if (type == 4.1) {
        cerr << "ERROR (demography file): Cannot specify migration to ancestral population:";
    }
    else if (type == 4.2) {
        cerr << "ERROR (demography file): Migration rate must be a percent and less than 1:";
    }
    else if (type == 4.3) {
        cerr << "ERROR (demography file): Must give rates for every specified generation. Amount of generations you have specified:";
    }
    else if (type == 5.0) {
        cerr << "ERROR (demography file): Must specify size for population:";
    }
    else if (type == 6.0) {
        cerr << "ERROR (demography file): Male generation 0 not made completeley of ancestral. Population number:";
    }
    else if (type == 6.1) {
        cerr << "ERROR (demography file): Female generation 0 not made completely of ancestral. Population number:";
    }
    else if (type == 6.2) {
        cerr << "ERROR (selection file): Undefined selection type. Type:";
    }
    else if (type == 6.3) {
        cerr << "ERROR (selection file): Undefined sex choice. Sex:";
    }
    else if (type == 6.4) {
        cerr << "ERROR (selection file): Specified a non-existent chromosome. Chromosome number:";
    }
    else if (type == 6.5) {
        cerr << "ERROR (selection file): Position specified is longer than the chromosome. Position:";
    }
    else if (type == 6.6) {
        cerr << "ERROR (selection file): For single locus selection, you must specify 3 genotypes. Number specified:";
    }
    else if (type == 6.7) {
        cerr << "ERROR (selection file): Selection coefficients must be less than 1. Coefficient specified:";
    }
    else if (type == 6.8) {
        cerr << "ERROR (selection file): For epistatic selection, you must specify 9 genotypes. Number specified:";
    }
    else if (type == 6.9) {
        cerr << "ERROR (selection file): For population selection, you must specify 3 genotypes per subpopulation. Number specified:";
    }
    else if (type == 7.0) {
        cerr << "ERROR (selection file): Mate choice can only apply to females or males. Sex specified:";
    }
    else if (type == 7.1) {
        cerr << "ERROR (selection file): For mate choice, you must specify 9 genotypes. Number specified:";
    }
    cerr << endl << endl << line << endl << endl;
    hit_error = true;
}
