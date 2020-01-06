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
#include "individual.h"
#include "subpopulation.h"

using namespace std;
const float EPSILON = .00001;



/* Population stored as a collection of subpopulations POPULAITONS dictated by the demography lookup table
 * DEMOGRAPHY and the current effective demography CURR_DEMOGRAPHY. The population stores parental individuals in 
 * a vector PARENTALS, the lengths of each chromosome in vector CHROMOSOME_LENGTHS, the optimal length of ancestry blocks
 * ANCESTRY_BLOCK_LENGTH, the portions of the chromosomes ANCESTRY_BLOCKS, and the ANCESTRY_BLOCK_NUMBER for garbage 
 * collection. A large vector of SELECTED_SITES serves as a concatenation of all selected sites between males and females
 * in all subpopulations. */
class population {

public:
    /// each subpopulation stored as an entry in this vector
    vector<subpopulation> populations ;
    const gsl_rng *rng;
    int num_sub; 
    int num_anc;
    /// demography functions read in / store
    vector<demography_parameters> demography ;
    int curr_demography;
    vector<output_parameters> output;
    int curr_output;
    /// initialize population ancestry
    int garbage_freq;
    vector<subpopulation> parental ;                 /// store ancestral individuals
    vector<float> ancestry_frequency ;               /// ancestry frequencies
    vector<float> chromosome_lengths;
    float ancestry_block_length;

    list<ancestry_block*> ancestry_blocks ;
    int ancestor_block_number ;                      /// always retain ancestor blocks, so this is the start position for the garbage collector
    map<int, map<int, float> > male_offspring;        // keep track of amount of offspring from each subpop- map<TO, map< FROM, AMOUNT > > 
    map<int, map<int, float> > female_offspring;
    vector< vector<float> > sites_to_track;          // master vector of sites to track on all chromosomes
    map< int, map<float, vector<float> > > ancestry_site_frequencies;         // map < CHROMOSOME, map< POSITION, VECTOR< ANCESTRAL FREQUENCIES> > 
    map<int, map<int, map<int, vector<individual*> > > > total_male_parents;         // a master map of offspring to parents in their own subpopulation
    map<int, map<int, map<int, vector<individual*> > > > total_fem_parents;          // map< TO map< FROM, map< OFFSPRING, *vector<PARENTS> > >
    bool track_lengths; // stats options
    bool track_freq;
    bool hermaphroditic;
    double male_recomb_scalar;

    void check_demography_file(cmd_line &options, string demography_file);
    void check_selection_file(cmd_line &options, string selection_file);
    void read_demography_file( string demography_file, cmd_line &options ) ;
    void read_selection_file( string selection_file ) ;
    void read_output_file(string output_file);
    void read_freq_file(string freq_file);
    void default_freq();
    void check_output_file();
    void compute_fitness() ;
    void initialize_ancestry () ;
    void add_mutations(individual old, individual &new_ind, int anc_type, bool male);
    void read_cmd_line(cmd_line &options);
    void initialize_demography(int g);
    void update_demography(int g);
    void print_demography();
    void migrate();
    void check_migrate();                           // make sure that correct migration numbers were computed
    void select_parents();
    individual* choose_dad(int from, individual* mom);
    void check_parents();                           // only used to debug and make sure subpopulation parent maps were created correctly
    void create_offspring(map<int, map<int, map<int, vector<individual*> > > > &parents, bool males, cmd_line options);
    vector<ancestry_block*> recombine(vector<float> &sites, vector<ancestry_block*> &c1, vector<ancestry_block*> &c2);
    ancestry_block* recombine_block(ancestry_block *b1, ancestry_block *b2, float &site);
    void add_offspring(); 
    void print_stats(int g);
    void print_stats_custom();
    void garbage_collect(); 

    population(const gsl_rng *r) {
    	rng = r; 
    }
} ;

/* Compute the fitness of the entire population. */
void population::compute_fitness () { 
    for (int k = 0; k < populations.size(); k++) {

        populations.at(k).compute_fitness( populations.at(k).females, populations.at(k).selection.female_sites, populations.at(k).selection.dmi_female,
                                            populations.at(k).selection.sl_female, ancestry_block_length, false );
        if (!hermaphroditic) {
            populations.at(k).compute_fitness(populations.at(k).males, populations.at(k).selection.male_sites, populations.at(k).selection.dmi_male, 
                                                populations.at(k).selection.sl_male, ancestry_block_length, true ) ;
        }
    }
}

/* Instantiate the population object with cmd_line OPTIONS. */
void population::read_cmd_line(cmd_line &options) {
    ancestry_block_length = options.ancestry_block_length;
    garbage_freq = options.garbage_freq;
    for (int i = 0; i < options.chrom_lengths.size(); i++) {
        chromosome_lengths.push_back(options.chrom_lengths.at(i));
    }
    curr_demography = 0;
    track_lengths = options.track_lengths;
    track_freq = options.track_freq;
    hermaphroditic = options.hermaphroditic;
    male_recomb_scalar = options.male_recomb_scalar;
}

/* Print the demography vector of population POP to ensure the external file was read in correctly. Used 
 * only for debugging purposes. */
void population::print_demography() {
    for (int i = 0; i < demography.size(); i++) {
        cout << demography.at(i).generation << " " << demography.at(i).pop_size_male[0] << " " << demography.at(i).pop_size_female[0] << " " << 
                demography.at(i).ancestral_rates_male[0][0] << " " << demography.at(i).ancestral_rates_female[0][0] << " " << 
                demography.at(i).ancestral_rates_male[0][1] << " " << demography.at(i).ancestral_rates_female[0][1] << 
                demography.at(i).migration_rates_male[0][0] << " " << demography.at(i).migration_rates_male[0][1] << endl;
    }
}


void population::update_demography(int g) {
    if (curr_demography + 1 < demography.size()) {
        int next_dem_gen = demography.at(curr_demography + 1).generation;
        if (g == next_dem_gen) {
            curr_demography += 1;
        }
    }
    demography.at(curr_demography).size = num_sub;
}

/* Check the input file DEMOGRAPHY_FILE for any formatting errors. */
void population::check_demography_file(cmd_line &options, string demography_file) {
    //read file
    map<int, bool> check_pop;
    ifstream demography_stream(demography_file);
    int num_gen = 0;
    set<int> sub_pop;
    num_sub = 0;
    num_anc = 0;

    string line;
    getline(demography_stream, line);
    stringstream line_stream;
    line_stream << line;
    string curr;
    line_stream >> curr; line_stream >> curr; line_stream >> curr;

    // Check that the first generation is 0
    line_stream >> curr;
    num_gen++;
    if (stoi(curr) != 0) {
        options.error(1.0, curr);
    }

    // Check that generations are in consecutive order and that there are no repeats
    int gen = 0;
    while (line_stream >> curr) {
        if (gen > stoi(curr)) {
            options.error(1.1, curr);
        } 
        else if (gen == stoi(curr)) {
            options.error(1.2, curr);
        }
        gen = stoi(curr);
        num_gen++;
    }

    while(getline(demography_stream, line)) {

        stringstream line_stream;
        string pop1;
        string pop2;
        string sex;
        string rate;

        line_stream << line;
        line_stream >> pop1 >> pop2 >> sex;
        sub_pop.insert(stoi(pop1));
        if (!isdigit(pop1[0])) {
            options.error(2.0, pop1);
        }
        if (stoi(pop1)+1 > num_sub) {
            num_sub = stoi(pop1) + 1;
        }
        if (sex != "M" && sex != "F" && sex != "A") {
            options.error(3.0, sex);
        }

        int check_num = 0;
        while (line_stream >> rate) {
            check_num++;
            if (pop1 == pop2) {
                if (stof(rate) < 1) {
                    options.error(4.0, rate);
                } else {
                    check_pop[stoi(pop1)] = true;
                }
            }
            if (pop1[0] == 'a') {
                options.error(4.1, pop1);
            }
            if (pop1 != pop2) {
                if (stof(rate) > 1) {
                    options.error(4.2, rate);
                }
            }
        }
        if (check_num != num_gen) {
            options.error(4.3, to_string(num_gen));
        }
        if (pop2[0] == 'a') {
            pop2.erase(pop2.begin());
            if (stoi(pop2) + 1 > num_anc) {
                num_anc = stoi(pop2) + 1;
            }
        }
    }
    for (auto it = sub_pop.begin(); it != sub_pop.end(); it++) {
        if (!check_pop[*it]) {
            options.error(5.0, to_string(*it));
        }
    }
}

/* Check the input file SELECTION FILE for any formatting errors. */
void population::check_selection_file(cmd_line &options, string selection_file) {
        
    ifstream selection_stream(selection_file);
    string line;
    
    while (getline(selection_stream, line)) {
        stringstream line_stream;

        line_stream << line;
        string type; string sex;
        line_stream >> type >> sex;

        if (type != "S" && type != "M" && type != "P" && type != "D") {
            options.error(6.2, type);
        }

        if (type == "S") {
            if (sex != "A" && sex != "M" && sex != "F") {
                options.error(6.3, sex);
            }
            string chrom; string pos;
            line_stream >> chrom >> pos;
            if (stoi(chrom) >= chromosome_lengths.size()) {
                options.error(6.4, chrom);
            }
            if (stoi(chrom) < chromosome_lengths.size()) {
                int c = atof(chrom.c_str());
                if (atof(pos.c_str()) > chromosome_lengths.at(c)) {
                    options.error(6.5, pos);
                }
            }
            vector<float> genotypes;
            string gen; 
            while (line_stream >> gen) {
                if (atof(gen.c_str()) > 1) {
                    options.error(6.7, gen);
                }
                genotypes.push_back(atof(gen.c_str()));
            }
            if (genotypes.size() != 3) {
                options.error(6.6, to_string(genotypes.size()));
            }
        } 

        else if (type == "E") {
            if (sex != "A" && sex != "M" && sex != "F") {
                options.error(6.3, sex);
            }
            string c1; string c2; string p1; string p2;
            line_stream >> c1 >> c2 >> p1 >> p2;
            if (stoi(c1) >= chromosome_lengths.size()) {
                options.error(6.4, c1);
            }
            if (stoi(c2) >= chromosome_lengths.size()) {
                options.error(6.4, c2);
            }

            if (stoi(c1) < chromosome_lengths.size()) {
                int c = atof(c1.c_str());
                if (atof(p1.c_str()) > chromosome_lengths.at(c)) {
                    options.error(6.5, p1);
                }
            }
            if (stoi(c2) < chromosome_lengths.size()) {
                int c = atof(c2.c_str());
                if (atof(p2.c_str()) > chromosome_lengths.at(c)) {
                    options.error(6.5, p2);
                }
            }

            string gen;
            vector<float> genotypes;
            while (line_stream >> gen) {
                if (atof(gen.c_str()) > 1) {
                    options.error(6.7, gen);
                }
                genotypes.push_back(atof(gen.c_str()));
            }
            if (genotypes.size() != 9) {
                options.error(6.8, to_string(genotypes.size()));
            }
        }

        else if (type == "P") {
            if (sex != "A" && sex != "M" && sex != "F") {
                options.error(6.3, sex);
            }
            string chrom; string pos;
            line_stream >> chrom >> pos;
            if (atof(chrom.c_str()) >= chromosome_lengths.size()) {
                options.error(6.4, chrom);
            }
            if (atof(chrom.c_str()) < chromosome_lengths.size()) {
                int c = atof(chrom.c_str());
                if (atof(pos.c_str()) > chromosome_lengths.at(c)) {
                    options.error(6.5, pos);
                }
            }

            string gen;
            vector<float> genotypes;
            while (line_stream >> gen) {
                if (atof(gen.c_str()) > 1) {
                    options.error(6.7, gen);
                }
                genotypes.push_back(atof(gen.c_str()));
            }
            if (genotypes.size() != num_sub * 3) {
                options.error(6.9, to_string(genotypes.size()));
            }
        }
        
        else if (type == "M") {
            if (sex != "M" && sex != "F") {
                options.error(7.0, sex);
            }
            string c1; string c2; string p1; string p2;
            line_stream >> c1 >> c2 >> p1 >> p2;
            if (stoi(c1) >= chromosome_lengths.size()) {
                options.error(6.4, c1);
            }
            if (stoi(c2) >= chromosome_lengths.size()) {
                options.error(6.4, c2);
            }

            if (stoi(c1) < chromosome_lengths.size()) {
                int c = atof(c1.c_str());
                if (atof(p1.c_str()) > chromosome_lengths.at(c)) {
                    options.error(6.5, p1);
                }
            }
            if (stoi(c2) < chromosome_lengths.size()) {
                int c = atof(c2.c_str());
                if (atof(p2.c_str()) > chromosome_lengths.at(c)) {
                    options.error(6.5, p2);
                }
            }
            string gen;
            vector<float> genotypes;
            while (line_stream >> gen) {
                if (atof(gen.c_str()) > 1) {
                    options.error(6.7, gen);
                }
                genotypes.push_back(atof(gen.c_str()));
            }
            if (genotypes.size() != 9) {
                options.error(7.1, to_string(genotypes.size()));
            }
        }

        if (options.hit_error) {
            exit(1);
        }
    }
}

void population::read_output_file(string output_file) {
    ifstream output_stream(output_file);
    string line;
    while (getline(output_stream, line)) {
        output_parameters new_output;
        stringstream line_stream;
        line_stream << line;
        line_stream >> new_output.gen >> new_output.subpop >> new_output.females >> new_output.males >> new_output.file;
        string type;
        output.push_back(new_output);
    }

    output_parameters new_output;
    new_output.gen = 100000000;
    output.push_back(new_output);
}


/* Create a look up table of sex specific migration rates and demographic composition of two subpopulations given
 * external file DEMOGRAPHY_FILE. Example line: 0 a0 A .5 .05 0 */
void population::read_demography_file ( string demography_file, cmd_line &options) {
    
    /// read file
    ifstream demography_stream ( demography_file) ;

    string line;
    //// start by reading in generations
    getline( demography_stream, line ) ;
    stringstream line_stream ;
    line_stream << line ;
    string trash ;
    line_stream >> trash ; line_stream >> trash ; line_stream >> trash ;        // First three items are trash
    demography_parameters new_demography_entry ;
    while ( line_stream >> new_demography_entry.generation ) {          // Store demography objects indexed by generation
        demography.push_back( new_demography_entry ) ;                  
    }
    for (int d = 0; d < demography.size(); d++) {
        for (int p = 0; p < num_sub; p++) {
            for (int m = 0; m < num_sub; m++) {
                demography.at(d).migration_rates_male[p][m] = 0;
                demography.at(d).migration_rates_female[p][m] = 0;
            }
            for (int a = 0; a < num_anc; a++) {
                demography.at(d).ancestral_rates_male[p][a] = 0;
                demography.at(d).ancestral_rates_female[p][a] = 0;
            }
        }
    }

    /// now read in all lines relevant to those generations
    while( getline(demography_stream, line) ) {
        
        stringstream line_stream;    
        int pop1 ;
        string pop2 ;
        line_stream << line ;           // Feed in line read to line_stream
        line_stream >> pop1 >> pop2 ;   // pop1 = population to; pop2 = population from
        /// read ancestral line -- iff pop2 is ancestral (e.g. a0/a1)
        if ( pop2[0] == 'a' ) { 
            int g = 0 ;
            string sex;
            line_stream >> sex;
            pop2.erase(pop2.begin()) ;
            if (sex == "M") {
                while ( g < demography.size() ) {       // fill up all generations with information
                    line_stream >> demography.at(g).ancestral_rates_male[pop1][stoi(pop2)] ; // Store proportion of ancestral population into subpopulation
                    g++;
                }
            } else if ( sex == "F") {           // define proportion for females
                while (g < demography.size()) {
                    line_stream >> demography.at(g).ancestral_rates_female[pop1][stoi(pop2)];
                    g++;
                }
            } else {                            // in the case that sex == "A"
                while (g < demography.size()) {
                    float proportion;
                    line_stream >> proportion;
                    demography.at(g).ancestral_rates_male[pop1][stoi(pop2)] = proportion;
                    demography.at(g).ancestral_rates_female[pop1][stoi(pop2)] = proportion;
                    g++;
                }
            }
        }
        /// read size line
        else if ( pop1 == stoi(pop2) ) {
            int g = 0 ;
            string sex;
            line_stream >> sex;
            if (sex == "M") {
                while (g < demography.size()) {
                    line_stream >> demography.at(g).pop_size_male[pop1];
                    g++ ;
                }
            } else if (sex == "F") {
                while (g < demography.size()) {
                    line_stream >> demography.at(g).pop_size_female[pop1];
                    g++;
                }
            } else {
                while (g < demography.size()) {
                    int size;
                    line_stream >> size;

                    demography.at(g).pop_size_female[pop1] = size;
                    demography.at(g).pop_size_male[pop1] = size;

                    g++;
                }
            }
        }
        
        // read migration line -- migration from one population to the other
        else {
            int g = 0;
            string sex;
            line_stream >> sex;
            if (sex == "M") {
                while ( g < demography.size() ) {
                    line_stream >> demography.at(g).migration_rates_male[pop1][stoi(pop2)] ;
                    g++ ;
                }
            } else if (sex == "F") {
                while (g < demography.size()) {
                    line_stream >> demography.at(g).migration_rates_female[pop1][stoi(pop2)];
                    g++;
                }
            } else {        // if sex == "A", then the rates apply to both males and females
                while(g < demography.size()) {
                    float rate;
                    line_stream >> rate;

                    demography.at(g).migration_rates_male[pop1][stoi(pop2)] = rate;
                    demography.at(g).migration_rates_female[pop1][stoi(pop2)] = rate;

                    g++;
                }
            }
        }
    }

    for (int k = 0; k < demography.size(); k++) {
        float keep_m;
        float keep_f;
        for (int j = 0; j < num_sub; j++) {
            keep_m = 1;
            for (int m = 0; m < demography.at(k).migration_rates_male[j].size(); m++) {
                keep_m -= demography.at(k).migration_rates_male[j][m];
            }
            for (int a = 0; a < demography.at(k).ancestral_rates_male[j].size(); a++) {
                keep_m -= demography.at(k).ancestral_rates_male[j][a];
            }
            demography.at(k).migration_rates_male[j][j] = keep_m;
        }
        for (int j = 0; j < num_sub; j++) {
            keep_f = 1;
            for (int m = 0; m < demography.at(k).migration_rates_female[j].size(); m++) {
                keep_f -= demography.at(k).migration_rates_female[j][m];
            }
            for (int a = 0; a < demography.at(k).ancestral_rates_female[j].size(); a++) {
                keep_f -= demography.at(k).ancestral_rates_female[j][a];
            }
            demography.at(k).migration_rates_female[j][j] = keep_f;
        }
    }
    if (!hermaphroditic) {
        for (int p = 0; p < demography.at(0).ancestral_rates_male.size(); p++) {
            float total_m = 0.0;
            for (int i = 0; i < demography.at(0).ancestral_rates_male.at(p).size(); i++) {
                total_m += (float) demography.at(0).ancestral_rates_male.at(p).at(i);
            }
            if (total_m - 1.0 > EPSILON) {
               options.error(6.0, to_string(p));
            }
        }
    }

    for (int p = 0; p < demography.at(0).ancestral_rates_female.size(); p++) {
        float total_f = 0.0;
        for (int i = 0; i < demography.at(0).ancestral_rates_female.at(p).size(); i++) {
            total_f += (float) demography.at(0).ancestral_rates_female.at(p).at(i);
        }

        if (total_f - 1.0 > EPSILON) {
            options.error(6.1, to_string(p));
        }
    }
    if (options.hit_error) {
        exit(1);
    }
}

void population::initialize_demography(int g) {
    populations.resize(num_sub);
    sites_to_track.resize(chromosome_lengths.size());
}

/* Feed in selection parameters to the subpopulations specific vectors from the external file SELECTION_FILE. */
void population::read_selection_file ( string selection_file ) {
    
    /// read file
    ifstream sel_stream ( selection_file.c_str() ) ;
    string line ;
    
    for (int k = 0; k < populations.size(); k++) {
        populations.at(k).selection.male_sites.resize(chromosome_lengths.size());
        populations.at(k).selection.female_sites.resize(chromosome_lengths.size());
    }


    while( getline( sel_stream, line ) ) {
        
        if (line.size() == 0) {
            cerr << "ERROR: CANNOT HAVE EMPTY SELECTION LINE... SKIPPING" << endl;
            continue;
        }

        /// first two columns
        char type ;
        char sex ;
        
        stringstream line_stream ;
        line_stream << line ;
        line_stream >> type >> sex ;
        
        /// DMI's
        if ( type == 'D' ) {
            
            //// read the epistasis
            epistatic_selection new_dmi ;
            new_dmi.chromosomes.resize(2) ;
            new_dmi.positions.resize(2) ;
            line_stream >> new_dmi.chromosomes.at(0) >> new_dmi.chromosomes.at(1) >> new_dmi.positions.at(0) >> new_dmi.positions.at(1);
            vector<float> selection_row(3) ;
            for ( int i = 0 ; i < 3 ; i ++ ) {
                line_stream >> selection_row.at(0) >> selection_row.at(1) >> selection_row.at(2) ;
                new_dmi.selection_coefficients.push_back( selection_row ) ;
            }
            /// store it in appropriate vectors for all subpopulations
            for ( int p = 0 ; p < populations.size() ; p ++ ) {
                if ( sex == 'M' || sex == 'A' ) {
                    populations.at(p).selection.dmi_male.push_back( new_dmi ) ;
                    populations.at(p).selection.male_sites.at(new_dmi.chromosomes.at(0)).push_back(new_dmi.positions.at(0)) ;
                    populations.at(p).selection.male_sites.at(new_dmi.chromosomes.at(1)).push_back(new_dmi.positions.at(1)) ;
                }
                if ( sex == 'F' || sex == 'A' ) {
                    populations.at(p).selection.dmi_female.push_back( new_dmi ) ;
                    populations.at(p).selection.female_sites.at(new_dmi.chromosomes.at(0)).push_back(new_dmi.positions.at(0)) ;
                    populations.at(p).selection.female_sites.at(new_dmi.chromosomes.at(1)).push_back(new_dmi.positions.at(1)) ;
                }

                sites_to_track.at(new_dmi.chromosomes.at(0)).push_back(new_dmi.positions.at(0));
                sites_to_track.at(new_dmi.chromosomes.at(1)).push_back(new_dmi.positions.at(1));
            }
        }

        else if (type == 'M') {
            mate_choice new_mc; new_mc.chromosomes.resize(2); new_mc.positions.resize(2);
            line_stream >> new_mc.chromosomes.at(0) >> new_mc.chromosomes.at(1) >> new_mc.positions.at(0) >> new_mc.positions.at(1);
            vector<float> selection_row(3);
            for (int i = 0; i < 3 ; i++) {
                line_stream >> selection_row.at(0) >> selection_row.at(1) >> selection_row.at(2);
                new_mc.selection_coefficients.push_back(selection_row);
            }
            for (int p = 0; p < populations.size(); p++) {
                populations.at(p).selection.female_mc.push_back(new_mc);
            }
            sites_to_track.at(new_mc.chromosomes.at(0)).push_back(new_mc.positions.at(0));
            sites_to_track.at(new_mc.chromosomes.at(1)).push_back(new_mc.positions.at(1));
        }
        
        else if ( type == 'S' ) {
            
            //// read in single locus information
            single_locus_selection new_sl ;
            float coeff1;
            float coeff2; 
            float coeff3;
            line_stream >> new_sl.chromosome >> new_sl.position ;
            
            /// next three positions should be selection parameters
            line_stream >> coeff1 >> coeff2 >> coeff3;
            new_sl.selection_coefficients.push_back(coeff1);
            new_sl.selection_coefficients.push_back(coeff2);
            new_sl.selection_coefficients.push_back(coeff3);
            
            /// store it in appropriate vectors
            for ( int p = 0 ; p < populations.size() ; p ++ ) {
                if ( sex == 'M' || sex == 'A' ) {
                    populations.at(p).selection.sl_male.push_back( new_sl ) ;
                    populations.at(p).selection.male_sites.at(new_sl.chromosome).push_back(new_sl.position) ;
                }
                if ( sex == 'F' || sex == 'A' ) {
                    populations.at(p).selection.sl_female.push_back( new_sl ) ;
                    populations.at(p).selection.female_sites.at(new_sl.chromosome).push_back(new_sl.position) ;
                }
                sites_to_track.at(new_sl.chromosome).push_back(new_sl.position);
            }
        }

        else if ( type == 'P' ) {
            single_locus_selection new_sl ;
            line_stream >> new_sl.chromosome >> new_sl.position ;
            
            for ( int p = 0 ; p < populations.size() ; p ++ ) {
                
                new_sl.selection_coefficients.clear() ; new_sl.selection_coefficients.resize(3) ;
                line_stream >> new_sl.selection_coefficients.at(0) >> new_sl.selection_coefficients.at(1) >> new_sl.selection_coefficients.at(2) ;
                
                /// if yes, store it in appropriate vectors
                if ( sex == 'M' || sex == 'A' ) {
                    populations.at(p).selection.sl_male.push_back( new_sl ) ;
                    populations.at(p).selection.male_sites.at(new_sl.chromosome).push_back(new_sl.position) ;
                }
                if ( sex == 'F' || sex == 'A' ) {
                    populations.at(p).selection.sl_female.push_back( new_sl ) ;
                    populations.at(p).selection.female_sites.at(new_sl.chromosome).push_back(new_sl.position) ;
                }
                sites_to_track.at(new_sl.chromosome).push_back(new_sl.position);
            }
        }
    }
    /// remove duplicate selection sites and only track one
    for ( int p = 0 ; p < populations.size() ; p ++ ) {
        for ( int c = 0 ; c < populations.at(p).selection.female_sites.size() ; c ++ ) {
            sort( populations.at(p).selection.female_sites.at(c).begin(), populations.at(p).selection.female_sites.at(c).end() );
            populations.at(p).selection.female_sites.at(c).erase( unique( populations.at(p).selection.female_sites.at(c).begin(),
                                                                populations.at(p).selection.female_sites.at(c).end() ),
                                                                populations.at(p).selection.female_sites.at(c).end() ) ;
        }
        for ( int c = 0 ; c < populations.at(p).selection.male_sites.size() ; c ++ ) {
            sort( populations.at(p).selection.male_sites.at(c).begin(), populations.at(p).selection.male_sites.at(c).end() );
            populations.at(p).selection.male_sites.at(c).erase( unique( populations.at(p).selection.male_sites.at(c).begin(),
                                                            populations.at(p).selection.male_sites.at(c).end() ),
                                                            populations.at(p).selection.male_sites.at(c).end() ) ;
        }
        for (int c = 0; c < sites_to_track.size(); c++) {
            sort (sites_to_track.at(c).begin(), sites_to_track.at(c).end());
            sites_to_track.at(c).erase(unique(sites_to_track.at(c).begin(), sites_to_track.at(c).end()), 
                                                sites_to_track.at(c).end());
        }
    }
}

/* Feed in the frequencies of selected sites for each ancestral population. */ 
void population::read_freq_file(string freq_file) {
    ifstream freq_stream ( freq_file);
    string line;
    float chrom;
    float site;
    float freq;
    while (getline( freq_stream, line )) {
        stringstream line_stream ;
        line_stream << line ;
        line_stream >> chrom;
        line_stream >> site;
        for (int s = 0; s < num_anc; s++) {
            line_stream >> freq;
            ancestry_site_frequencies[chrom][site].push_back(freq);
        }
    }
}

void population::default_freq() {
    for (int c = 0; c < sites_to_track.size(); c++) {
        for (int m = 0; m < sites_to_track.at(c).size(); m++) {
            ancestry_site_frequencies[c][sites_to_track.at(c).at(m)].push_back(0.0);
            ancestry_site_frequencies[c][sites_to_track.at(c).at(m)].push_back(1.0);
            for (int a = 2; a < num_anc; a++) {
                ancestry_site_frequencies[c][sites_to_track.at(c).at(m)].push_back(0.0);
            }
        }
    }
}

/* Begin the simulation by populating the population object with the parameters given by the user. */
void population::initialize_ancestry ( ) {

    parental.resize(num_anc);
    for (int i = 0; i < num_anc; i++) {
        individual female;
        female.chromosomes.resize(chromosome_lengths.size() * 2);
        /// create both by iterating through block structure
        for ( int r = 0 ; r < chromosome_lengths.size() ; r ++ ) {
            for ( int a = 0 ; a < chromosome_lengths.at(r) / ancestry_block_length ; a ++ ) {
                /// ancestry 0 block
                ancestry_block* anc_block = new ancestry_block ;
                pair<float, int> ancestry_switch = make_pair( a * ancestry_block_length, i) ;         // Make pair of (start_site, state)
                
		        anc_block->ancestry_switch.push_back( ancestry_switch ) ;
                
                ancestry_blocks.push_back( anc_block ) ;
                female.chromosomes.at(r*2).push_back( ancestry_blocks.back() ) ;        // Push back pointer to new ancestry block twice b/c diploid
                female.chromosomes.at(r*2+1).push_back( ancestry_blocks.back() ) ;
            }
        }
        parental.at(i).females.push_back(female) ;
    }
    
    /// remove second X chromosome if x is present and push each into male population vectors
    if (!hermaphroditic) {
        for (int i = 0; i < parental.size(); i++) {
            individual male;
            male = parental.at(i).females.at(0) ;
            male.chromosomes.pop_back() ;                   // Pop second x b/c males are XY
            parental.at(i).males.push_back(male) ;
        }
    }
    
    /// this stores the position of ancestral  ancestry blocks to prevent their garabage collection
    ancestor_block_number = ancestry_blocks.size() ;
    
    /// now store those individuals
    for ( int i = 0 ; i < populations.size() ; i ++ ) {
        // males and females at their respective population size
        populations.at(i).females.resize(demography.at(0).pop_size_female[i]) ;
        if (!hermaphroditic) {
            populations.at(i).males.resize(demography.at(0).pop_size_male[i]) ;
        }

        //// push proportion of each into population vectors
        int a = 0;
        int fem_ind = 0;
        while (a < num_anc) {
            if (a > 0) {
                fem_ind += demography.at(0).pop_size_female[i] * demography.at(0).ancestral_rates_female[i][a-1];
            }
            for (int f = fem_ind; f < demography.at(0).pop_size_female[i] * demography.at(0).ancestral_rates_female[i][a] + fem_ind; f++) {
                individual new_female ;
                add_mutations(parental.at(a).females.at(0), new_female, a, false);
                populations.at(i).females.at(f) = new_female;
            }
            a++;
        }
        
        /// now store male individuals
        a = 0;
        int male_ind = 0;
        if (!hermaphroditic) {
            while (a < num_anc) {
                if (a > 0) {
                    male_ind += demography.at(0).pop_size_male[i] * demography.at(0).ancestral_rates_male[i][a-1];
                }
                for (int m = male_ind; m < demography.at(0).pop_size_male[i] * demography.at(0).ancestral_rates_male[i][a] + male_ind; m++) {
                    individual new_male ;
                    add_mutations(parental.at(a).males.at(0), new_male, a, true);
                    populations.at(i).males.at(m) = new_male;
                }
                a++;
            }
        }
    }
    
    for (int p = 0; p < num_sub; p++) {
        if (!hermaphroditic) {
            populations.at(p).initialize_gsl(1, true);
            populations.at(p).initialize_gsl(1, false);
        } else {
            populations.at(p).initialize_gsl(1, false);
        }
    }
}

void population::add_mutations(individual old, individual &new_ind, int a, bool male) {
    
    new_ind.chromosomes.resize( old.chromosomes.size() ) ;
    
    for (int r = 0; r < sites_to_track.size(); r++) {
        
        for (int d = 0; d < 2; d++) {
            
            /// x chromosome / male check
            if (r*2+d == old.chromosomes.size() ) {
                continue;
            }
            
            /// now record mutations that are needed
            vector<vector<float> > sites_to_add( old.chromosomes.at(r*2+d).size() ) ;
            sites_to_add.resize( old.chromosomes.at(r*2+d).size() ) ;
            for (int m = 0; m < sites_to_track.at(r).size(); m++) {
                if (gsl_ran_flat(rng, 0, 1) < float(ancestry_site_frequencies[r][sites_to_track.at(r).at(m)].at(a))) {
                    int num = floor(sites_to_track.at(r).at(m) / ancestry_block_length);
                    sites_to_add[num].push_back( sites_to_track.at(r).at(m) ) ;
                }
            }
            
            /// now create new blocks with those mutations
            for ( int block = 0 ; block < sites_to_add.size() ; block ++ ) {
                
                /// block with no mutations are equal to ancestral
                if ( sites_to_add[block].size() == 0 ) {
                    new_ind.chromosomes.at(r*2+d).push_back( old.chromosomes.at(r*2+d).at(block) ) ;
                }
                
                else {
                    
                    ancestry_block *new_block = new ancestry_block ;
                    for (auto s = old.chromosomes.at(r*2+d).at(block)->ancestry_switch.begin(); s != old.chromosomes.at(r*2+d).at(block)->ancestry_switch.end(); s++) {
                        new_block->ancestry_switch.push_back(*s);
                    }
                
                    for ( int m = 0 ; m < sites_to_add.at(block).size() ; m ++ ) {
                        new_block->selected_mutations.push_back(sites_to_add.at(block).at(m)) ;
                    }
                    ancestry_blocks.push_back( new_block ) ;
                    new_ind.chromosomes.at(r*2+d).push_back( ancestry_blocks.back() ) ;
                }
            }
        }
    }
}

/* Compute how many individuals will move between subpopulations and populate parental vecotrs MALE_PARENTS and 
 * FEMALE_PARENTS. Each vector is indexed by the subpopulation and the amount of individuals coming from each 
 * population. */
void population::migrate() {

    for (int i = 0; i < populations.size(); i++) {      // iterate through each population TO
        for (int k = 0; k < demography.at(curr_demography).migration_rates_male[i].size(); k++) {  // iterate through each population FROM
            float mig = demography.at(curr_demography).pop_size_male[i] * demography.at(curr_demography).migration_rates_male[i][k];
            male_offspring[i][k] = mig;
        }
        for (int a = 0; a < demography.at(curr_demography).ancestral_rates_male[i].size(); a++) {
            float anc = demography.at(curr_demography).pop_size_male[i] * demography.at(curr_demography).ancestral_rates_male[i][a];
            populations.at(i).num_anc_male[a] = anc;
        }

        for (int k = 0; k < demography.at(curr_demography).migration_rates_female[i].size(); k++) {  // do the same for females
            float mig = demography.at(curr_demography).pop_size_female[i] * demography.at(curr_demography).migration_rates_female[i][k];
            female_offspring[i][k] = mig;
        }
        for (int a = 0; a < demography.at(curr_demography).ancestral_rates_female[i].size(); a++) {
            float anc = demography.at(curr_demography).pop_size_female[i] * demography.at(curr_demography).ancestral_rates_female[i][a];
            populations.at(i).num_anc_female[a] = anc;
        }
    }
}


/* Given a random number generator GENERATOR, calls subpopulation::select_parents() to feed into CREATE_OFFSPRING. */
void population::select_parents() {
    for (int to = 0; to < populations.size(); to++) {          // TO population
        for (int from = 0; from < populations.size(); from++) {      // FROM population
            total_male_parents[to][from].clear() ;
            total_fem_parents[to][from].clear() ;
            
            if (!hermaphroditic) {
                for ( int i = 0 ; i < male_offspring[to][from] ; i ++ ) {
                    total_male_parents[to][from][i].resize(2) ;
                    total_male_parents[to][from][i].at(1) = (&populations.at(from).females.at(gsl_ran_discrete(rng, populations.at(from).LT_F ))) ;
                    total_male_parents[to][from][i].at(0) = choose_dad(from, total_male_parents[to][from][i].at(1));
                }
            }
                
            for ( int i = 0 ; i < female_offspring[to][from] ; i ++ ) {
                total_fem_parents[to][from][i].resize(2) ;
                total_fem_parents[to][from][i].at(1) = (&populations.at(from).females.at(gsl_ran_discrete(rng, populations.at(from).LT_F ))) ;
                total_fem_parents[to][from][i].at(0) = choose_dad(from, total_fem_parents[to][from][i].at(1));
            }
        }
    }
}

individual* population::choose_dad(int from, individual* mom) {
    bool not_dad = true;
    int dad;

    // no mate choice, dad is just accpted
    if (populations.at(from).selection.female_mc.size() == 0) {
        if (!hermaphroditic) {
            dad = gsl_ran_discrete(rng, populations.at(from).LT_M);
            not_dad = false;
        } else {
            dad = gsl_ran_discrete(rng, populations.at(from).LT_F);
            not_dad = false;
        }
    }

    while (not_dad) {
        if (!hermaphroditic) {
            dad = gsl_ran_discrete(rng, populations.at(from).LT_M);
        } else {
            dad = gsl_ran_discrete(rng, populations.at(from).LT_F);
        }
        for ( int s = 0 ; s < populations.at(from).selection.female_mc.size() ; s ++ ) {      // iterate through all mate choice loci
            int dad_gen = 0 ;
            int mom_gen = 0 ;
            
            for ( int c = 0 ; c < 2 ; c ++ ) {          // iterate through each chromosome
                // First check if the specified chromosome is an x chrom
                if (!hermaphroditic) {
                    if (populations.at(from).selection.female_mc.at(s).chromosomes.at(1)*2+c == populations.at(from).males.at(dad).chromosomes.size()) {
                        dad_gen *= 2;
                    }
                    /// otherwise augment genotype
                    else {
                        auto it_d = find( populations.at(from).males.at(dad).chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(1)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(1)/ancestry_block_length))->selected_mutations.begin(), \
                                        populations.at(from).males.at(dad).chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(1)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(1)/ancestry_block_length))->selected_mutations.end(), \
                                        populations.at(from).selection.female_mc.at(s).positions.at(1)) ;
                        if ( it_d != populations.at(from).males.at(dad).chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(1)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(1)/ancestry_block_length))->selected_mutations.end()) {
                            dad_gen ++;
                        }
                    }
                } else {
                    auto it_d = find( populations.at(from).females.at(dad).chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(1)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(1)/ancestry_block_length))->selected_mutations.begin(), \
                                        populations.at(from).females.at(dad).chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(1)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(1)/ancestry_block_length))->selected_mutations.end(), \
                                        populations.at(from).selection.female_mc.at(s).positions.at(1)) ;
                    if ( it_d != populations.at(from).females.at(dad).chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(1)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(1)/ancestry_block_length))->selected_mutations.end()) {
                            dad_gen ++;
                    }
                }
                
                // Mom doesn't have the x chromosome issue, just check if it has the mate choice allele
                auto it_m = find( mom->chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(0)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(0)/ancestry_block_length))->selected_mutations.begin(), \
                                mom->chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(0)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(0)/ancestry_block_length))->selected_mutations.end(), \
                                populations.at(from).selection.female_mc.at(s).positions.at(0));
                if (it_m != mom->chromosomes.at(populations.at(from).selection.female_mc.at(s).chromosomes.at(0)*2+c).at(floor(populations.at(from).selection.female_mc.at(s).positions.at(0)/ancestry_block_length))->selected_mutations.end()) {
                    mom_gen ++;
                }
            }
            
            // Now compute the probability these two individuals can mate based on the genotype matrix
            float mc_prob = populations.at(from).selection.female_mc.at(s).selection_coefficients.at(mom_gen).at(dad_gen);
            if (gsl_ran_flat(rng, 0, 1) < mc_prob) {
                not_dad = false;
            } else {
                // even if a previous match was made, if this match is not good this cannot be the dad
                not_dad = true;
            }
        }
    }
    // return a pointer to an approved father
    if (!hermaphroditic) {
        return &populations.at(from).males.at(dad);
    }
    return &populations.at(from).females.at(dad);
}


void population::create_offspring(map<int, map<int, map<int, vector<individual*> > > > &parents, bool males, cmd_line options) {
    
    for (int p = 0; p < populations.size(); p++) {          // iterate through each subpopulation; p = to
        for (int k = 0; k < populations.size(); k++) {      // iterate through other subpopulations; k = from
            // Create each offspring from parents
            for (int i = 0; i < parents[p][k].size(); i++) {
                individual offspring;
                vector<individual*> par_vec = parents[p][k][i];
                for (int chrom = 0; chrom < chromosome_lengths.size(); chrom++) {
                    
                    // swap dad chromosome with probability 0.5 for mendelian inheritance
                    if (gsl_ran_flat(rng, 0, 1) < 0.5) {
                        if (chrom*2+1 < par_vec.at(0)->chromosomes.size()) {
                            swap(par_vec.at(0)->chromosomes.at(chrom*2), par_vec.at(0)->chromosomes.at(chrom*2+1));
                        }
                    }

                    // do the same swap with mom chromosome
                    if (gsl_ran_flat(rng, 0, 1) < 0.5) {
                        swap(par_vec.at(1)->chromosomes.at(chrom*2), par_vec.at(1)->chromosomes.at(chrom*2+1));
                    }
                    
                    // determine how many recombination events will take place
                    vector<float> dad_sites( gsl_ran_poisson(rng, (chromosome_lengths.at(chrom) * male_recomb_scalar)) ) ;
                    for (int c = dad_sites.size() - 1 ; c>= 0 ; c-- ) {
                        dad_sites.at(c) = gsl_ran_flat(rng, 0, chromosome_lengths.at(chrom)) ;
                        if ( gsl_ran_flat(rng,0,1) < options.gc_fraction ) {
			    dad_sites.push_back( dad_sites[c] + gsl_ran_exponential(rng, options.gc_rate ) ) ;
			}
                    }

                    vector<float> mom_sites( gsl_ran_poisson(rng, chromosome_lengths.at(chrom)) ) ;
                    for (int c = mom_sites.size() - 1 ; c >= 0 ; c-- ) {
                        mom_sites.at(c) = gsl_ran_flat(rng, 0, chromosome_lengths.at(chrom)) ;
			if ( gsl_ran_flat(rng,0,1) < options.gc_fraction ) {
                            mom_sites.push_back( mom_sites[c] + gsl_ran_exponential(rng, options.gc_rate ) ) ;
                        }
                    }

                    // no recombination for moms
                    if (mom_sites.size() == 0) {
                        offspring.chromosomes.push_back(par_vec.at(1)->chromosomes.at(chrom*2));
                    }
                    else {    // recombination events
                        offspring.chromosomes.push_back(recombine(mom_sites, par_vec.at(1)->chromosomes.at(chrom*2), par_vec.at(1)->chromosomes.at(chrom*2+1)));
                    }

                    
                    // if we're on the dad's x, give it to females only
                    if (chrom*2+1 == par_vec.at(0)->chromosomes.size()) {
                        if (!males) {
                            offspring.chromosomes.push_back(par_vec.at(0)->chromosomes.at(chrom*2));
                        }
                    }
                    // recombination for a non-x chromosomes
                    else if (dad_sites.size() != 0) {
                        offspring.chromosomes.push_back(recombine(dad_sites, par_vec.at(0)->chromosomes.at(chrom*2), par_vec.at(0)->chromosomes.at(chrom*2+1)));
                    // inherit non-recombinant autosome
                    } else {
                        offspring.chromosomes.push_back(par_vec.at(0)->chromosomes.at(chrom*2));
                    }
                }
                if (males) {
                    populations.at(p).new_males.push_back(offspring);
                } else {
                    populations.at(p).new_females.push_back(offspring);    
                }
            }
        }
    }
}

vector<ancestry_block*> population::recombine(vector<float> &sites, vector<ancestry_block*> &c1, vector<ancestry_block*> &c2) {
    sort( sites.begin(), sites.end() ) ;
    vector<ancestry_block*> child_chrom(c1.size());
    //// identify blocks to look at they are already sorted
    vector<int> block_sites ;
    for ( int s = 0 ; s < sites.size() ; s ++ ) {
        block_sites.push_back( floor( sites.at(s)/ancestry_block_length ) ) ;
    }
    block_sites.push_back( child_chrom.size() ) ;
    
    //// identify blocks before the first breakpoint and push them into the child chromosome
    int break_iterator = 0 ;
    int block_iterator = 0 ;
    while ( break_iterator < child_chrom.size() ) {
        
        //// inherit non-recombinant blocks
        while ( break_iterator < block_sites.at(block_iterator) ) {
            child_chrom.at(break_iterator) = c1.at(break_iterator) ;
            break_iterator ++ ;
        }

        /// now recombine block
        bool hit_before = false ;
        while ( break_iterator == block_sites.at(block_iterator) && block_sites.at(block_iterator) < child_chrom.size() ) {
            /// check if the two parents have identical blocks and just inherit if true
            if ( c1.at(break_iterator) == c2.at(break_iterator) ) {
                child_chrom.at(break_iterator) = c1.at(break_iterator) ;
            }
            else if ( hit_before == false ) {
                ancestry_blocks.push_back( recombine_block( c1.at(break_iterator), c2.at(break_iterator), sites.at(block_iterator) ) ) ;
                child_chrom.at(break_iterator) = ancestry_blocks.back() ;
                hit_before = true ;
            }
            else {
                ancestry_blocks.push_back( recombine_block( child_chrom.at(break_iterator), c2.at(break_iterator), sites.at(block_iterator) ) ) ;
                child_chrom.at(break_iterator) = ancestry_blocks.back() ;
            }
            swap( c1, c2 ) ;
            block_iterator ++ ;
        }
        break_iterator ++ ;
    }
    return child_chrom;
}

ancestry_block* population::recombine_block(ancestry_block *b1, ancestry_block *b2, float &site) {
    ancestry_block *child_block = new ancestry_block ;
    for (auto s = b1->ancestry_switch.begin(); s != b1->ancestry_switch.end(); s++) {
        if (s->first < site) {
            child_block->ancestry_switch.push_back(*s);
        } else {
            break;
        }
    }

    int current_track = child_block->ancestry_switch.back().second;
    list<pair<float, int> >::iterator insertion_site = child_block->ancestry_switch.end();
    for (auto s = b2->ancestry_switch.end(); s != b2->ancestry_switch.begin(); ) {
        s--;
        if (s->first > site) {
            child_block->ancestry_switch.insert(insertion_site, *s);
            insertion_site--;
        } else {
            if (s->second != current_track) {
                pair<float, int> new_switch = make_pair(site, s->second);
                child_block->ancestry_switch.insert(insertion_site, new_switch);
            }
            break;
        }
    }
    for ( int s = 0 ; s < b1->selected_mutations.size() ; s ++ ) {
        if ( b1->selected_mutations.at(s) < site ) {
            child_block->selected_mutations.push_back( b1->selected_mutations.at(s) ) ;
        }
        else {
            break ;
        }
    }
    for ( int s = 0 ; s < b2->selected_mutations.size() ; s ++ ) {
        if ( b2->selected_mutations.at(s) > site ) {
            child_block->selected_mutations.push_back( b2->selected_mutations.at(s) ) ;
        }
    }

    return child_block ;
}

void population::add_offspring() {
    for (int i = 0; i < populations.size(); i++) {
        for (int a = 0; a < populations.at(i).num_anc_male.size(); a++) {
            for (int s = 0; s < populations.at(i).num_anc_male[a]; s++) {
                individual new_male ;
                add_mutations(parental.at(a).males.at(0), new_male, a, true);
                populations.at(i).new_males.push_back(new_male);
            }
        }
        for (int a = 0; a < populations.at(i).num_anc_female.size(); a++) {
            for (int s = 0; s < populations.at(i).num_anc_female[a]; s++) {
                individual new_female ;
                add_mutations(parental.at(a).females.at(0), new_female, a, false);
                populations.at(i).new_females.push_back(new_female);
            }
        }
        swap(populations.at(i).males, populations.at(i).new_males);
        swap(populations.at(i).females, populations.at(i).new_females);
        populations.at(i).new_males.clear();
        populations.at(i).new_females.clear();
    }
    male_offspring.clear();
    female_offspring.clear();
}

void population::print_stats_custom() {
    ofstream output_stream(output.at(curr_output).file, ios::app);

    set<int> to_print;
    int count = 0;
    if (output.at(curr_output).females > 0) {
        while (to_print.size() < output.at(curr_output).females) {
            to_print.insert(gsl_ran_flat(rng, 0, populations.at(output.at(curr_output).subpop).females.size()));
        }
        for (auto it = to_print.begin(); it != to_print.end(); it++) {
            count++;
            populations.at(output.at(curr_output).subpop).females.at(*it).print(chromosome_lengths, output_stream, 0, count, output.at(curr_output).gen, output.at(curr_output).subpop);
        }
    }

    if (output.at(curr_output).males > 0) {
        count = 0;
        to_print.clear();
        while (to_print.size() < output.at(curr_output).males) {
            to_print.insert(gsl_ran_flat(rng, 0, populations.at(output.at(curr_output).subpop).males.size()));
        }
        for (auto it = to_print.begin(); it != to_print.end(); it++) {
            count++;
            populations.at(output.at(curr_output).subpop).males.at(*it).print(chromosome_lengths, output_stream, 1, count, output.at(curr_output).gen, output.at(curr_output).subpop);
        }
    }
}

void population::print_stats(int g) {
    // Average all tract lengths across first subpopulation males
    if (track_lengths) {
        map<char, float> avg;
        map<char, int> overall_total;
        map<char, vector<float> > overall_lengths;
        map<char, int> total;
        map<char, float> tlen;
        map<char, bool> hit;
        int total_fem = 0;
        for (int p = 0; p < populations.size(); p++) {
            for (int i = 0; i < populations.at(p).females.size(); i++) {
                total.clear();
                tlen.clear();
                hit.clear();
                char first_track;
                // Iterate through each ancestry block on a chromosome
                float current_start = populations.at(p).females.at(i).chromosomes.at(1).at(0)->ancestry_switch.front().first ;
                char current_track = populations.at(p).females.at(i).chromosomes.at(1).at(0)->ancestry_switch.front().second ;
                first_track = current_track;
                hit[first_track] = false;
                for ( int a = 0 ; a < populations.at(p).females.at(i).chromosomes.at(1).size() ; a ++ ) {
                    // iterate through each of the ancestry switches
                    for ( auto t = populations.at(p).females.at(i).chromosomes.at(1).at(a)->ancestry_switch.begin() ; t != populations.at(p).females.at(i).chromosomes.at(1).at(a)->ancestry_switch.end() ; t++ ) {
                        if (t->second != current_track) {
                            if (t->first != current_start) {
                               total[current_track]++;
                               overall_total[current_track]++;
                               tlen[current_track] += t->first - current_start;
                               overall_lengths[current_track].push_back(t->first - current_start);
                            }
                            hit[current_track] = true;
                            current_start = t->first;
                            current_track = t->second;
                        }
                    }
                }
                // If there aren't any tract changes, set tlen equal to the chromosome length
                if (!hit[first_track]) {
                    overall_total[first_track]++;
                    total[first_track]++;
                    tlen[first_track] += chromosome_lengths.at(0);
                    overall_lengths[first_track].push_back(chromosome_lengths.at(0));
                }
                // Add onto an existing average of all chromosome 0 tract lengths
                if (total['0'] > 0) {
                    avg['0'] += tlen['0'] / total['0'];
                }
                if (total['1'] > 0) {
                    avg['1'] += tlen['1'] / total['1'];
                }
            }
        }
        float avg0 = avg['0'] / overall_total['0'];
        float avg1 = avg['1'] / overall_total['1'];
        float var0 = 0;
        float var1 = 0;
        for (int i = 0; i < overall_lengths['0'].size(); i++) {
            var0 += (overall_lengths['0'].at(i) - avg0) * (overall_lengths['0'].at(i) - avg0);
        }
        var0 /= overall_total['0'];

        for (int i = 0; i < overall_lengths['1'].size(); i++) {
            var1 += (overall_lengths['1'].at(i) - avg1) * (overall_lengths['1'].at(i) - avg1);
        }
        var1 /= overall_total['1'];
        cout << g << "\t" << avg0 << "\t" << avg1 << "\t" << var0 << "\t" << var1 << "\t" << endl;
    }

    //track the selected site
    if (track_freq) {
        map<float, int> freq;
        int total = 0;
        for (int p = 0; p < populations.size(); p++) {      // iterate through all subpopulations
            for (int m = 0; m < populations.at(p).males.size(); m++) {      // iterate through male individuals first
                for (int c = 0; c < populations.at(p).males.at(m).chromosomes.size(); c++) {    // go through each individual's chromosome
                    for (int a = 0; a < populations.at(p).males.at(m).chromosomes.at(c).size(); a++) {      // access each ancestry block in the chromosome
                        for (int s = 0; s < populations.at(p).males.at(m).chromosomes.at(c).at(a)->selected_mutations.size(); s++) {
                            freq[populations.at(p).males.at(m).chromosomes.at(c).at(a)->selected_mutations.at(s)]++;        // keep a running tally of all selected sites totals          
                        }
                    }
                }
                total++;    // tally the number of individuals sampled
            }
            // do the same for females
            for (int f = 0; f < populations.at(p).females.size(); f++) {
                for (int c = 0; c < populations.at(p).females.at(f).chromosomes.size(); c++) {
                    for (int a = 0; a < populations.at(p).females.at(f).chromosomes.at(c).size(); a++) {
                        for (int s = 0; s < populations.at(p).females.at(f).chromosomes.at(c).at(a)->selected_mutations.size(); s++) {
                                freq[populations.at(p).females.at(f).chromosomes.at(c).at(a)->selected_mutations.at(s)]++;
                        }
                    }
                }
                total++;
            }
        }
        float site = sites_to_track.at(0).at(0);            // particular site we are interested in
        cout << (float) freq[site] / total << endl;         // compute a simple frequency average 
    }   
}

void population::garbage_collect() {
    
    // map to hold the ancestry blocks still in use
    map<ancestry_block*, bool> extant;

    // retain parental blocks
    // no need to record males, they're only a subset of the blocks present in females
    for ( int a = 0 ; a < parental.size() ; a ++ ) {
        for ( int c = 0 ; c < parental.at(a).females.at(0).chromosomes.size() ; c++ ) {
            for ( int b = 0 ; b < parental.at(a).females.at(0).chromosomes.at(c).size() ; b ++ ) {
                extant[parental.at(a).females.at(0).chromosomes.at(c).at(b)] = true ;
            }
        }
    }
    
    // retain used blockss
    for (int p = 0; p < populations.size(); p++) {
        
        for (int f = 0; f < populations.at(p).females.size(); f++) {
            for (int c = 0; c < populations.at(p).females.at(f).chromosomes.size(); c++) {
                for ( int a = 0 ; a < populations.at(p).females.at(f).chromosomes.at(c).size() ; a ++ ) {
                    extant[populations.at(p).females.at(f).chromosomes.at(c).at(a)] = true ;
                }
            }
        }

        for (int m = 0; m < populations.at(p).males.size(); m++) {
            for (int c = 0; c < populations.at(p).males.at(m).chromosomes.size(); c++) {
                for (int a = 0; a < populations.at(p).males.at(m).chromosomes.at(c).size(); a++) {
                    extant[populations.at(p).males.at(m).chromosomes.at(c).at(a)] = true;
                }
            }
        }
    }

    for ( auto a = ancestry_blocks.begin(); a != ancestry_blocks.end(); ) {
        if ( !extant[*a]) {
            delete *a ;
            a = ancestry_blocks.erase( a ) ;
        }
        else {
            a ++ ;
        }
    }
}

