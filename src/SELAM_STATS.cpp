/// SELAM STATS
/// THIS CODE IS DISTRIBUTED WITH SELAM AND WILL PRODUCE SOME BASIC STATISTICS FOR SELAM OUTPUT

/// VERSION 0.04

#include <algorithm>
#include <string.h> 
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <map>

using namespace std ;

// log\binom{n}{k}
static double lbinom(int n, int k)
{
    if (k == 0 || n == k) return 0;
    return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
static double hypergeo(int n11, int n1_, int n_1, int n)
{
    return exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
}

typedef struct {
    int n11, n1_, n_1, n;
    double p;
} hgacc_t;

// incremental version of hypergenometric distribution
static double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux) {
    if (n1_ || n_1 || n) {
        aux->n11 = n11; aux->n1_ = n1_; aux->n_1 = n_1; aux->n = n;
    } else { // then only n11 changed; the rest fixed
        if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
            if (n11 == aux->n11 + 1) { // incremental
                aux->p *= (double)(aux->n1_ - aux->n11) / n11
                * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
                aux->n11 = n11;
                return aux->p;
            }
            if (n11 == aux->n11 - 1) { // incremental
                aux->p *= (double)aux->n11 / (aux->n1_ - n11)
                * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
                aux->n11 = n11;
                return aux->p;
            }
        }
        aux->n11 = n11;
    }
    aux->p = hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
    return aux->p;
}

double kt_fisher_exact(int n11, int n12, int n21, int n22, double &_left, double &_right, double &two) {
    int i, j, max, min;
    double p, q, left, right;
    hgacc_t aux;
    int n1_, n_1, n;
    n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
    max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
    min = n1_ + n_1 - n;
    if (min < 0) min = 0; // min n11, for left tail
    two = 1 ;
    _left = _right = 1 ;
    if (min == max) return 1.; // no need to do test
    q = hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
    // left tail
    p = hypergeo_acc(min, 0, 0, 0, &aux);
    for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) // loop until underflow
        left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
    --i;
    if (p < 1.00000001 * q) left += p;
    else --i;

    // right tail
    p = hypergeo_acc(max, 0, 0, 0, &aux);
    for (right = 0., j = max - 1; p < 0.99999999 * q; --j) // loop until underflow
        right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
    ++j;
    if (p < 1.00000001 * q) right += p;
    else ++j;
    // two-tail
    two = left + right;
    cerr << two << endl ;
    if (two > 1.) two = 1.;
    // adjust left and right
    if (abs(i - n11) < abs(j - n11)) right = 1. - left + q;
    else left = 1.0 - right + q;
    _left = left; _right = right;
    return q;
}


class ancestry_tract {
public:
	float start;
	float end;
	int type;
} ;

class individual {
public:
	vector<vector<ancestry_tract> > chromosomes ;
    bool male ;
    int index ;
    int subpopulation ;
    int generation ;
} ;

class cmd_line {
public:
	char* input ;
    
    bool G ;
    char* genotype_positions ;
    bool g ;
    double genotype_window ;
    
    bool compute_fet ;
    bool ld ;
    bool lp ;
    bool lm ;
    bool lg ;
    bool lpm ;
    double linkage_window ;
    
	bool allele_frequency ;
	double allele_window ;
    
	void read_cmd_line(int argc, char *argv[]) ;
} ;

class tract_info {
public:
	map<int, double> chrom_lengths;
    vector<individual> sample ;
    
	void read_input(string input);
	void compute_allele_stats(double window);
    void compute_gametic_linkage_stats(double window, int parent, bool fet);
    void compute_composite_linkage_stats(double window ) ;
    void compute_interparental_ld( double window, bool fet ) ;
    void compute_ld_stats( double window, bool fet ) ;
    void output_genotypes( bool g, bool G, double window_size, char* genotype_positions ) ;
    
};

void cmd_line::read_cmd_line(int argc, char *argv[]) {

    input = NULL;
	allele_frequency = false;
	allele_window = 0.0;
    linkage_window = 0.0;
    compute_fet = false ;
    
    ld = false ;
    lg = false ;
    lp = false ;
    lm = false ;
    lpm = false ;
    g = false ;
    G = false ;

	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "--input") == 0 || strcmp(argv[i], "-i") == 0) {
			input = argv[++i];
		}
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--fet") == 0) {
            compute_fet = true;
        }
		if (strcmp(argv[i], "-a") == 0) {
			allele_window = atof(argv[++i]);
			allele_frequency = true;
		}
        if (strcmp(argv[i], "--ld") == 0) {
            linkage_window = atof(argv[++i]);
            ld = true;
        }
        if (strcmp(argv[i], "--lp") == 0) {
            linkage_window = atof(argv[++i]);
            lp = true;
        }
        if (strcmp(argv[i], "--lm") == 0) {
            linkage_window = atof(argv[++i]);
            lm = true;
        }
        if (strcmp(argv[i], "--lg") == 0) {
            linkage_window = atof(argv[++i]);
            lg = true;
        }
        if (strcmp(argv[i], "--lpm") == 0) {
            linkage_window = atof(argv[++i]);
            lpm = true;
        }
        if (strcmp(argv[i], "-g") == 0) {
            genotype_window = atof(argv[++i]);
            g = true;
        }
        if (strcmp(argv[i], "-G") == 0) {
            genotype_positions = argv[++i] ;
            G = true;
        }
	}
}

void tract_info::read_input(string input) {
    
    ifstream in(input);
	string line;
    
	while (getline(in, line)) {
        
		stringstream line_stream;
		line_stream << line;
		
        bool male ;
        int index ;
        int chrom ;
        int parent ;
        int subpopulation ;
        int generation ; 
        
        ancestry_tract new_tract ;
		line_stream >> generation >> subpopulation >> male >> index >> chrom >> parent >> new_tract.type >> new_tract.start >> new_tract.end ;
        
        if ( index > sample.size() ) {
            sample.resize(sample.size()+1) ;
            sample.back().index = index ;
            sample.back().generation = generation ;
            sample.back().subpopulation = subpopulation ;
            sample.back().male = male ;
        }
        
        if ( chrom * 2 + parent == sample.back().chromosomes.size() ) {
            sample.back().chromosomes.resize( sample.back().chromosomes.size() + 1 ) ;
            sample.back().chromosomes.back().push_back(new_tract) ;
        }
        else {
            sample.back().chromosomes.back().push_back(new_tract) ;
        }
        
        if (chrom_lengths[chrom]) {
            if ( chrom_lengths[chrom] < new_tract.end ) {
                chrom_lengths[chrom] = new_tract.end;
            }
        } else {
            chrom_lengths[chrom] = new_tract.end;
		}
	}
}

void tract_info::output_genotypes( bool g, bool G, double window_size, char* genotype_positions ) {
    
    if ( g ) {
        for ( auto c = chrom_lengths.begin() ; c != chrom_lengths.end() ; c ++ ) {
            cout << "##positions:\t" << c->first ;
            for ( double w = 0 ; w <= c->second ; w += window_size ) {
                cout << "\t" << w ;
            }
            cout << endl ;
        }
        
        for ( auto c = chrom_lengths.begin() ; c != chrom_lengths.end() ; c ++ ) {
            for ( int i = 0 ; i < sample.size() ; i ++ ) {
                for ( int p = 0 ; p < 2 ; p ++ ) {
                    /// will skip x chromosomes of males
                    if ( c->first*2+p == sample.at(i).chromosomes.size() ) {
                        continue ;
                    }
                    
                    cout << sample.at(i).generation << "\t" << sample.at(i).subpopulation << "\t" << sample.at(i).male << "\t" << sample.at(i).index << "\t" << c->first << "\t" << p ;
                    int tract_iterator = 0 ;
                    for ( double w = 0 ; w <= c->second ; w += window_size ) {
                        while ( tract_iterator < sample.at(i).chromosomes.at(c->first*2+p).size() - 1 && w > sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).end ) {
                            tract_iterator ++ ;
                        }
                        if ( w > sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).start && w <= sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).end ) {
                            cout << "\t" << sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).type ;
                        }
                    }
                    cout << endl ;
                }
            }
        }
    }
    
    if ( G ) {
        vector<vector<double> > positions ;
        ifstream IN( genotype_positions ) ;
        while( !IN.eof() ) {
            int chrom ;
            double position ;
            IN >> chrom >> position ;
            if ( chrom == positions.size() ) {
                positions.resize(chrom+1) ;
            }
            positions.at(chrom).push_back(position) ;
        }
        positions.back().pop_back() ;
        for ( int c = 0 ; c < positions.size() ; c ++ ) {
            sort(positions.at(c).begin(),positions.at(c).end()) ;
            cout << "##positions:\t" << c ;
            for ( int w = 0 ; w < positions.at(c).size() ; w ++ ) {
                cout << "\t" << positions.at(c).at(w) ;
            }
            cout << endl ;
        }
        for ( int i = 0 ; i < sample.size() ; i ++ ) {
            for ( int p = 0 ; p < 2 ; p ++ ) {
                for ( int c = 0 ; c < positions.size() ; c ++ ) {
                    
                    if ( c * 2 + p == sample.at(i).chromosomes.size() ) {
                        continue ;
                    }
                    
                    cout << sample.at(i).generation << "\t" << sample.at(i).subpopulation << "\t" << sample.at(i).male << "\t" << sample.at(i).index << "\t" << c << "\t" << p ;
                    int tract_iterator = 0 ;

                    for ( int w = 0 ; w < positions.at(c).size() ; w ++ ) {
                        while ( tract_iterator < sample.at(i).chromosomes.at(c*2+p).size() - 1 && positions.at(c).at(w) > sample.at(i).chromosomes.at(c*2+p).at(tract_iterator).end ) {
                            tract_iterator ++ ;
                        }
                        if ( positions.at(c).at(w) > sample.at(i).chromosomes.at(c*2+p).at(tract_iterator).start && positions.at(c).at(w) <= sample.at(i).chromosomes.at(c*2+p).at(tract_iterator).end ) {
                            cout << "\t" << sample.at(i).chromosomes.at(c*2+p).at(tract_iterator).type ;
                        }
                    }
                    cout << endl ;
                }
            }
        }
    }
}

void tract_info::compute_allele_stats(double window_size) {
    
    /// now iterate across all windows
    for ( auto c = chrom_lengths.begin() ; c != chrom_lengths.end() ; c ++ ) {
        
        class window_stats {
        public:
            vector<float> a0 ;
            vector<float> a1 ;
        } ;
        map<float,window_stats > stats ;
        
        for ( int i = 0 ; i < sample.size() ; i ++ ) {
            for ( int p = 0 ; p < 2 ; p ++ ) {
                
                /// will skip x chromosomes of males
                if ( c->first*2+p == sample.at(i).chromosomes.size() ) {
                    continue ;
                }
                
                int tract_iterator = 0 ;
                for ( double w = 0 ; w <= c->second ; w += window_size ) {
                    
                    while ( tract_iterator < sample.at(i).chromosomes.at(c->first*2+p).size() - 1 && w > sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).end ) {
                        tract_iterator ++ ;
                    }
                    
                    if ( w >= sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).start && w <= sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).end ) {
                        
                        if ( sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).type == 0 ) {
                            stats[w].a0.push_back( sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).end - sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).start ) ;
                        }
                        else {
                            stats[w].a1.push_back( sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).end - sample.at(i).chromosomes.at(c->first*2+p).at(tract_iterator).start ) ;
                        }
                    }
                }
            }
        }
        
        for ( auto w = stats.begin() ; w != stats.end() ; w ++ ) {
            
            float prop0 = w->second.a0.size() ;
            prop0 /= ( w->second.a0.size() + w->second.a1.size() ) ;
            
            float mean0 = 0 ;
            float mean1 = 0 ;
            for ( int i = 0 ; i < w->second.a0.size() ; i ++ ) {
                mean0 += w->second.a0.at(i) ;
            }
            for ( int i = 0 ; i < w->second.a1.size() ; i ++ ) {
                mean1 += w->second.a1.at(i) ;
            }
            
            mean0 /= w->second.a0.size() ;
            mean1 /= w->second.a1.size() ;
            
            float var0 = 0 ;
            float var1 = 0 ;
            for ( int i = 0 ; i < w->second.a0.size() ; i ++ ) {
                var0 += pow( mean0 - w->second.a0.at(i), 2 ) ;
            }
            for ( int i = 0 ; i < w->second.a1.size() ; i ++ ) {
                var1 += pow( mean1 - w->second.a1.at(i), 2 ) ;
            }
            
            var0 /= w->second.a0.size() ;
            var1 /= w->second.a1.size() ;
            
            cout << c->first << "\t" << w->first << "\t" ;
            cout << prop0 << "\t" << mean0 << "\t" ;
            cout << mean1 << "\t" << var0 << "\t" << var1 << endl ;
        }
    }
}

void tract_info::compute_ld_stats ( double window_size, bool fet ) {
    
    /// iterate across first chromosome
    for ( auto c1 = chrom_lengths.begin() ; c1 != chrom_lengths.end() ; c1 ++ ) {
        vector<vector<int> > iterator1 ( 2, vector<int> (sample.size(), 0 ) ) ;
        for ( double w1 = 0 ; w1 < c1->second ; w1 += window_size ) {
            
            /// iterate across first chromosome
            for ( auto c2 = c1 ; c2 != chrom_lengths.end() ; c2 ++ ) {
                vector<vector<int> > iterator2 ( 2, vector<int> (sample.size(), 0 ) ) ;
                
                double w2_start = 0 ;
                if ( c1 == c2 ) {
                    w2_start = w1 + window_size ;
                }
                
                for ( double w2 = w2_start ; w2 <= c2->second ; w2 += window_size ) {
                    
                    /// data structure for linkage stats
                    vector< vector< float > > ldm ( 2, vector<float> ( 2, 0 ) ) ;
                    float total = 0 ;
                    
                    /// iterate through all individuals
                    for ( int i = 0 ; i < sample.size() ; i ++ ) {
                        for ( int parent = 0 ; parent < 2 ; parent ++ ) {
                            
                            if ( c1->first*2+parent == sample.at(i).chromosomes.size() || c2->first*2+parent == sample.at(i).chromosomes.size() ) {
                                continue ;
                            }
                        
                            /// update iterators
                            while ( w1 > sample.at(i).chromosomes.at(c1->first*2+parent).at(iterator1.at(parent).at(i)).end ) {
                                iterator1.at(parent).at(i) ++ ;
                            }
                            while ( w2 > sample.at(i).chromosomes.at(c2->first*2+parent).at(iterator2.at(parent).at(i)).end ) {
                                iterator2.at(parent).at(i) ++ ;
                            }
                        
                            /// now update stats
                            total ++ ;
                            ldm[sample.at(i).chromosomes.at(c1->first*2+parent).at(iterator1.at(parent).at(i)).type][sample.at(i).chromosomes.at(c2->first*2+parent).at(iterator2.at(parent).at(i)).type] ++ ;
                        }
                    }
                    
                    double p, right, left ;
                    if ( fet ) {
                        double q = kt_fisher_exact( ldm.at(0).at(0), ldm.at(0).at(1), ldm.at(1).at(0), ldm.at(1).at(1), left, right, p ) ;
                    }
                    
                    for ( int i = 0 ; i < ldm.size() ; i ++ ) {
                        for ( int l = 0 ; l < ldm.at(i).size() ; l ++ ) {
                            ldm.at(i).at(l)/=total ;
                        }
                    }
                    
                    /// compute linkage stats
                    double d = ldm.at(0).at(0) - ( ldm.at(0).at(0)+ldm.at(0).at(1) )*(ldm.at(0).at(0)+ldm.at(1).at(0)) ;
                    
                    double r = d / sqrt( ( ( ldm.at(0).at(0)+ldm.at(0).at(1)) ) * ( 1 - ( ldm.at(0).at(0)+ldm.at(0).at(1) ) )
                                        * ( (ldm.at(0).at(0)+ldm.at(1).at(0)) )* ( 1- (ldm.at(0).at(0)+ldm.at(1).at(0)) ) ) ;
                    
                    cout << c1->first << "\t" << w1 << "\t" ;
                    cout << c2->first << "\t" << w2 << "\t" ;
                    cout << ldm.at(0).at(0) + ldm.at(0).at(1) << "\t" ;
                    cout << ldm.at(0).at(0) + ldm.at(1).at(0) << "\t" ;
                    cout << total << "\t" << d << "\t" << r ;
                    
                    if ( fet ) {
                        cout << "\t" << p  ;
                    }
                    cout << endl ;
                }
            }
        }
    }
}


void tract_info::compute_gametic_linkage_stats ( double window_size, int parent, bool fet ) {
    
    /// iterate across first chromosome
    for ( auto c1 = chrom_lengths.begin() ; c1 != chrom_lengths.end() ; c1 ++ ) {
        vector<int> iterator1 ( sample.size(), 0 ) ;
        for ( double w1 = 0 ; w1 < c1->second ; w1 += window_size ) {
            
            /// iterate across first chromosome
            for ( auto c2 = c1 ; c2 != chrom_lengths.end() ; c2 ++ ) {
                vector<int> iterator2 ( sample.size(), 0 ) ;
                
                double w2_start = 0 ;
                if ( c1 == c2 ) {
                    w2_start = w1 + window_size ;
                }
                
                for ( double w2 = w2_start ; w2 <= c2->second ; w2 += window_size ) {
                    
                    /// data structure for linkage stats
                    vector< vector< float > > ldm ( 2, vector<float> ( 2, 0 ) ) ;
                    float total = 0 ;
                    
                    /// iterate through all individuals
                    for ( int i = 0 ; i < sample.size() ; i ++ ) {
                        
                        if ( c1->first*2+parent == sample.at(i).chromosomes.size() || c2->first*2+parent == sample.at(i).chromosomes.size() ) {
                            continue ;
                        }
                        
                        /// update iterators
                        while ( w1 > sample.at(i).chromosomes.at(c1->first*2+parent).at(iterator1.at(i)).end ) {
                            iterator1.at(i) ++ ;
                        }
                        while ( w2 > sample.at(i).chromosomes.at(c2->first*2+parent).at(iterator2.at(i)).end ) {
                            iterator2.at(i) ++ ;
                        }
                        
                        /// now update stats
                        total ++ ;
                        ldm[sample.at(i).chromosomes.at(c1->first*2+parent).at(iterator1.at(i)).type][sample.at(i).chromosomes.at(c2->first*2+parent).at(iterator2.at(i)).type] ++ ;
                    }
                    
                    double p, right, left ;
                    if ( fet ) {
                        double q = kt_fisher_exact( ldm.at(0).at(0), ldm.at(0).at(1), ldm.at(1).at(0), ldm.at(1).at(1), left, right, p ) ;
                    }
                    
                    for ( int i = 0 ; i < ldm.size() ; i ++ ) {
                        for ( int l = 0 ; l < ldm.at(i).size() ; l ++ ) {
                            ldm.at(i).at(l)/=total ;
                        }
                    }
                    
                    /// compute linkage stats
                    double d = ldm.at(0).at(0) - ( ldm.at(0).at(0)+ldm.at(0).at(1) )*(ldm.at(0).at(0)+ldm.at(1).at(0)) ;
                    
                    double r = d / sqrt( ( ( ldm.at(0).at(0)+ldm.at(0).at(1)) ) * ( 1 - ( ldm.at(0).at(0)+ldm.at(0).at(1) ) )
                                        * ( (ldm.at(0).at(0)+ldm.at(1).at(0)) )* ( 1- (ldm.at(0).at(0)+ldm.at(1).at(0)) ) ) ;
                    
                    cout << c1->first << "\t" << w1 << "\t" ;
                    cout << c2->first << "\t" << w2 << "\t" ;
                    cout << ldm.at(0).at(0) + ldm.at(0).at(1) << "\t" ;
                    cout << ldm.at(0).at(0) + ldm.at(1).at(0) << "\t" ;
                    cout << total << "\t" << d << "\t" << r ;
                    
                    if ( fet ) {
                        cout << "\t" << p  ;
                    }
                    cout << endl ;
                }
            }
        }
    }
}

void tract_info::compute_composite_linkage_stats( double window_size ) {
    
    /// iterate across first chromosome
    for ( auto c1 = chrom_lengths.begin() ; c1 != chrom_lengths.end() ; c1 ++ ) {
        vector<int> iterator11 ( sample.size(), 0 ) ;
        vector<int> iterator12 ( sample.size(), 0 ) ;
        for ( double w1 = 0 ; w1 < c1->second ; w1 += window_size ) {
            
            /// iterate across first chromosome
            for ( auto c2 = c1 ; c2 != chrom_lengths.end() ; c2 ++ ) {
                vector<int> iterator21 ( sample.size(), 0 ) ;
                vector<int> iterator22 ( sample.size(), 0 ) ;
                
                double w2_start = 0 ;
                if ( c1 == c2 ) {
                    w2_start = w1 + window_size ;
                }
                
                for ( double w2 = w2_start ; w2 <= c2->second ; w2 += window_size ) {
                    
                    /// data structure for linkage stats
                    vector< vector< float > > ld ( 3, vector<float> ( 3, 0 ) ) ;
                    float total = 0 ;
                    
                    /// iterate through all individuals
                    for ( int i = 0 ; i < sample.size() ; i ++ ) {
                        
                        if ( c1->first*2+1 == sample.at(i).chromosomes.size() || c2->first*2+1 == sample.at(i).chromosomes.size() ) {
                            continue ;
                        }
                        
                        /// update iterators
                        while ( w1 > sample.at(i).chromosomes.at(c1->first*2).at(iterator11.at(i)).end ) {
                            iterator11.at(i) ++ ;
                        }
                        while ( w2 > sample.at(i).chromosomes.at(c2->first*2).at(iterator21.at(i)).end ) {
                            iterator21.at(i) ++ ;
                        }
                        while ( w1 > sample.at(i).chromosomes.at(c1->first*2+1).at(iterator12.at(i)).end ) {
                            iterator12.at(i) ++ ;
                        }
                        while ( w2 > sample.at(i).chromosomes.at(c2->first*2+1).at(iterator22.at(i)).end ) {
                            iterator22.at(i) ++ ;
                        }
                        
                        /// now update stats
                        total ++ ;
                        ld[sample.at(i).chromosomes.at(c1->first*2).at(iterator11.at(i)).type
                           +sample.at(i).chromosomes.at(c1->first*2+1).at(iterator12.at(i)).type]
                           [sample.at(i).chromosomes.at(c2->first*2).at(iterator21.at(i)).type
                            +sample.at(i).chromosomes.at(c2->first*2+1).at(iterator22.at(i)).type] ++ ;
                    }
                    
                    double Na1b1 = 2 * ld[0][0] + ld[0][1] + ld[1][0] + 0.5 * ld[1][1] ;
                    double pa1 = ( ld[0][0]+ld[0][1]+ld[0][2] )/total ;
                    double pb1 = ( ld[0][0]+ld[1][0]+ld[2][0] )/total ;
                    
                    double da1b1 = Na1b1/total - 2*pa1*pb1 ;
                    
                    cout << c1->first << "\t" << w1 << "\t" ;
                    cout << c2->first << "\t" << w2 << "\t" ;
                    cout << pa1 << "\t" << pb1 << "\t" ;
                    cout << total << "\t" << da1b1 << endl ;
                }
            }
        }
    }
}

void tract_info::compute_interparental_ld ( double window_size, bool fet ) {
    
    /// iterate across first chromosome
    for ( auto c1 = chrom_lengths.begin() ; c1 != chrom_lengths.end() ; c1 ++ ) {
        vector<int> iterator1 ( sample.size(), 0 ) ;
        for ( double w1 = 0 ; w1 < c1->second ; w1 += window_size ) {
            /// iterate across second chromosome
            for ( auto c2 = chrom_lengths.begin() ; c2 != chrom_lengths.end() ; c2 ++ ) {
                vector<int> iterator2 ( sample.size(), 0 ) ;
                for ( double w2 = 0 ; w2 <= c2->second ; w2 += window_size ) {
                    
                    vector< vector< float > > ld ( 2, vector<float> ( 2, 0 ) ) ;
                    float total = 0 ;
                    
                    /// iterate through all individuals
                    for ( int i = 0 ; i < sample.size() ; i ++ ) {
                        
                        if ( c1->first*2 == sample.at(i).chromosomes.size() || c2->first*2+1 == sample.at(i).chromosomes.size() ) {
                            continue ;
                        }
                        
                        /// update iterators
                        while ( w1 > sample.at(i).chromosomes.at(c1->first*2).at(iterator1.at(i)).end ) {
                            iterator1.at(i) ++ ;
                        }
                        while ( w2 > sample.at(i).chromosomes.at(c2->first*2+1).at(iterator2.at(i)).end ) {
                            iterator2.at(i) ++ ;
                        }
                        
                        /// now update stats
                        total ++ ;
                        ld[sample.at(i).chromosomes.at(c1->first*2).at(iterator1.at(i)).type][sample.at(i).chromosomes.at(c2->first*2+1).at(iterator2.at(i)).type] ++ ;
                    }
                    
                    double p, right, left ;
                    if ( fet ) {
                        double q = kt_fisher_exact( ld.at(0).at(0), ld.at(0).at(1), ld.at(1).at(0), ld.at(1).at(1), left, right, p ) ;
                    }
                    
                    for ( int i = 0 ; i < ld.size() ; i ++ ) {
                        for ( int l = 0 ; l < ld.at(i).size() ; l ++ ) {
                            ld.at(i).at(l)/=total ;
                        }
                    }
                    
                    /// compute linkage stats
                    double d = ld.at(0).at(0) - ( ld.at(0).at(0)+ld.at(0).at(1) )*(ld.at(0).at(0)+ld.at(1).at(0)) ;
                    
                    double r = d / sqrt( ( ( ld.at(0).at(0)+ld.at(0).at(1)) ) * ( 1 - ( ld.at(0).at(0)+ld.at(0).at(1) ) )
                                        * ( (ld.at(0).at(0)+ld.at(1).at(0)) )* ( 1- (ld.at(0).at(0)+ld.at(1).at(0)) ) ) ;
                    
                    cout << c1->first << "\t" << w1 << "\t" ;
                    cout << c2->first << "\t" << w2 << "\t" ;
                    cout << ld.at(0).at(0) + ld.at(0).at(1) << "\t" ;
                    cout << ld.at(0).at(0) + ld.at(1).at(0) << "\t" ;
                    cout << total << "\t" << d << "\t" << r ;
                    
                    if ( fet ) {
                        cout << "\t" << p  ;
                    }
                    cout << endl ;
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    
	cmd_line options ;
	options.read_cmd_line(argc, argv) ;
	tract_info info ;
    
    info.read_input( options.input ) ;
    
    if ( options.g || options.G ) {
        info.output_genotypes( options.g, options.G, options.genotype_window, options.genotype_positions ) ;
    }
    
	if ( options.allele_frequency ) {
        info.compute_allele_stats( options.allele_window ) ;
	}
    
    if ( options.lp || options.lm ) {
        info.compute_gametic_linkage_stats(options.linkage_window, options.lp, options.compute_fet ) ;
    }
    
    if ( options.lg ) {
        info.compute_composite_linkage_stats(options.linkage_window ) ;
    }
    
    if ( options.lpm ) {
        info.compute_interparental_ld( options.linkage_window, options.compute_fet ) ;
    }
    
    if ( options.ld ) {
        info.compute_ld_stats( options.linkage_window, options.compute_fet ) ;
    }

	return 0 ;
}


