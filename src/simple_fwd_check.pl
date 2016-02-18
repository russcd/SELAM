=begin

This script is written for the purpose of testing SELAM's results across relatively short time scales,
 where the coalescent is expected to fail. 

Results are generally identical provided enough distinct ancestry tracts are sampled (see Manual section 7.2), 
 and this script is substantially less efficeint than SELAM. 

=cut


use strict ; 
use warnings ; 
use Math::Random ;

my $n = 20000 ; 			### this is actually 2n
my $m = $ARGV[0] ; 
my $generations = $ARGV[1] ; 
my $sample_size = 200 ; 

my @pop_pos ; 
my @pop_type ;

foreach ( 1..$n*$m ) { 
	my @pos = ( 0 ) ; 
	my @type = ( 0 ) ; 
	push @pop_pos, \@pos ; 
	push @pop_type, \@type ; 
}

foreach ( $n*$m+1..$n ) { 
	my @pos = ( 0 ) ; 
	my @type = ( 1 ) ; 
	push @pop_pos, \@pos ; 
	push @pop_type, \@type ; 
}

foreach ( 1..$generations ) { 

	### create offspring
	my @new_pos ; 
	my @new_type ; 

	foreach ( 0..$n ) { 

		## id parents / recombination information
		my $parent = Math::Random::random_uniform_integer(1,1,($n/2))*2-1 ; 
		my $chrom = Math::Random::random_uniform_integer(1,0,1) ;
		my $breaks = Math::Random::random_poisson(1,1) ; 		
		my @sites = Math::Random::random_uniform($breaks,0,1) ; 
		push @sites, 1 ; 
		push @sites, 0 ; 
		@sites = sort {$a<=>$b} @sites ; 
		
		## new individual 
		my @pos ; 
		my @type ; 

		### now create chromsomes 
		foreach my $r ( 0..$#sites-1 ) { 
			foreach my $s ( 0..$#{ $pop_pos[$parent+$chrom] } ) { 
				if ( ${ $pop_pos[$parent+$chrom] }[$s] >= $sites[$r] ) {
					if ( ${ $pop_pos[$parent+$chrom] }[$s] < $sites[$r+1] ) {
						if ( $#type > -1 && $type[$#type] == ${ $pop_type[$parent+$chrom] }[$s] ) {
							next ; 
						}  
						push @pos, ${ $pop_pos[$parent+$chrom] }[$s] ;
						push @type, ${ $pop_type[$parent+$chrom] }[$s] ;
					}
					else {
						last ; 
					}
				}
				elsif ( $s == $#{ $pop_pos[$parent+$chrom] } || ${ $pop_pos[$parent+$chrom] }[$s+1] > $sites[$r] ) { 
					if ( $#type > -1 && $type[$#type] == ${ $pop_type[$parent+$chrom] }[$s] ) { 
						next ;
					}
					push @pos, $sites[$r] ; 
					push @type, ${ $pop_type[$parent+$chrom] }[$s] ;
				} 
			}
			if ( $chrom == 0 ) { 
				$chrom = 1 ;
			}
			else {
				$chrom = 0 ; 
			}
		} 

		push @new_pos, \@pos ; 
		push @new_type, \@type ; 
	}

	@pop_pos = @new_pos ; 
	@pop_type = @new_type ; 
}

## simply print the first chromosomes
foreach ( 0..$sample_size-1 ) { 
	foreach my $s ( 0..$#{ $pop_pos[$_] } ) { 
		print $_, "\t", ${ $pop_pos[$_] }[$s], "\t", ${ $pop_type[$_] }[$s], "\n" ;
	}
}
