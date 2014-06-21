#!/usr/bin/perl -w

#########################
#
# Copyright (c) 2013, SK Woolf <bcnskaa@gmail.com>.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#########################


#
# Usage: 
# perl process_vcf_to_freq.pl VCF-INFILE VCF-INFO-INFILE > OUT-FILE
#


# 1KGP VCF file
$fn = shift;
# 1KGP VCF info file, 
$fn_info = shift;


%POPS = ();
%INFO = ();
%GENDERS = ();
open(IN, "<", $fn_info) or die $!;
if(defined IN) {
	my $line_count = 0;
	while(<IN>) {
		chomp $_;
		
		if($line_count == 0) {
			# Header line
			# Do nothing
		} else {
			my @items = split /\t/, $_;
			
			$INFO{ $items[1] } = $items[0];
			$GENDERS{ $items[1] } = $items[4];
			
			#print qq{$items[4]\n};
			
			# Add a new population
			if(exists $POPS{$items[0]}) {
				$POPS{ $items[0] } = 0;
			}
		}
		
		$line_count++;
	}
	close(IN);
}

open(IN, "<", $fn) or die $!;
if(defined IN) {

	%IND_INFO = ();
	@IND_IDS = ();
	# Inidividual grouped under a same population
	%POP_IND = ();
	while(<IN>) {
		chomp $_;
		
		my @items = split /\t/, $_;
		
		my $chr = shift @items;
		my $pos = shift @items;
		my $ref = shift @items;
		
		# Process the header line
		if($_ =~ /^#/) {
			my $index = 0;
			foreach my $ind_id (@items) {
				push @IND_IDS, $ind_id;
				$IND_INFO{ $ind_id }{ "INDEX" } = $index;
				$IND_INFO{ $ind_id }{ "POPULATION" } = $INFO{ $ind_id };
				$IND_INFO{ $ind_id }{ "GENDER" } = $GENDERS{ $ind_id };
				push @{$POP_IND{$INFO{ $ind_id }}}, $ind_id;
				$index++;
			}	
		} else {	
			# Total genotype and allele counts
			my %GENO = ();
			my %ALLELES = ();
			
			# Population based 
			my %POP_GENO = ();
			my %POP_ALLELES = ();
			
			# Gender based
			my %GENDER_GENO = ();
			my %GENDER_ALLELES = ();
						
			
			my $index = 0;
			foreach my $cur_geno (@items) {
				chomp $cur_geno;
				my ($allele1, $allele2) = split /\//, $cur_geno;
				my $ind_id = $IND_IDS[$index];
				my $pop_id = $IND_INFO{ $ind_id }{ "POPULATION" };
				my $gender = $IND_INFO{ $ind_id }{ "GENDER" };
				# Population based
				if(exists $POP_ALLELES{ $pop_id }{ $allele1 }) {
					$POP_ALLELES{ $pop_id }{ $allele1 }++;
				} else {
					$POP_ALLELES{ $pop_id }{ $allele1 } = 1;
				}
				
				if(exists $POP_ALLELES{ $pop_id }{ $allele2 }) {
					$POP_ALLELES{ $pop_id }{ $allele2 }++;
				} else {
					$POP_ALLELES{ $pop_id }{ $allele2 } = 1;
				}		
				if(exists $POP_GENO{ $pop_id }{ $cur_geno }) {
					$POP_GENO{ $pop_id }{ $cur_geno }++;
				} else {
					$POP_GENO{ $pop_id }{ $cur_geno } = 1;
				}
				
							
				# Gender based
				if(exists $GENDER_ALLELES{ $gender }{ $allele1 }) {
					$GENDER_ALLELES{ $gender }{ $allele1 }++;
				} else {
					$GENDER_ALLELES{ $gender }{ $allele1 } = 1;
				}
				
				if(exists $GENDER_ALLELES{ $gender }{ $allele2 }) {
					$GENDER_ALLELES{ $gender }{ $allele2 }++;
				} else {
					$GENDER_ALLELES{ $gender }{ $allele2 } = 1;
				}		
				if(exists $GENDER_GENO{ $gender }{ $cur_geno }) {
					$GENDER_GENO{ $gender }{ $cur_geno }++;
				} else {
					$GENDER_GENO{ $gender }{ $cur_geno } = 1;
				}			
				

				# Total record
				if(exists $ALLELES{$allele1}) {
					$ALLELES{$allele1}++;
				} else {
					$ALLELES{$allele1} = 1;
				}
				if(exists $ALLELES{$allele2}) {
					$ALLELES{$allele2}++;
				} else {
					$ALLELES{$allele2} = 1;
				}

				if(exists $GENO{$cur_geno}) {
					$GENO{$cur_geno}++;
				} else {
					$GENO{$cur_geno} = 1;
				}
				
				# Next individual genotype
				$index++;
			}

			#$size = keys(%GENO);

			print qq{$chr\t$pos\t$ref\t};

			for my $allele (sort {$ALLELES{$b} <=> $ALLELES{$a}} (keys (%ALLELES))) {
				print $allele . ":" . $ALLELES{$allele} . ";";
			}
		
			print qq{\t};
			
			
			
				
			##### Total genotype and alleles
			my $total_allele_count = 0;
			my $maf_allele_count = 0;
			my $maf_freq = 0.0;
			for my $allele (sort {$ALLELES{$b} <=> $ALLELES{$a}} (keys (%ALLELES))) {
				#print $allele . ":" . $ALLELES{$allele} . ";";

				my $count = $ALLELES{$allele};

				$total_allele_count += $count;
				if($count > $maf_allele_count) {
 					$maf_allele_count = $count;
				}

                	}
                	$maf_freq = $maf_allele_count / $total_allele_count;
			print qq{$maf_freq};
			
			
			print qq{\t};			
			for my $geno (sort {$GENO{$b} <=> $GENO{$a}} (keys (%GENO))) {
				print $geno . ":" . $GENO{$geno} . ";";
			}
			
			
			##### Gender Alleles
			print qq{\t};		
			for my $my_gender (keys (%GENDER_ALLELES)) {
				print qq{$my_gender=};
				my $my_item_count = 0;
				for my $my_allele (keys (%{$GENDER_ALLELES{$my_gender}})) {
					my $my_count = $GENDER_ALLELES{$my_gender}{$my_allele};
					if($my_item_count > 0) {
						print qq{,};
					}
					
					print qq{$my_allele:$my_count};
					$my_item_count++;
				}
				
				print qq{;};
			}
						
			##### Gender Genotypes
			print qq{\t};		
			for my $my_gender (keys (%GENDER_GENO)) {
				print qq{$my_gender=};
				my $my_item_count = 0;
				for my $my_geno (keys (%{$GENDER_GENO{$my_gender}})) {
					my $my_count = $GENDER_GENO{$my_gender}{$my_geno};
					if($my_item_count > 0) {
						print qq{,};
					}
					
					print qq{$my_geno:$my_count};
					$my_item_count++;
				}
				
				print qq{;};
			}	
					
			##### Population Alleles		
			print qq{\t};		
			for my $my_pop (keys (%POP_ALLELES)) {
				print qq{$my_pop=};
				my $my_item_count = 0;
				for my $my_allele (keys (%{$POP_ALLELES{$my_pop}})) {
					my $my_count = $POP_ALLELES{$my_pop}{$my_allele};
					if($my_item_count > 0) {
						print qq{,};
					}
					
					print qq{$my_allele:$my_count};
					$my_item_count++;
				}
				
				print qq{;};
			}
			
			##### Population Genotypes	
			print qq{\t};			
			for my $my_pop (keys (%POP_GENO)) {
				print qq{$my_pop=};
				my $my_item_count = 0;
				for my $my_allele (keys (%{$POP_GENO{$my_pop}})) {
					my $my_count = $POP_GENO{$my_pop}{$my_allele};
					if($my_item_count > 0) {
						print qq{,};
					}
					
					print qq{$my_allele:$my_count};
					$my_item_count++;
				}
				
				print qq{;};
			}	
			
			
					
			print qq{\n};

		}
	}
	close IN;
}

# Sort based on allele frequency
sub sort_alele_freq {
	$ALLELES{$b} <=> $ALLELES{$a};
}

# Sort based on genotype frequency
sub sort_genotype_freq {
	$GENO{$b} <=> $GENO{$a};	
}
