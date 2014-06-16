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
process_fasta_c <- function(fasta_fn, header_id, base_to_be_counted) {
	
	dyn.load("process_fasta.so");
	#poss <- .Call("process_fasta_r", "/home/bcnskaa/share/db/gatk_resources/ucsc.hg19.fasta", "chr19", "C");
	poss <- .Call("process_fasta_r", fasta_fn, header_id, base_to_be_counted);
	
	return(poss);
}


process_fasta_bin_c <- function(fasta_fn, header_id, bin_size, base_to_be_counted) {
	
	dyn.load("process_fasta.so");
	#poss <- .Call("process_fasta_r", "/home/bcnskaa/share/db/gatk_resources/ucsc.hg19.fasta", "chr19", "C");
	poss <- .Call("process_fasta_r", fasta_fn, header_id, base_to_be_counted);
	
	return(poss);
}


get_subset <- function(fasta_fn, header_id, spos, epos)
{
	
	dyn.load("process_fasta.so");
	#poss <- .Call("process_fasta_r", "/home/bcnskaa/share/db/gatk_resources/ucsc.hg19.fasta", "chr19", "C");
	seq <- .Call("get_subseq_r", fasta_fn, header_id, as.integer(spos), as.integer(epos));

	return(seq);
}


test_code <- function() {
	test_loop <- 10000;
	pattern_len <- 11;
	nucleotide_bases <- get_nucleotide_bases();

	patterns <- generate_permutated_seq_patterns(pattern_len);
	
	fail_n <- 0;
	for(i in 1 : test_loop)
	{
		pattern <- paste(sample(nucleotide_bases, pattern_len, replace=T), sep="", collapse="");
		
		idx <- calculate_pattern_permutation_index(pattern);
		
		print(paste("Pattern=", pattern, " patterns[", idx, "]=", paste(patterns[idx,], sep="", collapse=""), sep=""));
		
		if(pattern != paste(patterns[idx,], sep="", collapse=""))
		{
			fail_n = fail_n + 1;
		}
	}
	
	print(paste("Test: ", fail_n, " out of ", test_loop, " failed.", sep=""));
}


# patterns <- generate_permutated_seq_patterns(3)
generate_permutated_seq_patterns <- function(pattern_len, N_included=T)
{
	#library(combinat);
	#patterns <- permn();
	
	library(gtools);
	
	nucleotide_bases <- get_nucleotide_bases(N_included);
	
	patterns <- permutations(length(nucleotide_bases), pattern_len, nucleotide_bases, repeats.allowed=T);
	
	return(patterns);
}


# Calculate the position of given pattern in a list of permutated sequence patterns
calculate_pattern_permutation_index <- function(pattern, N_included=T)
{
	nucleotide_bases <- get_nucleotide_bases(N_included);
	pattern_len <- nchar(pattern);
	base_rank <- length(nucleotide_bases);
	
	return(sum((base_rank ** seq(pattern_len - 1, 0)) * sapply(seq(1, pattern_len), FUN=function(i){get_base_value(substr(pattern, i, i), nucleotide_bases)})) + 1)
	# (1*(4^2)) + (2*(4^1)) + (1*(4^0)
	
}


get_nucleotide_bases <- function(N_included=T)
{
	if(N_included){
		return(c("A", "C", "G", "N", "T"));
	} else {
		return(c("A", "C", "G", "T"));
	}
}


get_base_value <- function(base, nucleotide_bases)
{
	return(which(nucleotide_bases == base) - 1);
}


# Calculate the heximal value
calculate_pattern_hex_code <- function(pattern, N_included=T)
{
	nucleotide_bases <- get_nucleotide_bases(N_included);
	
	pos <- 0
	hex_val <- 0
	for(i in nchar(pattern):1) 
	{
		base <- substr(pattern, i, i);
		hex_val = hex_val + (get_base_value(base, nucleotide_bases) * (10 ^ pos));
		pos = 1 + pos;
	}
	hex_val = hex_val + 1;
	
	return(hex_val)
}
