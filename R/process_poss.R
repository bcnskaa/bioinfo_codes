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
#source("process_poss.R")

# Process a list of dataframes, and convert them into poss objects based on given start and end position labels
# Assuming all the dataframes of list_df share the same spos_label_idx and epos_label_idx
process_list_to_poss <- function(list_df, spos_label_idx, epos_label_idx) {
	df_names <- names(list_df);
	
	print(paste("Converting ", length(df_names), " dataframe(s) into poss objects...",sep=""));
	
	list_df_poss <- list();
	for(df_name in df_names)
	{
		cur_df <- list_df[[df_name]];
		list_df_poss[[df_name]] <- process_poss_c(cur_df[,spos_label_idx], cur_df[,epos_label_idx]);
	}
	
	return(list_df_poss);
}


process_poss_c <- function(sposs, eposs) {
	
	dyn.load("process_poss.so");
	poss <- .Call("process_poss_r", as.integer(sposs), as.integer(eposs));
	
	dyn.unload("process_poss.so");
	return(poss);
}

#
# Construct a position-orientated list of values of length specified by len_m
#
# 
#res<-process_poss_vals_c(1000, c(100, 200, 400), c(130, 240, 450), c(1.1,2.1,0.2))
process_poss_vals_c <- function(sposs, eposs, vals, len_m)
{
	dyn.load("process_poss.so");
	print(paste("Max length: ", len_m));
	poss <- .Call("process_poss_vals_r", as.integer(sposs), as.integer(eposs), as.double(vals), as.integer(len_m));

	dyn.unload("process_poss.so");	
	return(poss);
}
	
#
#
#
convert_to_poss <- function(ranges) {
	return(process_poss_c(ranges$SPOS, ranges$EPOS));
}

#
#
#
convert_to_ranges <- function(ids, sposs, eposs)
{
	range_df <- data.frame(ID=ids, SPOS=sposs, EPOS=eposs);
	
	return(range_df);
}


#cut_poss_to_ranges <- function(ids, sposs, eposs, poss)
cut_poss_to_ranges <- function(poss, range_df)
{
	range_n <- nrow(range_df);
	print(paste("Number of ranges: ", range_n, sep=""));
	range_poss_vals <- rep(0.0, range_n);
	
	for(i in 1 : range_n)
	{
		range_poss_vals[i] <- sum(poss[poss >= sposs[i] & poss < eposs[i]]);
	}
	
	return(range_poss_vals);
}


# process_poss_vals_to_ranges_r(SEXP s, SEXP e, SEXP v, SEXP rs, SEXP re, SEXP len)
process_poss_vals_to_ranges <- function(sposs, eposs, vals, range_df, max_len)
{
	dyn.load("process_poss.so");
	print(paste("Max length: ", max_len));
	#bins <- .Call("process_poss_vals_r", as.integer(sposs), as.integer(eposs), as.double(vals), as.integer(range_df$START), as.integer(range_df$END), as.integer(max_len));
	bins <- .Call("process_poss_vals_to_ranges_r", as.integer(sposs), as.integer(eposs), as.double(vals), as.integer(range_df$START), as.integer(range_df$END), max_len);
	
	if(length(bins) != nrow(range_df))
	{
		print(paste("Something wrong with the data."));
	}
	
	dyn.unload("process_poss.so");
	return(bins);
}
