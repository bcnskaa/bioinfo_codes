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
# Constants
META_LIST_META_LABEL = "META";
META_LIST_DATA_LABEL = "DATA";
META_CLASS_LABEL_PATTERN_DELIMITOR = ";";


# Input: a list of dataframes with identical number and name of columns
# Output: a dataframe containing all elements of dataframes and a new column specifies name of the dataframe from which element is draw from
# Argument: df_list=a list of dataframes with identical number and name of columns
#	    new_class_label=the name of new column listing the source dataframes  
melt_df <- function(df_list, new_class_label="class") {
	if(length(unique(sapply(df_list, function(v) dim(v)[2]))) != 1)
	{
		print(paste("Input list contains dataframes of different columns' number/name, unable to melt the list."))
		return
	}
	
	class_labels <- names(df_list);
		
	melted_df <- data.frame();
	

	for(label in class_labels) {
		print(paste("Processing ", label, "...", sep=""));
		
		cur_df <- df_list[[label]];

		cur_df <- cbind(cur_df, X______X_____X=rep(label, nrow(cur_df)));
		if(dim(melted_df)[1] == 0)
		{
			melted_df <- cur_df;
		} else {
			melted_df <- rbind(melted_df, cur_df);
		}
	}

	colnames(melted_df)[colnames(melted_df) == "X______X_____X"] <- new_class_label;
	
	return(melted_df);
}



#
#
# 
process_df_with_label_vals <- function(dataframe, label_idx, label_vals, label_n_cutoff = 0)
{
	select_labels <- as.character(label_vals);
	
	label_n <- length(select_labels);
	
	# Create a new data.frame
	df <- list();
	
	for(i in 1 : label_n)
	{
		current_label <- select_labels[i];
		df[[current_label]] <- dataframe[which(dataframe[,label_idx] == current_label),];
	}
	
	return(df);
}

# Construct a dataframe based on label values in regular expression form
# Select items with labels started with "AluY", we should put "^AluY" in label_val_regex field
process_df_with_label_regex <- function(dataframe, label_idx, label_val_regex, label_n_cutoff = 0)
{
	#unique_labels <- as.character(unique(df[, label_idx]));
	counts <- table(dataframe[, label_idx]);
	counts <- counts[grep(label_val_regex, names(counts))];
	select_labels <- names(counts[counts >= label_n_cutoff]);
	
	return(process_df_with_label_vals(dataframe, label_idx, select_labels, labal_n_cutoff));
}


#
# construct a dataframe which contains a list of collections with unique label value
#
process_df <- function(dataframe, label_idx, label_n_cutoff = 0)
{
	counts <- table(dataframe[, label_idx])
	select_labels <- names(counts[counts >= label_n_cutoff]);
	
	return(process_df_with_label_vals(dataframe, label_idx, select_labels, labal_n_cutoff));
}



test_code_convert_meta_df_to_bin <- function() 
{
	sample_n <- 60000000;
	bin_size <- 2000000;
	label <- seq(1:sample_n);
	breaks <- c(seq(0,sample_n, by=bin_size));
	data <- sample(seq(0, 1, by=0.0001), sample_n, replace=T);
	poss_idx <- cut(label, breaks=breaks, labels=F);
	
	res <- sapply(seq(1:(length(breaks) - 1)), function(i) mean(data[poss_idx == i]));
	
}



#
factorize_meta_df <- function(meta_df_list, replace_original=T, enable_class_label_pattern=F, fixed_class_label_pattern=F)
{
	source("list_count.R");
	
	metum <- meta_df_list[[META_LIST_META_LABEL]];
	datum <- meta_df_list[[META_LIST_DATA_LABEL]];

	####################
	# bug here
	###################
	data_n <- length(datum);
	
	# Go through all dataset and see if they have subclass labels
	for(i in 1 : data_n)
	{
		meta <- metum[i, ];
		name <- as.character(meta$name);
		chrom_label <- as.character(meta$chrom_label);
		spos_label <- as.character(meta$spos_label);
		epos_label <- as.character(meta$epos_label);
		datatype <- as.character(meta$datatype);
		class_label <- as.character(meta$class_label);
		class_label_pattern <- as.character(meta$class_label_pattern);
		value_label <- as.character(meta$value_label);
		
		
		data <- datum[[name]];
		
		# Have subclasses
		if(nchar(class_label) != 0)
		{
			if(nchar(class_label_pattern) != 0)
			{
				class_label_patterns <- unlist(strsplit(class_label_pattern, META_CLASS_LABEL_PATTERN_DELIMITOR)); 
				
				for(pattern in class_label_patterns) 
				{
					print(paste("Pattern: ", pattern, sep=""));
					selected_idx <- grep(pattern, data[, class_label]);
					
					if(length(selected_idx) > 0)
					{
						df <- data[selected_idx,];
						df_meta <- meta;
						df_meta$name <- validate_colname(paste(name, "_", pattern, sep=""));
						df_meta$class_label <- "";
					
						meta_df_list$DATA[[df_meta$name]] <- df;
						meta_df_list$META <- rbind(meta_df_list$META, df_meta);
					} else {
						print("Pattern not found in data.");
					}
				}
				
			} else {
				label_counts <- list_count(data[, class_label]);
				print(paste(name, " has ", nrow(label_counts), " class labels", sep=""));
				label_n <- nrow(label_counts);
				
				# Go through individual label
				for(j in 1 : label_n)
				{
					label <- label_counts$LABEL[j];
					
					
					print(paste("Label: ", label, sep=""));
					
					selected_idx <- data[, class_label] == label;
					df <- data[selected_idx,];
					df_meta <- meta;
					df_meta$name <- validate_colname(paste(name, "_", label, sep=""));
					#df_meta$name <- gsub("[-]", "_", sub("[+]", "_", df_meta$name));
					df_meta$class_label <- "";
					
					meta_df_list$DATA[[df_meta$name]] <- df;
					meta_df_list$META <- rbind(meta_df_list$META, df_meta);
				}	
			}
			
			
			if(replace_original == T)
			{
				##########
				# bug
				##########
				meta_df_list$DATA[[name]] <- NULL;
				#meta_df_list$META <- meta_df_list$META[-i, ];
				meta_df_list$META <- meta_df_list$META[-which(meta_df_list$META$name == name), ];
			}
		}
	}
	
	print(paste("Name=", names(meta_df_list),sep=""))
	
	return(meta_df_list);
}

validate_colname <- function(colname) 
{
	return(sub("[-]", "_", sub("[+]", "_", colname)));
}



#
convert_meta_df_to_bin <- function(meta_df_list, bin_size, max_len)
{
	source("process_poss.R");
	
	metum <- meta_df_list[[META_LIST_META_LABEL]];
	datum <- meta_df_list[[META_LIST_DATA_LABEL]];
	data_n <- length(datum);
	
	bin_breaks <- seq(0, max_len, by=bin_size);
	bin_df <- data.frame();
	
	

	for(i in 1 : data_n)
	{
		meta <- metum[i, ];

		name <- as.character(meta$name);
		chrom_label <- as.character(meta$chrom_label);
		spos_label <- as.character(meta$spos_label);
		epos_label <- as.character(meta$epos_label);
		datatype <- as.character(meta$datatype);
		class_label <- as.character(meta$class_label);
		value_label <- as.character(meta$value_label);
	
		
		data <- datum[[name]];
		sposs <- as.numeric(data[, spos_label]);
		eposs <- as.numeric(data[, epos_label]);
		
		# Check if the chrom_label contains only one single value
		if(length(as.character(unique(data[, chrom_label]))) == 1 && nrow(data) > 0)
		{
			if(nchar(value_label) == 0)
			{
				if(datatype == "range")
				{
					print(paste("Processing range data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					poss <- process_poss_c(sposs, eposs);
					
				} else if(datatype == "mixed") {
					
					print(paste("Processing mixed data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					
					point_poss <- data[(eposs - sposs) <= 1, spos_label];
					
					print(paste("Number of point objects=", length(point_poss), sep=""))
					
					range_item_idx <- (eposs - sposs) > 1;
					range_poss <- process_poss_c(sposs[range_item_idx], eposs[range_item_idx]);
					
					poss <- c(point_poss, range_poss);
					
				} else if(datatype == "point") {
					
					print(paste("Processing point data ", name, " length=", nrow(data), sep=""));
					poss <- sposs;
					if(spos_label != epos_label)
					{
						poss <- c(poss, eposs);
					}	
				}
				
				if(dim(bin_df)[2] == 0)
				{
					bin_df <- data.frame(x_____X_____=as.numeric(table(cut(poss, breaks=bin_breaks, ordered_result=T))));
				} else {
					bin_df <- cbind(bin_df, data.frame(x_____X_____=as.numeric(table(cut(poss, breaks=bin_breaks, ordered_result=T)))));
				}
				colnames(bin_df)[which(colnames(bin_df) == "x_____X_____")] <- name;
				
			} else {  # This is a valued data
				
				data_values <- as.numeric(data[, value_label]);
						
				if(datatype == "range")
				{
					print(paste("Processing valued range data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					valued_poss <- process_poss_vals_c(sposs, eposs, data_values, max_len);
				
				} else if(datatype == "mixed") {
					
					print(paste("Processing valued mixed data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					
					range_item_idx <- (sposs - eposs) > 1;
					
					valued_poss <- process_poss_vals_c(sposs[range_item_idx], eposs[range_item_idx], data_values[range_item_idx], max_len);
					
					point_item_idx <- (eposs - sposs) <= 1;
					valued_poss[point_item_idx] <- data_values[point_item_idx];				
					
					
				} else if(datatype == "point") {
					print(paste("Processing valued point data ", name, " length=", nrow(data), sep=""));
					
					valued_poss <- rep(0.0, max_len);
					valued_poss[sposs] <- data_values;
					if(spos_label != epos_label)
					{
						valued_poss[eposs] <- data_values;
					}
				}
				
				
				valued_list <- c()
				#poss_idx <- cut(poss, breaks=bin_breaks, labels=F, ordered_result=T);
				for(i in 1 : (length(bin_breaks) - 1)) {
					bin_spos <- as.numeric(bin_breaks[i]);
					bin_epos <- as.numeric(bin_breaks[i + 1]);
					# print(paste("bin_spos=", bin_spos, ", bin_epos=", bin_epos, " mean=", mean(valued_poss[bin_spos:bin_epos]), sep=""))
					valued_list <- c(valued_list, mean(valued_poss[bin_spos:bin_epos]));
				}
				
				if(dim(bin_df)[2] == 0)
				{
					bin_df <- data.frame(x_____X_____=valued_list);
				} else {
					bin_df <- cbind(bin_df, data.frame(x_____X_____=valued_list));
				}
				colnames(bin_df)[which(colnames(bin_df) == "x_____X_____")] <- name;
			}

			
		} else {
			print(paste(name, ": values at the column,", chrom_label, "are inconsistent, data skipped.",  sep=" "));
		}
		
	}
	rownames(bin_df) <- bin_breaks[1:(length(bin_breaks)-1)];
	
	return(bin_df);
}



select_subset <- function(meta_df_list, meta_label, values, fixed=T, ignore_case=T)
{
	meta_df_list_subset <- list();
	
	
	#check if the meta label is available
	if(meta_label %in% colnames(meta_df_list$META))
	{
		if(fixed == T)
		{
			selected_idx <- which(meta_df_list$META[, meta_label] %in% values);
			
			meta_df_list_subset[[META_LIST_META_LABEL]] <- meta_df_list$META[selected_idx,];
			meta_df_list_subset[[META_LIST_DATA_LABEL]] <- meta_df_list$DATA[selected_idx];
		}
	} else {
		print(paste("Meta-Label, ", meta_label, " does not exist.", sep=""))
	}
	
	return(meta_df_list_subset);
}




convert_poss_list_to_bin <- function(poss_list, bin_size, max_len=integer(0))
{
	
}




convert_poss_list_to_bin <- function(poss_list, bin_size, max_len=integer(0))
{
	
}


# 
convert_poss_list_to_bin <- function(poss_list, bin_size, max_len=integer(0))
{
	data_n <- length(poss_list);
	
	print(paste("Begin converting with parameter: bin_size=", bin_size, " max-len=", max_len, sep=""));
	
	if(length(max_len) == 0)
	{
		max_len <- max(sapply(poss_list, max)) + bin_size;
	}
	
	print(paste("max=", max_len, " bin_size=", bin_size, sep=""))
	bin_breaks <- seq(0, max_len, by=bin_size);

	
	print(length(bin_breaks));
	
	
	#bin_list <- list()
	bin_df <- data.frame();
	for(i in 1 : data_n)
	{
		name <- names(poss_list)[i];
		print(paste("Converting ", name, " to a poss table...", sep=""));
		if(dim(bin_df)[2] == 0)
		{
			bin_df <- data.frame(x_____X_____=as.numeric(table(cut(poss_list[[i]], breaks=bin_breaks, ordered_result=T))));
		} else {
			bin_df <- cbind(bin_df, data.frame(x_____X_____=as.numeric(table(cut(poss_list[[i]], breaks=bin_breaks, ordered_result=T)))));
		}
		colnames(bin_df)[which(colnames(bin_df) == "x_____X_____")] <- name;
	}
	
	rownames(bin_df) <- format(bin_breaks[-length(bin_breaks)], scientific=F);
	
	return(bin_df);
}



# Provide a list of range objects, this function will cut data into corresponding range object
convert_df_to_poss_by_ranges <- function(meta_df_list, range_df, max_len, included_list=character(0))
{	
	source("process_poss.R");
	
	metum <- meta_df_list[[META_LIST_META_LABEL]];
	datum <- meta_df_list[[META_LIST_DATA_LABEL]];
	
	print(paste("Number of range object=", nrow(range_df), sep=""));
	
	
	labels <- names(datum);
	data_n <- length(labels);
	
	poss_list <- list();
	
	# Go through all the labels
	for(label in labels)
	{
		print(paste("Processing ", label, "...", sep=""));
		
		meta <- metum[metum$name == label,];
		data <- datum[[label]];
		
		# No data with a match label
		if(dim(data)[1] == 0)
			next;
		
		print(paste("size=", dim(meta)[2], sep=""))
		
		if(meta$datatype == "range")
		{
			values <- rep(1.0, nrow(data));
			if(nchar(meta$value_label) != 0) # It is a valued dataset
			{
				print(paste(">Valued range object identified.", sep=""));
				values <- data[,meta$value_label];
				#bins <- process_poss_vals_to_ranges(data[,meta$spos_label], data[,meta$epos_label], data[,meta$value_label], range_df, max_len);
				#poss <- process_poss_vals_c(data[,meta$spos_label], data[,meta$epos_label], data[,meta$value_label], max_len);
			}
			
			bins <- process_poss_vals_to_ranges(data[,meta$spos_label], data[,meta$epos_label], values, range_df, max_len);
			
			#poss_list[[label]] <- bins / (range_df$END - range_df$START);
			poss_list[[label]] <- bins;
		} else if(meta$datatype == "point") {
			values <- rep(1.0, nrow(data));
			if(nchar(meta$value_label) != 0) # It is a valued dataset
			{
				print(paste(">Valued range object identified.", sep=""));
				values <- data[,meta$value_label];
				#bins <- process_poss_vals_to_ranges(data[,meta$spos_label], data[,meta$epos_label], data[,meta$value_label], range_df, max_len);
				#poss <- process_poss_vals_c(data[,meta$spos_label], data[,meta$epos_label], data[,meta$value_label], max_len);
			}
			
			poss <- data[,meta$spos_label];
			bins <- process_poss_vals_to_ranges(poss, poss+1, values, range_df, max_len);
			
			poss_list[[label]] <- bins;
		}
		
	}
	
	return(poss_list);
}




convert_df_to_poss <- function(meta_list_df, merge_duplicates=F)
{
	return(convert_range_df_to_poss(meta_list_df, merge_duplicates));
}


META_LIST_META_LABEL = "META";
META_LIST_DATA_LABEL = "DATA";
# 
convert_range_df_to_poss <- function(meta_list_df, merge_duplicates=F)
{
	source("process_poss.R");
	
	metum <- meta_list_df[[META_LIST_META_LABEL]];
	datum <- meta_list_df[[META_LIST_DATA_LABEL]];
	
	
	poss_list = list();
	
	data_n <- length(datum);
	for(i in 1 : data_n)
	{
		meta <- metum[i, ];
		data <- datum[[i]];
		
		name <- as.character(meta$name);
		chrom_label <- as.character(meta$chrom_label);
		spos_label <- as.character(meta$spos_label);
		epos_label <- as.character(meta$epos_label);
		datatype <- as.character(meta$datatype);
		
		sposs <- as.numeric(data[, spos_label]);
		eposs <- as.numeric(data[, epos_label]);
		
		# Check if the chrom_label contains only one single value
		if(length(as.character(unique(data[, chrom_label]))) == 1)
		{
			if(datatype == "range")
			{
				print(paste("Processing range data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
				poss <- process_poss_c(sposs, eposs);
				print(paste("Number of positions generated: ", length(poss)));
			} else if(datatype == "mixed") {
				print(paste("Processing mixed data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
				point_poss <- data[(eposs - sposs) <= 1, spos_label];
				
				range_item_idx <- (eposs - sposs) > 1;
				range_poss <- process_poss_c(as.numeric(sposs[range_item_idx]), as.numeric(eposs[range_item_idx]));
				
				poss <- c(point_poss, range_poss);
				
			} else if(datatype == "point") {
				print(paste("Processing point data ", name, " length=", nrow(data), sep=""));
				poss <- sposs;
				if(spos_label != epos_label)
				{
					poss <- c(poss, eposs);
				}
			}
			poss_list[[name]] <- sort(poss);
		} else {
			print(paste(name, ": values at the column,", chrom_label, "are inconsistent, data skipped.",  sep=" "));
		}
	}
	
	return(poss_list);
}


get_ordered_label_meta_info_by_group_id <- function(meta_df_list) {
	group_labels <- sort(as.character(unique(meta_df_list$META$group)));
	

	ordered_labels <- character(0);
	ordered_group_labels <- character(0);
	
	for(i in 1 : length(group_labels)) {
		group_label <- group_labels[i];
		
		selected_labels <- meta_df_list$META$name[meta_df_list$META$group == group_label];
		selected_labels <- sort(selected_labels);
		
		ordered_labels <- c(ordered_labels, selected_labels);
		ordered_group_labels <- c(ordered_group_labels, rep(group_label, length(selected_labels)))
	}
	
	label_meta_info <- data.frame(LABEL=ordered_labels, GROUP=ordered_group_labels, stringsAsFactors=F);
	
	return(label_meta_info)
}


get_ordered_label_by_group_id <- function(meta_df_list) {
	group_labels <- sort(as.character(unique(meta_df_list$META$group)));
	
	
	ordered_labels <- character(0);
	
	for(i in 1 : length(group_labels)) {
		group_label <- group_labels[i];
		
		selected_labels <- meta_df_list$META$name[meta_df_list$META$group == group_label];
		selected_labels <- sort(selected_labels);
		
		ordered_labels <- c(ordered_labels, selected_labels);
	}
	
	return(ordered_labels)
}

# Construct a list of df from a list of filenames
construct_df <- function(fns_df, chr_id=character(0), quote="", header=T, sep="\t", comment_char="#", stringsAsFactors=F) {
	construct_meta_df(fns_df=fns_df, chr_id=chr_id, quote=quote, header=header, sep=sep, comment_char=comment_char, stringsAsFactors=stringsAsFactors);
}

construct_meta_df <- function(fns_df, chr_id=character(0), quote="", header=T, sep="\t", comment_char="#", stringsAsFactors=F) {
	fns_n <- nrow(fns_df);
	
	list_df = list();
	# data_label -> name
	list_df[["META"]] <- data.frame(group=character(0), name=character(0), datatype=character(0), class_label=character(0), class_label_pattern=character(0), value_label=character(0), chrom_label=character(0), spos_label=character(0), epos_label=character(0), description=character(0) , stringsAsFactors=stringsAsFactors);
	list_df[["DATA"]] <- list();
	
	
	for(i in 1 : fns_n)
	{
		fn <- fns_df$filename[i];
		chrom_label <- fns_df$chrom_label[i];
		spos_label <- fns_df$spos_label[i];
		epos_label <- fns_df$epos_label[i];
		
		name <- fns_df$name[i];
		group_label <- fns_df$group_label[i];
		class_label <- fns_df$class_label[i];
		class_label_pattern <- fns_df$class_label_pattern[i];
		value_label <- fns_df$value_label[i];
		data_type <- fns_df$data_type[i];
		desc <- fns_df$description[i];
		

		print(paste("construct_meta_df=>Scanning file ", fn, sep=""));
		
		
		list_df[["META"]] <- rbind(list_df[["META"]], data.frame(group=group_label, name=name, datatype=data_type, class_label=class_label, class_label_pattern=class_label_pattern, value_label=value_label, chrom_label=chrom_label, spos_label=spos_label, epos_label=epos_label, description=desc, stringsAsFactors=stringsAsFactors));
		
		
		if(length(chr_id)!=0)
		{
			df <- retrieve_subset(fn=fn, colname=chrom_label, selected_value=chr_id, comment_char="#", stringsAsFactors=stringsAsFactors);			
		} else {
			df <- read.table(fn, quote=quote, header=header, comment.char=comment_char, stringsAsFactors=stringsAsFactors);
		}

		if(nrow(df) > 0)
		{
			list_df[["DATA"]][[name]] <- df[sort(df[,spos_label], index.return=T)$ix,];
		}
	}
	
	return(list_df);
}


# Given a filename, and a given value of a column will be used to filter data of the file
retrieve_subset <- function(fn, colname, selected_value, quote="", header=T, sep="\t" , comment_char="#", stringsAsFactors=F)
{
	# path to python get_subset script
	path_to_script <- "/data/bcnskaa/UCSC/annotations/get_subset.py";
	filename_len <- 20;

	
	# generate a random file
	tmp_fn <- paste(getwd(), "/", paste(sample(c(LETTERS, letters, seq(0, 9)), filename_len, replace=T), sep="", collapse=""), ".temp", sep="");
	
	
	cmd <- paste("python3.3", path_to_script, fn, tmp_fn, colname, selected_value, sep=" ");
	print(paste("Extracting subset from ", fn, " (colname: ", colname, ", value: ", selected_value, ")...", sep=""));
	system(cmd);
	
	tmp_df <- read.table(tmp_fn, header=header, quote=quote, sep=sep, comment.char=comment_char, stringsAsFactors=stringsAsFactors);
	file.remove(tmp_fn);
	return(tmp_df);
}

#
import_df_from_multiple_files <- function(filenames, df_ids, colname=character(0), selected_value=character(0), sep="\t", header=T, quote="", comment_char="#", stringsAsFactors=F) 
{

	if(length(filenames) != length(df_ids))
	{
		print("The number of filenames and IDs are not consistent, abort now.");
		return;
	}
	
	df_list = list();
	
	for(i in 1 : length(filenames))
	{
		fn <- filenames[i];
		id <- df_ids[i];
			
		print(paste("Reading from", fn, "..."))
		
		if(file.exists(fn))
		{
			if(length(colname) == 1 & length(selected_value) == 1) {
				df <- retrieve_subset(fn, colname, selected_value, sep=sep, header=header, quote=quote, comment_char=comment_char, stringsAsFactors=stringsAsFactors);
			} else {
				df <- read.table(fn, sep=sep, header=header, quote=quote, comment.char=comment_char, stringsAsFactors=stringsAsFactors);
			}
			df_list[[id]] <- df;
		} else {
			print(paste("File, ", fn, " does not exist, skipped.", sep=""));
		}
	
	}
	return(df_list);
}


# df <- unlist_multiple_df(df_list)
unlist_multiple_df <- function(df_list, group_col_label="GROUP") {
	ids <- names(df_list);
	
	combined_df <- character(0);
	for(id in ids)
	{
		print(paste("Processing ", id, sep=""));
		df <- df_list[[id]];
		if(dim(df)[1] > 0)
		{
			df <- cbind(df, XXXXXXXX_XXXXXXX=rep(id, nrow(df)));
			
			if(is.null(dim(combined_df)))
			{
				combined_df <- df;
			} else {
				combined_df <- rbind(combined_df, df);
			}	
		}
	}
	colnames(combined_df)[which(colnames(combined_df) == "XXXXXXXX_XXXXXXX")] <- group_col_label;
	
	return(combined_df);
}



batch_files_conversion <- function() 
{
	path <- "."
	infile_pattern <- ".bb"
	outfile_pattern <- ".bed"
	converter_cmd <- "bigBedToBed"
	list_infiles <- list.files(path, pattern=infile_pattern);
	list_outfiles <- sub(infile_pattern, outfile_pattern, list_infiles);
	cmds <- paste(converter_cmd, list_infiles, list_outfiles);
	for(cmd in cmds) {
		print(paste("Running command:", cmd))
		system(cmd);
	}
	
	# dhc datasets specific
	cmds <- paste("cp", "dhcHumDerDenAncAllHighFreq.tbl", sub(".bed", ".tbl", list_outfiles));
	for(cmd in cmds) {
		print(paste("Running command:", cmd))
		system(cmd);
	}
	
	cmds <- paste("cp", list_outfiles, sub(".bed", ".txt", list_outfiles));
	for(cmd in cmds) {
		print(paste("Running command:", cmd))
		system(cmd);
	}
}



##########################
# Plot-related functions
##########################

# melt plot_mtx into linearized plot_data
create_plot_data_from_plot_mtx <- function(plot_mtx)
{
	row_n <- nrow(plot_mtx);
	col_n <- ncol(plot_mtx);
	
	col_names <- colnames(plot_mtx);
	row_names <- rownames(plot_mtx);
	
	plot_data <- data.frame(X=integer(0), Y=integer(0), Z=integer(0), stringsAsFactors=F);
	for(i in 1 : row_n)
	{
		plot_data <- rbind(plot_data, data.frame(X=col_names, Y=rep(row_names[i], col_n), Z=as.numeric(plot_mtx[i, ]), stringsAsFactors=F))
	}
	
	return(plot_data);
}


# Given an ordered list of labels, this function will sort the column order of the plot_data according to the given order list. 
sort_plot_data_by_group <- function(plot_data, ordered_labels)
{
	ordered_label_n <- length(ordered_labels);
	
	ordered_plot_data <- character(0);
	for(i in 1 : ordered_label_n) {
		label <- ordered_labels[i];
		
		
		coldata <- tryCatch(plot_data[, label], error=function(e) {cflag = F; print(paste("Label ", label, " does not existed in the plot data.", sep=""))});
		
		if(length(coldata) == 1)
		{
			next;
		}
		
		if(length(ordered_plot_data) == 0)
		{
			ordered_plot_data <- data.frame(X_______X_____X = coldata);
		} else {
			ordered_plot_data <- cbind(ordered_plot_data, data.frame(X_______X_____X = coldata));
		}
		
		
		colnames(ordered_plot_data)[colnames(ordered_plot_data) == "X_______X_____X"] <- label;
	}
	
	rownames(ordered_plot_data) <- rownames(plot_data);
	
	return(ordered_plot_data);
}


create_plot_data <- function(bin_list)
{
	data_n <- length(bin_list);
	labels <- names(bin_list);
	row_n <- as.numeric(unique(sapply(bin_list, function(v) length(v))))
	
	plot_data <- data.frame(matrix(0, row_n, data_n));
	colnames(plot_data) <- labels;
	
	
	for(label in labels) {
		print(label)
		bin_data <- bin_list[[label]];
		plot_data[, label] <- bin_data;
	}
	
	rownames(plot_data) <- rownames(bin_list);
	
	return(plot_data);
}



# focal_colname: the name of a column, on which the row of the plot data having a value on the focal column less than trim threshold will be trimmed
trim_plot_data <- function(plot_data, focal_colname=character(0), trim_threshold=0)
{
	row_n = dim(plot_data)[1];
	idx_to_be_trimmed <- c();
	
	trim_n <- 0;
	
	if(length(focal_colname) == 0)
	{
		for(i in 1 : row_n) 
		{
			total <- sum(plot_data[i,]);
			print(total)
			if(total <= trim_threshold)
			{
				idx_to_be_trimmed <- c(idx_to_be_trimmed, i);
				#plot_data <- plot_data[-i,];
				trim_n = trim_n + 1;
			}
		}
		
		plot_data <- plot_data[-idx_to_be_trimmed, ];
	} else {
		#idx_to_be_trimmed <- plot_data[,focal_colname] <= trim_threshold;		
		plot_data <- plot_data[plot_data[,focal_colname] > trim_threshold, ];
		
		trim_n <- row_n - dim(plot_data)[1];
	}
	
	
	
	print(paste(trim_n, " datapoint trimmed.", sep=""));
	
	return(plot_data);
}




multiplot <- function(..., plotlist=NULL, cols) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# Make the panel
	plotCols = cols                          # Number of columns of plots
	plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
	
	# Set up the page
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
	vplayout <- function(x, y)
		viewport(layout.pos.row = x, layout.pos.col = y)
	
	# Make each plot, in the correct location
	for (i in 1:numPlots) {
		curRow = ceiling(i/plotCols)
		curCol = (i-1) %% plotCols + 1
		print(plots[[i]], vp = vplayout(curRow, curCol ))
	}
	
}
