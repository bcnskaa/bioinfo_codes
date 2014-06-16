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

# Create an empty list of feature poss data




#process_DF_test_code("chr19", 1500000);
#process_DF_test_code("chr18", 1500000);
#process_DF_test_code("chrX", 1500000);
#process_DF_test_code("chr5", 1500000);
test_code_process_DF <- function(chr_id, bin_size) 
{
	
# Construct a list of df from a list of filenames
# fns_df: filename, group_label, data_label, data_type[point|range|mixed], chrom_label, class_label, description
gc();
	
source("process_DF.R");

fns <- read.table("/data/bcnskaa/UCSC/annotations/ucsc_annotations_for_R.lst", header=T, sep="\t", quote="", stringsAsFactors=F);
fns$filename <- paste("/data/bcnskaa/UCSC/annotations", "/", fns$filename, sep="");
#chr_id <- "chr13";
#bin_size <- 150000;
	
if(exists(commandArgs(TRUE)[1]) || !is.na(commandArgs(TRUE)[1]))
{
	chr_id <- commandArgs(TRUE)[1];
}	

if(exists(commandArgs(TRUE)[2]) || !is.na(commandArgs(TRUE)[2]))
{
	bin_size <- commandArgs(TRUE)[2];
}

print(paste("Processing chr_id=", chr_id, " and bin_size=", bin_size, sep=""))


source("chromosome_stats.R");
chr_len <- get_chr_len(chr_id);


meta_df_list <- construct_meta_df(fns, chr_id);


gc();

#poss_list <- convert_df_to_poss(meta_df_list);
#bin_list <- convert_poss_list_to_bin(poss_list, bin_size)
#bin_list2 <- convert_meta_df_to_bin(meta_df_list, bin_size);
	
factorized_meta_df_list <- factorize_meta_df(meta_df_list);
#source("process_DF.R");

fns <- read.table("/data/bcnskaa/UCSC/annotations/ucsc_annotations_for_R.lst", header=T, sep="\t", quote="", stringsAsFactors=F);
fns$filename <- paste("/data/bcnskaa/UCSC/annotations", "/", fns$filename, sep="");
#chr_id <- "chr13";
#bin_size <- 150000;

meta_df_list <- construct_meta_df(fns, chr_id);


gc();

#poss_list <- convert_df_to_poss(meta_df_list);
#bin_list <- convert_poss_list_to_bin(poss_list, bin_size)
#bin_list2 <- convert_meta_df_to_bin(meta_df_list, bin_size);

factorized_meta_df_list <- factorize_meta_df(meta_df_list);
factorized_meta_df_label_info <- get_ordered_label_meta_info_by_group_id(factorized_meta_df_list);

gc();

bin_list3 <- convert_meta_df_to_bin(factorized_meta_df_list, bin_size=bin_size, max_len=chr_len);


gc();

source("process_fasta.R")
# fasta sequence
fasta_fn <- paste("/data/bcnskaa/UCSC/", chr_id, ".fa", sep="");
N_poss <- process_fasta_c(fasta_fn, chr_id, "N");

gc();

N_poss_bin <- convert_poss_list_to_bin(list(N=N_poss), bin_size=bin_size, max_len=chr_len)
N_poss_bin <- (bin_size - N_poss_bin) / bin_size;
# Add the N_poss_bin and use it to filter bin having N-characters
bin_list3[["N"]] <- N_poss_bin$N;

gc();


source("plot_scatter_mtx.R");
#plot_scatter_mtx(bin_list3, png_outfile="plot_scatter_matrix.png", scatter_size=2);	
#plot_scatter_mtx(bin_list3, png_outfile="plot_scatter_matrix.png");

plot_data <- create_plot_data(bin_list3);
plot_data <- trim_plot_data(plot_data, focal_colname="N", trim_threshold=0.95);
tryCatch(plot_data <- trim_plot_data(plot_data, focal_colname="DecodeRecombRate_SexAverage"), error=function(e) print(e));


# N_poss column will not be included in the sorted plot_data
plot_data <- sort_plot_data_by_group(plot_data, factorized_meta_df_label_info$LABEL);

gc();

pdf_outfn <- paste("plot_scatter_matrix.", chr_id, ".", bin_size, ".all.pdf", sep="");
plot_scatter_mtx2(plot_data, pdf_outfile=pdf_outfn, scatter_size=150, plot_title=paste("Correlation plot of ", chr_id, " (", bin_size," bp)", sep=""));
system(paste("convert", pdf_outfn, sub(".pdf", ".png", pdf_outfn),sep=" "))
gc();



# Prepare a correlation matrix
cor_mtx <- cor(plot_data)
cor_plot_data <- create_plot_data_from_plot_mtx(cor_mtx)


source("plot_heat_map.R")

pdf_fn <- paste("correlation_map.", chr_id, ".pdf", sep="")
#plot_heatmap_with_group_info(cor_plot_data, plot_title=paste("Correlation Matrix: ", chr_id , sep=""), pdf_width=10, pdf_height=10, pdf_fn=pdf_fn, show_value=F);
plot_heatmap(cor_plot_data, xlabel_meta_info=factorized_meta_df_label_info, ylabel_meta_info=factorized_meta_df_label_info, plot_title=paste("Correlation Matrix: ", chr_id, " (Bin Size=", bin_size, " bp)", sep=""), pdf_width=12, pdf_height=12, pdf_fn=pdf_fn, show_value=F);
system(paste("convert", pdf_fn, sub(".pdf", ".png", pdf_fn)))

# Prepare a correlation matrix
pdf_fn <- paste("correlation_map.", chr_id, ".valued.pdf", sep="")
cor_plot_data$Z <- round(cor_plot_data$Z, digits=2)
#plot_heatmap_with_group_info(cor_plot_data, plot_title=paste("Correlation Matrix: ", chr_id , sep=""), pdf_width=20, pdf_height=20, pdf_fn=pdf_fn, show_value=T);
plot_heatmap(cor_plot_data, xlabel_meta_info=factorized_meta_df_label_info, ylabel_meta_info=factorized_meta_df_label_info, plot_title=paste("Correlation Matrix: ", chr_id, " (Bin Size=", bin_size, " bp)", sep=""), pdf_width=30, pdf_height=30, pdf_fn=pdf_fn, show_value=T);
system(paste("convert", pdf_fn, sub(".pdf", ".png", pdf_fn)))


pdf_fn <- paste("correlation_map.", chr_id, ".group-valued.pdf", sep="")
#plot_heatmap_with_group_info(cor_plot_data, plot_title=paste("Correlation Matrix: ", chr_id , sep=""), pdf_width=20, pdf_height=20, pdf_fn=pdf_fn, show_value=T);
plot_heatmap_with_group_info(cor_plot_data, xlabel_meta_info=factorized_meta_df_label_info, ylabel_meta_info=factorized_meta_df_label_info, show_value=T, plot_title=paste("Correlation plot of ", chr_id, "(bin size=", bin_size, " bp)", sep=""), title_x=-5, title_y=1.05, pdf_fn=pdf_fn, pdf_width=23, pdf_height=23)
system(paste("convert", pdf_fn, sub(".pdf", ".png", pdf_fn)))


labels <- colnames(plot_data)
data_n <- dim(plot_data)[2]
cor_pval_mtx <- matrix(0.0, data_n, data_n)
colnames(cor_pval_mtx) <- labels;
rownames(cor_pval_mtx) <- labels;


cor_df <- data.frame(labels=character(0), xlabel=character(0), ylabel=character(0), cor=numeric(0), pval=numeric(0), conf_low=numeric(0), conf_high=numeric(0), stringsAsFactors=F);
for(i in 1 : data_n)
{
	xlabel <- colnames(cor_pval_mtx)[i];
	for(j in 1 : data_n)
	{
		ylabel <- rownames(cor_pval_mtx)[j];
		res <- cor.test(plot_data[,which(colnames(plot_data) == xlabel)], plot_data[,which(colnames(plot_data) == ylabel)])
		#print(paste(xlabel, ylabel, "=", res$p.value, sep=""))
		cor_pval_mtx[i, j] <- res$p.value;
		
		cor_df <- rbind(cor_df, data.frame(labels=paste(xlabel, ":", ylabel, sep=""), xlabel=xlabel, ylabel=ylabel, cor=res$estimate, pval=res$p.value, conf_low=res$conf.int[1], conf_high=res$conf.int[2], stringsAsFactors=F));
	
	}
}
cor_pval_plot_data <- create_plot_data_from_plot_mtx(cor_pval_mtx)
cor_pval_plot_data$Z <- round(cor_pval_plot_data$Z, digits=3)

write.table(cor_df, paste("correlation_pvalue.", chr_id, ".", bin_size, ".txt", sep=""), row.names=F, quote=F, sep="\t")


pdf_fn <- paste("correlation_map.", chr_id, ".p_values.pdf", sep="")
#plot_heatmap_with_group_info(cor_plot_data, plot_title=paste("Correlation Matrix: ", chr_id , sep=""), pdf_width=20, pdf_height=20, pdf_fn=pdf_fn, show_value=T);
plot_heatmap_with_group_info(cor_pval_plot_data, xlabel_meta_info=factorized_meta_df_label_info, ylabel_meta_info=factorized_meta_df_label_info, show_value=T, plot_title=paste("Significance of Correlation plot of ", chr_id, " (bin size=", bin_size, " bp)", sep=""), title_x=-5, title_y=1.05, pdf_fn=pdf_fn, pdf_width=23, pdf_height=23, scale_colors=c("blue", "yellow"), scale_limits=c(0, 1))
system(paste("convert", pdf_fn, sub(".pdf", ".png", pdf_fn)))



rm(list=ls());
}


##########
# poss feature table
##########
#feature_poss[["COSMIC"]] <- cosmic_chr_bp;
#feature_poss[["1KGP_INDELS"]] <- indels_frq$POS;
#feature_poss[["1KGP_SNVS"]] <- snv_frq$POS;
#feature_poss[["DGV_BP"]] <- dgv_chr_bp;
#feature_poss[["DGV_BP_UNIQUE"]] <- dgv_chr_bp_unique;
#feature_poss[["DGV_BP_NR"]] <- dgv_chr_bp_nonredundant;
#feature_poss[["DGV_BP_GAIN"]] <- dgv_chr_bp_gain;
#feature_poss[["DGV_BP_GAINLOSS"]] <- dgv_chr_bp_gainloss;
#feature_poss[["DGV_BP_COMPLEX"]] <- dgv_chr_bp_complex;
#feature_poss[["DGV_BP_GAIN"]] <- dgv_chr_bp_gain;
#feature_poss[["DGV_BP_CNV"]] <- dgv_chr_bp_cnv;
#feature_poss[["DGV_BP_INSERTION"]] <- dgv_chr_bp_insertion;
#feature_poss[["DGV_BP_DELETION"]] <- dgv_chr_bp_deletion;

#########
#
# Construct feature_set based on unique labels
#
#########
# write.table(melt_isca, "/data/bcnskaa/UCSC/annotations/isca/iscaAllTissues.txt2", sep="\t", quote=F, row.names=F)

# Input: a list of dataframes with identical number and name of columns
# Output: a dataframe containing all elements of dataframes and a new column specifies name of the dataframe from which element is draw from
# Argument: df_list=a list of dataframes with identical number and name of columns
#			new_class_label=the name of new column listing the source dataframes  
melt_df <- function(df_list, new_class_label="class") {
	if(length(unique(sapply(df_list, function(v) dim(v)[2]))) != 1)
	{
		print(paste("Input list contains dataframes of different columns' number/name, unable to melt the list."))
		return
	}
	
	class_labels <- names(df_list);
		
	melted_df <- data.frame();
	
	# 
	for(label in class_labels) {
		print(paste("Processing ", label, "...", sep=""));
		
		cur_df <- df_list[[label]];
		#cur_df <- cbind(cur_df, class=rep(label, nrow(cur_df)));
		cur_df <- cbind(cur_df, X______X_____X=rep(label, nrow(cur_df)));
		if(dim(melted_df)[1] == 0)
		{
			melted_df <- cur_df;
		} else {
			melted_df <- rbind(melted_df, cur_df);
		}
	}
	#melted_df$X______X_____X <- as.character(melted_df$X______X_____X);
	
	#melted_df$class <- as.character(melted_df$class);
	colnames(melted_df)[colnames(melted_df) == "X______X_____X"] <- new_class_label;
	
	return(melted_df);
}


trim_df <- function(bin_list) 
{

}

#
#
# 
process_df_with_label_vals <- function(dataframe, label_idx, label_vals, label_n_cutoff = 0)
{
	#unique_labels <- as.character(unique(df[, label_idx]));
	#counts <- table(dataframe[, label_idx])
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
	#unique_labels <- as.character(unique(df[, label_idx]));
	counts <- table(dataframe[, label_idx])
	select_labels <- names(counts[counts >= label_n_cutoff]);
	
	return(process_df_with_label_vals(dataframe, label_idx, select_labels, labal_n_cutoff));
	
#	label_n <- length(select_labels);
#	# Create a new data.frame
#	df <- list();
#	
#	for(i in 1 : label_n)
#	{
#		current_label <- select_labels[i];
#		df[[current_label]] <- dataframe[which(dataframe[,label_idx] == current_label),];
#	}
#	
#	return(df);
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



# test <- factorize_meta_df(meta_df_list)
factorize_meta_df <- function(meta_df_list, replace_original=T, enable_class_label_pattern=F, fixed_class_label_pattern=F)
{
	source("list_count.R");
	
	metum <- meta_df_list[[META_LIST_META_LABEL]];
	datum <- meta_df_list[[META_LIST_DATA_LABEL]];
	####################
	# bug
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
						#df_meta$name <- paste(name, "_", pattern, sep="");
						#df_meta$name <- gsub("[-]", "_", sub("[+]", "_", df_meta$name));
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


# test <- meta_df_list$DATA$decodeRecombRate
# valued_poss <- process_poss_vals_c(test$chromStart, test$chromEnd, as.numeric(test$Rate),  max_len)
# xyplot(valued_list ~ bin_breaks[-length(bin_breaks)], data.frame(bin_breaks[-length(bin_breaks)], valued_list))
# xyplot(NucleaseAssessibleSitesPos ~ NucleaseAssessibleSitesNeg, data.frame(NucleaseAssessibleSitesPos=bin_list2$NucleaseAssessibleSitesPos, NucleaseAssessibleSitesNeg=bin_list2$NucleaseAssessibleSitesNeg))

#convert_meta_df_to_bin <- function(meta_df_list, bin_size, max_len=integer(0))
convert_meta_df_to_bin <- function(meta_df_list, bin_size, max_len)
{
	source("process_poss.R");
	
	metum <- meta_df_list[[META_LIST_META_LABEL]];
	datum <- meta_df_list[[META_LIST_DATA_LABEL]];
	data_n <- length(datum);
	
	
#	if(length(max_len) == 0)
#	{
#		max_len <- 0;
#		for(i in 1 : length(metum)) {
#			meta <- metum[i, ];
#			data <- datum[[i]];
#			
#			#print(meta$name)
#			#print(meta$epos_label)
#			# sapply(seq(1:length(meta_df_list$META$name)), function(i) max(meta_df_list$DATA[[i]][, meta_df_list$META$epos_label[i]]))
#			val <- max(data[,meta$epos_label]);
#			if(val > max_len)
#			{
#				max_len <- val;
#			}
#		}
#		max_len <- max_len + bin_size;
#	}
	
	#print(paste("Max length=", max_len,sep=""));
	
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
					#print(paste("Processing range data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					#poss <- process_poss_c(data[, spos_label], data[, epos_label]);
					
					print(paste("Processing range data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					poss <- process_poss_c(sposs, eposs);
					
				} else if(datatype == "mixed") {
					#print(paste("Processing mixed data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					#
					#point_poss <- data[(data[, epos_label] - data[, spos_label]) == 1, spos_label];
					#
					#range_item_idx <- (data[, epos_label] - data[, spos_label]) > 1;
					#range_poss <- process_poss_c(data[range_item_idx, spos_label], data[range_item_idx, epos_label]);
					#
					#poss <- c(point_poss, range_poss);
					
					
					print(paste("Processing mixed data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					
					point_poss <- data[(eposs - sposs) <= 1, spos_label];
					
					print(paste("Number of point objects=", length(point_poss), sep=""))
					
					range_item_idx <- (eposs - sposs) > 1;
					range_poss <- process_poss_c(sposs[range_item_idx], eposs[range_item_idx]);
					
					poss <- c(point_poss, range_poss);
					
				} else if(datatype == "point") {
					#print(paste("Processing point data ", name, " length=", nrow(data), sep=""));
					#poss <- data[, spos_label];
					#if(spos_label != epos_label)
					#{
					#	poss <- c(poss, data[, epos_label]);
					#}
					
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
					#print(paste("Processing valued mixed data ", name, " size=", nrow(data), "x", ncol(data), sep=""));
					#
					#range_item_idx <- (data[, epos_label] - data[, spos_label]) > 1;
					#
					#valued_poss <- process_poss_vals_c(as.numeric(data[range_item_idx, spos_label]), as.numeric(data[range_item_idx, epos_label]), as.numeric(data[range_item_idx, value_label]), max_len);
					#
					#point_item_idx <- (data[, epos_label] - data[, spos_label]) == 1;
					#valued_poss[point_item_idx] <- as.numeric((data[point_item_idx, value_label]));			
					
					
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
				
				#if(dim(bin_df)[2] == 0)
				#{
				#	bin_df <- data.frame(x_____X_____=as.numeric(table(cut(poss, breaks=bin_breaks, ordered_result=T))));
				#} else {
				#	bin_df <- cbind(bin_df, data.frame(x_____X_____=as.numeric(table(cut(poss, breaks=bin_breaks, ordered_result=T)))));
				#}
				#colnames(bin_df)[which(colnames(bin_df) == "x_____X_____")] <- name;
				
			}

			
		} else {
			print(paste(name, ": values at the column,", chrom_label, "are inconsistent, data skipped.",  sep=" "));
		}
		
	}
	#print(length(rownames(bin_df)));
	#print(length(bin_breaks));
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


# max_len <- 50000000
# poss_n <- 100000
# bin_size <- 250000
# poss <- sample(max_len, poss_n)
# poss <- sort(poss)
# ranges <- seq(1, max_len, by=bin_size)
# bin_list <- convert_poss_list_to_bin(poss_list, bin_size)
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
		#bin_list[[name]] <- table(cut(poss_list[[i]], breaks=bin_breaks, ordered_result=T));

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
			#} else { # It is not a valued dataset, just assign 1 as a default value
			#	print(paste(">Range object identified.", sep=""));
			#	bins <- process_poss_vals_to_ranges(data[,meta$spos_label], data[,meta$epos_label], rep(1.0, nrow(data)), range_df, max_len);
			#	#poss <- process_poss_c(data[,meta$spos_label], data[,meta$epos_label]);
			#}
			
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
			#poss <- c(data[,meta$spos_label], data[,meta$epos_label]);
			#values <- c(values, values);
			
			#bins <- process_poss_vals_to_ranges(data[,meta$spos_label], (data[,meta$spos_label] + 1), values, range_df, max_len);
			bins <- process_poss_vals_to_ranges(poss, poss+1, values, range_df, max_len);
			
			#poss_list[[label]] <- bins / (range_df$END - range_df$START);
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
# poss_list <- convert_df_to_poss(meta_df_list)
# merge_duplicate: merge the poss data if two lists have a same name
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

			#print(paste("Min pos=", min(data[, spos_label])));
			#print(paste("Max pos=", max(data[, spos_label])));
			#print(paste("Min pos=", min(data[, epos_label])));
			#print(paste("Max pos=", max(data[, epos_label])));
			#print(summary(data[, epos_label] - data[, spos_label]));
			#print(min(data[, epos_label] - data[, spos_label]));

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

# Fit distribution in R
# http://stats.stackexchange.com/questions/25568/estimating-the-distribution-from-data
# require(fitdistrplus)
#bin_list$refGene[bin_list$refGene > 0]
#tdat <- bin_list$refGene[bin_list$refGene > 0] - mean(bin_list$refGene[bin_list$refGene > 0])

#descdist(tdat, boot=100)



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
# fns_df: filename, group_label, data_label, data_type[point|range|mixed], chrom_label, class_label, description
# fns <- data.frame(filename=character(0), group_label=character(0), data_label=character(0), data_type=character(0), filtering_column=character(0), filtering_value=character(0), description=character(0));
# fns <- read.table("/data/bcnskaa/UCSC/annotations/ucsc_annotations_for_R.lst", header=T, sep="\t", stringsAsFactors=F)
# fns$filename <- paste("/data/bcnskaa/UCSC/annotations", "/", fns$filename, sep="")
# chr_id <- "chr19"
# res <- construct_df(fns, chr_id)
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
		#chrom_col_val <- fns_df$filtering_value[i];
		
		name <- fns_df$name[i];
		group_label <- fns_df$group_label[i];
		class_label <- fns_df$class_label[i];
		class_label_pattern <- fns_df$class_label_pattern[i];
		value_label <- fns_df$value_label[i];
		data_type <- fns_df$data_type[i];
		desc <- fns_df$description[i];
		

		print(paste("construct_meta_df=>Scanning file ", fn, sep=""));
		
		
		list_df[["META"]] <- rbind(list_df[["META"]], data.frame(group=group_label, name=name, datatype=data_type, class_label=class_label, class_label_pattern=class_label_pattern, value_label=value_label, chrom_label=chrom_label, spos_label=spos_label, epos_label=epos_label, description=desc, stringsAsFactors=stringsAsFactors));
		
		#print(paste("fn=", fn, ", chrom_label=", chrom_label, ", chr_id=", chr_id, sep=""));
		
		if(length(chr_id)!=0)
		{
			df <- retrieve_subset(fn=fn, colname=chrom_label, selected_value=chr_id, comment_char="#", stringsAsFactors=stringsAsFactors);			
		} else {
			df <- read.table(fn, quote=quote, header=header, comment.char=comment_char, stringsAsFactors=stringsAsFactors);
		}

		if(nrow(df) > 0)
		{
			#list_df[[name]] <- df;
			list_df[["DATA"]][[name]] <- df[sort(df[,spos_label], index.return=T)$ix,];
		}
	}
	
	return(list_df);
}


# Given a filename, and a given value of a column will be used to filter data of the file
retrieve_subset <- function(fn, colname, selected_value, quote="", header=T, sep="\t" , comment_char="#", stringsAsFactors=F) {
#retrieve_subset <- function(fn, colname, selected_value) {
#	sep="\t";
#	header=T;
#	quote="";

	# path to python get_subset script
	path_to_script <- "/data/bcnskaa/UCSC/annotations/get_subset.py";
	filename_len <- 20;

	
	# generate a random file
	tmp_fn <- paste(getwd(), "/", paste(sample(c(LETTERS, letters, seq(0, 9)), filename_len, replace=T), sep="", collapse=""), ".temp", sep="");
	
	
	cmd <- paste("python3.3", path_to_script, fn, tmp_fn, colname, selected_value, sep=" ");
	#print("fn:", fn);
	print(paste("Extracting subset from ", fn, " (colname: ", colname, ", value: ", selected_value, ")...", sep=""));
	system(cmd);
	
	tmp_df <- read.table(tmp_fn, header=header, quote=quote, sep=sep, comment.char=comment_char, stringsAsFactors=stringsAsFactors);
	file.remove(tmp_fn);
	return(tmp_df);
}
# system.time(rmsk <- read.delim2("rmsk.txt", sep="\t", quote="", header=T))
# system.time(rmsk_chr <- rmsk[rmsk$genoName == "chr19",])
# system.time(rmsk_chr2 <-retrieve_subset(fn="rmsk.txt", colname="genoName", selected_value="chr19"));
#system.time(rmsk_chr2 <-retrieve_subset("rmsk.txt", "genoName", "chr19"));


# fns <- list.files(".", pattern=".txt2")
# ids <- sub("isca", "", sub(".txt2", "", fns))
# colname <- "chrom"
# selected_value <- "chr19";
# df_list <- import_df_from_multiple_files(fns, ids, colname, selected_value);
# files_df: [filename]
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
	#list_outfiles <- list.files(".", pattern=".bed")
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
		#data <- data.frame(X=integer(col_n), Y=integer(col_n), Z=integer(col_n), stringsAsFactors=F);
		#for(j in 1 : col_n)
		#{
		#	data[j] <- c(rownames(plot_mtx[i,j]), colnames(plot_mtx[i,j]), plot_mtx[i,j]) 
		#}
		#plot_data <- rbind(plot_data,data);
		
		#plot_data <- rbind(plot_data, data.frame(X=colnames(plot_mtx[i,]),Y=rep(rownames(plot_mtx[i,]), col_n), Z=as.numeric(plot_mtx[i,]), stringsAsFactors=F))
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



