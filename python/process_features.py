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
import re
import math
import _thread
import json
import feature_descriptor_parser
import operator



global FEATURE_NAME_TAG, FEATURE_HEADER_TAG, FEATURE_CTX_TAG, FEATURE_LEN_TAG
FEATURE_NAME_TAG = "__FEATURE_NAME__"
FEATURE_HEADER_TAG = "__FEATURE_HEADER__"
FEATURE_CTX_TAG = "__FEATURE_CTX__"
FEATURE_LEN_TAG =  "__FEATURE_CTX_LEN__"




# read_feature_from_file: csv
def read_features_from_csv_file(feature_name, feature_ifn):
	features = {}
    	
	print("Reading features from ", feature_ifn, "...", sep="")
	features_data = open(feature_ifn).read()
	features_data = features_data.splitlines()
    
	features[FEATURE_NAME_TAG] = feature_name
	features[FEATURE_HEADER_TAG] = features_data[0].split()
	features[FEATURE_CTX_TAG] = []
	for i in range(1, len(features_data)):
		features[FEATURE_CTX_TAG].append(features_data[i].split("\t"))

	features[FEATURE_LEN_TAG] = len(features[FEATURE_CTX_TAG])
	return features





# read_feature_from_file: json
def read_features_from_json_file(feature_name, feature_ifn):
	features = {}
	print("reading features from ", feature_ifn, "...")
    
	features_data = open(feature_ifn).read()
	features_data = features.splitlines()
    
	with open(desc_fn) as json_data_src:
		json_data = json.load(json_data_src)    
        
	features[FEATURE_NAME_TAG] = feature_name
	features[FEATURE_HEADER_TAG] = json_data[FEATURE_HEADER_TAG]
	features[FEATURE_CTX_TAG] = json_data[FEATURE_CTX_TAG]
	features[FEATURE_LEN_TAG] = json_data[FEATURE_LEN_TAG]
	
	return features



#
def sort_features(features, class_label):
	class_label_idx = features[FEATURE_HEADER_TAG].index(class_label)
	
	sorted_features = features.sort(key=operator.itemgetter(class_label_idx))

	#feature_poss = [features[FEATURE_CTX_TAG][i][class_label_idx] for i in range(0, features[FEATURE_LEN_TAG])]
	
	return sorted_features


# 
def get_unique_labels(features, class_label):
	class_label_idx = features[FEATURE_HEADER_TAG].index(class_label)
	
	# initialize an empty list of labels
	unique_labels = []
	
	for i in range(0, features[FEATURE_LEN_TAG]):
		if not features[FEATURE_CTX_TAG][i][class_label_idx] in unique_labels:
			unique_labels.append(features[FEATURE_CTX_TAG][i][class_label_idx])
	
	unique_labels.sort()
	return unique_labels
	

	

def get_unique_label_count(features, class_label):
	class_label_idx = features[FEATURE_HEADER_TAG].index(class_label)	
	
	class_label_values = [features[FEATURE_CTX_TAG][i][class_label_idx] for i in range(0, features[FEATURE_LEN_TAG])]
	class_unique_labels = get_unique_labels(features, class_label)
	
	unique_label_count = {x:class_label_values.count(x) for i,x in enumerate(class_unique_labels)}

	return unique_label_count
	

# feature2csv
def features2csv():
	print("features2csv()")


# create_breaks
def create_breaks(spos, epos, break_size):
	print("create_breaks...")
	breaks = [[spos + (i * break_size), spos + ((i + 1) * break_size)]for i in range(0, (math.ceil((epos - spos) / break_size)))]
	return breaks



# cut_poss2bins 
def cut_poss2bins():
	print("cut_poss2bins()")



# process_feature
def process_range_features(features, spos_label, epos_label):
	
	features_header = features[FEATURE_HEADER_TAG]
	
	spos_label_idx = features_header.index(spos_label)
	epos_label_idx = features_header.index(epos_label)
	
	print("process_features(): spos_label_idx=", spos_label, "@col_idx=", spos_label_idx, " , epos_label_idx=", epos_label, "@col_idx=", epos_label_idx, sep="")
	
	poss = []
	for i in range(0, features[FEATURE_LEN_TAG]):
		#print("feature: ", features[FEATURE_CTX_TAG][i][spos_label_idx])
		poss.extend([i for i in range(int(features[FEATURE_CTX_TAG][i][spos_label_idx]), int(features[FEATURE_CTX_TAG][i][epos_label_idx]))])
	
	#print("Spos_label: ", spos_label, " Epos_label: ", epos_label, sep="")
	print("process_features() processd size: ", features[FEATURE_LEN_TAG], sep="")
	
	# sort the item		
	return poss


#  
def get_subclass_feature_re(features, class_label, subclass_value):
	print("get_subclass_feature_re(), class_label=", class_label, ", subclass_value=", subclass_value, sep="")
	class_label_idx = features[FEATURE_HEADER_TAG].index(class_label)
	class_discriminator = re.compile(subclass_value)
	
	subclass_features = {}
	subclass_features[FEATURE_HEADER_TAG] = features[FEATURE_HEADER_TAG]
	subclass_features[FEATURE_NAME_TAG] = features[FEATURE_NAME_TAG] + ":" + subclass_value
	subclass_features[FEATURE_CTX_TAG] = [features[FEATURE_CTX_TAG][i] for i in range(0, features[FEATURE_LEN_TAG]) if class_discriminator.match(features[FEATURE_CTX_TAG][i][class_label_idx])]	
	subclass_features[FEATURE_LEN_TAG] = len(subclass_features[FEATURE_CTX_TAG])
	
	return subclass_features



# poss2bin: assuming non-overlapping break region
def cut_poss2bins(poss, breaks):
	print("cut_poss2bin()")
	
	# Sort the poss
	poss.sort()		
	breaks_n = len(breaks)
	poss_n = len(poss)
	bins = [0 for i in range(0, breaks_n)]
	
	cur_break_idx = 0
	spos,epos = breaks[cur_break_idx]
	for i in range(0, poss_n):
		if poss[i] < epos:
			if poss[i] >= spos:
				bins[cur_break_idx] += 1
		else:
			if cur_break_idx < breaks_n:
				cur_break_idx += 1
				spos,epos = breaks[cur_break_idx]
			else:
				break
	return bins




# Print
def print_bins():
	print("print_bins()")



# 
def export_features(features, out_fn):
	print("export_features(): exporting features ", features[FEATURE_NAME_TAG], "(size=", features[FEATURE_LEN_TAG], ") to ", out_fn, "...", sep="")
	
	headers = "\t".join(features[FEATURE_HEADER_TAG])
	header_size = len(features[FEATURE_HEADER_TAG])
	
	#print("export_features(): header size=", header_size, sep="")
	
	with open(out_fn, "w") as f_out:
		f_out.write(headers)
		for i in range(0, features[FEATURE_LEN_TAG]):
			cur_line = ""
			for j in range(0, header_size):
				if j == 0:
					#f_out.write(features[FEATURE_CTX_TAG][i][j])
					cur_line = cur_line + features[FEATURE_CTX_TAG][i][j]
				else:
					#f_out.write("\t", features[FEATURE_CTX_TAG][i][j], sep="")
					cur_line = cur_line + "\t" + features[FEATURE_CTX_TAG][i][j]
				#print(cur_line)	
				
			cur_line = cur_line + "\n"		
			f_out.write(cur_line)
		
		


#  
def export_feature_table(features):
	print("export_feature_table")



# 
def print_feature_header(features):
	headers = features[FEATURE_HEADER_TAG]
	
	header_size = len(headers)
	print("Feature: ", features[FEATURE_NAME_TAG], sep="")
	for i in range(0, header_size):
		print(i, ":", headers[i], sep="")
		
		
# 
def print_feature(features, idx):
	headers = features[FEATURE_HEADER_TAG]
	
	header_size = len(headers)
	print("Feature: ", features[FEATURE_NAME_TAG], " at pos=", idx, sep="")
	
	for i in range(0, header_size):
		print(headers[i], ":", features[FEATURE_CTX_TAG][idx][i], sep="")
