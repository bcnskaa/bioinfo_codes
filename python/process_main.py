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
import _thread
import time
import traceback
import sys, getopt
import process_features
import re

import chrom_def
import feature_def
import process_ALL


# Initialize all the program parameters
def init():
	# Define some global variables
	global FEATURE_HEADER_TAG, FEATURE_LEN_TAG, FEATURE_CTX_TAG
	FEATURE_HEADER_TAG = feature_def.FEATURE_HEADER_TAG
	FEATURE_LEN_TAG = feature_def.FEATURE_LEN_TAG
	FEATURE_CTX_TAG = feature_def.FEATURE_CTX_TAG
	





# Main routine to start with
def main(argv):
	data = read_feature_from_file_ver2("test", argv[0])	
	
	#print("Feature length: ", data[FEATURE_LEN_TAG], sep="");
#	for header in data[FEATURE_HEADER_TAG]:
#		print()
#	for i in range(0, data[FEATURE_LEN_TAG]):
#		for header in data[HEADwER_TAG]:
	csv = features2csv_ver2(data)
	
	# subset data
	class_label = "repName"
	class_label_val = "AluY"
	#data_AluY = [data[i] for i in range(0, len(data) if data[class_label][i][0:len(class_label_val)] == class_label_val)]
	
	

	#breaks = process_features.create_breaks(0, max(poss), 100000)
	#bins = process_features.cut_poss2bins_ver3(poss, breaks)
	#bins = process_features.cut_poss2bins_ver4(poss, 100000)
	#print(bins)
	
	
	subclass_data = get_subclass_feature_re(data, "repName", "AluJo$")
	print("Subclass data size: ", subclass_data[FEATURE_LEN_TAG], sep="")	
	poss = process_features.process_ranged_features_ver2(subclass_data, "genoStart", "genoEnd")
	print("Number of poss: ", len(poss), sep="")
	bins = process_features.cut_poss2bins_ver4(poss, 100000)
	process_features.print_bins(bins)




if __name__ == '__main__':
	init()
	main(sys.argv[1:])
