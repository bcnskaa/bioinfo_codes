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
import json
import process_features
import process_fasta



# feature table
global FEATURE_TABLE_NAME, FEATURE_TABLE_BREAK_SIZE, FEATURE_TABLE_OUT_FILENAME, FEATURE_TABLE_OUT_FORMAT, FEATURE_NAME, FEATURE_DESCRIPTION, FEATURE_INPUT_FILENAME
# general feature
global FEATURE_LIST, FEATURE_HEADER, FEATURE_SORT, FEATURE_FORMAT, FEATURE_SUBCLASS_LABEL, FEATURE_TYPE,FEATURE_SUBCLASS_VALUE, FEATURE_HEADER_LABEL, FEATURE_HEADER_INDEX, FEATURE_FILTER, FEATURE_FILTER_LABEL, FEATURE_FILTER_LABEL_INDEX, FEATURE_FILTER_VALUE, FEATURE_MAPPER
# ranged feature
global FEATURE_HEADER_START_LABEL, FEATURE_HEADER_START_INDEX, FEATURE_HEADER_END_LABEL, FEATURE_HEADER_END_INDEX
# mapper
global FEATURE_MAPPER_HEADER, FEATURE_MAPPER_START_LABEL, FEATURE_MAPPER_START_INDEX, FEATURE_MAPPER_END_LABEL, FEATURE_MAPPER_END_INDEX


# Feature table
FEATURE_TABLE_NAME="feature_table_name"
FEATURE_TABLE_BREAK_SIZE="feature_table_break_size"
FEATURE_TABLE_OUT_FILENAME="feature_table_ofn"
FEATURE_TABLE_OUT_FORMAT="feature_table_ofn_format"
FEATURE_LIST="features"

# File formats
FILE_FORMART_FASTA="fasta"
FILE_FORMAT_CSV="csv"
FILE_FORMAT_JSON="json"
FILE_FORMAT_TXT="txt"
FILE_FORMAT_XML="xml"

# general feature
FEATURE_NAME="feature_name"
FEATURE_DESCRIPTION="feature_desc"
FEATURE_INPUT_FILENAME="feature_ifn"
FEATURE_INPUT_FORMAT="feature_ifn_format"

FEATURE_SORT="feature_sort"
FEATURE_TYPE="feature_type"
FEATURE_TYPE_POINT="point"
FEATURE_TYPE_RANGE="range"
FEATURE_TYPE_INTENSITY="intensity"
FEATURE_TYPE_SEQUENCE="sequence"
FEATURE_TYPE_MAPPER="mapper"
FEATURE_TYPE_FILTER="filter"

FEATURE_SUBCLASS_VALUE="feature_subclass_value"
FEATURE_HEADER="feature_header"  # True or False
FEATURE_HEADER_LABEL="feature_header_label"
FEATURE_HEADER_INDEX="feature_header_idx"

FEATURE_FILTER="feature_filter"	 # True or False
FEATURE_FILTER_LABEL="feature_filter_label"
FEATURE_FILTER_LABEL_INDEX="feature_filter_label_idx"
FEATURE_FILTER_VALUE="feature_filter_value"
FEATURE_FILTER_VALUE_FIXED="feature_filter_value_fixed"
FEATURE_FILTER_DELIMITER=","

FEATURE_SUBCLASS_LABEL="feature_subclass_label"
FEATURE_SUBCLASS_INDEX="feature_subclass_idx"
FEATURE_SUBCLASS_VALUES="feature_subclass_idx"

FEATURE_MAPPER="feature_mapper"

# ranged feature
FEATURE_HEADER_START_LABEL="feature_header_start_label"
FEATURE_HEADER_START_INDEX="feature_header_start_idx"
FEATURE_HEADER_END_LABEL="feature_header_end_label"
FEATURE_HEADER_END_INDEX="feature_header_end_idx"

# mapper
FEATURE_MAPPER_HEADER="feature_header"
FEATURE_MAPPER_START_LABEL="feature_header_start_label"
FEATURE_MAPPER_START_INDEX="feature_header_start_idx"
FEATURE_MAPPER_END_LABEL="feature_header_end_label"
FEATURE_MAPPER_END_INDEX="feature_header_end_idx"



# Class for holding and executing feature descriptors
class feature:
	FEATURE_CTX = "__FEATURE_CTX__"
	FEATURE_DESC = "__FEATURE_DESC__"
	
	# Constructor
	def __init__(self, features_json_data):

		feature_n = len(features_json_data)
		
		self.features = {}
		self.mappers = {}
	
		for i in range(0, feature_n):
			_feature = features_json_data[i]
			
			# Parsing a point or range feature
			if _feature[FEATURE_TYPE] == FEATURE_TYPE_POINT or _feature[FEATURE_TYPE] == FEATURE_TYPE_RANGE:
				# check if it has a legitimate definition 
				if self.validate_point_range_feature(_feature):
					self.features[_feature[FEATURE_NAME]] = _feature
				
					print("Feature type=", _feature[FEATURE_TYPE], sep="")	
					
					# Insert on 04/02/2014
					if(_feature[FEATURE_INPUT_FORMAT] == FILE_FORMAT_CSV):
						this_feature = process_features.read_features_from_csv_file(_feature[FEATURE_NAME], _feature[FEATURE_INPUT_FILENAME])
						spos_label = _feature[FEATURE_HEADER_START_LABEL]
						epos_label = _feature[FEATURE_HEADER_END_LABEL]
						
						print("Process range feature")
						#process_features.print_feature_header(this_feature)
						#process_features.print_feature(this_feature, 10)
						
						poss = process_features.process_range_features(this_feature, spos_label, epos_label)
						print("Number of poss=", len(poss), sep="")
						
						
						# Insert on 09/03/2014
						# Iterate the filters

						if self.get_bool_val(_feature[FEATURE_FILTER]):
							print("Filter: ", _feature[FEATURE_FILTER_LABEL], " Value: ", _feature[FEATURE_FILTER_VALUE], sep="")
						
							filter_n = self.get_filter_count(_feature[FEATURE_FILTER_LABEL], _feature[FEATURE_FILTER_VALUE])
							
							if filter_n == 1:
								subclass_features = process_features.get_subclass_feature_re(this_feature, _feature[FEATURE_FILTER_LABEL], _feature[FEATURE_FILTER_VALUE])
								print("Original feature size: ", this_feature[process_features.FEATURE_LEN_TAG], sep="")
								print("Subclass feature size (filter=", _feature[FEATURE_FILTER_VALUE], "@", _feature[FEATURE_FILTER_LABEL], "): ", subclass_features[process_features.FEATURE_LEN_TAG], sep="")
							elif filter_n > 1:
								filter_labels =  _feature[FEATURE_FILTER_LABEL].split(FEATURE_FILTER_DELIMITER)
								filter_values =  _feature[FEATURE_FILTER_VALUE].split(FEATURE_FILTER_DELIMITER)
								subclass_features = this_feature
								
								for filter_idx in range(0, filter_n):
									print("Iteration: ", filter_idx + 1, sep="")
									subclass_features = process_features.get_subclass_feature_re(subclass_features, filter_labels[filter_idx],filter_values[filter_idx])
								
									print("Original feature size: ", this_feature[process_features.FEATURE_LEN_TAG], sep="")
									print("Subclass feature size (filter=", filter_labels[filter_idx], "@", filter_values[filter_idx], "): ", subclass_features[process_features.FEATURE_LEN_TAG], sep="")
								
								#process_features.export_features(subclass_features, "test.out")
								
							else:
								self.report_error("Definition of filter " + _feature[FEATURE_NAME] + " is invalid, be skipped")
						else:
							print("No feature filter is identified.")
						
				else:
					print("Definition of feature ", _feature[FEATURE_NAME], " is invalid, be skipped", sep="")
					self.report_error("Definition of feature " + _feature[FEATURE_NAME] + " is invalid, be skipped")
			# Parsing a mapper feature			
			elif _feature[FEATURE_TYPE] == FEATURE_TYPE_MAPPER:
				# check if it has a legitimate definition 
				if self.validate_mapper_feature(_feature):
					self.features[_feature[FEATURE_NAME]] = _feature	
				else:
					print("Definition of mapper ", _feature[FEATURE_NAME], " is invalid, be skipped", sep="")	
					self.report_error("Definition of mapper " + _feature[FEATURE_NAME] + " is invalid, be skipped")
			elif _feature[FEATURE_TYPE] == FEATURE_TYPE_SEQUENCE:
				# check if the sequence feature has a legitimate definition 
				if self.validate_sequence_feature(_feature):
					self.features[_feature[FEATURE_NAME]] = {self.FEATURE_DESC:_feature, self.FEATURE_CTX:self.process_sequence_feature(_feature)}	
				else:
					print("Definition of sequence feature ", _feature[FEATURE_NAME], " is invalid, be skipped", sep="")
					self.report_error("Definition of sequence feature " + _feature[FEATURE_NAME] + " is invalid, be skipped")
			elif _feature[FEATURE_TYPE] == FEATURE_TYPE_INTENSITY:
				print("Validation of Intensity Feature: Not implemented.")
			
			elif _feature[FEATURE_TYPE] == FEATURE_TYPE_FILTER:
				print(feature_json_data[i][FEATURE_NAME], "is a filter")
			else:
				print("Unknown feature type: ", json_data[FEATURE_LIST][i][FEATURE_TYPE], sep="")
				self.report_error("Unknown feature type.")
		
		

	# Error mechanism
	def report_error(self, e_msg):
		raise Exception(e_msg)

	
	# Get the number of features processed by this parser
	def get_feature_n(self):
		return 0


	# Validate sequence feature
	def validate_sequence_feature(self, feature_json_data):
		status = True
		
		if not FEATURE_INPUT_FILENAME in feature_json_data:
			self.report_error("Feature input filename is not provided.")
			return False
		if not FEATURE_INPUT_FORMAT in feature_json_data:
			self.report_error("Feature input format is not specified.")
			return False
			
		return status	


	# Process sequence feature
	def process_sequence_feature(self, feature_json_data):
		fn = feature_json_data[FEATURE_INPUT_FILENAME]
		seq = process_fasta.read_fasta(fn)
		print("Seq: ", seq[process_fasta.FASTA_HEADER_TAG], ", len:", len(seq[process_fasta.FASTA_SEQ_TAG]), sep="")

		return seq

	
	# Validate mapper feature
	def validate_mapper_feature(self, feature_json_data):
		status = True
			
		if not FEATURE_HEADER in feature_json_data:
			return False
		else:
			s = feature_json_data[FEATURE_HEADER].upper()
			if s == "T" or s == "TRUE":
				if not (FEATURE_HEADER_START_LABEL in feature_json_data and FEATURE_HEADER_END_LABEL in feature_json_data):
					return False
			elif s == "F" or "FALSE":
				if not (FEATURE_HEADER_START_INDEX in feature_json_data and FEATURE_HEADER_END_INDEX in feature_json_data):
					return False
			else:
				return False
		return status
		
	
	# Validate point or range feature				
	def validate_point_range_feature(self, feature_json_data):
		status = True
		
		if not FEATURE_INPUT_FILENAME in feature_json_data:
			self.report_error("No " + FEATURE_INPUT_FILENAME + " defined.")
			return False
		
		if not FEATURE_INPUT_FORMAT in feature_json_data:
			self.report_error("No " + FEATURE_INPUT_FORMAT + " defined.")
			return False

		if not FEATURE_HEADER in feature_json_data:
			self.report_error("No " + FEATURE_HEADER + " defined.")
			return False
		else:
			# Testing point feature
			if feature_json_data[FEATURE_TYPE] == FEATURE_TYPE_POINT:
				s = feature_json_data[FEATURE_HEADER].upper()
				if s == "T" or s == "TRUE":
					if not FEATURE_HEADER_LABEL in feature_json_data:
						self.report_error("No " + FEATURE_HEADER_LABEL + " defined.")
						return False
				elif s == "F" or "FALSE":
					if not FEATURE_HEADER_INDEX in feature_json_data:
						self.report_error("No " + FEATURE_HEADER_INDEX + " defined.")
						return False
				else:
					self.report_error(FEATURE_HEADER + " is invalid.")
					return False
			
			# Testing range feature
			elif feature_json_data[FEATURE_TYPE] == FEATURE_TYPE_RANGE:
				s = feature_json_data[FEATURE_HEADER].upper()
				#if s == "T" or s == "TRUE":
				if self.get_bool_val(s):
					if not (FEATURE_HEADER_START_LABEL in feature_json_data and FEATURE_HEADER_END_LABEL in feature_json_data):
						self.report_error("No " + FEATURE_HEADER_START_LABEL + " defined.")
						return False
				#elif s == "F" or "FALSE":
				elif not self.get_bool_val(s):
					if not (FEATURE_HEADER_START_INDEX in feature_json_data and FEATURE_HEADER_END_INDEX in feature_json_data):
						self.report_error("No " + FEATURE_HEADER_START_INDEX + " defined.")
						return False
				else:
					self.report_error(FEATURE_HEADER + " is invalid.")
					return False	
		
		if FEATURE_FILTER in feature_json_data:
			s = feature_json_data[FEATURE_FILTER].upper()
			##if not (s == "T" or s == "TRUE" or s == "F" or s == "FALSE"):
			if not self.validate_bool_val(s):
				self.report_error(FEATURE_FILTER + " is invalid.")
				return False
			
			##if s == "T" or s == "TRUE":
			if self.get_bool_val(s):
				if not FEATURE_FILTER_VALUE in feature_json_data:
					return False
				
				##if feature_json_data[FEATURE_HEADER].upper() == "T" or feature_json_data[FEATURE_HEADER].upper() == "TRUE":
				if self.get_bool_val(feature_json_data[FEATURE_HEADER]):
					if not FEATURE_FILTER_LABEL in feature_json_data:
						return False
				else:
					if not FEATURE_FILTER_LABEL_INDEX in feature_json_data:
						return False
		
		return status	
		
			
	# Check if the val is a bool char or string or not none
	def validate_bool_val(self, val):
		if val is None:
			return False
			
		val = val.upper()
		
		return (val == "T" or val == "TRUE" or val == "F" or val == "FALSE")
	
	
	# Get the bool value
	def get_bool_val(self, val):
		if val is None:
			return False
	
		val = val.upper()
		if(val == "T" or val == "TRUE"):
			return True
		else:
			return False
		# return the number of filter in descriptor
		
		
	# Get the number of filter in the descriptor
	def get_filter_count(self, filter_labels, filter_vals):
		#labels = filter_labels.split(FEATURE_FILTER_DELIMITER)
		#values = filter_vals.split(FEATURE_FILTER_DELIMITER)
		
		if filter_labels.count(FEATURE_FILTER_DELIMITER) != filter_vals.count(FEATURE_FILTER_DELIMITER):
			print("Filter pair does not match. label=", filter_labels.count(FEATURE_FILTER_DELIMITER), ", value=", filter_vals.count(FEATURE_FILTER_DELIMITER), sep="")
			return -1
		else:
			return filter_labels.count(FEATURE_FILTER_DELIMITER) + 1





class feature_descriptor:
	# parse the feature descriptors
	def __init__(self, json_data):
		self.feature_table_name = json_data[FEATURE_TABLE_NAME]
		self.feature_break_size = json_data[FEATURE_TABLE_BREAK_SIZE]
		self.feature_table_ofn = json_data[FEATURE_TABLE_OUT_FILENAME]
		self.feature_table_format = json_data[FEATURE_TABLE_OUT_FORMAT]
		self.features = feature(json_data[FEATURE_LIST])
		self.feature_n = self.features.get_feature_n()
		
	def process(self):
		print("Processing descriptors...")
	



# Parser for working with feature descriptors
class feature_descriptor_parser:
	# constructor
	def __init__(self, desc_fn):
		print("Calling feature_descriptor_parser")
		
		self.descriptor_src = desc_fn
		self.feature_table_name = ""
		self.feature_break_size = ""
		self.feature_table_ofn = ""
		self.feature_table_format = ""
		self.feature_sort = False
		self.features = None
		self.parse(desc_fn)        


	# Validate descriptor header
	def validate_descriptor_header(self, desc_json_data):
		if not FEATURE_TABLE_NAME in desc_json_data:
			print("Feature table name is not set.")
			return False
		if not FEATURE_TABLE_OUT_FORMAT in desc_json_data:
			print("Feature table format is not specified.")
			return False
		if not FEATURE_LIST in desc_json_data:
			print("No feature is available. Abort now.")
			return False
		if not FEATURE_TABLE_OUT_FILENAME in desc_json_data:
			print("Feature table output filename is not specified.")
			return False
		if not FEATURE_TABLE_BREAK_SIZE in desc_json_data:
			print("Feature break size is not defined.")
			return False
		return True

 
	# read feature definition file
	def parse(self, desc_fn):
		print("Calling read_feature_descriptor(): ", desc_fn, sep="")

		with open(desc_fn) as json_data_src:
			json_data = json.load(json_data_src)

		if not self.validate_descriptor_header(json_data):
			print("Invalid descriptor header found.")
			return None

		self.feature_table_name= json_data[FEATURE_TABLE_NAME]
		self.feature_break_size=json_data[FEATURE_TABLE_BREAK_SIZE]
		self.feature_table_ofn=json_data[FEATURE_TABLE_OUT_FILENAME]
		self.feature_table_format=json_data[FEATURE_TABLE_OUT_FORMAT]



		if (json_data["feature_sort"].upper() == "TRUE" or json_data["feature_sort"].upper() == "T"):
			self.feature_sort=True
		else:
			self.feature_sort=False

		self.features=feature(json_data["features"])

		return(json_data)


	# export the processed file
	def dump_feature_descriptor(self, desc, desc_ofn):
		print("export_feature_descriptor to ", desc_ofn, sep="")



	# process the content of the descriptor
	def process_descriptor(self, descriptor):
		print("Processing descriptor")



				
	
		
