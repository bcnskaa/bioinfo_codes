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
import sys
import _thread
import time
import traceback
import math
import glob
# from os.path import isfile, join
from os import listdir
import os
import itertools
from numpy import *
import threading
from multiprocessing import Process, Queue


# #
# Todo: multithreading the worker code
#
# FASTA sequence object
class fasta_seq:
	header = ""
	seq = ""
	UNKNOWN_NUCLEOTIDE = "X"
	
	def __init__(self):
		self.header = ""
		self.seq = ""
	
	# get the length of sequence	
	def length(self):
		return len(self.seq)
	
	def name(self):
		return self.header
	
	def nucleotideAt(self, pos):
		if pos >= 0 and pos < self.length():
			return self.seq[pos]
		else:
			return self.INKNOWN_NUCLEOTIDE
	
	def subseqAt(self, spos, epos):
		return self.seq[spos:epos]
		 
	def index_of(self, pattern):
		pattern_indice = []
		pattern_len = len(pattern)
		search_bound_max = self.length() - pattern_len

		print("Searching pattern ", pattern, "...", sep="")
		
		for cur_pos in range(1, search_bound_max):
			subseq = self.subseqAt(cur_pos, cur_pos + pattern_len)
			if subseq == pattern:
				pattern_indice.append(cur_pos)
		
		return pattern_indice
		
	def nonoverlapping_index_of(self, pattern):
		pattern_indice = []
		pattern_len = len(pattern)
		
		seq_len = self.length() 
		
		search_bound_max = seq_len - pattern_len
		
		cur_pos = 0
		
		while(True):
			if cur_pos >= seq_len:
				break
			else:
				subseq = self.subseqAt(cur_pos, cur_pos + pattern_len)
				if subseq == pattern:
					pattern_indice.append(cur_pos)
					cur_pos += pattern_len
				else:
					cur_pos += 1
			
		return pattern_indice
		
	def read_fasta_seq_from_file(self, filename, fasta_seq_id):
	# Read from fasta sequences and select the sequence object specified by the obj_id				
	
		searched_object_n = 0
		cflag_found = False
		list_lines = []
		
		print("Searching ", filename, " for ", fasta_seq_id, "...", sep="")
		
		with open(filename, "r") as f_in:
			for line in f_in:
				# Remove space/newline character(s) from leading and tailing positions
				line.strip()
				
				if line.startswith(">"):
					id = line.replace(">", "")
					
					searched_object_n += 1
					
					# End of the selected object
					if cflag_found == True:
						break
						
					# Found the object, start to store the sequence
					if id == fasta_seq_id:
						print(fasta_seq_id, " found", sep="")
						cflag_found = True
						self.header = id
				else:
					# self.seq += line
					list_lines.append(line.rstrip('\n'))
					
		if len(list_lines) > 0:
			self.seq = "".join(list_lines)
			self.seq = self.seq.upper()




class indexed_fasta_seq:	
	fasta_seq = 0
	index_table = []
	index_word_size = 8
	
	patterns = []
	
	# Number of process to be created
	process_n = 8
	QueuePool = []
	ProcessPool = []
	ProcessIDs = []
	
	
	nucleotide_list = ["A", "C", "G", "N", "T"]
	base_rank = 0
	A_val = 0
	C_val = 1
	G_val = 2
	N_val = 3
	T_val = 4
	U_val = -1
	
	
	
	def __init__(self, filename, fasta_seq_id):
		self.index_table = []
		self.base_rank = len(self.nucleotide_list)
		self.fasta_seq = fasta_seq()
		self.fasta_seq.read_fasta_seq_from_file(filename, fasta_seq_id)
		# self.build_index()  
		self.build_index_by_processes()
	
	
	def read_fasta_seq_from_file(self, filename, fasta_seq_id): 
		self.fasta_seq.read_fasta_seq_from_file(filename, fasta_seq_id)
		# self.build_index()
		 
			  
	def show_overrepresent_pattern(self, top_n=10):
		pattern_occurrences = [len(v) for v in self.index_table]
		# lambda is a keyword for creating an inline function
		# enumerate(list) iterates elements of the given list
		# overrepresent_pattern_idx = max(enumerate(pattern_occurrences), key=lambda x: x[1])[0]
		# print("Over-represented pattern: ", "".join(self.patterns[overrepresent_pattern_idx]), " occurence=", pattern_occurrences[overrepresent_pattern_idx], sep="")
		
		sorted_patterns_with_high_occurrence = sorted(enumerate(pattern_occurrences), key=lambda x:x[1], reverse=True)
		
		print("Top ", top_n, " over-represented patterns with high hits:", sep="")
		[print("".join(self.patterns[sorted_patterns_with_high_occurrence[i][0]]), " occurrence=", pattern_occurrences[sorted_patterns_with_high_occurrence[i][0]], sep="") for i in range(top_n)]
	
		

	
	def build_index_by_processes(self):
		if(self.fasta_seq != 0):
			# Distribute the workload
			seq_len = self.fasta_seq.length()
					
			print("Generating pattern table of size ", self.index_word_size, "...", end="", sep="")
			
			self.patterns = self.generate_pattern_table(self.index_word_size)
			self.index_table = [[] for i in range(len(self.patterns))]
			
			print(len(self.patterns), " patterns generated.", sep="")
			
			
			numberOfBaseToBeCompletedByEachProcessed = round(seq_len / self.process_n)
			
			
			sposs = [i for i in range(0, seq_len, numberOfBaseToBeCompletedByEachProcessed)] 
			eposs = sposs[1:(len(sposs) - 1)]
			eposs.append(seq_len)
			
			
			[self.ProcessIDs.append("Process_"+str(i)) for i in range(self.process_n)]
			
			
			for i in range(self.process_n):
				pid = self.ProcessIDs[i]
				print("pid=", pid, " spos=", sposs[i], " epos=", eposs[i], sep="")
				# Queue pool
				q = Queue()
				p = Process(target=self.call_process_to_build_index, args=(sposs[i], eposs[i], pid, q,))
				self.QueuePool.append(q)
				self.ProcessPool.append(p)
				
			
			# Process Pool
			#self.ProcessPool = [Process(target=self.call_process_to_build_index, args=(self.ProcessIDs[i], sposs[i], eposs[i], self.QueuePool[i],)) for i in range(self.process_n)]
			
			[p.start() for p in self.ProcessPool]
			[p.join() for p in self.ProcessPool]
			
			

			
			print("Integerating results for ", len(self.patterns), " items...", sep="")
			
			# Combined index
			queue_res = [q.get() for q in self.QueuePool]

			half_progress = round(len(self.patterns) / 2)
			quater_progress = round(len(self.patterns) / 4)
			last_progress = round((4 * len(self.patterns)) / 3)
			
			print("0%...")
			for cur_idx in range(len(self.patterns)):
				print("idx=", cur_idx, sep="")
				[self.index_table[cur_idx].extend(queue_res[i][cur_idx]) for i in range(self.process_n)]
				#for i in range(self.process_n):
					
				if cur_idx == quater_progress:
					print("25%...")
				if cur_idx == half_progress:
					print("50%...")
				if cur_idx == last_progress:
					print("75%...")
			print("100% completed.")	
				
			print("Building of index table complete.")	
			
			self.ProcessPool = []
			self.ProcessIDs = []
			
		else:
			print("No sequence available, indexing aborted.")		




	def call_process_to_build_index(self, spos, epos, pid, queue):
		print(pid, "is about to start building index...")
		
		# Create an empty table for storing the position matching patterns
		my_index_table = [[] for i in range(len(self.patterns))]
		#print(len(self.patterns), " patterns generated.", sep="")
			
		print("Scanning between ", spos, " and ", epos, " and start to build index table.", sep="")
		
		count_n = 0
		for cur_pos in range(spos, epos):
			pattern = self.fasta_seq.subseqAt(cur_pos, cur_pos + self.index_word_size)
			idx = self.calculate_pattern_permutation_index(pattern)
			if(idx != self.U_val):
				my_index_table[idx].append(cur_pos)
			count_n += 1
		
		print(pid, ": number of data processed=", count_n, sep="")
		
		return my_index_table
	
	
	
	def build_index(self):
		if(self.fasta_seq != 0):
			seq_len = self.fasta_seq.length()
			max_bound = seq_len - self.index_word_size
			
			print("Generating pattern table of size ", self.index_word_size, "...", end="", sep="")
			self.patterns = self.generate_pattern_table(self.index_word_size)
			self.index_table = [[] for i in range(len(self.patterns))]
			#print(len(self.patterns), "patterns generated.", sep="")
			
			step_len = math.ceil(seq_len / 100)
			milestone_len = math.ceil(seq_len / 10)
			cur_progress = 0
			
			print("Scanning: ", max_bound, " windows", sep="")
			print("Building index table:")
			
			print("[", end="")
			
			for cur_pos in range(max_bound):
				pattern = self.fasta_seq.subseqAt(cur_pos, cur_pos + self.index_word_size)
				idx = self.calculate_pattern_permutation_index(pattern)
				if(idx != self.U_val):
					self.index_table[idx].append(cur_pos)
				
				# print(cur_pos, ".", sep="")
				
				if(cur_pos != 0 and cur_pos % step_len == 0):
					print(".", end="", sep="")
				if(cur_pos != 0 and cur_pos % milestone_len == 0):
					cur_progress += 10
					print("", cur_progress, "%", end="", sep="")
			print("]")
		else:
			print("No sequence available, indexing aborted.")
			
	
	def generate_pattern_table(self, pattern_len):
		return [p for p in itertools.product(self.nucleotide_list, repeat=pattern_len)]


	def calculate_pattern_permutation_index(self, pattern):
		base_values = [self.get_base_value(base) for base in pattern]
		
		# print("Base values: ", base_values, sep="")
		
		if(self.U_val in base_values):
			# print("Pattern contains unrecognized base character: ", pattern, sep="")
			return -1
		
		pos_values = list(range(0, len(pattern)))
		pos_values.reverse()
		
		# zip() function returns a list of tuples, where i-tuple contains i-elements from each of list
		idx = sum([b * (self.base_rank ** p) for b, p in zip(base_values, pos_values)])
		
		# print("Idx:", idx, sep="")
		
		return idx


	def get_base_value(self, base):
		# return nucleotide_list.index(base)
		if(base == "A"):
			return self.A_val
		elif(base == "C"):
			return self.C_val
		elif(base == "G"):
			return self.G_val
		elif(base == "T"):
			return self.T_val
		elif(base == "N"):
			return self.N_val
		else:
			return self.U_val


def test_code():
	sample_size = 100
	rand_size = [random.choice(range(100)) for i in range(1, sample_size)]
	a = [[random.choice(range(1000)) for j in range(n)] for n in rand_size]
	sorted_a = sorted(enumerate(a), key=lambda x:len(x[1]), reverse=True)
	
	

def main(argv):
	  
	if(len(argv) == 4):
		fasta_fn = argv[0]
		seq_obj_id = argv[1]
		seq_pattern = argv[2]
		idx_out_fn = argv[3]
		
		seq = fasta_seq()
		seq.read_fasta_seq_from_file(fasta_fn, seq_obj_id)
		print("Length of sequence: ", seq.length(), sep="")
		# print("Subseq: ", seq.subseqAt(0, 100))
		
		indexed_seq = indexed_fasta_seq(fasta_fn, seq_obj_id)
		
		#print("Calling threads...")
		# indexed_seq.startIndexing()
		
		idx = indexed_seq.calculate_pattern_permutation_index(seq_pattern)
		print("Index: ", idx, " pattern, ", seq_pattern, sep="")
		
		# print("Building index table...")
		# indexed_seq.build_index()
		indexed_seq.show_overrepresent_pattern()
		
		patterns = indexed_seq.generate_pattern_table(len(seq_pattern))
		patterns.sort()
		# print(patterns)
		print("Pattern at position ", idx, "=", patterns[idx], sep="")
		
		# patterns = seq.index_of(seq_pattern)
		# print("Number of pattern found:", len(patterns), sep="")
	else:
		print("Usage: FASTA-INFILE SEQ-OBJ-ID SEQ-PATTERN INDEX-OUTFILE")



# main 
if __name__ == "__main__":
	main(sys.argv[1:])

