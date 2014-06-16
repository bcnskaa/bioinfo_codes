/*
 * Copyright (c) 2013, SK Woolf <bcnskaa@gmail.com>.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/*
 * index.h
 *
 *  Created on: Nov 5, 2009
 *      Author: bcnskaa
 */

#ifndef INDEX_H_
#define INDEX_H_

extern const int DEFAULT_INDEX;

/**
 * Index table
 */
typedef struct {
	//int *sindex;		// Query Index, must be unique
	int *map_index;		// Subject index
	size_t mapping_num;
} INDEX;

extern INDEX* create_index(size_t);

/**
 * input sindex, return qindex
 */
extern int get_index(INDEX*, int);
extern void destroy_index(INDEX*);


extern int calculate_hash_index_value(char*, int);


typedef struct {
	char **patterns;
	size_t pattern_n;
	size_t **pattern_vals;
	size_t *pattern_val_lens;
	//size_t *pattern_val_buffer_sizes;
} PATTERN_TABLE;


extern PATTERN_TABLE* create_pattern_table(int pattern_len);
extern void destroy_pattern_table(PATTERN_TABLE* table);
extern PATTERN_TABLE* generate_pattern_profile(char* seq, int seq_len, int pattern_len, int overlapping_len);
extern char* get_pattern(int hash_idx, int pattern_len);
extern int get_hash_idx(char* pattern, int pattern_len);


#endif /* INDEX_H_ */
