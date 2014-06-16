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
 * index.c
 *
 *  Created on: Nov 16, 2009
 *      Author: bcnskaa
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <pthread.h>
//#include <cstring.h>

#include "index.h"
#include "threadpool.h"



const int BUFFER_SIZE = 100;

// gcc -o index index.c -lm

const int DEFAULT_INDEX = -1;


// A, T, G, C, N
const int HASH_VALUE_RANK = 5;
const int DEFAULT_PRECALCULATED_BASE_VALUES_SIZE = 20;


const char HASH_VALUE_A_CHAR = 'A';
const char HASH_VALUE_T_CHAR = 'T';
const char HASH_VALUE_G_CHAR = 'G';
const char HASH_VALUE_C_CHAR = 'C';
const char HASH_VALUE_N_CHAR = 'N';

const int HASH_VALUE_A = 0;
const int HASH_VALUE_C = 1;
const int HASH_VALUE_G = 2;
const int HASH_VALUE_T = 3;
const int HASH_VALUE_N = 4;
const int HASH_VALUE_ERROR = -1;



const int THREAD_N = 4;
pthread_mutex_t lock;
threadpool_t *pool;


void test_thread()
{
	pthread_mutex_init(&lock, NULL);

	// Create thread pool

}



//int PRECALCULATED_BASE_VALUES_SIZE = 20;
//extern size_t *PRECALCULATED_BASE_VALUES[PRECALCULATED_BASE_VALUES_SIZE];





int test_main(argc, argv)
	int argc;
	char **argv;
{
	PATTERN_TABLE *table;
	int i;

	printf("index: %d.\n", calculate_hash_index_value(argv[1], strlen(argv[1])));
	printf("index: %d.\n", calculate_hash_index_value("ATA", 3));

	printf("index=%d, pattern=%s\n", atoi(argv[2]), get_pattern(atoi(argv[2]),3));


	table = generate_pattern_profile(argv[3], strlen(argv[3]),  3, 2);

	for(i = 0; i < table->pattern_n; i++)
	{
		printf("%d %s: %d\n", i, table->patterns[i], table->pattern_val_lens[i]);
	}

	destroy_pattern_table(table); table = 0;



}

/**
 * Given a hash index and length of a pattern, this function will return the corresponding pattern
 */
char* get_pattern(hash_idx, pattern_len)
	int hash_idx;
	int pattern_len;
{
	char *pattern;
	int i, j;
	int remainder, base, denom;
	char base_list[] = {HASH_VALUE_A_CHAR, HASH_VALUE_T_CHAR, HASH_VALUE_G_CHAR, HASH_VALUE_C_CHAR, HASH_VALUE_N_CHAR};

	//printf("Pattern len=%d\n", pattern_len);

	pattern = (char*)calloc(pattern_len + 1, sizeof(char));
	memset(pattern, 0, pattern_len + 1);

	remainder = hash_idx;
	for(i = 0, j = pattern_len - 1; i < pattern_len; i++, j--)
	{
		denom = pow(HASH_VALUE_RANK, j);
		base = (int)floor(remainder / denom);
		remainder = remainder % denom;

		pattern[i] = base_list[base];
	}

	return pattern;
}


int get_hash_idx(char* pattern, int pattern_len)
{
	return calculate_hash_index_value(pattern, pattern_len);
}



size_t *PRECALCULATED_BASE_VALUES;
int PRECALCULATED_BASE_VALUES_SIZE = 0;
void initialize_precalculated_values()
{
	int i;
	int cardinality;

	cardinality = HASH_VALUE_RANK;

	PRECALCULATED_BASE_VALUES_SIZE = DEFAULT_PRECALCULATED_BASE_VALUES_SIZE;
	PRECALCULATED_BASE_VALUES = (size_t*)calloc(PRECALCULATED_BASE_VALUES_SIZE, sizeof(size_t));

	PRECALCULATED_BASE_VALUES[0] = 1;
	for(i = 1; i < PRECALCULATED_BASE_VALUES_SIZE; i++)
		PRECALCULATED_BASE_VALUES[i] = (size_t)cardinality * PRECALCULATED_BASE_VALUES[i - 1];


}

void deinitialize_precalculated_values()
{
	PRECALCULATED_BASE_VALUES_SIZE = 0;
	free(PRECALCULATED_BASE_VALUES); PRECALCULATED_BASE_VALUES = 0;
}

/**
 * Create a empty pattern table. The size of a table is calculated by the equation: HASH_VALUE_RANK ^ pattern_len.
 *
 */
PATTERN_TABLE* create_pattern_table(pattern_len)
	int pattern_len;
{
	PATTERN_TABLE* table;
	int i;

	if(PRECALCULATED_BASE_VALUES_SIZE == 0)
		initialize_precalculated_values();


	table = (PATTERN_TABLE*)calloc(1, sizeof(PATTERN_TABLE));

	// Calculate the permutation
	table->pattern_n = pow(HASH_VALUE_RANK, pattern_len);
	table->patterns = (char**)calloc(table->pattern_n, sizeof(char*));
	table->pattern_vals = (size_t**)calloc(table->pattern_n, sizeof(size_t*));
	table->pattern_val_lens = (size_t*)calloc(table->pattern_n, sizeof(size_t));
	//table->pattern_val_buffer_sizes = (size_t*)calloc(table->pattern_n, sizeof(size_t));

	for(i = 0; i < table->pattern_n; i++)
	{
		table->patterns[i] = get_pattern(i, pattern_len);
		//table->pattern_vals[i] = (size_t*)calloc(BUFFER_SIZE, sizeof(size_t));
		//table->pattern_vals[i] = 0;
		table->pattern_val_lens[i] = 0;

		//table->pattern_val_buffer_sizes[i] = BUFFER_SIZE;

		table->pattern_vals[i] = 0;
		//table->pattern_val_buffer_sizes[i] = 0;
	}


	return table;
}



void destroy_pattern_table(table)
	PATTERN_TABLE* table;
{
	int i;

	for(i = 0; i < table->pattern_n; i++)
	{
		free(table->pattern_vals[i]);
		table->pattern_vals[i] = 0;
		free(table->patterns[i]);
		table->patterns[i] = 0;
	}
	free(table->patterns);
	table->patterns = 0;
	free(table->pattern_vals);
	table->pattern_vals = 0;
	free(table->pattern_val_lens);
	table->pattern_val_lens = 0;
	//free(table->pattern_val_buffer_sizes);
	//table->pattern_val_buffer_sizes = 0;
	free(table);
	table = 0;


	if(PRECALCULATED_BASE_VALUES_SIZE != 0)
		deinitialize_precalculated_values();
}





PATTERN_TABLE* generate_pattern_profile(seq, seq_len, pattern_len, overlapping_len)
	char* seq;
	int seq_len;
	int pattern_len;
	int overlapping_len;
{
	PATTERN_TABLE* table;
	char pattern[pattern_len];
	int i, max_len, step, hash_idx;
	size_t *seq_profile;
	int *checkpoints;
	int cur_checkpoint;

	//int seq_len;

	table = create_pattern_table(pattern_len);

	// A profile of hash index of each position for linear time access
	seq_profile = (size_t*)calloc(seq_len, sizeof(size_t));
	memset(seq_profile, HASH_VALUE_ERROR, seq_len);


	//seq_len = strlen(seq);
	max_len = seq_len - pattern_len;
	step = (pattern_len - overlapping_len);
	for(i = 0; i <= max_len; i+= step)
	{
		strncpy(pattern, &seq[i], pattern_len);
		hash_idx = calculate_hash_index_value(pattern, pattern_len);


		seq_profile[i] = hash_idx;

		if(hash_idx == HASH_VALUE_ERROR)
		{
			printf("Unknown pattern: %s\n", pattern);
		} else {
			//printf("%d pattern=%s hash_idx=%d\n", i, pattern, hash_idx);


			table->pattern_val_lens[hash_idx] += 1;

//			if(table->pattern_val_lens[hash_idx] >= table->pattern_val_buffer_sizes[hash_idx])
//			{
//				table->pattern_val_buffer_sizes[hash_idx] += BUFFER_SIZE;
//				table->pattern_vals[hash_idx] = (size_t*)realloc(table->pattern_vals[hash_idx], table->pattern_val_buffer_sizes[hash_idx]);
//				printf("table->pattern_val_buffer_sizes[hash_idx]=%d, table->pattern_val_lens[hash_idx]=%d\n", table->pattern_val_buffer_sizes[hash_idx], table->pattern_val_lens[hash_idx]);
//
//			}

			//table->pattern_vals[hash_idx] = (size_t*)realloc(table->pattern_vals[hash_idx], table->pattern_val_lens[hash_idx]);
//			printf("d=%d\n", table->pattern_val_lens[hash_idx]);
//			table->pattern_vals[hash_idx][table->pattern_val_lens[hash_idx] - 1] = i;
		}
	}


	// checkpoints to store the current pattern_vals position
	checkpoints = (int*)calloc(table->pattern_n, sizeof(int));
	memset(checkpoints, 0, table->pattern_n);

	// Go through the pattern table and allocate the memory for the patten_len_vals
	for(i = 0; i < table->pattern_n; i+= step)
	{
		if(table->pattern_val_lens[i] > 0)
		{
			table->pattern_vals[i] = (size_t*)calloc(table->pattern_val_lens[i], sizeof(size_t));
			memset(table->pattern_vals[i], HASH_VALUE_ERROR, table->pattern_val_lens[i]);
 		}
	}


	// Assign the hash index into pattern_vals from seq_profile
	for(i = 0; i < seq_len; i++)
	{
		hash_idx = seq_profile[i];
		if(hash_idx != HASH_VALUE_ERROR)
		{
			cur_checkpoint = checkpoints[hash_idx];
			table->pattern_vals[hash_idx][cur_checkpoint] = i;
			checkpoints[hash_idx] = cur_checkpoint + 1;
		}
	}


	// Clean up
	free(seq_profile); seq_profile = 0;
	free(checkpoints); checkpoints = 0;


	return table;
}


//
//size_t PRECALCULATED_BASE_VALUES(k)
//	int k;
//{
//	size_t PRECALCULATED_BASE_VALUES[] = {1, HASH_VALUE_RANK, HASH_VALUE_RANK * HASH_VALUE_RANK, HASH_VALUE_RANK * HASH_VALUE_RANK * HASH_VALUE_RANK, HASH_VALUE_RANK * HASH_VALUE_RANK * HASH_VALUE_RANK};
//
//	return PRECALCULATED_BASE_VALUES[k];
//}

//
int calculate_hash_index_value(word, word_len)
	char *word;
	int word_len;
{
	int len;
	int i, j, k;
	int idx_value;
	char c;
	//int cardinality;
	int base_val;

	//len = strlen(word);
	len = word_len;

	//cardinality = HASH_VALUE_RANK;
	idx_value = 0;
	for(i = 0, j = len - 1; i < len; i++, j--)
	{
		c = word[i];

		base_val = PRECALCULATED_BASE_VALUES[j];
		/*
		base_val = 1;
		for(k = 0; k < j; k++)
			base_val *= cardinality;
		*/

		if(c == HASH_VALUE_A_CHAR)
			idx_value += (HASH_VALUE_A * base_val);
		else if(c == HASH_VALUE_T_CHAR)
			idx_value += (HASH_VALUE_T * base_val);
		else if(c == HASH_VALUE_G_CHAR)
			idx_value += (HASH_VALUE_G * base_val);
		else if(c == HASH_VALUE_C_CHAR)
			idx_value += (HASH_VALUE_C * base_val);
		else if(c == HASH_VALUE_N_CHAR)
			idx_value += (HASH_VALUE_N * base_val);
		else
			return HASH_VALUE_ERROR;

// switch case can't be used when
//		switch(c)
//		{
//		case 'A':
//			idx_value += (HASH_VALUE_A * (int)pow((double)cardinality, j));
//			break;
//		case 'T':
//			idx_value += (HASH_VALUE_T * (int)pow((double)cardinality, j));
//			break;
//		case 'G':
//			idx_value += (HASH_VALUE_G * (int)pow((double)cardinality, j));
//			break;
//		case 'C':
//			idx_value += (HASH_VALUE_C * (int)pow((double)cardinality, j));
//			break;
//		case 'N':
//			idx_value += (HASH_VALUE_N * (int)pow((double)cardinality, j));
//			break;
//		default:
//			return HASH_VALUE_ERROR;
//		}
	}

	return idx_value;
}



int calculate_hash_index_value_slow(word, word_len)
	char *word;
	int word_len;
{
	int len;
	int i, j;
	int idx_value;
	char c;
	int cardinality;

	//len = strlen(word);
	len = word_len;

	cardinality = HASH_VALUE_RANK;
	idx_value = 0;
	for(i = 0, j = len - 1; i < len; i++, j--)
	{
		c = word[i];

		if(c == HASH_VALUE_A_CHAR)
			idx_value += (HASH_VALUE_A * (int)pow((double)cardinality, j));
		else if(c == HASH_VALUE_T_CHAR)
			idx_value += (HASH_VALUE_T * (int)pow((double)cardinality, j));
		else if(c == HASH_VALUE_G_CHAR)
			idx_value += (HASH_VALUE_G * (int)pow((double)cardinality, j));
		else if(c == HASH_VALUE_C_CHAR)
			idx_value += (HASH_VALUE_C * (int)pow((double)cardinality, j));
		else if(c == HASH_VALUE_N_CHAR)
			idx_value += (HASH_VALUE_N * (int)pow((double)cardinality, j));
		else
			return HASH_VALUE_ERROR;

// switch case can't be used when
//		switch(c)
//		{
//		case 'A':
//			idx_value += (HASH_VALUE_A * (int)pow((double)cardinality, j));
//			break;
//		case 'T':
//			idx_value += (HASH_VALUE_T * (int)pow((double)cardinality, j));
//			break;
//		case 'G':
//			idx_value += (HASH_VALUE_G * (int)pow((double)cardinality, j));
//			break;
//		case 'C':
//			idx_value += (HASH_VALUE_C * (int)pow((double)cardinality, j));
//			break;
//		case 'N':
//			idx_value += (HASH_VALUE_N * (int)pow((double)cardinality, j));
//			break;
//		default:
//			return HASH_VALUE_ERROR;
//		}
	}

	return idx_value;
}


INDEX* create_index(size)
	size_t size;
{
	INDEX *idx;
	int i;

	idx = (INDEX*)calloc(1, sizeof(INDEX));
	idx->map_index = (int*)calloc(size, sizeof(int));
	//memset(idx->map_index, DEFAULT_INDEX, size);

	for(i = 0; i < size; i++) {
		idx->map_index[i] = DEFAULT_INDEX;
	}

	//idx->qindex = (int*)calloc(size, sizeof(int));
	idx->mapping_num = 0;

	return idx;
}

/**
 * input sindex, return qindex
 */
int get_index(idx, index)
	INDEX* idx;
	int index;
{
	int map_index;

	map_index = DEFAULT_INDEX;
	if(idx->mapping_num > index) {
		map_index = idx->map_index[index];
	}

	return map_index;
}

void destroy_index(idx)
	INDEX *idx;
{
	free(idx->map_index); idx->map_index = 0;
	//free(idx->sindex); idx->sindex = 0;
	free(idx); idx = 0;
}

