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
 * fasta_reader.h
 *
 *  Created on: Oct 30, 2009
 *      Author: bcnskaa
 */

#ifndef FASTA_READER_H_
#define FASTA_READER_H_

extern const int HEADER_NOT_FOUND;

typedef struct {
	char *header;
	char *seq;
	size_t seq_len;
} SeqObj;

extern size_t CHAR_PER_LINE;

// Input: filename, header
extern SeqObj* read_fasta_by_id(char*, char*);
extern void destroy_SeqObj(SeqObj*);

// Get the allele at the position specified
extern char charAt_SeqObj(SeqObj*, int);
//extern char get_allele_SeqObj(SeqObj*, int);
extern char* get_subseq_SeqObj(SeqObj*, int, int);


extern void print_fasta(char*, SeqObj*);
extern void display_fasta(SeqObj*, int, int);

// Convert the sequence into upper case
extern void seq_toupper(SeqObj*);

//extern int check_if_header_exist(char*, char*);
extern long int get_header_pos(char *infilename, char *header_id);
extern size_t get_sequence_length(char *infilename, char *header_id);
extern size_t get_sequence_length_with_pos(char *infilename, char *header_id, long int header_pos);


#endif /* FASTA_READER_H_ */
