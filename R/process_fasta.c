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
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdio.h>
#include <ctype.h>

#include "SeqTools/fasta_reader.h"

// http://adv-r.had.co.nz/C-interface.html
/**
 * 1. read fasta sequence
 * 2.
 *
 * Compilation: R CMD SHLIB process_fasta.c SeqTools/fasta_reader.c
 */
//SEXP process_fasta_r(SEXP fasta_fn, SEXP fasta_header, SEXP window_size, SEXP char_to_be_counted)
SEXP process_fasta_r(SEXP fasta_fn, SEXP fasta_header, SEXP char_to_be_counted)
{
	SEXP poss;
	char *fn, *header;
	SeqObj *so;
	char *token, *cur_pos;
	int i, j, token_n;
	int search_len;
	int match_n;
	//char *seq_uc;

	fn = CHAR(STRING_ELT(fasta_fn, 0));
	header = CHAR(STRING_ELT(fasta_header, 0));
	token = CHAR(STRING_ELT(char_to_be_counted, 0));

	printf("File name: %s\n", fn);

	so = read_fasta_by_id(fn, header);

	poss = 0;
	if(so != 0)
	{
		token_n = strlen(token);
		printf("Token length: %d, %s\n", token_n, token);

		search_len = (so->seq_len);
		match_n = 0;
		//so->seq = toupper(so->seq);

		//seq_uc = (char*)seq_toupper(so);
		seq_toupper(so);
		printf("Length of %s: %d bp\n", so->header, so->seq_len);
		if(token_n == 1)
		{
			for(i = 0; i < so->seq_len; i++)
			{
				if(so->seq[i] == token[0])
					match_n++;
				//if(seq_uc[i] == token[0])
				//	match_n++;
			}
		} else {
			for(i = 0, cur_pos = so->seq; i < search_len; i++, cur_pos++)
			{
				if(strncmp(cur_pos, token, token_n) == 0)
					match_n++;
			}
		}
		printf("There are %d matches of \"%s\" found in %s.\n", match_n, token, header);

		// Allocate memory
		if(match_n > 0)
		{
			PROTECT(poss = allocVector(INTSXP, match_n));

			j = 0;
			if(token_n == 1)
			{
				for(i = 0; i < so->seq_len; i++)
				{
					if(so->seq[i] == token[0])
					{
						INTEGER(poss)[j] = i;
						j++;
					}
				}
			} else {
				for(i = 0, cur_pos = so->seq; i < search_len; i++, cur_pos++)
				{
					if(strncmp(cur_pos, token, token_n) == 0)
					{
						INTEGER(poss)[j] = i;
						j++;
					}
				}
			}
		}
	} else {
		printf("No object named, %s, is found.\n", header);
	}


	UNPROTECT(1);

	if(so != 0)
	{
		destroy_SeqObj(so);
		//free(seq_uc); seq_uc = 0;
	}
	so = 0;

	return poss;
}


//
//SEXP get_subseqs_r(SEXP fasta_fn, SEXP fasta_header, SEXP ss, SEXP es)
//{
//	char *fn, *header;
//	int spos, epos;
//
//	fn = CHAR(STRING_ELT(fasta_fn, 0));
//	header = CHAR(STRING_ELT(fasta_header, 0));
//	spos = INTEGER(s);
//	epos = INTEGER(e);
//
//	return(get_subseq(fn, header, spos, epos));
//}


char *get_subseq(char *fasta_fn, char *fasta_header, int spos, int epos)
{
	SeqObj *so;
	char *seq;

	seq = 0;

	so = read_fasta_by_id(fasta_fn, fasta_header);

	if(so != 0)
	{
		seq = get_subseq_SeqObj(so, spos, epos);
		destroy_SeqObj(so);
		//free(seq_uc); seq_uc = 0;
	}
	so = 0;

	return(seq);
}



SEXP get_subseq_r(SEXP fasta_fn, SEXP fasta_header, SEXP s, SEXP e)
{
	char *fn, *header, *seq;
	int spos, epos;
	SEXP subseq;

	fn = CHAR(STRING_ELT(fasta_fn, 0));
	header = CHAR(STRING_ELT(fasta_header, 0));
	spos = INTEGER(s);
	epos = INTEGER(e);

	seq = get_subseq(fn, header, spos, epos);
	subseq = 0;
	if(seq != 0)
		subseq = mkChar(seq);

	return(subseq);
}

