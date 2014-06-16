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
 * fasta_reader.c
 *
 *  Created on: Oct 30, 2009
 *      Author: bcnskaa
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "fasta_reader.h"


size_t CHAR_PER_LINE = 50;
const int HEADER_NOT_FOUND = -1;


size_t buffer_size;
char *buffer;
int cflag;
int status;
int byte_read;



SeqObj* read_fasta_by_id_fast(infilename, header_id)
char *infilename;
char *header_id;
{
	int i;
	SeqObj *so;
	FILE *infile;
	char *ptr;
	long int header_pos;
	size_t seq_len;
	size_t cur_seq_pos;

	header_pos = get_header_pos(infilename, header_id);
	if(header_pos == HEADER_NOT_FOUND)
	{
		printf("%s cannot be found in %s, stop now.\n", header_id, infilename);
		return 0;
	}

	//printf("Header pos=%zd.\n", header_pos);
	seq_len = get_sequence_length_with_pos(infilename, header_id, header_pos);

	infile = (FILE*)fopen64(infilename, "r");
	cflag = 0;
	so = 0;
	if(infile != 0 && seq_len > 0)
	{
		buffer_size = 4096;
		buffer = (char*)calloc(buffer_size, sizeof(char));

		so = (SeqObj*)calloc(1, sizeof(SeqObj));
		so->header = (char*)strdup(header_id);
		so->seq_len = seq_len;
		so->seq = (char*)calloc(so->seq_len + 1, sizeof(char));
		memset(so->seq, '\0', so->seq_len + 1);
		cur_seq_pos = 0;

		// forward fstream to the sequence position
		fseek(infile, header_pos, SEEK_SET);

		while(feof(infile) == 0 && cur_seq_pos < seq_len)
		{
			// read line
			byte_read = getline(&buffer, &buffer_size, infile);

			//printf("read=%s.", buffer);
			if(byte_read > 1)
			{
				byte_read -= 1;
				for(i = 0; i < byte_read; i++)
				{
					//printf("cur_seq_pos=%d\n", cur_seq_pos);
					so->seq[cur_seq_pos] = buffer[i];
					cur_seq_pos++;
				}
			}
		}


	    free(buffer); buffer = 0;
	    buffer_size = 0;
	    fclose(infile); infile = 0;
	}


	return so;
}


SeqObj* read_fasta_by_id(infilename, header_id)
	char *infilename;
	char *header_id;
{
	int i, len, process_n;
	SeqObj *so;
	FILE *infile;
	char *ptr;
	long int header_pos;


	header_pos = get_header_pos(infilename, header_id);
	if(header_pos == HEADER_NOT_FOUND)
	{
		printf("%s cannot be found in %s, stop now.\n", header_id, infilename);
		return 0;
	}



	infile = (FILE*)fopen64(infilename, "r");

	printf("Reading from %s...\n", infilename);

	cflag = 0;
	process_n = 0;
	so = 0;
	if(infile != 0)
	{
		buffer_size = 4096;
		buffer = (char*)calloc(buffer_size, sizeof(char));


	    while(feof(infile) == 0 && cflag < 2)
	    {
	    	status = getline(&buffer, &buffer_size, infile);

	    	if(status > 1) {
		    	buffer[status - 1] = '\0'; buffer_size--; // Remove newline character


	            if(buffer[0] == '>') { // New header encounter

	            	ptr = buffer; ptr++; // exclude '>'

	            	if(strcmp(header_id, ptr) == 0) {
	            		printf("%s found\n", header_id);

	            		so = (SeqObj*)calloc(1, sizeof(SeqObj));
	            		so->header = (char*)strdup(ptr);
	            		so->seq = 0;
	            		so->seq_len = 0;


//	            		// Preestimate the size of the fasta object
//	            		filepos = ftell(infile);
//	            		 while(feof(infile) == 0 && cflag < 2)
//	            		 {
//	            			 status = getline(&buffer, &buffer_size, infile);
//	            			 if(statu)
//	            		 }



	            		cflag = 1;
	            	} else {


	            		if(cflag == 1)
	            			cflag = 2;
	            		else
	            			cflag = 0;
	            	}
	            	process_n++;
	            } else {
		    		if(cflag == 1) {
		    			len = status - 1;

		    			so->seq = (char*)realloc(so->seq, (so->seq_len + len) * sizeof(char));
		    			ptr = so->seq + so->seq_len;
		    			strncpy(ptr, buffer, len);
		    			so->seq_len += len;



		    			/*
		    			len = status - 1;
		    			//printf("Len: %d.\n", len);
		    			//so->seq = (char*)realloc(so->seq, (so->seq_len + strlen(buffer)) * sizeof(char));
		    			so->seq = (char*)realloc(so->seq, (so->seq_len + len + 1) * sizeof(char));
		    			ptr = so->seq + so->seq_len;
		    			//ptr = (char*)strdup(buffer);
		    			strncpy(ptr, buffer, len);
		    			*ptr++ = '\0';

		    			//printf("%s: %d.\n", buffer, strlen(buffer));
		    			//printf("%s.\n", so->seq);
		    			so->seq_len += len;
		    			*/
		    		}
	            }
	    	}
	    }

	    free(buffer); buffer = 0;
	    buffer_size = 0;
	    fclose(infile); infile = 0;

	    printf("Number of object processed: %d\n", process_n);
	}


	return so;
}



size_t get_sequence_length(infilename, header_id)
	char *infilename;
	char *header_id;
{
	long int header_pos;
	//FILE *infile;
	//size_t seq_len;


	//seq_len = 0;
	header_pos = get_header_pos(infilename, header_id);
	return get_sequence_length_with_pos(infilename, header_id, header_pos);


//	//printf("Header Position=%d.\n", header_pos);
//	if(header_pos != HEADER_NOT_FOUND)
//	{
//		buffer_size = 4096;
//		buffer = (char*)calloc(buffer_size, sizeof(char));
//
//		infile = (FILE*)fopen64(infilename, "r");
//
//		// forward fstream to the sequence position
//		fseek(infile, strlen(header_id) + 2, SEEK_SET);
//
//		while(feof(infile) == 0)
//		{
//			// read line
//			byte_read = getline(&buffer, &buffer_size, infile);
//
//			if(byte_read > 1)
//			{
//				buffer[byte_read - 1] = '\0';
//
//				printf("read=%s\n", buffer);
//
//				if(buffer[0] == '>') { // New header encounter
//					break;
//				} else {
//					seq_len += byte_read - 1;
//				}
//			}
//		}
//
//		// Clean up
//		free(buffer); buffer = 0;
//		fclose(infile); infile = 0;
//	}
//
//	return seq_len;
}

size_t get_sequence_length_with_pos(infilename, header_id, header_pos)
	char *infilename;
	char *header_id;
	long int header_pos;
{
	FILE *infile;
	size_t seq_len;


	seq_len = 0;
	//printf("Header Position=%d.\n", header_pos);
	if(header_pos != HEADER_NOT_FOUND)
	{
		buffer_size = 4096;
		buffer = (char*)calloc(buffer_size, sizeof(char));

		infile = (FILE*)fopen64(infilename, "r");

		// forward fstream to the sequence position
		//fseek(infile, header_pos + strlen(header_id) + 2, SEEK_SET);
		fseek(infile, header_pos, SEEK_SET);

		while(feof(infile) == 0)
		{
			// read line
			byte_read = getline(&buffer, &buffer_size, infile);

			if(byte_read > 1)
			{
				buffer[byte_read - 1] = '\0';

				//printf("read=%s\n", buffer);

				if(buffer[0] == '>') { // New header encounter
					break;
				} else {
					seq_len += byte_read - 1;
				}
			}
		}

		// Clean up
		free(buffer); buffer = 0;
		fclose(infile); infile = 0;
	}

	return seq_len;
}


//
//
//int check_if_header_exist(infilename, header_id)
//	char *infilename;
//	char *header_id;
//{
//	int i;
//	FILE *infile;
//	char *ptr;
//
//	infile = (FILE*)fopen64(infilename, "r");
//
//	cflag = 0;
//	if(infile != 0)
//	{
//		buffer_size = 4096;
//		buffer = (char*)calloc(buffer_size, sizeof(char));
//
//		while(feof(infile) == 0 && cflag < 2)
//		{
//			status = getline(&buffer, &buffer_size, infile);
//
//			if(status > 1)
//			{
//				buffer[status - 1] = '\0'; buffer_size--; // Remove newline character
//
//				if(buffer[0] == '>') { // New header encounter
//
//					ptr = buffer; ptr++; // exclude '>'
//
//					if(strcmp(header_id, ptr) == 0) {
//						cflag = 3;
//					}
//				}
//			}
//		}
//
//		free(buffer); buffer = 0;
//		fclose(infile); infile = 0;
//	}
//
//	if(cflag == 3)
//		return 1;
//	else
//		return 0;
//}


long int get_header_pos(infilename, header_id)
	char *infilename;
	char *header_id;
{
	int i;
	FILE *infile;
	char *ptr;
	long int pos;

	infile = (FILE*)fopen64(infilename, "r");

	pos = HEADER_NOT_FOUND;
	cflag = 0;

	if(infile != 0)
	{
		buffer_size = 4096;
		buffer = (char*)calloc(buffer_size, sizeof(char));

		while(feof(infile) == 0 && cflag < 2 && pos == HEADER_NOT_FOUND)
		{
			status = getline(&buffer, &buffer_size, infile);

			if(status > 1)
			{
				buffer[status - 1] = '\0'; buffer_size--; // Remove newline character

				if(buffer[0] == '>') { // New header encounter

					ptr = buffer; ptr++; // exclude '>'

					if(strcmp(header_id, ptr) == 0) {
						pos = ftell(infile);
					}
				}
			}
		}

		free(buffer); buffer = 0;
		fclose(infile); infile = 0;
	}

	return pos;
}




void print_fasta(outfilename, so)
	char *outfilename;
	SeqObj *so;
{
		FILE *outfile;
		int i, line_num, char_printed, seq_len;
		char *ptr;
		char buffer[CHAR_PER_LINE + 1];

		outfile = (FILE*)fopen(outfilename, "w");
		memset(buffer, '\0', CHAR_PER_LINE + 1);

		if(!so) {
			printf("Invalid SeqObj.\n");
			return;
		}

		seq_len = so->seq_len;

		if(outfile) {
			if((seq_len % CHAR_PER_LINE) == 0) {
				line_num = (int)(seq_len / CHAR_PER_LINE);
			} else {
				line_num = (int)((seq_len - (seq_len % CHAR_PER_LINE)) / CHAR_PER_LINE) + 1;
			}


			//printf("Line: %d", line_num);

			// Print header
			fputc('>', outfile);
			fputs(so->header, outfile);
			//fputc('\n', outfile);
			for(i = 0; i < so->seq_len; i++)
			{
				fputc(so->seq[i], outfile);
				//if(i != 0 && i % CHAR_PER_LINE == 0)
				//	fputc('\n', outfile);
			}

			/*
			ptr = so->seq;
			//for(i = 0, char_printed = 0; i < so->seq_len; i++) {
			for(i = 0, char_printed = 0; i <= line_num; i++) {
				strncpy(buffer, ptr, CHAR_PER_LINE);
				fputs(buffer, outfile);
				fputc('\n', outfile);

				ptr += CHAR_PER_LINE;

//				fputc(so->seq[i], outfile);
//				if(char_printed == CHAR_PER_LINE) {
//					char_printed = 0;
//					fputc('\n', outfile);
//				} else {
//					char_printed++;
//				}
			}
*/
			fclose(outfile); outfile = 0;
		}
}

void destroy_SeqObj(so)
	SeqObj *so;
{
	if(!so) return;

	free(so->header); so->header = 0;
	free(so->seq); so->seq = 0;
	free(so); so = 0;
}


char* get_subseq_SeqObj(so, spos, epos)
	SeqObj *so;
	int spos;
	int epos;
{
	char *seq;
	int i, j;
	size_t seq_len;

	seq = 0;
	if(spos >= 0 && epos < so->seq_len)
	{
		seq_len = (epos - spos) + 1;
		seq = (char*)calloc(seq_len, sizeof(char));

		for(j = 0, i = spos - 1; i < epos && i < so->seq_len; i++, j++)
			seq[j] = so->seq[i];
	}

	return seq;
}

char charAt_SeqObj(so, pos)
	SeqObj* so;
	int pos;
{
	char n;

	n = 0;
	if(pos >= 0 && pos < so->seq_len) {
		n = so->seq[pos - 1];
	}

	return n;
}

void display_fasta(so, spos, epos)
	SeqObj *so;
	int spos;
	int epos;
{
	char *seq;
	if((seq = get_subseq_SeqObj(so, spos, epos)) != 0) {
		printf("%s:%d..%d\n%s", so->header, spos, epos, seq);
	}
//	int i;
//	size_t seq_len;
//
//
//	if(strcmp(header, so->header) == 0)
//	{
//		printf("%s:%d..%d\n", header, spos, epos);
//		for(i = spos - 1; i < epos && i < so->seq_len; i++)
//		{
//			printf("%c", so->seq[i]);
//		}
//	}
}


//char* seq_toupper(so)
void seq_toupper(so)
	SeqObj* so;
{
	//int seq_n;
	int i;
	char c;
	//char *seq_uc;

	//seq_uc = (char*)calloc(so->seq_len, sizeof(char));

	for(i = 0; i < so->seq_len; i++)
	{
		c = so->seq[i];
		so->seq[i] = toupper(c);
	}

	//return i;
	//return seq_uc;
}

