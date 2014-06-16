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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// System Calls
#include <unistd.h>
#include <sys/stat.h> 
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <pthread.h>

#define IO_LARGE_BUFFER_SIZE 1000000000

/**
 * 2m29 
 * 1m16 mmap in count_read
 */

/**
 * 
 */
typedef struct {
	size_t read_length;
	int is_paired_end;
	int sequencing_platform;
	int quality_value_max;
	int quality_value_min;
} FASTQ_META;

typedef struct {
	char *header;	// header id
	char *seq;	// Sequence
	int *qv;	// Quality Value
} FASTQ_READ;


typedef struct {
	FASTQ_META *meta;
	FASTQ_READ **seq1;
	FASTQ_READ **seq2;
	size_t n;	
} FASTQ_PE;


/**
 * Phred quality value 
 * http://en.wikipedia.org/wiki/FASTQ_format
 *
 * 
 */

int PLATFORM_UNKNOWN = -1;
int PLATFORM_SANGER = 0;
int PLATFORM_SOLEXA = 1;
int PLATFORM_ILLUMINA13 = 2;
int PLATFORM_ILLUMINA15 = 3;
int PLATFORM_ILLUMINA18 = 4;

int SANGER_QV_OFFSET= '!';
int SANGER_QV_MIN = '!';
int SANGER_QV_MAX = 'I';

int SOLEXA_QV_OFFSET = '@';
int SOLEXA_QV_MIN = ';';
int SOLEXA_QV_MAX = 'h';
int ILLUMINA13_QV_OFFSET = '@';
int ILLUMINA13_QV_MIN = '@';
int ILLUMINA13_QV_MAX = 'h';
int ILLUMINA15_QV_OFFSET = '@';
int ILLUMINA15_QV_MIN = 'B';
int ILLUMINA15_QV_MAX = 'h';
int ILLUMINA18_QV_OFFSET = '!';
int ILLUMINA18_QV_MIN = '!';
int ILLUMINA18_QV_MAX = 'K';





/*
void print_PHRED_QV()
{
	int i;
	for(i = ILLUMINA_QV_OFFSET; i < (ILLUMINA_QV_OFFSET + 40); i++) {
		printf("%c: %d.\n", i, i);
	}
}
*/


int detect_platform(pe)
FASTQ_PE *pe;
{
	int i, j, isIdentified;
	int min, max;
	int sampling_n;
	char qv;
	
	sampling_n = 10000;
	min = 1000;
	max = -1;
	isIdentified = 0;
	for(i = 0; i < sampling_n && !isIdentified; i++) {
	//for(i = 0; i < pe->n && !isIdentified; i++) {
		//printf("H: %s\nS: %s\nQ: %s.\n Max: %d, Min: %d, l=%d\n", pe->seq1[i]->header,pe->seq1[i]->seq, pe->seq1[i]->qv, max, min, pe->meta->read_length);
		//printf("Max: %d, Min: %d\n", pe->seq1[i]->header,pe->seq1[i]->seq, pe->seq1[i]->qv, max, min);
		
		for(j = 0; j < pe->meta->read_length && !isIdentified; j++) {
			qv = pe->seq1[i]->qv[j];
			
			if(min > qv)
				min = qv;
			
			if(max < qv)
				max = qv;
		}
	}
	
	//printf("%c, %c\n", min, max);
	
	pe->meta->quality_value_max = max;
	pe->meta->quality_value_min = min;
	
	// Check SANGER
	if(min == SANGER_QV_MIN && max == SANGER_QV_MAX) {
		pe->meta->sequencing_platform = PLATFORM_SANGER;
		printf("PLATFORM_SANGER detected.\n");
		return PLATFORM_SANGER;
	} else if(min == SOLEXA_QV_MIN && max == SOLEXA_QV_MAX) {
		pe->meta->sequencing_platform = PLATFORM_SOLEXA;
		printf("PLATFORM_SOLEXA detected.\n");
		return PLATFORM_SOLEXA;
	} else if(min == ILLUMINA13_QV_MIN && max == ILLUMINA13_QV_MAX) {
		pe->meta->sequencing_platform = PLATFORM_ILLUMINA13;
		printf("PLATFORM_ILLUMINA13 detected.\n");
		return PLATFORM_ILLUMINA13;
	} else if(min == ILLUMINA15_QV_MIN && max == ILLUMINA15_QV_MAX) {
		pe->meta->sequencing_platform = PLATFORM_ILLUMINA15;
		printf("PLATFORM_ILLUMINA15 detected.\n");
		return PLATFORM_ILLUMINA15;			
	} else if(min == ILLUMINA18_QV_MIN && max == ILLUMINA18_QV_MAX) {
		pe->meta->sequencing_platform = PLATFORM_ILLUMINA18;
		printf("PLATFORM_ILLUMINA18 detected.\n");
		return PLATFORM_ILLUMINA18;	
	} else {
		pe->meta->sequencing_platform = PLATFORM_UNKNOWN;
		printf("PLATFORM_UNKNOWN detected.\n");
		return PLATFORM_UNKNOWN;
	}
}


/*
int indexing(char *) {
}
*/

void free_PE(pe)
FASTQ_PE *pe;
{
	int i;
	
	if(!pe) return;
	
	for(i = 0; i < pe->n; i++) {
		free(pe->seq1[i]->header); // pe->seq1[i]->header = 0;
		free(pe->seq1[i]->seq);
		free(pe->seq1[i]->qv);
		
		free(pe->seq2[i]->header); // pe->seq1[i]->header = 0;
		free(pe->seq2[i]->seq);
		free(pe->seq2[i]->qv);		
	}
	
	free(pe->seq1); pe->seq1 = 0;
	free(pe->seq2); pe->seq2 = 0;
	free(pe->meta); pe->meta = 0;
	free(pe); pe = 0;
}

int *convert_qv(sqv)
char *sqv; 
{
	int *iqv, i, len;
	
	len = strlen(sqv);
	iqv = (int*)calloc(sizeof(int), len);

	for(i = 0; i < len; i++) {
		iqv[i] = (int)(sqv[i]);
	}
	
	return iqv;
}


FASTQ_PE* process_mmap_fastq(fn1, fn2)
char* fn1;
char* fn2;
{
	FASTQ_PE *pe;
	int fastq_n;
	int n1, n2, i, j;
	char *fq1, *fq2, *fq1_line, *fq2_line;
	char *fq_start, *fq_end, temp[1000];
	int fd1, fd2;
	struct stat st1;
	struct stat st2;
	const char *NEW_LINE = "\n";
	pe = 0;
	fd1 = open64(fn1, O_RDONLY);
	fd2 = open64(fn2, O_RDONLY);	
	
	if(fd1 != -1 && fd2 != -1) {
    		if (fstat(fd1, &st1) != 0 || fstat(fd2, &st2) != 0)
			return 0;
        	
		printf("%s: %d, %s: %d.\n", fn1, st1.st_size, fn2, st2.st_size);
		
		
		// Memory mapping the two files
		fq1 = mmap(NULL, st1.st_size, PROT_READ, MAP_PRIVATE, fd1, 0);
		fq2 = mmap(NULL, st2.st_size, PROT_READ, MAP_PRIVATE, fd2, 0);
		
		// Close the file descriptor
		close(fd1);
		close(fd2);
		
		n1 = 0;
		i = 0;
		do {
			if(fq1[i] == '\n') n1++;
			i++;
		} while(i < st1.st_size);
		
		n2 = 0;
		i = 0;
		do {
			if(fq2[i] == '\n') n2++;
			i++;
		} while(i < st2.st_size);
		
		n1 = n1 / 4;
		n2 = n2 / 4;
		
		printf("Number of read: %s - %d, %s - %d.\n", fn1, n1, fn2, n2);
		if(n1 != n2 || n2 == 0 || n1 == 0 )
		{
			printf("Read file size not matched (%d, %d), abort.\n", n1, n2);
			return 0;
		}
		
		
		pe = (FASTQ_PE*)calloc(sizeof(FASTQ_PE), 1);
		pe->meta = (FASTQ_META*)calloc(sizeof(FASTQ_META), 1);
		pe->seq1 = 0;
		pe->seq2 = 0;
		pe->n = n1;
		
		pe->seq1 = (FASTQ_READ**)calloc(sizeof(FASTQ_READ*), n1);
		pe->seq2 = (FASTQ_READ**)calloc(sizeof(FASTQ_READ*), n2);

		printf("Preparing FASTQ object...");
		i = 0;			
		while(i < pe->n) {
			pe->seq1[i] = (FASTQ_READ*)calloc(sizeof(FASTQ_READ), 1);
			pe->seq2[i] = (FASTQ_READ*)calloc(sizeof(FASTQ_READ), 1);
			i++;	
		}
		printf("done.\n");
		
		//printf("Processing FASTQ object (%s)...\n");
		printf("Processing FASTQ object (%s)...", fn1);
		
		i = 0;
		
		fq_start = fq1;	
		while(i < pe->n && (fq_end = strchr(fq_start + 1, '\n')) != NULL) {
			
			//pe->seq1[i]->header = (char*)strndup(fq_start, 1 + (fq_end - fq_start));
			pe->seq1[i]->header = calloc(sizeof(char), (fq_end - fq_start) + 1);
			strncpy(pe->seq1[i]->header, fq_start, fq_end - fq_start);
			
			//pe->seq1[i]->header[fq_end - fq_start] = '\0';
			
			//printf("header: %s, %d - %d.\n", pe->seq1[i]->header, fq_start, fq_end);
			fq_start = fq_end + 1; fq_end = strchr(fq_start, '\n') ;
			//pe->seq1[i]->seq = (char*)strndup(fq_start, 1 + (fq_end - fq_start));
			//pe->seq1[i]->seq[fq_end - fq_start] = '\0';
			pe->seq1[i]->seq = calloc(sizeof(char), (fq_end - fq_start) + 1);
			strncpy(pe->seq1[i]->seq, fq_start, fq_end - fq_start);
			
			
			//printf("Sequence: %s.\n", pe->seq1[i]->seq);
			
			fq_start = fq_end + 1; fq_end = strchr(fq_start, '\n') ;
			fq_start = fq_end + 1; fq_end = strchr(fq_start, '\n') ;
			//pe->seq1[i]->qv = (char*)strndup(fq_start, 1 + (fq_end - fq_start));
			//pe->seq1[i]->qv[fq_end - fq_start] = '\0';
			//printf("QV: %s.\n", pe->seq1[i]->qv);
			
			
			//pe->seq1[i]->qv = calloc(sizeof(char), (fq_end - fq_start) + 1);
			//strncpy(pe->seq1[i]->qv, fq_start, fq_end - fq_start);
			
			strncpy(temp, fq_start, fq_end - fq_start);
			pe->seq1[i]->qv = convert_qv(temp);
			
			
			fq_start = fq_end + 1;
			i++;
		}
		
		
		/*
		fq1_line = strtok(fq1, NEW_LINE);
		while(fq1_line != NULL) {
			pe->seq1[i]->header = (char*)strdup(fq1_line);
			
			fq1_line = strtok(NULL, NEW_LINE);
			printf("header: %s.\n", fq1_line);
			
			pe->seq1[i]->seq = (char*)strdup(fq1_line);
			fq1_line = strtok(NULL, NEW_LINE);
			fq1_line = strtok(NULL, NEW_LINE);
			//pe->seq1[i]->qv = (char*)strdup(fq1_line);
			pe->seq1[i]->qv = convert_qv(fq1_line);
			
			i++;
			fq1_line = strtok(NULL, NEW_LINE);
		}
		*/
		
		printf("done (%d).\n", i);
		
		
		
		printf("Processing FASTQ object (%s)...", fn2);	
		
		i = 0;
		fq_start = fq2;	
		while(i < pe->n && (fq_end = strchr(fq_start + 1, '\n')) != NULL) {
			pe->seq2[i]->header = calloc(sizeof(char), (fq_end - fq_start) + 1);
			strncpy(pe->seq2[i]->header, fq_start, fq_end - fq_start);
			
			fq_start = fq_end + 1; fq_end = strchr(fq_start, '\n') ;
			
			pe->seq2[i]->seq = calloc(sizeof(char), (fq_end - fq_start) + 1);
			strncpy(pe->seq2[i]->seq, fq_start, fq_end - fq_start);
			
			fq_start = fq_end + 1; fq_end = strchr(fq_start, '\n') ;
			fq_start = fq_end + 1; fq_end = strchr(fq_start, '\n') ;
			
			//pe->seq2[i]->qv = calloc(sizeof(char), (fq_end - fq_start) + 1);
			//strncpy(pe->seq2[i]->qv, fq_start, fq_end - fq_start);
			
			strncpy(temp, fq_start, fq_end - fq_start);
			pe->seq2[i]->qv = convert_qv(temp);
			
			
			fq_start = fq_end + 1;
			i++;
		}
		
		/*
		fq2_line = strtok(fq2, NEW_LINE);
		while(fq2_line != NULL) {
			pe->seq2[i]->header = (char*)strdup(fq2_line);
			
			fq2_line = strtok(NULL, NEW_LINE);
			
			
			pe->seq2[i]->seq = (char*)strdup(fq2_line);
			fq2_line = strtok(NULL, NEW_LINE);
			fq2_line = strtok(NULL, NEW_LINE);
			//pe->seq2[i]->qv = (char*)strdup(fq2_line);
			pe->seq2[i]->qv = convert_qv(fq2_line);
			i++;
			fq2_line = strtok(NULL, NEW_LINE);
		}
		*/
		
		printf("done (%d).\n", i);

		// Set meta data
		pe->meta->is_paired_end = 1;
		pe->meta->read_length = strlen(pe->seq1[0]->seq);
		
		
		printf("Detecting sequencing platform...");
		detect_platform(pe);
		//printf("done.\n");
	
		munmap((void*)fq1, st1.st_size);
		munmap((void*)fq2, st2.st_size);
	}
	
	return pe;
		
}


//int process_fastq(fn)
FASTQ_PE* process_fastq(fn, fn2)
char* fn;
char* fn2;
{
	int fastq_n;
	int n1, n2;
	
	FILE *pfile, *pfile2;	
	
	FASTQ_PE *pe;
	
	//FASTQ_READ *read1, *read2;
	
	char *line_buffer, *token, *buffer;
	char *line_buffer2, *token2, *buffer2;	
	//char *header, *seq, *qv;
	//char *header2, *seq2, *qv2;
	//size_t line_count;
	size_t line_buffer_size, line_buffer_size2;
		
	
	
	// Initialize variables
	line_buffer_size = 0;
	line_buffer = 0;
	//line_count = 0;
	fastq_n = 0;



	// Estimate the number of reads
	n1 = get_read_count(fn) / 4;
	n2 = get_read_count(fn2) / 4;
		
	printf("Number of read: %s - %d, %s - %d.\n", fn, n1, fn2, n2);
	if(n1 != n2 || n2 == 0 || n1 == 0 )
	{
		printf("Read file size not matched (%d, %d), abort.\n", n1, n2);
		return 0;
	}
	

	
	// Make a connection to a local file
	//pfile = (FILE*)fopen64(fn, "r");
	//pfile2 = (FILE*)fopen64(fn2, "r");
	
	pfile = (FILE*)fopen64(fn, "r");
	pfile2 = (FILE*)fopen64(fn2, "r");
	setvbuf(pfile, NULL, _IOFBF, IO_LARGE_BUFFER_SIZE);
	setvbuf(pfile2, NULL, _IOFBF, IO_LARGE_BUFFER_SIZE);
	
	if(pfile && pfile2) {
		pe = (FASTQ_PE*)calloc(sizeof(FASTQ_PE), 1);
		pe->meta = (FASTQ_META*)calloc(sizeof(FASTQ_META), 1);
		pe->seq1 = 0;
		pe->seq2 = 0;
		pe->n = 0;
		
		pe->seq1 = (FASTQ_READ**)calloc(sizeof(FASTQ_READ*), n1);
		pe->seq2 = (FASTQ_READ**)calloc(sizeof(FASTQ_READ*), n2);
		
		// Read in header line of each read
		while(getline(&line_buffer, &line_buffer_size, pfile) > 0 && getline(&line_buffer2, &line_buffer_size2, pfile2) > 0) 
		{
			//line_count++;
			
			
			//pe->seq1 = (FASTQ_READ**)realloc(pe->seq1, sizeof(FASTQ_READ*) * (pe->n + 1));
			//pe->seq2 = (FASTQ_READ**)realloc(pe->seq2, sizeof(FASTQ_READ*) * (pe->n + 1));

			pe->seq1[pe->n] = (FASTQ_READ*)calloc(sizeof(FASTQ_READ), 1);
			pe->seq2[pe->n] = (FASTQ_READ*)calloc(sizeof(FASTQ_READ), 1);

			// Fill in the buffer
			pe->seq1[pe->n]->header = (char*)strdup(line_buffer);
			pe->seq2[pe->n]->header = (char*)strdup(line_buffer2);
			
			//token = (char*)strtok(line_buffer, "\t");
			
			// Seq
			getline(&line_buffer, &line_buffer_size, pfile);
			line_buffer[strlen(line_buffer) - 1] = 0;
			pe->seq1[pe->n]->seq = (char*)strdup(line_buffer);
			getline(&line_buffer2, &line_buffer_size2, pfile2);
			line_buffer2[strlen(line_buffer2) - 1] = 0;
			pe->seq2[pe->n]->seq = (char*)strdup(line_buffer2);
			
				
			// Wash line
			getline(&line_buffer, &line_buffer_size, pfile);
			getline(&line_buffer2, &line_buffer_size2, pfile2);
			
			// qv
			getline(&line_buffer, &line_buffer_size, pfile);
			//pe->seq1[pe->n]->qv = (char*)strdup(line_buffer);
			pe->seq1[pe->n]->qv = convert_qv(line_buffer);
			getline(&line_buffer2, &line_buffer_size2, pfile2);
			//pe->seq2[pe->n]->qv = (char*)strdup(line_buffer2);			
			pe->seq2[pe->n]->qv = convert_qv(line_buffer2);
			
			if(pe->n == 0)
			{
				pe->meta->is_paired_end = 1;
				pe->meta->read_length = strlen(pe->seq1[pe->n]->seq);
			}
			
			
			pe->n++;
			
			
			/*
			
			free(header); header = 0;
			free(seq); seq = 0;
			free(qv); qv = 0;
			
			free(header2); header2 = 0;
			free(seq2); seq2 = 0;
			free(qv2); qv2 = 0;
			*/				
		}
		
		detect_platform(pe);
		//pe->meta->sequencing_platform = detect_platform(pe);
		
		fclose(pfile);
		pfile = 0;
		
		fclose(pfile2);
		pfile2 = 0;
	}
	
	
	//return fastq_n;
	return pe;
}

/*
int PLATFORM_UNKNOWN = -1;
int PLATFORM_SANGER = 0;
int PLATFORM_SOLEXA = 1;
int PLATFORM_ILLUMINA13 = 2;
int PLATFORM_ILLUMINA15 = 3;
int PLATFORM_ILLUMINA18 = 4;
*/
void print_qv_table(pe, fn1, fn2) 
FASTQ_PE *pe;
char *fn1;
char *fn2;
{
	int i, j;
	FILE *pfile, *pfile2;
	
	pfile = (FILE*)fopen(fn1, "w");
	pfile2 = (FILE*)fopen(fn2, "w");

	for(i = 0; i < pe->meta->read_length; i++) {
		if(i ==0) {
			fprintf(pfile, "%d", i);
			fprintf(pfile2, "%d", i);
		} else {
			fprintf(pfile, "%d", i);
			fprintf(pfile2, "\t%d", i); 
		}
	}
	fputc('\n', pfile);
	fputc('\n', pfile2);
	
	
	if(pe->meta->sequencing_platform == PLATFORM_UNKNOWN || pe->meta->sequencing_platform == PLATFORM_SANGER || pe->meta->sequencing_platform == PLATFORM_ILLUMINA18) 
	{
		for(i = 0; i < pe->n; i++) {
			for(j = 0; j < pe->meta->read_length; j++) {
				if(j == 0) {
					fprintf(pfile, "%d", pe->seq1[i]->qv[j]);
					fprintf(pfile2, "%d", pe->seq2[i]->qv[j]);
				} else {
					fprintf(pfile, "\t%d", pe->seq1[i]->qv[j]);
					fprintf(pfile2, "\t%d", pe->seq2[i]->qv[j]);
				}
			}
			fputc('\n', pfile);
			fputc('\n', pfile2);
		}
	} else {
		for(i = 0; i < pe->n; i++) {
			for(j = 0; j < pe->meta->read_length; j++) {
				if(j == 0) {
					fprintf(pfile, "%d", pe->seq1[i]->qv[j] - SOLEXA_QV_OFFSET);
					fprintf(pfile2, "%d", pe->seq2[i]->qv[j] - SOLEXA_QV_OFFSET);
				} else {
					fprintf(pfile, "\t%d", pe->seq1[i]->qv[j] - SOLEXA_QV_OFFSET);
					fprintf(pfile2, "\t%d", pe->seq2[i]->qv[j] - SOLEXA_QV_OFFSET);
				}
			}
			fputc('\n', pfile);
			fputc('\n', pfile2);
		}
	}
	
	
	fclose(pfile); pfile = 0;
	fclose(pfile2); pfile2 = 0;
}

// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
/**
 * pthread properties
 * 
 * compiling command: gcc -D_REENTRANT read_mapper.c -o read_mapper -O3 -g -lpthread
 * 
 */
pthread_mutex_t mutex;
int current_processing_read_pos = 0;
int thread_n = 32;
/*
typedef struct {
	int tid;
	int pid;
	FASTQ_PE *pe;
	int **qv_table1;
	int **qv_table2;	
	int qv_level;
} THREAD_INFO;
*/

typedef struct {
	int tid;
	int pid;
	FASTQ_PE *pe;
	int **qv_table1;
	int **qv_table2;	
	int qv_level;
	int spos;
	int epos;
} THREAD_INFO;

void* pthread_process_qv(void *arg)
{
	int i, j, cur_pos, qv_idx;
	THREAD_INFO *info;
	FASTQ_PE *pe;
	i = 0;
	
	info = (THREAD_INFO*)arg;	
	pe = info->pe;
	
	printf("Running %d: %d.\n", info->tid, current_processing_read_pos);
	
	
	for(j = info->spos; j < info->epos; j++) {
		cur_pos = j;
		printf("Thread %d: processing position %d (%d-%d).\n", info->tid, cur_pos, info->spos, info->epos);
		for(i = 0; i < pe->n; i++)
		{
			qv_idx = (pe->seq1[i]->qv[cur_pos] - pe->meta->quality_value_min);
			
			if(qv_idx < 0 || qv_idx >= info->qv_level)
				printf("Abnormal at %d: %d.\n", i, qv_idx);
			
			(info->qv_table1[qv_idx][cur_pos])++;
			
			qv_idx = (pe->seq2[i]->qv[cur_pos] - pe->meta->quality_value_min);
			
			if(qv_idx < 0 || qv_idx >= info->qv_level)
				printf("Abnormal at %d: %d.\n", i, qv_idx);
			
			(info->qv_table2[qv_idx][cur_pos])++;		
			
		}	
	}
	pthread_exit((void*) 0);
		
	/*
	while(1) {
		pthread_mutex_lock(&mutex);
		if(current_processing_read_pos < pe->meta->read_length) {
			cur_pos = current_processing_read_pos;
			current_processing_read_pos++;
		} else {
			pthread_mutex_unlock(&mutex);	
			pthread_exit((void*) 0);
		}
		pthread_mutex_unlock(&mutex);
	
	
		printf("Thread %d: processing position %d.\n", info->tid, cur_pos);
		for(i = 0; i < pe->n; i++)
		{
			qv_idx = (pe->seq1[i]->qv[cur_pos] - pe->meta->quality_value_min);
			
			if(qv_idx < 0 || qv_idx >= info->qv_level)
				printf("Abnormal at %d: %d.\n", i, qv_idx);
			
			(info->qv_table1[qv_idx][cur_pos])++;
			
			qv_idx = (pe->seq2[i]->qv[cur_pos] - pe->meta->quality_value_min);
			
			if(qv_idx < 0 || qv_idx >= info->qv_level)
				printf("Abnormal at %d: %d.\n", i, qv_idx);
			
			(info->qv_table2[qv_idx][cur_pos])++;		
			
		}
	}
	*/
}

void pthread_print_qv_summary(pe, fn1, fn2) 
FASTQ_PE *pe;
char *fn1;
char *fn2;
{
	int i, j , qv_level, cur_pos, step, remain;
	int **qv_table1, **qv_table2;
	void *status;
	pthread_t threads[thread_n];
	pthread_attr_t attr;
	THREAD_INFO **info;
	FILE *pfile, *pfile2;
	
	pfile = (FILE*)fopen(fn1, "w");
	pfile2 = (FILE*)fopen(fn2, "w");

	// Print header
	fputs("QV\t", pfile);
	fputs("QV\t", pfile2);
	for(i = 0; i < pe->meta->read_length; i++) {
		if(i ==0) {
			fprintf(pfile, "%d", i);
			fprintf(pfile2, "%d", i);
		} else {
			fprintf(pfile, "\t%d", i);
			fprintf(pfile2, "\t%d", i); 
		}
	}
	fputc('\n', pfile);
	fputc('\n', pfile2);
	

	qv_level = (pe->meta->quality_value_max - pe->meta->quality_value_min) + 2;
	
	printf("QV level: %d.\n", qv_level);
	
	// initialize qv table
	qv_table1 = (int**)calloc(sizeof(int*), qv_level);
	for(i = 0; i < qv_level; i++)
		qv_table1[i] = (int*)calloc(sizeof(int), pe->meta->read_length);
	
	
	
	// initialize qv table
	qv_table2 = (int**)calloc(sizeof(int*), qv_level);
	for(i = 0; i < qv_level; i++)
		qv_table2[i] = (int*)calloc(sizeof(int), pe->meta->read_length);
	
	
	
	
	pthread_mutex_init(&mutex, NULL);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);


	cur_pos = 0; step = pe->meta->read_length / thread_n;
	remain = pe->meta->read_length % thread_n;
	printf("Preparing %d multithreading kernel executing %d (%d) positions each...\n", thread_n, step, pe->meta->read_length % thread_n);
	info = (THREAD_INFO**)calloc(sizeof(THREAD_INFO*), thread_n);
	for(i = 0; i < thread_n; i++) {
		info[i] = (THREAD_INFO*)calloc(sizeof(THREAD_INFO), 1);
		info[i]->tid = i;
		info[i]->pid = i;
		info[i]->pe = pe;
		info[i]->qv_table1 = qv_table1;
		info[i]->qv_table2 = qv_table2;		
		info[i]->qv_level = qv_level;
		info[i]->spos = cur_pos;
		if(remain> 0)
			info[i]->epos = cur_pos + step + 1;
		else 
			info[i]->epos = cur_pos + step;
		
		printf("Creating thread (id: %d)...\n", info[i]->tid);
		//pthread_create(&threads[i],  NULL, pthread_process_qv, (void*)info[i]);
		pthread_create(&threads[i], &attr, pthread_process_qv, (void*)info[i]);
		
		cur_pos = info[i]->epos;
		remain--;
	} 
	//printf("done.\n");
	
	
	(void)pthread_attr_destroy(&attr);
	
	for(i = 0; i < thread_n; i++) {
		pthread_join(threads[i], &status);
	} 
	
	/**
	 * print the table
	 */
	printf("Summarizing...");
	// Export the values
	for(i = 0; i < qv_level; i++) {
		for(j = 0; j < pe->meta->read_length; j++) {
			if(j == 0){
				fprintf(pfile, "%d\t%d", i, qv_table1[i][j]);
				fprintf(pfile2, "%d\t%d", i, qv_table2[i][j]);
			}else {
				fprintf(pfile, "\t%d", qv_table1[i][j]);
				fprintf(pfile2, "\t%d", qv_table2[i][j]);
			}
		}
		
		fputc('\n', pfile);
		fputc('\n', pfile2);
	}
	printf("done.\n");
	 
	 
	
	
		 
	
	for(i = 0; i < thread_n; i++) {
		free(info[i]); info[i] = 0;
	} 
	free(info); info = 0;


	
	for(i = 0; i < qv_level; i++) {
		free(qv_table1[i]);qv_table1[i] = 0;
	}
	free(qv_table1); qv_table1 = 0;
	
	for(i = 0; i < qv_level; i++) {
		free(qv_table2[i]);qv_table2[i] = 0;
	}
	free(qv_table2); qv_table2 = 0;	
	
	
	
	

	pthread_mutex_destroy(&mutex);
	//pthread_exit(NULL);	
}
// ************************************************************************************
// ************************************************************************************
// ************************************************************************************


void print_qv_summary(pe, fn1, fn2) 
FASTQ_PE *pe;
char *fn1;
char *fn2;
{
	int i, j, qv_idx, qv_level;
	
	int **qv_table1, **qv_table2;
	FILE *pfile, *pfile2;
	
	pfile = (FILE*)fopen(fn1, "w");
	pfile2 = (FILE*)fopen(fn2, "w");

	// Print header
	fputs("QV\t", pfile);
	fputs("QV\t", pfile2);
	for(i = 0; i < pe->meta->read_length; i++) {
		if(i ==0) {
			fprintf(pfile, "%d", i);
			fprintf(pfile2, "%d", i);
		} else {
			fprintf(pfile, "\t%d", i);
			fprintf(pfile2, "\t%d", i); 
		}
	}
	fputc('\n', pfile);
	fputc('\n', pfile2);
	
	
	qv_level = (pe->meta->quality_value_max - pe->meta->quality_value_min) + 2;
	
	printf("QV level: %d.\n", qv_level);
	
	// initialize qv table
	qv_table1 = (int**)calloc(sizeof(int*), qv_level);
	for(i = 0; i < qv_level; i++) {
		qv_table1[i] = (int*)calloc(sizeof(int), pe->meta->read_length);
		/*for(j = 0; j < pe->meta->read_length; j++) {
			qv_table1[i][j] = 0;
		}
		*/
	}
	
	
	
	// initialize qv table
	qv_table2 = (int**)calloc(sizeof(int*), qv_level);
	for(i = 0; i < qv_level; i++) {
		qv_table2[i] = (int*)calloc(sizeof(int), pe->meta->read_length);
		/*for(j = 0; j < pe->meta->read_length; j++) {
			qv_table2[i][j] = 0;
		}
		*/
	}
	
	
	printf("Summarizing...");
	
	
	// Summarize quality values over every position of every reads	
	for(i = 0; i < pe->n; i++) {
		for(j = 0; j < pe->meta->read_length; j++) {
			qv_idx = (pe->seq1[i]->qv[j] - pe->meta->quality_value_min);
			
			if(qv_idx < 0 || qv_idx >= qv_level)
				printf("Abnormal at %d: %d.\n", i, qv_idx);
			
			(qv_table1[qv_idx][j])++;
			qv_idx = (pe->seq2[i]->qv[j] - pe->meta->quality_value_min);
			if(qv_idx < 0 || qv_idx >= qv_level)
				printf("Abnormal at %d: %d.\n", i, qv_idx);	
			
			(qv_table2[qv_idx][j])++;
		}
	}
	
	
	printf("done.\n");
	
	// Export the values
	for(i = 0; i < qv_level; i++) {
		for(j = 0; j < pe->meta->read_length; j++) {
			if(j == 0){
				fprintf(pfile, "%d\t%d", i, qv_table1[i][j]);
				fprintf(pfile2, "%d\t%d", i, qv_table2[i][j]);
			}else {
				fprintf(pfile, "\t%d", qv_table1[i][j]);
				fprintf(pfile2, "\t%d", qv_table2[i][j]);
			}
		}
		
		fputc('\n', pfile);
		fputc('\n', pfile2);
	}
	
	
	
	
	for(i = 0; i < qv_level; i++) {
		free(qv_table1[i]);qv_table1[i] = 0;
	}
	free(qv_table1); qv_table1 = 0;
	
	for(i = 0; i < qv_level; i++) {
		free(qv_table2[i]);qv_table2[i] = 0;
	}
	free(qv_table2); qv_table2 = 0;	
	
	
	fclose(pfile); pfile = 0;
	fclose(pfile2); pfile2 = 0;
}


void print_qv(pe, idx)
FASTQ_PE *pe;
int idx;
{
	int i;
	
	for(i = 0; i < pe->meta->read_length; i++) {
		printf(" %d", pe->seq1[idx]->qv[i] - pe->meta->quality_value_min);
	}
	printf("\n");
}

int get_read_count(fn)
char *fn;
{
	int n, c, i;
	FILE *pfile;
	int fd;
	char *ptr;
	struct stat st; 

	n = -1;
	


	fd = open64(fn, O_RDONLY);
	
	
	if(fd != -1) {
    		if (fstat(fd, &st) != 0)
			return -1;
        	
		printf("Size: %d.\n", st.st_size);
		
		//setvbuf(pfile, NULL, _IOFBF, st.st_size);
		n = 0;
		i = 0;
		ptr = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
		close(fd);
		
		do {
			if(ptr[i] == '\n') n++;
			i++;
		} while(i < st.st_size);
		
		/*
		do {
			c = fgetc(pfile);
			if(c == '\n') n++;
		} while(c != EOF);
		*/
		munmap(ptr, st.st_size);
		
		//flush();	
		//fclose(pfile); pfile = 0;
	}
	
	return n;
}

int main(argc, argv)
	int argc;
	char** argv;
{
	FASTQ_PE *pe;
	char cmd[2048];

	
	if(argc != 5) {
		printf("Use: %s read1.fq read2.fq read1.qv read2.qv\n", argv[0]);
		exit(0);
	}



	printf("Processing %s and %s...\n", argv[1], argv[2]);
	//pe = process_fastq(argv[1], argv[2]);
	pe = process_mmap_fastq(argv[1], argv[2]);
	
	if(!pe) {
		printf("Problem on processing the input files. Abort now.\n");
		return -1;
	}
	
	printf("Number of reads: %d, read length: %d.\n", pe->n, pe->meta->read_length);
	
	//print_PHRED_QV();
	
	//printf("Exporting quality values to %s and %s...\n", argv[3], argv[4]);
	//print_qv_table(pe, argv[3], argv[4]);
	
	
	printf("Exporting summary of quality values to %s and %s...\n", argv[3], argv[4]);	
	//print_qv_summary(pe, argv[3], argv[4]);
	pthread_print_qv_summary(pe, argv[3], argv[4]);
	
	//print_qv(pe, 1000);
	
	printf("Clean up...");
	free_PE(pe); pe = 0;
	printf("done.");
	
	pthread_exit(NULL);	
	return 0;
} 

/*
 *
 ~/share/tools/axel-2.4/axel -n8 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG02545/sequence_read/ERR184137_1.filt.fastq.gz
 ~/share/tools/axel-2.4/axel -n8 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG02545/sequence_read/ERR184137_2.filt.fastq.gz
 ls -la
 folder: gpu10.ust.hk/data_local10/1KGP
 *
 */
