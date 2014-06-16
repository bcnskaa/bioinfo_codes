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
 *
 *
 *
 *
 * vcf_extract.c
 *
 *  vcf_extract: a small tool for extracting data from a data in VCF format
 *  Compilation:
 *
 *  Intel icc:
 *  /opt/intel/compilerpro-12.0.0.084/bin/intel64/icc -o vcf_extract vcf_extract.c -O3 -Wall -m64 -g -D_FILE_OFFSET_BITS=64
 *
 *	GNU gcc:
 *	gcc -o vcf_extract vcf_extract.c -O3 -Wall -m64 -g -D_FILE_OFFSET_BITS=64
 *
 *
 *	Modifications:
 *		1. Construct a SNP id list with chromosomal ID and position
 *		2. Assume the data of input VCF file is sorted based on their chromosomal coordinates, therefore the program
 *		   maintains two pointers, one for the current SNP id and other one is for the VCF data.
 *		   Assumption:
 *		   		1. Coordinates were sorted both for SNP data and VCF data
 *		   		2. Both SNP and VCF data are in the same build
 *
 *		   Algorithms:
 *		   		VCFMapping(VCFdata, SNPdata)
 *		   		 1. i <- 0, n <- size(SNPdata), m <- size(VCFdata)
 *		   		 2. curSNPdata <- SNPdata[0]
 *		   		 3. curSNPpos <- curSNPdata.pos
 *		   		 4. while curVCFdata <- VCFdata(0 to m - 1) and i < n
 *		   		 5. 	curVCFpos <- curVCFdata.pos
 *				 6.		cflag <- 0
 *				 7.		while curSNPpos < curVCFpos and i < n and cflag = 0
 *				 8.			if curVCFdata eq curSNPdata
 *				 9.				OUT <- curVCFData, cflag <- 1
 *				10.			i++
 *				11.			curSNPdata <- SNPdata[i]
 *				12.			curSNPpos <- curSNPdata.pos
 *
 *  Created on: Jan 26, 2011
 *  Last modified: April 4, 2013
 *  Author: bcnskaa
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>

#define __DEFAULT_IO_BUFFER__ 16000000000
#define __DEFAULT_BUFFER_SIZE__ 4092		/* Size of the buffer for reading infile */
#define __DEFAULT_CHROM_ID_SIZE__ 10
#define __DEFAULT_VCF_FIELD_SIZE__  25
#define __DEFAULT_VCF_DELIMITER__ '\t'
#define __DEFAULT_VCF_COMMENT_CHR__ '#'
#define __DEFAULT_NEWLINE_CHR__ '\n'


typedef struct {
	char *rsid;
	char *chrom;
	unsigned int pos;
	char *desc;
} SNPData;


typedef struct {
	int EXPORT_HEADER_LINE;		/* flag for exporting a VCF header lines */
	char *OUTFILENAME;			/* out file name */
	char *IDS_INFILENAME;		/* infile name for snp ids */
	char *CHROM;				/* chromosome id to be extracted */
} PARAMS;




/**
 *
 */
int vcf_get_indel(infilename)
	char *infilename;
{
	int export_count;
	int processed_count, read_size, buffer_size;
	FILE *in;
	//off64_t prev_pos, cur_pos;
	off_t prev_pos, cur_pos;
	char *token, *buf;
	char *cur_vcf, *snp_id;
	char *chrom, *cur_chrom;
	int isIndel;
	const char delimiters[] = "\t";

	export_count = 0;

	buffer_size = 0;
	buf = 0;

	processed_count = 0;
	if(in = fopen64(infilename, "r")) {
		while((read_size = getline(&buf, &buffer_size, in)) != -1)
		{
			isIndel = 0;

			if(read_size > 3)
			{
				cur_vcf = (char*)strdup(buf);

				if(buf[0] == '#') { // Comment line
					// do nothing
					printf("%s", buf);
				} else {
					//printf("Buf: %s.\n", buf);

					//token = (char*)strchr(buf, newline);

					// Chromosome
					//token[0] = '\0';

					// CHROM
					token = (char*)strtok(buf, "\t ");
					//printf("Token: %s\n", token);
					// POS
					token = (char*)strtok(NULL, "\t ");
					//printf("Token: %s\n", token);
					// ID
					token = (char*)strtok(NULL, "\t ");
					//printf("Token: %s\n", token);
					snp_id = (char*)strdup(token);

					// REF
					token = (char*)strtok(NULL, "\t ");
					if(strlen(token) > 1) {
						isIndel = 1;
					}

					// ALT
					token = (char*)strtok(NULL, "\t ");
					if(strlen(token) > 1) {
						isIndel = 1;
					}


					if(isIndel) {
						//printf("%s is a polymorphic site with indel variants.\n", snp_id);
						printf("%s", cur_vcf);
					}

					free(cur_vcf); free(buf); free(snp_id);
					buffer_size = 0; buf = 0; cur_vcf = 0; snp_id = 0;
				}
			}
		}
	}

	return export_count;
}



int vcf_get_fpos_by_chrom(infilename)
	char *infilename;
{
	int processed_count, read_size, buffer_size;
	FILE *in;
	//off64_t prev_pos, cur_pos;
	off_t prev_pos, cur_pos;
	char *token, *buf;
	char *chrom, *cur_chrom;

	processed_count = 0;
	if(in = fopen64(infilename, "r")) {
		chrom = 0;
		cur_chrom = 0;
		prev_pos = 0; cur_pos = 0;
		buf = 0;
		buffer_size = 0;
		while((read_size = getline(&buf, &buffer_size, in)) != -1)
		{
			if(buf[0] == '#') { // Comment line
				if(prev_pos == 0) {
					cur_pos = ftello64(in);
					prev_pos = cur_pos;
				} else {
					prev_pos = cur_pos;
					cur_pos = ftello64(in);
				}
			} else {
				token = strchr(buf, '\t');
				cur_chrom = (char*)strndup(buf, token - buf);


				prev_pos = cur_pos;
				cur_pos = ftello64(in);


				if(chrom == 0) {
					chrom = (char*)strdup(cur_chrom);
					printf("%s: %lld\n", cur_chrom, cur_pos);
				} else {
					if(strcmp(cur_chrom, chrom) != 0) {
						printf("%s: %lld\n", cur_chrom, prev_pos);
						free(chrom);
						chrom = cur_chrom;
					} else {
						free(cur_chrom); cur_chrom = 0;
					}
				}

			}
		}
	}
	return processed_count;
}


/**
 *		   Algorithms:
 *		   		VCFMapping(VCFdata, SNPdata)
 *		   		 1. i <- 0, n <- size(SNPdata), m <- size(VCFdata)
 *		   		 2. curSNPdata <- SNPdata[0]
 *		   		 3. curSNPpos <- curSNPdata.pos
 *		   		 4. while curVCFdata <- VCFdata(0 to m - 1) and i < n
 *		   		 5. 	curVCFpos <- curVCFdata.pos
 *				 6.		cflag <- 0
 *				 7.		while curSNPpos < curVCFpos and i < n and cflag = 0
 *				 8.			if curVCFdata eq curSNPdata
 *				 9.				OUT <- curVCFData, cflag <- 1
 *				10.			i++
 *				11.			curSNPdata <- SNPdata[i]
 *				12.			curSNPpos <- curSNPdata.pos
 */
int vcf_extract_by_snp_info(infilename, params)
	char *infilename;
	PARAMS *params;
{
	int extract_count, i;
	char **snp_ids;
	char **snp_chroms;
	int *snp_poss;
	int snp_count;
	const char delimiters[] = "\t";
	const char newline = __DEFAULT_NEWLINE_CHR__;
	char *token, *token2, *buf;
	int buffer_size, read_size;
	int cflag;
	FILE *in, *out;

	/* Transient state of current working SNP data */
	char *cur_snp_id, *cur_snp_chrom, *cur_vcf_id, *cur_vcf_chrom;
	int cur_snp_pos, cur_vcf_pos;

	/**
	 * Read in the SNP IDs
	 */
	snp_ids = 0;
	snp_chroms = 0;
	snp_poss = 0;
	snp_count = 0;
	token = 0;
	buf = 0;

	extract_count = -1;
	printf("Processing SNP data from %s...\n", params->IDS_INFILENAME);
	if(in = fopen64(params->IDS_INFILENAME, "r")) {
		while((read_size = getline(&buf, &buffer_size, in)) != -1)
		{
			if(read_size > 3) {  /* The name of a normal snp id should have more than 3 characters */
				snp_ids = (char**)realloc(snp_ids, (snp_count + 1) * sizeof(char*));
				snp_poss = (int*)realloc(snp_poss, (snp_count + 1) * sizeof(int));
				snp_chroms = (char**)realloc(snp_chroms, (snp_count + 1) * sizeof(char*));
				token = (char*)strchr(buf, newline);
				token[0] = '\0';

				token = (char*)strtok(buf, delimiters);
				printf("RS: %s.\n", token);
				snp_ids[snp_count] = (char*)strdup(token);
				token = (char*)strtok(NULL, delimiters);
				printf("chr: %s.\n", token);
				snp_chroms[snp_count] = (char*)strdup(token);
				token = (char*)strtok(NULL, delimiters);
				printf("pos: %s.\n", token);
				snp_poss[snp_count] = atoi(token);

				printf("ID[%d]: %s (%s:%d).\n", snp_count, snp_ids[snp_count], snp_chroms[snp_count], snp_poss[snp_count]);
				snp_count++;

				buffer_size = 0; buf = 0;
			}
		}



		free(buf); buf = 0;
		buffer_size = 0;

		fclose(in);
		in = 0;
	}
	printf("Total number of SNP IDs: %d.\n", snp_count);


	if(snp_count == 0) {
		return extract_count;
	}


	/**
	 * Process VCF data
	 *
	 * VCF format:
	 * 	CHROM\tPOS\tID
	 */
	if(in = fopen64(infilename, "r")) {
		i = 0;
		cur_snp_id = snp_ids[i];
		cur_snp_pos = snp_poss[i];
		cur_snp_chrom = snp_chroms[i];
		cur_vcf_chrom = (char*)calloc(__DEFAULT_VCF_FIELD_SIZE__, sizeof(char));
		cur_vcf_id = (char*)calloc(__DEFAULT_VCF_FIELD_SIZE__, sizeof(char));
		cur_vcf_pos = 0;
		out = fopen64(params->OUTFILENAME, "w");

		while((read_size = getline(&buf, &buffer_size, in)) != -1 && i < snp_count)
		{
			if(read_size > 3) {
				if(buf[0] == __DEFAULT_VCF_COMMENT_CHR__) {  /* The name of a normal snp id should have more than 3 characters */
					fputs(buf, out);
				} else {
					token = (char*)strchr(buf, __DEFAULT_VCF_DELIMITER__);
					strncpy(cur_vcf_chrom, buf, token - buf);

					// chrom id from VCF does not contain "chr" characters,
					token = cur_snp_chrom;
					cur_snp_chrom++; cur_snp_chrom++; cur_snp_chrom++;
					if(strcmp(cur_vcf_chrom, cur_snp_chrom) == 0) {

						cur_snp_chrom = token;

						// Chr
						token = (char*)strchr(buf, __DEFAULT_VCF_DELIMITER__);
						token++;
						// Pos
						token2 = (char*)strchr(token, __DEFAULT_VCF_DELIMITER__);
						strncpy(cur_vcf_id, token, token2 - token);

						cur_vcf_pos = atoi(cur_vcf_id);

						memset(cur_vcf_id, '\0', __DEFAULT_VCF_FIELD_SIZE__);

						token2++;
						token = (char*)strchr(token2, __DEFAULT_VCF_DELIMITER__);
						strncpy(cur_vcf_id, token2, token - token2);

						//printf("OK, they're same stuff. (%s->%s:%d vs %s->%s:%d)\n", cur_snp_id, cur_snp_chrom, cur_snp_pos, cur_vcf_id, cur_vcf_chrom, cur_vcf_pos);


						cflag = 0;
						while(cur_snp_pos <= cur_vcf_pos && i < snp_count && cflag == 0) {
							if(strcmp(cur_snp_id, cur_vcf_id) == 0) {
								printf("Extracting %d: %s...\n", extract_count, cur_vcf_id);
								fputs(buf, out);
								cflag = 1;
								extract_count++;
							}

							i++;

							if(i < snp_count) {
								cur_snp_id = snp_ids[i];
								cur_snp_pos = snp_poss[i];
								cur_snp_chrom = snp_chroms[i];
							}
						}
					} else {
						//printf("Oops, they're not same stuff. (%s vs %s)\n", cur_snp_chrom, cur_vcf_chrom);
					}

					//buffer_size = 0; buf = 0;
				}
			}

			memset(cur_vcf_id, '\0', __DEFAULT_VCF_FIELD_SIZE__);
			memset(cur_vcf_chrom, '\0', __DEFAULT_VCF_FIELD_SIZE__);
		}

		printf("%d data exported to %s.\n", extract_count, params->OUTFILENAME);

		free(cur_vcf_chrom); cur_vcf_chrom = 0;
		free(cur_vcf_id); cur_vcf_id = 0;

		fclose(in);
		in = 0;
		fclose(out);
		out = 0;

		free(buf); buf = 0;
		buffer_size = 0;

		cur_snp_id = 0;
		cur_snp_pos = 0;
		cur_snp_chrom = 0;
		token = 0;
		token2 = 0;
	}


	for(i = 0; i < snp_count; i++) {
		free(snp_ids[i]); snp_ids[i] = 0;
		free(snp_chroms[i]); snp_chroms[i] = 0;
	}
	free(snp_chroms); snp_chroms = 0;
	free(snp_ids); snp_ids = 0;
	free(snp_poss); snp_poss = 0;


	return extract_count;
}



int vcf_extract_by_snp_ids(infilename, params)
	char *infilename;
	PARAMS *params;
{
	int extract_count, cflag, i, idx, j;
	char **ids, **rec_ids;
	unsigned int id_count;

	FILE *in, *out;
	char *buf, *buf2, *token, *cur_id, *buffer;
	int buffer_size, read_size;
	int processed_count;

	/* Set buffer size */
	buffer_size = __DEFAULT_BUFFER_SIZE__;
	read_size = 0;
	token = 0; buf = 0; in = 0;
	extract_count = -1;


	/**
	 * Read in the SNP IDs
	 */
	ids = 0; id_count = 0;
	if(in = fopen64(params->IDS_INFILENAME, "r")) {
		while((read_size = getline(&buf, &buffer_size, in)) != -1)
		{
			if(read_size > 3) {  /* The name of a normal snp id should have more than 3 characters */
				ids = (char**)realloc(ids, (id_count + 1) * sizeof(char*));
				token = (char*)strchr(buf, '\n');
				token[0] = '\0';

				ids[id_count] = (char*)strdup(buf);
				printf("ID[%d]: %s.\n", id_count, ids[id_count]);
				id_count++;
			}
		}
	}
	printf("Total number of SNP IDs: %d.\n", id_count);


	if(id_count == 0) {
		return extract_count;
	}
	/**
	 * Open the large VCF file
	 */
	if(in = fopen64(infilename, "r")) {
		out = fopen64(params->OUTFILENAME, "w");
		extract_count = 0;
		processed_count = 0;
		//buffer = (char*)calloc(1024, sizeof(char));
		while((read_size = getline(&buf, &buffer_size, in)) != -1)
		{
			if(buf[0] == __DEFAULT_VCF_COMMENT_CHR__) {
				fputs(buf, out);
			} else {
				//token = buf;
				//(char*)strtok(token, "\t");
				token = strchr(buf, __DEFAULT_VCF_DELIMITER__);
				token++;
				token = strchr(token, __DEFAULT_VCF_DELIMITER__);
				token++;
				buf2 = token;
				token++;
				token = strchr(token, __DEFAULT_VCF_DELIMITER__);


				buffer = (char*)strndup(buf2, token - buf2);
				cflag = 0; i = 0;
				while(cflag == 0 && i < id_count) {
					cur_id = ids[i];

					//printf("Current ID: %s (%s).\n", cur_id, buffer);
					//if(strncmp(buf2, cur_id, token - buf2) == 0) {
					if(strcmp(buffer, cur_id) == 0) {
						printf("ID found: %s (%d).\n", cur_id, extract_count);
						fputs(buf, out);
						extract_count++;
						cflag = 1;

						/**
						 * Reconstruct the id list
						 */
						rec_ids = (char**)calloc(id_count - 1, sizeof(char*));
						idx = 0; j = 0;
						while(idx < id_count) {
							if(idx != i) {
								rec_ids[j] = (char*)strdup(ids[idx]);
								j++;
							}

							free(ids[idx]); ids[idx] = 0;
							idx++;
						}

						free(ids);
						ids = rec_ids;
						id_count--;
						rec_ids = 0;

						break;
					}
					i++;
				}
				free(buffer); buffer = 0;


				if(processed_count % 100000 == 0) {
					printf("Processed count: %d\n", processed_count);
				}

				processed_count++;
			}
		}
		fclose(in); in = 0;
		fclose(out); out = 0;
		free(buf); buf = 0;
	}


	return extract_count;
}


int vcf_get_by_chrom_range(infilename, chrom, spos, epos)
char *infilename;
char *chrom;
int spos;
int epos;
{
	size_t count, i;
	FILE *in;
	char *buf, *token, *token2, *buf2;
	size_t buffer_size, cflag;
	char tmp[1024];
	size_t pos, read_size;

	count = 0;

	buf = 0;
	token = 0;
	buf2 = 0;
	buffer_size = 0;



	printf("File: %s", infilename);

	if(in = (FILE*)fopen64(infilename, "r")) {

		cflag = 0;
		while((read_size = getline(&buf, &buffer_size, in)) != -1)
		{

			if(buf[0] == '#') { // Comment line
				printf("%s", buf);
			} else if(buf[0] != '\n') {
				token = strchr(buf, '\t');


				strncpy(tmp, buf, token - buf);

				//printf("%s %u %u\n", tmp, token, buf);


				if(strncmp(tmp, chrom, strlen(chrom)) == 0) {
					token2 = strchr(token + 1, '\t');

					strncpy(tmp, token, token2 - token);


					pos = atoi(tmp);
					//printf("%u.\n", pos);

					if(pos >= spos && pos <= epos) {
						//printf("%s, %d.\n", chrom, pos);
						printf("%s", buf);
					}
				}
			}

			memset(tmp, '\0', 1024);
		}

		fclose(in); in = 0;

		free(buf); buf = 0;
	}


	return count;
}


//int vcf_extract_by_chrom(infilename, params, ptr_func_extract)
//	char *infilename;
//	PARAMS *params;
//	int (*ptr_func_extract)(char*);
int vcf_extract_by_chrom(infilename, params)
		char *infilename;
		PARAMS *params;
{
	int extract_count;
	FILE *in, *out;
	char *buf, *token, *buf2;
	int buffer_size, read_size, cflag;
	char *IO_BUFFER;

	/* Set buffer size */
	buffer_size = __DEFAULT_BUFFER_SIZE__;
	read_size = 0;
	token = 0; buf = 0; in = 0;
	extract_count = -1;


	/**
	 * Open the large VCF file
	 */
	if(in = fopen64(infilename, "r")) {
		out = fopen64(params->OUTFILENAME, "w");
		extract_count = 0;
		cflag = 0;
		//buf2 = (char*)calloc(__DEFAULT_CHROM_ID_SIZE__, sizeof(char));
		//memset(buf2, '\0', __DEFAULT_CHROM_ID_SIZE__);

		IO_BUFFER = (char*)calloc(__DEFAULT_IO_BUFFER__, sizeof(char));
		setvbuf(in, IO_BUFFER, _IOFBF, (size_t)__DEFAULT_IO_BUFFER__);
		while((read_size = getline(&buf, &buffer_size, in)) != -1)
		{
			if(buf[0] == '#') { // Comment line
				fputs(buf, out);
			} else {
				//token = buf;
				//(char*)strtok(token, "\t");
				token = strchr(buf, '\t');
				//buf2 = (char*)strndup(buf, token - buf);
				//strncpy(buf2, buf, token - buf);

				//printf("Buf: %s.\n", buf2);
				if(strncmp(buf, params->CHROM, token - buf) == 0) {
				//if(strcmp(buf2, params->CHROM) == 0) {
					fputs(buf, out);
					extract_count++;
//					if(cflag == 0) {
//						printf("Stating extracting for %s...\n", params->CHROM);
//						cflag = 1;
//					}
				} else {
					if(cflag == 1) {
//						printf("Finished extracting.\n");
						break;
					}
				}

			}
		}
		free(IO_BUFFER); IO_BUFFER = 0;
		//free(buf2); buf2 = 0;
		fclose(in); in = 0;
		fclose(out); out = 0;
		free(buf); buf = 0;
	}

	return extract_count;
}

void print_usage() {
	printf("Usage\n");
	printf("\t(1): by-chrom VCF-INFILE-NAME OUTFILE-NAME CHROM\n");
	printf("\t(2): by-snp-id VCF-INFILE-NAME OUTFILE-NAME SNP-ID-INFILE-NAME\n");
	printf("\t(3): fpos-by-chrom VCF-INFILE-NAME\n");
	printf("\t(4): by-chrom-range VCF-INFILE-NAME OUTFILE-NAME CHR1:10-10000\n");
	printf("\t(5): indel VCF-INFILE-NAME OUTFILE-NAME\n");
	printf("\t(6): indel-w-threshold VCF-INFILE-NAME OUTFILE-NAME THRESHOLD\n");

	printf("Description\n");
	printf("\t(1) Extract all the vcf elements in a specified chromosome.\n");
	printf("\t(2) Format of SNP-ID-INFILE-NAME:");
	printf("\t\tSNP-ID\\tCHROM\\tPOS\n");
	printf("\t(3) Return the file position of the first occurence of data of given chromosome.\n");
	printf("\t(4) Extract all the vcf elements of a specific range of a chromosome.\n");
	printf("\t(5) Extract all the vcf elements with one or more indel variant (>1bp).\n");
	printf("\t(6) Extract all the vcf elements with one or more indel variant (>1bp) at least threshold.\n");


}

int main(argc, argv)
	int argc;
	char **argv;
{
	int status;
	PARAMS *params;

	//printf("%s %d.\n", argv[1], argc);

	if(argc > 5 || argc == 3 || argc == 4) {

	} else {
		//printf("Usage (1): by-chrom VCF-INFILE-NAME OUTFILE-NAME CHROM\n");
		//printf("Usage (2): by-snp-id VCF-INFILE-NAME OUTFILE-NAME SNP-ID-INFILE-NAME\n");
		print_usage();
		exit(0);
	}

	status = 1;


	if(strcmp(argv[1], "by-chrom") == 0) {
		params = (PARAMS*)calloc(1, sizeof(PARAMS));
		params->EXPORT_HEADER_LINE = 1;
		params->OUTFILENAME = argv[3];
		params->CHROM = argv[4];
		params->IDS_INFILENAME = 0;

		printf("Extracting VCF data from %s (chr%s)...\n", argv[2], params->CHROM);

		//status = vcf_extract_by_chrom(argv[2], params, 0);
		status = vcf_extract_by_chrom(argv[2], params);
		printf("Export: %d.\n", status);

		free(params); params = 0;
	} else if(strcmp(argv[1], "by-snp-id") == 0) {
		params = (PARAMS*)calloc(1, sizeof(PARAMS));
		params->EXPORT_HEADER_LINE = 1;
		params->OUTFILENAME = argv[3];
		params->IDS_INFILENAME = argv[4];
		params->CHROM = 0;

		//status = vcf_extract_by_snp_ids(argv[2], params, 0);
		status = vcf_extract_by_snp_info(argv[2], params);
		//printf("Export: %d.\n", status);

		free(params); params = 0;
	} else if(strcmp(argv[1], "fpos-by-chrom") == 0) {
		vcf_get_fpos_by_chrom(argv[2]);
	} else if(strcmp(argv[1], "by-chrom-range") == 0) {
		vcf_get_by_chrom_range(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));
	} else if(strcmp(argv[1], "indel") == 0) {
		//vcf_get_fpos_by_chrom(argv[2]);
		vcf_get_indel(argv[2]);
	} else {
		print_usage();
		exit(0);
	}

	return status;
}
