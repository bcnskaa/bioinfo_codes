/*
 *  Created on: Mar 3, 2012
 *      Author: bcnskaa
 */
#ifndef _VCF41_H_
#define _VCF41_H_

#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <vector>
#include <map>

using namespace std;

#define GENDER_UNKNOWN 0;
#define FEMALE 2;
#define MALE 1;

/**
 * Because of the std::to_string function, to compile the file including this header should use the flag -std=c++0x
 *
 */


/**
 * Population information for VCF data
 *
 * The infile should contain the following fields delimited by a tab character:
 * 1. Population ID (GBK, JPT, ...etc)
 * 2. Individual ID (HG00081... etc)
 * 3. SRA Accession number
 * 4. Gender (female or male)
 */
typedef struct {
	char *individual_id;
	int gender;
	char *population_id;
	int phenotype;
	char *fid;		// Family ID
	char *father;
	char *mother;
} VCF_INDIVIDUAL_INFO;


typedef struct {
	char *population_id; 				// population id
	vector<VCF_INDIVIDUAL_INFO*> individual_infos;

	/*vector<int> vcf_col_idx;			// column position of individual's genotype in the vcf data
	vector<string> vcf_col_names; 		// individual name
	vector<bool> gender;				// false - female, true - male
	*/
} VCF_POPULATION;


class VCF_POPULATION_INFO {
private:
	void process(const char* filename);
public:
	map<string, VCF_INDIVIDUAL_INFO*> map_individual2individualInfo;
	map<string, VCF_POPULATION*> map_population;		// column position of individual's genotype in the vcf data
	//map<string, VCF_POPULATION*> map_individual2population;				// tell which population individual belonging to
	map<string, string> map_individual2population;				// tell which population individual belonging to

	vector<string> individual_ids;
	vector<VCF_POPULATION*> populations;
	vector<string> population_ids;
	vector<VCF_INDIVIDUAL_INFO*> individual_infos;

	VCF_POPULATION_INFO(string infilename);
	VCF_POPULATION_INFO(const char *infilename);

	~VCF_POPULATION_INFO();

	//bool add_individual(char *individual_id, char *population_id, char *gender)
	bool add_individual(string ind_id, string pop_id, string gender);
	VCF_INDIVIDUAL_INFO* get_individual_info(string ind_id);
};





typedef struct {
	string *id;
	char *chrom;
	size_t pos;
	char *ref_allele;
	int ref_allele_n;
	char *alt_allele;
	int alt_allele_n;
	char *genotype;
	int size;
} VCF_DATA;

/**
 * Class for storing VCF data
 */
class VCF {
private:
	map<string, VCF_DATA*> map_site2vcf;
	map<string, int> map_individual2cindex;
	vector<string> individual_ids;
	vector<string> site_ids;
	vector<VCF_DATA*> vcf_data;

	/*  */
	int process_vcf(const char *infilename);
	int process_vcf_option; // 0 - do nothing, 1 - discard redundant, 2 - merge redundant


public:

	/**
	 * If nothing is read, the individual and site number will be zero
	 */
	VCF(string infilename);
	VCF(const char *infilename);
	VCF(const char *infilename, bool discardRedundant); // True: If more than one site having same chrom:pos value, only the first one will keep, False: try to merge the site
	~VCF();

	VCF_DATA* get_vcf(string site_id);
	VCF_DATA* get_vcf(int i);  // i = the index of site

	bool isReady();

	string get_ref_allele(string site_id);
	string get_alt_allele(string site_id);

	// i = site index
	string get_major_allele(int i);
	string get_major_allele(string site_id);
	string get_major_allele(VCF_POPULATION *population, string site_id);

	double get_ref_allele_freq(int i);
	double get_ref_allele_freq(string site_id);
	double get_ref_allele_freq(VCF_POPULATION *population, string site_id);

	/**
	 * type: 0 - VCF original coding, 1 - Allele
	 * Converted genotypes will be stored to the variable "genotypes"
	 */
	int convert_genotype(vector<string> &genotypes, VCF_DATA *vcf_data, int type);

	string get_site_id(int i);
	string get_individual_id(int i);

	vector<string> get_site_ids_by_range(size_t spos, size_t epos);
	vector<string> get_site_ids_by_idx(int i_start, int i_end);


	/* Get the number of site */
	size_t site_num();
	/* Get the number of individual in the dataset */
	size_t individual_num();


	bool is_individual_found(string individual_id);
	bool is_site_found(string site_id);

	/*  */
	int get_genotype(string site_id, string individual_id);
	int get_genotype(int site_idx, int ind_idx);
};

#endif
