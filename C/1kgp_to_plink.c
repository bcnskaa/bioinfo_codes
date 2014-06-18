/*
 * 1kgp_to_plink.cpp
 *
 *  Created on: Jun 3, 2013
 *      Author: bcnskaa
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <map>
#include <vector>
#include <string>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <fstream>



#include <boost/thread.hpp>




#include "vcf41.h"



//@headers = ();
//$REF_IDX = 3;
//$ALT_IDX = 4;
//$GENOTYPE_IDX = 9;
//$RSID_IDX = 2;
//$CHROM_IDX = 0;
//$POS_IDX = 1;
//@genotypes = ();
//$ind_n = 0;
//if(defined IN) {
//        #open PED_OUT, ">", $fn_ped;
//        open MAP_OUT, ">", $fn_map;
//        while(<IN>) {
//                chomp;
//
//                if(/^##/) { # Header line ignore.
//                } elsif(/^#C/) { # Column headers
//                        #print qq{$_\n};
//                        @headers = split /\t/, $_;
//                } else {
//                        #@items = split /\t/, $_;
//
//                        #push @genotypes, \@items;
//                        push @genotypes, [split /\t/];
//                }
//        }
//        close IN;
//
//        $ind_n = scalar @headers - $GENOTYPE_IDX;
//        print "Individual number: " . $ind_n . "\n";
//        print "Number of loci: " . scalar @genotypes . "\n";
//#       for(my $i = 0; $i < $ind_n; $i++){
//#               my $genotype_idx = $i + $GENOTYPE_IDX;
//#               print PED_OUT $headers[$genotype_idx];
//#               print "Processing " . $headers[$genotype_idx] . "\n";
//#
//#               my $ped_str = "";
//#               for(my $j = 0; $j < scalar @genotypes; $j++) {
//#                       my @snp_genotypes = @{$genotypes[$j]};
//#                       #print scalar @snp_genotypes . "\n";
//#               #foreach(@genotypes) {
//#                       #my @snp_genotypes = @$_;
//#                       my $genotype_str = $snp_genotypes[$genotype_idx];
//#                       #print substr($genotype_str, 0, 1) . "\n";
//#                       my $a1 = substr($genotype_str, 0, 1) + 1;
//#                       my $a2 = substr($genotype_str, 2, 1) + 1;
//#                       #print PED_OUT "  " . $a1 . " " . $a2;
//#                       $ped_str = $ped_str . "  " . $a1 . " " . $a2;
//#               }
//#               print PED_OUT $ped_str . "\n";
//#       }
//
//        my $unknown_id_n = 0;
//        for(my $j = 0; $j < scalar @genotypes; $j++) {
//                my @snp_genotypes = @{$genotypes[$j]};
//
//                my $chr = $snp_genotypes[$CHROM_IDX];
//                my $rs = $snp_genotypes[$RSID_IDX];
//                my $pos = $snp_genotypes[$POS_IDX];
//
//                if(length($rs) == 1) {
//                        $rs = "unknown_" . $unknown_id_n;
//                        $unknown_id_n++;
//                }
//                print MAP_OUT $chr . " " . $rs . " " . $pos . "\n";
//        }
//
//        close IN;
//        #close PED_OUT;
//        close MAP_OUT;
//        print scalar @genotypes . "\n";
//        @genotypes = ();
//}



typedef struct {
	string pop_id;
	vector<string> ind_ids;
	vector<string> site_ids;
	string ped_outfilename;
	string map_outfilename;
} JOB;

boost::mutex io_mutex;
vector<JOB*> job_list;
map<boost::thread::id, int> map_tid;
int THREAD_N = 8;
vector<boost::thread*> thread_pool;
bool isThreadEnabled = false;
bool isTerminated = false;


VCF *vcf;
VCF_POPULATION_INFO *vcf_info;


int export_data(VCF *vcf, VCF_POPULATION_INFO *vcf_info, vector<string> ind_ids, vector<string> site_ids, string ped_outfilename, string map_outfilename);


void thread_process()
{
	int tid = map_tid.find(boost::this_thread::get_id())->second;


	//boost::posix_time::seconds workTime(time);
	boost::mutex::scoped_lock lock(io_mutex);
	cout <<  "[Processing Core " << tid << "] " << "Start to run" << endl;
	lock.unlock();


	for(; job_list.size() > 0 && !isTerminated;)
	{
		boost::mutex::scoped_lock lock(io_mutex);
		JOB *job = job_list.back();
		cout <<  "[Processing Core " << tid << "] " << job->pop_id << ":" << job->ind_ids.size() << ":" << job->site_ids[0] << "-" << job->site_ids[job->site_ids.size() - 1] << endl;
		job_list.pop_back();
		lock.unlock();

		//
		export_data(vcf, vcf_info, job->ind_ids, job->site_ids, job->ped_outfilename, job->map_outfilename);
	}

}

void clear_thread_resource( )
{
	if(isThreadEnabled)
	{
		for(int i = 0; i < job_list.size(); i++)
		{
			JOB *job = job_list[i];
			job->ind_ids.clear();
			job->site_ids.clear();
			delete job;
		}
		job_list.clear();

		map_tid.clear();

		for(int i = 0; i < THREAD_N; i++)
		{
			boost::thread *thread = thread_pool[i];
			delete thread;
		}
		thread_pool.clear();
	}
}


/**
 *
 */
typedef struct {
	/*
	size_t spos;
	size_t epos;
	string block_ctx;
	size_t locus_n;
	*/
	string det_line;
	string block_line;
} HAPBLOCK;




/**
 *
 */
void reduce_map(vector<string> det_file_list, vector<string> block_file_list)
{
	vector<HAPBLOCK*> blocks;
	map<string, int> map_id2block;
	ifstream det_in, block_in;
	ofstream det_out, block_out;
	string det_in_line, block_in_line;


	if(det_file_list.size() != block_file_list.size()) {
		cout << "[reduce_map] det file and block file are not matched. Abort now."<< endl;
		return;
	}

	/**
	 * Assume the lists are sorted alphabetically
	 */
	for(int i = 0; i < det_file_list.size() ; i++) {
		det_in.open(det_file_list[i].c_str());		// with header
		block_in.open(block_file_list[i].c_str());	// withtout header line



		getline(det_in, det_in_line);				// flush the header line
		int j = 0;

		// Read all the data into a temp block
		vector<HAPBLOCK*> tmp_blocks;
		while(!det_in.eof() && !block_in.eof())
		{
			getline(det_in, det_in_line);
			getline(block_in, block_in_line);

			//
			HAPBLOCK *hb = new HAPBLOCK();

			hb->det_line = det_in_line;
			hb->block_line = block_in_line;

			tmp_blocks.push_back(hb);

			j++;
		}


		/**
		 *
		 */
		if(i == 0)	// The size of blocks is 0, simply add all the items inside tmp_blocks into blocks
		{
			//int j = 0;
			// for (vector<int>::iterator it = tmp_blocks.begin() ; it != tmp_blocks.end(); ++it)
			for(int j = 0; i < tmp_blocks.size(); j++)
			{
				blocks.push_back(tmp_blocks[j]);
				map_id2block.insert(pair<string, int>(tmp_blocks[j]->block_line, j));

				//j++;
			}
		} else {	//
			int j = map_id2block.find(tmp_blocks[0]->block_line) != map_id2block.end() ? map_id2block.find(tmp_blocks[0]->block_line)->second : -1;

			if(j != -1) {
				for (vector<HAPBLOCK*>::iterator it = tmp_blocks.begin() ; it != tmp_blocks.end(); ++it)
				{
					if(j < blocks.size()) {
						delete blocks[j];
						blocks[j] = *it;
					} else {
						map_id2block.insert(pair<string, int>(tmp_blocks[j]->block_line, j));
					}
					j++;
				}

			} else {
				cout << "Warning: No overlapping region exists for " << det_file_list[i] << " and " << block_file_list[i] << ", [block_line=" << tmp_blocks[0]->block_line << "]"<< endl;
			}

			tmp_blocks.clear();

		}

	}


	// clear the resource
	for (std::vector<HAPBLOCK*>::iterator it = blocks.begin() ; it != blocks.end(); ++it)
		delete *it;

	map_id2block.clear();
}





void proccess()
{

	if(isThreadEnabled)
	{
		// Prepare the thread pool
		for(int i = 0; i < THREAD_N; i++)
		{
			boost::thread *t = new boost::thread(thread_process);
			map_tid.insert(pair<boost::thread::id, int>(t->get_id(), i));
			thread_pool.push_back(t);
		}


		// Run all the job
		for(int i = 0; i < THREAD_N; i++)
			thread_pool[i]->join();
	}
}



/**
 *
 */
int export_data(VCF *vcf, VCF_POPULATION_INFO *vcf_info, vector<string> ind_ids, vector<string> site_ids, string ped_outfilename, string map_outfilename)
{
	ofstream ped_out, map_out;
	int export_n = 0;


	ped_out.open(ped_outfilename.c_str());
	for(int i = 0; i < ind_ids.size(); i++) {
		string ind_id = ind_ids[i];

		VCF_INDIVIDUAL_INFO* ind_info = vcf_info->get_individual_info(ind_id);

		if(ind_info == 0)
		{
			cout << ind_id << " does not have a valid meta data available, default missing values will be put." << endl;
			// Family ID\tIndividual ID\tGender Code (0: missing, 1: male, 2: female)\tPhenotype (0: unaffected, 1: affected, -9: unknown )
			ped_out << ind_id << "\t" << ind_id << "\t0" <<  "\t0";
		} else {
			if(ind_info->fid != 0)
				ped_out << ind_info->fid;
			else
				ped_out << ind_id;

			ped_out << "\t" << ind_id << "\t" << ind_info->gender << "\t0";
		}

		string line_buf;
		line_buf.clear();

		for(int j = 0; j < site_ids.size(); j++) {
			string site_id = site_ids[j];
			string ref_a, alt_a;

			char g1, g2;
			char geno = vcf->get_genotype(site_id, ind_id);

			//cout << geno << endl;

			if(geno == '0') {
				line_buf.append("  1 1");
				//ped_out << "  1 1";
				//g1 = '0'; g2 = '0';
			} else if(geno == '1') {
				line_buf.append("  1 2");
				//ped_out << "  1 2";
				//g1 = '0'; g2 = '1';
			} else if(geno == '2') {
				line_buf.append("  2 1");
				//ped_out << "  2 1";
				//g1 = '1'; g2 = '0';
			} else if(geno == '3') {
				line_buf.append("  2 2");
				//ped_out << "  2 2";
				//g1 = '1'; g2 = '1';
			} else { // Missing
				line_buf.append("  0 0");
				//ped_out << "  0 0";
				//g1 = '?'; g2 = '?';
			}
			export_n++;
		}
		ped_out << line_buf << endl;
		//ped_out << endl;
		ped_out.flush();
		line_buf.clear();
	}
	ped_out.close();


	map_out.open(map_outfilename.c_str());
	for(int i = 0; i < site_ids.size(); i++) {
		string site_id = site_ids[i];

		VCF_DATA* vcf_data = vcf->get_vcf(site_id);

		if(vcf_data)
			map_out <<  vcf_data->chrom << " "<< *vcf_data->id << " " << vcf_data->pos <<  endl;
		else
			cout << i << ") Problem with the site " << site_id << "." << endl;
	}
	map_out.close();

	return export_n;
}


 /**
  * Break the data into chunk and generate a run script
  *
  * map size defined in terms of number of sites
  */
int generate_split_map(VCF *vcf, VCF_POPULATION_INFO *vcf_info, int map_size, int overlapping_size, string outfile_prefix, string population_id)
{
	ofstream ped_out, map_out, run_out, list_out;
	vector<string> ind_ids, site_ids;
	int export_n = 0;

	string run_fn = outfile_prefix, ped_fn = outfile_prefix, map_fn = outfile_prefix, list_fn = outfile_prefix;

	cout << "[generate_split_map] Total number of site: " << vcf->site_num() << endl;
	cout << "[generate_split_map] Map size: " << map_size << " with overlapping size " << overlapping_size << endl;
	cout << "[generate_split_map] Number of map: " << (ceil(vcf->site_num() / (map_size - overlapping_size)) + 1) << endl;
	cout << "[generate_split_map] Creating command lines to " << run_fn << endl;


	if(population_id.empty())  { // Assume all individuals in the dataset are included
		// Create a list of individuals
		ind_ids.clear();
		for(int i = 0; i < vcf->individual_num(); i++) {
			ind_ids.push_back(vcf->get_individual_id(i));
		}
	} else {
		// Create a list of individuals
		ind_ids.clear();

		VCF_POPULATION *population = vcf_info->map_population.find(population_id)->second;
		if(population == 0)
		{
			cout << "[generate_split_map] No population included: " << outfile_prefix << ", abort now." << endl;
 	 		return 0;
		}

		for(int i = 0; i < population->individual_infos.size(); i++) {
			string ind_id = population->individual_infos[i]->individual_id;

			if(vcf->is_individual_found(ind_id))
				ind_ids.push_back(vcf->get_individual_id(i));
		}
	}

	// Check if there is any individual included
	if(ind_ids.size() > 0) {
		cout << "[generate_split_map] Number of individual: " << ind_ids.size() << endl;
	} else { // No individual included
		cout << "[generate_split_map] No individual included: " << outfile_prefix << ", abort now." << endl;
		return 0;
	}

 	run_fn.append(".run");
	run_out.open(run_fn);

	list_fn.append(".list");
 	list_out.open(list_fn);


 	// Create a list of site
 	for(int i = 0, j = 0; i < vcf->site_num() && j < vcf->site_num();) {
		if(i + map_size < vcf->site_num())
 			j = i + map_size;
 		else
 			j = vcf->site_num();

 		site_ids = vcf->get_site_ids_by_idx(i, j);


 		// Create a new files
 		ped_fn = outfile_prefix; ped_fn.append("."); ped_fn.append(std::to_string((long long int)i)); ped_fn.append("-"); ped_fn.append(std::to_string((long long int)j)); ped_fn.append(".ped");
 		map_fn = outfile_prefix; map_fn.append("."); map_fn.append(std::to_string((long long int)i)); map_fn.append("-"); map_fn.append(std::to_string((long long int)j)); map_fn.append(".map");

 		cout << "\t" << export_n  << ": creating " << ped_fn << " and " << map_fn << "..." << endl;
 		list_out << ped_fn << endl;

 		if(isThreadEnabled) {
 			JOB *job = new JOB();

 			for(int k = 0; k < ind_ids.size(); k++)
 				job->ind_ids.push_back(ind_ids[k]);
 			for(int k = 0; k < site_ids.size(); k++)
 				job->site_ids.push_back(site_ids[k]);

 			 job->ped_outfilename = ped_fn;
 			 job->map_outfilename = map_fn;
 			 job->pop_id = population_id;

 			 job_list.push_back(job);
 		} else {
 			// Create and export the file
 			export_data(vcf, vcf_info, ind_ids, site_ids, ped_fn, map_fn);
 		}

 		// Generate run command
 		run_out << "~/share/tools/plink/plink --ped " <<ped_fn << " --out " <<  ped_fn << " --map " <<  map_fn << "  --missing-phenotype 0 --no-parents --map3 --blocks --1 --noweb" << endl;

 		i = j - overlapping_size;
 		export_n++;
 	}

 	list_out.close();
 	run_out.close();

 	return export_n;
}



// int generate_split_map_pop(VCF *vcf, VCF_POPULATION_INFO *vcf_info, int map_size, int overlapping_size, string outfile_prefix)
// {
//	 int export_n = 0, n = 0;
//	 ofstream run_out;
//	 string run_fn = outfile_prefix;
//	 run_fn.append(".run");
//
//	 run_out.open(run_fn);
//
//	 for(int k = 0; k < vcf_info->population_ids.size(); k++) {
//		 string pop_id = vcf_info->population_ids[k];
//		 string pop_outfile_prefix= outfile_prefix;
//		 pop_outfile_prefix.append(".");
//		 pop_outfile_prefix.append(pop_id);
//
//		 n = generate_split_map(vcf, vcf_info, map_size, overlapping_size, pop_outfile_prefix,  pop_id);
//
//		 if(n > 0)
//			 run_out << "more " << pop_outfile_prefix<< ".run >> run.commands"  << endl;
//
//		 export_n += n;
//	 }
//
//	 run_out.close();
//
//	 cout << "Number of export: " << export_n << endl;
//
//	 return export_n;
// }




int generate_split_map_pop(VCF *vcf, VCF_POPULATION_INFO *vcf_info, int map_size, int overlapping_size, string outfile_prefix)
{
	int export_n = 0, n = 0;
	ofstream run_out;
	string run_fn = outfile_prefix;
	run_fn.append(".run");

	run_out.open(run_fn);


	for(int k = 0; k < vcf_info->population_ids.size(); k++) {
		string pop_id = vcf_info->population_ids[k];
		string pop_outfile_prefix= outfile_prefix;
		pop_outfile_prefix.append(".");
		pop_outfile_prefix.append(pop_id);

		n = generate_split_map(vcf, vcf_info, map_size, overlapping_size, pop_outfile_prefix,  pop_id);


		if(n > 0)
			run_out << "more " << pop_outfile_prefix<< ".run >> run.commands"  << endl;

		export_n += n;
	}

	run_out.close();

	cout << "Number of export: " << export_n << endl;

	return export_n;
}





 /**
  * Reduce map back
  */
 int reduce_map()
 {
 	return 0;
 }





//
//void generate_population()
//{
//
//	ofstream map_out, ped_out, run_out;
//
//	ped_out.open(ped_outfile.c_str());
//
//	for(int i = 0; i < vcf->individual_num(); i++) {
//		string ind_id = vcf->get_individual_id(i);
//		//ped_out << ind_id.c_str();
//
//		/*
//		if(i == 0)
//		cout << i << ") processing " << ind_id << "..." << endl;
//		 */
//		ind_info = vcf_info->get_individual_info(ind_id);
//
//		if(ind_info == 0)
//		{
//			cout << ind_id << " does not have a valid meta data available, default missing values will be put." << endl;
//			// Family ID\tIndividual ID\tGender Code (0: missing, 1: male, 2: female)\tPhenotype (0: unaffected, 1: affected, -9: unknown )
//			ped_out << ind_id << "\t" << ind_id << "\t0" <<  "\t0";
//		} else {
//			if(ind_info->fid != 0)
//				ped_out << ind_info->fid;
//			else
//				ped_out << ind_id;
//
//			ped_out << "\t" << ind_id;
//
//			ped_out << "\t" << ind_info->gender;
//
//			ped_out << "\t0";
//		}
//
//
//		for(int j = 0; j < vcf->site_num(); j++) {
//			string site_id = vcf->get_site_id(j);
//			string ref_a, alt_a;
//
//
//			/*
//			if(i == 0)
//				cout << site_id << endl;
//			 */
//
//			char g1, g2;
//			int geno = vcf->get_genotype(site_id, ind_id);
//
//			//ref_a = vcf.get_ref_allele(site_id);
//			//alt_a = vcf.get_ref_allele(site_id);
//
//			if(geno == '0') {
//				ped_out << "  1 1";
//				//g1 = '0'; g2 = '0';
//			} else if(geno == '1') {
//				ped_out << "  1 2";
//				//g1 = '0'; g2 = '1';
//			} else if(geno == '2') {
//				ped_out << "  2 1";
//				//g1 = '1'; g2 = '0';
//			} else if(geno == '3') {
//				ped_out << "  2 2";
//				//g1 = '1'; g2 = '1';
//			} else { // Missing
//				ped_out << "  0 0";
//				//g1 = '?'; g2 = '?';
//			}
//
//			//ped_out << g1 << " " << g2
//		}
//		ped_out << endl;
//		ped_out.flush();
//	}
//	ped_out.close();
//
//
//	cout << "Creating MAP file..." << endl;
//
//	map_out.open(map_outfile.c_str());
//	for(int i = 0; i < vcf->site_num(); i++) {
//		string site_id = vcf->get_site_id(i);
//
//		VCF_DATA* vcf_data = vcf->get_vcf(site_id);
//
//		if(vcf_data)
//			map_out <<  vcf_data->chrom << " "<< *vcf_data->id << " " << vcf_data->pos <<  endl;
//		else
//			cout << i << ") Problem with the site " << site_id << "." << endl;
//	}
//	map_out.close();
//
//
//	cout << "VCF size: " << vcf->site_num() << endl;
//
//
//
//
//
//	/**
//	 * Output ped file for each population
//	 */
//
//
//	run_out.open("run.pl");
//
//	vector<string> site_ids;
//	for(int i = 0; i < vcf->site_num(); i++) {
//		site_ids.push_back(vcf->get_site_id(i));
//	}
//
//
//	for(int k = 0; k < vcf_info->population_ids.size(); k++) {
//		VCF_POPULATION *pop = vcf_info->map_population.find(vcf_info->population_ids[k])->second;
//
//		cout << "Processing " << pop->population_id << "..." << endl;
//
//		int ind_n = 0;
//		vector<string> ind_ids;
//
//		for(int i = 0; i < pop->individual_infos.size(); i++)
//		{
//			string ind_id = pop->individual_infos[i]->individual_id;
//			if(vcf->is_individual_found(ind_id)) {
//				ind_ids.push_back(ind_id);
//				ind_n++;
//			}
//		}
//
//		if(ind_n > 0) {
//			cout << k << ": Creating PED file for the population " << pop->population_id << endl;
//
//			string pop_ped_fn = ped_outfile;
//			pop_ped_fn.append(".");
//			pop_ped_fn.append(pop->population_id);
//			pop_ped_fn.append(".ped");
//			string pop_map_fn = ped_outfile;
//			pop_map_fn.append(".");
//			pop_map_fn.append(pop->population_id);
//			pop_map_fn.append(".map");
//
//
//			run_out << "~/share/tools/plink/plink --ped " << pop_ped_fn << " --out " <<  pop_ped_fn << " --map " <<  pop_map_fn << "  --missing-phenotype 0 --no-parents --map3 --blocks --1 &" << endl;
//			int export_n = export_data(vcf, vcf_info, ind_ids, site_ids, pop_ped_fn, pop_map_fn);
//			cout << "Number of exported for " << pop->population_id << ": " << export_n << endl;
//		} else {
//			cout << k << ": No VCF data for the population " << pop->population_id << endl;
//		}
//
//		ind_ids.clear();
//	}
//
//}


 void print_usage()
 {
		cout << "Usage: 1kgp_to_plink THREAD_N VCF_DATA VCF_INFO 0|1" << endl;
		cout << "       1kgp_to_plink THREAD_N VCF_DATA VCF_INFO 0|1 MAP_SIZE OVERLAP_SIZE" << endl;
 }


// plink --file mydata --r2 --missing-phenotype 0 --no-fid --no-parents --no-sex --no-pheno --map3 --ld-window-r2
// plink --file mydata --map3 --no-parent --no-fid




// inside the vcfparser folder, run: svn co svn+ssh://143.89.24.58/home/bcnskaa/svn/vcfparser
// Compilation: g++ -o 1kgp_to_plink vcfparser/1kgp_to_plink.cpp vcfparser/vcf41_parser.cpp -O3  -std=c++0x -lm
 // g++ -o 1kgp_to_plink vcfparser/1kgp_to_plink.cpp vcfparser/vcf41_parser.cpp -O3 -I /home/bcnskaa/share/lib/boost-1.53.0/include -lboost_thread -L/home/bcnskaa/share/lib/boost-1.53.0/lib -std=c++0x -lm

// Plink to estimate haplotype block:  ~/share/tools/plink/plink --out ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.GABRB2.vcf.ped --ped ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.GABRB2.vcf.ped --map ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.GABRB2.vcf.map --missing-phenotype 0 --no-parents --map3 --blocks --1
// Plink to estimate haplotype freq:  ~/share/tools/plink/plink --out ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.GABRB2.vcf.ped --ped ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.GABRB2.vcf.ped --map ALL.chr5.phase1_release_v3.20101123.snps_indels_svs.genotypes.GABRB2.vcf.map --missing-phenotype 0 --no-parents --map3 --blocks --1

// option: --1 is added to inform plink that phenotype is coded in 0/1 fashion.
int main(int argc, char **argv)
{
	if(argc != 5 && argc != 7) {
		print_usage();
		return 0;
	}


	THREAD_N = atoi(argv[1]);
	string vcf_infile = argv[2];
	string vcf_info_infile = argv[3];
	int pop_export_flag = atoi(argv[4]);
	int map_size = atoi(argv[5]);
	int overlap_size = atoi(argv[6]);

	string map_outfile = vcf_infile;
	map_outfile.append(".map");
	string ped_outfile = vcf_infile;
	ped_outfile.append(".ped");

	//VCF_POPULATION_INFO *vcf_info;
	VCF_INDIVIDUAL_INFO* ind_inf;

	if(THREAD_N > 1 && THREAD_N <= 32)
		isThreadEnabled = true;
	else
		isThreadEnabled = false;


	cout << "Reading from " << vcf_info_infile << "..."<< endl;
	vcf_info = new VCF_POPULATION_INFO(vcf_info_infile);
	cout << vcf_info->individual_ids.size() << " individuals added." << endl;


	cout << "Reading from " << vcf_infile << "..."<< endl;
	vcf  = new VCF(vcf_infile);
	//VCF *vcf = new VCF(argv[1]);

	if(vcf == NULL || !vcf->isReady())
		return 0;



	if(argc == 7)
	{
		int map_n = 0;
		if(pop_export_flag == 1) { // Group individuals ethnically
			map_n = generate_split_map_pop(vcf, vcf_info,  map_size, overlap_size, vcf_infile);
			proccess();
			//cout << "Number of map generated: " << map_n << endl;
		} else if(pop_export_flag == 0) {  // No population divided
			map_n = generate_split_map(vcf, vcf_info,  map_size, overlap_size, vcf_infile, "");
			proccess();
		} else {
			print_usage();
		}
		cout << "Number of map generated: " << map_n << endl;
	} else {

	}



	clear_thread_resource();

	delete vcf;
	delete vcf_info;


	return 0;
}
