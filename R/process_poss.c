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

/**
 * 1. read fasta sequence
 * 2.
 *
 * Compilation: R CMD SHLIB process_poss.c
 */
SEXP process_poss_r(SEXP s, SEXP e)
{
	SEXP poss, Rdim;
	int sposs_n, eposs_n;
	int *sposs, *eposs;
	int spos, epos;
	int i, j, n;
	int cur_pos;

	sposs = INTEGER(s);
	eposs = INTEGER(e);

	//Rdim = getAttrib(s, R_DimSymbol);
	//sposs_n = INTEGER(Rdim)[0];
	sposs_n = length(s);
	eposs_n = length(e);
	n = 0;


	if(eposs_n != sposs_n)
	{
		Rprintf("sposs and eposs have different size.");
	}


	Rprintf("Number of regions: %d\n", sposs_n);

	for(i = 0; i < sposs_n; i++)
	{
		//spos = INTEGER(s)[i];
		//epos = INTEGER(e)[i];
		spos = sposs[i];
		epos = eposs[i];
		n += epos - spos;
	}

	PROTECT(poss = allocVector(INTSXP, n));

	// Reset the array into zeros
	for(i = 0; i < n; i++)
		INTEGER(poss)[i] = 0;

	//cur_pos = INTEGER(poss);
	cur_pos = 0;
	for(i = 0; i < sposs_n; i++)
	{
		//INTEGER(poss)[i] = i;
		//spos = INTEGER(s)[i];
		//epos = INTEGER(e)[i];
		spos = sposs[i];
		epos = eposs[i];

		//Rprintf("Processing: %d-%d\n", spos, epos);
		for(j = spos; j < epos; j++, cur_pos++)
		{
			INTEGER(poss)[cur_pos] = j;
		}
	}

	UNPROTECT(1);

	return poss;
}

//
//SEXP process_poss_r(SEXP s, SEXP e)
//{
//	SEXP poss, Rdim;
//	int sposs_n, eposs_n;
//	int *sposs, *eposs;
//	int spos, epos;
//	int i, j, n;
//	int cur_pos;
//
//	sposs = INTEGER(s);
//	eposs = INTEGER(e);
//
//	//Rdim = getAttrib(s, R_DimSymbol);
//	//sposs_n = INTEGER(Rdim)[0];
//	sposs_n = length(s);
//	eposs_n = length(e);
//	n = 0;
//
//
//	if(eposs_n != sposs_n)
//	{
//		Rprintf("sposs and eposs have different size.");
//	}
//
//
//	Rprintf("Number of regions: %d\n", sposs_n);
//
//	for(i = 0; i < sposs_n; i++)
//	{
//		//spos = INTEGER(s)[i];
//		//epos = INTEGER(e)[i];
//		spos = sposs[i];
//		epos = eposs[i];
//		n += epos - spos;
//	}
//
//	PROTECT(poss = allocVector(INTSXP, n));
//
//	// Reset the array into zeros
//	for(i = 0; i < n; i++)
//		INTEGER(poss)[i] = 0;
//
//	//cur_pos = INTEGER(poss);
//	cur_pos = 0;
//	for(i = 0; i < sposs_n; i++)
//	{
//		//INTEGER(poss)[i] = i;
//		//spos = INTEGER(s)[i];
//		//epos = INTEGER(e)[i];
//		spos = sposs[i];
//		epos = eposs[i];
//
//		//Rprintf("Processing: %d-%d\n", spos, epos);
//		for(j = spos; j < epos; j++, cur_pos++)
//		{
//			INTEGER(poss)[cur_pos] = j;
//		}
//	}
//
//	UNPROTECT(1);
//
//	return poss;
//}
//




//SEXP process_poss_vals_r(SEXP s, SEXP e, SEXP v, int max_pos)
SEXP process_poss_vals_r(SEXP s, SEXP e, SEXP v, SEXP len)
{
	SEXP sexp_vals, Rdim;
	int sposs_n, eposs_n, vals_n, max_len;
	int *sposs, *eposs;
	double *vals;
	int spos, epos;
	int i, j, n;
	int cur_pos;
	double val;

	sposs = INTEGER(s);
	eposs = INTEGER(e);
	vals = REAL(v);
	max_len = INTEGER(len)[0];

	//Rdim = getAttrib(s, R_DimSymbol);
	//sposs_n = INTEGER(Rdim)[0];
	sposs_n = length(s);
	eposs_n = length(e);
	vals_n = length(v);
	n = 0;


	if(eposs_n != sposs_n && vals_n != sposs_n)
	{
		Rprintf("sposs and eposs, or poss and values have different size.");
		return(0);
	}

	// Check if the max
	/*
	if(eposs_n != sposs_n && vals_n != sposs_n)
	{
			Rprintf("sposs and eposs, or poss and values have different size.");
			return(0);
	}
	*/


	Rprintf("Number of regions: %d\n", sposs_n);
	Rprintf("Max length: %d\n", max_len);

	/*
	for(i = 0; i < sposs_n; i++)
	{
		//spos = INTEGER(s)[i];
		//epos = INTEGER(e)[i];
		spos = sposs[i];
		epos = eposs[i];
		n += epos - spos;
	}
	*/


	PROTECT(sexp_vals = allocVector(REALSXP, max_len));

	//Rprintf("list intialized.\n");

	// initialize vector
	for(i = 0; i < max_len; i++)
	{
		REAL(sexp_vals)[i] = 0.0;
	}

	//cur_pos = INTEGER(poss);
	//cur_pos = 0;
	for(i = 0; i < sposs_n; i++)
	{
		//REAL(sexp_vals)[i] = i;
		//spos = INTEGER(s)[i];
		//epos = INTEGER(e)[i];
		spos = sposs[i] - 1;
		epos = eposs[i] - 1;
		val = vals[i];

		//Rprintf("Spos=%d, Epos=%d, Val=%f.\n", spos, epos, val);

		//Rprintf("Processing: %d-%d\n", spos, epos);
		for(j = spos; j < epos && j < max_len; j++)
		{
			REAL(sexp_vals)[j] = val;
		}
	}

//
//	for(i = 0; i < max_len; i+=1000) {
//		Rprintf("%d: %f\n", i, REAL(sexp_vals)[i]);
//	}

	UNPROTECT(1);

	return sexp_vals;
}




double* process_poss_vals_c(int *s, int *e, double *v, int n, int len)
{
	int i, j, spos, epos, processed_n;
	double *vals, cur_val;

//Rprintf("Len: %d\n", len);

	processed_n = 0;

	vals = (double*)calloc(len, sizeof(double));
	for(i = 0; i < len; i++)
		vals[i] = 0.0;


	for(i = 0; i < n; i++)
	{
		spos = s[i] - 1;
		epos = e[i] - 1;
		cur_val = v[i];

//Rprintf("spos=%d epos=%d\n", spos, epos);

		for(j = spos; j < epos; j++)
		{
			if(vals[j] < cur_val)
				vals[j] = cur_val;
			processed_n++;
		}
	}

	Rprintf("Processed=%d\n", processed_n);

	return (vals);
}





//SEXP process_poss_vals_r(SEXP s, SEXP e, SEXP v, int max_pos)
SEXP process_poss_vals_to_ranges_r(SEXP s, SEXP e, SEXP v, SEXP rs, SEXP re, SEXP len)
{
	SEXP sexp_vals, sexp_rvals, Rdim;
	int sposs_n, eposs_n, vals_n, max_len, rsposs_n, reposs_n;
	int *sposs, *eposs, *rsposs, *reposs;
	double *vals, *poss;
	int spos, epos, rspos, repos;
	int i, j, n;
	int cur_pos;
	double val;


	//Rdim = getAttrib(s, R_DimSymbol);
	//sposs_n = INTEGER(Rdim)[0];
	sposs = INTEGER(s);
	eposs = INTEGER(e);
	sposs_n = length(s);
	eposs_n = length(e);

	// Range definition
	rsposs = INTEGER(rs);
	reposs = INTEGER(re);
	rsposs_n = length(rs);
	reposs_n = length(re);

	vals = REAL(v);
	vals_n = length(v);
	n = 0;

	max_len = INTEGER(len)[0];

	//
	if(eposs_n != sposs_n && vals_n != sposs_n)
	{
		Rprintf("sposs and eposs, or poss and values have different size.");
		return(0);
	}

	Rprintf("Number of range object: %d\n", rsposs_n);
	Rprintf("Number of data regions: %d\n", sposs_n);
	Rprintf("Max length: %d\n", max_len);


	poss = process_poss_vals_c(sposs, eposs, vals, sposs_n, max_len);


//Rprintf("Passed\n", max_len);


	// Convert the range or value object into poss objects
//	sexp_vals = process_poss_vals_r(s, e, v, len);
//	if(sexp_vals == 0)
//	{
//		Rprintf("Problem on constructing a poss list.");
//		return(0);
//	}
//
//	// Convert the SEXP object into double
//	vals = REAL(sexp_vals);

	// Prepare memory for storing range values
	PROTECT(sexp_rvals = allocVector(REALSXP, rsposs_n));


	// Reset the memory
	for(i = 0; i < rsposs_n; i++)
	{
		REAL(sexp_rvals)[i] = 0.0;

		rspos = rsposs[i] - 1;
		repos = reposs[i] - 1;

//Rprintf("rspos=%d repos=%d\n", rspos, repos);

		for(j = rspos; j < repos; j++)
		{
//			if(j >= 58662131 && j < 58662148)
//			{
//				if(poss[j] > 0)
//				{
//					Rprintf("%d: %f, %f\n", j, poss[j], REAL(sexp_rvals)[i] );
//
//				}
//			}
			REAL(sexp_rvals)[i] += poss[j];
		}
	}

	free(poss); poss = 0;

	UNPROTECT(1);

	return sexp_rvals;
}



/*
int* process_poss(int n, int *sposs, int *eposs)
{
	int i;

	for()
	{

	}
	for(i in 1:nrow(all_rmsk)) {
		r <- all_rmsk[i, ];
		if(exists("my_poss")) {
			my_poss <- c(my_poss, process_rmsk_2_poss(r));
		} else {
			my_poss <- process_rmsk_2_poss(r);
		}
	}
	my_poss <- unique(my_poss);
	return(my_poss);
}
*/
