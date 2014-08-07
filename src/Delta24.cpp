#include <vector>
#include <list>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

#include <cmath>
#include <math.h>
#include <algorithm>
#include <numeric>

//needed to read bam files.
#include "SamFile.h"
#include "SamFlag.h"
#include "SamValidation.h"

//The likelihood functions
#include "PyRoeNewton.h"

//memory managment routines
#include <unistd.h>
#include <sys/sysinfo.h> 

#include "matHash.hpp"
#include "bam24.hpp"

using namespace std;

// use count_t instead of uint, so we might be able to switch types in a future version
typedef unsigned int count_t;

/**
 * @brief Maps nucleotides to a two bit code.
 */
inline static count_t char2uint( const char& c){
	count_t ret = 0;
	switch( c){
		case 'A': ret = 0; break;
		case 'C': ret = 1; break;
		case 'G': ret = 2; break;
		case 'T': ret = 3; break;
	}
	return ret;
}

/**
 * This function does some sorting and counting.
 */
matHash* make_sorted_count ( size_t distance, const mappedReads_t::const_iterator begin, const mappedReads_t::const_iterator end){

	matHash *matCounts = new matHash();

	mappedReads_t::const_iterator I = begin;
	mappedReads_t::const_iterator J = begin;
	advance(J, distance);

	// Iterate over all positions, with `I` and `J` being `distance` positions apart.
	for(; J < end; I++, J++){
		if( I->empty() == true || J->empty() == true ){
			continue;
		}

		array<count_t, 24> count;
		count.fill(0);
		count_t countA[4] = {0}, countB[4] = {0};

		array<count_t, 4> sortA, sortB;
		iota( sortA.begin(), sortA.end(), 0);
		iota( sortB.begin(), sortB.end(), 0);

		for( const auto& it: *I ){
			countA[char2uint(it.second)]++;
		}

		for( const auto& it: *J ){
			countB[char2uint(it.second)]++;
		}

		/* For a reason I dont know yet, we need to sort the nucleotides by frequency.
		 * The sortA/B arrays represent the sorted sequence in a fast lookup table.
		 * Here the sorting is implemented using `std::sort` and the new and shiny
		 * lambdas.
		 * Sort is descending!
		 */
		auto cmpA = [&countA]( count_t a, count_t b){
			return countA[a] > countA[b];
		};
		sort( sortA.begin(), sortA.end(), cmpA);
		sort( sortB.begin(), sortB.end(), [&countB](count_t a, count_t b){
			return countB[a] > countB[b];
		});

		auto ii = I->cbegin();
		auto ie = I->cend();

		auto ji = J->cbegin();
		auto je = J->cend();

		// I dont know, what these next lines are for.
		while( ii != ie && ji != je ){
			if( ii->first < ji->first ){
				count[sortA[ char2uint(ii->second) ] ]++;
				++ii;
			} else if ( ii->first > ji->first ){
				count[sortB[ char2uint(ji->second)] + 4]++;
				++ji;
			} else {
				size_t offset = sortB[char2uint(ji->second)] + sortA[char2uint(ii->second)] * 4;
				count[8 + offset]++;
				++ii;
				++ji;
			}
		}

		for(; ii != ie; ++ii){
			count[ sortA[char2uint(ii->second)]]++;
		}

		for(; ji != je; ++ji){
			count[ sortB[char2uint(ji->second)] + 4]++;
		}

		matCounts->inc(count);
	}

	return matCounts;
}

void make_four_count ( matHash matCounts, const mappedReads_t::const_iterator begin, const mappedReads_t::const_iterator end ){
	
	// Iterate over all elements.
	for (auto i = begin; i != end; ++i){
		// Deal with positions without mapped reads.
		array<count_t, 24> count;
		count.fill(0);

		if( i->empty() == true ) {
			matCounts.inc(count);
			continue;
		}

		// Count the number of mapped reads at this position.
		for( auto j = i->begin(); j != i->end(); ++j){
			count[ char2uint(j->second) ]++;
		}

		matCounts.inc(count);
	}
}

void setcoef(float *coef, float *parms){
	//These are a set of coeficents that appear numerous times in the likelihood calculations. They are computed here for efficency.
	coef[1]=log(0.333333333333333*parms[1]*(-parms[1] + 1));
	coef[2]=log(0.0555555555555556*pow(parms[1],2) + 0.166666666666667*parms[1]*(-parms[1] + 1));
	coef[3]=log(0.111111111111111*pow(parms[1],2) );
	coef[4]=log(pow((-parms[1] + 1),2));
	coef[5]=log(0.166666666666667*parms[1]*(-parms[1] + 1) + 0.5*pow((-parms[1] + 1),2));
	coef[6]=log(0.0555555555555556*pow(parms[1], 2) + 0.5*pow(1-parms[1], 2));
	coef[7]=log(-0.333333333333333*parms[1] + 0.5);
	coef[8]=log(0.333333333333333*parms[1]);
	coef[9]=log(-parms[1] + 1);
	coef[10]=pow(1-parms[0],2)+parms[2]*(1-parms[0])*parms[0];
	coef[11]=2*(1-parms[2])*(1-parms[0])*parms[0];
	coef[12]=pow(parms[0],2)+parms[2]*(1-parms[0])*parms[0];
	coef[13]=(-2*pow(parms[0],2)+2*parms[0]);
	coef[14]=(-pow(parms[0],2)+parms[0]);
}

void compute( char* filename, size_t start, size_t stop, size_t inc ){
	float parms[4] = {0.01, 0.01, 0.0, 0.0};
	float coef[15] = {0};
	float R[2];
	float &pi = parms[0];
	float &eps = parms[1];

	cout << start << ", " << stop << endl;
	cout << "Starting main loop.\n";

	auto matCounts = matHash();
	auto foobar = bam24(filename);

	make_four_count( matCounts, foobar->begin(), foobar->end() );

	// FIXME: Hm, 30? Why not 42?
	parms[3] = 30.0;
	R[0] = 100;
	R[1] = 100;
	
	// some loop
	while ( (fabsf(R[0])+fabsf(R[1]) > 0.00001f )|| isnan(R[0]) || isnan(R[1]) ){
		setcoef(coef, parms);
		float iJ[2][2], J[2][2];

		J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
		R[0] = R[1] = 0.0;

		for( auto it : matCounts){
			float *X = it.second;
			float C = X[24];
			J[0][0] += J00(parms, X, coef) * C;
			J[0][1] += J01(parms, X, coef) * C;
			J[1][0] += J10(parms, X, coef) * C;
			J[1][1] += J11(parms, X, coef) * C;
			R[0] += R0(parms, X, coef) * C;
			R[1] += R1(parms, X, coef) * C;
		}

		// invert J
		auto detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

		iJ[0][0] = 1/detJ * J[1][1];
		iJ[0][1] =-1/detJ * J[0][1];
		iJ[1][0] =-1/detJ * J[1][0];
		iJ[1][1] = 1/detJ * J[0][0];

		// FIXME: Is it just me, or do the following lines look like a bug?
		R[0] = R[0] * iJ[0][0] + R[1] * iJ[0][1];
		R[1] = R[0] * iJ[1][0] + R[1] * iJ[1][1];

		if( pi > R[0]) {
			pi-= R[0];
		} else {
			pi/= 2.0;
		}

		if( eps > R[1]) {
			eps-= R[1];
		} else {
			eps/= 2.0;
		}
	}

	cout << "Pi=" << pi << ", Epsilon=" << eps  << ", R=" << fabs(R[0])+fabs(R[1]) << endl;
	cout << "R(" << R[0] << "," << R[1] << ")" << endl;

	matCounts.clear();

	#pragma omp parallel for
	for (size_t D = start; D < stop; D += inc){
		auto matCounts = *make_sorted_count ( D, foobar->begin(), foobar->end());

		float lnL_0 = 0;
		float lnL_1 = 0;
		float D_1 = pi;
		float D_0 = pi;

		setcoef(coef, parms);
		D_1 = D_0 / 2;
		parms[2] = D_1;

		for( const auto& it: matCounts){
			float *X = it.second;
			lnL_0 += F0(parms, X, coef) * X[24];
		}

		setcoef(coef, parms);
		for( auto it: matCounts){
			float *X = it.second;
			lnL_1 += F0(parms, X, coef) * X[24];
		}

		int passes = 0;
		while (fabs(D_0 - D_1) > 0.0000001 && passes <= 15){
			float temp = parms[2] - lnL_1 * (D_1 - D_0) / (lnL_1 - lnL_0);

			lnL_0 = lnL_1;
			D_0 = D_1;
			D_1 = parms[2] = temp;

			lnL_1 = 0;
			setcoef(coef, parms);

			for( auto it: matCounts){
				float *X = it.second;
				lnL_1 += F0(parms, X, coef) * X[24];
			}

			passes++;
		}

		if (passes > 15){
			cerr << "Failure to converge\n";
		}

		cout << "D=" << D << ", Delta=" << D_1 << ", Pi=" << pi << ", Epsilon=" << eps << endl;
		matCounts.clear();
	}
}

int main (int argc, char**argv){

	if( argc != 5 ){
		printf("usage: %s <bam file> <start distance> <stop distance> <increment>\n", argv[0]);
		exit(0);
	}

	size_t start = atoi(argv[2]);
	size_t stop  = atoi(argv[3]);
	size_t inc   = atoi(argv[4]);

	printf("bam file:%s start distance:%lu stop distance:%lu increment:%lu\n", argv[1], start, stop, inc);

	compute(argv[1], start, stop, inc);	
}
