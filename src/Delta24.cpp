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

typedef unsigned int count_t;

inline static count_t char2uint( const char c){
	count_t ret = 0;
	switch( c){
		case 'A': ret = 0; break;
		case 'C': ret = 1; break;
		case 'G': ret = 2; break;
		case 'T': ret = 3; break;
	}
	return ret;
}

matHash* make_sorted_count ( size_t distance, const mappedReads_t::const_iterator begin, const mappedReads_t::const_iterator end){

	matHash *matCounts = new matHash();

	mappedReads_t::const_iterator I = begin;
	mappedReads_t::const_iterator J = begin;
	advance(J, distance);

	//size_t length = projection->size();
	//for (size_t i = 0; i < length - distance; ++i){
	for(; J < end; I++, J++){
		if( !*I || !*J || !(*I)->size() || !(*J)->size() ){
			continue;
		}

		count_t count[24] = {0};
		count_t countA[4] = {0}, countB[4] = {0};

		array<count_t, 4> sortA, sortB;

		for( size_t k=0; k<4; k++){
			sortA[k] = sortB[k] = k;
		}

		// FIXME: The following line is fugly. Fix that type!
		for( auto it: **I ){
			countA[char2uint(it.second)]++;
		}

		for( auto it: **J ){
			countB[char2uint(it.second)]++;
		}

		auto cmpA = [&countA]( count_t a, count_t b){
			return countA[a] > countA[b];
		};
		sort( sortA.begin(), sortA.end(), cmpA);
		sort( sortB.begin(), sortB.end(), [&countB](count_t a, count_t b){
			return countB[a] > countB[b];
		});

		assert(countA[sortA[1]] >= countA[sortA[2]]);

		auto ii = (*I)->begin();
		auto ie = (*I)->end();

		auto ji = (*J)->begin();
		auto je = (*J)->end();

		while( ii != ie && ji != je ){
			if( ii->first < ji->first ){
				count[sortA[ char2uint(ii->second) ] ]++;
				ii++;
			} else if ( ii->first > ji->first ){
				count[sortB[ char2uint(ji->second)] + 4]++;
				ji++;
			} else {
				size_t offset = sortB[char2uint(ji->second)] + sortA[char2uint(ii->second)] * 4;
				count[8 + offset]++;
				ii++;
				ji++;
			}
		}

		for(; ii != ie; ii++){
			count[ sortA[char2uint(ii->second)]]++;
		}

		for(; ji != je; ji++){
			count[ sortB[char2uint(ji->second)] + 4]++;
		}

		matCounts->inc(count);
	}

	return matCounts;
}

float make_four_count ( matHash matCounts, const mappedReads_t::const_iterator begin, const mappedReads_t::const_iterator end ){
	unsigned int count[24];
	for (auto i = begin; i != end; ++i){
		memset(count, 0, sizeof(unsigned int)*24);

		if( !*i || !(*i)->size()) {
			matCounts.inc(count);
			continue;
		}

		for( auto j = (*i)->begin(); j != (*i)->end(); j++){
			count[ char2uint(j->second) ]++;
		}

		matCounts.inc(count);
	}

	return 30.0;
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

int main (int argc, char**argv){

	if(argc!=5){
		printf("usage: %s <bam file> <start distance> <stop distance> <increment>\n", argv[0]);
		exit(0);
	}

	int start=atoi(argv[2]);
	int stop=atoi(argv[3]);
	int inc=atoi(argv[4]);

	printf("bam file:%s start distance:%d stop distance:%d increment:%d\n", argv[1], start, stop, inc);


	float D_0, D_1, lnL_0, lnL_1, T;
	float parms[4]={0.01, 0.01, 0.0, 0.0};
	float coef[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float iJ[2][2], J[2][2], R[2];

	float *X, C;

	cout << start << ", " << stop << endl;
	cout << "Starting main loop.\n";

	R[0]=100;
	R[1]=100;

	auto matCounts = matHash();
	auto foobar = bam24(argv[1]);

	make_four_count( matCounts, foobar->begin(), foobar->end() );

	parms[3] = 30.0;
	
	while ((fabs(R[0])+fabs(R[1])>0.00001 )|| isnan(R[0]) || isnan(R[1]) ){
		setcoef(coef, parms);

		J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
		R[0] = R[1] = 0.0;

		for( auto it : matCounts){
			X = it.second;
			C = X[24];
			J[0][0] += J00(parms, X, coef) * C;
			J[0][1] += J01(parms, X, coef) * C;
			J[1][0] += J10(parms, X, coef) * C;
			J[1][1] += J11(parms, X, coef) * C;
			R[0] += (R0(parms, X, coef) ) * C;
			R[1] += (R1(parms, X, coef) ) * C;
		}

		auto detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

		iJ[0][0] = 1/detJ * J[1][1];
		iJ[0][1] =-1/detJ * J[0][1];
		iJ[1][0] =-1/detJ * J[1][0];
		iJ[1][1] = 1/detJ * J[0][0];

		R[0] = R[0] * iJ[0][0] + R[1] * iJ[0][1];
		R[1] = R[0] * iJ[1][0] + R[1] * iJ[1][1];

		if( parms[0] > R[0]) {
			parms[0]-= R[0];
		} else {
			parms[0]/= 2.0;
		}

		if( parms[1] > R[1]) {
			parms[1]-= R[1];
		} else {
			parms[1]/= 2.0;
		}
	}

	cout << "Pi=" << parms[0] << ", Epsilon=" << parms[1]  << ", R=" << fabs(R[0])+fabs(R[1]) << endl;

	D_0 = parms[0];
	D_1 = parms[0];

	for (int D = start; D < stop; D += inc){
		auto matCounts = *make_sorted_count ( D, foobar->begin(), foobar->end());

		lnL_0 = 0;

		if ( D_1 < 1 && D_1 > 0 ) {
			parms[2] = D_1;
		} else {
			parms[2] = parms[0];
			D_0 = parms[0];
		}

		setcoef(coef, parms);

		for( auto it: matCounts){
			X = it.second;
			C = X[24];
			lnL_0 += F0(parms, X, coef) * C;
		}

		lnL_1=0;
		D_1=D_0/2;
		parms[2]=D_1;
		setcoef(coef, parms);
		for( auto it: matCounts){
			X = it.second;
			C = X[24];
			lnL_1 += F0(parms, X, coef) * C;
		}

		int passes = 0;
		while (fabs(D_0 - D_1) > 0.0000001 && passes <= 15){
			T = parms[2]-lnL_1*(D_1-D_0)/(lnL_1-lnL_0);
			lnL_0 = lnL_1;
			D_0 = D_1;
			D_1 = T;
			parms[2] = T;
			lnL_1 = 0;
			setcoef(coef, parms);

			for( auto it: matCounts){
				X = it.second;
				C = X[24];
				lnL_1 += F0(parms, X, coef) * C;
			}

			passes++;
		}

		if (passes > 15){
			cerr << "Failure to converge\n";
		}

		cout << "D=" << D << ", Delta=" << D_1 << ", Pi=" << parms[0] << ", Epsilon=" << parms[1] << endl;
		matCounts.clear();
	}
}
