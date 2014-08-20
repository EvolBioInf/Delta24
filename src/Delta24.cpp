#include <vector>
#include <list>
#include <unordered_map>

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

#include "bam24.hpp"
#include "dml.hpp"

using namespace std;

// use count_t instead of uint, so we might be able to switch types in a future version
typedef int count_t;

/* uncomment for unordered hash map */

template <class T>
inline void hash_combine(std::size_t& seed, const T& v){
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// A stupid hash function. Not sure if there might be a better idea.
struct matHashFn {
	std::size_t operator()( std::array<count_t, 24> key) const {
		std::size_t ret = 0;
		for( std::size_t i=0; i<24; i++){
			hash_combine(ret, key[i]);
		}
		return ret;
	}
};

typedef unordered_map< array<count_t, 24>, count_t, matHashFn> map_t;

void inc( map_t& map, const array<count_t, 24>& key ){
	auto elem = map.find(key);
	if( elem != map.end()){
		elem->second++;
	} else {
		map.insert(map_t::value_type(key,1));
	}
}

/**
 * This function does some sorting and counting.
 */
map_t make_sorted_count ( size_t distance, const mappedReads_t::const_iterator begin, const mappedReads_t::const_iterator end){

	map_t matCounts{};

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
			countA[ it.getCode() ]++;
		}

		for( const auto& it: *J ){
			countB[ it.getCode() ]++;
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
			if( ii->getSeqID() < ji->getSeqID() ){
				count[sortA[ ii->getCode() ] ]++;
				++ii;
			} else if ( ii->getSeqID() > ji->getSeqID() ){
				count[sortB[ ji->getCode() ] + 4]++;
				++ji;
			} else {
				size_t offset = sortB[ ji->getCode() ] + sortA[ii->getCode()] * 4;
				count[8 + offset]++;
				++ii;
				++ji;
			}
		}

		for(; ii != ie; ++ii){
			count[ sortA[ ii->getCode() ]]++;
		}

		for(; ji != je; ++ji){
			count[ sortB[ ji->getCode() ] + 4]++;
		}

		inc( matCounts, count);
	}

	return matCounts;
}

map_t make_four_count ( const mappedReads_t::const_iterator begin, const mappedReads_t::const_iterator end ){
	
	map_t matCounts{};

	// Iterate over all elements.
	for (auto i = begin; i != end; ++i){
		// Deal with positions without mapped reads.
		array<count_t, 24> count;
		count.fill(0);

		if( i->empty() == true ) {
			inc( matCounts, count);
			continue;
		}

		// Count the number of mapped reads at this position.
		for( auto j: *i){
			count[ j.getCode() ]++;
		}

		inc( matCounts, count);
	}

	return matCounts;
}

void setcoef(float *coef, float pi, float eps, float delta){
	//These are a set of coeficents that appear numerous times in the likelihood calculations. They are computed here for efficency.
	coef[1]=log(0.333333333333333*eps*(-eps + 1));
	coef[2]=log(0.0555555555555556*pow(eps,2) + 0.166666666666667*eps*(-eps + 1));
	coef[3]=log(0.111111111111111*pow(eps,2) );
	coef[4]=log(pow((-eps + 1),2));
	coef[5]=log(0.166666666666667*eps*(-eps + 1) + 0.5*pow((-eps + 1),2));
	coef[6]=log(0.0555555555555556*pow(eps, 2) + 0.5*pow(1-eps, 2));
	coef[7]=log(-0.333333333333333*eps + 0.5);
	coef[8]=log(0.333333333333333*eps);
	coef[9]=log(-eps + 1);
	coef[10]=pow(1-pi,2)+delta*(1-pi)*pi;
	coef[11]=2*(1-delta)*(1-pi)*pi;
	coef[12]=pow(pi,2)+delta*(1-pi)*pi;
	coef[13]=(-2*pow(pi,2)+2*pi);
	coef[14]=(-pow(pi,2)+pi);
}

void compute( char* filename, size_t start, size_t stop ){
	float parms[4] = {0.01, 0.01, 0.0, 0.0};
	float coef[15] = {0};
	float R[2];
	float &pi = parms[0];
	float &eps = parms[1];

	cout << start << ", " << stop << endl;
	cout << "Starting main loop.\n";

	mappedReads_t foobar = bam24(filename);

	auto matCounts = make_four_count( foobar.begin(), foobar.end() );

	R[0] = 100;
	R[1] = 100;
	
	// some loop
	while ( (fabsf(R[0])+fabsf(R[1]) > 0.00001f )|| isnan(R[0]) || isnan(R[1]) ){
		setcoef(coef, pi, eps, 0.0);
		float iJ[2][2], J[2][2];

		J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
		R[0] = R[1] = 0.0;

		for( auto it : matCounts){
			float X[25];
			for(int i=0;i<24;i++){
				X[i] = static_cast<float>(it.first[i]);
			}
			int C = it.second;
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

		// R_ = R*iJ
		auto oldR0 = R[0];
		R[0] = oldR0 * iJ[0][0] + R[1] * iJ[0][1];
		R[1] = oldR0 * iJ[1][0] + R[1] * iJ[1][1];

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

	std::vector<double> delta( stop - start );

	#pragma omp parallel for
	for (size_t D = start; D < stop; D++){
		map_t matCounts = make_sorted_count( D, foobar.begin(), foobar.end());

		float dML_prev = 0;
		float dML_curr = 0;
		float &D_curr = parms[2];
		float D_prev = pi;

		D_curr = pi;

		setcoef( coef, pi, eps, D_curr);
		
		vector<dml_s> partial;

		for (auto i = matCounts.begin(); i != matCounts.end(); ++i){
			partial.push_back(dml_init(parms, i->first, i->second, coef ));
		}

		matCounts.clear();

		for( const auto &it: partial){
			dML_prev += dml_comp(it, coef);
		}

		D_curr = D_prev / 2;
		setcoef( coef, pi, eps, D_curr);

		for( auto it: partial){
			dML_curr += dml_comp(it, coef);
		}

		/**
		 * The next loop is a Newton iteration, or so I guess: From the (imaginary)
		 * function ML(D) we need the maximum. Standard school maths tells us to
		 * look at ML'(D)=0. I this code ML'(D)=\sum_X F0(X,D,…). 
		 *
		 * The process used here can be expressed as:
		 *
		 *  D_{n+1} = D_n - ML'(D) / ML''(D)
		 *
		 * Since ML'(D) is already pretty complicated and we don't want to have yet 
		 * another partial derivative ML''(D) is simply computed via the difference
		 * quotient (slope):
		 *
		 *  ML''(D_n) ≃ (ML'(D_n) - ML'(D_{n-1})) / (D_n - D_{n-1})
		 *
		 * Variable Names:
		 *  D_curr = D_n
		 *  D_prev = D_{n-1}
		 *  ML'(D_n) = dML_curr
		 *  ML'(D_{n-1}) = dML_prev
		 *
		 */

		int passes = 0;
		while (fabs(D_prev - D_curr) > 0.0000001 && passes <= 15){
			double slope = (dML_curr - dML_prev) / (D_curr - D_prev);
			double temp = D_curr - dML_curr / slope;

			dML_prev = dML_curr;
			D_prev = D_curr;
			D_curr = temp;

			dML_curr = 0;
			setcoef( coef, pi, eps, D_curr);

			for( auto it: partial){
				dML_curr += dml_comp(it, coef);
			}

			passes++;
		}

		if (passes > 15){
			//cerr << "Failure to converge\n";
			delta[D-start] = -42.0;
		} else {
			delta[D-start] = D_curr;
		}
	}

	for( uint i=0; i<( stop-start ); i++){
		if( delta[i] == -42.0){
			cout << "D=" << i+start << " Delta=" << "NC" << endl;
		} else {
			cout << "D=" << i+start << " Delta=" << delta[i] << endl;
		}
	}
}

int main (int argc, char**argv){

	if( argc < 4 ){
		printf("usage: %s <bam file> <start distance> <stop distance> \n", argv[0]);
		exit(0);
	}

	size_t start = atoi(argv[2]);
	size_t stop  = atoi(argv[3]);
	printf("bam file:%s start distance:%lu stop distance:%lu \n", argv[1], start, stop);

	compute(argv[1], start, stop);
}
