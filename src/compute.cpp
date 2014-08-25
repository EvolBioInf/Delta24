#include <list>
#include <unordered_map>

#include <iostream>


#include <assert.h>
#include <string>

#include <algorithm>
#include <numeric>

//The likelihood functions
#include "PyRoeNewton.h"

#include "mapNucl.hpp"
#include "dml.hpp"

using namespace std;

// use count_t instead of uint, so we might be able to switch types in a future version
typedef int count_t;
typedef array<count_t, 24> count_24mer_t;

template <class T>
inline void hash_combine( size_t& seed, const T& v){
    hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// A stupid hash function. Not sure if there might be a better idea.
struct matHashFn {
	size_t operator()( const array<count_t, 24>& key) const {
		size_t ret = 0;
		for( size_t i=0; i<24; i++){
			hash_combine(ret, key[i]);
		}
		return ret;
	}
};


typedef unordered_map< count_24mer_t, count_t, matHashFn> count_map_t;

/**
 * This function does some sorting and counting.
 */
void make_sorted_count ( count_map_t& countMap, size_t distance, const mapped_nucl_t::const_iterator begin, const mapped_nucl_t::const_iterator end){

	mapped_nucl_t::const_iterator I = begin;
	mapped_nucl_t::const_iterator J = begin;
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
			if( ii->getReadID() < ji->getReadID() ){
				count[sortA[ ii->getCode() ] ]++;
				++ii;
			} else if ( ii->getReadID() > ji->getReadID() ){
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

		countMap[count]++;
	}
}

/** @brief Count the ACGT-patterns at each postion.
 *
 * For the computation of pi and epsilon the full 24-mer count is not yet
 * needed. Instead a simple cumulative count of the ACGT-patterns is provided.
 *
 * @param begin - An iterator to the first mapped nucleotide.
 * @param end - An iterator to the last mapper nucleotide.
 * @returns A counting map of the patterns.
 */
count_map_t make_four_count ( const mapped_nucl_t::const_iterator& begin, const mapped_nucl_t::const_iterator& end ){
	
	count_map_t countMap{};

	// Iterate over all positions.
	for (auto i = begin; i != end; ++i){
		array<count_t, 24> pattern;
		pattern.fill(0);

		// Deal with positions without mapped reads.
		if( i->empty() == true ) {
			countMap[pattern]++;
			continue;
		}

		// Create the ACGT-Pattern from this position.
		for( const auto& j: *i){
			pattern[ j.getCode() ]++;
		}

		countMap[pattern]++;
	}

	return countMap;
}

/** @brief Cache some coefficients.
 *
 * These are a set of coefficents that appear numerous times in the likelihood 
 * calculations. They are computed here for efficency.
 *
 * @param coef - Output parameter of the coefficient array.
 * @param pi - The mutation rate.
 * @param eps - The epsilon parameter.
 * @paran delta - The delta parameter.
 */
void setcoef(double *coef, double pi, double eps, double delta){
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

/** @brief Compute pi and epsilon.
 *
 * @param begin - An iterator to the first mapped nucleotide.
 * @param end - An iterator to the last mapper nucleotide.
 * @returns A double pair: The first is pi, the second is epsilon.
 */
pair<double,double> compute_pi_eps ( const mapped_nucl_t::const_iterator& begin, const mapped_nucl_t::const_iterator& end){
	double parms[4] = {0.01, 0.01, 0.0, 0.0};
	double coef[15] = {0};
	double R[2] = { 100, 100};
	double &pi = parms[0];
	double &eps = parms[1];

	auto patternMap = make_four_count( begin, end);

	// Resolve ambiguity.
	using std::isnan;

	// some loop
	while ( isnan(R[0]) || isnan(R[1]) || fabs(R[0]) + fabs(R[1]) > 0.00001){
		setcoef(coef, pi, eps, 0.0);
		double iJ[2][2], J[2][2];

		J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
		R[0] = R[1] = 0.0;

		// iterate over all patterns
		for( const auto& it : patternMap){
			// convert to legacy interface
			double X[25] = {
				static_cast<double>(it.first[0]),
				static_cast<double>(it.first[1]),
				static_cast<double>(it.first[2]),
				static_cast<double>(it.first[3]),
				0.0
			};
			double C = static_cast<double>(it.second);
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

	// TODO: Find a better way to print information.
	cout << "Pi=" << pi << ", Epsilon=" << eps << endl;

	return {pi,eps};
}

void compute( const char* filename, size_t start, size_t stop, size_t lumping ){

	mapped_nucl_t mappedNucls = mapNucl(filename);

	auto pieps = compute_pi_eps( mappedNucls.cbegin(), mappedNucls.cend());
	double pi = pieps.first;
	double eps = pieps.second;

	double coef[15] = {0};

	size_t length = (stop-start)/lumping + 1;
	std::vector<double> delta( length );

	#pragma omp parallel for
	for (size_t i = 0; i < length; ++i){
		size_t D = i * lumping + start;
		count_map_t countMap {};

		for(uint l=0; l<lumping; l++){
			make_sorted_count( countMap, D+l, mappedNucls.begin(), mappedNucls.end());
		}

		double parms[4] = {pi, eps, 0.01, 0.0};
		double dML_prev = 0;
		double dML_curr = 0;
		double &D_curr = parms[2];
		double D_prev = pi;

		D_curr = pi;

		setcoef( coef, pi, eps, D_curr);
		
		vector<dml_s> partial;

		for (auto i = countMap.begin(); i != countMap.end(); ++i){
			partial.push_back(dml_init(parms, i->first, i->second, coef ));
		}

		countMap.clear();

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

		if (passes > 15 || D_curr < -1.0 || D_curr > 1.0){
			//cerr << "Failure to converge\n";
			delta.at(i) = -42.0;
		} else {
			delta.at(i) = D_curr;
		}
	}

	for( uint i=0; i< length; i++){
		if( delta.at(i) == -42.0){
			cout << "D=" << i*lumping+start << " Delta=" << "NC" << endl;
		} else {
			cout << "D=" << i*lumping+start << " Delta=" << delta[i] << endl;
		}
	}
}