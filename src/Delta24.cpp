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

int names_allocated=0;

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

template <class Type>
class matFile{
	private:

	Type *block;		//the array which stores read and call information.
	unsigned int sites;	//total number of sites covered.
	unsigned int allocated;	//total number of 'Type' allocated. 
	char tempname[100];	//name of temporary file (used for automatic saves if memory runs low). 
	Type *end_;		//the first unallocated address after the array 
	Type *ptr;		//the pointer to the address 'currently' being read.

	public:
	float coverage;

	matFile (){
		/*basic constructor, just's sets everything to zeros.*/
		sprintf(tempname, "tempfile_%d.bin", names_allocated);
		names_allocated++;
				sites=0;                                               
				allocated=0;                                               
		block=NULL;
		end_=NULL;
		ptr=NULL;
	};

	void write(){
		ofstream tempfile(tempname, ios::out | ios::binary);
		tempfile.write((char *)&sites, sizeof(unsigned int) );
		tempfile.write((char *)&allocated, sizeof(unsigned int) );
		tempfile.write((char *)block, sizeof(Type)*allocated);
		tempfile.close();
		if (block!=NULL) delete block;
		block=NULL;
		end_=NULL;
		ptr=NULL;
		sites=0;
		allocated=0;
	};

	void close(){
		if (block!=NULL) delete block;
		block=NULL;
		end_=NULL;
		ptr=NULL;
		sites=0;
		allocated=0;
	};

	void read(){
		ifstream tempfile(tempname, ios::in | ios::binary);
		tempfile.read((char*)&sites, sizeof(unsigned int));
		tempfile.read((char*)&allocated, sizeof(unsigned int));
		if (block!=NULL) delete block;
		block=NULL;
		block=new Type [allocated];
		ptr=block;
		end_=block+allocated-1;
		tempfile.read((char*)block, sizeof(Type)*allocated);
		tempfile.close();
	};

	matFile (vector <Type> *a, unsigned int size_){

		sprintf(tempname, "tempfile_%d.bin", names_allocated);
		names_allocated++;

		/*This is the constructor. vector *a is the structure in which data about base calls is temporarily storred so that the matFile is not constantly being updated. 
		a matFile only stores information from a single contig/scaffold. size_ 

		The basic strucutre of the file is a long array of 'Type's arranged in the following way:

		[N0][C0][C1]..[CN0][N1][C0][C1]..[CN1]

		where N0 is the number of 'C' values that follow the value N0, N1 is the number of 'C' values that follow N1 and so on.
		each 'C' is a combined read number and nucleotide call. The two lowest value bits specify the nucleotide read, and the remainder of the bits specify the read.
	
		The total amount of memory used by this structure should be one 'Type' for each nucleotide called within a read, and one 'Type' for each site in the reference genome.
		
		*/
		unsigned int total_size=0; 						//
		for (int x=0; x<size_; x++) total_size+=a[x].size();			//
		sites=size_;
		allocated=sites+total_size;					//
		block=new Type [allocated];

		Type *site=block;							//
		for (int x=0; x<sites; x++){						//
			*site=a[x].size();
			site++;
			for (int y=0; y<a[x].size(); y++){ *site=a[x][y]; site++; }
		};
		end_=site;								//
	};
	void add_reads(vector <Type> *a, unsigned int new_sites){
		/*This function alocates a new array larger array, copies over the old values and then adds the new values.
		Currently the old array data is copied over the slowest way possible, but I am too lazy to fix it. */

		unsigned int total_size=allocated;

		for (int x=0; x<new_sites; x++){
			total_size+=a[x].size();
		};

		Type *new_block = new Type [total_size+new_sites-sites];
		allocated = total_size+(new_sites-sites);

		//TODO I should check to make sure this alocation worked.

		Type *new_site=new_block;
		Type *old_site=block;

		unsigned int old_size;
		for (int x=0; x<sites; x++){
			old_size=*old_site;
			*new_site=*old_site+a[x].size();
			new_site++;
			old_site++;
			for (int y=0; y<old_size; y++){
				*new_site=*old_site;
				new_site++;
				old_site++;
			}; 
			for (int y=0; y<a[x].size(); y++){*new_site=a[x][y]; new_site++;};
			a[x].clear();
		};

		for (int x=sites; x<new_sites; x++){
			*new_site=a[x].size();
			new_site++;
			for (int y=0; y<a[x].size(); y++){ *new_site=a[x][y]; new_site++;}
			a[x].clear();
		};
		
		// Look, I even delete the old block so we don't have a memory leak! Wow, what skill! 
		if  (block!=NULL) delete block;

		block=new_block;
		end_=new_site;
		sites=new_sites;
		coverage=float(allocated-sites)/float(sites);

	};

	Type *operator[](int idx) {
		//The basic operator A[X] returns a pointer to the Xth site in the file (i.e NX, not CX).  
		ptr=block;
		for (int x=0; x<idx; x++){
			ptr+=(*ptr+1);
		};
		return ptr;
	};

	class iterator {
		//The iterator for a matFile.  
		private:
		friend class matFile;

		public: 

		//The very basic operators, nothing special here
			Type* ptr_;
		iterator(Type * other){	ptr_=other;};
		iterator(){};	
		Type* operator->() { return ptr_; }
		const Type* operator->() const { return ptr_; }
		Type& operator*() { return *ptr_; }
		const Type& operator*() const { return *ptr_; }
		bool operator==(const iterator & other) const {return ptr_ == other.ptr_; }
		bool operator!=(const iterator & other) const {return ptr_ != other.ptr_; }

		iterator& operator=(iterator other){
			ptr_=other.ptr_;
			return *this;
		}

		//prefix incrementor. 
		iterator & operator++(){
			ptr_=ptr_+(*ptr_+1);
			return *this;
		};

		//I really can't remember why I thought I needed this....
		int to_int (){
			return int(ptr_);
		};
	
		//This returns the adress of C0 and N1 in the array [N0][C0]...[CN0][N1][
		Type *inner_begin() { return ptr_+1; }
		Type *inner_end() { return (ptr_+(*ptr_)+1); }
		Type inner_size() { return *ptr_; }
	};
	//This returns the adress of the begining and end of block.
	iterator begin() { return iterator(block); }
	iterator end() { return iterator(end_); }
};

//reads a SamRecord into a character array. The buffer filled in this function is used by a switch/case statement in the function "read" to construct a matFile. 

vector<char>* getCalls( SamRecord &record, size_t size ){
	std::vector<char> *v = new vector<char>(size);
	int get = 0;
	for( size_t i=0; i< size; i++){
		get = record.getCigarInfo()->getQueryIndex(i);
		(*v)[i] = (get >= 0 ? record.getSequence(i) : '*');
	}
	return v;
}

//returns the sum of all the different reads at a specific site, and allows us to calculate read depth.
template< typename T>
T sum (T* X){
	T ret = (T)0;
	for( size_t i= 0; i< 24; i++ ){
		ret += X[i];
	}
	return ret;
}

//A global. Can easily just be moved into main()
matHash matCounts;

int uint_cmp(const void *a, const void *b){
	const unsigned int *ia = (const unsigned int *)a;
	const unsigned int *ib = (const unsigned int *)b;
	return *ia - *ib;
}

int uint_cmp_desc(const void *a, const void *b){
	return - uint_cmp(a,b);
}

matHash* make_sorted_count ( size_t distance, vector<vector<entry>*>* projection){
	matHash *matCounts = new matHash();

	size_t length = projection->size();
	for (size_t i = 0; i < length - distance; ++i){
		size_t j = i + distance;

		if( (*projection)[i] == NULL || (*projection)[j] == NULL){
			continue;
		}

		count_t count[24] = {0};
		count_t countA[4] = {0}, countB[4] = {0};

		array<count_t, 4> sortA, sortB;

		for( size_t k=0; k<4; k++){
			sortA[k] = sortB[k] = k;
		}

		// FIXME: The following line is fugly. Fix that type!
		for( auto& it: *((*projection)[i]) ){
			countA[char2uint(it.second)]++;
		}

		for( auto& it: *((*projection)[j]) ){
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

		auto ii = ((*projection)[i])->begin();
		auto ie = ((*projection)[i])->end();

		auto ji = ((*projection)[j])->begin();
		auto je = ((*projection)[j])->end();

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

void make_sorted_count (int dist, matFile <unsigned int> *file, int LS, unsigned int MIN, unsigned int MAX){
	/*this function totals up the number of occrance of A,C,G and T at site A and at site B, sorts them, and then makes the full 24 count*/

	matFile<unsigned int>::iterator A;	//A pointer to site A.
	matFile<unsigned int>::iterator B;	//A pointer to site B.
	matFile<unsigned int>::iterator end;	//A pointer to the end of the matFile.
	unsigned int *a, *b, *aend, *bend;	//Pointers to calls at sites A and B.
	unsigned int mask=3;			//The mask to extract the call (i.e A,C,G or T) from the read number. 

	int X=0;				//An int to keep track of the number of sites read for debuging purposes. 
	long unsigned int All=0;		//An int to store the number of reads covering a site. Probably won't show up in a final version.

	unsigned int count[24], countA[4]={0,0,0,0}, countB[4]={0,0,0,0}, sortA[4]={0,1,2,3}, sortB[4]={0,1,2,3};
	for (int f=0; f<LS; f++){

		//thing thing thing
		file[f].read();
		A=file[f].begin();
		B=file[f][dist];
		end=file[f].end();
		while (B!=end){
			/*obtain seperate mono-nucleotide counts for sites A nd B .*/

			//Clear ALL the DATA!
			memset(count, 0, sizeof(int)*24);
			memset(countA, 0, sizeof(int)*4);
			memset(countB, 0, sizeof(int)*4);

			sortA[0]=0; sortA[1]=1; sortA[2]=2; sortA[3]=3;
			sortB[0]=0; sortB[1]=1; sortB[2]=2; sortB[3]=3;

			a=A.inner_begin();
			b=B.inner_begin();
			aend=A.inner_end();
			bend=B.inner_end();


			for (; a!=aend; a++){
				countA[ ((*a)&mask)] += 1;
			}
			for (; b!=bend; b++){
				countB[ ((*b)&mask)] += 1;
			}


			/*An explicit sorting network to make 0 largest and 3 smallest count*/
			{
				if (countA[0]<countA[2]){
					swap(countA[0], countA[2]);
					swap(sortA[0], sortA[2]);
				}
				if(countA[1]<countA[3]){
					swap(countA[1], countA[3]);
					swap(sortA[1], sortA[3]);
				}
				if(countA[2]<countA[3]){
					swap(countA[2], countA[3]);
					swap(sortA[2], sortA[3]);
				}
				if(countA[0]<countA[3]){
					swap(countA[0], countA[3]);
					swap(sortA[0], sortA[3]);
				}
				if(countA[1]<countA[2]){
					swap(countA[1], countA[2]);
					swap(sortA[1], sortA[2]);
				}


				if (countB[0]<countB[2]){
					swap(countB[0], countB[2]);
					swap(sortB[0], sortB[2]);
				}
				if(countB[1]<countB[3]){
					swap(countB[1], countB[3]);
					swap(sortB[1], sortB[3]);
				}
				if(countB[2]<countB[3]){
					swap(countB[2], countB[3]);
					swap(sortB[2], sortB[3]);
				}
				if(countB[0]<countB[3]){
					swap(countB[0], countB[3]);
					swap(sortB[0], sortB[3]);
				}
				if(countB[1]<countB[2]){
					swap(countB[1], countB[2]);
					swap(sortB[1], sortB[2]);
				}
			}
			/*end sort. Get DI-nucleotide counts sorted. */
			a=A.inner_begin();
			b=B.inner_begin();
			while(a!=aend && b!=bend){
				if ( ((*a)>>2) < ((*b)>>2)){
					count[sortA[((*a)&mask)]] += 1;
					a++;
				}
				else if ( ((*a)>>2)>((*b)>>2)){
					count[sortB[((*b)&mask)]+4] += 1;
					b++;
				}
				else {
					count[8+sortB[((*b)&mask)]+sortA[((*a)&mask)]*4] += 1;
					a++;
					b++;
				}
			}
			while (a!=aend){
				count[ sortA[((*a)&mask)]] += 1;
				a++;
			}
			while (b!=bend){
				count[ sortB[((*b)&mask)]+4] += 1;
				b++;
			}

			All=sum(count);
			if (All>MIN && All<MAX) matCounts.inc(count);
			++A;
			++B;
			X++;
		}

		file[f].close();
	}
}

float make_four_count ( matHash matCounts, vector<vector<entry>*>* map){
	unsigned int count[24];
	for (auto i = map->begin(); i != map->end(); ++i){
		memset(count, 0, sizeof(unsigned int)*24);

		if( !*i) {
			matCounts.inc(count);
			continue;
		}

		for( auto j = (*i)->begin(); j != (*i)->end(); j++){
			size_t offset = 0;
			switch( j->second){
				case 'A': offset = 0; break;
				case 'C': offset = 1; break;
				case 'G': offset = 2; break;
				case 'T': offset = 3; break;
			}
			count[ offset ]++;
		}

		matCounts.inc(count);
	}

	return 30.0;
}

matFile <unsigned int> *read (char *filename, int LS){
	//This is the basic function to set up an array of matFiles from a bam file. The argument LS determines the number of scaffolds to be used from the bam file.

	//Junk to set up and read the bam file.
	BamIndex bamIndex;

	SamFile samIn;
	SamFileHeader samHeader;
	SamRecord samRecord, next_samRecord;

	if(!samIn.OpenForRead(filename))
	{ 
		std::cout << "Failure to open " << filename << std::endl; 
	}
	if(!samIn.ReadHeader(samHeader)){
		std::cout << "Failure to read header on " << filename << std::endl; 
	};

	SamHeaderRecord *headerRecord;

	vector <int> scaffolds;
	vector <int> test (1,10);

	//Checks the headers to figure out the length of each scaffold.

	//AVALIBLE MEMORY
	// size_t tLS;

	for (int x=0; x<LS; x++){
		headerRecord=samHeader.getNextSQRecord();
		scaffolds.push_back(atoi(headerRecord->getTagValue("LN") ) );
		//?
	}
	size_t slice_start=0, slice_stop=0;

	int ReadIndex=0;
	int start0, end0, start1, end1;
	int S0, S1;

	matFile <unsigned int> *myFile;
	myFile=new matFile <unsigned int> [LS];

	for (int y=0; y<LS; y++) myFile[y].write();
	for (int y=0; y<LS; y++) myFile[y].close();

	struct sysinfo info;
	unsigned int vector_allocated=0;

	sysinfo(&info);	

	unsigned long memory_cap=(info.freeram+info.bufferram-(long int)(1073741824)/info.mem_unit*2)/sizeof(unsigned int)*info.mem_unit;

	cout << "memory_cap: " << memory_cap << endl;

	while (slice_stop<LS){
		size_t tLN=0;
		slice_start=slice_stop;
		while (tLN<memory_cap) {
			tLN+=scaffolds[slice_stop];
			slice_stop++;
			if (slice_stop==LS) break;
		};
		cout << "Reading slice " << slice_start << ", " << slice_stop << endl;
		vector <unsigned int> **very_large_array=new vector <unsigned int> *[scaffolds.size()];
	
		for (int y=slice_start; y<slice_stop; y++){
			very_large_array[y]=new vector <unsigned int> [scaffolds[y]];
		}

		while(samIn.ReadRecord(samHeader, samRecord) ){
			ReadIndex+=1;
			S0=samRecord.getReferenceID();

			if (S0<slice_stop && S0>=slice_start) {
				start0=samRecord.get0BasedPosition();
				end0=samRecord.get0BasedAlignmentEnd();
				int length = samRecord.getReadLength();

				/*
				if(SamFlag::isProperPair(samRecord.getFlag() ) ) {
					cout << "NUNUNU" << endl;
					samIn.ReadRecord(samHeader, next_samRecord);
					S1=next_samRecord.getReferenceID();

					if (S1<slice_stop && S1>=slice_start) {
						auto v = getCalls(next_samRecord, end1-start1);
						start1 = next_samRecord.get0BasedPosition();
						end1 = next_samRecord.get0BasedAlignmentEnd();

						int x = start1;
						for( auto it: *v){
							size_t offset = 0;
							switch(it){
								case 'A': offset = 0; break;
								case 'C': offset = 1; break;
								case 'G': offset = 2; break;
								case 'T': offset = 3; break;
							}

							very_large_array[S1][x].push_back( (ReadIndex << 2) + offset );
							x++;
						}

						delete v;
					}
				} */

				auto v = getCalls( samRecord, length);

				for( int x = 0; x < v->size(); x++){
					vector_allocated++;
					size_t offset = 0;
					switch( (*v)[x]){
						case 'A': offset = 0; break;
						case 'C': offset = 1; break;
						case 'G': offset = 2; break;
						case 'T': offset = 3; break;
					}
					very_large_array[S0][x + start0].push_back( (ReadIndex << 2) + offset );
				}

				delete v;
			}
		}
		//The slice is read, let's check to see if we should save the vectors (to free up some memory), then set up the matFiles...
		cout << "Finished reading slice...\n";
		for (int y=0; y<scaffolds.size(); y++){
			myFile[y].read();
			myFile[y].add_reads(very_large_array[y], scaffolds[y]);
			cout << "Scaffold " << y << " : " << myFile[y].coverage << endl;
			delete [] very_large_array[y];
			myFile[y].write();
			myFile[y].close();
		}	
		//and let's delete the buffer that we won't need anymore.
		delete [] very_large_array;
	};
	return myFile;
};

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
};

void
printCount( matHash matCounts ){
	ofstream derp("matCounts.tmp", ios::out);
		for( auto it : matCounts){
			for( size_t i = 0; i<25; i++){
				derp << it.second[i] << " ";
			}
			derp << endl;
		}
	derp.close();
}

int main (int argc, char**argv){

	//options : Flat (i.e. no paired end info)
	//

	if(argc!=8){
		printf("usage: %s <bam file> <number of scaffolds> <start distance> <stop distance> <increment> <min> <max>\n", argv[0]);
		exit(0);
	};

	int LS=atoi(argv[2]); 
	int start=atoi(argv[3]);
	int stop=atoi(argv[4]);
	int inc=atoi(argv[5]);
	int min=atoi(argv[6]);
	int max=atoi(argv[7]);

	printf("bam file:%s number of scaffolds:%d start distance:%d stop distance:%d increment:%d min:%d max:%d\n", argv[1], LS, start, stop, inc, min, max);

	cout << min << ", " << max << endl;


	float D_0, D_1, lnL_0, lnL_1, T;
	float parms[4]={0.01, 0.01, 0.0, 0.0};
	float coef[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float iJ[2][2], J[2][2], R[2];

	matFile <unsigned int> *data=read(argv[1], LS);
	maptype::iterator it, end;

	float *X, C;

	cout << start << ", " << stop << endl;
	cout << "Starting main loop.\n";

	R[0]=100;
	R[1]=100;

	matCounts = matHash();

	auto foobar = bam24(argv[1]);

	make_four_count( matCounts, foobar );

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
		//matCounts.clear();
		//make_sorted_count(D, data, LS, min, max);
		auto matCounts = *make_sorted_count ( D, foobar);
		printCount(matCounts);

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
