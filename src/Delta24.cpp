#include <vector>
#include <list>

#include <stdio.h>
#include <stdlib.h>

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

//used to initialize random number generater. Probably won't exist in a final program.
#include <time.h>

using namespace std;

int names_allocated=0;

class matHash{
	private:
	float **tree;
	float *fptr;
	list <float *> values;
	list <float *> keys; // keys is never read!
	bool alloc;
	unsigned int MAX;

	public:
	matHash(){
	};

	void init(unsigned int _MAX){
		MAX = _MAX;
		tree= (float**) malloc(MAX * sizeof(float **));
		for(size_t i=0; i<MAX; i++){
			tree[i] = NULL;
		}
	};

	void inc(unsigned int *key){
		float **tree_ptr = tree;
		unsigned int *key_ptr = key;

		bool alloc = false;
		for (int x=0; x< 24; x++){
			if((*key_ptr) >= MAX){
				cout << "Whoops, you lied to me!\n";
			}

			if( tree_ptr[*key_ptr] == NULL){
				tree_ptr+=*key_ptr;
				*tree_ptr=(float *) calloc (MAX, sizeof(float **)); //This is probably a horible way to do this...
				keys.push_back( (*tree_ptr) );
				tree_ptr=(float **)(*tree_ptr);
				alloc=true;
			}
			else {
				tree_ptr=(float **)(tree_ptr[*key_ptr]);
			}
			key_ptr++;
		}
		if (alloc){
			float *fptr = new float[25];
			values.push_back(fptr);
			*tree_ptr = fptr;
			key_ptr = key;

			for (int x=0; x<24; x++){
				fptr[x] = float(key_ptr[x]);
			}

			fptr[24] = 1.0;
		}
		else{
			(*tree_ptr)[24]++;
		};
	}

	void clear(void){
		cout << "Clearing\n";

		for( auto it: values){
			free(it);
		}
		values.clear();

		for( auto it: keys){
			free(it);
		}
		keys.clear();
		cout << "Done\n";

		for( size_t i=0; i< MAX; i++){
			tree[i] = NULL;
		}
	};

	list <float *>::iterator begin() { return values.begin(); }
	list <float *>::iterator end() { return values.end(); }	
};

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
		try {block=new Type [allocated];}
		catch (std::bad_alloc& ba){std::cerr << "insufficent memory to read in buffer : " << ba.what() << '\n';}
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
		try {block=new Type [allocated];}
		catch (std::bad_alloc& ba){std::cerr << "insufficent memory to read in buffer : " << ba.what() << '\n';}
				Type *site=block;							//
				for (int x=0; x<sites; x++){						//
					*site=a[x].size();
						site++;
					for (int y=0; y<a[x].size(); y++){ *site=a[x][y]; site++; };
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

		Type *new_block;
		unsigned long MB=1048576;
		try {new_block=new Type [total_size+new_sites-sites];}
		catch (std::bad_alloc& ba){std::cerr << "insufficent memory to create file " << tempname << ". Requested " << (total_size+new_sites)*sizeof(Type)/MB << " MB.\n";}
		allocated=total_size+(new_sites-sites);

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
int getCalls(SamRecord &this_record, char *buffer, int size){

		int get=0;
	if (size>200) size=200;
		for (int x=0;x<size; x++){
				get=this_record.getCigarInfo()->getQueryIndex(x);
		if (get>0){
					buffer[x]=this_record.getSequence(get);
		}
		else{
			buffer[x]='*';
		}
		}
		return 1;
};

//returns the sum of all the different reads at a specific site, and allows us to calculate read depth.
float sum(float *X){
	float ret = 0.0;
	for( size_t i= 0; i< 24; i++ ){
		ret += X[i];
	}
	return ret;
}

//let's overload it a bit.
unsigned int sum(unsigned int *X){
	auto ret = 0;
	for( size_t i= 0; i< 24; i++ ){
		ret += X[i];
	}
	return ret;
}

//A global. Can easily just be moved into main()
map <string, float *> counts;
matHash matCounts;

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

			while(a!=aend && b!=bend){
					countA[((*a)&mask)]+=1;
					a++;
					countB[((*b)&mask)]+=1;
					b++;
			};
			while (a!=aend){
					countA[ ((*a)&mask)]+=1;
					a++;
			};
			while (b!=bend){
					countB[ ((*b)&mask)]+=1;
					b++;
			}
			a=A.inner_begin();
			b=B.inner_begin();

			/*An explicit sorting network to make 0 largest and 3 smallest count*/
		
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
			
			/*end sort. Get DI-nucleotide counts sorted. */
		
			while(a!=aend && b!=bend){
					if ( ((*a)>>2) < ((*b)>>2)){
							count[sortA[((*a)&mask)]]+=1;
							a++;
					}
					else if ( ((*a)>>2)>((*b)>>2)){
							count[sortB[((*b)&mask)]+4]+=1;
							b++;
					}
					else {
							count[8+sortB[((*b)&mask)]+sortA[((*a)&mask)]*4]+=1;
							a++;
							b++;
					};
			};
			while (a!=aend){
					count[ sortA[((*a)&mask)]]+=1;
					a++;
			};
			while (b!=bend){
					count[ sortB[((*b)&mask)]+4]+=1;
					b++;
			}

			char key[200];
			All=sum(count);
			if (All>MIN && All<MAX) matCounts.inc(count);
			++A;
			++B;
			X++;
			};
		file[f].close();
	};
//        delete counts["0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"];
//        counts.erase("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0");
};


float make_four_count (matFile <unsigned int> *file, int LS, unsigned int MIN, unsigned int MAX){
	//returns an estimate of average read depth, ignoring coverage 0 sites.

		matFile<unsigned int>::iterator A;
		matFile<unsigned int>::iterator end;
		unsigned int *a, *aend;
		unsigned int count[24];
		unsigned int mask=3;
		int X=0;
		long unsigned int Di=0, All=0;
	
	for (int f=0; f<LS; f++){
		float depth=0, sites=0;
		file[f].read();
			A=file[f].begin();
			end=file[f].end();
			while (A!=end){
			//cout << int(end.ptr_)-int(A.ptr_) << endl;
					memset(count, 0, sizeof(unsigned int)*24);
					a=A.inner_begin();
					aend=A.inner_end();
					while (a!=aend){
							count[ ((*a)&mask)]+=1;
							a++;
					};

					All=(count[0]+count[1]+count[2]+count[3]);
					if (All>MIN && All<MAX) matCounts.inc(count);
					++A;
					X++;
			};
		file[f].close();
	}
		//delete counts["0,0,0,0"];
		//counts.erase("0,0,0,0");
	return 30;
};


matFile <unsigned int> *read (char *filename, int LS){
	//This is the basic function to set up an array of matFiles from a bam file. The argument LS determines the number of scaffolds to be used from the bam file.

	//Junk to set up and read the bam file.
	uint64_t start, stop;
	BamIndex bamIndex;
	SamStatus::Status returnStatus = SamStatus::SUCCESS;

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
	size_t usage=test.capacity()*sizeof(int) + sizeof(test);

	//Checks the headers to figure out the length of each scaffold.

	//AVALIBLE MEMORY
	size_t tLS;

	for (int x=0; x<LS; x++){
		headerRecord=samHeader.getNextSQRecord();
		scaffolds.push_back(atoi(headerRecord->getTagValue("LN") ) );
		//?
	}
	size_t slice_start=0, slice_stop=0;

	int ReadIndex=0;
	int start0, end0, start1, end1;
	int S0, S1, N0, N1, M, Match, Tries;
	char buffer[200];

	matFile <unsigned int> *myFile;
	myFile=new matFile <unsigned int> [LS];

	for (int y=0; y<LS; y++) myFile[y].write();
	for (int y=0; y<LS; y++) myFile[y].close();

	struct sysinfo info;
	unsigned int vector_allocated=0;

	sysinfo(&info);	

	unsigned long memory_cap=(info.freeram+info.bufferram-(long int)(1073741824)/info.mem_unit*2)/sizeof(unsigned int)*info.mem_unit;

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
			if (vector_allocated>100000000){
				cout << "ding!\n";
							for (int y=0; y<scaffolds.size(); y++){
								myFile[y].read();
									myFile[y].add_reads(very_large_array[y], scaffolds[y]);
									myFile[y].write();
									myFile[y].close();
							}
				vector_allocated=0;
						}

			if (S0<slice_stop && S0>=slice_start) {
				start0=samRecord.get0BasedPosition();
				end0=samRecord.get0BasedAlignmentEnd();
				if(SamFlag::isProperPair(samRecord.getFlag() ) ) {
					//cout << "found mate pair\n";
						samIn.ReadRecord(samHeader, next_samRecord);
					S1=next_samRecord.getReferenceID();
							//if(not(next_samRecord.getCigarInfo()->hasIndel() ) ){
					if (S1<slice_stop && S1>=slice_start) {
						getCalls(next_samRecord, buffer, end1-start1);
						start1=next_samRecord.get0BasedPosition();
						end1=next_samRecord.get0BasedAlignmentEnd();
						if (end1-start1<200){
							for(int x=start1; x<end1; x++) {
									vector_allocated++;
									switch (buffer[x-start1]){
										case 'A':
											very_large_array[S1][x].push_back((ReadIndex<<2));
										break;
										case 'C':
											very_large_array[S1][x].push_back((ReadIndex<<2)+1);
										break;
										case 'G':
											very_large_array[S1][x].push_back((ReadIndex<<2)+2);
										break;
										case 'T':
											very_large_array[S1][x].push_back((ReadIndex<<2)+3);
										break;
										default:
											break;
									};
							}
						}
					}
					//}
				};
				getCalls(samRecord, buffer, end0-start0);
				if (end0-start0<200){
					for(int x=start0; x<end0; x++){
							vector_allocated++;
							switch (buffer[x-start0]){
								case 'A':
									very_large_array[S0][x].push_back((ReadIndex<<2));
								break;
								case 'C':
									very_large_array[S0][x].push_back((ReadIndex<<2)+1);
								break;
								case 'G':
									very_large_array[S0][x].push_back((ReadIndex<<2)+2);
									break;
								case 'T':
									very_large_array[S0][x].push_back((ReadIndex<<2)+3);
								break;
								default:
									break;
							}

					}
				}
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

void parse (float *datum, const char *read){
	//The data generated by the make_?_count functions is stored in a dictionary that pairs a string representation of the observation
	//with a count of the number of times an observation has been seen. This function  takes the string representation of the observation
	//and turns it back into an array of floats. 
	const char *cptr=read;
	char buffer[200];
	memset(buffer, 0, sizeof(char)*200); 
	int x=0, y=0;
	while (*cptr!=0){
		if (*cptr==','){
			datum[x]=atoi(buffer);
			memset(buffer, 0, sizeof(char)*200); 
			x++;
			y=0;
		}
		else{
			buffer[y]=*cptr;
			y++;
		}
		cptr++;
	}
		datum[x]=atoi(buffer);
};

float fRand(){
	//Returns a random float between 0 and 1.
	//I don't think I use this currently.
	double f = (double)rand() / RAND_MAX;
	return (float)(f);
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
};

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


	float H[2]={0,0}, D_0, D_1, lnL_0, lnL_1, T;
	float parms[4]={0.01, 0.01, 0.0, 0.0};
	float coef[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float lnL=0, iJ[2][2], J[2][2], R[2];

	matFile <unsigned int> *data=read(argv[1], LS);
	list <float*>::iterator it, end;

	float *X, C;

	cout << start << ", " << stop << endl;
	cout << "Starting main loop.\n";

	R[0]=100;
	R[1]=100;

	matCounts.init(max);
	parms[3]=make_four_count(data, LS, min, max);
	
	while ((fabs(R[0])+fabs(R[1])>0.00001 )|| isnan(R[0]) || isnan(R[1]) ){
		setcoef(coef, parms);
				memset(J[0], 0, sizeof(float)*2);
				memset(J[1], 0, sizeof(float)*2);
				memset(R, 0, sizeof(float)*2);
				it=matCounts.begin();
				end=matCounts.end();
				while (it!=end ){
			X=*it;
			C=X[24];
						J[0][0]+=J00(parms, X, coef)*C;
						J[0][1]+=J01(parms, X, coef)*C;
						J[1][0]+=J10(parms, X, coef)*C;
						J[1][1]+=J11(parms, X, coef)*C;
						R[0]+=(R0(parms, X, coef) )*C;
						R[1]+=(R1(parms, X, coef) )*C;
						++it;
				};

		iJ[0][0]=1/(J[0][0]*J[1][1]-J[0][1]*J[1][0])*J[1][1];
		iJ[0][1]=-1/(J[0][0]*J[1][1]-J[0][1]*J[1][0])*J[0][1];
		iJ[1][0]=-1/(J[0][0]*J[1][1]-J[0][1]*J[1][0])*J[1][0];
		iJ[1][1]=1/(J[0][0]*J[1][1]-J[0][1]*J[1][0])*J[0][0];

		R[0]=(R[0]*iJ[0][0]+R[1]*iJ[0][1]);
		R[1]=(R[0]*iJ[1][0]+R[1]*iJ[1][1]);

		if (parms[0]>R[0]) parms[0]-=R[0];
		else parms[0]/=2.0;
		if (parms[1]>R[1]) parms[1]-=R[1];
		else parms[1]/=2.0;
		};
	cout << "Pi=" << parms[0] << ", Epsilon=" << parms[1]  << ", R=" << fabs(R[0])+fabs(R[1]) << endl;

	D_0=parms[0];
	D_1=parms[0];

	for (int D=start; D<stop; D+=inc){
		matCounts.clear();
		make_sorted_count(D, data, LS, min, max);
		
		lnL_0=0;
		if ( D_1<1 && D_1>0 ) parms[2]=D_1;
		else {parms[2]=parms[0]; D_0=parms[0];}

				it=matCounts.begin();
				end=matCounts.end();
		setcoef(coef, parms);

				while (it!=end ){
			X=*it;
			C=X[24];
					lnL_0+=F0(parms, X, coef)*C;
					++it;
				};
				lnL_1=0;
				D_1=D_0/2;
				parms[2]=D_1;
				it=matCounts.begin();
		setcoef(coef, parms);
				while (it!=end ){
			X=*it;
			C=X[24];
			lnL_1+=F0(parms, X, coef)*C;
			++it;
				};
		int inc=0;
		while (fabs(D_0-D_1)>0.0000001 ){
			if (inc>15){
				cerr << "Failure to converge\n";
				break;
			}
			inc++;
			T=parms[2]-lnL_1*(D_1-D_0)/(lnL_1-lnL_0);
			lnL_0=lnL_1;
			D_0=D_1;
			D_1=T;
			parms[2]=T;
			lnL_1=0;
			it=matCounts.begin();
			setcoef(coef, parms);

			while (it!=end ){
				X=*it;
				C=X[24];
				lnL_1+=F0(parms, X, coef)*C;
				++it;
			};
		}
		cout << "D=" << D << ", Delta=" << D_1 << ", Pi=" << parms[0] << ", Epsilon=" << parms[1] << endl;
	};
};
