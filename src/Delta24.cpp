#include <vector>
#include <iostream>
#include <cstdio>

#include "compute.hpp"


int main (int argc, char**argv){

	if( argc < 4 ){
		printf("usage: %s <bam file> <start distance> <stop distance> [lumping]\n", argv[0]);
		exit(0);
	}

	size_t start = atoi(argv[2]);
	size_t stop  = atoi(argv[3]);
	size_t lumping = argc >= 5 ? atoi(argv[4]) : 1;
	printf("bam file:%s start distance:%lu stop distance:%lu lumping:%lu\n", argv[1], start, stop, lumping);

	compute(argv[1], start, stop, lumping);
}
