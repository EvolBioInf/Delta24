#include <vector>
#include <iostream>
#include <cstdio>

#include <getopt.h>

#include "compute.hpp"

using namespace std;

void version();
void usage();

/** @brief The main function.
 *
 * Parses all commandline arguments and starts the computation with the correct
 * parameters.
 *
 * @returns 0 iff successful
 */
int main (int argc, char**argv){

	static struct option long_options[] = {
		{"version", no_argument, NULL, 'v'},
		{"help", no_argument, NULL, 'h'},
		{"min", required_argument, NULL, 'm'},
		{"max", required_argument, NULL, 'M'},
		{"lumping", required_argument, NULL, 'l'},
		{NULL,0,NULL,0}
	};

	size_t start = 1;
	size_t stop = 5;
	size_t lumping = 1;
	const char *filename = nullptr;

	while( 1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "hm:M:l:v", long_options, &option_index);

		if( c == -1){
			break;
		}

		switch(c){
			case 0: break;
			case 'h': usage(); break;
			case 'v': version(); break;
			case 'm':
				start = stoull(optarg);
				break;

			case 'M':
				stop = stoull(optarg);
				break;

			case 'l':
				lumping = stoull(optarg);
				break;

			case '?': // intentional fallthrough
			default:
				usage();
				break;
		}
	}

	if( stop < start){
		cerr << "The argument to --max should be greater than the one to --min. Aborting." << endl;
		exit(1);
	}

	if( optind == argc ){
		cerr << "Missing filename. Aborting." << endl;
		exit(1);
	} else {
		filename = argv[optind];
	}

	compute( filename, start, stop, lumping);
}

/** @brief Prints the version, copyright and license. Does not return. */
void version(){
	const char str[] = {
		"delta24 v0.1\n"
	};

	cout << str;
	exit(0);
}

/** @brief Prints usage instructions. Does not return. */
void usage(){
	const char str[] = {
		"usage: delta24 [OPTIONS] FILE\n"
		"\tFILE can be either a BAM or SAM file.\n"
		"Available Options:\n"
		"  -m, --min <num>     The smallest distance to check (inclusive); default: 1\n"
		"  -M, --max <num>     The biggest distance to check (inclusive); default: 5\n"
		"  -l, --lumping <num> Lumping; default: 1\n"
		"  -h, --help          Print this help.\n"
		"  -v, --version       Print the version information.\n"
	};

	cout << str;
	exit(0);
}
