#ifndef _BAM24_HPP_
#define _BAM24_HPP_ 1

#include <utility>
#include <vector>
#include <list>

typedef std::pair< std::size_t, char> seqNuc_t;

typedef std::vector<std::list<seqNuc_t>> mappedReads_t;

mappedReads_t* bam24( char * filename);

#endif
