#ifndef _BAM24_HPP_
#define _BAM24_HPP_ 1

#include <utility>
#include <vector>

typedef std::pair< std::size_t, char> seqNuc_t;

typedef std::vector<std::vector<seqNuc_t>*> mappedReads_t;

typedef std::pair< std::size_t, char> entry;

mappedReads_t* bam24( char * filename);

#endif
