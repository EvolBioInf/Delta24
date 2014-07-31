#ifndef _BAM24_HPP_
#define _BAM24_HPP_ 1

#include <utility>
#include <vector>

typedef std::pair< std::size_t, char> entry;

std::vector<std::vector<entry>*>* bam24( char * filename);

#endif
