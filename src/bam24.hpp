#ifndef _BAM24_HPP_
#define _BAM24_HPP_ 1

#include <utility>
#include <vector>
#include <forward_list>

#include "nucl.hpp"

typedef std::vector<std::vector<Nucl>> mappedReads_t;

mappedReads_t* bam24( char * filename);

#endif
