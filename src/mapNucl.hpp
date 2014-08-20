#ifndef _MAPNUCL_HPP_
#define _MAPNUCL_HPP_ 1

#include <utility>
#include <vector>

#include "nucl.hpp"

typedef std::vector<std::vector<Nucl>> mapped_nucl_t;

mapped_nucl_t mapNucl( const char * filename);

#endif
