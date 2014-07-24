#ifndef __MATHASH_H__
#define __MATHASH_H__
#include <list>
#include <map>
#include <unordered_map>
#include <array>
#include <utility>

/* uncomment for unordered hash map 

// A stupid hash function. Not sure if there might be a better idea.
struct matHashFn {
	std::size_t operator()( std::array<unsigned int, 24> key) const {
		std::size_t ret = 0;
		for( std::size_t i=0; i<24; i++){
			ret += key[i] << i;
		}
		return std::hash<size_t>()(ret);
	}
};

typedef std::unordered_map< std::array<unsigned int, 24>, float*, matHashFn> maptype;
 */

// This is the type for a ordered map.
typedef std::map< std::array<unsigned int, 24>, float*> maptype;

class matHash {
private:
	maptype* hmap;
public:
	matHash();
	void init( std::size_t MAX);

	void inc( unsigned int *key);

	void clear();

	maptype::iterator begin();
	maptype::iterator end();
};

#endif /* __MATHASH_H__ */
