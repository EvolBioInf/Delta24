#ifndef _NUCL_HPP_
#define _NUCL_HPP_

#include <vector>

using std::size_t;

class Nucl {
private:
	size_t data;
public:
	Nucl(size_t seqID, char n);
	Nucl();
	static size_t char2code( char c);
	static char code2char( size_t d);
	bool operator< (const Nucl&) const;
	bool operator> (const Nucl&) const;

	size_t getCode() const{
		return data & 0x3UL;
	}

	size_t getSeqID() const {
		return data & (~0UL ^ 0x3UL);
	}
};

#endif // _NUCL_HPP_
