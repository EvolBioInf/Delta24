#ifndef _NUCL_HPP_
#define _NUCL_HPP_

#include <vector>

class Nucl {
private:
	std::size_t data;
public:
	Nucl(std::size_t seqID, char n);
	Nucl();
	static std::size_t char2code( char c);
	static char code2char( std::size_t d);
	std::size_t getCode() const;
	std::size_t getSeqID() const;
	bool operator< (const Nucl&) const;
	bool operator> (const Nucl&) const;
};

#endif // _NUCL_HPP_
