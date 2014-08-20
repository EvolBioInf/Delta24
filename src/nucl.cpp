#include "nucl.hpp"

using namespace std;

size_t Nucl::char2code( char c){
	size_t ret = 0;
	switch( c){
		case 'A': ret = 0; break;
		case 'C': ret = 1; break;
		case 'G': ret = 2; break;
		case 'T': ret = 3; break;
	}
	return ret;
}


char Nucl::code2char( std::size_t d){
	char c = 0;
	switch( d & 0x3UL){
		case 0 : c = 'A'; break;
		case 1 : c = 'C'; break;
		case 2 : c = 'G'; break;
		case 3 : c = 'T'; break;
	}
	return c;
}

Nucl::Nucl(size_t readID, char c): data(readID<<2) {
	data |= char2code(c);
}

Nucl::Nucl(): data(0){};

bool Nucl::operator< (const Nucl& b) const {
	return getReadID() < b.getReadID();
}

bool Nucl::operator> (const Nucl& b) const {
	return getReadID() > b.getReadID();
}
