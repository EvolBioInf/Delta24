#c++ main.cpp -I ../statgen/lib/include/ -L ../statgen/lib/ -lStatGen -lz

CXXFLAGS=-I ../lib/libStatGen/include -L ../lib/libStatGen/ -lStatGen -lz -g -std=c++0x #-lcrypto

CC=c++

Delta24: Delta24.cpp
	$(CC) Delta24.cpp $(CXXFLAGS)-O3 -o Delta24

clean:
	rm Delta24

