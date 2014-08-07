#include <tuple>
#include <stdlib.h>
#include <array>
#include <unordered_map>
#include <iostream>

#include "matHash.hpp"

using namespace std;

/* Simple standard constructor */
matHash::matHash(){
	this->hmap = new maptype;
}

/**
 * Increase the count of a pattern
 */
void matHash::inc( unsigned int *key){
	// derive the key for the map
	array<unsigned int, 24> K;
	for (size_t i = 0; i < 24; ++i){
		K[i] = key[i];
	}

	// check if the key is already in the map
	auto elem = this->hmap->find(K);
	if( elem != this->hmap->end() ){
		// increase the count
		elem->second[24]++;
	} else {
		// create a new value to be saved in the map
		auto value = (float*) new float[25];
		for (size_t i = 0; i < 24; ++i){
			value[i] = float(key[i]);
		}
		value[24] = 1.0;

		// insert the new value
		this->hmap->insert(maptype::value_type(K, value ));
	}
}

void matHash::inc( const array<unsigned int, 24>& key){
	// check if the key is already in the map
	auto elem = this->hmap->find(key);
	if( elem != this->hmap->end() ){
		// increase the count
		elem->second[24]++;
	} else {
		// create a new value to be saved in the map
		auto value = (float*) new float[25];
		for (size_t i = 0; i < 24; ++i){
			value[i] = float(key[i]);
		}
		value[24] = 1.0;

		// insert the new value
		this->hmap->insert(maptype::value_type(key, value ));
	}
}

/**
 * Get rid of all data
 */
void matHash::clear(){
	auto it = this->hmap->begin();
	for(; it != this->hmap->end(); it++){
		delete (*it).second;
	} 
	
	this->hmap->clear();
}

maptype::iterator matHash::begin(){
	return this->hmap->begin();
}

maptype::iterator matHash::end(){
	return this->hmap->end();
}
