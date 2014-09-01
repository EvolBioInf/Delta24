#include <vector>
#include <string>
#include <functional>
#include <stdexcept>

#include "SamFile.h"
#include "SamFlag.h"
#include "SamValidation.h"

#include "mapNucl.hpp"

using namespace std;

const unsigned char MIN_QUAL = '!';

/** @brief Maps nucleotides to their corresponding position on the
 * reference sequence.
 *
 * This function reads a BAM/SAM file and maps the reads contained
 * within. The resulting datatype is pos -> [(readID, Nuc)] in SML
 * notation. Basically, its a vector with the length of the reference
 * sequence. At each position there is another vector of all reads 
 * which have a nucleotide mapping to that position.
 *
 * @param filename - The path to the BAM/SAM file.
 */
mapped_nucl_t mapNucl( const char * filename){
	SamFile file(ErrorHandler::RETURN);
	SamFileHeader header;

	if( !file.OpenForRead(filename) || !file.ReadHeader(header) ){
		throw runtime_error("Failed to open input BAM file.");
	}

	// Get the Reference Sequence (SQ)
	// ATM only the first SQ is read.
	SamHeaderRecord *header_record = header.getNextSQRecord();
	ssize_t ref_length = atoi(header_record->getTagValue("LN"));

	mapped_nucl_t ret { static_cast<size_t>(ref_length), mapped_nucl_t::value_type()};

	SamRecord record, next_record;

	// Get the next read.
	while( file.ReadRecord(header, record)){
		// filter for unmapped reads and bad alignments
		auto flag = record.getFlag();
		bool requirements = SamFlag::isMapped(flag) && !(flag & SamFlag::SECONDARY_ALIGNMENT) &&
			!SamFlag::isQCFailure(flag) && !SamFlag::isDuplicate(flag) &&
			!(flag & 0x800 /* supplementary alignment */) &&
			(SamFlag::isPaired(flag) ? SamFlag::isProperPair(flag) : true);

		if( !requirements ) continue;

		/* For the computation we need to take track of nucleotides that come from
		 * the same read. Since two reads can be pairs, but reside at different
		 * places within a BAM file, this read information must be carried along
		 * through the whole computation.
		 * Here I use a hash over the read name (a string) to compute an ID (64bit).
		 */
		string str (record.getReadName());
		std::hash<string> str_hash;
		size_t readID = str_hash(str);

		ssize_t start = record.get0BasedPosition();
		ssize_t length = record.getReadLength();

		// skip broken records.
		if( start < 0 || length < 0 || strcmp(record.getReferenceName(), "*") == 0 ){
			continue;	
		}

		str = record.getQuality();

		// iterate over the read.
		for( ssize_t i = 0; i < length; i++){
			int pos_in_ref = record.getCigarInfo()->getRefPosition(i, start);

			// skip indels
			if( pos_in_ref < 0 || pos_in_ref >= ref_length){
				continue;
			}

			// skip low quality nucleotides
			if( str[i] < MIN_QUAL){
				continue;
			}

			// add the current nucleotide to the map.
			Nucl p {readID, record.getSequence(i)};

			ret[pos_in_ref].push_back(p);
		}
	}

	file.Close();

	return ret;
}
