#include "BamReader.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>
#include <boost/iostreams/filter/gzip.hpp>
#include <memory>
#include "ContigsHg38.hpp"

using namespace std;
using namespace boost::iostreams;

BamReader::BamReader() : m_useReverseReads(false), m_done(false), m_position(-1) {}
BamReader::BamReader(const string &bamFileName, const bool useReverseReads) : m_useReverseReads(useReverseReads), m_done(false), m_position(-1) {
	setBamFile(bamFileName);
}

void BamReader::setBamFile(const string &bamFileName) {
	unique_ptr<ifstream> filePtr(new ifstream(bamFileName, ios_base::in | ios_base::binary));
	m_inputFilePtr = move(filePtr);
	m_bamStream.push(gzip_decompressor());
	m_bamStream.push(*m_inputFilePtr);

	// Make sure the magic code matches
    char bamCode[5];
    m_bamStream.sgetn(bamCode, 4);
    bamCode[4] = 0;
    if (strcmp(bamCode, "BAM\1") != 0) {
    	cerr << "Invalid BAM file: " << bamFileName << endl;
    }

    // Get through the header data so we are ready to get reads
    int32_t l_text;
    m_bamStream.sgetn((char*)&l_text, 4);
    for (int j = 0; j < l_text; ++j) m_bamStream.snextc();

    // see how many reference contigs we have
    int32_t n_ref;
    m_bamStream.sgetn((char*)&n_ref, 4);
    if (n_ref != ContigsHg38::count) {
		cerr << "Mismatch in reference contig count! Expected " << ContigsHg38::count << " but got " << n_ref << endl;
	}

	// compute offsets for all the contigs
    char buff[100];
    buff[99] = 0;
    for (int j = 0; j < n_ref; ++j) {

    	// read the data for this reference contig
    	int32_t l_name;
    	m_bamStream.sgetn((char*)&l_name, 4);
    	char* name = new char[l_name];
    	m_bamStream.sgetn(buff, min(l_name, 99));
    	string refName = buff;
    	int32_t l_ref;
    	m_bamStream.sgetn((char*)&l_ref, 4);

    	// Make sure we match what we expect
    	if (l_ref != ContigsHg38::sizes[j] || refName.compare(ContigsHg38::names[j]) != 0) {
    		cerr << "Mismatch in contigs! Expected " << j << " " << ContigsHg38::names[j] << "(" << ContigsHg38::sizes[j] <<
    			") but got " << refName << "(" << l_ref << ")" << endl;
    	}

    	// save the offset for this contig
    	m_refOffsets.push_back(accumulate(ContigsHg38::sizes, ContigsHg38::sizes+j, static_cast<int64_t>(0)));
    }
    advance();
}

int64_t BamReader::position() const {
	return m_position;
}

void BamReader::advance() {
	try {

		while (!m_done) {

			// See if we have reached the end of the file
			if (m_inputFilePtr->eof()) {
		    	m_done = true;
		    	m_position = -1;
		    	m_bamStream.reset();
		    	return;
		    }

		    int32_t block_size;
		    m_bamStream.sgetn((char*)&block_size, 4);

		    // get the reference contig this read maps to
		    int32_t refID;
		    m_bamStream.sgetn((char*)&refID, 4);

	    	// get the read position
		    int32_t pos;
		    m_bamStream.sgetn((char*)&pos, 4);
		    if (refID != -1) m_position = static_cast<int64_t>(pos) + m_refOffsets[refID];

		    // skip the bin_mq_nl field
		    for (int k = 0; k < 4; ++k) m_bamStream.snextc();

		   	// see if we are reverse complemented
		    uint32_t flag_nc;
		    m_bamStream.sgetn((char*)&flag_nc, 4);
		    const bool forward = (flag_nc & REVERSE_FLAG) == 0;

		    // skip the rest of the entry
		    for (int k = 0; k < block_size-16; ++k) {
		    	if (m_inputFilePtr->eof()) break;
		    	m_bamStream.snextc();
		    }

		    // break if we found a read in the right direction
		    if (refID != -1 && ((forward && !m_useReverseReads) || (!forward && m_useReverseReads))) return;
		}

	// get the gzip error at the end of the file I can't figure out
	} catch (boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::iostreams::gzip_error> > &e) {
		m_done = true;
    	m_position = -1;
    	m_bamStream.reset();
	}
}
