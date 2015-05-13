#ifndef BAMREADER_H
#define BAMREADER_H

#include <boost/iostreams/filtering_streambuf.hpp>
#include <fstream>
#include <memory>
#include <vector>

class BamReader {
private:
	std::unique_ptr<std::ifstream> m_inputFilePtr;
	boost::iostreams::filtering_streambuf<boost::iostreams::input> m_bamStream;
	bool m_useReverseReads;
	bool m_done;
	int64_t m_position;
	std::vector<int64_t> m_refOffsets;

public:
	const uint32_t REVERSE_FLAG = 16 << 16;
	
	BamReader();
	BamReader(const std::string &bamFileName, const bool useReverseReads);

	void setUseReverseReads(const bool useReverseReads);
	void setBamFile(const std::string &bamFileName);

	int64_t position() const; // returns -1 when we are done
	void advance();
};

#endif