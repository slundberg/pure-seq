#ifndef FEATUREBLOCKS_H
#define FEATUREBLOCKS_H

#include <armadillo>
#include <vector>
#include "BamReader.hpp"

class FeatureBlocks {
private:
	std::vector<arma::mat> m_blocks;
	std::vector<BamReader> &m_controlReaders;
	BamReader &m_targetReader;

	bool m_done;
	int64_t m_position;
	const int m_blockSize;
	const int m_blockWidth;

	int m_prevInd;
	int m_currInd;
	int m_nextInd;

	void fillBlock(const int index);
	void computeCurrBlockFeatures();
	void buildSmoothedFeature(int srcCol, int destCol, int windowSize);

public:

	FeatureBlocks(BamReader &targetReader, std::vector<BamReader> &controlReaders, const int blockSize);

	int blockWidth() const { return m_blockWidth; };
	int blockSize() const { return m_blockSize; };
	bool done() const { return m_done; };
	void advance();
	const arma::mat &block() const { return m_blocks[m_currInd]; };
};

#endif