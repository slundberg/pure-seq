#include "FeatureBlocks.hpp"

using namespace std;

FeatureBlocks::FeatureBlocks(BamReader &targetReader, vector<BamReader> &controlReaders, const int blockSize) : 
		m_targetReader(targetReader), m_controlReaders(controlReaders), m_position(0), m_blockSize(blockSize),
		m_blockWidth(controlReaders.size()*4 + 2) {

	// create the prev, curr, and next blocks
	m_blocks.push_back(arma::zeros<arma::mat>(m_blockSize, m_blockWidth));
	m_blocks.push_back(arma::zeros<arma::mat>(m_blockSize, m_blockWidth));
	m_blocks.push_back(arma::zeros<arma::mat>(m_blockSize, m_blockWidth));

	fillBlock(1);
	fillBlock(2);
	

	m_prevInd = 0;
	m_currInd = 1;
	m_nextInd = 2;

	computeCurrBlockFeatures();
}

void FeatureBlocks::advance() {

	// advance the block cache indexes
	m_prevInd = (m_prevInd + 1)%3;
	m_currInd = (m_currInd + 1)%3;
	m_nextInd = (m_nextInd + 1)%3;

	// fill the next block with new values
	fillBlock(m_nextInd);
	computeCurrBlockFeatures();
}

void FeatureBlocks::fillBlock(const int index)
{
	arma::mat &block = m_blocks[index];
	block.zeros();
	//cout << "1 A" << endl;
	m_done = true;

	// Fill in the target column
	while (m_targetReader.position() != -1 && m_targetReader.position() < m_position + m_blockSize) {
		block(m_targetReader.position() - m_position,0) += 1;
		//cout << "target pos " << m_targetReader.position() << endl;
		m_targetReader.advance();
		m_done = false;
	}
	//cout << "1 B " << block.col(0) << endl << index << endl;
	// Fill in the control columns
	for (int i = 0; i < m_controlReaders.size(); ++i) {
		BamReader &reader = m_controlReaders[i];

		while (reader.position() != -1 && reader.position() < m_position + m_blockSize) {
			block(reader.position() - m_position,i+1) += 1;
			reader.advance();
			m_done = false;
		}
	}

	// See if we are really done or just found a blank block
	if (m_done) {
		m_done = m_done && m_targetReader.position() == -1;
		for (int i = 0; i < m_controlReaders.size(); ++i) {
			m_done = m_done && m_controlReaders[i].position() == -1;
		}
	}

	m_position += m_blockSize;
}

void FeatureBlocks::computeCurrBlockFeatures() {

	// build smoothed features
	for (int i = 0; i < m_controlReaders.size(); ++i) {
		buildSmoothedFeature(i+1, i+1+m_controlReaders.size(), 10);
		buildSmoothedFeature(i+1, i+1+2*m_controlReaders.size(), 100);
		buildSmoothedFeature(i+1, i+1+3*m_controlReaders.size(), 1000);
	}

	// build the constant feature
	arma::mat &block = m_blocks[m_currInd];
	for (int i = 0; i < m_blockSize; ++i) block(i, m_blockWidth-1) = 1;
}

void FeatureBlocks::buildSmoothedFeature(int srcCol, int destCol, int windowSize) {
	arma::mat &prevBlock = m_blocks[m_prevInd];
	arma::mat &currBlock = m_blocks[m_currInd];
	arma::mat &nextBlock = m_blocks[m_nextInd];
    double windowSum = 0;

    for (int i = -2*windowSize; i < m_blockSize; ++i) {
        
        // add new positions to the window
        if (i+windowSize < 0) {
        	//cout << "adding " << i+windowSize+m_blockSize << " prev" << endl;
            windowSum += prevBlock(i+windowSize+m_blockSize,srcCol);
        } else if (i+windowSize < m_blockSize) {
            windowSum += currBlock(i+windowSize,srcCol);
            //cout << "adding " << i+windowSize << " curr" << endl;
        } else {
            windowSum += nextBlock(i+windowSize-m_blockSize,srcCol);
            //cout << "adding " << i+windowSize-m_blockSize << " next" << endl;
        }

        
        if (i >= 0) {

        	// store the current window average
        	currBlock(i,destCol) = windowSum/(2*windowSize+1);

        	// remove old positions
            if (i-windowSize < 0) {
            	//cout << "removing " << i-windowSize+m_blockSize << " prev " << -windowSize << " " << i << endl;
                windowSum -= prevBlock(i-windowSize+m_blockSize,srcCol);
            } else if (i-windowSize < m_blockSize) {
            	//cout << "removing " << i-windowSize << " curr" << endl;
                windowSum -= currBlock(i-windowSize,srcCol);
            } else {
            	//cout << "removing " << i-windowSize-m_blockSize << " next" << endl;
                windowSum -= nextBlock(i-windowSize-m_blockSize,srcCol);
            }
        }
        
        
    }
    //cout << "done buildSmoothedFeature" << endl;
}