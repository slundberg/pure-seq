#include <memory>
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <tclap/CmdLine.h>
#include <armadillo>

#include "BamReader.hpp"
#include "FeatureBlocks.hpp"

int main(int argc, char **argv) {
	using namespace std;
	using namespace boost::iostreams;
	try {
		// Read in the command line arguments
		TCLAP::CmdLine cmd("Builds XtX data matrix.", ' ', "0.1");
		TCLAP::ValueArg<std::string> targetArg("t","target","Target",true,"homer","target BAM file");
	    cmd.add(targetArg);
		TCLAP::UnlabeledMultiArg<string> controlArg("controls","Control BAM files",true,"control BAM file");
	    cmd.add(controlArg);
	    cmd.parse(argc, argv);
	    vector<string> controlFileNames = controlArg.getValue();
	    BamReader targetReader(targetArg.getValue(), false);
	    
	    // open the positiom streams to the bam files
	    vector<BamReader> controlReaders(controlFileNames.size());
	    for (int i = 0; i < controlReaders.size(); ++i) {
	    	controlReaders[i].setBamFile(controlFileNames[i]);
	    }

	    FeatureBlocks featureBlocks(targetReader, controlReaders, 1000000);

	    arma::mat C(featureBlocks.blockWidth(), featureBlocks.blockWidth());
	    C.zeros();

	    int count = 0;
	    while (!featureBlocks.done()) {
	    	const arma::mat &block = featureBlocks.block();
	    	C = C + block.t()*block;
	    	featureBlocks.advance();
	    	if (count % 100 == 0) {
	    		cout << "finished block " << count << endl;
	    		cout << C << endl;
	    	}
	    	count += 1;
	    }

	    cout << "done!" << endl;
	    
		return 0;

	} catch (TCLAP::ArgException &e) {
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
	}
}