/*
 * hitsio: reading and writing hits files
 */

#include <boost/serialization/array_wrapper.hpp>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <math.h>
#include <string>
#include <string.h>
#include <fstream>
#include <assert.h>
#include <algorithm>
#include <time.h>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

#define LOG(ARG) cerr << __FILE__ << ":" << __LINE__ << " " << (ARG) << endl

#define MMSEQ_HEADER "MMSEQ_HITSFILE"

class HitsfileWriter {
private:
	boost::iostreams::filtering_ostream ofs;
	int hitsfileSchema;
	vector<string> transcriptName;
	map<string, double> transcriptEffectiveLength;
	map<string, int> transcriptTrueLength;
  map<string, int> transcriptToIndex;
	map<string, vector<string> > geneIsoforms;
	string currentGeneName;
	vector<vector<string> > identicalTranscripts;
	string currentReadName;
	vector<string> currentReadTranscripts;
	// internal variables of the string delta encoding macros:
	size_t nMatchesBeg;
    size_t nMatchesEnd;
	string deltaBuffer;
public:
	HitsfileWriter(string argHitsfileFormat);
	void addTranscriptMetaData(string iTranscriptName, double iTranscriptEffectiveLength, int iTranscriptTrueLength);
	void addGeneIsoformRecord(string geneName);
	void addTranscriptToGeneIsoformRecord(string transcriptName);
	void addIdenticalTranscriptsRecord();
	void addTranscriptToIdenticalTranscriptsRecord(string transcriptName);
	void writeHeaderSchema0();
	void writeHeaderSchema1();
	void writeHeader();
	void addReadMapRecord(string readName);
	void addTranscriptToReadMapRecord(string transcriptName);
	void writeReadMapRecordSchema0();
	void writeReadMapRecordSchema1();
	void writeReadMapRecord();
};

class HitsfileReader {
private:
	vector<string> headerTranscriptName;
	boost::iostreams::filtering_istream ifs;
	ifstream ifsFile;
	int hitsfileSchema;
	int countReadMapRecord;
	int nextInt;
	string nextStr;
	// internal variables of the string delta encoding macros:
	size_t nMatchesBeg;
    size_t nMatchesEnd;
	string deltaBuffer;
public:
	HitsfileReader(string fileName);
	//~HitsfileReader();
	void readHeaderSchema0(vector<string> *transcriptName, map<string, double > *transcriptEffectiveLength, map<string,int> *transcriptTrueLength,
			map<string, vector<string> > *geneIsoforms, vector<vector<string> > *identicalTranscripts);
	void readHeaderSchema1(vector<string> *transcriptName, map<string, double> *transcriptEffectiveLength, map<string, int> *transcriptTrueLength,
			map<string, vector<string> > *geneIsoforms, vector<vector<string> > *identicalTranscripts);
	void readHeader(vector<string> *transcriptName, map<string, double> *transcriptEffectiveLength, map<string, int> *transcriptTrueLength,
			map<string, vector<string> > *geneIsoforms, vector<vector<string> > *identicalTranscripts);
	bool readReadMapRecordReadIDSchema0(string &readID);
	bool readReadMapRecordReadIDSchema1(string &readID);
	bool readReadMapRecordReadID(string &readID);
	bool readReadMapRecordTranscriptIDSchema0(string &transcriptID);
	bool readReadMapRecordTranscriptIDSchema1(string &transcriptID);
	bool readReadMapRecordTranscriptID(string &transcriptID);
};
