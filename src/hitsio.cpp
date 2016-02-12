/*
 * hitsio: reading and writing hits files
 */

#include "hitsio.hpp"

template <typename T>
std::string to_string(T const& value) {
  stringstream sstr;
  sstr << value;
  return sstr.str();
}

double string_to_double( const std::string& s ) {
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return 0;
  return x;
} 

#define WRITE_UINT(STREAM, VAL) { \
	/* rationale for writing numeric values as uint32_t: size_t is usually defined as uint32_t */ \
	if (VAL > numeric_limits<uint32_t>::max()) { \
		cerr << "Numeric value overflow when writing the hits file" << endl; \
		exit(1); \
	} \
	uint32_t BUF = (VAL); \
	(STREAM).write((char*)&(BUF), sizeof(uint32_t)); \
};

#define READ_UINT(STREAM, VAL) { \
	(STREAM).read((char*)&(VAL), sizeof(uint32_t)); \
};

#define WRITE_UINT_SMALL(STREAM, VAL) { \
	if ((VAL) < 255) { \
		unsigned char BUF = (unsigned char)(VAL); \
		(STREAM).write((char*)&(BUF), sizeof(unsigned char)); \
	} else { \
		unsigned char BUF = 255; \
		(STREAM).write((char*)&(BUF), sizeof(unsigned char)); \
		WRITE_UINT(STREAM, VAL); \
	} \
};

#define READ_UINT_SMALL(STREAM, VAL) { \
	unsigned char BUF; \
	(STREAM).read((char*)&(BUF), sizeof(unsigned char)); \
	if (BUF == 255) { \
		READ_UINT(STREAM, VAL); \
	} else { \
		VAL = BUF; \
	}; \
};

#define WRITE_STR(STREAM, STR) { \
	(STREAM) << (STR) << endl; \
};

#define READ_STR(STREAM, VAR) { \
	getline((STREAM), (VAR)); \
};

/*
WRITE_STR_DELTA and READ_STR_DELTA

Write/read (STR) to (STREAM) using a simple delta encoding against the string deltaBuffer:
1) count number matching characters in (STR) and deltaBuffer
2) replace the matching characters by their counts

Both assume the existence of the following variables (declared as member variables of HitsReader/HitsWriter):
    size_t nMatchesBeg;
    size_t nMatchesEnd;
	string deltaBuffer = "";
*/
#define WRITE_STR_DELTA(STREAM, STR) { \
	nMatchesBeg = 0; \
	nMatchesEnd = 0; \
	/* count number of common characters in the beginning */ \
	while((nMatchesBeg < min(deltaBuffer.size(), (STR).size())) && \
		(deltaBuffer[nMatchesBeg] == (STR)[nMatchesBeg])) { nMatchesBeg++; } \
	/* count number of common characters in the end */ \
	while(((nMatchesBeg + nMatchesEnd) < min(deltaBuffer.size(), (STR).size())) && \
		(deltaBuffer[deltaBuffer.size() - 1 - nMatchesEnd] == \
			(STR[(STR).size() - 1 - nMatchesEnd]))) { nMatchesEnd++; } \
	/* if no similarities: write the whole read_id as normal */ \
	if ((nMatchesBeg == 0) && (nMatchesEnd == 0)) { WRITE_STR(STREAM, STR) } \
	/* else: write an empty string followed by the delta-encoded string */ \
	else { \
		WRITE_STR(STREAM, ""); \
		WRITE_UINT_SMALL(STREAM, nMatchesBeg); \
		/*WRITE_UINT(STREAM, nMatchesBeg);*/ \
		WRITE_STR(STREAM, (STR).substr(nMatchesBeg, (STR).size() - nMatchesBeg - nMatchesEnd)); \
		WRITE_UINT_SMALL(STREAM, nMatchesEnd); \
		/*WRITE_UINT(STREAM, nMatchesEnd);*/ \
	} \
	deltaBuffer = STR; \
};

#define READ_STR_DELTA(STREAM, VAR) { \
	READ_STR((STREAM), (VAR)); \
	if ((VAR) == "") { \
		READ_UINT_SMALL((STREAM), nMatchesBeg); \
		/*READ_UINT((STREAM), nMatchesBeg);*/ \
		READ_STR((STREAM), VAR); \
		READ_UINT_SMALL((STREAM), nMatchesEnd); \
		/*READ_UINT((STREAM), nMatchesEnd);*/ \
		if (!(STREAM.eof())) { \
			(VAR) = deltaBuffer.substr(0, nMatchesBeg) + (VAR) + \
				deltaBuffer.substr(deltaBuffer.size() - nMatchesEnd, nMatchesEnd); \
		} \
	} \
	deltaBuffer = (VAR); \
};

HitsfileWriter::HitsfileWriter(string argHitsfileFormat) {
	if (argHitsfileFormat[0] == 't') {
		hitsfileSchema = 0;
	} else if (argHitsfileFormat[0] == 'b') {
		hitsfileSchema = 1;
	// no file format specified => select default (binary)
	} else {
		hitsfileSchema = 1;
	}
	// if output is a binary format: initialise compression
	if (hitsfileSchema == 1) { ofs.push(boost::iostreams::zlib_compressor(boost::iostreams::zlib::best_speed)); }
	ofs.push(cout);
}

void HitsfileWriter::addTranscriptMetaData(string argTranscriptName, double argTranscriptEffectiveLength, int argTranscriptTrueLength) {
	transcriptName.push_back(argTranscriptName);
	transcriptEffectiveLength.insert(make_pair(argTranscriptName, argTranscriptEffectiveLength));
	transcriptTrueLength.insert(make_pair(argTranscriptName, argTranscriptTrueLength));
	// The "straightforward way":
	//     transcriptToIndex[argTranscriptName] = transcriptToIndex.size();
	// does not seem to work: at least in one implementation, transcriptToIndex.size() is evaluated after
	// the new element is added to transcriptToIndex.
	size_t transcriptToIndexSize = transcriptToIndex.size();
	transcriptToIndex[argTranscriptName] = transcriptToIndexSize;
}

void HitsfileWriter::addGeneIsoformRecord(string geneName) {
	vector<string> currentGeneIsoforms;
    geneIsoforms.insert(make_pair(geneName, currentGeneIsoforms));
	currentGeneName = geneName;
}

void HitsfileWriter::addTranscriptToGeneIsoformRecord(string transcriptName) {
    (geneIsoforms[currentGeneName]).push_back(transcriptName);
}

void HitsfileWriter::addIdenticalTranscriptsRecord() {
	vector<string> currentIdenticalTranscripts;
	identicalTranscripts.push_back(currentIdenticalTranscripts);
}

void HitsfileWriter::addTranscriptToIdenticalTranscriptsRecord(string transcriptName) {
	identicalTranscripts[identicalTranscripts.size() - 1].push_back(transcriptName);
}

void HitsfileWriter::writeHeaderSchema0() {
	for(int i = 0; i < transcriptName.size(); i++) {
		ofs << "@TranscriptMetaData\t" << transcriptName[i] << "\t" << transcriptEffectiveLength[ transcriptName[i] ] 
        << "\t" << transcriptTrueLength[ transcriptName[i] ] << endl;
	}
	for(map<string, vector<string> >::iterator i = geneIsoforms.begin(); i != geneIsoforms.end(); i++) {
		ofs << "@GeneIsoforms\t" << i->first;
		for(vector<string>::iterator j = (i->second).begin(); j != (i->second).end(); j++)
			ofs  << "\t" << (*j);
		ofs << endl;
	}
	for(vector<vector<string> >::iterator i = identicalTranscripts.begin(); i != identicalTranscripts.end(); i++) {
		ofs << "@IdenticalTranscripts";
		for(vector<string>::iterator j = i->begin(); j != i->end(); j++) {
			ofs << "\t" << (*j);
		}
		ofs << endl;
	}
}

void HitsfileWriter::writeReadMapRecordSchema0() {
	ofs << ">" << currentReadName << endl;
	for(int i = 0; i < currentReadTranscripts.size(); i++) {
		ofs << currentReadTranscripts[i] << endl;
	}
}

void HitsfileWriter::writeHeaderSchema1() {
	ofs << MMSEQ_HEADER << endl;
	WRITE_UINT(ofs, hitsfileSchema);
	WRITE_UINT(ofs, transcriptName.size());
	for(int i = 0; i < transcriptName.size(); i++) {
		WRITE_STR(ofs, transcriptName[i]);
    WRITE_STR(ofs, to_string(transcriptEffectiveLength[transcriptName[i]]));
		WRITE_UINT(ofs, transcriptTrueLength[transcriptName[i] ])
	}
	WRITE_UINT(ofs, geneIsoforms.size());
	for(map<string, vector<string> >::iterator i = geneIsoforms.begin(); i != geneIsoforms.end(); i++) {
		WRITE_STR(ofs, i->first);
		WRITE_UINT(ofs, (i->second).size());
		for(vector<string>::iterator j = (i->second).begin(); j != (i->second).end(); j++) {
			WRITE_STR(ofs, *j);
		}
	}
	WRITE_UINT(ofs, identicalTranscripts.size());
	for(vector<vector<string> >::iterator i = identicalTranscripts.begin(); i != identicalTranscripts.end(); i++) {
		WRITE_UINT(ofs, i->size());
		for(vector<string>::iterator j = i->begin(); j != i->end(); j++) {
			WRITE_STR(ofs, *j);
		}
	}
}

void HitsfileWriter::writeHeader() {
	if (hitsfileSchema == 0) {
		writeHeaderSchema0();
	} else if (hitsfileSchema == 1) {
		writeHeaderSchema1();
	}
}

void HitsfileWriter::addReadMapRecord(string readName) {
	currentReadName = readName;
	currentReadTranscripts.clear();
}

void HitsfileWriter::addTranscriptToReadMapRecord(string transcriptName) {
	currentReadTranscripts.push_back(transcriptName);
}

void HitsfileWriter::writeReadMapRecordSchema1() {
	WRITE_STR_DELTA(ofs, currentReadName);
	WRITE_UINT(ofs, currentReadTranscripts.size());
	size_t matchIndex;
	for(int i = 0; i < currentReadTranscripts.size(); i++) {
		matchIndex = transcriptToIndex[ currentReadTranscripts[i] ];
		WRITE_UINT(ofs, matchIndex);
	}
}

void HitsfileWriter::writeReadMapRecord() {
	if (hitsfileSchema == 0) {
		writeReadMapRecordSchema0();
	} else if (hitsfileSchema == 1) {
		writeReadMapRecordSchema1();
	}
}

HitsfileReader::HitsfileReader(string fileName) {
	ifsFile.open(fileName.c_str(), ios_base::in | ios_base::binary);
	if(!ifsFile.good()) {
		cerr << "Error reading hits file \"" << fileName << "\".\n";
		exit(1);
	}
	// "lazy but sufficient" way to check if the file is gz-compressed
	//if (ifsFile.peek() == 0x1F) { ifs.push(boost::iostreams::gzip_decompressor()); }
	if (ifsFile.peek() == 0x78) { ifs.push(boost::iostreams::zlib_decompressor()); }
	ifs.push(ifsFile);

	// detect schema of input file by looking at the first line of the file
	string line;
	getline(ifs, line);
	stringstream tokens(line);
	string token1;
	tokens >> token1;
	if (token1 == "@TranscriptMetaData") {
		hitsfileSchema = 0;
	} else {
		READ_UINT(ifs, hitsfileSchema);
	}
	if ((hitsfileSchema < 0) || (hitsfileSchema > 4)) {
		cerr << "Input file \"" << fileName << "\" does not seem to be a hits file.\n";
		exit(1);
	}
	// ifs.seekg(ios::beg)
    // "brute force replacement" for the previous line, which does not seem to have the intended effect
    ifsFile.close();
    ifs.reset();
    ifsFile.open(fileName.c_str(), ios_base::in | ios_base::binary);
    //if (ifsFile.peek() == 0x1F) { ifs.push(boost::iostreams::gzip_decompressor()); }
	if (ifsFile.peek() == 0x78) { ifs.push(boost::iostreams::zlib_decompressor()); }
    ifs.push(ifsFile);
}

void HitsfileReader::readHeaderSchema0(vector<string> *transcriptName, map<string, double> *transcriptEffectiveLength,
    map<string, int> *transcriptTrueLength,
		map<string, vector<string> > *geneIsoforms, vector<vector<string> > *identicalTranscripts) {
	string line;
	vector<string> iReadMatches;
	while (ifs && (ifs.peek() != '>')) {
		getline(ifs, line);
		stringstream tokens(line);
		string token1;
		tokens >> token1;
		if (token1 == "@TranscriptMetaData") {
			string iTranscriptName;
			double iTranscriptEffectiveLength;
			int iTranscriptTrueLength;
			tokens >> iTranscriptName >> iTranscriptEffectiveLength >> iTranscriptTrueLength;
			transcriptName->push_back(iTranscriptName);
			transcriptEffectiveLength->insert(make_pair(iTranscriptName, iTranscriptEffectiveLength));
			transcriptTrueLength->insert(make_pair(iTranscriptName, iTranscriptTrueLength));
		} else if (token1 == "@GeneIsoforms") {
			string geneID;
			string iTranscriptID;
	        vector<string> transcriptID;
	        tokens >> geneID;
			tokens >> iTranscriptID;
	        while (tokens) {
				transcriptID.push_back(iTranscriptID);
				tokens >> iTranscriptID;
			}
	        geneIsoforms->insert(make_pair(geneID, transcriptID));
		} else if (token1 == "@IdenticalTranscripts") {
			string iTranscriptID;
	        vector<string> transcriptID;
			tokens >> iTranscriptID;
	        while (tokens) {
				transcriptID.push_back(iTranscriptID);
				tokens >> iTranscriptID;
			}
	        identicalTranscripts->push_back(transcriptID);
		} else {
			cerr << "Hits file looks malformed.\n";
			exit(1);
		}
	}
}

bool HitsfileReader::readReadMapRecordReadIDSchema0(string &readID) {
	if (ifs.eof()) { return false; }
	getline(ifs, readID);
	assert(readID[0] == '>');
	readID = readID.substr(1);
  if (ifs.peek() == std::ifstream::traits_type::eof()) {
    cerr << "Warning: read record without any mapping transcripts found"
         << " at the end of the hits file. The hits file may be corrupted.\n";
    return false;
  }
	return true;
}

bool HitsfileReader::readReadMapRecordTranscriptIDSchema0(string &transcriptID) {
	if ((ifs.eof()) || (ifs.peek() == '>')) { return false; }
	return(static_cast<bool>(getline(ifs, transcriptID)));
}

void HitsfileReader::readHeaderSchema1(vector<string> *transcriptName, map<string, double> *transcriptEffectiveLength,
    map<string, int> *transcriptTrueLength,
		map<string, vector<string> > *geneIsoforms, vector<vector<string> > *identicalTranscripts) {
	string schema_id;
	getline(ifs, schema_id);
	READ_UINT(ifs, hitsfileSchema);
	int nextInt;
	size_t nextStrLen;
	string nextStr;
	string nextStr2;
	/* read @TranscriptMetaData */
	READ_UINT(ifs, nextInt); int nTranscriptName = nextInt;
	for(int i = 0; i < nTranscriptName; i++) {
		READ_STR(ifs, nextStr);
		transcriptName->push_back(nextStr);
		headerTranscriptName.push_back(nextStr);
		READ_STR(ifs, nextStr2);
		transcriptEffectiveLength->insert(make_pair(nextStr, string_to_double(nextStr2)));
		READ_UINT(ifs, nextInt);
		transcriptTrueLength->insert(make_pair(nextStr, nextInt));
	}
	/* read @GeneIsoforms */
	READ_UINT(ifs, nextInt);
	int nGeneIsoforms = nextInt;
	for (int i = 0; i < nGeneIsoforms; i++) {
		READ_STR(ifs, nextStr);
		string geneID = nextStr;
		READ_UINT(ifs, nextInt);
		int nGeneTranscripts = nextInt;
		vector<string> transcriptID;
		for(int j = 0; j < nGeneTranscripts; j++) {
			READ_STR(ifs, nextStr);
			transcriptID.push_back(nextStr);
		}
        geneIsoforms->insert(make_pair(geneID, transcriptID));
	}
	/* read @IdenticalTranscripts */ \
	READ_UINT(ifs, nextInt);;
	int nIdenticalTranscripts = nextInt;
	for (int i = 0; i < nIdenticalTranscripts; i++) {
		READ_UINT(ifs, nextInt);
		int nTranscripts = nextInt;
        vector<string> transcriptID;
        for(int j = 0; j < nTranscripts; j++) {
        	READ_STR(ifs, nextStr);
			transcriptID.push_back(nextStr);
		}
        identicalTranscripts->push_back(transcriptID);
	}
}

void HitsfileReader::readHeader(vector<string> *transcriptName, map<string, double> *transcriptEffectiveLength,
    map<string, int> *transcriptTrueLength,
		map<string, vector<string> > *geneIsoforms, vector<vector<string> > *identicalTranscripts) {
	if (hitsfileSchema == 0) {
		readHeaderSchema0(transcriptName, transcriptEffectiveLength, transcriptTrueLength, geneIsoforms, identicalTranscripts);
	} else if (hitsfileSchema == 1) {
		readHeaderSchema1(transcriptName, transcriptEffectiveLength, transcriptTrueLength, geneIsoforms, identicalTranscripts);
	} else {
		cerr << "We should never get to this state!\n";
		exit(1);
	}
}

bool HitsfileReader::readReadMapRecordReadIDSchema1(string &readID) {
	READ_STR_DELTA(ifs, nextStr);
	if (ifs.eof()) { return false; }
	readID = nextStr;
	READ_UINT(ifs, nextInt);
	countReadMapRecord = nextInt;
	return true;
}

bool HitsfileReader::readReadMapRecordReadID(string &readID) {
	if (hitsfileSchema == 0) {
		return (readReadMapRecordReadIDSchema0(readID));
	} else if (hitsfileSchema == 1) {
		return (readReadMapRecordReadIDSchema1(readID));
	} else {
		cerr << "We should never get to this state!\n";
		exit(1);
	}
}

bool HitsfileReader::readReadMapRecordTranscriptIDSchema1(string &transcriptID) {
	if (ifs.eof() or (countReadMapRecord == 0)) { return false; }
	READ_UINT(ifs, nextInt);
	transcriptID = headerTranscriptName[nextInt];
	countReadMapRecord--;
	return true;
}

bool HitsfileReader::readReadMapRecordTranscriptID(string &readID) {
	if (hitsfileSchema == 0) {
		return (readReadMapRecordTranscriptIDSchema0(readID));
	} else if (hitsfileSchema == 1) {
		return (readReadMapRecordTranscriptIDSchema1(readID));
	}
}
