/*
 * hitstools: command-line utility for inspecting and converting hitsfiles
 */

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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "hitsio.hpp"

using namespace std;

void printUsage(char *bin, ostream& out) {
	out << "Utilities for handling hitsfile (.hits) and bitsfile (.bits, `binary hitsfile`) formats." << endl
		<< endl
		<< "View a summary of a hits file:" << endl
		<< "  hitstools inspect in.hits" << endl
    << "Output header:" <<endl
    << " hitstools header in.hits" << endl
		<< "Convert from text format to binary format:" << endl
		<< "  hitstools b in_text.hits > out_binary.hits" << endl
		<< "Convert from binary format to text format:" << endl
		<< "  hitstools t in_binary.hits > out_text.hits" << endl
		<< endl;
}

void convert(char *hits_file, string hitsfileFormat) {
	vector<string> transcriptName;
	map<string, double> transcriptEffectiveLength;
	map<string, int> transcriptTrueLength;
	map<string, vector<string> > geneIsoforms;
	vector<vector<string> > identicalTranscripts;
	HitsfileReader hitsfileReader(hits_file);
	hitsfileReader.readHeader(&transcriptName, &transcriptEffectiveLength, &transcriptTrueLength, &geneIsoforms, &identicalTranscripts);
	HitsfileWriter hitsfileWriter(hitsfileFormat);
	for(int i = 0; i < transcriptName.size(); i++) {
		hitsfileWriter.addTranscriptMetaData(transcriptName[i], transcriptEffectiveLength[ transcriptName[i] ],
      transcriptTrueLength[ transcriptName[i] ]);
	}
	for(map<string, vector<string> >::iterator i = geneIsoforms.begin(); i != geneIsoforms.end(); i++) {
		hitsfileWriter.addGeneIsoformRecord( i->first );
		for(vector<string>::iterator j = (i->second).begin(); j != (i->second).end(); j++)
			hitsfileWriter.addTranscriptToGeneIsoformRecord(*j);
	}
	for(vector<vector<string> >::iterator i = identicalTranscripts.begin(); i != identicalTranscripts.end(); i++) {
		hitsfileWriter.addIdenticalTranscriptsRecord();
		for(vector<string>::iterator j = i->begin(); j != i->end(); j++) {
			hitsfileWriter.addTranscriptToIdenticalTranscriptsRecord(*j);
		}
	}
	hitsfileWriter.writeHeader();
	string readID, transcriptID;
	while(hitsfileReader.readReadMapRecordReadID(readID)) {
		hitsfileWriter.addReadMapRecord(readID);
		while (hitsfileReader.readReadMapRecordTranscriptID(transcriptID)) {
			hitsfileWriter.addTranscriptToReadMapRecord(transcriptID);
		}
		hitsfileWriter.writeReadMapRecord();
	}
}

#define NINSPECT 3

int main(int argc, char **argv) {
	vector<string> transcriptName;
	map<string, double> transcriptEffectiveLength;
	map<string, int> transcriptTrueLength;
	map<string, vector<string> > geneIsoforms;
	vector<vector<string> > identicalTranscripts;

	if (argc != 3) {
		printUsage(argv[0], cerr);
		exit(1);
	}
	if (strcmp(argv[1], "inspect") == 0) {
		HitsfileReader hitsfileReader(argv[2]);
		hitsfileReader.readHeader(&transcriptName, &transcriptEffectiveLength, &transcriptTrueLength, &geneIsoforms, &identicalTranscripts);
		cout << "Hits file: " << argv[2] << endl;
		cout << "Total number of transcripts: " << transcriptName.size() << endl;
		int count = 0;
		for(int i = 0; i < transcriptName.size(); i++) {
			count++; if (count > NINSPECT) { cout << "..." << endl; break; }
			cout << "@TranscriptMetaData\t" << transcriptName[i] << "\t" << transcriptEffectiveLength[ transcriptName[i] ]
           << "\t" << transcriptTrueLength[ transcriptName[i] ] << endl;
		}
		cout << "Total number of genes: " << geneIsoforms.size() << endl;
		count = 0;
		for(map<string, vector<string> >::iterator i = geneIsoforms.begin(); i != geneIsoforms.end(); i++) {
			count++; if (count > NINSPECT) { cout << "..." << endl; break; }
			cout << "@GeneIsoforms\t" << i->first;
			for(vector<string>::iterator j = (i->second).begin(); j != (i->second).end(); j++)
				cout  << "\t" << (*j);
			cout << endl;
		}
		cout << "Total number of identical transcript sets: " << identicalTranscripts.size() << endl;
		count = 0;
		for(vector<vector<string> >::iterator i = identicalTranscripts.begin(); i != identicalTranscripts.end(); i++) {
			count++; if (count > NINSPECT) { cout << "..." << endl; break; }
			cout << "@IdenticalTranscripts";
			for(vector<string>::iterator j = i->begin(); j != i->end(); j++) {
				cout << "\t" << (*j);
			}
			cout << endl;
		}
		cout << "Reads:" << endl;
		count = 0;
		string readID, transcriptID;
		while(hitsfileReader.readReadMapRecordReadID(readID)) {
			count++;
			if (count == 4) { cout << "..." << endl; }
			if (count >= 4) { while (hitsfileReader.readReadMapRecordTranscriptID(transcriptID)) {}; continue; }
			cout << '>' << readID << endl;
			while (hitsfileReader.readReadMapRecordTranscriptID(transcriptID)) {
				cout << transcriptID << endl;
			}
		}
		cout << "Total number of reads: " << count << endl;
	} else if(strcmp(argv[1], "header") == 0 ) {
		HitsfileReader hitsfileReader(argv[2]);
		hitsfileReader.readHeader(&transcriptName, &transcriptEffectiveLength, &transcriptTrueLength, &geneIsoforms, &identicalTranscripts);
    for(int i = 0; i < transcriptName.size(); i++) {
      cout << "@TranscriptMetaData\t" << transcriptName[i] << "\t" << transcriptEffectiveLength[ transcriptName[i] ]
                 << "\t" << transcriptTrueLength[ transcriptName[i] ] << endl;
    }
    for(map<string, vector<string> >::iterator i = geneIsoforms.begin(); i != geneIsoforms.end(); i++) {
      cout << "@GeneIsoforms\t" << i->first;
      for(vector<string>::iterator j = (i->second).begin(); j != (i->second).end(); j++)
        cout  << "\t" << (*j);
      cout << endl;
    }
    for(vector<vector<string> >::iterator i = identicalTranscripts.begin(); i != identicalTranscripts.end(); i++) {
      cout << "@IdenticalTranscripts";
      for(vector<string>::iterator j = i->begin(); j != i->end(); j++) {
        cout << "\t" << (*j);
      }
      cout << endl;
    }
  } else if ((strcmp(argv[1], "t") == 0) || (strcmp(argv[1], "b") == 0)) {
		convert(argv[2], argv[1]);
	} else {
		printUsage(argv[0], cerr);
		exit(1);
	}
}
