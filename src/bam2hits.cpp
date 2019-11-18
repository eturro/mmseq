/*    Copyright 2011-2012 Ernest Turro
**
**    This program is free software; you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation; either version 2 of the License, or
**    (at your option) any later version.
**
**    This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program; if not, write to the Free Software
**    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/

#include "htslib/sam.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <limits.h>
#include <stdlib.h>
#include <libgen.h>
#include <math.h>
#include <map>
#include <sys/stat.h> 
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_cdf.h>

#include "hitsio.hpp"

extern "C" {
  #include "fasta.h"
}
#ifdef _OPENMP
  #include <omp.h>
  #define OMP_GET_MAX_THREADS omp_get_max_threads()
  #define OMP_GET_THREAD_NUM omp_get_thread_num()
#else
  #define OMP_GET_MAX_THREADS 1
  #define OMP_GET_THREAD_NUM 0
#endif // _OPENMP

#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

#ifndef VERSION
  #define VERSION ?
#endif

#define TOPTHRES 2.5 // use top 2.5% of isizes to set default deviation threshold
#define PAIRSFORISIZEFILTER 2000000 // use 2M read pairs to set isize expected and threshold
#define WINDOWFORISIZEMODE 3 // smooth isize frequencies in windows between min isize and TOPTHRES (must be odd number)
#define MINDEVIATION 20 // deviation less than 20bp is surprisingly tight.. issue warning and cancel filter
#define HISTHEIGHT 80
#define NCENTRALBINS 20 // get 20 bins
#define SDSFORFILTER 1.96 // for central 95% of a Gaussian (must be consistent with TOPTHRES)
#define MINEFFLEN 0.1 // hard lower bound on possible transcript effective length

using namespace std;

class FastaRecord {
  string _transcript;
  string _gene;
  string _seq;
  double _effectiveLength;
  public:
  FastaRecord(string transcript, string gene, string seq);
  string transcript() { return _transcript; }
  string gene() { return _gene; }
  string seq() { return _seq; }
  void effectiveLength(double l) { _effectiveLength=l; }
  double effectiveLength() { return _effectiveLength; }
  string getItem(char what) { // 'g' for gene(), 's' for seq()
    switch(what) {
      case 'g':
        return(this->gene());
      case 's':
        return(this->seq());
      default: 
        cerr << "Error: FastaRecord::"<<what << " does not exist.\n";
        exit(1);
    }
  }
};

FastaRecord::FastaRecord(string transcript, string gene, string seq) {
  _transcript=transcript;
  _gene=gene;
  _seq=seq;
  _effectiveLength = seq.size();
}

inline bool compareFastasByGene(FastaRecord a, FastaRecord b) {
  return (a.gene()<b.gene());
}
inline bool compareFastasBySeq(FastaRecord a, FastaRecord b) {
  return (a.seq()<b.seq());
}

void outputMergeHeaders(vector<FastaRecord> &sortedfastas, const char what, const char *atname, bool moreThanOne, HitsfileWriter &hitsfileWriter) {
  //ostream& out) {
  vector<FastaRecord> curSet;
  curSet.push_back(sortedfastas[0]);
  for(int i=0; i < sortedfastas.size()-1; i++) {
    if(sortedfastas[i+1].getItem(what)==curSet[0].getItem(what)) {
      curSet.push_back(sortedfastas[i+1]);
      continue;
    }
    if((moreThanOne && curSet.size()>1) || !moreThanOne) {
      //out << atname << "\t";
      //if(what=='g') out << curSet[0].getItem(what) << "\t";
      if (what == 'g') { hitsfileWriter.addGeneIsoformRecord( curSet[0].getItem(what) ); }
      if (what == 's') { hitsfileWriter.addIdenticalTranscriptsRecord(); }
      for(int j=0; j < curSet.size()-1; j++) {
        //out << curSet[j].transcript() << "\t";
        if (what == 'g') { hitsfileWriter.addTranscriptToGeneIsoformRecord( curSet[j].transcript() ); }
        if (what == 's') { hitsfileWriter.addTranscriptToIdenticalTranscriptsRecord( curSet[j].transcript() ); }
      }
      //out << curSet[curSet.size()-1].transcript() << endl;
      if (what == 'g') { hitsfileWriter.addTranscriptToGeneIsoformRecord( curSet[curSet.size()-1].transcript() ); }
      if (what == 's') { hitsfileWriter.addTranscriptToIdenticalTranscriptsRecord( curSet[curSet.size()-1].transcript() ); }
    }
    curSet.clear();
    curSet.push_back(sortedfastas[i+1]);
  }
  if((moreThanOne && curSet.size()>1) || !moreThanOne) {
    //out << atname << "\t";
    //if(what=='g') out << sortedfastas[sortedfastas.size()-1].getItem(what) << "\t";
    if (what == 'g') { hitsfileWriter.addGeneIsoformRecord( sortedfastas[sortedfastas.size()-1].getItem(what) ); }
    if (what == 's') { hitsfileWriter.addIdenticalTranscriptsRecord(); }

  for(int j=0; j < curSet.size()-1; j++) {
      //out << curSet[j].transcript() << "\t";
      if (what == 'g') { hitsfileWriter.addTranscriptToGeneIsoformRecord( curSet[j].transcript() ); }
      if (what == 's') { hitsfileWriter.addTranscriptToIdenticalTranscriptsRecord( curSet[j].transcript() ); }
  }
    //out << curSet[curSet.size()-1].transcript() << endl;
    if (what == 'g') { hitsfileWriter.addTranscriptToGeneIsoformRecord( curSet[curSet.size()-1].transcript() ); }
    if (what == 's') { hitsfileWriter.addTranscriptToIdenticalTranscriptsRecord( curSet[curSet.size()-1].transcript() ); }
  }
}

class Alignment {
  string _read;
  string _transcript;
  int _isize;
  int _nm;
  int _pos;
  int _mpos;
  public:
  Alignment(string r, string t, int i, int n, int pos, int mpos) ;
  string read() { return _read; }
  string transcript() { return _transcript; }
  void transcript(string t) { _transcript=t; }
  int isize() { return _isize; }
  void isize(int n){ _isize=n; }
  int NM() { return _nm; }
  void NM(int n) { _nm=n; }
  int pos() { return _pos; }
  void pos(int p){ _pos=p; }
  int mpos() { return _mpos; }
  void mpos(int p){ _mpos=p; }
};

Alignment::Alignment(string r, string t, int i, int n, int pos, int mpos) {
  _read=r;
  _transcript=t;
  _isize=i;
  _nm=n;
  _pos=pos;
  _mpos=mpos;
}


vector<Alignment> bestStratumFilter(vector<Alignment> & alns) {
  int minNM=INT_MAX;
  for(vector<Alignment>::iterator it= alns.begin(); it !=  alns.end(); it++) {
    if(minNM>it->NM()) minNM=it->NM();
  }
  vector<Alignment> res;
  for(vector<Alignment>::iterator it= alns.begin(); it !=  alns.end(); it++) {
    if(it->NM()==minNM) res.push_back(*it);
  }

  return(res);
}

vector<Alignment> repetitiveTranscriptFilter(vector<Alignment> & alns, int expected) {
  map<string, int> t2dev;
  map<string, int> t2ind;
  vector<int> torem;
  int i=0;
  for(vector<Alignment>::iterator it= alns.begin(); it !=  alns.end(); ++it) {
    if(t2dev.count(it->transcript())==0) {
      t2dev[it->transcript()]=abs(it->isize()-expected);
      t2ind[it->transcript()]=i;
    } else {
      if(abs(it->isize()-expected) < t2dev[it->transcript()]) {
        torem.push_back(t2ind[it->transcript()]);
        t2dev[it->transcript()]=abs(it->isize()-expected);
        t2ind[it->transcript()]=i;
      } else {
        torem.push_back(i);
      }
    }
    i++;
  }
  if(torem.size()==0) return(alns);

  vector<Alignment> res;
  sort(torem.begin(), torem.end());
  i=0;
  for(int j=0; j < alns.size(); j++) {
    if(i < torem.size() && j==torem[i]) {
      i++;
    } else {
      res.push_back(alns[j]);
    }
  }
  return(res);
}

// Warning: assumes each alignment contains a different transcript (i.e. repetitiveTranscriptFilter has been called)
vector<Alignment> isizeDeviationFilter(vector<Alignment> & alns, int expected, int threshold) {
  int minDeviation=INT_MAX;
  for(vector<Alignment>::iterator it= alns.begin(); it !=  alns.end(); it++) {
    int deviation=abs(it->isize()-expected);
    if(minDeviation > deviation) minDeviation = deviation;
  }
  
  if(minDeviation > threshold) return(alns);
  vector<Alignment> res;
  for(vector<Alignment>::iterator it= alns.begin(); it !=  alns.end(); it++) {
    int deviation=abs(it->isize()-expected);
    if(deviation <= threshold) res.push_back(*it);
  }
  return(res);
}

bool compareTranscripts(Alignment a, Alignment b) {
  return(a.transcript() < b.transcript());
}
bool sameTranscript(Alignment a, Alignment b) {
  return(a.transcript() == b.transcript());
}

void outputHits(vector<Alignment> alns, HitsfileWriter &hitsfileWriter) {//ostream& out) {
  //out << ">" << alns[0].read() << endl;
  hitsfileWriter.addReadMapRecord(alns[0].read());
  for(vector<Alignment>::iterator it= alns.begin(); it !=  alns.end(); it++) {
    //out << it->transcript() /*<<  "\t" << it->isize() << "\t" <<  it->NM() */ << endl;
    hitsfileWriter.addTranscriptToReadMapRecord(it->transcript());
  }
  hitsfileWriter.writeReadMapRecord();
}


string stripHapSuffix(string s) {
  if(s.rfind("_A") == s.size()-2 || s.rfind("_B") == s.size()-2 ) {
    return(s.substr(0, s.size()-2));
  } else {
    return(s);
  }
}

void flushOutput(vector<Alignment> &currentTranscripts, bool useIsizeFilter, int & expected, int &threshold, 
        bool mergeHaps, int& stratarem, int &reptranrem, int &isizerem, int &alignmentsretained, int &readsretained, HitsfileWriter &hitsfileWriter) {
    //, ostream& out) {
  int before=currentTranscripts.size();
  currentTranscripts=bestStratumFilter(currentTranscripts);
  stratarem+=(before-currentTranscripts.size());
  before=currentTranscripts.size();
  currentTranscripts=repetitiveTranscriptFilter(currentTranscripts, expected);
  reptranrem+=(before-currentTranscripts.size());
  if(useIsizeFilter) {
    before=currentTranscripts.size();
    currentTranscripts=isizeDeviationFilter(currentTranscripts, expected, threshold);
    isizerem+=(before-currentTranscripts.size());
  }

  if(mergeHaps) {
    for(vector<Alignment>::iterator it=currentTranscripts.begin(); it != currentTranscripts.end(); it++) {
      it->transcript(stripHapSuffix(it->transcript()));
    }
  }

  sort(currentTranscripts.begin(), currentTranscripts.end(), compareTranscripts);
  currentTranscripts.erase(
    unique(currentTranscripts.begin(), currentTranscripts.end(), sameTranscript), currentTranscripts.end()
  );

  alignmentsretained+= currentTranscripts.size();
  readsretained+=1;
  outputHits(currentTranscripts, hitsfileWriter);
}


vector<string> extractElements(const char *str, boost::regex re, vector<int> backrefs) {
  boost::cmatch matches;
  if(! boost::regex_match(str, matches, re)) {
    cerr << "Error: entry did not match regular expression:\n";
    cerr << "\tentry: " << str << endl;
    cerr << "\tregex: " << re << endl;
    cerr << "\tHints: use bowtie's --fullref option on FASTAs not following the Ensembl convention" << endl;
    cerr << "\t       use bam2hits's -m tg_regexp t_ind g_ind option to match transcript and gene IDs" << endl;
    exit(1);
  }
  vector<string> res;
  for(vector<int>::iterator it=backrefs.begin(); it != backrefs.end(); it++) {
    res.push_back(string(matches[*it].first,matches[*it].second));
  }
  return(res);
}

void checkSortOrder(bam_hdr_t *bam_head, const char *sortorder, ostream &out) {
  boost::cmatch ma;
  if(boost::regex_match(bam_head->text, ma, boost::regex(".*^@HD.*SO:(.[^\n]*)\n.*") )) {
    if(string(ma[1].first, ma[1].second) != sortorder) {
      out << "\nWarning: BAM header sort order '" << string(ma[1].first, ma[1].second)
           << "', should be '" << sortorder << "'.\n";
    }
  }
}

bool isFile(const char *path) {
  struct stat s;
  if( stat(path,&s) == 0 ) {
    if( S_ISREG(s.st_mode) ) return true;
    else return false;
  }
  return false;
}

void tokenise(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void printUsage(char *bin, ostream& out) {
  out << "Usage: bam2hits [-v] [-t] [-m tg_regexp t_ind g_ind] [-i expected_isize sd_of_isizes] [-k] [-b] [-c alt_lengths] [-u alt_lengths] transcript_fasta sorted_bam > hits_file" << endl
      << "Generate hits file for use with `mmseq`." << endl
      << "" << endl
      << "Mandatory arguments:" << endl
      << "  transcript_fasta:    FASTA file used as reference to produce SAM file." << endl
      << "  sorted_bam:          BAM file sorted by read name containing alignments to reference transcripts." << endl
      << "" << endl
      << "Option -v:             print version and quit." << endl
      << "" << endl
      << "Option -t:             output hits file in plain text instead of compressed binary format." << endl
      << "" << endl
      << "Arguments to option -m for use with non-Ensembl transcript FASTA files:" << endl
      << "  tg_regexp:           regular expression matching FASTA entry names, where pairs of brackets" << endl
      << "                       are used to capture transcript and gene IDs. E.g.: \"(\\S+).*gene:(\\S+).*\"" << endl
      << "  t_ind:               index of bracket pair that captures the transcript ID. E.g.: 1." << endl
      << "  g_ind:               index of bracket pair that captures the gene ID. E.g.: 2." << endl
      << "" << endl
      << "Arguments to option -i for calculating effective transcript lengths and filtering very unlikely alignments:" << endl
      << "  expected_isize:      expected insert size (e.g. 180). Default: for paired-end reads, the mode obtained using the" << endl
      << "                       first 2M single-mapping proper pairs, smoothing the histogram in windows of 3;" << endl
      << "                       for single-end reads, 180." << endl
      << "  sd_of_isizes:        standard deviation of insert sizes (e.g. 25). Default: for paired-end reads, 1/" << SDSFORFILTER << " times" << endl
      << "                       the distance between the mode and the top " << TOPTHRES << "% insert sizes in the first 2M single-mapping proper pairs;" << endl
      << "                       for single-end reads, 25." << endl
      << "" << endl
      << "Option -k:             for paired-end reads, disable insert size deviation filter, which removes alignments with" << endl
      << "                       insert sizes beyond " << SDSFORFILTER << " sd of the mean if alternative alignments exist" << endl
      << "                       with insert sizes within " << SDSFORFILTER << " sd of the mean." << endl
      << "" << endl
      << "Option -b:             merge hits to _A and _B arbitrary-phase haplo-isoforms. This is usually done if the only" << endl
      << "                       reason for creating haplo-isoform sequences is to reduce allelic mapping biases." << endl
      << "" << endl
      << "Argument to option -c for adjusting transcript lengths for non-uniformity with 'mseq' R package (required):" << endl
      << "  alt_lengths:         name of file in which to save the adjusted lengths (can be reused in future with -u)" << endl
      << "" << endl
      << "Argument to option -u for correcting for non-uniformity:" << endl
      << "  alt_lengths:         name of file containing adjusted transcript lengths" << endl
      << "" << endl;
}


int main(int argc, char *argv[]) {
  // Set max_threads here because on Mac OS X omp_get_max_threads() erroneously returns 1 in parallel regions
  int max_threads=OMP_GET_MAX_THREADS;

  cerr.setf(ios_base::fixed, ios_base::floatfield);
  cerr.precision(2);
  // max # BAM lines to read to get counts to estimate non-uniformity adjustment
  const int MAX_BAM_COUNT_ITERS=10000000;
  // max # transcripts to write out counts for if -c set
  const int MAX_T_COUNTS=500;
  // # mseq coefficients before/after start pos
  const int MSEQ_BEFORE=20;
  const int MSEQ_AFTER=21;

  // set some defaults
  boost::regex re("(\\S+).*gene:(\\S+).*");
  vector<int> tgBackrefs;
  tgBackrefs.push_back(1);
  tgBackrefs.push_back(2);
  bool customRef=false;
  
  bool useIsizeFilter=true;
  bool estimateIsizeDistroParams=true;
  int expected=180;
  int sd=25;
  int threshold=sd*SDSFORFILTER;

  bool mergeHaps=false;

  bool mseqCorrection=false;
  string adjLenFile;

  int readLength= -1;

  // general purpose string
  string str;

  // parse arguments
  vector<string> arguments;
  for(int i=1; i < argc; i++) arguments.push_back(string(argv[i]));

  if(find(arguments.begin(), arguments.end(), string("-v")) != arguments.end()) {
    cerr << "bam2hits-" << QUOTE(VERSION) << endl;
    exit(1);
  }

  string hitsfileFormat = "";

  while(arguments.size() != 2) {
    if(arguments.size() >0 && arguments[0]=="-m") {
      customRef=true; arguments.erase(arguments.begin());
      re.assign(arguments[0]); arguments.erase(arguments.begin());
      tgBackrefs[0] = boost::lexical_cast<int>(arguments[0]); arguments.erase(arguments.begin());
      tgBackrefs[1] = boost::lexical_cast<int>(arguments[0]); arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-i") {
      estimateIsizeDistroParams=false; arguments.erase(arguments.begin());
      expected = boost::lexical_cast<int>(arguments[0]); arguments.erase(arguments.begin());
      sd = boost::lexical_cast<int>(arguments[0]); arguments.erase(arguments.begin());
      threshold=sd*SDSFORFILTER;
    } else if(arguments.size() > 0 && arguments[0]=="-k") {
      useIsizeFilter=false; arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-b") {
      mergeHaps=true; arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-c") {
      mseqCorrection=true; arguments.erase(arguments.begin());
      adjLenFile = arguments[0]; arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-u") {
      arguments.erase(arguments.begin());
      adjLenFile = arguments[0]; arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-t") {
      hitsfileFormat = "t"; arguments.erase(arguments.begin());
    } else {
      if(arguments.size() != 2) {
        printUsage(argv[0], cerr);
        exit(1);
      }
    }
  }

  if (arguments.size() != 2) {
    printUsage(argv[0], cerr);
    exit(1);
  }

  char *ff = new char[arguments[0].size() + 1];
  char *bf = new char[arguments[1].size() + 1];
  strcpy(ff,arguments[0].c_str());
  strcpy(bf,arguments[1].c_str());

  if(expected <=0 || sd <= 0) {
    cerr << "Error: expected_isize and sd_of_isizes must be > 0." << endl; exit(1);
  }
  if(!isFile(ff)) {
    cerr << "Error: " << ff << " is not a file." << endl; exit(1);
  }
  if(!isFile(bf)) {
    cerr << "Error: " << bf << " is not a file." << endl; exit(1);
  }
  if(!mseqCorrection && adjLenFile.size()>0 && !isFile(adjLenFile.c_str())) {
    cerr << "Error: " << adjLenFile << " is not a file." << endl; exit(1);
  }
  
  // Parse FASTA and output @headers
  cerr << "Parsing FASTA file...";
  vector<FastaRecord> fastas;
  map<string, int> t2ind;
  {
    FASTAFILE *ffp;
    char *seq;
    char *name;
    int L;
    vector<string> tg;
    ffp = OpenFASTA(ff);
    while(ReadFASTA(ffp, &seq, &name, &L)) {
      tg = extractElements(name, re, tgBackrefs);
      if(mergeHaps) {
        int _B = tg[0].rfind("_B");
        if(_B != tg[0].size()-2) {
          int _A = tg[0].rfind("_A");
          if(_A == tg[0].size()-2) {
            _A = tg[1].rfind("_A");
            if(_A != tg[1].size()-2) {
              cerr << "Error. _A transcript does not belong to _A gene.\n";
              exit(1);
            }
            tg[0]=stripHapSuffix(tg[0]);
            tg[1]=stripHapSuffix(tg[1]);
          }
          fastas.push_back(FastaRecord(tg[0], tg[1], string(seq)));
          t2ind[tg[0]] = fastas.size()-1;
        }
      } else {
        fastas.push_back(FastaRecord(tg[0], tg[1], string(seq)));
        t2ind[tg[0]] = fastas.size()-1;
      }
      free(seq);
      free(name);
    }
    CloseFASTA(ffp);
  }
  cerr << "done (" << fastas.size() << " transcripts).\n"; 

  // Check that all the SN names in the header are also in the FASTA
  samFile *bamIn;
  if ((bamIn = sam_open(bf, "r")) == 0) {
    fprintf(stderr, "Failed to open BAM file %s\n", bf);
    exit(1);  
  }

  bam_hdr_t *bam_head = sam_hdr_read(bamIn);
  
  vector<string> tokens;
  tokenise(string(bam_head->text), tokens, "\n");
  for(vector<string>::iterator it=tokens.begin(); it != tokens.end(); it++) {
    vector<string> tokens2;
    tokenise(*it, tokens2, "\t");
    if(tokens2[0]=="@SQ") {
      for(int i=1; i < tokens2.size(); i++) {
        if(tokens2[i].find("SN:")==0) {
          tokens2[i]=tokens2[i].substr(3);
          str = customRef ? extractElements(tokens2[i].c_str(), re, tgBackrefs)[0].c_str() : tokens2[i];
          if(mergeHaps) {
            str = stripHapSuffix(str);
          }
          if(t2ind.find(str)==t2ind.end()) {
            cerr << "Error: transcript '" << str
                 << "' appears in the BAM header but it is not present in the transcript FASTA file "
                 << "(a different FASTA file must have been used for alignment).\n";
            exit(1);
          }
        }
      }
    }
  }

  // Set read length using first alignment in BAM file
  // Also determine if this is paired-end data
  bool pairedend=false;
  {
    bam1_t *line = bam_init1();
    if(sam_read1(bamIn, bam_head, line)>=0) {
      if(readLength == -1) {
        readLength=line->core.l_qseq;
        cerr << "Getting read length from BAM...";
        cerr << readLength << "bp.\n";
      }
      if( (line->core.flag & 0x1) == 0x1) {
        cerr << "Reads are paired-end." << endl;
        pairedend=true;
      } else {
        cerr << "Reads are single-end." << endl;
        estimateIsizeDistroParams=false;
        useIsizeFilter=false;
      }
    } else {
      cerr << "Error: no alignments in BAM file.\n";
      exit(1);
    }
  }
  bam_hdr_destroy(bam_head);
  sam_close(bamIn);

  if(mseqCorrection) {// correct for non-uniform coverage

    cerr << "Parsing BAM file...";
    
    vector<vector< int> > rawCounts(fastas.size());
    vector<double> transcriptCounts(fastas.size());
    for(int i=0; i < fastas.size(); i++) {
      rawCounts[i].resize(fastas[i].seq().size(), 0);
    }

    if ((bamIn = sam_open(bf, "r")) == 0) {
      fprintf(stderr, "Failed to open BAM file %s\n", bf);
      exit(1);  
    }

    bam_hdr_t *bam_head = sam_hdr_read(bamIn);
    checkSortOrder(bam_head, "queryname", cerr);

    bam1_t *line = bam_init1();
    int j=0;
    while(sam_read1(bamIn, bam_head, line)>=0) {
      if( (line->core.flag & 0xC) == 0 ) {
        str = customRef ? extractElements(bam_head->target_name[line->core.tid], re, tgBackrefs)[0].c_str() : string(bam_head->target_name[line->core.tid]).c_str();
        if(mergeHaps) {
          str = stripHapSuffix(str);
        }
        if(t2ind.find(str.c_str())==t2ind.end()) {
          cerr << "Error: transcript '" << str.c_str() 
               << "' has alignments but it is not present in the transcript FASTA file "
               << "(a different FASTA file must have been used for alignment).\n";
          exit(1);
        }
        rawCounts[t2ind[str.c_str()]][line->core.pos]++;
        transcriptCounts[t2ind[str.c_str()]]++;
      }
      j++;
      if(j>MAX_BAM_COUNT_ITERS) break;
    }
    bam_hdr_destroy(bam_head);
    sam_close(bamIn);
    cerr << "done.\n";

    for(int i=0;i < fastas.size(); i++) transcriptCounts[i] /= (double)fastas[i].seq().size();
    vector<double> tc2 = transcriptCounts;
    sort(tc2.rbegin(), tc2.rend());
    double thres=tc2[MAX_T_COUNTS]; 
    vector<string> seenGenes;
    cerr << "Mean reads per nucleotide of " << MAX_T_COUNTS << "th transcript: " << thres << endl;
  
    char cd[L_tmpnam+1];
    if(tmpnam (cd) == NULL) {
      cerr << "Could not create temporary file name.\n";
      exit(1);
    }
    ofstream ofs(cd);
    ofs << "index\ttag\tseq\tcount\n";
    for(int i=0; i < fastas.size(); i++) {
      if(transcriptCounts[i]>=thres && find(seenGenes.begin(), seenGenes.end(),fastas[i].gene()) == seenGenes.end()) {
        seenGenes.push_back(fastas[i].gene());
        for(int j =0; j < fastas[i].seq().size(); j++) {
          int tag = (j <= MSEQ_BEFORE || j >= fastas[i].seq().size()-MSEQ_AFTER) ? -2 : 0;
          ofs << fastas[i].transcript() << "\t" << tag << "\t" << fastas[i].seq()[j] << "\t" << rawCounts[i][j] << endl;
        }
      }
    }
    ofs << flush;

    cerr << "Calculating 'mseq' coefficients...\n";
    char cf[L_tmpnam+1];
    if(tmpnam (cf) == NULL) {
      cerr << "Could not create temporary file name.\n";
      exit(1);
    }
    
    char runmseq[L_tmpnam+1];
    if(tmpnam (runmseq) == NULL) {
      cerr << "Could not create temporary file name.\n";
      exit(1);
    }
    
    ofstream ofs2(runmseq);
    ofs2 << "library(methods)\n"
	<< "library(mseq)\n"
	<< "# need --args infile outfile\n"
	<< "MSEQ_BEFORE=20\n"
	<< "MSEQ_AFTER=21\n"

	<< "args=commandArgs()\n"
	<< "if(args[length(args)-2] != \"--args\")\n"
	<< "  stop(\"mmseq.R requires infile and outfile arguments.\n\")\n"

	<< "infile=args[length(args)-1]\n"
	<< "outfile=args[length(args)]\n"

	<< "my_counts <- read.delim(infile, header=T)\n"

  << "if('N' %in% unique(my_counts$seq)) {\n"
  << "  write(\"Warning: observed Ns in counts file:\", stderr())\n"
  << "  t <- table(my_counts$seq)\n"
  << "  write(paste(c(\"\\t\",names(t)), collapse=\"\\t\"), stderr())\n"
  << "  write(paste(c(\"Counts\t\",t), collapse=\"\\t\"), stderr())\n"
  << "  cat(\"Randomly changing Ns to (A,T,C,G).\")\n"
  << "  my_counts$seq[my_counts$seq==\"N\"] <- c('A','T','C','G')[floor(runif(t[names(t)==\"N\"],1,5))]\n"
  << "}\n"

	<< "for(t in unique(my_counts$index)) {\n"
	<< "  if(sum(my_counts$count[my_counts$index==t & my_counts$tag == 0])==0) {\n"
	<< "    my_counts <- my_counts[my_counts$index != t,]\n"
	<< "  }\n"
	<< "}\n"

	<< "data <- expData(my_counts, MSEQ_BEFORE, MSEQ_AFTER)\n"
	<< "data.glm <- iterGlm(data)\n"

	<< "write.table(data.glm$coefficients, file=outfile, quote=FALSE, col.names=FALSE)";
      ofs2 << flush;

#ifdef DEBUG
    cerr << "DEBUG: saving runmseq.R, countFile.txt and coefFile.txt to current working directory\n";
    system((string("cp ") + runmseq + " " + "runmseq.R").c_str()); 
    system((string("cp ") + cd + " " + "countFile.txt").c_str()); 
    system((string("cp ") + cf + " " +  "coefFile.txt").c_str()); 
#endif

    if(system((string("Rscript --vanilla ") + runmseq + " --args " + cd + " " + cf 
            + " 1>&2").c_str()) != 0) {
      cerr << "Error running Rscript --vanilla " << runmseq << " --args " << cd << " " << cf 
            << " 1>&2" << endl;
      exit(1);
    }
    cerr << "done.\n";

    ofs.close();
    remove(cd);
    ofs2.close();
    remove(runmseq);

    cerr << "Adjusting lengths to account for non-uniformity and saving to " << adjLenFile << "...\n";
    cerr << "  (this make take a long time, but if non-uniformity biases are the same for several" << endl
         << "   BAMs, then there is no need to do this more than once (use -u " << adjLenFile << " next time))\n";

    map<int, map<char,double> > coefs;
    for(int i = -MSEQ_BEFORE; i < MSEQ_AFTER; i++) {
      coefs[i]['A'] = 0.0;
      coefs[i]['T'] = 0.0;
      coefs[i]['C'] = 0.0;
      coefs[i]['G'] = 0.0;
    }
    ifstream ifs(cf);
    double alphacoef;
    if(ifs.is_open()) {
      string str;
      boost::cmatch ma;
      if(ifs.good()) {
        ifs >> str;
        ifs >> alphacoef;
      }
      while(ifs.good()) {
        ifs >> str;
        if(boost::regex_match(str.c_str(),ma, boost::regex("^pM([0-9]+)([ATCG]).*"))) { // rear coefficient
          ifs >> coefs[-atoi(string(ma[1].first, ma[1].second).c_str())][string(ma[2].first, ma[2].second).c_str()[0]];
        } else if(boost::regex_match(str.c_str(), ma, boost::regex("^p([0-9]+)([ATCG]).*"))){ // foward coefficient
          ifs >> coefs[atoi(string(ma[1].first, ma[1].second).c_str())][string(ma[2].first, ma[2].second).c_str()[0]];
        } else {
          cerr << "Incorrectly formatted 'mseq' coefficients file " << cf << endl;
        }
      }
    } else {
      cerr << "Could not open 'mseq' coefficients file " << cf << endl;
      exit(1);
    }

    ifs.close();
    remove(cf);

    double adjlen;
    double sum;
    int f=0;
    int nf=fastas.size();
    int basethread = -1;
    #pragma omp parallel for schedule(static) private(sum,adjlen) firstprivate(f) shared(basethread)
    for(int i=0; i < nf; i++) {
      if(i==0) basethread=OMP_GET_THREAD_NUM;
      adjlen=MSEQ_BEFORE+MSEQ_AFTER-1;
      if(fastas[i].seq().length() <= adjlen) { // can't adjust, keep length as is
        fastas[i].effectiveLength(fastas[i].seq().length());
        continue;
      }
      for(int j=MSEQ_BEFORE + 1; j <= fastas[i].seq().length()-MSEQ_AFTER; j++) {
        sum=0;
        for(int k= -MSEQ_BEFORE; k < MSEQ_AFTER ; k++) {
          sum += coefs[k][fastas[i].seq()[j+k]];
        }
        adjlen+= exp(sum+alphacoef);
      }
      fastas[i].effectiveLength(max(adjlen-readLength+(pairedend ? 0 : 1), MINEFFLEN));
      f++;
      if(OMP_GET_THREAD_NUM==basethread) cerr << "Progress " << 100.0*f*max_threads/nf << "%\r";
    }

    ofs.open(adjLenFile.c_str()); 
    for(vector<FastaRecord>::iterator it= fastas.begin(); it != fastas.end(); it++) {
      ofs << it->transcript() << "\t" << it->effectiveLength() << endl;
    }
    ofs.close();

    cerr << "done.           \n";
  }
  if(!mseqCorrection && adjLenFile.size()>0) {
    cerr << "Reading adjusted lengths from " << adjLenFile << "...";
    ifstream ifs(adjLenFile.c_str());
    if(ifs.is_open()) {
      string str;
      double adj;
      while(ifs.good()) {
        ifs >> str;
        ifs >> adj;
        if(t2ind.count(str)==0) {
          if(t2ind.count(str+"_A") > 0) {
            str=str+"_A";
          } else if(t2ind.count(str+"_B") > 0) {
            str=str+"_B";
          } else {
            cerr << "Error: transcript names in adjusted transcript lengths file do not match names in FASTA file.\n";
            cerr << "       " <<  str << " not found in FASTA." << endl;
            exit(1);
          }
        }
        fastas[t2ind[str]].effectiveLength(adj);
      }
    } else {
      cerr << "Could not open adjusted transcript lengths file " << adjLenFile << endl;
      exit(1);
    }
    ifs.close();
    cerr << "done.\n";
  }

  bool bam_eof=false;
  bam1_t *line;
  int i=0;

  // Use up to first 2M properly aligned pairs with unambiguous insert sizes in BAM file to calculate insert size distribution
  vector <int> testIsizes;

  if(pairedend) {
    if ((bamIn = sam_open(bf, "r")) == 0) {
      fprintf(stderr, "Failed to open BAM file %s\n", bf);  
      exit(1);  
    }

    string readname_prev="";
    int isize_prev = -1;
    bool same_isize=true;

    if(estimateIsizeDistroParams) {
      cerr << "Estimating insert size distribution parameters ";
    } else {
      cerr << "Printing insert size histogram "; 
    }
    cerr << "using first " << PAIRSFORISIZEFILTER << " proper read pairs with unambiguous insert isizes...\n";
    // get through the header
    bam_head = sam_hdr_read(bamIn);

    line = bam_init1();

    while(!bam_eof && testIsizes.size() <= PAIRSFORISIZEFILTER) {
      if(i==0 && sam_read1(bamIn, bam_head, line)<0) break;

      // forget about insert size filtering if not paired-end
      if( (line->core.flag & 0x1) == 0 ) {
        cerr << "aborting, as this is a single-end BAM file.\n";
        useIsizeFilter=false;
        break;
      }
      // skip if the query sequence itself is unmapped or the mate is unmapped 
      if( (line->core.flag & 0xC) > 0 ) {
        if(sam_read1(bamIn, bam_head, line)<0) bam_eof=true;
        continue;
      }
      if(line->core.isize > 0) {
        if(readname_prev != "" && string(bam_get_qname(line)) != readname_prev && same_isize) {
          testIsizes.push_back(isize_prev);
          same_isize=true;
        }
        if(string(bam_get_qname(line)) == readname_prev && isize_prev != line->core.isize) {
          same_isize=false;
        }
        if(string(bam_get_qname(line)) != readname_prev && !same_isize) {
          same_isize=true;
        }

        isize_prev=line->core.isize;
        readname_prev=string(bam_get_qname(line));
      }
      if(sam_read1(bamIn, bam_head, line)<0) {
        bam_eof=true;
        break;
      }
    }

    bam_hdr_destroy(bam_head);
    sam_close(bamIn);
    sort(testIsizes.begin(), testIsizes.end());
    int estthreshold = testIsizes[floor(testIsizes.size() * (100.0-TOPTHRES) / 100.0)];
    vector<int> frequencies;
    vector<int> values;
    int total=0;
    // Not the most efficient algorithm but should be quick in practice
    for(int j=testIsizes[0]; j <= estthreshold; j++) {
      values.push_back(j);
      frequencies.push_back(count(testIsizes.begin(), testIsizes.end(), j));
      total+=frequencies[frequencies.size()-1];
    }

    int estexpected=0;
    int maxsmooth=0;
    for(int j=1; j < frequencies.size()-WINDOWFORISIZEMODE + 1; j++) {
      int smoothsum=0;
      for(int k=j - (WINDOWFORISIZEMODE-1)/2; k <= j + (WINDOWFORISIZEMODE-1)/2; k++) {
        smoothsum += frequencies[k];
      }
      if(maxsmooth < smoothsum) {
        maxsmooth=smoothsum;
        estexpected=values[j];
      }
    }
    estthreshold = estthreshold - estexpected;
    if(estimateIsizeDistroParams && useIsizeFilter) {
      expected=estexpected;
      threshold=estthreshold;
      sd = estthreshold/SDSFORFILTER;
    }

    // Print isize histogram
    vector<int> heights, labels;
    heights.push_back(0);
    int maxheight=0;
    int isizehistbins = 2*estthreshold/NCENTRALBINS;
    if(isizehistbins==0) {
      cerr << "Error: insert size distribution does not have a central mode - check alignment parameters.\n";
      exit(1);
    }
    int curdiv=testIsizes[0]/isizehistbins;
    labels.push_back(curdiv * isizehistbins);
    for(int i=0; i < testIsizes.size(); i++) {
      if(testIsizes[i]/isizehistbins == curdiv) {
        heights[heights.size()-1] += 1;
      } else {
        curdiv = testIsizes[i]/isizehistbins;
        heights.push_back(0);
        labels.push_back(curdiv*isizehistbins);
        if(heights[heights.size()-2] > maxheight) maxheight= heights[heights.size()-2];
      }
    }
    int minbin= -1, maxbin=0;
    for(int i=0; i < heights.size(); i++) {
      if((HISTHEIGHT * heights[i]/maxheight) > 0) {
        if(minbin == -1 ) minbin=i;
        if(maxbin <= i) maxbin=i;
      }
    }
    for(int i=minbin; i <= maxbin; i++) {
      cerr << right << setw(4) << labels[i] << "|";
      for(int j=0; j < HISTHEIGHT * heights[i]/maxheight;j++) cerr << "-";
      cerr << endl;
    }
  }
  testIsizes.clear();

  cerr << "Expected insert size: " << expected << "bp";
  if(!estimateIsizeDistroParams) cerr << " (provided by user)";
  cerr << endl;
  cerr << "Standard deviation: " << sd << "bp";
  if(!estimateIsizeDistroParams) cerr << " (provided by user)";
  cerr << endl;

  if(useIsizeFilter) {
    cerr << "Insert size deviation filter thresholds: " << expected << " +/- " << threshold << "bp" << endl;
  }


  if(threshold < MINDEVIATION) {
    cerr << "Warning: thresholds less than " << MINDEVIATION << "bp from expected insert size, which seems very low. "
         << "Check your insert size distribution.\n"
         << "         Continuing without using insert size deviation filter.\n";
    useIsizeFilter=false;
  } 

  if(!mseqCorrection && adjLenFile.size()==0) {
    cerr << "Adjusting transcript lengths for insert size distribution...";
    double adjlen;
    for(vector<FastaRecord>::iterator it = fastas.begin(); it != fastas.end(); it++) {
      adjlen=0;
      for(int fl=1; fl <= (it->seq()).size(); fl++) {
        adjlen += (gsl_cdf_gaussian_P(fl+.5-expected, sd) - gsl_cdf_gaussian_P(fl-.5- expected,sd))* ((it->seq()).size()- fl + 1);
      }
      it->effectiveLength(max(adjlen, MINEFFLEN));
    }
    cerr << "done.\n";
  }

  // initialize HitsfileWriter to output to screen
  HitsfileWriter hitsfileWriter(hitsfileFormat);

  for(vector<FastaRecord>::iterator it = fastas.begin(); it != fastas.end(); it++) {
    if(it->effectiveLength() <= 0) {
      cerr << "Error: transcript effective lengths must be positive." << endl;
      exit(1);
    }
    hitsfileWriter.addTranscriptMetaData(it->transcript(), it->effectiveLength(), (it->seq()).size());
    // cout << "@TranscriptMetaData\t" << it->transcript() << "\t" << max(it->effectiveLength()-readLength+(pairedend ? 0 : 1),1.0) << "\n";
  }

  sort(fastas.begin(),fastas.end(),compareFastasByGene);
  outputMergeHeaders(fastas, 'g', "@GeneIsoforms", false, hitsfileWriter);
  sort(fastas.begin(),fastas.end(),compareFastasBySeq);
  outputMergeHeaders(fastas, 's', "@IdenticalTranscripts", true, hitsfileWriter);
  hitsfileWriter.writeHeader();


  // Parse _pre-sorted_ BAM file and output hits
  cerr << "Parsing BAM file...";

  if ((bamIn = sam_open(bf, "r")) == 0) {
    fprintf(stderr, "Failed to open BAM file %s\n", bf);  
    exit(1);  
  }

  bam_head = sam_hdr_read(bamIn);

//  boost::cmatch matches;
  boost::cmatch ma;
  if(boost::regex_match(bam_head->text, ma, boost::regex(".*^@HD.*SO:(.[^\n]*)\n.*") )) {
    if(string(ma[1].first, ma[1].second) != "queryname") {
      cerr << "\nWarning: BAM header sort order '" << string(ma[1].first, ma[1].second)
           << "', should be 'queryname'.\n";
    }
  }

  line = bam_init1();
  i=0;
  int stratarem=0, reptranrem=0, isizerem=0, readsretained=0, alignmentsretained=0;
  string currentRead;
  bool foundMate=false;
  vector<Alignment> currentTranscripts;
  bam_eof=false;

  while(!bam_eof) {
    if(i==0 && sam_read1(bamIn, bam_head, line)<0) break;
    
    // skip if the query sequence itself is unmapped or the mate is unmapped 
    if( (line->core.flag & 0xC) > 0 ) {
      if(sam_read1(bamIn, bam_head, line)<0) bam_eof=true;
      continue;
    }

    currentRead = string(bam_get_qname(line));
//    cerr << "currentRead: " << currentRead << " " << string(bam_head->target_name[line->core.tid]) << endl;

    currentTranscripts.clear(); 
    str = customRef ? extractElements(bam_head->target_name[line->core.tid], re, tgBackrefs)[0] : string(bam_head->target_name[line->core.tid]);
  //  cerr << "str: " << str << endl;
    if(mergeHaps) {
      str = stripHapSuffix(str);
    }
    Alignment aln(currentRead, str, abs(line->core.isize), bam_aux_get(line,"NM") == NULL ? 0 : bam_aux2i(bam_aux_get(line,"NM")), line->core.pos, line->core.mpos);
    currentTranscripts.push_back(aln);

//    if(sam_read1(bamIn, bam_head, line)<0) {
//      bam_eof=true;
//      flushOutput(currentTranscripts, useIsizeFilter, expected, threshold, mergeHaps,
//                  stratarem, reptranrem, isizerem, alignmentsretained, readsretained, hitsfileWriter);//cout);
//      break;
//    }
    currentRead = string(bam_get_qname(line)); i++;
//    cerr << "currentRead2: " << currentRead << " " << string(bam_head->target_name[line->core.tid]) << endl;

    while(currentRead == currentTranscripts.back().read()) {
      // if we found a mate in paired-end data, merely increase NM for that alignment
      foundMate=false;

      if(pairedend) {
        // if the transcript name is the same and pos1=mpos2 and pos2=mpos1
        for(vector<Alignment>::iterator ait= currentTranscripts.begin(); ait != currentTranscripts.end(); ait++) {
          str = customRef ? extractElements(bam_head->target_name[line->core.tid], re, tgBackrefs)[0] : string(bam_head->target_name[line->core.tid]);
          if(mergeHaps) {
            str = stripHapSuffix(str);
          }
     //     cerr << ait->transcript() << "," << str << "; " << ait->NM() << "+" << bam_aux2i(bam_aux_get(line,"NM")) << endl;
          if(ait->transcript() == str &&           
            ait->pos() == line->core.mpos && ait->mpos() == line->core.pos) {
            int f=ait->NM() + (bam_aux_get(line,"NM") == NULL ? 0 : bam_aux2i(bam_aux_get(line,"NM")));
      //      cerr << "Found Mate. Setting to " << f << endl;
            ait->NM(f);
            foundMate=true;
       //     cerr << "Found Mate. Total for " << ait->transcript() << ": " << ait->NM() << "\n";
            break;
          }
        }
      }

      // otherwise add as separate alignments
      if(!foundMate) {
        str = customRef ? extractElements(bam_head->target_name[line->core.tid], re, tgBackrefs)[0] : string(bam_head->target_name[line->core.tid]);
        if(mergeHaps) {
          str = stripHapSuffix(str);
        }
        Alignment aln(currentRead, str, abs(line->core.isize), bam_aux_get(line,"NM") == NULL ? 0 : bam_aux2i(bam_aux_get(line,"NM")), line->core.pos, line->core.mpos);
        currentTranscripts.push_back(aln);
      }

      if(sam_read1(bamIn, bam_head, line)<0) {
        bam_eof=true; 
        flushOutput(currentTranscripts, useIsizeFilter, expected, threshold, mergeHaps,
                  stratarem, reptranrem, isizerem, alignmentsretained, readsretained, hitsfileWriter);//cout);
        break;
      }
      currentRead = string(bam_get_qname(line)); i++;
    }
    
    flushOutput(currentTranscripts, useIsizeFilter, expected, threshold, mergeHaps,
                stratarem, reptranrem, isizerem, alignmentsretained, readsretained, hitsfileWriter);  
    cout << flush;
    i++;
  }
  bam_hdr_destroy(bam_head);
  sam_close(bamIn);
  cerr << "done.\n"
       << "Alignments removed by mismatch stratum filter: " << stratarem << "." << endl
       << "Alignments removed by repetitive transcript filter: " << reptranrem << "." << endl
       << "Alignments removed by isize deviation filter: " << isizerem << "." << endl
       << "Alignments retained: " << alignmentsretained << "." << endl;
  if(pairedend) {
    cerr << "Read pairs retained: " << readsretained << "." << endl;
  } else {
    cerr << "Reads retained: " << readsretained << "." << endl;
  }

  delete [] ff;
  delete [] bf;

  cerr << "Now run `mmseq`." <<  endl;
  return 0;
}
