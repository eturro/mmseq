#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string.h>

extern "C" {
    #include "fasta.h"
}

using namespace std;

static const int FASTAWIDTH=60;

void printUsage(char *bin, ostream& out) {
  out << "Usage: extract_transcripts fasta_file ensembl_gtf_file" << endl
      << endl;
}

class GffRecord {
  string _chr;
  string _strand;
  vector<int> _starts; // 0-based, inclusive
  vector<int> _ends; // 0-based, inclusive
  map<string,string> _features;
  public:
  GffRecord() {_chr=""; };
  string chr() {return(_chr);}
  void chr(string chr) {_chr=chr; };
  void add_feature(string feature_type, string feature_name){
    _features[feature_type]=feature_name;
  }
  map<string,string> features() { return(_features); }
  void strand(string strand) {_strand=strand;}
  string strand() { return(_strand);}
  void clear() {
    _starts.clear(); _ends.clear(); _features.clear();
  }
  vector<int> * starts() { return(&_starts);}
  vector<int> * ends() { return(&_ends);}
  void add_start(int n) { _starts.push_back(n);}
  void add_end(int n) { _ends.push_back(n);}
};

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

void strip(string & str, string strip=" ") { 
  str.erase(0, str.find_first_not_of(strip));
  str.erase(str.find_last_not_of(strip)+1);
}

string revcomp(string str) {
  string res;
  for(int i=str.size()-1; i >= 0; i--) {
    if(str[i]=='A') res += 'T';
    if(str[i]=='T') res += 'A';
    if(str[i]=='G') res += 'C';
    if(str[i]=='C') res += 'G';
  }
  return(res);
}


int main(int argc, char** argv) {
  // parse arguments
  if(argc != 3) {
    cerr << "Error: mandatory arguments missing.\n";
    printUsage(argv[0], cerr);
    exit(1);
  }

  vector<string> arguments;
  for(int i=1; i < argc; i++) arguments.push_back(string(argv[i]));
  char *ff = new char[arguments[0].size() + 1];
  char *gf = new char[arguments[1].size() + 1];
  strcpy(ff,arguments[0].c_str());
  strcpy(gf,arguments[1].c_str());


  cerr << "Parsing Ensembl GTF file...";

  fstream ifs(gf);
  char buf16384[16384];
  string str;


  map<string, vector<GffRecord> > gff_records;

  int line=0;
  GffRecord current_rec;
  while(true) {
    ifs.getline(buf16384,16384);
    if(ifs.eof()) { break; }
    str=buf16384;
    vector<string> tokens;
    tokenise(str, tokens, "\t");
    if(tokens[2].compare("exon")==0) {
      vector<string> tokens2;
      string transcript_id;
      string gene_id;
      tokenise(tokens[8], tokens2, ";");
      for(int i=0; i < tokens2.size(); i++) {
        strip(tokens2[i]);
        vector<string> tokens3;
        tokenise(tokens2[i], tokens3, " ");
        if(tokens3[0] == "transcript_id") {
          strip(tokens3[1], "\"");
          transcript_id = tokens3[1];
        } else if(tokens3[0]=="gene_id") {
          strip(tokens3[1], "\"");
          gene_id=tokens3[1];
        }
      }

      if(current_rec.features()["transcript_id"] != transcript_id) {
        if(current_rec.chr() != "") {
          gff_records[current_rec.chr()].push_back(current_rec);  
          current_rec.clear();
        }
        current_rec.chr(tokens[0]);
        current_rec.strand(tokens[6]);
        current_rec.add_feature("gene_id", gene_id);
        current_rec.add_feature("transcript_id", transcript_id);
      } 
      current_rec.add_start(atoi(tokens[3].c_str())-1);
      current_rec.add_end(atoi(tokens[4].c_str())-1);
    }
  }
  gff_records[current_rec.chr()].push_back(current_rec);
  ifs.close();

  cerr << "done.\n";

  FASTAFILE *ffp;
  char *seq;
  char *name;
  int L;
  ffp = OpenFASTA(ff);

  cerr << "Parsing DNA FASTA file..." << endl;
  while(ReadFASTA(ffp, &seq, &name, &L)) {
    str=name;
    vector<string> tokens;
    tokenise(str,tokens, " ");
    string chr=str.substr(0,str.find_first_of(" "));
    str=seq;
    int chr_start=0; // get offset from name (e.g. for Y chromosome, which starts in pos 2649521)
    int chr_end=str.size()-1;
    for(int i=0; i < tokens.size(); i++) {
      if(tokens[i].find("chromosome:") != string::npos) { // otherwise assume no offset
        vector<string> tokens2;
        tokenise(tokens[i], tokens2, ":");
        chr_start=atoi(tokens2[3].c_str())-1;
        chr_end=atoi(tokens2[4].c_str()) -1;
      }
    }
    cerr << "  Chromosome " << string(name).substr(0,string(name).find_first_of(" ")) 
        << " (start=" << chr_start << ", end=" <<  chr_end << ")\n";
    vector<GffRecord>::iterator it;
    if((gff_records.count(chr)==0)) {
      cerr << "  WARNING: no GTF records found for chromosome " << chr << endl;
      continue;
    }
    for(it = (gff_records[chr]).begin(); it != (gff_records[chr]).end(); it++) {
      sort(it->starts()->begin(), it->starts()->end());
      sort(it->ends()->begin(), it->ends()->end());
      if((*it->starts())[0] < chr_start || (*it->ends())[(*it->ends()).size()-1] > chr_end) {
        cerr << "\nWarning: skipping " << it->features()["transcript_id"]
             << " because it is not fully embedded in " << chr << "[" << chr_start << ", "
             << chr_end << "].\n";
        continue;
      }

      string stitch;
      for(int i=0; i < it->starts()->size(); i++) {
        stitch = stitch + str.substr((*it->starts())[i],  (*it->ends())[i] - (*it->starts())[i] + 1);
      }
      if(it->strand()=="-") {
        stitch = revcomp(stitch);
      }
      cout << ">" << it->features()["transcript_id"] << " gene:" << it->features()["gene_id"] << endl;
      int faline=0;
      while(faline < stitch.size()) {
        cout << stitch.substr(faline, FASTAWIDTH) << endl;
        faline += FASTAWIDTH;
      }

    }

    free(seq);
    free(name);
  }
  CloseFASTA(ffp);
  cerr << "Parsing DNA FASTA file...done.                                           \n";

  return 0;
}



