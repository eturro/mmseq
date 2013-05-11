#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string.h>
#include <sstream>

using namespace std;

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
  out << "Usage: offsetGTF ensembl_gtf_file indel_vcf_file" << endl
      << endl;
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
  char *gf = new char[arguments[0].size() + 1];
  char *vf = new char[arguments[1].size() + 1];
  strcpy(gf,arguments[0].c_str());
  strcpy(vf,arguments[1].c_str());

  cerr << "Parsing Ensembl GTF file...";

  map<string, vector<int> > gff_positions;
  map<string, bool> isstart; // start position => true, end position=>false;

  fstream ifs(gf);
  string str;

  while(true) {
    getline(ifs,str);
    if(ifs.eof()) { break; }
    if(str[0]=='#') continue;
    vector<string> tokens;
    tokenise(str, tokens, "\t");
    if(tokens.size() < 3) {
      cerr << "Error: GTF doesn't have enough fields" << endl;
      cerr << str<< endl;
      exit(1);
    }
    if(tokens[2] != "exon") {
      continue;
    }
    if(gff_positions.count(tokens[0]) > 0) {
      gff_positions[tokens[0]].push_back(atoi(tokens[3].c_str()));
      gff_positions[tokens[0]].push_back(atoi(tokens[4].c_str()));
    } else {
      gff_positions[tokens[0]]=vector<int>(1, atoi(tokens[3].c_str()));
      gff_positions[tokens[0]].push_back(atoi(tokens[4].c_str()));
    }
  }
  ifs.close();
  ifs.clear();

  cerr << "done.\n";

  cerr << "Sorting GTF positions...";

  std::vector<int>::iterator vit;
  for(map<string, vector<int> >::iterator mit = gff_positions.begin(); mit != gff_positions.end(); mit++) {
    sort((mit->second).begin(), (mit->second).end());
    vit = unique ((mit->second).begin(), (mit->second).end());
    (mit->second).resize( distance((mit->second).begin(), vit) );
  }
  cerr << "done.\n";

  map<string, int> gff_newleftpositions; // Here string is in the form chr:pos
  map<string, int> gff_newrightpositions; // Here string is in the form chr:pos

  cerr << "Parsing VCF file and mapping positions...\n";

  ifs.open(vf);

  int offset=0;
  string chrom="";
  string prevchrom="";
  int pos=1;
  int gff_pos= -1;
  int gff_ind = 0;
  vector<string> tokens;
  while(true) {
    getline(ifs,str);
    if(ifs.eof()) { 
      int curoff=offset;
      while(gff_ind < gff_positions[chrom].size()) {
        if ( tokens[4].size() < tokens[3].size() && pos  < gff_positions[chrom][gff_ind] &&
              pos + tokens[3].size() - tokens[4].size() >= gff_positions[chrom][gff_ind]) {

          ostringstream ss;
          ss << gff_positions[chrom][gff_ind];
          gff_newleftpositions[chrom + ":" + ss.str()] = pos + 1 +  tokens[3].size() - tokens[4].size() + offset;
          gff_newrightpositions[chrom + ":" + ss.str()] = pos + tokens[3].size() - tokens[4].size() + offset;
          gff_ind++;
        } else {
          ostringstream ss;
          ss << gff_positions[chrom][gff_ind];
          gff_newleftpositions[chrom + ":" + ss.str()] = gff_positions[chrom][gff_ind] + offset;
          gff_newrightpositions[chrom + ":" + ss.str()] = gff_newleftpositions[chrom + ":" + ss.str()];
          gff_ind++;
        }
      }
      break;
    }
    if(str[0]=='#') continue;
    tokens.clear();
    tokenise(str, tokens, "\t");
    if(tokens.size() < 4) {
      cerr << "Error: VCF doesn't have enough fields" << endl;
      cerr << str << endl;
      exit(1);
    }
    chrom = tokens[0];
    if(gff_positions.count(chrom)==0) {
      cerr << "Warning: chromosome " << chrom << " appears in the VCF but not in the GTF.\n";
      gff_ind=0;
      continue;
    }
    pos = atoi(tokens[1].c_str());
    if(prevchrom != chrom && prevchrom != "") { // We've changed chromosome on the VCF, finish up with the old chromosome on the GTF
      while(gff_ind < gff_positions[prevchrom].size()) {
        ostringstream ss;
        ss << gff_positions[prevchrom][gff_ind];
        gff_newleftpositions[prevchrom + ":" + ss.str()] = gff_positions[prevchrom][gff_ind] + offset;
        gff_newrightpositions[prevchrom + ":" + ss.str()] = gff_newleftpositions[prevchrom + ":" + ss.str()];
        gff_ind++;
      }
      offset=0;
      gff_ind=0;
    }

    prevchrom=chrom;
    if ( tokens[4].size() < tokens[3].size() && pos < gff_positions[chrom][gff_ind] &&
          pos + tokens[3].size() - tokens[4].size() >= gff_positions[chrom][gff_ind]) {

      ostringstream ss;
      ss << gff_positions[chrom][gff_ind];
      gff_newleftpositions[chrom + ":" + ss.str()] = pos + 1 + offset;
      gff_newrightpositions[chrom + ":" + ss.str()] = pos + offset;
      gff_ind++;
      offset += tokens[4].size() - tokens[3].size();
      continue;
    }
    if( pos < gff_positions[chrom][gff_ind]) {
      offset += tokens[4].size() - tokens[3].size();
      continue;
    }
    while(gff_ind < gff_positions[chrom].size() && gff_positions[chrom][gff_ind] <= pos) {
      ostringstream ss;
      ss << gff_positions[chrom][gff_ind];
      gff_newleftpositions[chrom + ":" + ss.str()] = gff_positions[chrom][gff_ind] + offset;
      gff_newrightpositions[chrom + ":" + ss.str()] = gff_newleftpositions[chrom + ":" + ss.str()];
      gff_ind++;
    }
    offset += tokens[4].size() - tokens[3].size();
  }
  ifs.close();
  ifs.clear();

  cerr << "done.\n";

  cerr << "Outputting new Ensembl GTF file (only exon records)...";

  ifs.open(gf);

  while(true) {
    getline(ifs,str);
    if(ifs.eof()) { break; }
    if(str[0]=='#') {
      cout << str << endl;
      continue;
    }
    tokens.clear();
    tokenise(str, tokens, "\t");
    if(tokens.size() < 3) {
      cerr << "Error: GTF doesn't have enough fields" << endl;
      cerr << str<< endl;
      exit(1);
    }
    if(tokens[2] != "exon") {
      continue;
    }
    if(gff_newleftpositions.count(tokens[0] + ":" + tokens[3])>0 && gff_newrightpositions.count(tokens[0] + ":" + tokens[4])>0 && gff_newrightpositions[tokens[0] + ":" + tokens[4]] < gff_newleftpositions[tokens[0] + ":" + tokens[3]]) {
        cerr << "Skip " << tokens[0] << ":" << tokens[3] << "-"<<tokens[4] << "\n";
      continue;
    }
    for(int i=0; i < 3; i++) { 
      cout << tokens[i] << "\t";
    }
    if(gff_newleftpositions.count(tokens[0] + ":" + tokens[3])==0) {
      cerr << "Warning: couldn't find a new position for " << tokens[0] << ":" << tokens[3] << endl;
      cout << tokens[3] << "\t";
    } else {
      cout << gff_newleftpositions[tokens[0] + ":" + tokens[3]] << "\t";
    }
    if(gff_newrightpositions.count(tokens[0] + ":" + tokens[4])==0) {
      cerr << "Warning: couldn't find a new position for " << tokens[0] << ":" << tokens[4] << endl;
      cout << tokens[4];
      if(5==tokens.size()) cout << "\n";
      else cout << "\t";
    } else {
      cout << gff_newrightpositions[tokens[0] + ":" + tokens[4]];
      if(5==tokens.size()) cout << "\n";
      else cout << "\t";
    }
    for(int i=5; i < tokens.size(); i++) { 
      cout << tokens[i];
      if(i < tokens.size()-1) {
        cout << "\t";
      } else {
        cout << endl;
      }
    }
  }
  ifs.close();
  ifs.clear();

  cerr << "done.\n";

  return 0;
}



