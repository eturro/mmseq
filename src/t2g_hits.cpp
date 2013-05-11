/*    Copyright 2011 Ernest Turro
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
#include <vector>
#include <iterator>

using namespace std;

#define HITCHAR '>'
#define COMMENTCHAR '@'

void printUsage(const char *bin, ostream& out) {
  out << "Usage: t2g_hits hits_file > gene_hits_file" << endl
      << endl
      << "Mandatory arguments:"
      << endl
      << "  hits_file          hits file generated with `bam2hits`\n" 
      << endl;
}

int main(int argc, char **argv) {
  if(argc != 2) {
    printUsage("t2g_hits", cerr);
    exit(1);
  }

  string hits_file = string(argv[1]);

  // START helpers
  ifstream ifs;
  string str;
  vector<string>::iterator strit;
  char buf16384[16384];
  // END helpers

  // Open hits file for reading
  ifs.open(hits_file.c_str());
  if(!ifs.good()) {
    cerr << "Error reading hits file \"" << hits_file << "\".\n";
    exit(1);
  }


  // Read meta-data at top of hits file 
  // Expecting several lines of @TranscriptMetaData ID length 
  // followed by several lines of @GeneIsoforms GID TID1 TID2 ...
  // followed optionally by several lines of @IdenticalTranscripts TID1 TID2 ...
  // first line after header lines should contain >read_id and is saved in str 
  
  map<string, string> t2g; 
  {
    while(true) {
      ifs.getline(buf16384,16384);
      if(ifs.eof()) {
        cerr << "Error. Malformed hits file.\n";
        exit(1);
      }
      str=buf16384;
      istringstream iss(str);
      vector<string> tokens;
      copy(istream_iterator<string>(iss),
        istream_iterator<string>(),
        back_inserter<vector<string> >(tokens));
      if(tokens[0].compare("@GeneIsoforms")==0) {
        for(int i=2; i < tokens.size(); i++) {
          t2g[tokens[i]] = tokens[1];
        }
      } else if(tokens[0][0] == HITCHAR) {
        break;
      } else {
        if(tokens[0].compare("@TranscriptMetaData") !=0 && tokens[0].compare("@IdenticalTranscripts") != 0) {
          cerr << "Error. Malformed hits file.\n";
          exit(1);
        }
      }
    }
  }

  while(!ifs.eof()) {
    if(str[0]==HITCHAR) {
      cout << str << endl;
      vector<string> comb;
      while(!(ifs >> str).eof() && str[0] != HITCHAR) {
        if(find(comb.begin(), comb.end(), t2g[str]) == comb.end()) comb.push_back(t2g[str]); 
      }
      sort(comb.begin(), comb.end());
      unique(comb.begin(), comb.end());
      for(strit=comb.begin(); strit != comb.end(); strit++) {
        cout << *strit << endl;
      }
    }
  }
  ifs.close(); ifs.clear();

  return 0;
}
