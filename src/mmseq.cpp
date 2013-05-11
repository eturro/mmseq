/*    Copyright 2010, 2011 Ernest Turro
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
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#ifdef _OPENMP
  #include <omp.h>
  #define OMP_GET_MAX_THREADS omp_get_max_threads()
  #define OMP_GET_THREAD_NUM omp_get_thread_num()
#else
  #define OMP_GET_MAX_THREADS 1
  #define OMP_GET_THREAD_NUM 0
#endif // _OPENMP

#include "sokal.hh"

#include "hitsio.hpp"

#include "uh.hh"

#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

#ifndef VERSION
  #define VERSION ?
#endif

#define HITCHAR '>'
#define COMMENTCHAR '@'

namespace ublas = boost::numeric::ublas;
using namespace std;
using ublas::prod;
using ublas::inner_prod;
using ublas::matrix_row;
using ublas::matrix_column;

typedef ublas::compressed_matrix<bool> boolMat;
typedef boolMat::iterator1 boolMatIt1;
typedef boolMat::iterator2 boolMatIt2;
typedef ublas::compressed_matrix<int> intMat;
typedef intMat::iterator1 intMatIt1;
typedef intMat::iterator2 intMatIt2;
typedef ublas::compressed_matrix<double> doubleMat;
typedef intMat::iterator1 doubleMatIt1;
typedef intMat::iterator2 doubleMatIt2;

bool secondRowGreater(matrix_row<boolMat> a, matrix_row<boolMat> b) {
  for(matrix_row<boolMat>::iterator ita = a.begin(), itb = b.begin(); ; ita++, itb++) {
    bool enda = ita==a.end();
    bool endb = itb==b.end();
    if(enda && endb) return false;
    if(enda && !endb) return true;
    if(endb && !enda) return false;
    
    if(ita.index() < itb.index()) return false;
    if(ita.index() > itb.index()) return true;
  }
  cerr << "Error in secondRowGreater()\n";
  exit(1);
}

boolMat M(0,0); // GLOBAL M matrix
bool compareRowsByInt(int i, int j) {
  return(secondRowGreater(matrix_row<boolMat>(M,j), matrix_row<boolMat>(M,i)));
}

int powerof2 (unsigned int x) {
  while (((x & 1) == 0) && x > 1) x >>= 1;
  return (x == 1);
}

vector<double> meanvar(vector<double> v, double floor=numeric_limits<double>::min(), double ceil=numeric_limits<double>::max()) {
  double sum=0;
  double sum2=0;
  for(int i=0;i < v.size(); i++) {
    if(v[i] < floor) v[i] = floor;
    if(v[i] > ceil) v[i] = ceil;
    sum += v[i];
    sum2 += v[i]*v[i];
  }
  vector<double> res;
  res.push_back(sum/v.size());
  res.push_back((sum2-sum*sum/v.size())/(v.size()-1.0));
  return(res);
}


bool equalRows(matrix_row<boolMat> a, matrix_row<boolMat> b) {
  for(matrix_row<boolMat>::iterator ita = a.begin(), itb = b.begin(); ; ita++, itb++) {
    bool enda = ita==a.end();
    bool endb = itb==b.end();
    if(enda && endb) return true;
    if(enda && !endb) return false;
    if(endb && !enda) return false;
    if(ita.index() != itb.index()) return false;
  }
  cerr << "Error in equalRows()\n";
  exit(1);
}

void printUsage(char *bin, ostream& out) {
  out << "Usage: mmseq [OPTIONS...] hits_file output_base" << endl
      << endl
      << "Mandatory arguments:"
      << endl
      << "  hits_file          hits file generated with `bam2hits`\n" 
      << "  output_base        base name for output files"
      << endl << endl
      << "Optional arguments:\n"
      << "  -alpha FLOAT       value of alpha in Gamma prior for mu (default: 0.1)" << endl
      << "  -beta FLOAT        value of beta in Gamma prior for mu (default: 0.1)" << endl
      << "  -max_em_iter INT   maximum number of EM iterations (default: 1000)" << endl
      << "  -epsilon FLOAT     minimum loglik ratio between successive EM iterations (default: 0.1)" << endl
      << "  -gibbs_iter INT    number of Gibbs iterations (default: 16384)" << endl
      << "  -gibbs_ss INT      subsampling interval for Gibbs output (default: gibbs_iter/1024)" << endl
      << "  -seed INT          seed for the PRNG in thread 0 (default: 1234)" << endl
      << "  -debug             output additional diagnostic files" << endl
      << "  -help              print this help message" << endl
      << "  -version           print the version" << endl
      << endl;
}

int main(int argc, char **argv) {
  // Set max_threads here because on Mac OS X omp_get_max_threads() erroneously returns 1 in parallel regions
  int max_threads=OMP_GET_MAX_THREADS;

  // DEFAULT PARAMETER VALUES
  double alpha=0.1; 
  double beta=0.1;

  int max_em_iter=1000;
  double epsilon=0.1;

  int gibbs_iter=16384;
  int trace_length=1024;
  int gibbs_ss=gibbs_iter/trace_length;

  int seed=1234;

  bool debug=false;

  // parse arguments
  vector<string> arguments;
  for(int i=1; i < argc; i++) arguments.push_back(string(argv[i]));

  while(true) {
    if(arguments.size() >0 && arguments[0]=="-alpha") {
      arguments.erase(arguments.begin());
      alpha=strtod(arguments[0].c_str(), NULL);
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-beta") {
      arguments.erase(arguments.begin());
      beta=strtod(arguments[0].c_str(), NULL);
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && arguments[0]=="-max_em_iter") {
      arguments.erase(arguments.begin());
      max_em_iter=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && arguments[0]=="-epsilon") {
      arguments.erase(arguments.begin());
      epsilon=strtod(arguments[0].c_str(), NULL);
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && arguments[0]=="-gibbs_iter") {
      arguments.erase(arguments.begin());
      gibbs_iter=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && arguments[0]=="-gibbs_ss") {
      arguments.erase(arguments.begin());
      gibbs_ss=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && arguments[0]=="-seed") {
      arguments.erase(arguments.begin());
      seed=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && arguments[0]=="-debug") {
      debug=true;
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && ( arguments[0]=="-h" || arguments[0]=="--help" || arguments[0]=="-help")) {
      cerr << "Calculate mmseq expression estimates.\n";
      printUsage(argv[0], cerr);
      exit(1);
    } else if(arguments.size() > 0 && ( arguments[0]=="-v" || arguments[0]=="--version" || arguments[0]=="-version")) {
      cerr << "mmseq-" << QUOTE(VERSION) << endl;
      exit(1);
    } else {
      if(arguments.size()==2) {
        break;
      } else {
        if(arguments.size()>0 && arguments[0][0]=='-') {
          cerr << "Error: unrecognised option " << arguments[0] << ".\n";
        } else {
          cerr << "Error: mandatory arguments missing.\n";
        }
        printUsage(argv[0], cerr);
        exit(1);
      }
    }
  }

  if( gibbs_iter % gibbs_ss != 0 ) {
    cerr << "Error: gibbs_iter must be divisible by gibbs_ss.\n";
    printUsage(argv[0],cerr);
    exit(1);
  }

  gibbs_ss=gibbs_iter/trace_length;

  if( gibbs_iter <= 0 || trace_length <= 0 ) {
    cerr << "Error: no. of iteratons or trace length <= 0. Possible integer overflow - is gibbs_iter too high?\n";
    printUsage(argv[0],cerr);
    exit(1);
  }

  if(!powerof2(trace_length)) {
    cerr << "Error: gibbs_iter/gibbs_ss must be a power of 2.\n";
    printUsage(argv[0],cerr);
    exit(1);
  }

  string hits_file = arguments[0];
  string output_base(arguments[1]);

  // START helpers
  ifstream ifs;
  ofstream ofs;
  boost::iostreams::filtering_stream<boost::iostreams::output> gofs;
  gofs.push(boost::iostreams::gzip_compressor());
  string str, str2;
  vector<string>::iterator strit;
  char buf16384[16384];
  // END helpers

  // Open hits file for reading
  HitsfileReader hitsfileReader(hits_file);

  cout << "Running mmseq with parameters:\n"
       << "  alpha:         " << alpha << endl
       << "  beta:          " << beta << endl
       << "  max_em_iter:   " << max_em_iter << endl
       << "  epsilon:       " << epsilon << endl
       << "  gibbs_iter:    " << gibbs_iter << endl
       << "  gibbs_ss:      " << gibbs_ss << endl
       << "  seed[0]:       " << seed << endl
       << "  debug:         " << debug << endl
       << "  threads:       " << max_threads << endl;


  // Read meta-data at top of hits file 
  // Expecting several lines of @TranscriptMetaData ID length 
  // followed by several lines of @GeneIsoforms GID TID1 TID2 ...
  // followed optionally by several lines of @IdenticalTranscripts TID1 TID2 ...
  // first line after header lines should contain >read_id and is saved in str 
  
  map<string, double> sidLen; 
  map<string, int> sidSeqLen; 
  vector<string> transcriptList;
//  map<string, bool> sidRemcand; 

  map<string, vector<string> > gene2transcripts;
  vector< vector<string> > identical_transcripts;

  hitsfileReader.readHeader(&transcriptList, &sidLen, &sidSeqLen, &gene2transcripts, &identical_transcripts);

  map<string, string> transcript2gene;
  map<string, int> gene2index;
  {
    int g=0;
    vector<string> transcriptListGI;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin(); git != gene2transcripts.end(); git++) {
      gene2index[git->first] = g;
      for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
        if(transcript2gene.count(*strit) > 0) {
          cerr << "Error: transcripts must be nested within genes in GeneIsoforms metadata.\n";
          exit(1);
        }
        transcriptListGI.push_back(*strit);
        transcript2gene[*strit] = git->first;
      }
      g++;
    }
    vector<string> transcriptList2 = transcriptList;
    vector<string> transcriptListGI2 = transcriptListGI;
    sort(transcriptList2.begin(), transcriptList2.end());
    sort(transcriptListGI2.begin(), transcriptListGI2.end());
    strit = unique (transcriptList2.begin(), transcriptList2.end());
    strit = unique (transcriptListGI2.begin(), transcriptListGI2.end());
    if(transcriptList2.size() != transcriptList.size()) {
      cerr << "Error: duplicate transcripts in @TranscriptMetaData entries.\n";
      exit(1);
    }
    if(transcriptListGI2.size() != transcriptListGI.size()) {
      cerr << "Error: duplicate transcripts in @GeneIsoforms entries.\n";
      exit(1);
    }
    for(strit = transcriptList.begin(); strit != transcriptList.end(); strit++) {
      if(transcript2gene.count(*strit) ==0) {
        cerr << "Error: " << *strit << " does not belong to a gene in the @GeneIsoforms header entries.\n";
        exit(1);
      }
    }
  }

  // START Initialise variables M,  m, n, k, mu, l
  map<string, int> sidIndex;
  map<int, string> indexSid;
  map<vector<int>, int > indexComb;
  vector<int> k; 

  int m=0; // counts distinct combinations encountered
  int n=0; // counts distinct transcripts encountered

  vector<int> doublehits;
  bool doublehit;

  int numbermappedreads=0;
  
  while(hitsfileReader.readReadMapRecordReadID(str)) { // read *next read id* for *current read* to *str*
    doublehit=false;
    // start reading mappings for the next read->list(transcripts) section
    //if(str[0]==HITCHAR) {
    //cout << "read_id: " << str << endl;
    numbermappedreads++;
    vector<int> comb;
    while(hitsfileReader.readReadMapRecordTranscriptID(str)) { // read *next transcript id* for *current read* to *str*
      if(sidIndex.count(str)==0) { sidIndex[str] = n; indexSid[n]=str;n++; doublehits.push_back(0); }
      if(find(comb.begin(), comb.end(), sidIndex[str]) == comb.end()) comb.push_back(sidIndex[str]); 
      else {
//        cout << "Double hit: " << str << " " << indexSid[sidIndex[str]] << endl;
        doublehits[sidIndex[str]]++;
        doublehit=true;
      }
      // TODO: think about what to do when a read hits the same transcript in more than one place.
    }
    sort(comb.begin(), comb.end());
    if(indexComb.count(comb) == 0) {
      cout << "Found " << n << " transcripts in " << m << " transcript combinations.\r";
      //	for(vector<int>::iterator intit=comb.begin(); intit != comb.end(); intit++) cout << *intit << " ";
      //	cout << endl;
      indexComb[comb] = m;
      m++;
      k.resize(m);
      /* Here I'd rather do M.resize(c,t,true), but this is not implemented in Boost yet */
      // Dynamic matrix
      int factor1 = 1; int factor2 = 1;
      if(M.size1() < m ) factor1=2;
      if(M.size2() < n ) factor2=2;
      if(factor1==2 || factor2==2 ) {
        boolMat temp(M);
        int s1=factor1==2?m*factor1: M.size1();
        int s2=factor2==2?n*factor2: M.size2();
        M.resize(s1,s2,false);
        //cout <<"Now: M.size1()=="<<M.size1()<<", M.size2() == " << M.size2() << endl << flush;
        for(boolMatIt1 i1 = temp.begin1(); i1 != temp.end1(); ++i1) {
          for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
            M(i2.index1(),i2.index2())= *i2;
          }
        }
      }
      for(vector<int>::iterator intit=comb.begin(); intit != comb.end(); intit++)
        M(m-1, *intit) = 1;
    }
    k[indexComb[comb]]++;
  }
  cout << endl;

/*  cout << "Checking if sufficient memory (" << sizeof(double)*(n + identical_transcripts.size() + gene2transcripts.size())* trace_length *1/1048576<< "MB) is available...";
  try {
    double *temp=new double[(n + identical_transcripts.size() + gene2transcripts.size())* trace_length];
    delete [] temp;
  } catch(...) {
    cerr << "Insufficient memory available. Aborting." << endl;
    exit(1);
  }
  cout << "done.\n";
*/

  // Resize M to get rid of excessive rows/columns from dynamic bit above
  boolMat temp(M);
  M.resize(m,n,false);
  for(boolMatIt1 i1 = temp.begin1(); i1 != temp.end1(); ++i1) {
    for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
      M(i2.index1(),i2.index2())= *i2;
    }
  }

  time_t before,after;

/*
  // Absent transcript removal step, based on M and removal candidate status

  cout << "Transposing M (this may take a while)..." << flush;
  time(&before);
  M = trans(M);
  time(&after); 
  cout << "done (" << difftime(after,before) << " seconds).\n";

  cout << "Removing absent transcripts..." << flush;

  // Sort rows
  {
    vector<int> Mord(n,0);
    for(int i=1; i < n; i++) Mord[i] = Mord[i-1]+1;
    sort(Mord.begin(), Mord.end(), compareRowsByInt);
    map<int, string> indexSid_ = indexSid;
    for(int t=0; t < n; t++) indexSid[Mord[t]] = indexSid_[t];
    boolMat newMat(n,m);
    for(int t=0; t < n; t++) {
      matrix_row<boolMat> mr(M, t);
      for(matrix_row<boolMat>::iterator it = mr.begin(); it != mr.end(); ++it) {
        newMat(Mord[t],it.index()) = 1;
      }
    }
    M = newMat;
  }

  // Remove removal candidate iff it has duplicates and not all its duplicates are also removal candidates
  vector<int> removal_indices;
  for(int t=0; t < n; t++) {
    if(sidRemcand[indexSid[t]]) {
      vector<int> dup_indices(t,1);
      int start=t,end=t;
      while(start > 0 && equalRows(matrix_row<boolMat>(M,t), matrix_row<boolMat>(M, start-1))) start--;
      while(end < n-2 && equalRows(matrix_row<boolMat>(M,t), matrix_row<boolMat>(M, end+1))) end++; 
      if(end-start > 1) {
        bool allremcands=true;
        vector<int> potential_indices;
        for(int v=start; v <= end; v++) {
          if(!sidRemcand[indexSid[v]]) {
            allremcands=false;
          } else {
            potential_indices.push_back(v);
          }
        }
        if(!allremcands) removal_indices.insert(removal_indices.end(), potential_indices.begin(), potential_indices.end());
      }
    }
  }
  sort(removal_indices.begin(), removal_indices.end());
  removal_indices.resize(unique(removal_indices.begin(), removal_indices.end()) - removal_indices.begin());
#ifdef DEBUG
  ofs.open("removedTranscripts.txt");
  for(vector<int>::iterator it=removal_indices.begin(); it != removal_indices.end(); it++) {
    ofs << indexSid[*it] << " " << *it << endl;
  }
  ofs.close(); ofs.clear();
#endif

//  ofs.open("Mcols_before.txt");
//  for(int t=0; t < n; t++) {
//    ofs << indexSid[t] << "(" << sidRemcand[indexSid[t]] << ","
//    << !(find(removal_indices.begin(),removal_indices.end(),t)==removal_indices.end())
//    << ") ";
//    matrix_row<boolMat> mr(M, t);
//    for (matrix_row<boolMat>::iterator it = mr.begin(); it != mr.end(); ++it) {
//      ofs << it.index() << " ";
//    }
//    ofs << endl;
//  }
//  ofs.close();ofs.clear();
//
  if(removal_indices.size()>0) {
    // Shift indexSid
    int offset=0;
    for(int i=0; i < n; i++) {
      if(i==removal_indices[offset]) offset++;
      else {
        indexSid[i-offset]=indexSid[i];
        doublehits[i-offset]=doublehits[i];
      }
    }
    for(int i=n-removal_indices.size(); i <n; i++) indexSid.erase(i);

    // Update sidIndex just in case we need it
    sidIndex.clear();
    for(int i=0; i < indexSid.size(); i++) sidIndex[indexSid[i]] = i;
     
    // Reconstruct Mt without removed transcripts
    boolMat newMat(n-removal_indices.size(),m);
    int remind=0;
    int currow=0;
    for(int t=0; t < n; t++) {
      if(remind == removal_indices.size() || t != removal_indices[remind]) {
        matrix_row<boolMat> mr(M, t);
        for (matrix_row<boolMat>::iterator it = mr.begin(); it != mr.end(); ++it) {
          newMat(currow,it.index()) = 1;
        }
        currow++;
      } else {
        remind++;
      }
    }
    M = newMat;
    n=n-removal_indices.size();
  }

  cout << "done\n";
*/  
  cout << "Transposing M (this may take a while)..." << flush;
  time(&before);
  boolMat Mt = trans(M);
  time(&after);
  cout << "done (" << difftime(after,before) << " seconds).\n";

  vector<int> Mrowsum(m,0);
  vector<int> Mcolsum(n,0);
  for(boolMatIt1 i1 = M.begin1(); i1 != M.end1(); ++i1) {
    for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
      Mrowsum[i2.index1()]++;
      Mcolsum[i2.index2()]++;
    }
  }

  ///////////////
  
  vector<double> l(n);
  for(unsigned t=0; t < n; t++)  {
    if(indexSid.count(t)==0) {
      cerr << "Error: transcript with index " << indexSid[t] << " has no name.\n";
      exit(1);
    }
    if(sidLen.count(indexSid[t])==0) {
      cerr << "Error: transcript '" << indexSid[t] << "' has no length.\n";
      exit(1);
    }
    l[t] = (double)sidLen[indexSid[t]] * (double)numbermappedreads/1000000000.0; // apply this constant to get RPKM
    if(l[t] <= 0) { 
      cerr << "Error: transcript '" << indexSid[t] << "' has a length of zero.\n";
      exit(1);
    }
  }

  vector<vector<int> > counts_shared; // counts shared with different sized sets of transcripts

  for(int t=0; t < n; t++){
    vector<int> counts(100,0); // bin together sets >= 100
    counts_shared.push_back(counts);
  }
  // Set starting mu's
  ublas::vector<double> mu(n);
  for(int t=0; t < n; t++) mu[t] = 0.0;
  #pragma omp parallel
  {
    boolMatIt2 i1 = M.begin2();
    bool firstloop=true;
    int row;
    #pragma omp for schedule(static) 
    for(int t=0; t < n; t++){
      if(firstloop) {
        advance(i1,t);
        firstloop=false;
      }
      for (boolMatIt1 i2 = i1.begin(); i2 != i1.end(); ++i2) {
        row=i2.index1();
        mu[t] += (double)k[row]/Mrowsum[row];
        counts_shared[t][min(Mrowsum[row], 100)-1] += k[row];
      }
      mu[t] /= l[t];
      i1++;
    }
  }

  // END Initialise variables M, m, n, k, mu, l

  //// count unique hits for sets of transcripts (identical and genes)
  vector<int> identical_unique_hits;
  vector<int> gene_unique_hits;
  {
    // identical
    cerr << "Counting unique hits to sets of identical transcripts...";
    boolMat T(n, identical_transcripts.size());
    {
      int v=0;
      for(vector< vector<string> >::iterator idit = identical_transcripts.begin();
        idit != identical_transcripts.end(); idit++) {
        for(strit = idit->begin(); strit != idit->end(); strit++) {
          if(sidIndex.count(*strit)>0) {
            T(sidIndex[*strit], v) = 1;
          }
        }
        v++;
      }
    }
    identical_unique_hits = uh(M, T, k);
    // genes
    cerr << "done." << endl;
    cerr << "Counting unique hits to genes...";
    boolMat G(n, gene2transcripts.size());
    {
      int v=0;
      for(map<string, vector<string> >::iterator git = gene2transcripts.begin();
        git != gene2transcripts.end(); git++) {
        for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
          if(sidIndex.count(*strit)>0) {
            G(sidIndex[*strit], v) = 1;
          }
        }
        v++;
      }
    }
    gene_unique_hits = uh(M,G,k);
    cerr << "done." << endl;
  }

  ofs.open((output_base + ".k").c_str());
  for(int i=0;i < m; i++) ofs << k[i] <<  endl;
  ofs.close(); ofs.clear();

  ofs.open((output_base+".M").c_str());
  ofs << "#";
  for(int t=0;t < n; t++) ofs << "\t" << indexSid[t];
  ofs << endl;
  for (boolMatIt1 i1 = M.begin1(); i1 != M.end1(); ++i1) {
    for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
      ofs << i2.index1() << "\t" << i2.index2() << endl;
    }
  }
  ofs.close(); ofs.clear();
  
  if(debug) {
    // Output shared counts histograms
    ofs.open((output_base + ".sharedcounts").c_str());
    for(strit = transcriptList.begin(); strit != transcriptList.end(); strit++) {
      ofs << *strit << "\t";
      for(int i=0; i < 100; i++) {
        if(sidIndex.count(*strit) > 0 ) {
          ofs << counts_shared[sidIndex[*strit]][i] << "\t";
        } else {
          ofs << "0\t";
        }
      }
      ofs << endl;
    }
    ofs.close();ofs.clear();

    ofs.open((output_base + ".doublehits").c_str());
    for(int i=0;i < n; i++) ofs << doublehits[i] <<  endl;
    ofs.close(); ofs.clear();

    ofs.open((output_base+".Mt-nodups").c_str());
    vector<int> dupTs;
    ofstream ofs2((output_base + ".dupIDs").c_str());
    for (boolMatIt1 i1 = Mt.begin1(); i1 != Mt.end1(); ++i1) {
      if(i1.index1()>0 && equalRows(matrix_row<boolMat>(Mt, i1.index1()), matrix_row<boolMat>(Mt,i1.index1()-1))) {
        ofs2 << indexSid[i1.index1()] << endl;
      } else {
        for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
          ofs << i2.index1() << "\t" << i2.index2() << endl;
        }
      }
    }
    ofs.close(); ofs.clear();
    ofs2.close();ofs2.clear();
  } // end debug

  if(debug) { 
    ofs.open((output_base + ".trace_em.gz").c_str(), ios_base::out | ios_base::binary | ios_base::trunc);
    gofs.push(ofs);
    for(unsigned t=0; t < n; t++)
      gofs << indexSid[t] << " ";
    gofs << endl;
  }

  ublas::vector<double> mu_temp(n); // build up next iteration's mu before overwriting mu

  double loglik=0; // log likelihood

  #pragma omp parallel
  {
    #pragma omp for schedule(static) reduction(+:loglik) 
    for(int i=0; i < m; i++) {
      matrix_row<boolMat > mr (M, i);
      loglik += k[i]*log(inner_prod(mr, mu));
    }
    #pragma omp for schedule(static) reduction(+:loglik) 
    for(int t=0; t < n; t++) loglik -= mu[t]*l[t];
  }
  double loglik_temp=0; // next iteration's log likelihood
  double llr=epsilon+1; // log likelihood ratio between successive iterations

  int iter=0;
  cout.precision(5);
  cout.setf(ios::fixed,ios::floatfield);
  while(iter < max_em_iter && llr > epsilon) {
    cout << "EM iteration " << iter << flush;

    if(debug) {
      for(unsigned t=0; t < n; t++)
        gofs << mu[t] << " "; 
      gofs << endl;
    }

    loglik_temp=0;

    #pragma omp parallel 
    {
      int threadid=OMP_GET_THREAD_NUM;

      boolMatIt1 i1 = Mt.begin1();
      bool firstloop=true;
      int row;

      #pragma omp for schedule(static) 
      for(int t=0; t < n; t++) {
       if(firstloop) {
          advance(i1,t);
          firstloop=false;
        }
        double sum=0;
        for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
          row=i2.index2();
          matrix_row<boolMat > mr (M, row);
          sum += k[row]/inner_prod(mr, mu);
        }
        mu_temp[t] = mu[t]*sum/l[t];
        i1++;
      }
    
      #pragma omp for schedule(static) reduction(+:loglik_temp) 
      for(int i=0; i < m; i++) {
        matrix_row<boolMat > mr (M, i);
        loglik_temp += k[i]*log(inner_prod(mr, mu_temp));
      }
      #pragma omp for schedule(static) reduction(-:loglik_temp) 
      for(int t=0; t < n; t++) loglik_temp -= mu_temp[t]*l[t];
    }
    mu = mu_temp;
    llr = loglik_temp - loglik;
    loglik = loglik_temp;
        
    cout << ", log likelihood ratio: " << llr << "            \r";
   
    iter++;
  }

  if(debug) {
    gofs.pop();
    ofs.close(); ofs.clear();
  }

  cout <<  endl;

  ublas::vector<double> mu_em = mu;

  // gibbs trace file
  char *tgt_file = new char[(output_base + ".trace_gibbs.gz").length()+1];
  strcpy(tgt_file, (output_base + ".trace_gibbs.gz").c_str());
  ofs.open(tgt_file, ios_base::out | ios_base::binary | ios_base::trunc);
  gofs.push(ofs);
  double * mu_trace = new double[n * trace_length];

  for(unsigned t=0; t < n; t++)
    gofs << indexSid[t] << " "; 
  gofs << endl;

//  gsl_rng * rg = gsl_rng_alloc (gsl_rng_mt19937); //random number generator
  gsl_rng * rg[max_threads];
  for(int i=0; i < max_threads; i++) {
    rg[i] =  gsl_rng_alloc (gsl_rng_mt19937); 
    gsl_rng_set(rg[i], seed+i);
  }
  

  // set up X matrix with non-zero's in appropriate cells
  intMat X(m,n);
  for(boolMatIt1 i1 = M.begin1(); i1 != M.end1(); ++i1) {
    for(boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
      X(i2.index1(), i2.index2()) = 1;
    }
  }

  int *Xcolsum = new int[n];
  int * Xcolsums = new int[n*max_threads]; // each thread calculates partial Xcolsums
  for(iter=0; iter < gibbs_iter; iter++) {
    cout << "Gibbs iteration " << iter << "       \r";
    
    memset(Xcolsum,0,n*sizeof(int));
    memset(Xcolsums,0, n*max_threads*sizeof(int));

    #pragma omp parallel
    {
      int threadid=OMP_GET_THREAD_NUM;

      // Update X's
      intMatIt1 i1 = X.begin1();
      bool firstloop=true;
      #pragma omp for schedule(static) 
      for(int i=0; i < m; i++) {
        if(firstloop) {
          advance(i1,i);
          firstloop=false;
        }

        unsigned int x[Mrowsum[i1.index1()]];
        double p[Mrowsum[i1.index1()]];
        int count=0;

        for(intMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
          p[count] = mu[i2.index2()];
          count++;
        }
        
        gsl_ran_multinomial(rg[threadid], Mrowsum[i1.index1()], k[i1.index1()], p, x);

        count=0;
        for(intMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
          X(i2.index1(), i2.index2()) = x[count];
        //  #pragma omp atomic
        //  myXcolsum[i2.index2()] += x[count];
          Xcolsums[i2.index2()+n*threadid] += x[count];
          count++;
        }
        i1++;
      }

      #pragma omp barrier
      
      // sum up threads' partial Xcolsums
      #pragma omp for schedule(static)
      for(int t=0; t < n ; t++) {
        for(int i=0; i < max_threads; i++) Xcolsum[t] += Xcolsums[t+i*n];
      }

      #pragma omp barrier
      // Now Xcolsum is ready

      // Update mu's
      #pragma omp for schedule(static)
      for(int t=0; t < n ; t++) {
        mu[t] = gsl_ran_gamma(rg[threadid], alpha + Xcolsum[t] , 1.0/(beta + l[t]));
      }
    }

    if(iter % gibbs_ss == 0) {
      for(unsigned t=0; t < n; t++) {
        gofs << mu[t] << " ";
        mu_trace[t*trace_length + iter/gibbs_ss] = mu[t];
      }
      gofs << endl;
    }
  }

  delete [] Xcolsum;
  delete [] Xcolsums;

  cout << endl;
  gofs.pop();
  ofs.close(); ofs.clear();

  cout << "Amalgamating transcripts and calculating summary statistics..." << flush;

  // We have to deal with three types of situations:
  // 1: All components have hits (add MCMC traces)
  // 2: No components have hits (simulate from Gammas and add up)
  // 3: Some but not all components have hits (use MCMC traces and supplement with simulated Gamma variates)


  // For sets of identical transcripts, cases 1 and 2 apply
  // Deal with case 1 here. Case 2 is dealt with when outputting the tables
  // add up identical transcripts and isoforms within genes
  double * mu_trace_identical = new double[identical_transcripts.size() * trace_length];

  for(int i=0; i < identical_transcripts.size() * trace_length; i++) mu_trace_identical[i]=0;
  {
    int t=0;
    for(vector< vector<string> >::iterator idit = identical_transcripts.begin();
      idit != identical_transcripts.end(); idit++) {
      for(strit = idit->begin(); strit != idit->end(); strit++) {
        if(sidIndex.find(*strit) != sidIndex.end()) {
          for(int i=0; i < trace_length; i++) {
            mu_trace_identical[t*trace_length + i] += mu_trace[sidIndex[*strit]*trace_length + i]; 
          }
        }
      }
      t++;
    }
  }

  // For isoforms within genes, cases 1, 2 and 3 can apply
  // Here, deal with cases 1, 2 and 3
  // Record transcript traces for cases 2 and 3 in mu_trace_simu for proportion calcs later
  map<string, vector<double> > mu_trace_simu;
  double * mu_trace_gene = new double[gene2transcripts.size() * trace_length];
  for(int i=0; i < gene2transcripts.size() * trace_length; i++) mu_trace_gene[i]=0;
  {
    int g=0;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin();
      git != gene2transcripts.end(); git++) {
      for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
        if(sidIndex.find(*strit) != sidIndex.end()) {
          for(int i=0; i < trace_length; i++) {
            mu_trace_gene[g*trace_length + i] += mu_trace[sidIndex[*strit]*trace_length + i]; 
          }
        } else { // No hits, so add through simulation
          vector<double> temp;
          for(int i=0; i < trace_length; i++) {
            temp.push_back(gsl_ran_gamma(rg[0], alpha, 1.0/(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0)));
            mu_trace_gene[g*trace_length + i] += temp.back();
          }
          mu_trace_simu[*strit] = temp; // record trace for proportion estimates
        }
      }
      g++;
    }
  }

  // calculate proportion of expression for transcripts within genes
  double * prop_trace = new double[n * trace_length];
  map<string, vector<double> > prop_trace_simu;
  for(int t=0; t < n; t++) prop_trace[t]=NAN;
  {
    int g=0;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin();
      git != gene2transcripts.end(); git++) {
      for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
        if(sidIndex.find(*strit) != sidIndex.end()) {
          for(int i=0; i < trace_length; i++) {
            prop_trace[sidIndex[*strit]*trace_length + i] =
              mu_trace[sidIndex[*strit]*trace_length + i]/mu_trace_gene[g*trace_length + i];
          }
        } else {
          vector<double> temp;
          for(int i=0; i < trace_length; i++) {
            temp.push_back(mu_trace_simu[*strit][i]/mu_trace_gene[g*trace_length + i]);
          }
          prop_trace_simu[*strit] = temp;
        }
      }
      g++;
    }
  }  

  // output the traces for transcripts with zero unique hits
  //char *tgh_file = new char[(output_base + ".h0.trace_gibbs").length()+1];
  //strcpy(tgh_file, (output_base + ".h0.trace_gibbs").c_str());
/*
  ofs.open(tgh_file);
  for(strit = transcriptList.begin(); strit != transcriptList.end(); strit++) {
    if(sidIndex.count(*strit) > 0 && counts_shared[sidIndex[*strit]][0]==0) {
      ofs << *strit << " ";
    }
  }
  ofs << endl;
  for(int i=0; i < trace_length; i++) {
    for(strit = transcriptList.begin(); strit != transcriptList.end(); strit++) {
      if(sidIndex.count(*strit) > 0 && counts_shared[sidIndex[*strit]][0]==0) {
        ofs << mu_trace[sidIndex[*strit]*trace_length + i] << " ";
      }
    }
    ofs << endl;
  }
  ofs.close(); ofs.clear();
*/
  // output the traces for transcripts, sum of identical transcripts and sum of transcripts belonging to the same gene
  // and proportion of each transcript in each gene
  char *tgi_file = new char[(output_base + ".identical.trace_gibbs.gz").length()+1];
  strcpy(tgi_file, (output_base + ".identical.trace_gibbs.gz").c_str());
  char *tgg_file = new char[(output_base + ".gene.trace_gibbs.gz").length()+1];
  strcpy(tgg_file, (output_base + ".gene.trace_gibbs.gz").c_str());
  char *tgp_file = new char[(output_base + ".prop.trace_gibbs.gz").length()+1];
  strcpy(tgp_file, (output_base + ".prop.trace_gibbs.gz").c_str());

  ofs.open(tgi_file, ios_base::out | ios_base::binary | ios_base::trunc);
  gofs.push(ofs);
  {
    int t=0;
    for(vector< vector<string> >::iterator idit = identical_transcripts.begin(); idit != identical_transcripts.end(); idit++) {
      if(isfinite(log(mu_trace_identical[t*trace_length])) != 0) {
        for(strit = idit->begin(); strit != idit->end(); strit++) {
          gofs << *strit;
          if(strit->compare(idit->back())!=0) gofs << "+";
        }
        gofs << " ";
      }
      t++;
    }
  }
  gofs << endl;
  for(int i=0; i < trace_length; i++) {
    int t=0;
    for(vector< vector<string> >::iterator idit = identical_transcripts.begin(); idit != identical_transcripts.end(); idit++) {
      if(isfinite(log(mu_trace_identical[t*trace_length])) != 0) {
        gofs << mu_trace_identical[t*trace_length + i] << " ";
      }
      t++;
    }
    gofs << endl;
  }
  gofs.pop();
  ofs.close(); ofs.clear();

  ofs.open(tgg_file, ios_base::out | ios_base::binary | ios_base::trunc);
  gofs.push(ofs);
  {
    int g=0;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin(); git != gene2transcripts.end(); git++) {
      if(isfinite(log(mu_trace_gene[g*trace_length])) != 0) {
        gofs << git->first  << " ";
      }
      g++;
    }
  }
  gofs << endl;
  for(int i=0; i < trace_length; i++) {
    int g=0;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin(); git != gene2transcripts.end(); git++) {
      if(isfinite(log(mu_trace_gene[g*trace_length])) != 0) {
        gofs << mu_trace_gene[g*trace_length + i] << " ";
      }
      g++;
    }
    gofs << endl;
  }
  gofs.pop();
  ofs.close(); ofs.clear();


  ofs.open(tgp_file, ios_base::out | ios_base::binary | ios_base::trunc);
  gofs.push(ofs);

  for(unsigned t=0; t < n; t++) gofs << indexSid[t] << " ";
  gofs << endl;

  for(int i=0; i < trace_length; i++) {
    for(int t=0; t < n; t++) {
      gofs << prop_trace[t*trace_length + i] << " ";
    }
    gofs << endl;
  }
  gofs.pop();
  ofs.close(); ofs.clear();

  // get means of logged posterior samples
  // Also do this for simulated transcript traces
  vector<double> meanmu(n,0);
  vector<double> meanmu_identical(identical_transcripts.size(), 0);
  vector<double> meanmu_gene(gene2transcripts.size(), 0);
//  map<string, double> meanmu_simu;

  // log mu traces to make them more normal before calculating standard errors
  for(int i=0; i < trace_length; i++) {
    for(int t=0; t < n; t++) {
      mu_trace[t*trace_length + i] = log(mu_trace[t*trace_length + i] );
      meanmu[t] +=  mu_trace[t*trace_length + i] ;
    }
    for(int t=0; t < identical_transcripts.size() ; t++) {
      mu_trace_identical[t*trace_length + i] = log(mu_trace_identical[t*trace_length + i] );
      meanmu_identical[t] += mu_trace_identical[t*trace_length + i];
    }
    for(int t=0; t < gene2transcripts.size(); t++){
      mu_trace_gene[t*trace_length + i] = log(mu_trace_gene[t*trace_length + i] );
      meanmu_gene[t] += mu_trace_gene[t*trace_length + i];
    }
//    double temp=0.0;
//    for(map<string, vector<double> >::iterator tit = prop_trace_simu.begin();
//        tit != prop_trace_simu.end(); tit++) {
//      if(meanmu_simu.count(tit->first)<1) meanmu_simu[tit->first] = 0.0;
//      mu_trace_simu[tit->first][i] = log(mu_trace_simu[tit->first][i]);
//      meanmu_simu[tit->first] += mu_trace_simu[tit->first][i];
//    }
  }

  for(int t=0; t < n; t++) meanmu[t] /= trace_length;
  for(int i=0; i < identical_transcripts.size() ; i++) meanmu_identical[i] /= trace_length; 
  for(int i=0; i < gene2transcripts.size(); i++) meanmu_gene[i]/= trace_length;
/*  for(map<string, vector<double> >::iterator tit = prop_trace_simu.begin();
      tit != prop_trace_simu.end(); tit++) {
      meanmu_simu[tit->first] /= trace_length;
  }*/

  // now means and traces are on log scale.

  // Calculate proportion summaries for transcripts with hits
  vector<double> meanprop(n,0);
  vector<double> meanprobitprop(n,0);
  vector<double> ssprobitprop(n,0);
  vector<double> sdprobitprop(n,0);
  for(int i=0; i < trace_length; i++) {
    double temp;
    for(int t=0; t < n; t++) {
      meanprop[t] += prop_trace[t*trace_length + i] ;
      if(gene2transcripts.count(transcript2gene[indexSid[t]]) == 0) {
        cerr << "Error:" << transcript2gene[indexSid[t]] 
             << " not a key in gene2transcripts[" << indexSid[t] << "] (" << t << ")" << endl;
      }
      if(gene2transcripts[transcript2gene[indexSid[t]]].size() > 1) {
        temp = gsl_cdf_ugaussian_Pinv(min(max(prop_trace[t*trace_length + i], 0.000000001), 0.999999999));
      } else {
        temp = numeric_limits<double>::infinity();
      }
      meanprobitprop[t] += temp;
      ssprobitprop[t] += temp*temp;
    }
  }
  for(int t=0; t < n; t++) meanprop[t] /= trace_length;

  for(int t=0; t < n; t++) {
    sdprobitprop[t] = sqrt((ssprobitprop[t] - meanprobitprop[t]*meanprobitprop[t]/trace_length)/(trace_length - 1.0));
  }
  for(int t=0; t < n; t++) {
    meanprobitprop[t] /= trace_length;
  }

  // Calculate proportion summaries for transcripts with no hits
  map<string, double> meanprop_simu;
  map<string, double> meanprobitprop_simu;
  map<string, double> ssprobitprop_simu;
  map<string, double> sdprobitprop_simu;

  for(map<string, vector<double> >::iterator tit = prop_trace_simu.begin();
      tit != prop_trace_simu.end(); tit++) {
    meanprop_simu[tit->first]=0.0;
    meanprobitprop_simu[tit->first]=0.0;
    ssprobitprop_simu[tit->first]=0.0;
    sdprobitprop_simu[tit->first]=0.0;
    double temp;
    for(int i=0; i < trace_length; i++) {
      meanprop_simu[tit->first] += prop_trace_simu[tit->first][i];
      if(gene2transcripts.count(transcript2gene[tit->first]) == 0) {
        cerr << "Error:" << transcript2gene[tit->first] << " not a key in gene2transcripts[" << tit->first << "]" << endl;
      } 
      if(gene2transcripts[transcript2gene[tit->first]].size() > 1) {
        temp = gsl_cdf_ugaussian_Pinv(min(max(prop_trace_simu[tit->first][i], 0.000000001), 0.999999999));
      } else {
        temp = numeric_limits<double>::infinity();
      }
    meanprobitprop_simu[tit->first] += temp;
    ssprobitprop_simu[tit->first] += temp*temp;
    }
  }
  for(map<string, vector<double> >::iterator tit = prop_trace_simu.begin();
    tit != prop_trace_simu.end(); tit++) {
    meanprop_simu[tit->first] /= trace_length;
  }
  for(map<string, vector<double> >::iterator tit = prop_trace_simu.begin();
    tit != prop_trace_simu.end(); tit++) {
    sdprobitprop_simu[tit->first] = sqrt((ssprobitprop_simu[tit->first] - meanprobitprop_simu[tit->first]*meanprobitprop_simu[tit->first]/trace_length)/(trace_length - 1.0));
  }
  for(map<string, vector<double> >::iterator tit = prop_trace_simu.begin();
    tit != prop_trace_simu.end(); tit++) {
    meanprobitprop_simu[tit->first] /= trace_length;
  }

  // calculate sd, mumcse and iact for transcripts with hits, identical sets and genes
  double mumcse[n];
  double iact[n];
  double sd[n];
  for(int t=0; t < n; t++) {
    double var=0;
    double tau=0;
    int m=0;
    int len=trace_length;
    if(sokal(&len,&mu_trace[t*trace_length],&var,&tau,&m)!=0) {
      mumcse[t] = trace_length;
      iact[t] = NAN;
    } else {
      mumcse[t] = sqrt(tau*var/trace_length);
      iact[t] = tau;
    }
    sd[t]=sqrt(var);
  }
  delete [] mu_trace;

  double mumcse_identical[identical_transcripts.size()];
  double iact_identical[identical_transcripts.size()];
  double sd_identical[identical_transcripts.size()];
  for(int t=0; t < identical_transcripts.size(); t++) {
    double var=0;
    double tau=0;
    int m=0;
    int len=trace_length;
    if(sokal(&len,&mu_trace_identical[t*trace_length],&var,&tau,&m)!=0) {
      mumcse_identical[t] = trace_length;
      iact_identical[t] = NAN;
    } else {
      mumcse_identical[t] = sqrt(tau*var/trace_length);
      iact_identical[t] = tau;
    }
    sd_identical[t]=sqrt(var);
  }
  delete [] mu_trace_identical;

  double mumcse_gene[gene2transcripts.size()];
  double iact_gene[gene2transcripts.size()];
  double sd_gene[gene2transcripts.size()];
  for(int g=0; g < gene2transcripts.size(); g++) {
    double var=0;
    double tau=0;
    int m=0;
    int len=trace_length;
    if(sokal(&len,&mu_trace_gene[g*trace_length],&var,&tau,&m)!=0) {
      mumcse_gene[g] = trace_length;
      iact_gene[g] = NAN;
    } else {
      mumcse_gene[g] = sqrt(tau*var/(double)trace_length);
      iact_gene[g] = tau;
    }
    sd_gene[g] = sqrt(var);
  }
  delete [] mu_trace_gene;


  // calculate sd, mumcse and iact for transcripts without hits
  map<string, double> meanprop_map;
  map<string, double> meanprobitprop_map;
  map<string, double> sdprobitprop_map;


  double digalpha= gsl_sf_psi(alpha);
  double sqrtpolygalpha= sqrt(gsl_sf_psi_n(1,alpha));

  // Calc gene-level weighted sum of transcript lengths
  vector<double> gene_lengths(gene2transcripts.size(), 0.0);
  {
    int g=0;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin(); git != gene2transcripts.end(); git++) {
      if(isfinite(meanmu_gene[g]) != 0) {
        double sum=0;
        for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
          if(sidIndex.find(*strit) != sidIndex.end()) {
            gene_lengths[g] += sidLen[*strit] * exp(meanmu[sidIndex[*strit]]);
            sum += exp(meanmu[sidIndex[*strit]]);
          } else {
            gene_lengths[g] += sidLen[*strit] * exp((digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0)));
            sum += exp((digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0)));
          }
        }
        gene_lengths[g] /= sum;
      }
      g++;
    }
  }


/*
  //// Calculate proportions for transcripts with zero hits
  map<string, double> meanprop_map;
  map<string, double> meanprobitprop_map;
  map<string, double> sdprobitprop_map;
  map<string, double> meangene;
  map<string, double> sdgene;
  map<string, double> lengene;

  {
    int g=0;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin();
      git != gene2transcripts.end(); git++) {
      if(!isfinite(meanmu_gene[g])) {
      // Calculate E and sd of log(mu_a+mu_b+...) by MCMC (can't be done analytically)
        vector<double> gtrace(trace_length, 0.0);
        map<string, vector <double> > ttrace;
        for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
          vector<double> vec(trace_length, 0.0);
          ttrace[*strit] = vec;
        }

        for(int i=0; i < trace_length; i++) {
          for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
            ttrace[*strit][i] = gsl_ran_gamma(rg[0], alpha, 1.0/(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0));
            gtrace[i] += ttrace[*strit][i];
          }
        }
        
        // Get weighted length
        double gefflen=0;
        double esum=0;
        for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
          gefflen += sidLen[*strit] * exp((digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0)));  
          esum +=  exp((digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0)));
        }
        gefflen /= esum;
        lengene[git->first]=gefflen;
        
        // Get props
        vector<double> ptrace(trace_length);
        vector<double> psum(trace_length, 0.0);
        for(strit = (git->second).begin(); strit != (git->second).end(); strit++) {
          for(int i=0; i < trace_length; i++) {
            ptrace[i] = ttrace[*strit][i]/gtrace[i];
          }
          vector<double> mv = meanvar(ptrace);
          meanprop_map[*strit] = mv[0];

          for(int i=0; i < trace_length; i++) {
            ptrace[i] = gsl_cdf_ugaussian_Pinv(ptrace[i]);
          }
          // Now on probit scale
          mv = meanvar(ptrace);
          meanprobitprop_map[*strit] = mv[0];
          sdprobitprop_map[*strit] = sqrt(mv[1]);
        }

        for(int i=0; i < trace_length; i++) {
          gtrace[i]=log(gtrace[i]);
        }

        vector<double> mv=meanvar(gtrace);
        meangene[git->first]=mv[0];
        sdgene[git->first]=sqrt(mv[1]);

      } // Otherwise gene has hits and if component transcripts have no hits, should've taken care of them above
      g++;
    }
  }
*/
  ofs.open((output_base + ".mmseq").c_str());
  ofs << "# Mapped fragments: " << numbermappedreads << endl;
  ofs << "feature_id\tlog_mu\tsd\tmcse\tiact\teffective_length\ttrue_length\tunique_hits\tmean_proportion\tmean_probit_proportion\tsd_probit_proportion\tlog_mu_em\tobserved\tntranscripts\n";
  for(strit = transcriptList.begin(); strit != transcriptList.end(); strit++) {
    if(sidIndex.count(*strit) > 0 ) {
      if(gene2transcripts.count(transcript2gene[*strit])==0) {
        cerr << "Error: can't find transcript2gene[" << *strit << "]" << endl;
      }
      ofs << *strit << "\t" << meanmu[sidIndex[*strit]] << "\t" 
          << sd[sidIndex[*strit]] << "\t"
          << mumcse[sidIndex[*strit]] << "\t" << iact[sidIndex[*strit]] << "\t" << sidLen[*strit]  << "\t"
          << sidSeqLen[*strit]  << "\t" 
          << counts_shared[sidIndex[*strit]][0] << "\t"
          << meanprop[sidIndex[*strit]] << "\t"
          << meanprobitprop[sidIndex[*strit]] << "\t"
          << sdprobitprop[sidIndex[*strit]] << "\t"
          << log(mu_em[sidIndex[*strit]]) << "\t"
          << "1" << "\t" << gene2transcripts[transcript2gene[*strit]].size() << endl;
    } else {
      ofs << *strit << "\t" << digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0) << "\t" << sqrtpolygalpha << "\t" 
          << "0" << "\t"
          << "1" << "\t" << sidLen[*strit] << "\t" << sidSeqLen[*strit] << "\t" << 0 << "\t" << meanprop_simu[*strit] << "\t" 
          << meanprobitprop_simu[*strit] << "\t" << sdprobitprop_simu[*strit] << "\t" << "NA" << "\t" << "0" << "\t"
          << gene2transcripts[transcript2gene[*strit]].size() << endl;
   /*   if(meanprop_map.count(*strit)>0) { // gene has hits but gene does
        ofs << *strit << "\t" << digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0) << "\t" << sqrtpolygalpha << "\t" << "0" << "\t"
            << "NA" << "\t" << sidLen[*strit] << "\t" << sidSeqLen[*strit] << "\t" << 0 << "\t" << meanprop_map[*strit] << "\t" 
            << meanprobitprop_map[*strit] << "\t" << sdprobitprop[*strit] << "\t" << "NA" << endl;
      } else { // neither transcript nor gene have hits
        ofs << *strit << "\t" << digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0) << "\t" << sqrtpolygalpha << "\t" << "0" << "\t"
            << "NA" << "\t" << sidLen[*strit] << "\t" << sidSeqLen[*strit] << "\t" << 0 << "\t" << meanprop_map[*strit] << "\t" 
            << meanprobitprop_map[*strit] << "\t" << sdprobitprop_map[*strit] << "\t" << "NA" << endl; */
    }
  }
  ofs.close();ofs.clear();
  
  ofs.open((output_base + ".identical.mmseq").c_str());
  ofs << "# Mapped fragments: " << numbermappedreads << endl;
  ofs << "feature_id\tlog_mu\tsd\tmcse\tiact\teffective_length\ttrue_length\tunique_hits\tobserved\tntranscripts\n";
  {
    int t=0;
    for(vector< vector<string> >::iterator idit = identical_transcripts.begin();
      idit != identical_transcripts.end(); idit++) {
      if(isfinite(meanmu_identical[t])) {
        for(strit = idit->begin(); strit != idit->end(); strit++) {
          ofs << *strit;
          if(strit->compare(idit->back())!=0) ofs << "+";
          else ofs << "\t" << meanmu_identical[t] << "\t"
                   << sd_identical[t] << "\t" << mumcse_identical[t]
                   << "\t" << iact_identical[t] << "\t" << sidLen[*(idit->begin())] << "\t"
                   << sidSeqLen[*(idit->begin())] << "\t" << identical_unique_hits[t] << "\t" << "1" << "\t" << idit->size() << endl;
        }
      } else {
        for(strit = idit->begin(); strit != idit->end(); strit++) {
          ofs << *strit;
          if(strit->compare(idit->back())!=0) ofs << "+";
          else ofs << "\t" << log(idit->size()) + digalpha - log(beta + sidLen[*strit] * (double)numbermappedreads/1000000000.0) << "\t" 
                   << sqrtpolygalpha << "\t"
                   << "0" << "\t" << "NA" << "\t" << sidLen[*(idit->begin())] << "\t"
                   << sidSeqLen[*(idit->begin())] << "\t" << 0 << "\t" << "0" << "\t" << idit->size() << endl;
        }
      }
      t++;
    }
  }
  ofs.close();ofs.clear();

  ofs.open((output_base + ".gene.mmseq").c_str());
  ofs << "# Mapped fragments: " << numbermappedreads << endl;
  ofs << "feature_id\tlog_mu\tsd\tmcse\tiact\teffective_length\ttrue_length\tunique_hits\tntranscripts\tobserved\n";
  {
    int g=0;
    for(map<string, vector<string> >::iterator git = gene2transcripts.begin();
      git != gene2transcripts.end(); git++) {
      bool obs=false;
      for(vector<string>::iterator strit = (git->second).begin(); strit != (git->second).end(); strit++) {
        if(sidIndex.count(*strit) > 0) {
          obs=true;
          break;
        }
      }
      if(obs) {
        ofs << git->first  << "\t" << meanmu_gene[g] << "\t" 
            << sd_gene[g] << "\t" << mumcse_gene[g] << "\t" << iact_gene[g] << "\t" << gene_lengths[g] <<  "\t"
            << "NA" << "\t" << gene_unique_hits[g] << "\t"
            << (git->second).size() << "\t" << "1" << endl;
      } else {
          ofs << git->first <<  "\t" << meanmu_gene[g] << "\t" << sd_gene[g] << "\t"
            << sd_gene[g]/sqrt(trace_length) << "\t" << 1 << "\t" << gene_lengths[g] << "\t" 
            << "NA" << "\t" << "0" << "\t" << (git->second).size() << "\t" << "0" << endl;
      }
      g++;
    }
  }
  ofs.close();ofs.clear();

  cout << "done." << endl;

  for(int i=0; i < max_threads; i++) gsl_rng_free (rg[i]);

  cout << "Output files: " << endl
       << "  " << output_base << ".mmseq" << endl
       << "  " << output_base << ".identical.mmseq" << endl
       << "  " << output_base << ".gene.mmseq" << endl;
  cout << "  " << output_base << ".M" << endl
       << "  " << output_base << ".k" << endl << endl;
//       << "  " << output_base << ".h0.trace_gibbs" << endl;
  cout << "  " << output_base << ".trace_gibbs.gz" << endl
       << "  " << output_base << ".identical.trace_gibbs.gz" << endl
       << "  " << output_base << ".gene.trace_gibbs.gz" << endl
       << "  " << output_base << ".prop.trace_gibbs.gz" << endl << endl;


  if(debug) {
    cout << endl
         << "  " << output_base << ".trace_em.gz" << endl
         << "  " << output_base << ".sharedcounts" << endl
         << "  " << output_base << ".Mt-nodups" << endl
         << "  " << output_base << ".doublehits" << endl
         << "  " << output_base << ".dupIDs" << endl;
  }

/*
cout << "Outputting information matrix.. ";

  ofs.open("infomat.txt");

  ublas::vector<double> gamma(m);
  for(int i=0; i < m; i++) {
    matrix_row<boolMat > mr (M, i);
    gamma[i] = inner_prod(mr, mu);
  }

  doubleMat J(n,n);

  for(int t1=0; t1 < n; t1++) {
    for(int t2=0; t2 < n; t2++) {
      double sum=0;
      for(int i=0; i < m; i++) {
        if(M(i,t1)==1 && M(i,t2)==1) {
          sum += k[i]/(gamma[i]*gamma[i]);
        }
      }
      J(t1,t2) = sum;
      ofs << J(t1,t2) << " ";
    }
    ofs << endl;
  }
         
  ofs.close();
  cout << "done.\n";
*/
  return 0;
}
