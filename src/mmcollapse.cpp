/*    Copyright 2012 Ernest Turro
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
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <armadillo>
#include "sokal.hh"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "uh.hh"
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

#define TRACELEN 1024
#define IACTTHRES 1.1
#define SDPENALTY 0.0

using ublas::prod;
using ublas::inner_prod;
using ublas::matrix_row;
using ublas::matrix_column;

using namespace std;
using namespace arma;

namespace ublas = boost::numeric::ublas;
typedef ublas::compressed_matrix<bool> boolMat;
typedef boolMat::iterator1 boolMatIt1;
typedef boolMat::iterator2 boolMatIt2;

bool strorder ( const pair<string,int>& l, const pair<string,int>& r){
  return l.first < r.first;
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

bool stringCompare( const string &left, const string &right ){
   for( string::const_iterator lit = left.begin(), rit = right.begin(); lit != left.end() && rit != right.end(); ++lit, ++rit )
      if( tolower( *lit ) < tolower( *rit ) )
         return true;
      else if( tolower( *lit ) > tolower( *rit ) )
         return false;
   if( left.size() < right.size() )
      return true;
   return false;
}

void printUsage(char *bin, ostream& out) {
      out << "Usage: mmcollapse [-thres FLOAT] basename1 [basename2...]" << endl;
      out << "       -thres FLOAT   stopping threshold as -percentile of the maximum (over transcripts)" << endl
          << "                      mean (over samples) correlation (default 97.5)" << endl;
}

// find 0-unique-hit transcripts with SD > max(SD of transcripts with 1 unique hit and IACT < thres)
// including sets of identical transcripts
void get_unidentifiable_transcripts(const string filename,
    vector<string> & candidates, vector<string> & zeros, vector<string> & toremove,
    map<string, map<string, double> > & zeros_feature2iact,
    map<string, map<string, double> > & zeros_feature2sd,
    map<string, map<string, double> > & zeros_feature2efflen,
    map<string, int> & numbermappedreads,
    vector<string> & all_features,
    map<string, bool> & isIdent) {
  vector<double> sd;
  vector<double> iact;
  vector<string> tokens, tokens2;
  ifstream ifs;
  string str;
  bool retallfeatures = false;
  if(all_features.size()==0) retallfeatures =true;

  ifs.open((filename + ".identical.mmseq").c_str());
  if(!ifs.good()) {
    cerr << "Error: cannot open " << filename << ".identical.mmseq" << endl;
    exit(1);
  }

  getline(ifs,str);
  // skip comments
  while(ifs.good() && str[0]=='#') {
    getline(ifs,str);
  }
  if(!ifs.good()) {
    cerr << "Error: truncated *identical.mmseq file?\n";
    exit(1);
  }

  tokens.clear();
  tokenise(str, tokens, "\t");
  int ic = -1; // feature ID column
  int ac = -1; // iact column
  int uh = -1; // unique_hits column
  int ob = -1; // observed column
  int si = -1; // sd column
  int lm = -1; // log_mu column
  int el = -1; // effective_length column
  for(int i=0; i < tokens.size(); i++) {
    if(tokens[i]=="feature_id") {
      ic=i;
    } else if(tokens[i]=="iact") {
      ac=i;
    } else if(tokens[i]=="unique_hits") {
      uh=i;
    } else if(tokens[i]=="observed") {
      ob=i;
    } else if(tokens[i]=="sd") {
      si=i;
    } else if(tokens[i]=="log_mu") {
      lm=i;
    } else if(tokens[i]=="effective_length") {
      el=i;
    }
  }
  if(ic == -1) {
    cerr << "Error: " << filename << ".identical.mmseq file does not contain \"feature_id\" column.\n";
    exit(1);
  }
  if(ac == -1) {
    cerr << "Error: " << filename << ".identical.mmseq file does not contain \"iact\" column.\n";
    exit(1);
  }
  if(uh == -1) {
    cerr << "Error: " << filename << ".identical.mmseq file does not contain \"unique_hits\" column.\n";
    exit(1);
  }
  if(ob == -1) {
    cerr << "Error: " << filename << ".identical.mmseq file does not contain \"observed\" column.\n";
    exit(1);
  }
  if(si == -1) {
    cerr << "Error: " << filename << ".identical.mmseq file does not contain \"sd\" column.\n";
    exit(1);
  }
  if(lm == -1) {
    cerr << "Error: " << filename << ".identical.mmseq file does not contain \"log_mu\" column.\n";
    exit(1);
  }
  if(el == -1) {
    cerr << "Error: " << filename << ".identical.mmseq file does not contain \"effective_length\" column.\n";
    exit(1);
  }


  while(!ifs.eof()) {
    getline(ifs,str);
    tokens.clear();
    tokenise(str, tokens, "\t");
    if(tokens.size()==0) break;
    candidates.push_back(tokens[ic]);
    tokens2.clear();
    tokenise(tokens[ic], tokens2, "+");
    for(int i=0; i < tokens2.size(); i++) isIdent[tokens2[i]]=1;
    if(retallfeatures) all_features.push_back(tokens[ic]);
    if(tokens[ob]=="0") {
      sd.push_back(math::inf());
      iact.push_back(math::inf());
      zeros.push_back(tokens[ic]);
      zeros_feature2iact[filename][tokens[ic]]=atof(tokens[ac].c_str());
      zeros_feature2sd[filename][tokens[ic]]=atof(tokens[si].c_str());
      zeros_feature2efflen[filename][tokens[ic]]=atof(tokens[el].c_str());
    } else {
      sd.push_back(atof(tokens[si].c_str()));
      iact.push_back(atof(tokens[ac].c_str()));
    }
  }
  ifs.close(); ifs.clear();

  ifs.open((filename + ".mmseq").c_str());
  if(!ifs.good()) {
    cerr << "Error: cannot open " << filename << ".mmseq" << endl;
    exit(1);
  }

  getline(ifs,str);
  // skip comments
  while(ifs.good() && str[0]=='#') {
    getline(ifs,str);
    if(str.find("Mapped fragments") != string::npos) {
      tokens.clear();
      tokenise(str, tokens, " ");
      numbermappedreads[filename] = atof(tokens.back().c_str());
    }
  }
  if(!ifs.good()) {
    cerr << "Error: truncated *mmseq file?\n";
    exit(1);
  }

  tokens.clear();
  tokenise(str, tokens, "\t");
  ic = -1; // feature ID column
  int uc = -1; // unique hits column
  ac = -1; // iact column
  ob = -1; // observed column
  si = -1; // sd column
  lm = -1; // log_mu column
  el = -1; // effective_length column
  for(int i=0; i < tokens.size(); i++) {
    if(tokens[i]=="unique_hits") {
      uc=i;
    } else if(tokens[i]=="feature_id") {
      ic=i; 
    } else if(tokens[i]=="iact") {
      ac=i;
    } else if(tokens[i]=="observed") {
      ob=i;
    } else if(tokens[i]=="sd") {
      si=i;
    } else if(tokens[i]=="log_mu") {
      lm=i;
    } else if(tokens[i]=="effective_length") {
      el=i;
    }
  }
  if(uc == -1) {
    cerr << "Error: " << filename << ".mmseq file does not contain \"feature_id\" column.\n";
    exit(1);
  }
  if(uc == -1) {
    cerr << "Error: " << filename << ".mmseq file does not contain \"unique_hits\" column.\n";
    exit(1);
  }
  if(ac == -1) {
    cerr << "Error: " << filename << ".mmseq file does not contain \"iact\" column.\n";
    exit(1);
  }
  if(ob == -1) {
    cerr << "Error: " << filename << ".mmseq file does not contain \"observed\" column.\n";
    exit(1);
  }
  if(si == -1) {
    cerr << "Error: " << filename << ".mmseq file does not contain \"sd\" column.\n";
    exit(1);
  }
  if(lm == -1) {
    cerr << "Error: " << filename << ".mmseq file does not contain \"log_mu\" column.\n";
    exit(1);
  }
  if(el == -1) {
    cerr << "Error: " << filename << ".mmseq file does not contain \"effective_length\" column.\n";
    exit(1);
  }

  double max_h1_sd=0;

  int v=0;
  while(!ifs.eof()) {
    getline(ifs,str);
    tokens.clear();
    tokenise(str, tokens, "\t");
    if(tokens.size()==0) break;
    if(isIdent.count(tokens[ic])==0) {
      if(retallfeatures) all_features.push_back(tokens[ic]);
    }
    if(tokens[uc]=="0" && isIdent.count(tokens[ic])==0) {
      candidates.push_back(tokens[ic]);
      if(tokens[ob]=="0") {
        sd.push_back(math::inf());
        iact.push_back(math::inf());
        zeros.push_back(tokens[ic]);
        zeros_feature2iact[filename][tokens[ic]]=atof(tokens[ac].c_str());
        zeros_feature2sd[filename][tokens[ic]]=atof(tokens[si].c_str());
        zeros_feature2efflen[filename][tokens[ic]]=atof(tokens[el].c_str());
      } else {
        sd.push_back(atof(tokens[si].c_str()));
        iact.push_back(atof(tokens[ac].c_str()));
      }
    } else if(tokens[uc]=="1") {
      if(atof(tokens[si].c_str()) > max_h1_sd && atof(tokens[ac].c_str()) < IACTTHRES) {
        max_h1_sd=atof(tokens[si].c_str()) ;
      }
      toremove.push_back(tokens[ic]);
    } else {
      toremove.push_back(tokens[ic]);
    }
    v++;
  }

  cerr << "SD thres=" << max_h1_sd << endl;

  for(int i= candidates.size()-1; i >=0; i--) {
    if(sd[i] < max_h1_sd || isIdent.count(candidates[i])>0 ) {
      toremove.push_back(candidates[i]);
      candidates.erase(candidates.begin()+i);
    }
  }
}

/* void get_identical_sets(string filename, vector< vector<string> > & ident, vector< double> & ident_mcse) {
  vector<string> tokens;
  vector<string> tokens2;
  ifstream ifs;
  string str;
  ifs.open((filename + ".identical.mmseq").c_str());
  if(!ifs.good()) {
    cerr << "Error: cannot open " << filename << ".identical.mmseq" << endl;
    exit(1);
  }

  getline(ifs,str);
  // skip comments
  while(ifs.good() && str[0]=='#') {
    getline(ifs,str);
  }
  if(!ifs.good()) {
    cerr << "Error: truncated *identical.mmseq file?\n";
    exit(1);
  }

  int ic = 0; // feature ID column
  int mc = -1; // mcse column
  tokens.clear();
  tokenise(str, tokens, "\t");
  for(int i=0; i < tokens.size(); i++) {
    if(tokens[i]=="mcse") {
      mc=i;
    }
  }
  if(mc == -1) {
    cerr << "Error: *.identical.mmseq file does not contain \"mcse\" column.\n";
    exit(1);
  }

  while(!ifs.eof()) {
    getline(ifs,str);
    tokens.clear();
    tokenise(str, tokens, "\t");
    ident_mcse.push_back(atof(tokens[mc].c_str()));
    if(tokens.size()==0) break;
    tokens2.clear();
    tokenise(tokens[ic], tokens2, "+");
    ident.push_back(tokens2);
  }
}*/

// Collapse sets of transcripts given by sorted! indexes, updating R, candidates and cand2ind
void collapse(cube &R, vector<string> & candidates, map<string, int> & cand2ind, vector<int> indexes) {
  // Update variances and covariances
  for(int s=0; s < R.n_slices; s++) {
    for(int t=0; t < candidates.size(); t++) {
      for(int j=1; j < indexes.size(); j++) {
        if(t !=indexes[0]) {
          R(indexes[0], t, s) += R(t,indexes[j], s);
          R(t, indexes[0], s) += R(t,indexes[j], s);
        } else {
          R(indexes[0],t,s) += R(indexes[j],indexes[j],s);
          R(indexes[0],t,s) += R(indexes[0],indexes[j],s);
          R(indexes[0],t,s) += R(indexes[j],indexes[0],s);
        }
      }
    }
    // Set duplicate rows and cols to 1's to avoid recolapsing
    for(int j=1; j < indexes.size(); j++) {
      //R.slice(s).col(indexes[j]) = R.slice(s).col(indexes[0]);
      //R.slice(s).row(indexes[j]) = R.slice(s).row(indexes[0]);
      for(int i=0; i < R.n_cols; i++) {
        R(indexes[j],i,s)=math::nan();
        R(i,indexes[j],s)=math::nan();
      }
    }
  }

  // Sort new name
  vector<string> tset;
  for(int i=0; i < indexes.size(); i++) {
    tset.push_back(candidates[indexes[i]]);
  }
  sort(tset.begin(), tset.end());
  string newname = tset[0];
  for(int i=1; i < tset.size(); i++) {
    newname += "*" + tset[i];
  }
  cand2ind[newname] = indexes[0];
  candidates[indexes[0]] = newname;

  for(int i=1; i < indexes.size(); i++) {
    cand2ind.erase(candidates[indexes[i]]);
    candidates[indexes[i]] = "NA";
  }
}

map<int,int> join_traces(mat &M, vector<string> &allnames, vector<string> candidates) {
  map<string, int> tmap;
  for(int t=0; t < allnames.size(); t++) tmap[allnames[t]]=t;
  sort(candidates.begin(), candidates.end());
  vector<string>::iterator strit = unique(candidates.begin(), candidates.end());
  candidates.resize(strit-candidates.begin());

  vector<string> tokens;
  vector<int> toshed;
  for(int c=0; c < candidates.size(); c++) {
    tokens.clear();
    tokenise(candidates[c], tokens, "*");
    vector<int> indexes;
    for(int t=0; t < tokens.size(); t++) {
      if(tmap.count(tokens[t]) < 1) {
        cerr << "No trace for " << tokens[t] << " and trying to collapse.\n";
        exit(1);
      }
      indexes.push_back(tmap[tokens[t]]);
    }
    sort(indexes.begin(),indexes.end());
    for(int t=1; t < indexes.size(); t++) {
      for(int i=0; i < TRACELEN; i++) {
        M(i, indexes[0]) += M(i, indexes[t]);
      }
    }
    for(int t=1; t < indexes.size(); t++) {
      toshed.push_back(indexes[t]);
    }
    allnames[indexes[0]] = candidates[c];
  }

  map<int,int> shed;
  for(int i=0; i < toshed.size(); i++) {
    shed[toshed[i]]=1;
  }
  return(shed);
}

// set mean corrs in V
void mean_corrs(const cube &R, const umat &S, mat &V, mat &W, vector<int> ts, 
  vector<int> bl = vector<int>(), double sdpenalty=0) {
  for(vector<int>::iterator it=ts.begin(); it != ts.end(); it++) {
    for(int v= 0; v < R.n_rows; v++) {
      colvec r = (colvec)R.subcube(*it,v,0,*it,v,R.n_slices-1);
      r = r/sqrt((colvec)R.subcube(*it,*it,0,*it,*it,R.n_slices-1));
      r = r/sqrt((colvec)R.subcube(v,v,0,v,v,R.n_slices-1));
      urowvec u(S.row(*it) % S.row(v));
      // where u is zero, make sure r is zero (might not be if there was /0)
      for(int i=0; i < u.size(); i++) {
        if(u(i)==0) {
          r(i)=0;
        }
      }
      V(*it,v) = V(v,*it) = conv_to<double>::from(
        u * r/(u*ones<vec>(R.n_slices)));
      if(R.n_slices>1) {
        W(*it,v) = W(v,*it) = sqrt(conv_to<double>::from((u*ones<vec>(R.n_slices)/(u*ones<vec>(R.n_slices)-1)) * (
          u * (r % r)/(u*ones<vec>(R.n_slices)) - V(*it,v)*V(*it,v))));
        if(!is_finite(W(*it,v))) {
          W(*it,v) = W(v,*it) = 0;
        }
      } else {
        W(*it,v) = W(v,*it) = 0.0;
      }
      V(*it,v)=V(v,*it) = V(*it,v)+sdpenalty*W(*it,v);
    }
  }
}


void get_corrs(cube &R, mat &M, const vector<string> & basenames, map<string, int> & cand2ind) {

  cerr << "Reading posterior traces and computing variance-covariance matrices using " << min(OMP_GET_MAX_THREADS, (int)basenames.size()) << " thread(s)" << endl
        << "this step requires up to " << setprecision(2) << M.n_cols*M.n_cols*7.450581e-09 << " GB per thread; re-run with a lower value for OMP_NUM_THREADS to use less memory but more CPU...\n";
  #pragma omp parallel for shared(R, M) schedule(static)
  for(int s=0; s < basenames.size(); s++) {
    vector<string> sources;
    vector<string> tokens;
    ifstream ifs;
    boost::iostreams::filtering_stream<boost::iostreams::input> gifs;
    gifs.push(boost::iostreams::gzip_decompressor());
    string str;
    sources.push_back(".");
    sources.push_back(".identical.");
    mat myM(M);
    #pragma omp critical
    {
    for(vector<string>::iterator strit=sources.begin(); strit < sources.end(); strit++) {
      cerr << "\t" << basenames[s] << *strit << "trace_gibbs.gz" << endl;
      ifs.open((basenames[s] + *strit + "trace_gibbs.gz").c_str(), ios_base::in | ios_base::binary);
      gifs.push(ifs);
      if(!ifs.good()) {
        cerr << "Error: couldn't open " << basenames[s] << *strit << "trace_gibbs.gz.\n";
        exit(1);
      }
      getline(gifs,str); // get IDs
      tokens.clear();
      tokenise(str,tokens," ");
      for(int i=0; i < TRACELEN; i++) {
        for(int t=0; t < tokens.size(); t++) {
          gifs >> str;
          if(cand2ind.count(tokens[t]) > 0 ) {
            if(str=="NA") myM(i,cand2ind[tokens[t]])=math::nan();
            myM(i,cand2ind[tokens[t]])=atof(str.c_str());
          }
        }
      }
      gifs.pop();
      ifs.close();ifs.clear();
    }
    }
    R.slice(s) = cov(myM);
    // Set nans to 0
    for(cube::slice_iterator it=R.begin_slice(s); it != R.end_slice(s); it++) {
      if(! is_finite(*it)) *it = 0;
    }
  }
}


int main(int argc, char **argv) {

  gsl_rng * rg = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(rg, 13837);

  bool debug=false;

  // DEFAULT PARAMETER VALUES
  double alpha=0.1;
  double beta=0.1;

  vector<string> tokens, tokens2, tokens3; // buffer
  vector<string> strvec; // buffer
  vector<string>::iterator strit, strit2, strit3;
  string str; // buffer
  ifstream ifs;
  boost::iostreams::filtering_stream<boost::iostreams::input> gifs;
  gifs.push(boost::iostreams::gzip_decompressor());
  ofstream ofs,ofs2;

  vector<string> basenames;

  double stoppingthreshold=0.975;

  // parse arguments
  for(int i=1; i < argc; i++) {
    if(string(argv[i])=="-thres" && argc> i+1) {
      stoppingthreshold=strtod(argv[i+1], NULL)/100.0;
      if(stoppingthreshold <= 0 || stoppingthreshold > 1) {
        printUsage(argv[0],cerr);
      }
      i++;
    } else if(string(argv[i]).find("-") == 0) {
      printUsage(argv[0],cerr);
      exit(1);
    } else {
      basenames.push_back(string(argv[i]));
    }
  }

  if(basenames.size() < 1) {
    printUsage(argv[0], cerr);
    exit(1);
  }

  cerr << "Stopping threshold: -" << stoppingthreshold*100 << "th percentile of the distribution of correlations\n";


  // Find "unidentifiable" candidate transcripts
  // Later remove transcripts with zero hits across all samples
  // or which are identifiable in at least one sample
  vector<string> candidates;
  vector<string> zeros;
  vector<string> all_features;
  map<string, map<string, double> > zeros_feature2iact;
  map<string, map<string, double> > zeros_feature2sd;
  map<string, map<string, double> > zeros_feature2efflen;
  map<string, int> numbermappedreads;
  vector<string> toremove;
  map<string, bool> isIdent;
  cerr << basenames[0] << " ";
  get_unidentifiable_transcripts(basenames[0], candidates, zeros, toremove, zeros_feature2iact, zeros_feature2sd, zeros_feature2efflen, numbermappedreads, all_features, isIdent);
  vector< vector<string> > all_zeros;
  all_zeros.push_back(zeros);
  for(int i=1; i < basenames.size(); i++) {
    vector<string> cands;
    vector<string> zh;
    cerr << basenames[i] << " ";
    get_unidentifiable_transcripts(basenames[i], cands, zh, toremove, zeros_feature2iact, zeros_feature2sd, zeros_feature2efflen, numbermappedreads, all_features, isIdent);

    all_zeros.push_back(zh);
    vector<string> newcands(candidates.size() + cands.size());
    sort(candidates.begin(), candidates.end());
    sort(cands.begin(), cands.end());
    strit=set_union(candidates.begin(), candidates.end(), cands.begin(), cands.end(), newcands.begin());
    candidates=newcands;
    candidates.resize(strit-newcands.begin());

    vector<string> newzh(zeros.size() + zh.size());
    sort(zeros.begin(), zeros.end());
    sort(zh.begin(), zh.end());
    strit=set_intersection(zeros.begin(), zeros.end(), zh.begin(), zh.end(), newzh.begin());
    zeros=newzh;
    zeros.resize(strit-newzh.begin());
  }

  { // Remove transcripts absent from all samples
    vector<string> newcandidates(candidates.size());
    sort(candidates.begin(), candidates.end());
    sort(zeros.begin(), zeros.end());
    strit=set_difference(candidates.begin(), candidates.end(), zeros.begin(), zeros.end(), newcandidates.begin());
    candidates=newcandidates;
    candidates.resize(strit-newcandidates.begin());
  }

  { // Remove transcripts with at least one hit or MCSE < thres in at least one sample
    vector<string> newcandidates(candidates.size());
    sort(candidates.begin(), candidates.end());
    sort(toremove.begin(), toremove.end());
    toremove.erase(unique(toremove.begin(), toremove.end()), toremove.end());
    strit=set_difference(candidates.begin(), candidates.end(), toremove.begin(), toremove.end(), newcandidates.begin());
    candidates=newcandidates;
    candidates.resize(strit-newcandidates.begin());
  }

  cerr << candidates.size() << " transcripts or sets of identical transcripts are unidentifiable in all samples.\n";
  if(debug) {
    ofs.open("candidates");
    for(strit=candidates.begin(); strit != candidates.end(); strit++) {
      ofs << *strit << endl;
    }
    ofs.close();ofs.clear();
  }

  // create hash of candidates to indexes
  map<string, int> cand2ind;
  for(strit=candidates.begin(); strit != candidates.end(); strit++) {
    cand2ind[*strit]=int(strit-candidates.begin());
  }

  // create binary T x S matrix telling us whether transcript t was observed in sample s
  umat S(candidates.size(), basenames.size());
  S.fill(1);
  for(int s=0; s < basenames.size(); s++) {
    for(strit=all_zeros[s].begin(); strit != all_zeros[s].end(); strit++) {
      if(cand2ind.count(*strit) > 0) {
        S(cand2ind[*strit], s) = 0;
      }
    }
  }

  mat M(TRACELEN, candidates.size());
  cube R(candidates.size(), candidates.size(), basenames.size());
  get_corrs(R, M, basenames, cand2ind);

  cerr << "Calculating mean and s.d. of correlations...";
  mat V(candidates.size(), candidates.size()); // mean correlations
  mat W(candidates.size(), candidates.size()); // sds
  vector<int> ts;
  for(int i=0; i < candidates.size(); i++) ts.push_back(i);

  mean_corrs(R, S, V, W, ts, vector<int>(), SDPENALTY);
  cerr << "done.\n";
 
  cerr << "Threshold for mean anti-correlation: ";
  vector<double> maxcorrs, mincorrs;
  for(map<string,int>::iterator it=cand2ind.begin(); it != cand2ind.end(); it++) {
    double m = -1;
    double mm = 1;
    for(map<string,int>::iterator it2=cand2ind.begin(); it2 != cand2ind.end(); it2++) {
      if( it->second != it2->second && V(it->second,it2->second) > m) {
        m = V(it->second,it2->second);
      }
      if( it->second != it2->second && V(it->second,it2->second) < mm) {
        mm = V(it->second,it2->second);
      }
    }
    maxcorrs.push_back(m);
    mincorrs.push_back(mm);
  }
  sort(maxcorrs.begin(), maxcorrs.end());
  sort(mincorrs.begin(), mincorrs.end());
  if(debug) {
    ofs.open("maxcorrs");
    for(int i=0 ;i < maxcorrs.size(); i++) {
      ofs << maxcorrs[i] << " ";
    }
    ofs << endl;
    ofs.close(); ofs.clear();
    ofs.open("mincorrs");
    for(int i=0 ;i < mincorrs.size(); i++) {
      ofs << mincorrs[i] << " ";
    }
    ofs << endl;
    ofs.close(); ofs.clear();

  }


  double thr= -maxcorrs[floor(maxcorrs.size() * stoppingthreshold)];
  cerr << thr << endl;

  uword minrow, mincol;
  double mincor=V.min( minrow, mincol );
  vec onesb = ones<vec>(basenames.size());
  if(debug) {
    ofs.open("minavecortrace");
  }
  int nco=0;
  vector<int> blacklist;
  while(mincor < thr) {
    cerr << "Collapsing unidentifiable transcripts based on mean anti-correlations...(min mean cor="
         << left << setw(10) << mincor << ")\r";
    vector<int> inds;
    inds.push_back(minrow);
    inds.push_back(mincol);
    sort(inds.begin(), inds.end());
    if(debug) {
      ofs << candidates[inds[0]] << " " << candidates[inds[1]] << " " << mincor << endl;
    }
    blacklist.push_back(inds[1]);
    collapse(R, candidates, cand2ind, inds);

    mean_corrs(R, S, V, W, inds, vector<int> (&inds[1], &inds[1]), SDPENALTY);
    mincor=V.min( minrow, mincol );
    if(find(blacklist.begin(), blacklist.end(), minrow)!=blacklist.end()) {
      cerr << "WOOPS " << minrow <<  " " << mincol << " " << V(minrow,mincol) << endl;
    }
    if(find(blacklist.begin(), blacklist.end(), mincol)!=blacklist.end()) {
      cerr << "BOOPS " << minrow <<  " " << mincol << " " << V(minrow,mincol) << endl;
    }


    if(debug) {
      maxcorrs.clear();
      mincorrs.clear();
      for(map<string,int>::iterator it=cand2ind.begin(); it != cand2ind.end(); it++) {
        double m = -1;
        double mm = 1;
        for(map<string,int>::iterator it2=cand2ind.begin(); it2 != cand2ind.end(); it2++) {
          if( it->second != it2->second && V(it->second,it2->second) > m) {
            m = V(it->second,it2->second);
          }
          if( it->second != it2->second && V(it->second,it2->second) < mm) {
            mm = V(it->second,it2->second);
          }
        }
        maxcorrs.push_back(m);
        mincorrs.push_back(mm);
      }
      sort(maxcorrs.begin(), maxcorrs.end());
      sort(mincorrs.begin(), mincorrs.end());
      if(debug) {
        ofs2.open("maxcorrs", std::ofstream::app);
        for(int i=0 ;i < maxcorrs.size(); i++) {
          ofs2 << maxcorrs[i] << " ";
        }
        ofs2 << endl;
        ofs2.close(); ofs2.clear();
        ofs2.open("mincorrs", std::ofstream::app);
        for(int i=0 ;i < mincorrs.size(); i++) {
          ofs2 << mincorrs[i] << " ";
        }
        ofs2 << endl;
        ofs2.close(); ofs2.clear();
      }
    }



    nco++;
  }
  cerr << "Collapsing unidentifiable transcripts based on mean anti-correlations..." 
       << "done (" << nco << left << setw(20) << " iterations)." << endl; 
  if(debug) {
    ofs.close(); ofs.clear();
  }


  for(int s=0; s < basenames.size(); s++) {
    vector<string> ids;
    map<string, bool> gottrace;
    cerr << "Reading posterior traces from \"" << basenames[s] << ".trace_gibbs.gz\"...";
    ifs.open((basenames[s] + ".trace_gibbs.gz").c_str(), ios_base::in | ios_base::binary);
    gifs.push(ifs);
    if(!ifs.good()) {
      cerr << "Error: couldn't open " << basenames[s] << ".trace_gibbs.gz.\n";
      exit(1);
    }
    getline(gifs,str); // get IDs
    tokens.clear();
    tokenise(str,tokens," ");
    int isident=0;
    for(strit = tokens.begin(); strit != tokens.end(); strit++) {
      if(isIdent.count(*strit) < 1) {
        ids.push_back(*strit);
        gottrace[*strit]=true;
      }
    }

    M.set_size(TRACELEN, ids.size());

    for(int i=0; i < TRACELEN; i++) {
      int v=0;
      for(int t=0; t < tokens.size(); t++) {
        gifs >> str;
        if(isIdent.count(tokens[t]) < 1) {
          M(i,v)=atof(str.c_str());
          v++;
        }
      }
    }
    gifs.pop();
    ifs.close();ifs.clear();
    cerr << "done.\n";
    cerr << "Reading posterior traces from \"" << basenames[s] << ".identical.trace_gibbs.gz\"...";
    ifs.open((basenames[s] + ".identical.trace_gibbs.gz").c_str(), ios_base::in | ios_base::binary);
    gifs.push(ifs);
    if(!ifs.good()) {
      cerr << "Error: couldn't open " << basenames[s] << ".identical.trace_gibbs.gz.\n";
      exit(1);
    }
    getline(gifs,str); // get IDs
    tokens2.clear();
    tokenise(str,tokens2," ");
    int firstset=ids.size();
    for(strit = tokens2.begin(); strit != tokens2.end(); strit++) {
      ids.push_back(*strit);
      gottrace[*strit]=true;
    }

    M.resize(TRACELEN, M.n_cols + tokens2.size());

    for(int i=0; i < TRACELEN; i++) {
      for(int t=0; t < tokens2.size(); t++) {
        gifs >> str;
        M(i,firstset + t)=atof(str.c_str());
      }
    }
    gifs.pop();
    ifs.close();ifs.clear();
    cerr << "done.\n";

    cerr << "Simulating traces for transcripts/sets of identical transcripts with no hits...\n";
    int addcol=0;
    for(strit = all_features.begin(); strit != all_features.end(); strit++) {
      if(gottrace.count(*strit) < 1) {
        addcol++;
      }
    }
    M.resize(TRACELEN, M.n_cols + addcol);

    int t=0;
    for(strit = all_features.begin(); strit != all_features.end(); strit++) {
      if(gottrace.count(*strit) < 1) {
        ids.push_back(*strit);
        for(int i=0; i < TRACELEN; i++) {
          if(zeros_feature2efflen[basenames[s]].count(*strit) < 1) {
            cerr << "Error: couldn't get effective length for " << *strit << endl;
            exit(1);
          }
          M(i,firstset + tokens2.size() + t) = 
                gsl_ran_gamma(rg, alpha, 1.0/(beta + zeros_feature2efflen[basenames[s]][*strit] * (double)numbermappedreads[basenames[s]]/1000000000.0));
        }
        t++;
      }
    }

    cerr << "done.\n";

    cerr << "\tSaving new table \"" << basenames[s] << ".collapsed.mmseq\".\n";

    vector<string> forcollapsing;
    for(strit=candidates.begin(); strit != candidates.end(); strit++) {
      if(strit->find("*") != string::npos) {
        forcollapsing.push_back(*strit);
      }
    }

    map<int,int> shed = join_traces(M, ids, forcollapsing);

    M=log(M);
    double mumcse[M.n_cols];
    double iact[M.n_cols];
    double mu_trace[M.n_rows];
    double sd[M.n_rows];
    for(int t=0; t < M.n_cols; t++) {
      if(shed.count(t) >0) continue;
      memcpy(mu_trace, &conv_to<std::vector<double> >::from(M.col(t))[0], sizeof( double ) * TRACELEN );
      double var=0;
      double tau=0;
      int m=0;
      int len=TRACELEN;
      if(sokal(&len,&mu_trace[0],&var,&tau,&m)!=0) {
        mumcse[t] = TRACELEN;
        iact[t] = NAN;
      } else {
        mumcse[t] = sqrt(tau*var/TRACELEN);
        iact[t] = tau;
      }
      sd[t]=sqrt(var);
    }

    map<string, int> feature2ind;
    // Calculate unique hits //
    boolMat MM(1,1);
    ifs.open((basenames[s] + ".M").c_str());
    getline(ifs,str);
    if(!ifs.good() || str[0]!='#') {
      cerr << "Error reading " << basenames[s] << ".M file." << endl;
      exit(1);
    }
    tokens2.clear();
    tokenise(str, tokens2, "\t");
    // tokens2[0] should contain "#"
    for(int i=1; i < tokens2.size(); i++) {
      feature2ind[tokens2[i]]=i-1;
    }

    int i, j, maxi=0, maxj=0;
    while(true) {
      ifs >> i;
      ifs >> j;
      if(ifs.eof()) break;
      if(maxi < i) maxi=i;
      if(maxj < j) maxj=j;

      int factor1 = 1; int factor2 = 1;
      if(MM.size1() < i+1 ) factor1=2;
      if(MM.size2() < j+1 ) factor2=2;
      if(factor1==2 || factor2==2 ) {
        boolMat temp(MM);
        int s1=factor1==2?(i+1)*factor1: MM.size1();
        int s2=factor2==2?(j+1)*factor2: MM.size2();
        MM.resize(s1,s2,false);
        for(boolMatIt1 i1 = temp.begin1(); i1 != temp.end1(); ++i1) {
          for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
            MM(i2.index1(),i2.index2())= *i2;
          }
        }
      }
      MM(i,j)=1;
    }
    boolMat temp(MM);
    MM.resize(maxi+1,maxj+1,false);
    for(boolMatIt1 i1 = temp.begin1(); i1 != temp.end1(); ++i1) {
      for (boolMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2) {
        MM(i2.index1(),i2.index2())= *i2;
      }
    }
    ifs.close(); ifs.clear();

    vector<int> k;
    ifs.open((basenames[s] + ".k").c_str());
    ifs >> str;
    while(!ifs.eof()) {
      k.push_back(atoi(str.c_str()));
      ifs >> str;
    }
    ifs.close(); ifs.clear();

    boolMat T(feature2ind.size(), ids.size());
    for(int t=0; t < ids.size(); t++) {
      if(shed.count(t)>0) continue;
      tokens2.clear();
      tokenise(ids[t], tokens2, "+");
      for(strit2 = tokens2.begin(); strit2 != tokens2.end(); strit2++) {
        tokens3.clear();
        tokenise(*strit2, tokens3, "*");
        for(strit3 = tokens3.begin(); strit3 != tokens3.end(); strit3++) {
          if(feature2ind.count(*strit3)>0 ) {
            T(feature2ind[*strit3],t)=1;
          }
        }
      }
    }
    vector<int> unique_hits =uh(MM,T,k);
    // // // //

    ofs.open((basenames[s] + ".collapsed.mmseq").c_str());
    if(!ofs.good()) {
      cerr << "Error: cannot open " << basenames[s] << ".collapsed.mmseq for writing" << endl;
      exit(1);
    }
    ifs.open((basenames[s] + ".mmseq").c_str());
    if(!ifs.good()) {
      cerr << "Error: cannot open " << basenames[s] << ".mmseq" << endl;
      exit(1);
    }
    
    getline(ifs,str);
    // skip comments
    while(ifs.good() && str[0]=='#') {
      ofs << str << endl;
      getline(ifs,str);
    }
    ifs.close(); ifs.clear();

    /*

    map<string, bool> iszero;

    vector<string> all_features;
    for(int t = 0; t < M.n_cols; t++) {
      all_features.push_back(tokens[t]);
    }
    for(strit = all_zeros[s].begin(); strit != all_zeros[s].end(); strit++) {
      all_features.push_back(*strit);
      iszero[*strit]=true;
    }
    sort(all_features.begin(), all_features.end());

    cerr << "All features " << all_features.size() << endl;

    ofs << "feature_id\tlog_mu\tmcse\tweighted_length\tiact\tunique_hits\n";
    for(strit = all_features.begin(); strit != all_features.end(); strit++) {
      if(iszero.count(*strit)<1) {
        if(shed.count(feature2ind[*strit]) > 0) continue;
        double thismu=conv_to<double>::from(M.col(feature2ind[*strit]).t() * ones<colvec>(TRACELEN) / TRACELEN);
        if(isfinite(thismu)) {
          ofs << tokens[feature2ind[*strit]] << "\t" << thismu << "\t"
            << mumcse[feature2ind[*strit]] << "\t" << "NA" << "\t" << iact[feature2ind[*strit]] << "\t" << unique_hits[feature2ind[*strit]] << endl;
        } else {
          cerr << "shouldn't happen\n";
          exit(1);
        }
      } else {
          ofs << *strit << "\t" << "MU" << "\t"
              << zeros_feature2mcse[basenames[s]][*strit] << "\t" << "NA" 
              << zeros_feature2iact[basenames[s]][*strit] << "\t" << "0" << endl;
      }
    }
*/

    vector<pair<string,int> > strind;
    for(int i=0; i < ids.size(); i++) {
      pair<string, int> pa;
      pa.first=ids[i];
      pa.second=i;
      strind.push_back(pa);
    }
    sort(strind.begin(), strind.end(), strorder);

    ofs << "feature_id\tlog_mu\tsd\tmcse\teffective_length\tiact\tunique_hits\n";
    for(int v = 0; v < M.n_cols; v++) {
      int t=strind[v].second;
      if(shed.count(t) > 0) continue;
      double thismu=conv_to<double>::from(M.col(t).t() * ones<colvec>(TRACELEN) / TRACELEN);
      if(isfinite(thismu)) {
        ofs << ids[t] << "\t" << thismu << "\t"
          << sd[t] << "\t" << mumcse[t] << "\t" << "NA" << "\t" << iact[t] << "\t" << unique_hits[t] << endl;
      } else {
        cerr << "shouldn't happen\n";
        exit(1);
      }
    }

    ofs.close(); ofs.clear();
  }

  if(debug) {
    cerr << "Trace file: minavecortrace\n";
    cerr << "Candidates file: candidates\n";
  }

/*  for(int t=0; t < candidates.size(); t++) {
    cout << candidates[t] << " ";
  }
  cout << endl;
  for(int t=0; t < candidates.size(); t++) {
    for(int v=0; v < candidates.size(); v++) {
      cout << V(t,v) << " ";
    }
    cout << endl;
  }
*/
  return 0;
}



