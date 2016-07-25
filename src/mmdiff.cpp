/*    Copyright 2012 Ernest Turro and Will J Astle
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
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string>
#include <string.h>
#include <fstream>
#include <assert.h>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <armadillo>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <unistd.h>
#include <sys/stat.h>

#include "bms.hpp"

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

#define OUTLEN 1024
#define MAXBATCHES 8192
#define LOTSOFFEATURES 500000

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

bool endsWith(const string& a, const string& b) {
  if (b.size() > a.size()) return false;
  return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}


void parse_mmseq(vector<string> filenames, vector<string> & features, Mat<double> & y, Mat<double> & e,
                 Mat<double> & uh, int range_start, int range_end, double trace_len, bool useprops) {
  ifstream ifs;
  string str;
  vector<string> tokens;
  int feature_ind=-1, lmg_ind=-1, sd_ind=-1, mcse_ind=-1, iact_ind=-1, mp_ind=-1, sp_ind=-1, uh_ind=-1;
  y.resize(LOTSOFFEATURES,filenames.size());
  e.resize(LOTSOFFEATURES, filenames.size());
  uh.resize(LOTSOFFEATURES, filenames.size());

  cerr << "Parsing ";
  for(int i=0; i < filenames.size(); i++) {
    cerr << filenames[i] << " ";
    ifs.open(filenames[i].c_str());
    if(!ifs.good()) {
      cerr << "Error: couldn't open " << filenames[i] << endl;
      exit(1);
    }
    getline(ifs,str);
    while(str[0]=='#') getline(ifs,str);
    tokens.clear();
    tokenise(str, tokens, "\t");
    for(int j=0; j < tokens.size(); j++) {
      if(tokens[j]=="feature_id") feature_ind=j;
      if(tokens[j]=="log_mu") lmg_ind=j;
      if(tokens[j]=="sd") sd_ind=j;
      if(tokens[j]=="mcse") mcse_ind=j;
      if(tokens[j]=="iact") iact_ind=j;
      if(tokens[j]=="mean_probit_proportion") mp_ind=j;
      if(tokens[j]=="sd_probit_proportion") sp_ind=j;
      if(tokens[j]=="unique_hits") uh_ind=j;
    }
    if(useprops) {
      if((feature_ind+1)*(mp_ind+1)*(sp_ind+1)*(uh_ind+1) == 0) {
        cerr << "Error: input tables must have feature_id, mean_probit_proportion, sd_probit_proportion and unique_hits columns.\n";
        exit(1);
      }
    } else {
      if((feature_ind+1)*(lmg_ind+1)*(sd_ind+1)*(uh_ind+1) == 0) {
        cerr << "Error: input tables must have feature_id, log_mu, sd and unique_hits columns.\n";
        exit(1);
      }
    }
    int k=0;
    while(true) {
      getline(ifs,str);
      if(ifs.eof()) break;
      tokens.clear();
      tokenise(str, tokens, "\t");
      if(k + 1 > features.size()) {
        features.push_back(tokens[feature_ind]);
      } else {
        if(features[k] != tokens[feature_ind]) {
          cerr << "Error: features across files do not match (" << k << ","
               << tokens[feature_ind] << "," << features[k] << ")\n";
          exit(1);
        }
      }
      if(k+1 > y.n_rows) {
        y.resize((k+1)*2,filenames.size());
        e.resize((k+1)*2,filenames.size());
        uh.resize((k+1)*2,filenames.size());
      }
      if(!useprops && tokens[lmg_ind]=="NA" || useprops && tokens[mp_ind]=="NA") {
        cerr << "Error: encountered NA" << endl;
        exit(1);
      } else {
        if(!useprops) {
          y(k,i) = atof(tokens[lmg_ind].c_str());
          e(k,i) = atof(tokens[sd_ind].c_str());
          uh(k,i) = atof(tokens[uh_ind].c_str());
        } else {
          y(k,i) = atof(tokens[mp_ind].c_str());
          e(k,i) = atof(tokens[sp_ind].c_str());
          uh(k,i) = atof(tokens[uh_ind].c_str());
        }
      }
      k++;
    }
    y.resize(k, filenames.size());
    e.resize(k, filenames.size());
    uh.resize(k, filenames.size());
    ifs.close(); ifs.clear();
  }
  if(range_start >=0 && range_end > range_start && range_end < y.n_rows) {
    y = y.submat(range_start, 0, range_end, y.n_cols-1);
    e = e.submat(range_start, 0, range_end, e.n_cols-1);
    uh = uh.submat(range_start, 0, range_end, uh.n_cols-1);
    vector<string> features2;
    for(int i=range_start; i <= range_end; i++) features2.push_back(features[i]);
    features=features2;
  }
  if(useprops) {
    Mat<double> newy=y;
    Mat<double> newe=e;
    Mat<double> newuh=uh;
    vector<string> newfeatures=features;

    for(int i=y.n_rows-1; i >= 0; i--) {
      if(isinf(y(i,0))) {
        newy.shed_row(i);
        newe.shed_row(i);
        newuh.shed_row(i);
        newfeatures.erase(newfeatures.begin()+i);
      } 
    }

    y = newy;
    e = newe;
    uh = newuh;
    features = newfeatures;
    cerr << endl << "Kept " << y.n_rows << " transcripts belonging to multi-isoform genes";
  } else {
    cerr << endl << "Analysing " << y.n_rows << " features";
  }
  cerr << endl;
}

// DESeq style
// y[,i] <- y[,i] - median(y[,i] - rowMeans(y))
// apply only using rows where fraction of samples having uh>=1 is greater than uhfrac
void apply_normalisation(vector<string> filenames, Mat<double> & y, Mat<double> & uh, double uhfrac) {
  vector<int> use; // use only these features
  for(int i=0; i < uh.n_rows; i++) {
    int sum=0;
    for(int j=0; j < uh.n_cols; j++) {
      if(uh(i,j) > 0) sum++;
    }
    if((double)sum/(double)uh.n_cols >= uhfrac) {
      use.push_back(i);
    }
  }

  if(use.size() < 100) {
    cerr << "Warning: fewer than 100 features found for normalisation. Skipping.\n";
    return;
  }

  cerr << "Using " << use.size() << "/" << y.n_rows << " features for normalisation.\n";

  vector<double> rowMeans;
  for(int i=0; i < use.size(); i++) {
    double sum=0;
    for(int j=0; j < y.n_cols; j++) {
      sum+= y(use[i], j);
    }
    rowMeans.push_back(sum/(double)y.n_cols);
  }

  cerr << "Log scale normalisation factors:\n";
  for(int sample=0; sample < y.n_cols; sample++) {
    vector<double> ydiff;
    for(int i=0; i < use.size(); i++) {
      ydiff.push_back(y(use[i], sample) - rowMeans[i]);
    }
    sort(ydiff.begin(), ydiff.end());
    cerr << "\t" << filenames[sample]  << "\t" << ydiff[ydiff.size()/2] << endl;
    for(int i=0; i < y.n_rows; i++) {
      y(i, sample) = y(i, sample) - ydiff[ydiff.size()/2];
    }
  }
}

void apply_permutation(Mat<double> & y, Mat<double> & e) {
  cerr << "Permuting input data...";
  Mat<double> newy = y;
  Mat<double> newe = e;
  vector<int> randinds;
  for(int j=0; j< y.n_cols; j++) randinds.push_back(j);
  for(int i=0; i < y.n_rows; i++) {
    random_shuffle(randinds.begin(), randinds.end());
    for(int j=0; j < y.n_cols; j++) {
      newy(i,j) = y(i, randinds[j]);
      newe(i,j) = e(i, randinds[j]);
    }
  }
  y=newy;
  e=newe;
  cerr << "done\n";
}

void parse_matrices(string file, Mat<double> & M, Mat<double> & P0, Mat<double> & P1, Mat<int> & C, int nrows) {
  ifstream ifs;
  ifs.open(file.c_str());
  if(!ifs.good()) {
    cerr << "Error: couldn't open " << file << endl;
    exit(1);
  }
  string str;
  vector<string> tokens;
  int m=-1,i=0;
  int P0reducedcols=0;
  int P1reducedcols=0;
  while(true) {
    bool spacer=false;
    do {
      if(!ifs.eof()) {
        getline(ifs,str);
      } else {
        str="";
      }
      str=str.substr(0, str.find_first_of("#"));
      tokens.clear();
      tokenise(str, tokens, " ");
      if(tokens.size() < 1) spacer=true;
    } while(!ifs.eof() && tokens.size() < 1);
    if(str.size()==0 && ifs.eof()) break;
    if(spacer) {
      m++;
      i=0;
    }
    if(m==-1) m=0;

    if(m==0) M.resize(nrows,tokens.size());
    if(m==1) C.resize(nrows,tokens.size());
    if(m==2) P0.resize(nrows,tokens.size());
    if(m==3) P1.resize(nrows,tokens.size());
    if(i >= nrows) {
      cerr << "Error: number of rows of matrices greater than number of samples.\n";
      exit(1);
    }
    for(int t=0; t < tokens.size(); t++) {
      if(m==0) M(i,t) = atof(tokens[t].c_str());
      if(m==1) C(i,t) = atoi(tokens[t].c_str());
      if(m==2) {
        for(int j=0; j < nrows; j++) {
          if(i > max(C.col(0))) { 
              cerr << "Error: more distinct rows of P than classes for model 0 (" << i << " > " << max(C.col(0)) << ").\n";
              exit(1);
          }
          if(C(j,0)==i) {
            P0(j,t)= atof(tokens[t].c_str());
          }
        }
        P0reducedcols = max(P0reducedcols, i);
      }
      if(m==3) {
        for(int j=0; j < nrows; j++) {
          if(i > max(C.col(1))) { 
              cerr << "Error: more distinct rows of P than classes for model 1 (" << i << " > " << max(C.col(1)) << ").\n";
              exit(1);
          }
          if(C(j,1)==i) {
            P1(j,t)= atof(tokens[t].c_str());
          }
        }
        P1reducedcols = max(P1reducedcols, i);
      }
    }
    i++;
  }
  if(C.min() != 0) {
    cerr << "Error: need at least one class in each model labelled 0." << endl;
    exit(1);
  }
  if(max(C.col(0)) != P0reducedcols) {
    cerr << "Error: number of classes does not correspond to number of disinct rows of P for model 0.\n";
    exit(1);
  }
  if(max(C.col(1)) != P1reducedcols) {
    cerr << "Error: number of classes does not correspond to number of disinct rows of P for  model 1.\n";
    exit(1);
  }

}

void printUsage(char *bin, ostream& out) {
  out << "Usage: mmdiff [OPTIONS...] [-de n1 n2 ... nC | -m matrices_file] mmseq_file1 mmseq_file2... > out.mmdiff" << endl
      << "       matrices_file contains M P0 P1 each separated by an empty line" << endl
      << endl;
  out << "Mandatory arguments:" << endl
      << "  ONE OF:" << endl
      << "  -de INT INT...    simple differential expression between several groups of samples, where" << endl
      << "                    each INT corresponds to a grouping of MMSEQ files into one condition" << endl
      << "  -m STRING         path to matrices file specifying the two models to compare" << endl;
  out << "Optional arguments:" << endl
      << "  -tracedir STRING  directory in which to save MCMC traces (default: "" (do not write traces))" << endl
      << "  -useprops         run on isoform/gene proportions instead of expression" << endl
      << "  -permute          run on permuted dataset; combine with non-permuted results to obtain q-values" << endl
      << "  -p FLOAT          prior probability of the second model (default: 0.1)" << endl
      << "  -d FLOAT          d hyperparameter (default: 1.4)" << endl 
      << "  -s FLOAT          s hyperparameter (default: 2.0)" << endl 
      << "  -l INT            length of MCMC trace used to produce the input estimates (default: 1024)" << endl
      << "  -fixalpha         fix alpha=0 and do not update (default: estimate alpha)" << endl
      << "  -nonorm           do not normalise the input data" << endl
      << "  -pdash FLOAT      initial value of pdash for improved mixing of gamma (stays fixed if -notune is set) (default: 0.5)" << endl
      << "  -notune           do not tune pdash to improve mixing for gamma" << endl
      << "  -uhfrac FLOAT     if normalising, use features which have at least one unique hit" << endl
      << "                    in at least uhfrac of the samples (default: max(0.2, (N - floor(N^2/160))/N))" << endl
      << "  -burnin INT       burnin iterations (default: 8192)" << endl
      << "  -iter INT         MCMC iterations (default: 16384)" << endl
      << "  -seed INT         seed for PRNG (default: 1234)" << endl
      << "  -range INT INT    select features indexed within range (default: all)" << endl;
}

int main(int argc, char** argv) {
  string matrices_file="";
  string tracedir="";
  double p=0.1;
  double d=1.4;
  double s=2.0;
  int l=1024;
  int burnin=8192;
  int mcmciters=16384;
  int seed=1234;
  int range_start=-1;
  int range_end=-1;
  bool useprops=false;
  bool fixalpha=false;
  bool normalise=true;
  bool customuhfrac=false;
  double uhfrac = -1;
  bool permute=false;
  double pdash=0.5;
  int batchlen=128;
  int tuneiters=batchlen;
  int numbatches=1;
  vector <int> simple_de;
  int ss=0;

  vector <string> arguments;

  for(int i=1; i < argc; i++) arguments.push_back(string(argv[i]));

  while(true) {
    if(arguments.size() >0 && arguments[0]=="-tracedir") {
      arguments.erase(arguments.begin());
      tracedir=arguments[0];
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-m") {
      arguments.erase(arguments.begin());
      for(int i=0; i < arguments.size(); i++) {
        if(arguments[i].find("-")==0) {
          cerr << "Error: optional arguments must be specified before -de or -m." << endl << endl;
          printUsage(argv[0], cerr);
          exit(1);
        }
      }
      matrices_file=arguments[0];
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-de") {
      arguments.erase(arguments.begin());
      for(int i=0; i < arguments.size(); i++) {
        if(arguments[i].find("-")==0) {
          cerr << "Error: optional arguments must be specified before -de or -m." << endl << endl;
          printUsage(argv[0], cerr);
          exit(1);
        }
      }
      cerr << "Number of samples in each group:";
      while(ss < arguments.size()) {
        simple_de.push_back(atoi(arguments[0].c_str()));
        cerr << " " << simple_de.back();
        arguments.erase(arguments.begin());
        if(simple_de.back() < 1) {
          cerr << endl << "Error: each grouping must contain at least one sample" << endl;
          exit(1);
        }

        ss+=simple_de.back();
      }
      cerr << endl;
      if(ss != arguments.size()) {
        cerr << "Error: total number of samples specified with -de must equal number of MMSEQ files" << endl;
        exit(1);
      }
    } else if(arguments.size() >0 && arguments[0]=="-useprops") {
      arguments.erase(arguments.begin());
      useprops=true;
    } else if(arguments.size() >0 && arguments[0]=="-p") {
      arguments.erase(arguments.begin());
      p=strtod(arguments[0].c_str(), NULL);
      if(p < 0 || p > 1) { 
        cerr << "Error: p must be between 0 and 1.\n";
        exit(1);
      }
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-s") {
      arguments.erase(arguments.begin());
      s=strtod(arguments[0].c_str(), NULL);
      if(s<=0) {
        cerr << "Error: s must be positive.\n";
        exit(1);
      }
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-d") {
      arguments.erase(arguments.begin());
      d=strtod(arguments[0].c_str(), NULL);
      if(d<=0) {
        cerr << "Error: d must be positive.\n";
        exit(1);
      }
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-pdash") {
      arguments.erase(arguments.begin());
      pdash=strtod(arguments[0].c_str(), NULL);
      if(pdash < 0 || pdash > 1) { 
        cerr << "Error: pdash must be between 0 and 1.\n";
        exit(1);
      }
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-fixalpha") {
      arguments.erase(arguments.begin());
      fixalpha=true;
    } else if(arguments.size() >0 && arguments[0]=="-l") {
      arguments.erase(arguments.begin());
      l=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-burnin") {
      arguments.erase(arguments.begin());
      burnin=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-iter") {
      arguments.erase(arguments.begin());
      mcmciters=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-seed") {
      arguments.erase(arguments.begin());
      seed=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-nonorm") {
      arguments.erase(arguments.begin());
      normalise=false;
    } else if(arguments.size() >0 && arguments[0]=="-notune") {
      arguments.erase(arguments.begin());
      tuneiters=0;
    } else if(arguments.size() >0 && arguments[0]=="-permute") {
      arguments.erase(arguments.begin());
      permute=true;
    } else if(arguments.size() >0 && arguments[0]=="-uhfrac") {
      arguments.erase(arguments.begin());
      customuhfrac=true;
      uhfrac=atof(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() >0 && arguments[0]=="-range") {
      arguments.erase(arguments.begin());
      range_start=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
      range_end=atoi(arguments[0].c_str());
      arguments.erase(arguments.begin());
    } else if(arguments.size() > 0 && ( arguments[0]=="-h" || arguments[0]=="--help" || arguments[0]=="-help")) {
      cerr << "Bayesian model selection for RNA-seq expression estimates.\n";
      printUsage(argv[0], cerr);
      exit(1);
    } else if(arguments.size() > 0 && ( arguments[0]=="-v" || arguments[0]=="--version" || arguments[0]=="-version")) {
      cerr << "mmdiff-" << QUOTE(VERSION) << endl;
      exit(1);
    } else {
      if(arguments.size()>0 && arguments[0][0]=='-') {
        cerr << "Error: unrecognised option " << arguments[0] << ".\n";
        printUsage(argv[0], cerr);
        exit(1);
      } else if(arguments.size() <= 2) {
        cerr << "Error: mandatory arguments missing.\n";
        printUsage(argv[0], cerr);
        exit(1);
      } else {
        break;
      }
    }
  }

  if(burnin <= 0 || mcmciters <= 0) {
    cerr << "Error: negative burnin and iter parameters.\n";
    printUsage(argv[0], cerr);
    exit(1);
  }
  if(burnin % OUTLEN != 0 || mcmciters % OUTLEN != 0) {
    cerr << "Error: burnin and iter parameters must be multiples of " << OUTLEN << endl;
    printUsage(argv[0], cerr);
    exit(1);
  }

  if(useprops) {
    cerr << "Using proportions, therefore disabling normalisation.\n";
    normalise=false;
  }

  if(matrices_file=="" && simple_de.size()==0) {
    cerr << "Error: either -de or -m must be specified" << endl;
    exit(1);
  }

  if(matrices_file=="" && simple_de.size()==1) {
    cerr << "Error: -de requires at least two groupings" << endl;
    exit(1);
  }

  vector <string> filenames;
  for(int i=0; i < arguments.size(); i++) {
    filenames.push_back(arguments[i]);
  }

  if(customuhfrac && (uhfrac > 1 || uhfrac < 1.0/(double)filenames.size())) {
    cerr << "Error: uhfrac must be <= 1 and >= 1/N." << endl;
    exit(1);
  }

  if(!customuhfrac) {
    uhfrac = max(0.2, (double)(filenames.size() - filenames.size()*filenames.size()/160)/(double)filenames.size());
  }
  if(normalise) {
    cerr << "Min unique hits fraction for normalisation: " << uhfrac << endl;
  }

  if(tracedir.size() > 0) {
    mkdir(tracedir.c_str(), 0755);
  }

  if(tracedir.size() > 0) {
    if(::access(tracedir.c_str(), F_OK | R_OK | W_OK | X_OK) == -1) {
      cerr << "Error: can't write to trace directory " << tracedir << "." << endl;
      exit(1);
    }
  }


  int max_threads=OMP_GET_MAX_THREADS;

  cerr << "No. threads: " << max_threads << endl;

  gsl_rng * rg[max_threads];
  for(int i=0; i < max_threads; i++) {
    rg[i] =  gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set(rg[i], seed+i);
  }

  vector<string> features;
  Mat<double> y;
  Mat<double> e;
  Mat<double> uh;

  Mat<double> M;
  Mat<double> P0;
  Mat<double> P1;
  Mat<int> C; // classes

  parse_mmseq(filenames, features, y, e, uh, range_start, range_end, l, useprops);
  if(normalise) apply_normalisation(filenames, y, uh, uhfrac);
  if(permute) apply_permutation(y, e);
  if(matrices_file=="") {
    M.zeros(ss, 1);
    C.zeros(ss, 2);
    P0.zeros(ss,1);
    if(simple_de.size()>2) {
      P1.zeros(ss,simple_de.size());
    } else {
      P1.zeros(ss,1);
    }
    int k=0;
    int l=0;
    for(int i=0; i < simple_de.size(); i++) {
      for(int j=0; j < simple_de[i]; j++) {
        C(k,1)=l;
        P0(k,0)=1.0;
        if(simple_de.size()>2) {
          P1(k,i)=1.0;
        } else {
          if(i==0) {
            P1(k,0)=.5;
          } else {
            P1(k,0)=-.5;
          }
        }
        k++;
      }
      l++;
    }
  } else {
    parse_matrices(matrices_file, M, P0, P1, C, y.n_cols);
  }

  if(!(filenames.size() == M.n_rows && M.n_rows == P0.n_rows && P0.n_rows == P1.n_rows)) {
    cerr << "Error: Rows in M, P0, P1 and number of mmseq files must match." << endl;
    exit(1);
  }

  BMS mcmc = BMS(&y, e, M, P0, P1, C, pdash, d, s, tracedir, rg, max_threads, fixalpha);

  cerr << "Design matrix for model 0 ([";
  if(!fixalpha) cerr << "1";
  if(!fixalpha && (!mcmc.Misnil() || !mcmc.Pisnil(0))) cerr << "|";
  if(!mcmc.Misnil()) cerr << "M";
  if(!mcmc.Pisnil(0) && (!fixalpha || !mcmc.Misnil())) cerr << "|";
  if(!mcmc.Pisnil(0)) cerr << "P0";
  cerr << "]):\n";
  Mat<double> Z= ones(P0.n_rows, fixalpha ? 0 : 1);
  if(!mcmc.Misnil()) Z.insert_cols(Z.n_cols,M);

  if(det(trans(Z)*Z)==0) {
    cerr << "Error: collinearity in combined matrix of intercept and covariates for model 0" << endl;
    exit(1);
  }

  if(!mcmc.Pisnil(0)) {
    Z.insert_cols(Z.n_cols,P0);
    if(det(trans(P0)*P0)==0) {
      cerr << "Error: collinearity in matrix P0" << endl;
      exit(1);
    }
  }

  cerr << Z;

  if(det(trans(Z)*Z)==0) {
    cerr << "Warning: collinearity in full matrix of predictors for model 0" << endl;
  }


  cerr << "Design matrix for model 1 ([";
  if(!fixalpha) cerr << "1";
  if(!fixalpha && (!mcmc.Misnil() || !mcmc.Pisnil(1))) cerr << "|";
  if(!mcmc.Misnil()) cerr << "M";
  if(!mcmc.Pisnil(1) && (!fixalpha || !mcmc.Misnil())) cerr << "|";
  if(!mcmc.Pisnil(1)) cerr << "P0";
  cerr << "]):\n";
  Z= ones(P1.n_rows, fixalpha ? 0 : 1);
  if(!mcmc.Misnil()) Z.insert_cols(Z.n_cols,M);

  if(det(trans(Z)*Z)==0) {
    cerr << "Error: collinearity in combined matrix of intercept and covariates for model 1" << endl;
    exit(1);
  }

  if(!mcmc.Pisnil(1)) {
    Z.insert_cols(Z.n_cols,P1);
    if(det(trans(P1)*P1)==0) {
      cerr << "Error: collinearity in matrix P1" << endl;
      exit(1);
    }
  }
  cerr << Z;

  if(det(trans(Z)*Z)==0) {
    cerr << "Warning: collinearity in full matrix of predictors for model 1" << endl;
  }

  if(mcmc.Misnil()) 
    cerr << "Note: no betas\n";
  if(mcmc.Pisnil(0)) 
    cerr << "Note: no etas in model 0\n";
  if(mcmc.Pisnil(1)) 
    cerr << "Note: no etas in model 1\n";

  int subsample_burnin=burnin / OUTLEN;
  int subsample_iter=mcmciters / OUTLEN;

  bool inburnin=true;

  inburnin=true;

#ifdef DEBUG
  int start,end, talpha=0, teta=0, tlambda=0, tsigmasq=0, trho=0, tgamma=0;
#endif

  for(int iter=0; (inburnin && iter < burnin) || (tuneiters > 0 && mcmc.stilltotune() > 0 && numbatches < MAXBATCHES) || iter < (tuneiters + mcmciters); iter++) {
    if(!inburnin && iter < tuneiters && iter > 0 && iter % batchlen == 0) mcmc.printtune(batchlen);
    if(!inburnin && iter+1==tuneiters && mcmc.stilltotune() > 0 && numbatches != MAXBATCHES) { tuneiters += batchlen; numbatches++; }

    if(inburnin && iter % subsample_burnin ==0) 
      cerr << "BURNIN ITER " << iter << "\r";
    if(!inburnin && iter % subsample_iter == 0) {
      if(iter < tuneiters) {
        cerr << "TUNING BATCH " << numbatches << " (" << mcmc.stilltotune() << " left)";
      } else {
        cerr << "TRACE ITER " << iter - tuneiters << " (sampling after " << numbatches << " tuning batches)";
      }
      cerr << "\r";
    }

    if(inburnin && iter == OUTLEN/10 || !inburnin && iter==tuneiters) {
      mcmc.record();
    }

#pragma omp parallel for schedule(static) shared(mcmc)
    for(int feature=0; feature < features.size(); feature++) {
      if(!inburnin && iter < tuneiters && iter > 0) {
        if(mcmc.tuned(feature)) continue;
        else if(iter % batchlen == 0) {
          mcmc.tunep(feature, iter, batchlen);
        }
      }

#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) start=clock();
#endif
      mcmc.update_alpha(feature,0, inburnin);
      mcmc.update_alpha(feature,1, inburnin);
#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) talpha += clock()-start;
#endif

      mcmc.update_beta(feature,0, inburnin);
      mcmc.update_beta(feature,1, inburnin);

#ifdef DEBUG
        if(OMP_GET_THREAD_NUM==0) start=clock();
#endif
        mcmc.update_eta(feature,0, inburnin);
        mcmc.update_eta(feature,1, inburnin);
#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) teta += clock()-start;
#endif

#ifdef DEBUG
        if(OMP_GET_THREAD_NUM==0) start=clock();
#endif
        mcmc.update_lambda(feature,0, inburnin);
#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) tlambda += clock()-start;
#endif
        mcmc.update_lambda(feature,1, inburnin);

#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) start=clock();
#endif
      mcmc.update_sigmasq_RW(feature,0, inburnin, !inburnin);
      mcmc.update_sigmasq_RW(feature,1, inburnin, !inburnin);
#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) tsigmasq += clock()-start;
#endif

#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) start=clock();
#endif
      mcmc.update_rho(feature,0, inburnin);
      mcmc.update_rho(feature,1, inburnin);
#ifdef DEBUG
      if(OMP_GET_THREAD_NUM==0) trho += clock()-start;
#endif

      if(!inburnin) {
#ifdef DEBUG
        if(OMP_GET_THREAD_NUM==0) start=clock();
#endif
        mcmc.update_gamma(feature);
#ifdef DEBUG
        if(OMP_GET_THREAD_NUM==0) tgamma += clock()-start;
#endif
      }
    }

    if(mcmc.stilltotune() > 0) {
      mcmc.updatetuned();
    }

    mcmc.tick();

    if(inburnin && iter % subsample_burnin ==0 || !inburnin && iter >= tuneiters && (iter - tuneiters) % subsample_iter == 0) {
      mcmc.print(!inburnin);
    }

    if(inburnin && iter==burnin-1) {
      mcmc.close_streams();

      cerr << "\nSetting pseudopriors...";
//      mcmc.set_pseudoprior_alpha(OUTLEN, OUTLEN/10);
      mcmc.set_pseudoprior_alpha();
      mcmc.set_pseudoprior_beta();
//      mcmc.set_pseudoprior_eta(OUTLEN, OUTLEN/10);
      mcmc.set_pseudoprior_eta();
//      mcmc.set_pseudoprior_lambda(OUTLEN, OUTLEN/10);
      mcmc.set_pseudoprior_lambda();
//      mcmc.set_pseudoprior_sigmasq(OUTLEN, OUTLEN/10);
      mcmc.set_pseudoprior_sigmasq();
//      mcmc.set_pseudoprior_rho(OUTLEN, OUTLEN/10);
      mcmc.set_pseudoprior_rho();
      
      mcmc.initialise_streams();
//      mcmc.copy_pseudo();
      inburnin=false;
      iter=-1;
      mcmc.reset();
      mcmc.print_pseudo();
      cerr << "done.\n";
    }
  }
  cerr << "\nDONE MCMC\n";

#ifdef DEBUG
  cerr << "talpha "  << talpha/1000000.0 << endl;
  cerr << "teta "  << teta/1000000.0 << endl;
  cerr << "tlambda "  << teta/1000000.0 << endl;
  cerr << "tsigmasq "  << tsigmasq/1000000.0 << endl;
  cerr << "trho "  << trho/1000000.0 << endl;
  cerr << "tgamma "  << tgamma/1000000.0 << endl;
#endif

  double postlogodds, BF;
  cout << "#prior_probability=" << p << endl;
  cout << "feature_id\tbayes_factor\tposterior_probability\t";

  for(int model=0; model < 2; model++) {
    if(!fixalpha) {
      cout << "alpha" << model << "\t";
    }
    if(!mcmc.Misnil()) {
      for(int l=0; l < M.n_cols; l++) {
        cout << "beta" << model << "_" << l << "\t";
      }
    }
    if(!mcmc.Pisnil(model)) {
      for(int l=0; l < mcmc.pc(model); l++) {
        cout << "eta" << model << "_" << l << "\t";
      }
    }
  }
  vector<string> samplenames = filenames;

  for(int f=0; f < samplenames.size(); f++) {
    if(endsWith(filenames[f],".mmseq")) {
      int found = filenames[f].find_last_of(".");
      int found2 = filenames[f].find_last_of("/");
      found2 = found2==string::npos ? -1 : found2;
      samplenames[f]=filenames[f].substr(found2+1,found - found2 -1);
    }
    cout << "mu_" << samplenames[f] << "\t";
  }
  for(int f=0; f < samplenames.size(); f++) {
    cout << "sd_" << samplenames[f];
    if(f < samplenames.size()-1) cout << "\t";
    else cout << "\n";
  }


  for(int feature=0; feature < features.size(); feature++) {
    BF = mcmc.gammamean(feature, true)/(1.0-mcmc.gammamean(feature, true)) * (1.0-mcmc.getp(feature))/mcmc.getp(feature);
    postlogodds = log(BF) + log(p) -log1p(-p);
    double pp = 1.0/(1.0+exp(-postlogodds));
    if(BF >= DBL_MAX) pp=1.0;
    cout << features[feature] << "\t" << BF << "\t" << pp << "\t";
    for(int model=0; model < 2; model++) {
      if(!fixalpha) {
        cout << mcmc.alphamean(model, feature) << "\t";
      }
      if(!mcmc.Misnil()) {
        for(int l=0; l < M.n_cols; l++) {
          cout << mcmc.betamean(model, l, feature) << "\t";
        }
      }
      if(!mcmc.Pisnil(model)) {
        for(int l=0; l < mcmc.pc(model); l++) {
          cout << mcmc.etamean(model, l, feature) << "\t";
        }
      }
    }
    for(int f=0; f < samplenames.size(); f++) {
      cout << y(feature,f) << "\t";
    }
    for(int f=0; f < samplenames.size(); f++) {
      cout << e(feature,f);
      if(f < samplenames.size()-1) cout << "\t";
      else cout << "\n";
    }
  }

  if(tracedir.size() > 0) {
    for(int model=0; model < 2; model++) {
      stringstream ss;
      ss << "tracedir" << "/sigar" << model << ".txt";
      ofstream sar((ss.str()).c_str());
      for(int feature=0; feature < features.size(); feature++) {
        for(int c =0; c < mcmc.nc(model); c++) {
          sar << mcmc.sigar(model, feature,c) << " ";
        }
        sar << endl;
      }
      sar.close();
    }
  }

  return 0;
}


