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

#ifndef _BMS_H
#define _BMS_H

#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <math.h>
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
#include <gsl/gsl_sf_gamma.h>

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

using namespace std;
using namespace arma;


class BMS {
  private:
  vector<Col<double> > alpha;
  vector<Col<double> > alphaS;
  vector<Col<double> > alphaSS;
  vector<Col<int> > alphaN;
  vector<Col<double> > A;
  vector<Col<double> > Valpha;
  double v_alpha, v_alpharoot;
  vector<Mat<double> > D;
  vector<Mat<double> > S_inv;
  vector<Mat<double> > F;
  vector<Mat<double> > Veta;
  vector<Mat<double> > B;
  vector<Mat<double> > Vbeta;
  vector<Mat<double> > J;
  vector<Mat<double> > L;
  vector<Col<double> > Q;
  vector<Col<double> > R;
  double v_beta, v_betaroot;
  double v_gamma, v_gammaroot;
  double s;
  double d;
  vector<double> logitp;
  vector<double> LOsum;
  vector<bool> istuned;
  Col<int> gamma;
  vector< Mat<double> > beta;
  vector< Mat<double> > betaS;
  vector< Mat<double> > betaSS;
//  Cube<double> Vbeta;
//  Cube<double> VbetaCholeskyTrans;
  vector<Mat<double> > eta;
  vector<Mat<double> > etaS;
  vector<Mat<double> > etaSS;
  vector<Mat<int> > etaN;
  vector< Mat<double> > lambda;
  vector< Mat<double> > lambda_invS;
  vector< Mat<double> > lambda_loginvS;
  vector< Mat<double> > nu;
  vector< Mat<double> > sigmasq;
  vector< Mat<double> > sigmasq_invS;
  vector< Mat<double> > sigmasq_loginvS;
  vector< Col<double> > rho;
  vector< Col<double> > rhoS;
  vector< Col<double> > rho_logS;
  double k, q, r;
  int stilltotune_;
  vector<int> justtuned;
  vector<Col<int> > n; // observations in each condition

  vector<Col<double> > Mb, Mb1; // store M*beta, one for each thread

  vector<Col<double> >Peta; // store P*eta, one for each thread

  const Mat<double> *y;
  Mat<double> esq;

  Col<double> gammasum;
  int runlen;

  vector< Mat<double> > _sigacc; // accepted sigmas
  vector< Mat<double> > _sigprop; // proposed sigmas

  Mat<double> M;
  vector< Mat<double> > P;
  Mat<int> C;
  vector<int> nclasses; // number of classes for each model

  bool Mnil;
  vector<bool> Pnil;
  bool fixalpha;

  static const double g=2.0;
  gsl_rng ** rg;
  int max_threads;

  // ofstreams
  ofstream *ofs;
  int nstreams;
  vector<string> params;
  map<string, int> pind;
  string tracedir;
  const string burnin_suffix;

  Col<int> m2v(Mat<double> P);
  void tokenise(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ");
  double rnd_inv_gauss(gsl_rng * rg, double mu, double lambda);
  vector<double> meanvar(Col<double> v);
  double var_from_S_SS(double x, double xx, int n);
  Col<int> hadprod(Col<int> a, Col<int> b);

  ofstream *lpof;

  bool recording; // whether or not we are recording S and SS for burnin estimates

  public:
  BMS(const Mat<double> * _y, const Mat<double> & _e, Mat<double> _M, Mat<double> _P0, Mat<double> _P1, Mat<int> _C,
    double pdash, double d, double s, string _tracedir, gsl_rng ** _rg, const int _max_threads, bool _fixalpha);
   
  ~BMS();

  void initialise_streams(string suffix=string(""));

  void close_streams();

  int nfeatures();

  void print(bool also_gamma);
  void printtune(int batchlen);

  void set_pseudoprior_alpha();
  void set_pseudoprior_beta();
  void set_pseudoprior_eta();
  void set_pseudoprior_sigmasq();
  void set_pseudoprior_lambda();
  void set_pseudoprior_rho();

  void print_pseudo();

  double getp(int feature);
  void tunep(int feature, int k, int batchlen);
  bool tuned(int feature);
  int stilltotune();
  void updatetuned();

  void copy_pseudo();

  void update_alpha(int feature, int model, bool fit);

  void update_beta(int feature, int model, bool fit);

  void update_eta(int feature, int model, bool fit);

  void update_lambda(int feature, int model, bool fit);

  void update_sigmasq(int feature, int model, bool fit, bool rec_ar=true);
  void update_sigmasq_RW(int feature, int model, bool fit, bool rec_ar=true);

  void update_rho(int feature, int model, bool fit);

  void update_gamma(int feature);

  double log_target_posterior(int feature, int model);
  double log_target_pseudo(int feature, int model);

  double gammamean(int feature);
  double alphamean(int model, int feature);
  double etamean(int model, int cov, int feature);
  int nc(int model);
  int pc(int model);

  double sigar(int model, int feature, int condition);

  void reset();

  void tick();
  void record();
};

#endif // _BMS_H
