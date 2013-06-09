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

// Bayesian Model Selection for RNA-seq


#include "bms.hpp"

#define SQRT2 1.4142135623730950488016887242096980785696718753769480731
#define LOGIT07 0.8472979

double log_inv_gamma_pdf( double x, double shape, double scale) {
  return(-gsl_sf_lngamma(shape) + shape*log(scale) + (-shape-1.0) * log(x) -scale/x);
}


double inv_gamma_pdf( double x, double shape, double scale) {
  return(1.0/gsl_sf_gamma(shape)*pow(scale, shape)*pow(x,-shape-1.0) *exp(-scale/x));
}

Col<int> BMS::m2v(Mat<double> P) {
  Col<int> res(P.n_rows);
  for(int i=0; i < P.n_cols; i++) {
    for(int j=0; j < P.n_rows; j++) {
      if(P(j,i)>0.0000001) res(j) = i;
    }
  }
  return(res);
}

BMS::~BMS() {
  if(tracedir.size() > 0) {
    delete [] ofs;
  }
}

void BMS::close_streams() {
  if(tracedir.size() >0 ) {
    for(int i=0; i < nstreams; i++) {
      ofs[i].close();
      ofs[i].clear();
    }
  }
}

double BMS::rnd_inv_gauss(gsl_rng * rg, double mu, double lambda)
{
  // using algorithm on Wikipedia due to Michael, Schucany and Haas
  double v=gsl_ran_ugaussian(rg);
  double y = v*v;
  double x = mu + 0.5*((mu*mu*y)/lambda - (mu/lambda)* sqrt(4.0*mu*lambda*y + mu*mu*y*y));
 
  double test = gsl_rng_uniform(rg);  // sample from a uniform distribution between 0 and 1
  if (test > (mu)/(mu + x))
    return (mu*mu)/x;
  else
    return x;
}

vector<double> BMS::meanvar(Col<double> v) {
  double sum=0;
  double sum2=0;
  for(int i=0;i < v.n_elem; i++) {
    sum += v(i);
    sum2 += v(i)*v(i);
  }
  vector<double> res;
  res.push_back(sum/v.n_elem);
  res.push_back((sum2-sum*sum/(double)v.n_elem)/(v.n_elem-1.0));
  return(res);
}

double BMS::var_from_S_SS(double sum, double sum2, int n) {
  return((sum2-sum*sum/(double)n)/(n-1.0));
}


void BMS::print(bool also_gamma) {
  if(tracedir.size()==0) return;
  for(int feature=0; feature < y->n_rows; feature++) {
    for(int model=0; model < 2; model++) {
      if(model==0) ofs[pind["alpha0"]] << alpha[model](feature) << " ";
      else ofs[pind["alpha1"]] << alpha[model](feature) << " ";
      for(int i=0; i < M.n_cols; i++) {
        stringstream ss;
        ss << "beta" << model << "_" << i;
        ofs[pind[ss.str()]] << beta[model](i, feature) << " ";
      }
      for(int l=0; l < P[model].n_cols; l++) {
        stringstream ss;
        ss << "eta" << model << "_" << l;
        ofs[pind[ss.str()]] << eta[model](l, feature) << " ";
        ss.str("");
        ss << "lambda" << model << "_" << l;
        ofs[pind[ss.str()]] << lambda[model](feature, l) << " ";
        ss.str("");
      }
      for(int c=0; c < nclasses[model]; c++) {
        stringstream ss;
        ss << "sigmasq" << model << "_" << c;
        ofs[pind[ss.str()]] << sigmasq[model](c, feature) << " ";
      }
      stringstream ss;
      ss << "rho" << model;
      ofs[pind[ss.str()]] << rho[model](feature) << " ";
    }
    if(also_gamma) {
      ofs[pind["gamma"]] << gamma(feature) << " ";
    }
  }
  for(int i=0; i < nstreams; i++) {
    if( ( params[i]=="gamma" && !also_gamma ) || params[i]=="logitp") continue;
    ofs[i] << endl;
  }
}

void BMS::printtune(int batchlen) {
  if(tracedir.size()==0) return;
  for(int feature=0; feature < y->n_rows; feature++) {
    ofs[pind["meanLO"]] << LOsum[feature]/(double)batchlen << " ";
    ofs[pind["logitp"]] << logitp[feature] << " ";
  }
  ofs[pind["logitp"]] << endl;
  ofs[pind["meanLO"]] << endl;
}

int BMS::stilltotune() {
  return(stilltotune_);
}

/*
void BMS::set_pseudoprior_alpha(int len, int discard) {
  ifstream *ifs = new ifstream[nstreams];
  ifs[pind["alpha0"]].open((tracedir + "/alpha0" + burnin_suffix).c_str());
  ifs[pind["alpha1"]].open((tracedir + "/alpha1" + burnin_suffix).c_str());
  Cube<double> temp(len-discard, y->n_rows, 2); // read in for each of the two models
  double foo;
  for(int i=0; i < len; i++) {
    for(int j=0; j < y->n_rows; j++) {
      if(i < discard) {
        ifs[pind["alpha0"]] >> foo;
        ifs[pind["alpha1"]] >> foo;
      } else {
        ifs[pind["alpha0"]] >> temp(i-discard,j, 0);
        ifs[pind["alpha1"]] >> temp(i-discard,j, 1);
      }
    }
  }
  ifs[pind["alpha0"]].close();
  ifs[pind["alpha0"]].clear();
  ifs[pind["alpha1"]].close();
  ifs[pind["alpha1"]].clear();
  delete [] ifs;
  for(int model=0; model < 2; model++) {
    for(int j=0; j < y->n_rows; j++) {
      vector<double> res=meanvar(temp.slice(model).col(j));
      Valpha[model](j) = res[1];
      A[model](j) = res[0];
      if(model==1 && j==1) {
        cerr << res[0] << " " << res[1] << endl;
      }
//      cerr << " alpha pseudoprior[" << model << "](" <<  j << ") ~ N(A=" << res[0] << ", Valpha= " << res[1] << ")" << endl;
    }
  }
}

*/

void BMS::set_pseudoprior_alpha() {
  if(fixalpha) return;
  for(int model=0; model < 2; model++) {
    for(int j=0; j < y->n_rows; j++) {
      A[model](j) = alphaS[model](j)/(double)runlen;
      Valpha[model](j) = var_from_S_SS(alphaS[model](j), alphaSS[model](j), runlen); 
    }
  }
}

/*
void BMS::set_pseudoprior_beta(int len, int discard) {
  stringstream ss[2 * M.n_cols];
  for(int b=0; b < M.n_cols; b++) {
    ss[b] << "beta0_" << b;
  }
  for(int b=0; b < M.n_cols; b++) {
    ss[M.n_cols + b] << "beta1_" << b;
  }

  ifstream *ifs = new ifstream[nstreams];
  for(int b=0; b < M.n_cols; b++) {
    ifs[pind[ss[b].str()]].open((tracedir + "/" + ss[b].str() + burnin_suffix).c_str());
  }
  for(int b=0; b < M.n_cols; b++) {
    ifs[pind[ss[M.n_cols+b].str()]].open((tracedir + "/" + ss[M.n_cols + b].str() + burnin_suffix).c_str());
  }
  Cube<double> temp(len-discard, y->n_rows, 2*M.n_cols); // read in for each of the two models
  double foo;
  for(int i=0; i < len; i++) {
    for(int j=0; j < y->n_rows; j++) {
      for(int b=0; b < M.n_cols; b++) {
        if(i < discard) {
          ifs[pind[ss[b].str()]] >> foo;
        } else {
          ifs[pind[ss[b].str()]] >> temp(i-discard,j, b);
        }
      }
      for(int b=0; b < M.n_cols; b++) {
        if(i < discard) {
          ifs[pind[ss[P[0].n_cols+b].str()]] >> foo;
        } else {
          ifs[pind[ss[P[0].n_cols+b].str()]] >> temp(i-discard,j, M.n_cols + b);
        }
      }
    }
  }
  for(int b=0; b < M.n_cols; b++) {
    ifs[pind[ss[b].str()]].close();
    ifs[pind[ss[b].str()]].clear();
  }
  for(int b=0; b < M.n_cols; b++) {
    ifs[pind[ss[M.n_cols+b].str()]].close();
    ifs[pind[ss[M.n_cols+b].str()]].clear();
  }
  delete [] ifs;

  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) {
      for(int b=0; b < M.n_cols; b++) {
        vector<double> res=meanvar(temp.slice(P[0].n_cols*model + b).col(j));
        B[model](j, b) = res[0];
        Vbeta[model](j, b) = res[1];
      }
    }
  }
}
*/

void BMS::set_pseudoprior_beta() {
  if(Mnil) return;
  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) {
      for(int b=0; b < M.n_cols; b++) {
        B[model](j, b) = betaS[model](b, j)/(double)runlen;
        Vbeta[model](j, b) = var_from_S_SS(betaS[model](b, j), betaSS[model](b, j), runlen);
      }
    }
  }
}

/*
void BMS::set_pseudoprior_eta(int len, int discard) {
  stringstream ss[P[0].n_cols + P[1].n_cols];
  for(int c=0; c < P[0].n_cols; c++) {
    ss[c] << "eta0_" << c;
  }
  for(int c=0; c < P[1].n_cols; c++) {
    ss[P[0].n_cols + c] << "eta1_" << c;
  }

  ifstream *ifs = new ifstream[nstreams];
  for(int c=0; c < P[0].n_cols; c++) {
    ifs[pind[ss[c].str()]].open((tracedir + "/" + ss[c].str() + burnin_suffix).c_str());
  }
  for(int c=0; c < P[1].n_cols; c++) {
    ifs[pind[ss[P[0].n_cols+c].str()]].open((tracedir + "/" + ss[P[0].n_cols + c].str() + burnin_suffix).c_str());
  }
  Cube<double> temp(len-discard, y->n_rows, P[0].n_cols + P[1].n_cols); // read in for each of the two models
  double foo;
  for(int i=0; i < len; i++) {
    for(int j=0; j < y->n_rows; j++) {
      for(int c=0; c < P[0].n_cols; c++) {
        if(i < discard) {
          ifs[pind[ss[c].str()]] >> foo;
        } else {
          ifs[pind[ss[c].str()]] >> temp(i-discard,j, c);
        }
      }
      for(int c=0; c < P[1].n_cols; c++) {
        if(i < discard) {
          ifs[pind[ss[P[0].n_cols+c].str()]] >> foo;
        } else {
          ifs[pind[ss[P[0].n_cols+c].str()]] >> temp(i-discard,j, P[0].n_cols + c);
        }
      }
    }
  }
  for(int c=0; c < P[0].n_cols; c++) {
    ifs[pind[ss[c].str()]].close();
    ifs[pind[ss[c].str()]].clear();
  }
  for(int c=0; c < P[1].n_cols; c++) {
    ifs[pind[ss[P[0].n_cols+c].str()]].close();
    ifs[pind[ss[P[0].n_cols+c].str()]].clear();
  }
  delete [] ifs;

  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) {
      for(int c=0; c < P[model].n_cols; c++) {
        if(Pnil[model]) continue;
        vector<double> res=meanvar(temp.slice(P[0].n_cols*model + c).col(j));
        if(res[1] < 0.00001) {
          cerr << "Warning: no variance in burnin for eta "
          << model << " " << c << " " << j << ".\n";
          cerr << temp.slice(P[0].n_cols*model + c)(0,0) << endl;
          Veta[model](j, c) = 0.00001;
        } else {
          Veta[model](j, c) = res[1];
        }
        F[model](j, c) = res[0];
//        cerr << " eta pseudoprior[" <<  model << "](" << j << ", " << c << ") ~ N(F=" << res[0] << ", Veta=" << res[1] << ")" << endl;
      }
    }
  }
}
*/

void BMS::set_pseudoprior_eta() {
  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) {
      for(int c=0; c < P[model].n_cols; c++) {
        if(Pnil[model]) continue;
        F[model](j, c) = etaS[model](c, j)/(double)runlen;
        Veta[model](j, c) = var_from_S_SS(etaS[model](c,j), etaSS[model](c,j), runlen);
      }
    }
  }
}

/*
void BMS::set_pseudoprior_sigmasq(int len, int discard) {
  stringstream ss[nclasses[0] + nclasses[1]];
  for(int c=0; c < nclasses[0]; c++) {
    ss[c] << "sigmasq0_" << c;
  }
  for(int c=0; c < nclasses[1]; c++) {
    ss[nclasses[0] + c] << "sigmasq1_" << c;
  }

  ifstream *ifs = new ifstream[nstreams];
  for(int c=0; c < nclasses[0]; c++) {
    ifs[pind[ss[c].str()]].open((tracedir + "/" + ss[c].str() + burnin_suffix).c_str());
  }
  for(int c=0; c < nclasses[1]; c++) {
    ifs[pind[ss[nclasses[0]+c].str()]].open((tracedir + "/" + ss[nclasses[0] + c].str() + burnin_suffix).c_str());
  }
  Cube<double> temp(len-discard, y->n_rows, nclasses[0] + nclasses[1]); // read in for each of the two models
  Cube<double> temp2(len-discard, y->n_rows, nclasses[0] + nclasses[1]); // read in for each of the two models
  double foo;
  for(int i=0; i < len; i++) {
    for(int j=0; j < y->n_rows; j++) {
      for(int c=0; c < nclasses[0]; c++) {
        if(i < discard) {
          ifs[pind[ss[c].str()]] >> foo;
        } else {
          ifs[pind[ss[c].str()]] >> temp(i-discard,j, c);
          temp(i-discard,j,c) = 1.0/temp(i-discard,j,c) ; // Gamma distro scale
          temp2(i-discard,j,c) = log(temp(i-discard,j, c));
        }
      }
      for(int c=0; c < nclasses[1]; c++) {
        if(i < discard) {
          ifs[pind[ss[nclasses[0]+c].str()]] >> foo;
        } else {
          ifs[pind[ss[nclasses[0]+c].str()]] >> temp(i-discard,j, nclasses[0] + c);
          temp(i-discard,j, nclasses[0] + c) = 1.0/temp(i-discard,j, nclasses[0] + c) ;
          temp2(i-discard,j, nclasses[0] + c) = log(temp(i-discard,j, nclasses[0] + c));
        }
      }
    }
  }
  for(int c=0; c < nclasses[0]; c++) {
    ifs[pind[ss[c].str()]].close();
    ifs[pind[ss[c].str()]].clear();
  }
  for(int c=0; c < nclasses[1]; c++) {
    ifs[pind[ss[nclasses[0]+c].str()]].close();
    ifs[pind[ss[nclasses[0]+c].str()]].clear();
  }
  delete [] ifs;

  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) { 
      for(int c=0; c < nclasses[model]; c++) {
        vector<double> res=meanvar(temp.slice(nclasses[0]*model + c).col(j));
        if(res[1] < 0.0000001) {
          cerr << "Error: no variance in burnin for sigmasq (" << model << ", " << c << ", " << j << ") " << res[0] << ", " << res[1] << ".\n";
          exit(1);
        }
        vector<double> res2=meanvar(temp2.slice(nclasses[0]*model + c).col(j));
        double _s  = log(res[0]) - res2[0];
        J[model](j,c) = (3.0 - _s + sqrt(pow(_s-3.0, 2.0)+24.0*_s))/(12.0*_s);
        L[model](j,c) = J[model](j,c)/res[0];
//        cerr << " sigmasq pseudoprior[" <<  model << "](" << j << ", " << c << ") ~  1/Gamma(J=" << J[model](j,c)  << ", 1/L=" << 1.0/L[model](j,c) << "(scale))" << endl;
      }
    }
  }
}
*/

void BMS::set_pseudoprior_sigmasq() {
  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) { 
      for(int c=0; c < nclasses[model]; c++) {
        double res= sigmasq_invS[model](c, j)/(double)runlen;
        double res2= sigmasq_loginvS[model](c, j)/(double)runlen;
        double _s  = log(res) - res2;
        J[model](j,c) = (3.0 - _s + sqrt(pow(_s-3.0, 2.0)+24.0*_s))/(12.0*_s);
        L[model](j,c) = J[model](j,c)/res;
      }
    }
  }
}

/*
void BMS::set_pseudoprior_lambda(int len, int discard) {
  stringstream ss[P[0].n_cols + P[1].n_cols];
  for(int c=0; c < P[0].n_cols; c++) {
    ss[c] << "lambda0_" << c;
  }
  for(int c=0; c < P[1].n_cols; c++) {
    ss[P[0].n_cols + c] << "lambda1_" << c;
  }

  ifstream *ifs = new ifstream[nstreams];
  for(int c=0; c < P[0].n_cols; c++) {
    ifs[pind[ss[c].str()]].open((tracedir + "/" + ss[c].str() + burnin_suffix).c_str());
  }
  for(int c=0; c < P[1].n_cols; c++) {
    ifs[pind[ss[P[0].n_cols+c].str()]].open((tracedir + "/" + ss[P[0].n_cols + c].str() + burnin_suffix).c_str());
  }
  Cube<double> temp(len-discard, y->n_rows, P[0].n_cols + P[1].n_cols); // read in for each of the two models
  double foo;
  for(int i=0; i < len; i++) {
    for(int j=0; j < y->n_rows; j++) {
      for(int c=0; c < P[0].n_cols; c++) {
        if(i < discard) {
          ifs[pind[ss[c].str()]] >> foo;
        } else {
          ifs[pind[ss[c].str()]] >> temp(i-discard,j, c);
        }
      }
      for(int c=0; c < P[1].n_cols; c++) {
        if(i < discard) {
          ifs[pind[ss[P[0].n_cols+c].str()]] >> foo;
        }  else {
          ifs[pind[ss[P[0].n_cols+c].str()]] >> temp(i-discard,j, P[0].n_cols + c);
        }
      }
    }
  }
  for(int c=0; c < P[0].n_cols; c++) {
    ifs[pind[ss[c].str()]].close();
    ifs[pind[ss[c].str()]].clear();
  }
  for(int c=0; c < P[1].n_cols; c++) {
    ifs[pind[ss[P[0].n_cols+c].str()]].close();
    ifs[pind[ss[P[0].n_cols+c].str()]].clear();
  }
  delete [] ifs;

  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) {
      for(int l=0; l < P[model].n_cols; l++) {
        if(Pnil[model]) continue;
        vector<double> res=meanvar(temp.slice(P[0].n_cols*model + l).col(j));
        if(res[1] < 0.00001) {
          cerr << "Warning: no variance in burnin for lambda " 
               << model << " " << l << " " << j << ".\n";
        } 
        S_inv[model](j, l) = res[0];
//        cerr << " lambda pseudoprior[" <<  model << "](" << j << ", " << c << ") ~ Exp(S=" << 1.0/res[0] << ")" << endl;
      }
    }
  }
}
*/

/*void BMS::set_pseudoprior_lambda() {
  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) {
      for(int l=0; l < P[model].n_cols; l++) {
        if(Pnil[model]) continue;
        S_inv[model](j, l) = lambdaS[model](j, l)/(double)runlen;
      }
    }
  }
}*/

void BMS::set_pseudoprior_lambda() {
  for(int j=0; j < y->n_rows; j++) {
    for(int model=0; model < 2; model++) { 
      for(int l=0; l < P[model].n_cols; l++) {
        if(Pnil[model]) continue;
        double res= lambda_invS[model](j, l)/(double)runlen;
        double res2= lambda_loginvS[model](j, l)/(double)runlen;
        double _s  = log(res) - res2;
        D[model](j,l) = (3.0 - _s + sqrt(pow(_s-3.0, 2.0)+24.0*_s))/(12.0*_s);
        S_inv[model](j,l) = res/D[model](j,l);
      }
    }
  }
}

/*
void BMS::set_pseudoprior_rho(int len, int discard) {
  ifstream *ifs = new ifstream[nstreams];
  ifs[pind["rho0"]].open((tracedir + "/rho0" + burnin_suffix).c_str());
  ifs[pind["rho1"]].open((tracedir + "/rho1" + burnin_suffix).c_str());
  Cube<double> temp(len-discard, y->n_rows, 2); // read in for each of the two models
  Cube<double> temp2(len-discard, y->n_rows, 2); // read in for each of the two models
  double foo;
  for(int i=0; i < len; i++) {
    for(int j=0; j < y->n_rows; j++) {
      if(i < discard) {
        ifs[pind["rho0"]] >> foo;
        ifs[pind["rho1"]] >> foo;
      } else {
        ifs[pind["rho0"]] >> temp(i-discard,j, 0);
        ifs[pind["rho1"]] >> temp(i-discard,j, 1);
        temp2(i-discard,j,0)=log(temp(i-discard,j,0));
        temp2(i-discard,j,1)=log(temp(i-discard,j,1));
      }
    }
  }
  ifs[pind["rho0"]].close();
  ifs[pind["rho0"]].clear();
  ifs[pind["rho1"]].close();
  ifs[pind["rho1"]].clear();
  delete [] ifs;
  for(int model=0; model < 2; model++) {
    for(int j=0; j < y->n_rows; j++) {
      vector<double> res=meanvar(temp.slice(model).col(j));
      vector<double> res2=meanvar(temp2.slice(model).col(j));
      double _s  = log(res[0]) - res2[0];
      Q[model](j) = (3.0 - _s + sqrt(pow(_s-3.0, 2.0)+24.0*_s))/(12.0*_s);
      R[model](j) = Q[model](j)/res[0];
//      cerr << " rho pseudoprior[" <<  model << "](" << j << ") ~ Gamma(Q=" << Q[model](j) << ", R=" << R[model](j) << "(rate))" << endl;
    }
  }
}
*/

void BMS::set_pseudoprior_rho() {
  for(int model=0; model < 2; model++) {
    for(int j=0; j < y->n_rows; j++) {
      double res=rhoS[model](j)/(double)runlen;
      double res2=rho_logS[model](j)/(double)runlen;
      double _s = log(res) - res2;
      Q[model](j) = (3.0 - _s + sqrt(pow(_s-3.0, 2.0)+24.0*_s))/(12.0*_s);
      R[model](j) = Q[model](j)/res;
    }
  }
}

// probably doesn't work with non-standard Ps
void BMS::copy_pseudo() { // debug function: copy pseudo params from model 0 to model 1; warning: P0 must be equal to P1
  if(P[0].n_cols != P[1].n_cols) {
    cerr << "Error: copy_pseudo requires P0 = P1.\n";
    exit(1);
  }
  umat res = (P[0] == P[1]);
  res = res.t() * res;
  res = res.t() * res;
  if( conv_to<int>::from(res) == 0) {
    cerr << "Error: copy_pseudo requires P0 = P1.\n";
    exit(1);
  }
  for(int j=0; j < y->n_rows; j++) {
    Q[1](j) = Q[0](j);
    R[1](j) = R[0](j);
    for(int c=0; c < P[0].n_cols; c++) {
      S_inv[1](j, c)=S_inv[0](j, c);
      J[1](j,c) = J[0](j,c);
      L[1](j,c) = L[0](j,c);
      F[1](j, c) = F[0](j, c);
      Veta[1](j, c) = Veta[0](j, c);
    }
    Valpha[1](j) = Valpha[0](j);
    A[1](j) = A[0](j);
  }
}

void BMS::update_alpha(int feature, int model, bool fit) {
  if(fixalpha) return;
  if(fit || gamma(feature)==model) {
    Mb[OMP_GET_THREAD_NUM]=M*beta[model].col(feature);
    Peta[OMP_GET_THREAD_NUM]=P[model]*eta[model].col(feature);
    double V=0;
    for(int j=0; j < y->n_cols; j++) {
      V += 1.0/(esq(feature,j) + sigmasq[model](C(j,model), feature));
    }
    V += 1.0/v_alpha;
    V = 1.0/V;

    double sum=0;
    for(int i=0; i < y->n_cols; i++) {
      sum += ((*y)(feature,i) - Mb[OMP_GET_THREAD_NUM](i)
        - Peta[OMP_GET_THREAD_NUM](i))/(esq(feature,i) + sigmasq[model](C(i,model), feature));
    }
    alpha[model](feature) = gsl_ran_gaussian(rg[OMP_GET_THREAD_NUM], sqrt(V)) + V*sum;
  } else {
    alpha[model](feature) = gsl_ran_gaussian(rg[OMP_GET_THREAD_NUM], sqrt(Valpha[model](feature))) + A[model](feature);
  } 
  if(recording) {
    if(fit || gamma(feature)==model) {
      alphaS[model](feature) += alpha[model](feature);
      alphaSS[model](feature) += alpha[model](feature)*alpha[model](feature);
      alphaN[model](feature) += 1;
    }
  }
}

void BMS::update_beta(int feature, int model, bool fit) {
  if(Mnil) return;
  if(fit || gamma(feature)==model) {
    Peta[OMP_GET_THREAD_NUM]=P[model]*eta[model].col(feature);
    for(int j=0; j < y->n_cols; j++) {
      invE[OMP_GET_THREAD_NUM](j,j)=1.0/esq(feature, j) + 1.0/sigmasq[model](C(j,model), feature);
    }

    V[OMP_GET_THREAD_NUM] = inv(M.t()*invE[OMP_GET_THREAD_NUM]*M + eye(M.n_cols, M.n_cols)/v_beta);
    Vchol[OMP_GET_THREAD_NUM] = trans(chol(V[OMP_GET_THREAD_NUM]));

    for(int i=0; i < M.n_cols; i++) {
      cholSim[OMP_GET_THREAD_NUM](i) = gsl_ran_gaussian(rg[OMP_GET_THREAD_NUM], 1);
    }
    beta[model].col(feature) = (Vchol[OMP_GET_THREAD_NUM]*cholSim[OMP_GET_THREAD_NUM] +
                       (V[OMP_GET_THREAD_NUM]*M.t()*invE[OMP_GET_THREAD_NUM])*
                       (y->row(feature).t() - alpha[model](feature)*ones(y->n_cols)
                         - Peta[OMP_GET_THREAD_NUM])).t();
  } else {
    for(int i=0; i < M.n_cols; i++) {
      beta[model](i, feature) =  gsl_ran_gaussian(rg[OMP_GET_THREAD_NUM], sqrt(Vbeta[model](feature,i))) + B[model](feature,i);
    }
  }
  if(recording) {
    for(int i=0; i < M.n_cols; i++) {
      betaS[model](i, feature) += beta[model](i, feature);
      betaSS[model](i, feature) += beta[model](i, feature)*beta[model](i, feature);
    }
  }
}

void BMS::update_eta(int feature, int model, bool fit) {
  if(Pnil[model]) return;
  if(fit || gamma(feature)==model) {
    Mb[OMP_GET_THREAD_NUM]=M*beta[model].col(feature);
    Peta[OMP_GET_THREAD_NUM]=P[model]*eta[model].col(feature);
    Col<double> V(P[model].n_cols);
    V.fill(0.0);

    for(int l=0; l < P[model].n_cols; l++) {
      V(l) = 1.0/lambda[model](feature,l);
      for(int i=0; i < y->n_cols; i++) {
        V(l) += P[model](i,l)*P[model](i,l)/(esq(feature,i) + sigmasq[model](C(i,model), feature));
      }
    }
    for(int l=0; l < P[model].n_cols; l++) {
      V(l) = 1.0/V(l);
    }
    Col<double> sum(P[model].n_cols);
    sum.fill(0.0);
    Col<double>  Pev;
    for(int l=0; l < P[model].n_cols; l++) {
      Pev = Peta[OMP_GET_THREAD_NUM] - P[model].col(l)*eta[model](l, feature);
      for(int i=0; i < y->n_cols; i++) {
        sum(l) += P[model](i,l) * ((*y)(feature,i) - Mb[OMP_GET_THREAD_NUM](i)
                - alpha[model](feature) - Pev(i) )/
                  (esq(feature,i) + sigmasq[model](C(i,model), feature));
      }
    }
    for(int l=0; l < P[model].n_cols; l++) {
      eta[model](l, feature) = gsl_ran_gaussian(rg[OMP_GET_THREAD_NUM], sqrt(V(l))) + V(l) * sum(l);
    }
  } else {
    for(int l=0; l < P[model].n_cols; l++) {
      eta[model](l, feature) = gsl_ran_gaussian(rg[OMP_GET_THREAD_NUM], sqrt(Veta[model](feature, l))) + F[model](feature, l);
    }
  }
  if(recording) {
    for(int l=0; l < P[model].n_cols; l++) {
      if(fit || gamma(feature)==model) {
        etaS[model](l,feature) += eta[model](l,feature);
        etaSS[model](l,feature) += eta[model](l,feature)*eta[model](l,feature);
        etaN[model](l,feature) += 1;
      }
    }
  }
}

void BMS::update_lambda(int feature, int model, bool fit) {
  if(Pnil[model]) return;
  if(fit || gamma(feature)==model) {
    for(int l=0; l < P[model].n_cols; l++) {
      //lambda[model](feature, l) = 1.0 / rnd_inv_gauss(rg[OMP_GET_THREAD_NUM], SQRT2*sqrt(s)/abs(eta[model](l, feature)), 2.0*s);
      lambda[model](feature, l) = 1.0 / gsl_ran_gamma(rg[OMP_GET_THREAD_NUM], d+.5, 1.0/( s + .5*eta[model](l, feature)*eta[model](l, feature)));
    }
  } else {
    for(int l=0; l < P[model].n_cols; l++) {
      //lambda[model](feature, l) = gsl_ran_exponential(rg[OMP_GET_THREAD_NUM], S_inv[model](feature,l));
      lambda[model](feature, l) = 1.0 / gsl_ran_gamma(rg[OMP_GET_THREAD_NUM], D[model](feature,l), S_inv[model](feature,l));
    }
  }
  if(recording) {
    double temp;
    for(int l=0; l < P[model].n_cols; l++) {
      temp = 1.0/lambda[model](feature,l);
      lambda_invS[model](feature,l) += temp;
      lambda_loginvS[model](feature,l) += log(temp);
    }
  }
}

void BMS::update_sigmasq(int feature, int model, bool fit, bool rec_ar) {
  if(fit || gamma(feature)==model) {
    Mb[OMP_GET_THREAD_NUM]=M*beta[model].col(feature);
    Col<double> sum(nclasses[model]);
    sum.fill(0.0);
    Col<double> logsqrtratio(nclasses[model]);
    logsqrtratio.fill(0.0);
    Col<double> sigmasqprop(nclasses[model]);

    double temp;
    for(int i=0; i < y->n_cols; i++) {
      temp = (*y)(feature,i) - alpha[model](feature) - Mb[OMP_GET_THREAD_NUM](i)
                   - eta[model](C(i,model), feature);
      sum(C(i,model)) += temp*temp - esq(feature,i);
    }
    
    for(int c=0; c < nclasses[model]; c++) {
      sigmasqprop(c) = 1.0/gsl_ran_gamma(rg[OMP_GET_THREAD_NUM], n[model](c)/2.0 + k/2.0,
        1.0/(max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0));
    }

    for(int i=0; i < y->n_cols; i++) {
      logsqrtratio(C(i,model)) += 
          log(esq(feature,i) + sigmasq[model](C(i,model), feature)) - log(sigmasq[model](C(i,model), feature)) - 
          log(esq(feature,i) + sigmasqprop(C(i,model))) + log(sigmasqprop(C(i,model)));
    }

    for(int c=0; c < nclasses[model]; c++) {
      logsqrtratio(c) = .5*logsqrtratio(c);
      double sumcur=0.0, sumprop=0.0, logar;

      if(sum(c)>0) {
        for(int i=0; i < y->n_cols; i++) {
          if(C(i,model)==c) {
            sumcur += esq(feature,i)/sigmasq[model](C(i,model), feature) * 
                      (pow(((*y)(feature,i) - alpha[model](feature) - Mb[OMP_GET_THREAD_NUM](i)
                        - eta[model](C(i,model), feature)), 2)/(esq(feature,i) + sigmasq[model](C(i,model), feature) ) - 1.0);
            sumprop += esq(feature,i)/sigmasqprop(C(i,model)) * 
                      (pow(((*y)(feature,i) - alpha[model](feature) - Mb[OMP_GET_THREAD_NUM](i)
                        - eta[model](C(i,model), feature)), 2)/(esq(feature,i) + sigmasqprop(C(i,model)) ) - 1.0);
          }
        }
        logar=logsqrtratio(c) + .5 *(sumprop - sumcur);
      } else {
        for(int i=0; i < y->n_cols; i++) {
          if(C(i,model)==c) {
            sumcur +=  
                      pow(((*y)(feature,i) - alpha[model](feature) - Mb[OMP_GET_THREAD_NUM](i)
                        - eta[model](C(i,model), feature)), 2)/(esq(feature,i) + sigmasq[model](C(i,model), feature) );
            sumprop+=  
                      pow(((*y)(feature,i) - alpha[model](feature) - Mb[OMP_GET_THREAD_NUM](i)
                        - eta[model](C(i,model), feature)), 2)/(esq(feature,i) + sigmasqprop(C(i,model)) );
          }
        }
        logar=logsqrtratio(c) + .5 *(sumcur - sumprop);
      }

      if(rec_ar) _sigprop[model](feature,c) += 1;
      if(logar > log(gsl_rng_uniform(rg[OMP_GET_THREAD_NUM])) ) {
        if(c==1 && model==1 && feature==367) {
          cerr << "YES " << sigmasq[model](c, feature) << " --> " << sigmasqprop(c) << "\n";
          double olds= sigmasq[model](c,feature);
          cerr << "      " << log_target_posterior(feature, model) << " --> ";
          sigmasq[model](c,feature)=sigmasqprop(c);
          cerr << log_target_posterior(feature, model) << endl;
          cerr << "Proposal params shape, scale: " << n[model](c)/2.0 + k/2.0 << ", " << max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0 << endl;
          cerr << log_inv_gamma_pdf(olds, n[model](c)/2.0 + k/2.0, max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0) << " <-- "
               << log_inv_gamma_pdf(sigmasqprop(c), n[model](c)/2.0 + k/2.0, max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0) << endl;
          cerr << "c(";
          double maxs=-1;
          double ml = -100000000000;
          for(double s=0.001; s <= 0.2501; s+=0.001) {
            sigmasq[model](c, feature)=s;
            double lp=log_target_posterior(feature, model);
            if(ml < lp) {
              ml = lp;
              maxs=s;
            }
            cerr << lp ;
            if(s < 0.250) cerr << ",";
          }
          cerr << ")" << endl;
          cerr << "Max LP " << ml << " at " << maxs << endl;
        }
        sigmasq[model](c, feature) = sigmasqprop(c);
        if(rec_ar) _sigacc[model](feature,c) += 1;
      } else {
       if(c==1 && model==1 && feature==367) {
          cerr << "NO to " << sigmasq[model](c,feature) << " --> " << sigmasqprop(c) << "\n";
          double olds= sigmasq[model](c,feature);
          cerr << "      " << log_target_posterior(feature, model) << " --> ";
          sigmasq[model](c,feature)=sigmasqprop(c);
          cerr << log_target_posterior(feature, model) << endl;

          cerr << "Proposal params shape, scale: " << n[model](c)/2.0 + k/2.0 << ", " << max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0 << endl;

          cerr << log_inv_gamma_pdf(olds, n[model](c)/2.0 + k/2.0, max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0) << " <-- "
               << log_inv_gamma_pdf(sigmasqprop(c), n[model](c)/2.0 + k/2.0, max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0) << endl;
          cerr << "c(";
          double maxs=-1;
          double ml = -100000000000;
          for(double s=0.001; s < 0.2501; s+=0.001) {
            sigmasq[model](c, feature)=s;
            double lp=log_target_posterior(feature, model);
            if(ml < lp) {
              ml = lp;
              maxs=s;
            }
            cerr << lp ;
            if(s < 0.250) cerr << ",";
          }
          cerr << ")" << endl;
          maxs= -1;
          ml = -1000000000000;
          
          cerr << "c(";
          for(double s=0.001; s < 0.2501; s+=0.001) {
            sigmasq[model](c, feature)=s;
            double lp=log_inv_gamma_pdf(s, n[model](c)/2.0 + k/2.0, max(0.0, .5 * sum(c)) + k*rho[model](feature)/2.0) ;
            if(ml < lp) {
              ml = lp;
              maxs=s;
            }
            cerr << lp ;
            if(s < 0.250) cerr << ",";
          }
          cerr << ")" << endl;
       

          cerr << "Max LP " << ml << " at " << maxs << endl;
          sigmasq[model](c, feature)=olds;
        } 
      } 
    }
  } else {
    for(int c=0; c < nclasses[model]; c++) {
      sigmasq[model](c, feature) = 1.0/gsl_ran_gamma(rg[OMP_GET_THREAD_NUM], J[model](feature, c), 1.0/L[model](feature, c) );
    }
  }
}

void BMS::update_sigmasq_RW(int feature, int model, bool fit, bool rec_ar) {
  if(fit || gamma(feature)==model) {
    Mb[OMP_GET_THREAD_NUM]=M*beta[model].col(feature);
    Peta[OMP_GET_THREAD_NUM]=P[model]*eta[model].col(feature);
    Col<double> logsigmasq(nclasses[model]);
    Col<double> sigmasqprop(nclasses[model]);
    Col<double> logsigmasqprop(nclasses[model]);

    for(int c=0; c < nclasses[model]; c++) {
      logsigmasq(c) = log(sigmasq[model](c, feature));
      logsigmasqprop(c) = gsl_ran_gaussian(rg[OMP_GET_THREAD_NUM], g) + logsigmasq(c);
      sigmasqprop(c) = exp(logsigmasqprop(c));
    }

    double temp;
    Col<double> sum(nclasses[model]);
    sum.fill(0.0);
    for(int i=0; i < y->n_cols; i++) {
      temp = (*y)(feature,i) - alpha[model](feature) - Mb[OMP_GET_THREAD_NUM](i)
                   - Peta[OMP_GET_THREAD_NUM](i);
      sum(C(i,model)) += log(esq(feature,i) + sigmasqprop(C(i,model))) - log(esq(feature,i) + sigmasq[model](C(i,model), feature)) +
            temp*temp*
              (1.0/(esq(feature,i) + sigmasqprop(C(i,model))) - 1.0/(esq(feature,i) + sigmasq[model](C(i,model),feature)));
    }

    for(int c=0; c < nclasses[model]; c++) {
      double logar = -.5*sum(c) - .5*k*rho[model](feature)*(1.0/sigmasqprop(c) - 1.0/sigmasq[model](c, feature))
              - .5*k*(logsigmasqprop(c) - logsigmasq(c));

      if(rec_ar) _sigprop[model](feature,c) += 1;
      if(logar > log(gsl_rng_uniform(rg[OMP_GET_THREAD_NUM])) ) {
        sigmasq[model](c, feature) = sigmasqprop(c);
        if(rec_ar) _sigacc[model](feature,c) += 1;
      } 
    }
  } else {
    for(int c=0; c < nclasses[model]; c++) {
      sigmasq[model](c, feature) = 1.0/gsl_ran_gamma(rg[OMP_GET_THREAD_NUM], J[model](feature, c), 1.0/L[model](feature, c) );
    }
  }
  if(recording) {
    double temp;
    for(int c=0; c < nclasses[model]; c++) {
      temp = 1.0/sigmasq[model](c, feature);
      sigmasq_invS[model](c, feature) += temp;
      sigmasq_loginvS[model](c, feature) += log(temp);
    }
  }
}


void BMS::update_rho(int feature, int model, bool fit) {
  if(fit || gamma(feature)==model) {
    double sum=0.0;
    for(int c=0; c < nclasses[model]; c++) {
      sum += 1.0/sigmasq[model](c, feature);
    }
    rho[model](feature)=gsl_ran_gamma(rg[OMP_GET_THREAD_NUM], nclasses[model]*.5*k + q, 1.0/(r + .5*k*sum));
  } else {
    rho[model](feature)=gsl_ran_gamma(rg[OMP_GET_THREAD_NUM], Q[model](feature), 1.0/R[model](feature));
  }
  if(recording) {
    rhoS[model](feature) += rho[model](feature);
    rho_logS[model](feature) += log(rho[model](feature));
  }
}

double BMS::log_target_posterior(int feature, int model) {
  double res=0.0, sum=0.0, sum2=0.0;
  Mb[OMP_GET_THREAD_NUM]=M*beta[model].col(feature);
  Peta[OMP_GET_THREAD_NUM]=P[model]*eta[model].col(feature);
  for(int i=0; i < y->n_cols; i++) {
    sum += log((esq(feature,i) + sigmasq[model](C(i,model), feature)));
    sum2 += (pow((*y)(feature,i) - alpha[model](feature) - Mb[OMP_GET_THREAD_NUM](i) -
                Peta[OMP_GET_THREAD_NUM](i), 2))/(esq(feature,i) + sigmasq[model](C(i,model), feature));
  }
  res += -.5 * sum -.5*sum2;
  res += -.5 * alpha[model](feature)*alpha[model](feature)/v_alpha;
  sum=0;
  for(int j=0; j < M.n_cols; j++) {
    sum += beta[model](j, feature)*beta[model](j, feature);
  }
  res += -.5/v_beta * sum; 
  sum=0;
  if(!Pnil[model]) {
    for(int l=0; l < P[model].n_cols; l++) {  
      sum += -.5 * eta[model](l, feature)*eta[model](l, feature)/lambda[model](feature, l)
             -(1.5+d)*log(lambda[model](feature,l)) + d*log(s) - gsl_sf_lngamma(d) - s/lambda[model](feature,l);
    }
  }
  for(int c=0; c < nclasses[model]; c++) {  
    sum += k/2.0*log(rho[model](feature)) - k*rho[model](feature)/(2.0*sigmasq[model](c, feature))
          -(1.0+k/2.0)*log(sigmasq[model](c, feature));
  }
  res += sum;
  res += k*nclasses[model]/2.0* log(k/2.0) - nclasses[model] * gsl_sf_lngamma(k/2.0);
  res += (q-1.0)*log(rho[model](feature)) - r*rho[model](feature); 
  return(res);
}

double BMS::log_target_pseudo(int feature, int model) {
  double res=0.0, sum=0.0;
  if(!fixalpha) {
    res += -.5 * log(Valpha[model](feature)) - .5*pow(alpha[model](feature) - A[model](feature), 2)/Valpha[model](feature);
  }
  sum=0;
  if(!Pnil[model]) {
    for(int l=0; l < P[model].n_cols; l++) {
      sum += -.5*log(Veta[model](feature,l)) - .5*pow((eta[model](l,feature) - F[model](feature, l)),2)/Veta[model](feature, l);
      sum += D[model](feature, l)*log(1.0/S_inv[model](feature, l)) - (D[model](feature, l)+1)*log(lambda[model](feature, l));
      sum += - gsl_sf_lngamma(D[model](feature, l)) - 1.0/(S_inv[model](feature,l)*lambda[model](feature, l));
    }
  }
  for(int c=0; c < nclasses[model]; c++) {
    sum += J[model](feature, c)*log(L[model](feature,c)) - gsl_sf_lngamma(J[model](feature,c));
    sum += (-J[model](feature, c) -1.0) * log(sigmasq[model](c, feature)) - L[model](feature,c)/sigmasq[model](c, feature);
  }
  res += sum;

  res += Q[model](feature) *log(R[model](feature)) - gsl_sf_lngamma(Q[model](feature));
  res += (Q[model](feature) -1.0)*log(rho[model](feature)) - R[model](feature)*rho[model](feature);
  
  return(res);
} 

void BMS::update_gamma(int feature) {
  double LO =  log_target_posterior(feature,1) + log_target_pseudo(feature, 0)
                   - log_target_posterior(feature,0) - log_target_pseudo(feature, 1)
                   + logitp[feature];
 
  double x=gsl_rng_uniform(rg[OMP_GET_THREAD_NUM]);
  x=log(x)-log(1-x);
  LOsum[feature]+=LO;
  if( x < LO) {
    gamma(feature)=1;
  } else {
    gamma(feature)=0;
  }
  
  if(recording) {
    gammasum(feature) += gamma(feature);
  }
}

double BMS::gammamean(int feature) {
  double res=gammasum(feature)/(double)runlen;
  return(res);
}

int BMS::nc(int model) {
  return(nclasses[model]);
}

int BMS::pc(int model) {
  return(P[model].n_cols);
}

double BMS::alphamean(int model, int feature) {
  double res=alphaS[model](feature)/(double)alphaN[model](feature);
  return(res);
}

double BMS::etamean(int model, int cov, int feature) {
  double res=etaS[model](cov, feature)/(double)etaN[model](cov, feature);
  return(res);
}

double BMS::sigar(int model, int feature, int condition) {
  return(_sigacc[model](feature, condition)/_sigprop[model](feature,condition));
}

void BMS::reset() {
  gammasum.fill(0.0);
  runlen=0;

  for(int model=0; model < 2; model++) {
    alphaS[model].fill(0.0);
    alphaSS[model].fill(0.0);
    alphaN[model].fill(0.0);
    betaS[model].fill(0.0);
    betaSS[model].fill(0.0);
    etaS[model].fill(0.0);
    etaSS[model].fill(0.0);
    etaN[model].fill(0.0);
    lambda_invS[model].fill(0.0);
    lambda_loginvS[model].fill(0.0);
    sigmasq_invS[model].fill(0.0);
    sigmasq_loginvS[model].fill(0.0);
    rhoS[model].fill(0.0);
    rho_logS[model].fill(0.0);
  }
  recording=false;
}

void BMS::tick() { // count iters
  if(recording)
    runlen++;
}

void BMS::record() { // start using burnin iters for pseudoprior
  recording=true;
}

void BMS::print_pseudo() {
  if(tracedir.size()==0) return;
  string file=(tracedir + "/pseudo");

  ofstream ofp(file.c_str());
  for(int model=0; model < 2; model++) {
    ofp << "A" << model << "\tValpha" << model << "\t";
    for(int l=0; l < P[model].n_cols; l++) {
      ofp << "F" << model << "_" << l << "\tVeta" << model << "_" << l << "\t";
      ofp << "S" << model << "_" << l << "\t";
    }
    for(int c=0; c < nclasses[model]; c++) {
      ofp << "J" << model << "_" << c << "\tL" << model << "_" << c << "\t";
    }
    ofp << "Q" << model << "\tR" << model << "\t";
  }
  ofp << endl;

  for(int i=0; i < y->n_rows; i++) {
    for(int model=0; model < 2; model++) {
      ofp << A[0](i,0) << "\t" << Valpha[0](i,0) << "\t";
      for(int l=0; l < P[model].n_cols; l++) {
        ofp << F[model](i, l) << "\t" << Veta[model](i,l) << "\t"
            << 1.0/S_inv[model](i, l) << "\t";
      }
      for(int c=0; c < nclasses[model]; c++) {
        ofp << J[model](i,c) << "\t" << L[model](i,c) << "\t";
      }
      ofp << Q[model](i) << "\t" << R[model](i) << "\t";
    }
    ofp << endl;
  }
}


BMS::BMS(const Mat<double> * _y, const Mat<double> & _e, Mat<double> _M, Mat<double> _P0, Mat<double> _P1, Mat<int> _C, 
    double pdash, double _d, double _s, string _tracedir, gsl_rng ** _rg, const int _max_threads, bool _fixalpha) : burnin_suffix("-burnin") {
  y=_y;
  M=_M;
  P.push_back(_P0);
  P.push_back(_P1);
  C=_C;
  rg=_rg;
  max_threads=_max_threads;
  esq=pow(_e,2);
  gamma.set_size(y->n_rows);
  gamma.fill(0);
  v_alpha=25;
  v_alpharoot=sqrt(v_alpha);
  d=_d;
  s=_s; // lower means more shrinkage of etas towards small values
  tracedir=_tracedir;
  fixalpha=_fixalpha;
  if(fixalpha) {
    cerr << "Fixing alpha=0.\n";
  }
  logitp.resize(y->n_rows, log(pdash) - log(1.0-pdash));
  istuned.resize(y->n_rows, false);
  LOsum.resize(y->n_rows, 0.0);
  stilltotune_=y->n_rows;
  justtuned.resize(OMP_GET_MAX_THREADS,0);
  v_beta=pow(10,3);
  v_betaroot=sqrt(v_beta);
  k=4.0;
  q=1.2; // control shrinkage of sigma^2
  r=2.0; // control shrinkage of sigma^2
 
  Pnil.resize(2);
  Mnil=M.n_cols==1 && M.max() - M.min() < 0.00001 ? true : false;
  Pnil[0]=P[0].n_cols==1 && P[0].max() - P[0].min() < 0.00001 ? true : false;
  Pnil[1]=P[1].n_cols==1 && P[1].max() - P[1].min() < 0.00001 ? true : false;
  if(M.n_cols==1 && M.max() - M.min() < 0.00001) Mnil=true;
  if(Mnil) cerr << "Skipping betas." << endl;
  if(Pnil[0]) cerr << "Skipping etas for model 0." << endl;
  if(Pnil[1]) cerr << "Skipping etas for model 1." << endl;


  Mat<double> Z = ones(y->n_cols);
  if(!Mnil) Z.insert_cols(1, M);
  if(det(trans(Z)*Z)==0) {
    cerr << "Error: singular covariate matrix " << endl << Z;
    exit(1);
  }

  gammasum.set_size(y->n_rows);
  gammasum.fill(0.0);
  runlen=0;
  recording=false;

  for(int i=0; i < OMP_GET_MAX_THREADS; i++) {
    Col<double> c;
    Mb.push_back(c);
    Mb1.push_back(c);
    Peta.push_back(c);
    cholSim.push_back(c);
    Mb[i].set_size(y->n_rows);
    Mb1[i].set_size(y->n_rows);
    Peta[i].set_size(y->n_rows);
    cholSim[i].set_size(M.n_cols);

    Mat<double> m;
    V.push_back(m);
    Vchol.push_back(m);
    invE.push_back(m);
    V[i].set_size(M.n_cols, M.n_cols);
    Vchol[i].set_size(M.n_cols, M.n_cols);
    invE[i].set_size(y->n_cols, y->n_cols);
  }

  // nclasses
  nclasses.resize(2);
  nclasses[0]=nclasses[1]=0;
  for(int model=0; model < 2; model++) {
    map<int, bool> seen;
    vector<int> _n;
    for(int i=0; i < y->n_cols; i++) {
      if(seen.find(C(i,model)) == seen.end()) {
        seen[C(i,model)]=true;
        nclasses[model]++;
        if(C(i,model)+1 > _n.size()) {
          _n.resize(C(i,model)+1, 0);
        }
        _n[C(i,model)]++;
      }
    }
    n.push_back(Col<int>(_n));
  }

  // alpha
  for(int model=0; model < 2; model++) {
    alpha.push_back(Col<double>(y->n_rows));
    alpha[model].fill(0.0);
    alphaS.push_back(Col<double>(y->n_rows));
    alphaS[model].fill(0.0);
    alphaSS.push_back(Col<double>(y->n_rows));
    alphaSS[model].fill(0.0);
    alphaN.push_back(Col<int>(y->n_rows));
    alphaN[model].fill(0.0);
    Valpha.push_back(Col<double>(y->n_rows));
    Valpha[model].fill(v_alpha);
    A.push_back(Col<double>(y->n_rows));
    A[model].fill(0.0);
  }

  // beta
  for(int model=0; model < 2; model++) {
    beta.push_back(Mat<double>(M.n_cols, y->n_rows));
    beta[model].fill(0.0);
    betaS.push_back(Mat<double>(M.n_cols, y->n_rows));
    betaS[model].fill(0.0);
    betaSS.push_back(Mat<double>(M.n_cols, y->n_rows));
    betaSS[model].fill(0.0);
  }

  // eta, lambda, sigmasq
  for(int model=0; model < 2; model++) {
    eta.push_back(Mat<double>(P[model].n_cols, y->n_rows));
    eta[model].fill(0.0);
    etaS.push_back(Mat<double>(P[model].n_cols, y->n_rows));
    etaS[model].fill(0.0);
    etaSS.push_back(Mat<double>(P[model].n_cols, y->n_rows));
    etaSS[model].fill(0.0);
    etaN.push_back(Mat<int>(P[model].n_cols, y->n_rows));
    etaN[model].fill(0.0);
    lambda.push_back(Mat<double>(y->n_rows, P[model].n_cols));
    lambda[model].fill(s/(d-1.0));
    lambda_invS.push_back(Mat<double>(y->n_rows, P[model].n_cols));
    lambda_invS[model].fill(0.0);
    lambda_loginvS.push_back(Mat<double>(y->n_rows, P[model].n_cols));
    lambda_loginvS[model].fill(0.0);
    sigmasq.push_back(Mat<double>(nclasses[model], y->n_rows));
    sigmasq[model].fill(0.5);
    sigmasq_invS.push_back(Mat<double>(nclasses[model], y->n_rows));
    sigmasq_invS[model].fill(0.0);
    sigmasq_loginvS.push_back(Mat<double>(nclasses[model], y->n_rows));
    sigmasq_loginvS[model].fill(0.0);
    _sigacc.push_back(Mat<double>(y->n_rows, nclasses[model]));
    _sigacc[model].fill(0.0);
    _sigprop.push_back(Mat<double>(y->n_rows, nclasses[model]));
    _sigprop[model].fill(0.0);

    J.push_back(Mat<double>(y->n_rows, nclasses[model]));
    J[model].fill(2.0);
    L.push_back(Mat<double>(y->n_rows, nclasses[model]));
    L[model].fill(0.5);
    F.push_back(Mat<double>(y->n_rows, P[model].n_cols));
    F[model].fill(0.0);
    Veta.push_back(Mat<double>(y->n_rows, P[model].n_cols));
    Veta[model].fill(1.0);
    B.push_back(Mat<double>(y->n_rows, M.n_cols));
    B[model].fill(0.0);
    Vbeta.push_back(Mat<double>(y->n_rows, M.n_cols));
    Vbeta[model].fill(1.0);
    D.push_back(Mat<double>(y->n_rows, P[model].n_cols));
    D[model].fill(0.0);
    S_inv.push_back(Mat<double>(y->n_rows, P[model].n_cols));
    S_inv[model].fill(1.0/s);
  }

  // rho
  for(int model=0; model < 2; model++) {
    rho.push_back(Col<double>(y->n_rows));
    rho[model].fill(.2);
    rhoS.push_back(Col<double>(y->n_rows));
    rhoS[model].fill(0.0);
    rho_logS.push_back(Col<double>(y->n_rows));
    rho_logS[model].fill(0.0);
    Q.push_back(Col<double>(y->n_rows));
    Q[model].fill(2.0);
    R.push_back(Col<double>(y->n_rows));
    R[model].fill(10.0);
  }

  this->initialise_streams(string("-burnin"));
}

double BMS::getp(int feature) {
  if(logitp[feature] >0) {
    return 1.0/(1.0+exp(-logitp[feature]));
  } else {
    return exp(logitp[feature])/(1.0+exp(logitp[feature]));
  }
}

void BMS::tunep(int feature, int k, int batchlen) {
  double mp = LOsum[feature]/(double)batchlen;
  if(mp > -LOGIT07 && mp < LOGIT07) {
    istuned[feature]=true;
    justtuned[OMP_GET_THREAD_NUM]++;
    return;
  }
  if(mp > 0) {
    logitp[feature] -= 1.0/(sqrt(2+(k/batchlen)));
  } else { 
    logitp[feature] += 1.0/(sqrt(2+(k/batchlen)));
  }
//  ofs[pind["gammap"]] << psum[feature]/(double)batchlen << " ";
  LOsum[feature]=0;
//  ofs[pind["p"]] << p[feature] << " ";
//  ofs[pind["p"]] << endl;
//  ofs[pind["gammap"]] << endl;
}

void BMS::updatetuned() {
  for(int i=0; i < OMP_GET_MAX_THREADS; i++) {
    stilltotune_ = stilltotune_ - justtuned[i];
    justtuned[i]=0;
  }
}

bool BMS::tuned(int feature) {
  return(istuned[feature]);
}

void BMS::initialise_streams(string suffix){
  if(tracedir.size()==0) return; // Only output if we have a directory for the traces
  nstreams=2+ 2*M.n_cols + 2*P[0].n_cols + 2*P[1].n_cols + nclasses[0] + nclasses[1] + 3 + 1 + 1;
  params.resize(nstreams);
  int j=0;
  for(int i=0; i < 2; i++) {
    stringstream ss;
    ss << "alpha" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int i=0; i < M.n_cols; i++) {
    stringstream ss;
    ss << "beta0_" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int i=0; i < M.n_cols; i++) {
    stringstream ss;
    ss << "beta1_" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int i=0; i < P[0].n_cols; i++) {
    stringstream ss;
    ss << "eta0_" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int i=0; i < P[1].n_cols; i++) {
    stringstream ss;
    ss << "eta1_" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int i=0; i < P[0].n_cols; i++) {
    stringstream ss;
    ss << "lambda0_" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int i=0; i < P[1].n_cols; i++) {
    stringstream ss;
    ss << "lambda1_" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int c=0; c < nclasses[0]; c++) {
    stringstream ss;
    ss << "sigmasq0_" << c;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int c=0; c < nclasses[1]; c++) {
    stringstream ss;
    ss << "sigmasq1_" << c;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  for(int i=0; i < 2; i++) {
    stringstream ss;
    ss << "rho" << i;
    params[j]=(ss.str()).c_str();pind[(ss.str()).c_str()]=j;
    j++;
  }
  params[j]="gamma"; pind["gamma"]=j;
  j++;
  params[j]="logitp"; pind["logitp"]=j;
  j++;
  params[j]="meanLO"; pind["meanLO"]=j;

  ofs = new ofstream[nstreams];
  for(int i=0; i < nstreams; i++) {
    ofs[i].open((tracedir + "/" + params[i] + suffix).c_str(), ios_base::trunc);
  }
}

