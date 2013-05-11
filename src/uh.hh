#ifndef _UH_H
#define _UH_H

#include <omp.h>
#include <vector>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
namespace ublas = boost::numeric::ublas;
typedef ublas::compressed_matrix<bool> boolMat;
typedef boolMat::const_iterator1 boolMatConstIt1;
typedef boolMat::const_iterator2 boolMatConstIt2;

std::vector<int> uh(const boolMat & M, const boolMat & T, const std::vector<int> & k);

#endif // _UH_H
