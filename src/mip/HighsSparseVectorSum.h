/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SPARSE_VECTOR_SUM_H_
#define HIGHS_SPARSE_VECTOR_SUM_H_

#include <algorithm>
#include <cassert>
#include <vector>

#include "util/HighsCD0uble.h"
#include "util/HighsInt.h"

class HighsSparseVectorSum {
 public:
  std::vector<uint8_t> nonzeroflag;
  std::vector<HighsCD0uble> values;
  std::vector<HighsInt> nonzeroinds;
  HighsSparseVectorSum() = default;

  HighsSparseVectorSum(HighsInt dimension) { setDimension(dimension); }

  void setDimension(HighsInt dimension) {
    values.resize(dimension);
    nonzeroflag.resize(dimension);
    nonzeroinds.reserve(dimension);
  }

  void add(HighsInt index, HighsFloat value) {
    assert(index >= 0 && index < (HighsInt)nonzeroflag.size());
    if (nonzeroflag[index]) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroflag[index] = true;
      nonzeroinds.push_back(index);
    }
  }

  void add(HighsInt index, HighsCD0uble value) {
    if (nonzeroflag[index]) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroflag[index] = true;
      nonzeroinds.push_back(index);
    }
  }

  void set(HighsInt index, HighsFloat value) {
    values[index] = value;
    nonzeroinds.push_back(index);
  }

  void set(HighsInt index, HighsCD0uble value) {
    values[index] = value;
    nonzeroinds.push_back(index);
  }

  void chgValue(HighsInt index, HighsFloat val) { values[index] = val; }
  void chgValue(HighsInt index, HighsCD0uble val) { values[index] = val; }

  const std::vector<HighsInt>& getNonzeros() const { return nonzeroinds; }

  HighsFloat getValue(HighsInt index) const { return HighsFloat(values[index]); }

  void clear() {
    for (HighsInt i : nonzeroinds) nonzeroflag[i] = false;

    nonzeroinds.clear();
  }

  template <typename Pred>
  HighsInt partition(Pred&& pred) {
    return std::partition(nonzeroinds.begin(), nonzeroinds.end(), pred) -
           nonzeroinds.begin();
  }

  template <typename IsZero>
  void cleanup(IsZero&& isZero) {
    HighsInt numNz = nonzeroinds.size();

    for (HighsInt i = numNz - 1; i >= 0; --i) {
      HighsInt pos = nonzeroinds[i];
      HighsFloat val = HighsFloat(values[pos]);

      if (isZero(pos, val)) {
        values[pos] = 0.0;
        nonzeroflag[pos] = 0;
        --numNz;
        std::swap(nonzeroinds[numNz], nonzeroinds[i]);
      }
    }

    nonzeroinds.resize(numNz);
  }
};

#endif
