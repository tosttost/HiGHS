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
/**@file presolve/HAggregator.h
 * @brief
 */
#ifndef PRESOLVE_HAGGREGATOR_H_
#define PRESOLVE_HAGGREGATOR_H_
#include <cassert>
#include <cmath>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "util/HighsCD0uble.h"
#include "util/HighsHash.h"

namespace presolve {

class HAggregator {
  // triplet storage
  std::vector<HighsFloat> Avalue;
  std::vector<HighsInt> Arow;
  std::vector<HighsInt> Acol;

  // linked list links for column based links for each nonzero
  std::vector<HighsInt> colhead;
  std::vector<HighsInt> Anext;
  std::vector<HighsInt> Aprev;

  // splay tree links for row based iteration and lookup
  std::vector<HighsInt> rowroot;
  std::vector<HighsInt> ARleft;
  std::vector<HighsInt> ARright;

  std::vector<HighsInt> rowsize;
  std::vector<HighsInt> colsize;

  std::vector<HighsInt> iterstack;
  std::vector<HighsInt> rowpositions;
  std::unordered_map<HighsInt, int> fillinCache;
  std::vector<HighsInt> impliedLbRow;
  std::vector<HighsInt> impliedUbRow;
  std::vector<HighsFloat> col_numerics_threshold;
  // priority queue to reuse free slots
  std::priority_queue<HighsInt, std::vector<HighsInt>, std::greater<HighsInt>>
      freeslots;

  // vectors holding row activities
  std::vector<HighsCD0uble> minact;
  std::vector<HighsCD0uble> maxact;
  std::vector<HighsInt> ninfmin;
  std::vector<HighsInt> ninfmax;

  struct ImpliedFreeVarReduction {
    HighsInt row;
    HighsInt col;
    HighsInt rowlen;
    HighsInt collen;
    HighsInt stackpos;
    HighsFloat eqrhs;
    HighsFloat colcost;
    HighsFloat substcoef;
  };

 public:
  struct PostsolveStack {
    friend class HAggregator;

   private:
    std::vector<std::pair<HighsInt, HighsFloat>> reductionValues;
    std::vector<ImpliedFreeVarReduction> reductionStack;

   public:
    void undo(HighsSolution& solution, HighsBasis& basis) const;

    void undo(std::vector<HighsInt>& colFlag, std::vector<HighsInt>& rowFlag,
              std::vector<HighsFloat>& col_value, std::vector<HighsFloat>& col_dual,
              std::vector<HighsFloat>& row_dual,
              std::vector<HighsBasisStatus>& col_status,
              std::vector<HighsBasisStatus>& row_status) const;

    void undo(std::vector<HighsInt>& colFlag, std::vector<HighsInt>& rowFlag,
              std::vector<HighsFloat>& colvalue) const;

    void clear() {
      reductionStack.clear();
      reductionValues.clear();
    }

    void unsetFlags(std::vector<HighsInt>& rowFlag,
                    std::vector<HighsInt>& colFlag) const {
      for (const ImpliedFreeVarReduction& reduction : reductionStack) {
        rowFlag[reduction.row] = 0;
        colFlag[reduction.col] = 0;
      }
    }

    bool empty() const { return reductionStack.empty(); }
  };

 private:
  // set with equation rows and a vector to access there iterator positions in
  // the set by index
  std::set<std::pair<HighsInt, int>> equations;
  std::vector<std::set<std::pair<HighsInt, int>>::iterator> eqiters;

  // settings used for substitution behavior
  HighsFloat drop_tolerance;
  HighsFloat bound_tolerance;
  HighsFloat markowitz_tol;
  HighsInt maxfillin;

  // references to row and column information. Row and objective information is
  // updated in the aggregator
  std::vector<HighsFloat>& rowLower;
  std::vector<HighsFloat>& rowUpper;
  std::vector<HighsFloat>& colCost;
  HighsFloat& objOffset;
  const std::vector<HighsVarType>& integrality;
  const std::vector<HighsFloat>& colLower;
  const std::vector<HighsFloat>& colUpper;

  void link(HighsInt pos);

  void unlink(HighsInt pos);

  void dropIfZero(HighsInt pos);

  HighsFloat getImpliedLb(HighsInt row, HighsInt col);

  HighsFloat getImpliedUb(HighsInt row, HighsInt col);

  bool isImpliedFree(HighsInt col);

  void computeActivities(HighsInt row);

  bool checkFillin(HighsInt row, HighsInt col);

  void substitute(PostsolveStack& postsolveStack, HighsInt row, HighsInt col);

#ifndef NDEBUG
  void debugPrintRow(HighsInt row);

  void debugPrintSubMatrix(HighsInt row, HighsInt col);
#endif

  template <typename Func>
  void loopRow(HighsInt row, Func&& func) {
    HighsInt current = rowroot[row];

    while (true) {
      while (current != -1) {
        iterstack.push_back(current);
        current = ARleft[current];
      }

      if (iterstack.empty()) return;

      if (func(iterstack.back())) {
        iterstack.clear();
        return;
      }

      current = ARright[iterstack.back()];
      iterstack.pop_back();
    }
  }

  HighsInt countFillin(HighsInt row);

  void storeRowPositions(HighsInt pos);

  HighsInt findNonzero(HighsInt row, HighsInt col);

 public:
  HAggregator(std::vector<HighsFloat>& rowLower, std::vector<HighsFloat>& rowUpper,
              std::vector<HighsFloat>& colCost, HighsFloat& objOffset,
              const std::vector<HighsVarType>& integrality,
              const std::vector<HighsFloat>& colLower,
              const std::vector<HighsFloat>& colUpper);

  void setMaxFillin(HighsInt maxfillin) { this->maxfillin = maxfillin; }

  void setDropTolerance(HighsFloat drop_tolerance) {
    this->drop_tolerance = drop_tolerance;
  }

  void setBoundTolerance(HighsFloat bound_tolerance) {
    this->bound_tolerance = bound_tolerance;
  }

  void setMarkowitzTolerance(HighsFloat markowitz_tol) {
    this->markowitz_tol = markowitz_tol;
  }

  HighsInt numNonzeros() const { return int(Avalue.size() - freeslots.size()); }

  void addNonzero(HighsInt row, HighsInt col, HighsFloat val);

  void fromCSC(const std::vector<HighsFloat>& Aval,
               const std::vector<HighsInt>& Aindex,
               const std::vector<HighsInt>& Astart);

  void fromDynamicCSC(const std::vector<HighsFloat>& Aval,
                      const std::vector<HighsInt>& Aindex,
                      const std::vector<HighsInt>& Astart,
                      const std::vector<HighsInt>& Aend,
                      const std::vector<HighsInt>& rowFlag,
                      const std::vector<HighsInt>& colFlag);

  void fromCSR(const std::vector<HighsFloat>& ARval,
               const std::vector<HighsInt>& ARindex,
               const std::vector<HighsInt>& ARstart);

  void toCSC(std::vector<HighsFloat>& Aval, std::vector<HighsInt>& Aindex,
             std::vector<HighsInt>& Astart);

  void toCSR(std::vector<HighsFloat>& ARval, std::vector<HighsInt>& ARindex,
             std::vector<HighsInt>& ARstart);

  PostsolveStack run();

  void substitute(HighsInt substcol, HighsInt staycol, HighsFloat offset,
                  HighsFloat scale);

  void removeFixedCol(HighsInt col);

  void removeRow(HighsInt row);

  void removeRedundantRows(std::vector<uint8_t>& rowdeleted);

  HighsInt strengthenInequalities();
};

}  // namespace presolve
#endif
