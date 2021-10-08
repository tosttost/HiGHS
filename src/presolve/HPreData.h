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
/**@file presolve/HPreData.h
 * @brief
 */
#ifndef PRESOLVE_HPREDATA_H_
#define PRESOLVE_HPREDATA_H_

#include <cstring>
#include <list>
#include <stack>
#include <utility>
#include <vector>

#include "lp_data/HConst.h"
#include "test/KktCh2.h"

using std::pair;
using std::stack;
using std::string;
using std::vector;

namespace presolve {
struct change {
  HighsInt type;
  HighsInt row;
  HighsInt col;
};

class HPreData {
 public:
  HPreData();
  virtual ~HPreData() = default;

  // Model data
  HighsInt numCol;
  HighsInt numRow;
  HighsInt numRowOriginal;
  HighsInt numColOriginal;
  HighsInt numTot;

  vector<HighsInt> Astart;
  vector<HighsInt> Aindex;
  vector<HighsFloat> Avalue;
  vector<HighsFloat> colCost;
  vector<HighsFloat> colLower;
  vector<HighsFloat> colUpper;
  vector<HighsFloat> rowLower;
  vector<HighsFloat> rowUpper;
  vector<HighsVarType> integrality;

  // during postsolve hold the reduced solution, then at the end of postsolve
  // they hold the recovered. passed to dev kkt checker.
  vector<HighsFloat> colValue;
  vector<HighsFloat> colDual;
  vector<HighsFloat> rowValue;
  vector<HighsFloat> rowDual;

  // Row wise copy of matrix.
  vector<HighsInt> ARstart;
  vector<HighsInt> ARindex;
  vector<HighsFloat> ARvalue;

  vector<HighsInt> Aend;

  // Solution
  // The first numColOriginal elements are the primal variables, slacks after
  vector<HighsFloat> valuePrimal;
  vector<HighsFloat> valueColDual;
  vector<HighsFloat> valueRowDual;

  vector<HighsInt> nzCol;  // nonzeros in columns and rows
  vector<HighsInt> nzRow;
  vector<HighsInt> flagCol;
  vector<HighsInt> flagRow;

  const bool use_simplex_basis_logic = false;  // true;//
  vector<HighsInt> nonbasicFlag;

  // Record of whether a column or row is basic or nonbasic
  vector<HighsBasisStatus> col_status;
  vector<HighsBasisStatus> row_status;

  vector<HighsFloat> colCostAtEl;

  void makeARCopy();
  void makeACopy();
  HighsFloat getaij(HighsInt i, HighsInt j);
  bool isZeroA(HighsInt i, HighsInt j);
  HighsFloat getRowValue(HighsInt i);

  stack<HighsFloat> postValue;

  // to match reduced solution to original
  vector<HighsInt> rIndex;
  vector<HighsInt> cIndex;

  dev_kkt_check::KktChStep chk2;

  stack<change> chng;
  stack<pair<HighsInt, vector<HighsFloat>>> oldBounds;  //(j, l, u)
};

struct MainLoop {
  HighsInt rows;
  HighsInt cols;
  HighsInt nnz;
};

struct DevStats {
  HighsInt n_loops = 0;
  std::vector<MainLoop> loops;
};

struct PresolveStats {
  DevStats dev;

  HighsInt n_rows_removed = 0;
  HighsInt n_cols_removed = 0;
  HighsInt n_nnz_removed = 0;
};

void initPresolve(PresolveStats& stats);

}  // namespace presolve

#endif /* PRESOLVE_HPREDATA_H_ */
