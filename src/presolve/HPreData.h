/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPreData.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_HPREDATA_H_
#define PRESOLVE_HPREDATA_H_

#include <cstring>
#include <list>
#include <stack>
#include <utility>
#include <vector>

#include "lp_data/HConst.h"
#include "test/KktChStep.h"

using std::pair;
using std::stack;
using std::string;
using std::vector;

namespace presolve {
struct change {
  int type;
  int row;
  int col;
};

class HPreData {
 public:
  HPreData();

  // Model data
  int numCol;
  int numRow;
  int numRowOriginal;
  int numColOriginal;
  int numTot;

  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;

  vector<double> colValue;
  vector<double> colDual;
  vector<double> rowValue;
  vector<double> rowDual;

  // Row wise copy of matrix.
  vector<int> ARstart;
  vector<int> ARindex;
  vector<double> ARvalue;

  vector<int> Aend;

  // Solution
  // The first numColOriginal elements are the primal variables, slacks after
  vector<double> valuePrimal;
  vector<double> valueColDual;
  vector<double> valueRowDual;

  vector<int> nzCol;  // nonzeros in columns and rows
  vector<int> nzRow;
  vector<int> flagCol;
  vector<int> flagRow;

  const bool use_simplex_basis_logic = false;  // true;//
  vector<int> nonbasicFlag;

  // Record of whether a column or row is basic or nonbasic
  vector<HighsBasisStatus> col_status;
  vector<HighsBasisStatus> row_status;

  vector<double> colCostAtEl;
  vector<double> rowLowerAtEl;
  vector<double> rowUpperAtEl;

  void print(int k);
  void printAR(int i);
  void makeARCopy();
  void makeACopy();
  double getaij(int i, int j);
  bool isZeroA(int i, int j);
  void printSolution();
  double getRowValue(int i);

  stack<double> postValue;

  // to match reduced solution to original
  vector<int> rIndex;
  vector<int> cIndex;

  KktChStep chk;

  stack<change> chng;
  stack<pair<int, vector<double>>> oldBounds;  //(j, l, u)
};

struct MainLoop {
  int rows;
  int cols;
  int nnz;
};

struct DevStats {
  int n_loops = 0;
  std::vector<MainLoop> loops;
};

struct PresolveStats {
  DevStats dev;

  int n_rows_removed = 0;
  int n_cols_removed = 0;
  int n_nnz_removed = 0;
};

void initPresolve(PresolveStats& stats);

class PresolveList {
 public:
  void setSize(const int size) { index.assign(size, -1); }

  inline int front() {
    assert(lst.size() > 0);
    return lst[0];
  }
  inline void push_back(const int S) {
    lst.push_back(S);
    index[S] = lst.size() - 1;
  }

  inline void remove(const int S) {
    // blows up assert(index.size() > S && S > 0);
    if (S < 0) return;
    // blows up assert(index[S] >= 0 && index[S] < lst.size());
    if (index[S] < 0) return;
    assert(lst[index[S]] == S);
    lst[index[S]] = lst[lst.size() - 1];
    assert(lst[lst.size() - 1] >= 0 && lst[lst.size() - 1] < index.size());
    index[lst[lst.size() - 1]] = index[S];
    index[S] = -1;

    lst.pop_back();
  }

  inline int length() { return lst.size(); }

  bool empty() {
    if (lst.size() == 0) return true;
    return false;
  }

  std::vector<int> lst;
 private:
  std::vector<int> index;
};

}  // namespace presolve

#endif /* PRESOLVE_HPREDATA_H_ */
