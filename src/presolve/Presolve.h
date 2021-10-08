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
/**@file presolve/HPresolve.h
 * @brief
 */
#ifndef PRESOLVE_PRESOLVE_H_
#define PRESOLVE_PRESOLVE_H_

#include <list>
#include <map>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"
#include "presolve/HAggregator.h"
#include "presolve/HPreData.h"
#include "presolve/PresolveAnalysis.h"
#include "test/DevKkt.h"

using std::list;
using std::string;

enum class HighsPostsolveStatus {
  kNotPresolved = -1,
  kReducedSolutionDimenionsError,
  kSolutionRecovered,
  kBasisError
};

namespace presolve {

enum class Presolver {
  kMainEmpty,
  kMainRowSingletons,
  kMainForcing,
  kMainColSingletons,
  kMainD0ublet0nEq,
  kMainDominatedCols,
  kMainSingletonsOnly,
  kMainMipDualFixing,
};

const std::map<Presolver, std::string> kPresolverNames{
    {Presolver::kMainEmpty, "Empty & fixed ()"},
    {Presolver::kMainRowSingletons, "Row singletons ()"},
    {Presolver::kMainForcing, "Forcing rows ()"},
    {Presolver::kMainColSingletons, "Col singletons ()"},
    {Presolver::kMainD0ublet0nEq, "D0ublet0n eq ()"},
    {Presolver::kMainDominatedCols, "Dominated Cols()"},
    {Presolver::kMainSingletonsOnly, "Singletons only()"},
    {Presolver::kMainMipDualFixing, "Dual fixing ()"}};

class Presolve : public HPreData {
 public:
  Presolve(HighsTimer& timer_ref) : timer(timer_ref) {}
  virtual ~Presolve() {}

  HighsPresolveStatus presolve();
  HighsPostsolveStatus postsolve(const HighsSolution& reduced_solution,
                                 const HighsBasis& reduced_basis,
                                 HighsSolution& recovered_solution,
                                 HighsBasis& recovered_basis);

  HighsPostsolveStatus primalPostsolve(
      const std::vector<HighsFloat>& reduced_solution,
      HighsSolution& recovered_solution);

  void setNumericalTolerances();
  void load(const HighsLp& lp, bool mip = false);
  // todo: clear the public from below.
  string modelName;

  // Options
  std::vector<Presolver> order;

  struct AggregatorCall {
    HAggregator::PostsolveStack postsolveStack;
    std::vector<HighsFloat> colCostAtCall;
    std::vector<HighsFloat> ARvalueAtCall;
    std::vector<HighsInt> ARindexAtCall;
    std::vector<HighsInt> ARstartAtCall;
  };

  std::vector<AggregatorCall> aggregatorStack;

  HighsInt max_iterations = 0;

  void setTimeLimit(const HighsFloat limit) {
    assert(limit < inf && limit > 0);
    timer.time_limit = limit;
  }

  HighsInt iPrint = 0;
  HighsLogOptions log_options;
  HighsFloat objShift;

 private:
  HighsInt iKKTcheck = 0;
  HighsInt presolve(HighsInt print);

  const bool report_postsolve = false;

  void initializeVectors();
  void runAggregator();
  void runPropagator();
  void detectImpliedIntegers();
  void applyMipDualFixing();
  void setProblemStatus(const HighsInt s);
  void reportTimes();

  // new bounds on primal variables for implied free detection
  vector<HighsFloat> implColLower;
  vector<HighsFloat> implColUpper;
  vector<HighsInt> implColLowerRowIndex;
  vector<HighsInt> implColUpperRowIndex;

  vector<HighsInt> implRowDualLowerSingColRowIndex;
  vector<HighsInt> implRowDualUpperSingColRowIndex;

  // new bounds on row duals y_i
  vector<HighsFloat> implRowDualLower;
  vector<HighsFloat> implRowDualUpper;

  vector<HighsFloat> implColDualLower;
  vector<HighsFloat> implColDualUpper;
  vector<HighsFloat> implRowValueLower;
  vector<HighsFloat> implRowValueUpper;

  PresolveTimer timer;  // holds enum for main presolve rules

  enum Stat {
    kUnset = 0,
    kInfeasible = 1,
    kUnboundedOrInfeasible = 2,
    kOptimal = 4,
    kReduced = 5,
    kTimeout = 6,
  };

 private:
  bool mip;
  bool hasChange = true;
  HighsInt status = Stat::kUnset;

  list<HighsInt> singRow;  // singleton rows
  list<HighsInt> singCol;  // singleton columns

  // original data
 public:
  vector<HighsFloat> colCostOriginal;

 private:
  vector<HighsFloat> rowLowerOriginal;
  vector<HighsFloat> rowUpperOriginal;
  vector<HighsFloat> colLowerOriginal;
  vector<HighsFloat> colUpperOriginal;

  // functions
  void setPrimalValue(const HighsInt j, const HighsFloat value);
  void checkForChanges(HighsInt iteration);
  void resizeProblem();
  void resizeImpliedBounds();

  // easy transformations
  void removeFixedCol(HighsInt j);
  void removeEmpty();
  void removeFixed();
  void removeEmptyRow(HighsInt i);
  void removeEmptyColumn(HighsInt j);
  void removeRow(HighsInt i);
  void checkBoundsAreConsistent();

  // singleton rows
  void removeRowSingletons();
  HighsInt getSingRowElementIndexInAR(HighsInt i);
  HighsInt getSingColElementIndexInA(HighsInt j);

  // forcing constraints
  void removeForcingConstraints();
  pair<HighsFloat, HighsFloat> getImpliedRowBounds(HighsInt row);
  void setVariablesToBoundForForcingRow(const HighsInt row, const bool isLower);
  void dominatedConstraintProcedure(const HighsInt i, const HighsFloat g,
                                    const HighsFloat h);

  // HighsFloatton equations
  void removeD0ublet0nEquations();
  pair<HighsInt, HighsInt> getXYD0ublet0nEquations(const HighsInt row);
  void processRowD0ublet0nEquation(const HighsInt row, const HighsInt x,
                                   const HighsInt y, const HighsFloat akx,
                                   const HighsFloat aky, const HighsFloat b);
  pair<HighsFloat, HighsFloat> getNewBoundsD0ublet0nConstraint(HighsInt row,
                                                       HighsInt col, HighsInt j,
                                                       HighsFloat aik, HighsFloat aij);
  void UpdateMatrixCoeffD0ublet0nEquationXzero(
      const HighsInt i, const HighsInt x, const HighsInt y, const HighsFloat aiy,
      const HighsFloat akx, const HighsFloat aky);
  void UpdateMatrixCoeffD0ublet0nEquationXnonZero(
      const HighsInt i, const HighsInt x, const HighsInt y, const HighsFloat aiy,
      const HighsFloat akx, const HighsFloat aky);

  // column singletons
  void removeColumnSingletons();
  bool removeIfImpliedFree(HighsInt col, HighsInt i, HighsInt k);
  void removeFreeColumnSingleton(const HighsInt col, const HighsInt row,
                                 const HighsInt k);
  void removeZeroCostColumnSingleton(const HighsInt col, const HighsInt row,
                                     const HighsInt k);
  bool removeColumnSingletonInD0ublet0nInequality(const HighsInt col,
                                                  const HighsInt i,
                                                  const HighsInt k);
  void removeSecondColumnSingletonInD0ublet0nRow(const HighsInt j,
                                                 const HighsInt i);
  pair<HighsFloat, HighsFloat> getBoundsImpliedFree(HighsFloat lowInit, HighsFloat uppInit,
                                            const HighsInt col,
                                            const HighsInt i, const HighsInt k);
  void removeImpliedFreeColumn(const HighsInt col, const HighsInt i,
                               const HighsInt k);

  // dominated columns
  void removeDominatedColumns();
  void rowDualBoundsDominatedColumns();
  pair<HighsFloat, HighsFloat> getImpliedColumnBounds(HighsInt j);
  void removeIfWeaklyDominated(const HighsInt j, const HighsFloat d,
                               const HighsFloat e);

  //    void findDuplicateRows();
  //    void findDuplicateColumns();
  //    void removeDuplicateRows(HighsInt i, HighsInt k, HighsFloat v);
  //    HighsInt makeCheckForDuplicateRows(HighsInt k, HighsInt i,
  //    vector<HighsFloat>& coeff, vector<HighsInt>& colIndex, HighsFloat v, HighsInt
  //    whichIsFirst); void removeDuplicateColumns(HighsInt j,HighsInt k, HighsFloat
  //    v); bool checkDuplicateRows(HighsInt i, HighsInt k) ;
  //	  bool checkDuplicateColumns(HighsInt i, HighsInt k) ;

  // old or test
  // void updateRemovedColRow(HighsInt dim);
  // void updateRowsByNZ();
  void testAnAR(HighsInt post);

  void countRemovedRows(PresolveRule rule);
  void countRemovedCols(PresolveRule rule);

  HighsFloat tol = 0.0000001;
  const HighsFloat default_primal_feasiblility_tolerance = 1e-7;
  const HighsFloat default_dual_feasiblility_tolerance = 1e-7;
  const HighsFloat default_small_matrix_value = 1e-9;
  HighsFloat inconsistent_bounds_tolerance;
  HighsFloat fixed_column_tolerance;
  HighsFloat HighsFloatton_equation_bound_tolerance;
  HighsFloat HighsFloatton_inequality_bound_tolerance;
  HighsFloat presolve_small_matrix_value;
  HighsFloat empty_row_bound_tolerance;
  HighsFloat dominated_column_tolerance;
  HighsFloat weakly_dominated_column_tolerance;

  // postsolve
  bool noPostSolve = false;

  void addChange(const PresolveRule type, const HighsInt row,
                 const HighsInt col);
  void fillStackRowBounds(const HighsInt col);
  void setKKTcheckerData();

  void getBoundOnLByZj(const HighsInt row, const HighsInt j, HighsFloat* lo,
                       HighsFloat* up, const HighsFloat colLow, const HighsFloat colUpp);
  HighsFloat getRowDualPost(const HighsInt row, const HighsInt col);
  HighsFloat getColumnDualPost(const HighsInt col);
  void roundIntegerBounds(HighsInt col);
  string getDualsForcingRow(const HighsInt row, vector<HighsInt>& fRjs);
  void getDualsSingletonRow(const HighsInt row, const HighsInt col);
  void getDualsD0ublet0nEquation(const HighsInt row, const HighsInt col);
  void recordCounts(const string fileName);
  void trimA();

  void setBasisElement(const change c);

  // test basis matrix singularity
  //
  // public:
  //	vector<HighsInt> nbffull;
  //	vector<HighsInt> bindfull;
  //	void cmpNBF(HighsInt row, HighsInt col);
  //	void setNBFfullproblem(vector<HighsInt>& nbfFull, vector<HighsInt>&
  // bnFull); 	HighsInt testBasisMatrixSingularity();
  //

  // Dev presolve
  // April 2020
  void reportDevMainLoop();
  void reportDevMidMainLoop();
  PresolveStats stats;
  HighsInt runPresolvers(const std::vector<Presolver>& order);

  void checkKkt(const bool final = false);
  dev_kkt_check::State initState(const bool intermediate = false);

  void caseTwoSingletonsD0ublet0nInequality(const HighsInt row,
                                            const HighsInt x, const HighsInt y);

  // August 2020
  void removeSingletonsOnly();
  bool isKnapsack(const HighsInt col) const;
  void removeKnapsack(const HighsInt col);
};

}  // namespace presolve

#endif /* PRESOLVE_HPRESOLVE_H_ */
