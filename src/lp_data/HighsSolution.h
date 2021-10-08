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
/**@file lp_data/HighsSolution.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSSOLUTION_H_
#define LP_DATA_HIGHSSOLUTION_H_

#include <string>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsStatus.h"
#include "model/HighsModel.h"

class HighsLp;
struct IpxSolution;
class HighsOptions;

using std::string;

struct HighsPrimalDualErrors {
  HighsInt num_nonzero_basic_duals;
  HighsInt num_large_nonzero_basic_duals;
  HighsFloat max_nonzero_basic_dual;
  HighsFloat sum_nonzero_basic_duals;
  HighsInt num_off_bound_nonbasic;
  HighsFloat max_off_bound_nonbasic;
  HighsFloat sum_off_bound_nonbasic;
  HighsInt num_primal_residual;
  HighsFloat max_primal_residual;
  HighsFloat sum_primal_residual;
  HighsInt num_dual_residual;
  HighsFloat max_dual_residual;
  HighsFloat sum_dual_residual;
};

void getKktFailures(const HighsOptions& options, const HighsModel& model,
                    const HighsSolution& solution, const HighsBasis& basis,
                    HighsInfo& highs_info);

void getKktFailures(const HighsOptions& options, const HighsModel& model,
                    const HighsSolution& solution, const HighsBasis& basis,
                    HighsInfo& highs_info,
                    HighsPrimalDualErrors& primal_dual_errors,
                    const bool get_residuals = false);

void getLpKktFailures(const HighsOptions& options, const HighsLp& lp,
                      const HighsSolution& solution, const HighsBasis& basis,
                      HighsInfo& highs_info);

void getLpKktFailures(const HighsOptions& options, const HighsLp& lp,
                      const HighsSolution& solution, const HighsBasis& basis,
                      HighsInfo& highs_info,
                      HighsPrimalDualErrors& primal_dual_errors,
                      const bool get_residuals = false);

void getKktFailures(const HighsOptions& options, const HighsLp& lp,
                    const std::vector<HighsFloat>& gradient,
                    const HighsSolution& solution, const HighsBasis& basis,
                    HighsInfo& highs_info,
                    HighsPrimalDualErrors& primal_dual_errors,
                    const bool get_residuals = false);

void getVariableKktFailures(const HighsFloat primal_feasibility_tolerance,
                            const HighsFloat dual_feasibility_tolerance,
                            const HighsFloat lower, const HighsFloat upper,
                            const HighsFloat value, const HighsFloat dual,
                            HighsBasisStatus* status_pointer,
                            HighsFloat& primal_infeasibility,
                            HighsFloat& dual_infeasibility, HighsFloat& value_residual);

HighsFloat computeObjectiveValue(const HighsLp& lp, const HighsSolution& solution);

void refineBasis(const HighsLp& lp, const HighsSolution& solution,
                 HighsBasis& basis);

HighsStatus ipxSolutionToHighsSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<HighsFloat>& rhs, const std::vector<char>& constraint_type,
    const HighsInt ipx_num_col, const HighsInt ipx_num_row,
    const std::vector<HighsFloat>& ipx_x, const std::vector<HighsFloat>& ipx_slack_vars,
    // const std::vector<HighsFloat>& ipx_y,
    HighsSolution& highs_solution);
HighsStatus ipxBasicSolutionToHighsBasicSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<HighsFloat>& rhs, const std::vector<char>& constraint_type,
    const IpxSolution& ipx_solution, HighsBasis& highs_basis,
    HighsSolution& highs_solution);

void resetModelStatusAndHighsInfo(HighsLpSolverObject& solver_object);
void resetModelStatusAndHighsInfo(HighsModelStatus& model_status,
                                  HighsInfo& highs_info);
bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis);

bool isPrimalSolutionRightSize(const HighsLp& lp,
                               const HighsSolution& solution);
bool isDualSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);
bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);
bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis);

void clearPrimalSolutionUtil(HighsSolution& solution);
void clearDualSolutionUtil(HighsSolution& solution);
void clearSolutionUtil(HighsSolution& solution);
void clearBasisUtil(HighsBasis& solution);
#endif  // LP_DATA_HIGHSSOLUTION_H_
