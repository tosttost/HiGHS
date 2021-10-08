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
/**@file simplex/HighsSimplexAnalysis.h
 * @brief Analyse simplex iterations, both for run-time control and data
 * gathering
 */
#ifndef SIMPLEX_HIGHSSIMPLEXANALYSIS_H_
#define SIMPLEX_HIGHSSIMPLEXANALYSIS_H_

#include <cassert>
#include <memory>
#include <sstream>

#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "simplex/HVector.h"
#include "simplex/SimplexConst.h"
#include "util/HighsInt.h"
#include "util/HighsTimer.h"
#include "util/HighsUtils.h"

//#ifdef OPENMP
//#include "omp.h"
//#endif

enum TRAN_STAGE {
  TRAN_STAGE_FTRAN_LOWER = 0,
  TRAN_STAGE_FTRAN_UPPER_FT,
  TRAN_STAGE_FTRAN_UPPER,
  TRAN_STAGE_BTRAN_UPPER,
  TRAN_STAGE_BTRAN_UPPER_FT,
  TRAN_STAGE_BTRAN_LOWER,
  NUM_TRAN_STAGE_TYPE,
};

struct TranStageAnalysis {
  std::string name_;
  HighsScatterData rhs_density_;
  HighsInt num_decision_;
  HighsInt num_wrong_original_sparse_decision_;
  HighsInt num_wrong_original_hyper_decision_;
  HighsInt num_wrong_new_sparse_decision_;
  HighsInt num_wrong_new_hyper_decision_;
};

const HighsInt kAnIterTraceMaxNumRec = 20;
const HighsLogType kIterationReportLogType = HighsLogType::kVerbose;

/**
 * @brief Analyse simplex iterations, both for run-time control and data
 * gathering
 */
class HighsSimplexAnalysis {
 public:
  HighsSimplexAnalysis() {}
  // Pointer to timer
  HighsTimer* timer_;

  void setup(const std::string lp_name, const HighsLp& lp,
             const HighsOptions& options,
             const HighsInt simplex_iteration_count);
  void messaging(const HighsLogOptions& log_options_);
  void iterationReport();
  void invertReport();
  void invertReport(const bool header);
  void userInvertReport(const bool force);
  void userInvertReport(const bool header, const bool force);
  bool predictEndDensity(const HighsInt tran_stage_id,
                         const HighsFloat start_density, HighsFloat& end_density);
  void afterTranStage(const HighsInt tran_stage_id, const HighsFloat start_density,
                      const HighsFloat end_density, const HighsFloat historical_density,
                      const HighsFloat predicted_end_density,
                      const bool use_solve_sparse_original_HFactor_logic,
                      const bool use_solve_sparse_new_HFactor_logic);

  void simplexTimerStart(const HighsInt simplex_clock,
                         const HighsInt thread_id = 0);
  void simplexTimerStop(const HighsInt simplex_clock,
                        const HighsInt thread_id = 0);
  bool simplexTimerRunning(const HighsInt simplex_clock,
                           const HighsInt thread_id = 0);
  HighsInt simplexTimerNumCall(const HighsInt simplex_clock,
                               const HighsInt thread_id = 0);
  HighsFloat simplexTimerRead(const HighsInt simplex_clock,
                          const HighsInt thread_id = 0);

  HighsTimerClock* getThreadFactorTimerClockPointer();

  const std::vector<HighsTimerClock>& getThreadSimplexTimerClocks() {
    return thread_simplex_clocks;
  }
  HighsTimerClock* getThreadSimplexTimerClockPtr(HighsInt i) {
    assert(i >= 0 && i < (HighsInt)thread_simplex_clocks.size());
    return &thread_simplex_clocks[i];
  }

  const std::vector<HighsTimerClock>& getThreadFactorTimerClocks() {
    return thread_factor_clocks;
  }
  HighsTimerClock* getThreadFactorTimerClockPtr(HighsInt i) {
    assert(i >= 0 && i < (HighsInt)thread_factor_clocks.size());
    return &thread_factor_clocks[i];
  }

  void iterationRecord();
  void iterationRecordMajor();
  void operationRecordBefore(const HighsInt operation_type,
                             const HVector& vector,
                             const HighsFloat historical_density);
  void operationRecordBefore(const HighsInt operation_type,
                             const HighsInt current_count,
                             const HighsFloat historical_density);
  void operationRecordAfter(const HighsInt operation_type,
                            const HVector& vector);
  void operationRecordAfter(const HighsInt operation_type,
                            const HighsInt result_count);
  void summaryReport();
  void summaryReportFactor();
  void reportSimplexTimer();
  void reportFactorTimer();
  void updateInvertFormData(const HFactor& factor);
  void reportInvertFormData();

  // Control methods to be moved to HEkkControl
  void dualSteepestEdgeWeightError(const HighsFloat computed_edge_weight,
                                   const HighsFloat updated_edge_weight);
  //  bool switchToDevex();

  std::vector<HighsTimerClock> thread_simplex_clocks;
  std::vector<HighsTimerClock> thread_factor_clocks;
  HighsTimerClock* pointer_serial_factor_clocks;

  // Local copies of LP data
  HighsInt numRow;
  HighsInt numCol;
  HighsInt numTot;
  std::string model_name_;
  std::string lp_name_;

  // Local copies of IO data
  HighsLogOptions log_options;

  // Interpreted shortcuts from bit settings in highs_analysis_level
  bool analyse_lp_data;
  bool analyse_simplex_summary_data;
  bool analyse_simplex_runtime_data;
  bool analyse_simplex_time;
  bool analyse_factor_data;
  bool analyse_factor_time;
  bool analyse_simplex_data;

  // Control parameters moving to info
  //  bool allow_dual_steepest_edge_to_devex_switch;
  //  HighsFloat dual_steepest_edge_weight_log_error_threshold;

  // Local copies of simplex data for reporting
  HighsInt simplex_strategy = 0;
  DualEdgeWeightMode edge_weight_mode = DualEdgeWeightMode::kSteepestEdge;
  HighsInt solve_phase = 0;
  HighsInt simplex_iteration_count = 0;
  HighsInt devex_iteration_count = 0;
  HighsInt pivotal_row_index = 0;
  HighsInt leaving_variable = 0;
  HighsInt entering_variable = 0;
  HighsInt rebuild_reason = 0;
  std::string rebuild_reason_string = "";
  HighsFloat reduced_rhs_value = 0;
  HighsFloat reduced_cost_value = 0;
  HighsFloat edge_weight = 0;
  HighsFloat primal_delta = 0;
  HighsFloat primal_step = 0;
  HighsFloat dual_step = 0;
  HighsFloat pivot_value_from_column = 0;
  HighsFloat pivot_value_from_row = 0;
  HighsFloat factor_pivot_threshold = 0;
  HighsFloat numerical_trouble = 0;
  HighsFloat objective_value = 0;
  HighsInt num_primal_infeasibility = 0;
  HighsInt num_dual_infeasibility = 0;
  HighsFloat sum_primal_infeasibility = 0;
  HighsFloat sum_dual_infeasibility = 0;
  // This triple is an original infeasiblility record, so it includes max,
  // but it's only used for reporting
  HighsInt num_dual_phase_1_lp_dual_infeasibility = 0;
  HighsFloat max_dual_phase_1_lp_dual_infeasibility = 0;
  HighsFloat sum_dual_phase_1_lp_dual_infeasibility = 0;
  HighsInt num_devex_framework = 0;
  HighsFloat col_aq_density;
  HighsFloat row_ep_density;
  HighsFloat row_ap_density;
  HighsFloat row_DSE_density;
  HighsFloat col_basic_feasibility_change_density;
  HighsFloat row_basic_feasibility_change_density;
  HighsFloat col_BFRT_density;
  HighsFloat primal_col_density;
  HighsFloat dual_col_density;
  HighsInt num_costly_DSE_iteration;
  HighsFloat costly_DSE_measure;

  // Local copies of parallel simplex data for reporting
  HighsInt multi_iteration_count = 0;
  HighsInt multi_chosen = 0;
  HighsInt multi_finished = 0;
  HighsInt min_threads = 0;
  HighsInt num_threads = 0;
  HighsInt max_threads = 0;

  // Unused
  //  HighsInt multi_num = 0; // Useless
  //  HighsFloat basis_condition = 0; // Maybe useful

  // Records of how pivotal row PRICE was done
  HighsInt num_col_price = 0;
  HighsInt num_row_price = 0;
  HighsInt num_row_price_with_switch = 0;

  HighsValueDistribution before_ftran_upper_sparse_density;
  HighsValueDistribution ftran_upper_sparse_density;
  HighsValueDistribution before_ftran_upper_hyper_density;
  HighsValueDistribution ftran_upper_hyper_density;
  HighsValueDistribution cost_perturbation1_distribution;
  HighsValueDistribution cost_perturbation2_distribution;
  HighsValueDistribution cleanup_dual_change_distribution;
  HighsValueDistribution cleanup_primal_step_distribution;
  HighsValueDistribution cleanup_dual_step_distribution;
  HighsValueDistribution cleanup_primal_change_distribution;

  HighsInt num_primal_cycling_detections = 0;
  HighsInt num_dual_cycling_detections = 0;

  HighsInt num_quad_chuzc = 0;
  HighsInt num_heap_chuzc = 0;
  HighsFloat sum_heap_chuzc_size = 0;
  HighsInt max_heap_chuzc_size = 0;

  HighsInt num_correct_dual_primal_flip = 0;
  HighsFloat min_correct_dual_primal_flip_dual_infeasibility = kHighsInf;
  HighsFloat max_correct_dual_primal_flip = 0;
  HighsInt num_correct_dual_cost_shift = 0;
  HighsFloat max_correct_dual_cost_shift_dual_infeasibility = 0;
  HighsFloat max_correct_dual_cost_shift = 0;
  HighsInt net_num_single_cost_shift = 0;
  HighsInt num_single_cost_shift = 0;
  HighsFloat max_single_cost_shift = 0;
  HighsFloat sum_single_cost_shift = 0;

  // Tolerances for analysis of TRAN stages - could be needed for
  // control if this is ever used again!
  vector<HighsFloat> original_start_density_tolerance;
  vector<HighsFloat> new_start_density_tolerance;
  vector<HighsFloat> historical_density_tolerance;
  vector<HighsFloat> predicted_density_tolerance;
  vector<TranStageAnalysis> tran_stage;

  std::unique_ptr<std::stringstream> analysis_log;

 private:
  void iterationReport(const bool header);
  void reportAlgorithmPhase(const bool header);
  void reportIterationObjective(const bool header);
  void reportInfeasibility(const bool header);
  void reportThreads(const bool header);
  void reportMulti(const bool header);
  void reportOneDensity(const HighsFloat density);
  void printOneDensity(const HighsFloat density);
  void reportDensity(const bool header);
  void reportInvert(const bool header);
  //  void reportCondition(const bool header);
  void reportIterationData(const bool header);
  void reportRunTime(const bool header, const HighsFloat run_time);
  void reportFreeListSize(const bool header);
  HighsInt intLog10(const HighsFloat v);
  bool dualAlgorithm();

  //  HighsFloat AnIterCostlyDseFq;  //!< Frequency of iterations when DSE is costly
  //  HighsFloat AnIterCostlyDseMeasure;

  HighsInt num_dual_steepest_edge_weight_check = 0;
  HighsInt num_dual_steepest_edge_weight_reject = 0;
  HighsInt num_wrong_low_dual_steepest_edge_weight = 0;
  HighsInt num_wrong_high_dual_steepest_edge_weight = 0;
  HighsFloat average_frequency_low_dual_steepest_edge_weight = 0;
  HighsFloat average_frequency_high_dual_steepest_edge_weight = 0;
  HighsFloat average_log_low_dual_steepest_edge_weight_error = 0;
  HighsFloat average_log_high_dual_steepest_edge_weight_error = 0;
  HighsFloat max_average_frequency_low_dual_steepest_edge_weight = 0;
  HighsFloat max_average_frequency_high_dual_steepest_edge_weight = 0;
  HighsFloat max_sum_average_frequency_extreme_dual_steepest_edge_weight = 0;
  HighsFloat max_average_log_low_dual_steepest_edge_weight_error = 0;
  HighsFloat max_average_log_high_dual_steepest_edge_weight_error = 0;
  HighsFloat max_sum_average_log_extreme_dual_steepest_edge_weight_error = 0;

  HighsInt num_invert_report_since_last_header = -1;
  HighsInt num_iteration_report_since_last_header = -1;
  HighsFloat last_user_log_time = -kHighsInf;
  HighsFloat delta_user_log_time = 1e0;

  HighsFloat average_num_threads;
  HighsFloat average_fraction_of_possible_minor_iterations_performed;
  HighsInt sum_multi_chosen = 0;
  HighsInt sum_multi_finished = 0;

  // Analysis of INVERT form
  HighsInt num_invert = 0;
  HighsInt num_kernel = 0;
  HighsInt num_major_kernel = 0;
  HighsFloat max_kernel_dim = 0;
  HighsFloat sum_kernel_dim = 0;
  HighsFloat running_average_kernel_dim = 0;
  HighsFloat sum_invert_fill_factor = 0;
  HighsFloat sum_kernel_fill_factor = 0;
  HighsFloat sum_major_kernel_fill_factor = 0;
  HighsFloat running_average_invert_fill_factor = 1;
  HighsFloat running_average_kernel_fill_factor = 1;
  HighsFloat running_average_major_kernel_fill_factor = 1;

  HighsInt AnIterIt0 = 0;
  HighsInt AnIterPrevIt;

  // Major operation analysis struct
  struct AnIterOpRec {
    HighsFloat AnIterOpHyperCANCEL;
    HighsFloat AnIterOpHyperTRAN;
    HighsInt AnIterOpRsDim;
    HighsInt AnIterOpNumCa;
    HighsInt AnIterOpNumHyperOp;
    HighsInt AnIterOpNumHyperRs;
    HighsFloat AnIterOpSumLog10RsDensity;
    HighsInt AnIterOpRsMxNNZ;
    std::string AnIterOpName;
    HighsValueDistribution AnIterOp_density;
  };
  AnIterOpRec AnIterOp[kNumSimplexNlaOperation];

  struct AnIterTraceRec {
    HighsFloat AnIterTraceTime;
    HighsFloat AnIterTraceMulti;
    HighsFloat AnIterTraceDensity[kNumSimplexNlaOperation];
    HighsFloat AnIterTraceCostlyDse;
    HighsInt AnIterTraceIter;
    HighsInt AnIterTrace_dual_edge_weight_mode;
  };

  HighsInt AnIterTraceNumRec;
  HighsInt AnIterTraceIterDl;
  AnIterTraceRec AnIterTrace[1 + kAnIterTraceMaxNumRec + 1];

  HighsInt AnIterNumInvert[kRebuildReasonCount];
  HighsInt AnIterNumEdWtIt[(HighsInt)DualEdgeWeightMode::kCount];

  HighsValueDistribution primal_step_distribution;
  HighsValueDistribution dual_step_distribution;
  HighsValueDistribution simplex_pivot_distribution;
  HighsValueDistribution numerical_trouble_distribution;
  HighsValueDistribution factor_pivot_threshold_distribution;
};

#endif /* SIMPLEX_HIGHSSIMPLEXANALYSIS_H_ */
