#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
//#include "lp_data/HConst.h"

const bool dev_run = true;

TEST_CASE("Presolve", "[highs_test_presolve]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  const HighsModel& presolved_model = highs.getPresolvedModel();
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kNotReduced);

  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  highs.readModel(model_file);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kReduced);
  REQUIRE((presolved_model.lp_.numCol_ > 0 && presolved_model.lp_.numRow_ > 0));

  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  SpecialLps special_lps;

  special_lps.scipLpi3Lp(lp, require_model_status);
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kInfeasible);
  REQUIRE(highs.getModelStatus() == require_model_status);

  special_lps.distillationLp(lp, require_model_status, optimal_objective);
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(lp.equalButForNames(presolved_model.lp_));
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kNotReduced);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kNotset);

  special_lps.primalDualInfeasible1Lp(lp, require_model_status);
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kInfeasible);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
}
