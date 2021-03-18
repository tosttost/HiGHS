#ifndef __SRC_LIB_FEASIBILITYHIGHS_HPP__
#define __SRC_LIB_FEASIBILITYHIGHS_HPP__

#include "Highs.h"

#include "feasibility.hpp"

CrashSolution computestartingpoint(Instance& instance) {
   // compute initial feasible point
   Highs highs;
	
	// set HiGHS to be silent
	highs.setHighsOutput(NULL);
	highs.setHighsLogfile(NULL);


	HighsLp lp;
	lp.Aindex_ = *((std::vector<int>*)&instance.A.mat.index);
	lp.Astart_ = *((std::vector<int>*)&instance.A.mat.start);
	lp.Avalue_ = instance.A.mat.value;
   lp.colCost_.assign(instance.num_var, 0.0);
	// lp.colCost_ = instance.c.value;
	lp.colLower_ = instance.var_lo;
	lp.colUpper_ = instance.var_up;
	lp.rowLower_ = instance.con_lo;
	lp.rowUpper_ = instance.con_up;
	lp.numCol_ = instance.num_var;
	lp.numRow_ = instance.num_con;
   
	highs.passModel(lp);
	highs.run();

   HighsModelStatus phase1stat = highs.getModelStatus();
   if (phase1stat == HighsModelStatus::PRIMAL_INFEASIBLE) {
      printf("QP infeasible\n");
      exit(0);
   }

	HighsSolution sol = highs.getSolution();
	HighsBasis bas = highs.getBasis();

	Vector x0(instance.num_var);
	Vector ra(instance.num_con);
	for (unsigned int i=0; i<x0.dim; i++) {
		if (fabs(sol.col_value[i]) > 10E-5) {
			x0.value[i] = sol.col_value[i];
			x0.index[x0.num_nz++] = i;
		}
	}

	for (unsigned int i=0; i<ra.dim; i++) {
		if (fabs(sol.row_value[i]) > 10E-5) {
			ra.value[i] = sol.row_value[i];
			ra.index[ra.num_nz++] = i;
		}
	}

	std::vector<unsigned int> initialactive;
   std::vector<unsigned int> initialinactive;
	std::vector<BasisStatus> atlower;
	for (int i=0; i<bas.row_status.size(); i++) {
		if (bas.row_status[i] == HighsBasisStatus::LOWER) {
			initialactive.push_back(i);
			atlower.push_back(BasisStatus::ActiveAtLower);
		} else if (bas.row_status[i] == HighsBasisStatus::UPPER) {
			initialactive.push_back(i);
			atlower.push_back(BasisStatus::ActiveAtUpper);
		} else if (bas.row_status[i] != HighsBasisStatus::BASIC) {
         printf("row %u nonbasic\n", i);
         initialinactive.push_back(instance.num_con + i);
      } else {
         assert(bas.row_status[i] == HighsBasisStatus::BASIC);
      }
	}

	for (int i=0; i<bas.col_status.size(); i++) {
		if (bas.col_status[i] == HighsBasisStatus::LOWER) {
			initialactive.push_back(i + instance.num_con);
			atlower.push_back(BasisStatus::ActiveAtLower);
		} else if (bas.col_status[i] == HighsBasisStatus::UPPER) {
			initialactive.push_back(i + instance.num_con);
			atlower.push_back(BasisStatus::ActiveAtUpper);
		} else if (bas.col_status[i] == HighsBasisStatus::ZERO) {
         // printf("col %u free and set to 0 %d\n", i, (int)bas.col_status[i]);
         initialinactive.push_back(instance.num_con + i);
      } else if (bas.col_status[i] != HighsBasisStatus::BASIC) {
         printf("Column %d basis stus %d\n", i, (int)bas.col_status[i]);
      } else {
         assert(bas.col_status[i] == HighsBasisStatus::BASIC);
      }
	}


   assert(initialactive.size() + initialinactive.size() == instance.num_var);

   CrashSolution crash(instance.num_var, instance.num_con);
   crash.rowstatus = atlower;
   crash.active = initialactive;
   crash.inactive = initialinactive;
	crash.primal = x0;
	crash.rowact = ra;

   return crash;
}

#endif
