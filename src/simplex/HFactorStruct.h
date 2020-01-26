/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HFactorStruct.h
 * @brief Basis matrix factorization, update and solves for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef HFACTORSTRUCT_H_
#define HFACTORSTRUCT_H_

struct UseSolveSparse {
  double historical_density;
  double start_density;
  double predicted_end_density;
  double end_density;
  bool use_original_logic;
  bool logic;
  bool logic_original;
  bool logic_new;
};  
  
#endif /* HFACTORSTRUCT_H_ */
