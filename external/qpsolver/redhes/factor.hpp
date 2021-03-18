#ifndef __SRC_LIB_FACTOR_HPP__
#define __SRC_LIB_FACTOR_HPP__

#include <vector>
#include <cassert>
#include "math.h"

#include "matrix.hpp"
#include "runtime.hpp"

class CholeskyFactor {
private:
   bool uptodate = false;

   Runtime& runtime;

   std::vector<std::vector<double>> orig;

   std::vector<std::vector<double>> L;

   // solve L * w = rhs
   Vector solveL(const Vector& rhs) {
      Vector res = rhs;

      for (int r=0; r<res.dim; r++) {
         double sum = rhs.value[r];

         for (int j=0; j<r; j++) {
            sum -= res.value[j]*L[r][j];
         }

         res.value[r] = sum / L[r][r];
      }

      return res;
   }

   // solve L' u = v
   Vector solveLT(const Vector& rhs) {
      Vector res = rhs;

      for (int i=rhs.dim-1; i>=0; i--) {
         double sum = 0.0;
         for (int j=rhs.dim-1; j>i; j--) {
            sum += res.value[j] * L[j][i];
         } 
         res.value[i] = (rhs.value[i] - sum) / L[i][i];
      }

      return res;
   }

   void recompute(Matrix& Z) {
      // M = (Q' * Z)' * Z
      Matrix m = Z.tran_mat(runtime.instance.Q).mat_mat(Z);


      orig.assign(m.mat.num_col, std::vector<double>(m.mat.num_col, 0.0));

      L.assign(m.mat.num_col, std::vector<double>(m.mat.num_col, 0.0));
 
      for (unsigned int i=0; i<m.mat.num_col; i++) {
         for (int j=m.mat.start[i]; j<m.mat.start[i+1]; j++) {
            int row = m.mat.index[j];
            orig[row][i] = m.mat.value[j];
         }
      }

      for (int col = 0; col < orig.size(); col++) { 
        for (int row = 0; row <= col; row++) { 
            double sum = 0; 
            if (row == col) { 
               for (int k = 0; k < row; k++) 
                    sum += L[row][k] * L[row][k];
               L[row][row] = sqrt(orig[row][row] - sum); 
            } else { 
               for (int k = 0; k < row; k++) 
                  sum += (L[col][k] * L[row][k]); 
               L[col][row] = (orig[col][row] - sum) / L[row][row]; 
            } 
         } 
      } 

      uptodate = true;
   }

public:
   CholeskyFactor(Runtime& rt) : runtime(rt) {
      uptodate = false;
   }
 
   Vector solve(Matrix& Z, const Vector& rhs) {
      if (!uptodate) {
         recompute(Z);
      }

      Vector result = solveLT(solveL(rhs));

      result.resparsify();

      return result;
   }

   void expand(Matrix& Z, Vector& yp) {
      if (!uptodate) {
         return;
      }

      // m = Z' * G * yp
      Vector m = Z.vec_mat(runtime.instance.Q.vec_mat(yp));
      // mu = yp' * G * yp
      double mu = runtime.instance.Q.vec_mat(yp) * yp;
      // solve L * l = m
      Vector l = solveL(m);
      l.resparsify();
      // lambda = mu - l'l
      double lambda = mu - l.norm2();

      assert(lambda > 0);
      
      std::vector<double> newlrow = l.value;
      newlrow.push_back(sqrt(lambda));
      L.push_back(newlrow);

      uptodate = true;
   }

   void reduce() {
      uptodate = false;

      // TODO
   }

   void report(std::string name = "") {
      printf("%s\n", name.c_str());
      for (unsigned int i=0; i<L.size(); i++) {
         for (unsigned int j=0; j<=i; j++) {
            printf("%lf ", L[i][j]);
         }
         printf("\n");
      }
   }

};

#endif
