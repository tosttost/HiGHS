#include "Highs.h"

#include <stdio.h>

int main(void) {
   Highs highs;

   highs.readFromFile(std::string("diet.lp"));

   highs.run();

   HighsSolution sol = highs.getSolution();
   for (int i=0; i< sol.col_value.size(); i++) {
      printf("x%d =  %lf\n", i, sol.col_value[i]);
   }
   printf("obj: %lf\n", highs.getObjectiveValue());
}