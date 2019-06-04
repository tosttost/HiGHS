#include "Highs.h"

#include <stdio.h>

int main(void) {
   int idx[] = {0, 1};
   double row1[] = {2, 2};
   double row2[] = {-2, 2};
   double row3[] = {4, 2};

   HighsLp lp;
   Highs highs;
   highs.initializeLp(lp);
   highs.changeObjectiveSense(objSense::OBJSENSE_MAXIMIZE);
   highs.addCol(-1.0, 0, highs.options_.infinite_bound, 0, NULL, NULL);
   highs.addCol(2.0, 0, highs.options_.infinite_bound, 0, NULL, NULL);
   highs.addRow(3.0, highs.options_.infinite_bound, 2, idx, row1);
   highs.addRow(-highs.options_.infinite_bound, 3.0, 2, idx, row2);
   highs.addRow(-highs.options_.infinite_bound, 18.5, 2, idx, row3);

   highs.run();

   HighsSolution sol = highs.getSolution();
   for (int i=0; i< sol.col_value.size(); i++) {
      printf("x%d =  %lf\n", i, sol.col_value[i]);
   }
   printf("obj: %lf\n", highs.getObjectiveValue());
}