#include "Highs.h"

#include <stdio.h>

int main(void) {
   HighsLp lp;
   HighsOptions options;
   options.filename = "qap10.mps";
   
   Filereader* reader = Filereader::getFilereader("qap10.mps");
   reader->readModelFromFile(options, lp);

   Highs highs;

   highs.options_ = options;
   highs.initializeLp(lp);
   highs.run();

   HighsSolution sol = highs.getSolution();
   for (int i=0; i< sol.col_value.size(); i++) {
      printf("x%d =  %lf\n", i, sol.col_value[i]);
   }
   printf("obj: %lf\n", highs.getObjectiveValue());
}