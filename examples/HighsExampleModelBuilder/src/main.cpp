#include "Highs.h"

#include <stdio.h>

class FoodItem {
  public: 
   std::string commodity;
   std::string unit;
   double price;
   double calories;
   double protein;
   double calcium;
   double iron;
   double vit_a;
   double vit_b1;
   double vit_b2;
   double niacin;
   double vit_c;
};

class Nutricient {
 public:
   std::string name;
   double requirement;
};

static std::vector<FoodItem*> fooditems;
static std::vector<Nutricient*> nutricients;

void readDiet(std::string filename) {
   char buf[256];
   FILE* file = fopen(filename.c_str(), "r");
   assert(file != NULL);

   char* ret = fgets(buf, 255, file);
   while (strlen(buf) > 1 && ret != NULL) {
      Nutricient* nut = new Nutricient();
      nut->name = std::string(strtok(buf, ",")); 
      sscanf(strtok(NULL, ","), "%lf", &nut->requirement);
      
      nutricients.push_back(nut);

      ret = fgets(buf, 255, file);
   }

   fclose(file);
}

void readFoodItems(std::string filename) {

   char buf[256];
   FILE* file = fopen(filename.c_str(), "r");
   assert(file != NULL);

   char* ret = fgets(buf, 255, file);
   while (strlen(buf) > 1 && ret != NULL) {
      FoodItem* fi = new FoodItem();
      fi->commodity = std::string(strtok(buf, ","));      
      fi->unit = std::string(strtok(NULL, ","));
      sscanf(strtok(NULL, ","), "%lf", &fi->price);
      sscanf(strtok(NULL, ","), "%lf", &fi->calories);
      sscanf(strtok(NULL, ","), "%lf", &fi->protein);
      sscanf(strtok(NULL, ","), "%lf", &fi->calcium);
      sscanf(strtok(NULL, ","), "%lf", &fi->iron);
      sscanf(strtok(NULL, ","), "%lf", &fi->vit_a);
      sscanf(strtok(NULL, ","), "%lf", &fi->vit_b1);
      sscanf(strtok(NULL, ","), "%lf", &fi->vit_b2);
      sscanf(strtok(NULL, ","), "%lf", &fi->niacin);
      sscanf(strtok(NULL, ","), "%lf", &fi->vit_c);
      std::cout << "Food Item: " << fi->commodity << std::endl;
      
      fooditems.push_back(fi);

      ret = fgets(buf, 255, file);
   }
   
   fclose(file);
}

double getNutricientValue(FoodItem* fi, int nid) {
   switch(nid) {
      case 0:
         return fi->calories;
      case 1:
         return fi->protein;
      case 2:
         return fi->calcium;
      case 3:
         return fi->iron;
      case 4:
         return fi->vit_a;
      case 5:
         return fi->vit_b1;
      case 6:
         return fi->vit_b2;
      case 7: 
         return fi->niacin;
      case 8:
         return fi->vit_c;
      default:
         return 0;
   }
}

void dietproblem() {
   int nfoods = fooditems.size();
   int nnutricients = nutricients.size(); 

   HighsModelBuilder builder;
   HighsVar** vars = new HighsVar*[nfoods];

   Highs highs;
   HighsOptions options = highs.options_;

   for (int i=0; i<nfoods; i++) {
      builder.HighsCreateVar(NULL, 0.0, options.infinite_bound, 1.0, HighsVarType::CONT, &vars[i]);
   }

   for (int i=0; i<nnutricients; i++) {
      HighsLinearCons* cons;
      builder.HighsCreateLinearCons(NULL, nutricients[i]->requirement, options.infinite_bound, &cons);
      for (int j=0; j<nfoods; j++) {
         HighsLinearConsCoef* coef;
         builder.HighsCreateLinearConsCoef(vars[j], getNutricientValue(fooditems[j], i), &coef);
         builder.HighsAddLinearConsCoefToCons(cons, coef);
      }
   }

   builder.objSense = 1; 

   HighsLp lp;
   builder.HighsBuildTechnicalModel(&lp);

   highs.initializeLp(lp);
   highs.run();
   
   HighsSolution sol = highs.getSolution();
   for (int i=0; i<sol.col_value.size(); i++) {
      if (sol.col_value[i] > 10E-5){
         std::cout << "Buy " << sol.col_value[i] << " of " << fooditems[i]->commodity << std::endl;
      }
   }

   std::cout << "Daily cost: " << highs.getObjectiveValue() << std::endl;
   std::cout << "Annual cost: " << highs.getObjectiveValue()*365 << std::endl;
}

int main() {
   readDiet("nutrdat2.csv");
   readFoodItems("nutrdat1.csv");
   dietproblem();
}
