#pragma once
#ifndef __SRC_LIB_BASIS_HPP__
#define __SRC_LIB_BASIS_HPP__

#include <map>
#include <vector>

#include "simplex/HFactor.h"

#include "pricing/pricing.hpp"
#include "runtime.hpp"
#include "instance.hpp"
#include "snippets.hpp"

#include "simplex/HVector.h"

inline HVector vec2hvec(const Vector& vec) {
   HVector hvec;
   hvec.setup(vec.dim);
   for (int i = 0; i < vec.num_nz; i++) {
      hvec.index[i] = vec.index[i];
      hvec.array[vec.index[i]] = vec.value[vec.index[i]];
   }
   hvec.count = vec.num_nz;
   hvec.packFlag = true;
   return hvec;
}

inline Vector& hvec2vec(const HVector& hvec, Vector& target) {
   target.reset();
   for (int i = 0; i < hvec.count; i++) {
      target.index[i] = hvec.index[i];
      target.value[target.index[i]] = hvec.array[hvec.index[i]];
   }
   target.num_nz = hvec.count;
   return target;
}

inline Vector hvec2vec(const HVector& hvec) {
   Vector vec(hvec.size);

   return hvec2vec(hvec, vec);
}

enum class BasisStatus {
   Default,
   ActiveAtLower = 1,
   ActiveAtUpper,
   ActiveAtZero,
   Inactive
};

class Basis {
   Runtime& runtime;
   HFactor basisfactor;
   unsigned int updatessinceinvert = 0;

   MatrixBase Atran;

   // indices of active constraints in basis
	std::vector<unsigned int> activeconstraintidx;

	// ids of constraints that are in the basis but not active
	// I need to extract those columns to get Z
	std::vector<unsigned int> nonactiveconstraintsidx;

	// ids of constraints that are in the basis 
	// std::vector<int> baseindex;
   int* baseindex;

   std::map<unsigned int, BasisStatus> basisstatus;

   // index i: -1 if constraint not in basis, [0, num_var] if constraint in basis (active or not)
	std::vector<int> constraintindexinbasisfactor;

   void build();
   void rebuild();

   // buffer to avoid recreating vectors
   Vector buffer_column_aq;
   Vector buffer_row_ep;

public:
   Basis(Runtime& rt, std::vector<unsigned int> active, std::vector<BasisStatus> atlower, std::vector<unsigned int> inactive);

   unsigned int getnupdatessinceinvert() {
      return updatessinceinvert;
   }

   unsigned int getnumactive() const { return activeconstraintidx.size(); };

   unsigned int getnuminactive() const { return nonactiveconstraintsidx.size(); };

   const std::vector<unsigned int>& getactive() const { return activeconstraintidx; };

   const std::vector<unsigned int>& getinactive() const { return nonactiveconstraintsidx; };

   const std::vector<int>& getindexinfactor() const { return constraintindexinbasisfactor; };

   BasisStatus getstatus(unsigned int conid) { return basisstatus[conid]; };

   void report();

   // move that constraint into V section basis (will correspond to Nullspace from now on)
   void deactivate(unsigned int conid);

   void activate(Runtime& rt, unsigned int conid, BasisStatus atlower, unsigned int nonactivetoremove, Pricing* pricing);

   void updatebasis(Runtime& rt, int newactivecon, int droppedcon, Pricing* pricing);

   Vector btran(const Vector& rhs) const;

   Vector ftran(const Vector& rhs) const;

   Vector& btran(const Vector& rhs, Vector& target) const;

   Vector& ftran(const Vector& rhs, Vector& target) const;

   Vector recomputex(const Instance& inst);

   void write(std::string filename);
};

#endif
