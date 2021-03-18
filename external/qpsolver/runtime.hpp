#ifndef __SRC_LIB_RUNTIME_HPP__
#define __SRC_LIB_RUNTIME_HPP__

#include "eventhandler.hpp"
#include "instance.hpp"
#include "settings.hpp"
#include "statistics.hpp"

struct Runtime {
   Instance instance;
   Settings settings;
   Statistics statistics;

   Eventhandler<Runtime&> endofiterationevent;

   Vector primal;
   Vector rowactivity;
   // Vector dual;

   Runtime(Instance& inst) : instance(inst), primal(Vector(instance.num_var)), rowactivity(Vector(instance.num_con)) {

   }
};

#endif
