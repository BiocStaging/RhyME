// Single TMB compilation unit dispatching package model templates.
#define TMB_LIB_INIT R_init_RhyME

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Winfinite-recursion"
#elif defined(__clang__)
#pragma clang diagnostic ignored "-Wignored-attributes"
#endif

#include <TMB.hpp>
#include "beta_glmm.h"
#include "beta_fe.h"

template <class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if (model == "beta_glmm") {
    return beta_glmm(this);
  } else if (model == "beta_fe") {
    return beta_fe(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
