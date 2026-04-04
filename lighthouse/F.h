#ifndef LIGHTHOUSE_F_H
#define LIGHTHOUSE_F_H

// Griewank's lighthouse function
// lighthouse at unit distance from quay wall
#include "types.h"
#include <cmath>

template <typename T>
void F(const X_t<T>& x, Y_t<T>& y) {
  T v=x(1)*x(2);
  v=tan(v);
  T w=x(0)-v;
  y(0)=v/w;
  y(1)=y(0)*x(0);
}

#endif

