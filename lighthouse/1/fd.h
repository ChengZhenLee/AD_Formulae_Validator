/*
 * Jacobian by central finite differences
 */
#ifndef F_H
#define F_H

#include "F.h"

template<typename T>
void F_x(X_t<T> x, Y_X_t<T>& y_x) {
  // activate x and y
  Y_t<T> y_up, y_down;
  for (int i=0;i<n;++i) {
    // store original value of x(i)
    T x_i=x(i);
    // determine suitable perturbation of x(i) 
    T h=sqrt(std::numeric_limits<T>::epsilon());
    h= (x(i)!=0) ? h*fabs(x(i)) : h;
    // perturb x(i) upwards
    x(i)+=h;
    F(x,y_up);
    // perturb x(i) downwards
    x(i)=x_i-h;
    F(x,y_down);
    // reset x(i)
    x(i)=x_i;
    // approximate i-th column of the Jacobian by central finite difference
    for (int j=0;j<m;++j) 
      y_x(j)(i)=(y_up(j)-y_down(j))/(h+h);
  }
}	

#endif