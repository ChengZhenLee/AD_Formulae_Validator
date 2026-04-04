// Hessian by vector tangent of vector tangent AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void vtvt_F_xx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx) {
  // activate x and y
  X_t<T_t<T_t<T,n>,n>> x; Y_t<T_t<T_t<T,n>,n>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value().value()=x_values(i);
    // set X^{(1)}
    x(i).tangent(i).value()=1;
    // set X^{(2)}
    x(i).value().tangent(i)=1;
  }
  // run overloaded F
  F(x,y);
  // extract Hessian from Y^{(1,2)}
  for (int j=0;j<m;++j) 
    for (int n1=0;n1<n;++n1) 
      for (int n2=0;n2<n;++n2) 
        y_xx(j)(n1)(n2)=y(j).tangent(n1).tangent(n2);

  // extract Jacobian from Y^{(1)}
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      y_x(j)(i) = y(j).tangent(i).value();
    }
  }
}	