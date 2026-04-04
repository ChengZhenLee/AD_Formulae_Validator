// Hessian by vector tangent of vector tangent of vector tangent AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xxx(const X_t<T>& x_values, Y_XXX_t<T>& y_xxx) {
  // activate x and y
  X_t<T_t<T_t<T_t<T,n>,n>,n>> x; Y_t<T_t<T_t<T_t<T,n>,n>,n>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value().value()=x_values(i);
    // set X^{(1)}
    x(i).tangent(i).value().value()=1;
    // set X^{(2)}
    x(i).value().tangent(i).value()=1;
    // set X^{(3)}
    x(i).value().value().tangent(i)=1;
  }
  // run overloaded F
  F(x,y);
  // extract F^{[3]} from Y^{(1,2,3)}
  for (int j=0;j<m;++j) 
    for (int n1=0;n1<n;++n1) 
      for (int n2=0;n2<n;++n2) 
        for (int n3=0;n3<n;++n3) 
          y_xxx(j)(n1)(n2)(n3)=y(j).tangent(n1).tangent(n2).tangent(n3);
}	
