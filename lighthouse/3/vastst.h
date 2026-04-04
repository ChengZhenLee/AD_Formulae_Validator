// Third derivative by vector adjoint of scalar tangent of scalar adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xxx(const X_t<T>& x_values, Y_XXX_t<T>& y_xxx) {
  for (int n2=0;n2<n;++n2) {
    for (int n1=0;n1<n;++n1) {
      // activate x and y
      X_t<T_t<T_t<A_t<T,n>,1>,1>> x; Y_t<T_t<T_t<A_t<T,n>,1>,1>> y;
      for (int i=0;i<n;++i) {
        // set x
        x(i).value().value().value()=x_values(i);
        // register x with A_t<T,n>::tape
        x(i).value().value().register_input(); 
      }
      // set X^{(1)}
      x(n1).tangent().value().value()=1;
      // set X^{(2)}
      x(n2).value().tangent().value()=1;
      // run overloaded F
      F(x,y);
      // allocate vector of adjoints of A_t<T,n>::tape
      A_t<T,n>::tape::init_adjoints();
      // set Y^{(1,2)}_{(3)}
      for (int m3=0;m3<m;++m3) 
        y(m3).tangent().tangent().adjoint(m3)=1;
      // interpret A_t<T,n>::tape
      A_t<T,n>::tape::interpret();
      // extract third derivative from X_{(3)}
      for (int i3=0;i3<n;++i3) 
        for (int m3=0;m3<m;++m3) 
          y_xxx(m3)(n1)(n2)(i3)=x(i3).value().value().adjoint(m3);
      A_t<T,n>::tape::reset();
    }	
  }
}
