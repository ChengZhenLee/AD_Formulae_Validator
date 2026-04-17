// Third derivative by vector tangent of vector adjoint of scalar adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xxx(const X_t<T>& x_values, Y_XXX_t<T>& y_xxx) {
  for (int i1=0;i1<n;++i1) {
    // activate x and y
    X_t<T_t<A_t<T_t<T,n>,m>,1>> x; Y_t<T_t<A_t<T_t<T,n>,m>,1>> y;
    for (int i=0;i<n;++i) {
      // set x
      x(i).value().value().value()=x_values(i);
      // register x with A_t<T_t<T,n>,n>::tape
      x(i).value().register_input(); 
    }
    // set X^{(1)}
    x(i1).tangent().value().value()=1;
    // set X^{(3)}
    for (int i3=0;i3<n;++i3) 
      x(i3).value().value().tangent(i3)=1;
    // run overloaded F
    F(x,y);
    // allocate vector of adjoints of A_t<T_t<T,n>,m>::tape
    A_t<T_t<T,n>,m>::tape::init_adjoints();
    // set Y^{(1)}_{(2)}
    for (int j=0;j<m;++j) 
      y(j).tangent().adjoint(j).value()=1;
    // interpret A_t<T_t<T,n>,m>::tape
    A_t<T_t<T,n>,m>::tape::interpret();
    // extract third derivative from X^{(3)}_{(2)}
    for (int j=0;j<m;++j) 
      for (int n3=0;n3<n;++n3) 
        for (int n2=0;n2<n;++n2) 
          y_xxx(j)(i1)(n2)(n3)=x(n2).value().adjoint(j).tangent(n3);
    // reset tapes
    A_t<T_t<T,n>,m>::tape::reset();
  }
}