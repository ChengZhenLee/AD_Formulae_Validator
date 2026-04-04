// Third derivative by vector tangent of vector adjoint of scalar adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xxx(const X_t<T>& x_values, Y_XXX_t<T>& y_xxx) {
  for (int m1=0;m1<m;++m1) {
    // activate x and y
    X_t<A_t<A_t<T_t<T,n>,n>,1>> x; Y_t<A_t<A_t<T_t<T,n>,n>,1>> y;
    for (int i=0;i<n;++i) {
      // set x
      x(i).value().value().value()=x_values(i);
      // register x with A_t<A_t<T_t<T,n>,m>,1>::tape
      x(i).register_input(); 
      // register x with A_t<T_t<T,n>,m>::tape
      x(i).value().register_input(); 
      // set X^{(3)}
      x(i).value().value().tangent(i)=1;
    }
    // run overloaded F
    F(x,y);
    // allocate vector of adjoints of A_t<A_t<T_t<T,n>,m>,1>::tape
    A_t<A_t<T_t<T,n>,n>,1>::tape::init_adjoints();
    // set Y_{(1)}
    y(m1).adjoint().value().value()=1;
    // interpret A_t<A_t<T_t<T,n>,m>,1>::tape
    A_t<A_t<T_t<T,n>,n>,1>::tape::interpret();
    // allocate vector of adjoints of A_t<T_t<T,n>,m>::tape
    A_t<T_t<T,n>,n>::tape::init_adjoints();
    // set X_{(1,2)}
    for (int m2=0;m2<n;++m2) 
      x(m2).adjoint().adjoint(m2).value()=1;
    // interpret A_t<T_t<T,n>,m>::tape
    A_t<T_t<T,n>,n>::tape::interpret();
    // extract third derivative from X_{(3)}
    for (int m2=0;m2<n;++m2) 
      for (int i2=0;i2<n;++i2) 
        for (int n3=0;n3<n;++n3) 
          y_xxx(m1)(i2)(m2)(n3)=x(i2).value().adjoint(m2).tangent(n3);
    // reset tapes
    A_t<T_t<T,n>,n>::tape::reset();
    A_t<A_t<T_t<T,n>,n>,1>::tape::reset();
  }
}
