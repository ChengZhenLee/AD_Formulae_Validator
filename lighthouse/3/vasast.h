// Third derivative by vector adjoint of vector tangent of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xxx(const X_t<T>& x_values, Y_XXX_t<T>& y_xxx) {
  for (int m2=0;m2<m;++m2) {
    for (int n1=0;n1<n;++n1) {
      // activate x and y
      X_t<T_t<A_t<A_t<T,n>,1>,1>> x; Y_t<T_t<A_t<A_t<T,n>,1>,1>> y;
      for (int i=0;i<n;++i) {
        // set x
        x(i).value().value().value()=x_values(i);
        // register x with A_t<A_t<T,n>,1>::tape
        x(i).value().register_input(); 
        // register x with A_t<T,n>::tape
        x(i).value().value().register_input(); 
      }
      // set X^{(2)}
      x(n1).tangent().value().value()=1;
      // run overloaded F
      F(x,y);
      // allocate vector of adjoints of A_t<A_t<T,n>,1>::tape
      A_t<A_t<T,n>,1>::tape::init_adjoints();
      // set Y^{(1)}_{(2)}
      y(m2).tangent().adjoint().value()=1;
      // interpret A_t<A_t<T,n>,1>::tape
      A_t<A_t<T,n>,1>::tape::interpret();
      // allocate vector of adjoints of A_t<T,n>::tape
      A_t<T,n>::tape::init_adjoints();
      // set X_{(2,3)}
      for (int m3=0;m3<n;++m3) 
        x(m3).value().adjoint().adjoint(m3)=1;
      A_t<T,n>::tape::interpret();
      // extract third derivative from X_{(3)}
      for (int i1=0;i1<n;++i1) 
        for (int m3=0;m3<n;++m3) 
          y_xxx(m2)(i1)(n1)(m3)=x(i1).value().value().adjoint(m3);
      A_t<T,n>::tape::reset();
      A_t<A_t<T,n>,1>::tape::reset();
    }	
  }
}
