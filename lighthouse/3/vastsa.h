// Third derivative by vector adjoint of vector tangent of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xxx(const X_t<T>& x_values, Y_XXX_t<T>& y_xxx) {
  for (int m1=0;m1<m;++m1) {
    for (int n2=0;n2<n;++n2) {
      // activate x and y
      X_t<A_t<T_t<A_t<T,n>,1>,1>> x; Y_t<A_t<T_t<A_t<T,n>,1>,1>> y;
      for (int i=0;i<n;++i) {
        // set x
        x(i).value().value().value()=x_values(i);
        // register x with A_t<T_t<A_t<T,n>,1>,1>::tape
        x(i).register_input(); 
        // register x with A_t<T,n>::tape
        x(i).value().value().register_input(); 
      }
      // set X^{(2)}
      x(n2).value().tangent().value()=1;
      // run overloaded F
      F(x,y);
      // allocate vector of adjoints of A_t<T_t<A_t<T,n>,n>,m>::tape
      A_t<T_t<A_t<T,n>,1>,1>::tape::init_adjoints();
      // set Y_{(1)}
      y(m1).adjoint().value().value()=1;
      // interpret A_t<T_t<A_t<T,n>,n>,m>::tape
      A_t<T_t<A_t<T,n>,1>,1>::tape::interpret();
      // set X_{(1,3)}^{(2)}
      A_t<T,n>::tape::init_adjoints();
      for (int m3=0;m3<n;++m3) 
        x(m3).adjoint().tangent().adjoint(m3)=1;
      A_t<T,n>::tape::interpret();
      // set X_{(1,3)}^{(2)}
      // extract third derivative from X_{(3)}
      for (int i1=0;i1<n;++i1) 
        for (int m3=0;m3<n;++m3) 
          y_xxx(m1)(i1)(n2)(m3)=x(i1).value().value().adjoint(m3);
      A_t<T,n>::tape::reset();
      A_t<T_t<A_t<T,n>,1>,1>::tape::reset();
    }	
  }
}
