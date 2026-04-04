// Third derivative by vector adjoint of vector adjoint of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xxx(const X_t<T>& x_values, Y_XXX_t<T>& y_xxx) {
  for (int m1=0;m1<m;++m1) {
    for (int m2=0;m2<n;++m2) {
      // activate x and y
      X_t<A_t<A_t<A_t<T,n>,1>,1>> x; Y_t<A_t<A_t<A_t<T,n>,1>,1>> y;
      for (int i=0;i<n;++i) {
        // set x
        x(i).value().value().value()=x_values(i);
        // register x with A_t<A_t<A_t<T,n>,1>,1>::tape
        x(i).register_input(); 
        // register x with A_t<A_t<T,n>,1>::tape
        x(i).value().register_input(); 
        // register x with A_t<T,n>::tape
        x(i).value().value().register_input(); 
      }
      // run overloaded F
      F(x,y);
      // allocate vector of adjoints of A_t<A_t<A_t<T,n>,1>,1>::tape
      A_t<A_t<A_t<T,n>,1>,1>::tape::init_adjoints();
      // set Y_{(1)}
      y(m1).adjoint().value().value()=1;
      // interpret A_t<A_t<A_t<T,n>,1>,1>::tape
      A_t<A_t<A_t<T,n>,1>,1>::tape::interpret();
      // allocate vector of adjoints of A_t<A_t<T,n>,1>::tape
      A_t<A_t<T,n>,1>::tape::init_adjoints();
      // set X_{(1,2)}
      x(m2).adjoint().adjoint().value()=1;
      // interpret A_t<A_t<T,n>,1>::tape
      A_t<A_t<T,n>,1>::tape::interpret();
      // allocate vector of adjoints of A_t<T,n>::tape
      A_t<T,n>::tape::init_adjoints();
      // set X_{(1,2)}
      for (int m3=0;m3<n;++m3) 
        x(m3).value().adjoint().adjoint(m3)=1;
      // interpret A_t<T,n>::tape
      A_t<T,n>::tape::interpret();
      // extract third derivative from X_{(3)}
      for (int m3=0;m3<n;++m3) 
        for (int i3=0;i3<n;++i3) 
          y_xxx(m1)(i3)(m2)(m3)=x(i3).value().value().adjoint(m3);
      // reset tapes
      A_t<T,n>::tape::reset();
      A_t<A_t<T,n>,1>::tape::reset();
      A_t<A_t<A_t<T,n>,1>,1>::tape::reset();
    }	
  }
}
