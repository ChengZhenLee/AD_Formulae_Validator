// Hessian by vector tangent of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xx(const X_t<T>& x_values, Y_XX_t<T>& y_xx) {
  // activate x and y
  X_t<A_t<T_t<T,n>,m>> x; Y_t<A_t<T_t<T,n>,m>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value()=x_values(i);
    // set X^{(2)}
    x(i).value().tangent(i)=1;
    // register x with A_t<T_t<T,n>,m>::tape
    x(i).register_input(); 
  }
  // run overloaded F
  F(x,y);
  // allocate vector of adjoints of A_t<T_t<T,n>,m>::tape
  A_t<T_t<T,n>,m>::tape::init_adjoints();
  // set Y_{(1)}
  for (int m1=0;m1<m;++m1) 
    y(m1).adjoint(m1).value()=1;
  // interpret A_t<T_t<T,n>,m>::tape
  A_t<T_t<T,n>,m>::tape::interpret();
  // extract Hessian from X_{(1)}^{(2)}
  for (int m1=0;m1<m;++m1) 
    for (int i1=0;i1<n;++i1) 
      for (int n2=0;n2<n;++n2) 
        y_xx(m1)(i1)(n2)=x(i1).adjoint(m1).tangent(n2);
}	
