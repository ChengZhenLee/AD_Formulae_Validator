// Hessian by vector adjoint of vector tangent AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xx(const X_t<T>& x_values, Y_XX_t<T>& y_xx) {
  for (int n1=0;n1<n;++n1) { // shared trailing index
    // activate x and y
    X_t<T_t<A_t<T,m>,n>> x; Y_t<T_t<A_t<T,m>,n>> y;
    for (int i=0;i<n;++i) {
      // set x
      x(i).value().value()=x_values(i);
      // register x with A_t<T,m>::tape
      x(i).value().register_input();
    }
    // set X^{(1)}
    x(n1).tangent(n1).value()=1;
    // run overloaded F
    F(x,y);
    // allocate vector of adjoints of A_t<T,m>::tape
    A_t<T,m>::tape::init_adjoints();
    // set Y^{(1)}_{(2)}
    for (int m2=0;m2<m;++m2)
      y(m2).tangent(n1).adjoint(m2)=1; 
    // interpret A_t<T,m>::tape
    A_t<T,m>::tape::interpret();
    // extract (m x n)-slice of Hessian from X_{(2)}
    for (int m2=0;m2<m;++m2)
      for (int i2=0;i2<n;++i2)
        y_xx(m2)(i2)(n1)=x(i2).value().adjoint(m2);
    // plot dag of F^{(1)}
    A_t<T,m>::tape::todot("g.dot");
    // reset tape for next n1
    A_t<T,m>::tape::reset();
  }
}	
