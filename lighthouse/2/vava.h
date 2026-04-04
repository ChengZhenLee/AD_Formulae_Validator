// Hessian by vector adjoint of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void F_xx(const X_t<T>& x_values, Y_XX_t<T>& y_xx) {
  for (int m1=0;m1<m;++m1) { // shared trailing index
    // activate x and y
    X_t<A_t<A_t<T,n>,m>> x; Y_t<A_t<A_t<T,n>,m>> y;
    for (int i=0;i<n;++i) {
      // set x
      x(i).value().value()=x_values(i);
      // register x with A_t<T,n>::tape
      x(i).value().register_input(); 
      // register x with A_t<A_t<T,n>,m>::tape
      x(i).register_input(); 
    }
    // run overloaded F
    F(x,y);
    // allocate vector of adjoints of A_t<A_t<T,n>,m>::tape
    A_t<A_t<T,n>,m>::tape::init_adjoints();
    // set Y_{(1)}
    y(m1).adjoint(m1).value()=1;
    // interpret A_t<A_t<T,n>,m>::tape
    A_t<A_t<T,n>,m>::tape::interpret();
    // allocate vector of adjoints of A_t<T,n>::tape
    A_t<T,n>::tape::init_adjoints();
    // set X_{(1,2)}
    for (int m2=0;m2<n;++m2) 
      x(m2).adjoint(m1).adjoint(m2)=1;
    // interpret A_t<T,n>::tape
    A_t<T,n>::tape::interpret();
    // extract (n x n)-slice of Hessian from X_{(2)}
    for (int m2=0;m2<n;++m2) 
      for (int i2=0;i2<n;++i2) 
        y_xx(m1)(m2)(i2)=x(i2).value().adjoint(m2);
    // plot dag of F
    A_t<A_t<T,n>,m>::tape::todot("g1.dot");
    // plot dag of F_{(1)}
    A_t<T,n>::tape::todot("g2.dot");
    // reset both tapes for next m1
    A_t<A_t<T,n>,m>::tape::reset();
    A_t<T,n>::tape::reset();
  }
}	
