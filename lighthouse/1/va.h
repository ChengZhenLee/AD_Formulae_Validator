// Jacobian by vector adjoint AD

#include "F.h"
#include "ad.h"
#include <limits>


template<typename T, int size>
using A_t=ad::adjoint_t<T,size>;

template<typename T>
void va_F_x(const X_t<T>& x_values, Y_X_t<T>& y_x) {
  // activate x and y
  X_t<A_t<T,m>> x; Y_t<A_t<T,m>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value()=x_values(i);
    // register x with  A_t<T,m>::tape
    x(i).register_input(); 
  }
  // run overloaded F
  F(x,y);
  // allocate vector of adjoints of A_t<T,m>::tape
  A_t<T,m>::tape::init_adjoints();
  // set Y_{(1)}
  for (int j=0;j<m;++j) 
    y(j).adjoint(j)=1;
  // interpret A_t<T,m>::tape
  A_t<T,m>::tape::interpret();
  // extract Jacobian from X_{(1)}
  for (int j=0;j<m;++j) 
    for (int i=0;i<n;++i) 
      y_x(j)(i)=x(i).adjoint(j);

  A_t<T,m>::tape::reset();
}	


template<typename T>
void AD_F_x(const X_t<T>& x_values, Y_t<T>& y_1, X_t<T>& x_1) {
  X_t<A_t<T,m>> x; Y_t<A_t<T,m>> y;
  for (int i=0;i<n;++i) {
    x(i).value()=x_values(i);
    x(i).register_input(); 
  }

  F(x,y);

  A_t<T,m>::tape::init_adjoints();

  for (int j=0;j<m;++j) 
    y(j).adjoint(j)=y_1(j);

  A_t<T,m>::tape::interpret();

  for (int i=0;i<n;++i) {
    x_1(i) = 0;
    for (int j=0;j<m;j++) {
      x_1(i)+=x(i).adjoint(j);
    }
  }

  A_t<T,m>::tape::reset();
}


template<typename T>
void Formula_F_x(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_t<T>& y_1, X_t<T>& x_1) {
  for (int i = 0; i < n; i++) {
    x_1(i) = 0;
    for (int j = 0; j < m; j++) {
      x_1(i) += y_x(j)(i) * y_1(j);
    }
  }
}


template<typename T>
bool Validate_va() {
  X_t<T> x_values = X_t<T>::Random();

  T tol = std::sqrt(std::numeric_limits<T>::epsilon());

  Y_X_t<T> y_x;
  X_t<T> AD_x_1;
  X_t<T> Formula_x_1;
  Y_t<T> y_1 = Y_t<T>::Random();

  // Populate y_x
  va_F_x(x_values, y_x);

  // Run the AD version
  AD_F_x(x_values, y_1, AD_x_1);

  // Run the Formula version
  Formula_F_x(x_values,y_x, y_1, Formula_x_1);

  // Print the results
  std::cout << "AD X_{(1)}\n" << AD_x_1 << "\n\n";
  std::cout << "Formula X_{(1)\n" << Formula_x_1 << "\n\n";

  // Validate
  for (int i = 0; i < n; i++) {
    T diff = std::abs(Formula_x_1(i) - AD_x_1(i));
    if (diff > tol) {
      std::cout << "Validation failed at index " << i << " Diff: " << diff << "\n";
      return false;
    }
  }
  return true;
}