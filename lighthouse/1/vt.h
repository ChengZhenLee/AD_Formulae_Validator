// Jacobian by vector tangent AD

#include "F.h"
#include "ad.h"

template<typename T, int size>
using T_t=ad::tangent_t<T,size>;

template<typename T>
void vt_F_x(const X_t<T>& x_values, Y_X_t<T>& y_x) {
  // activate x and y
  X_t<T_t<T,n>> x; Y_t<T_t<T,n>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value()=x_values(i);
    // set X^{(1)}
    x(i).tangent(i)=1;
  }
  // run overloaded F
  F(x,y);
  // extract Jacobian from Y^{(1)}
  for (int j=0;j<m;++j) 
    for (int i=0;i<n;++i) 
      y_x(j)(i)=y(j).tangent(i);
}	


template<typename T>
void AD_F_x(const X_t<T>& x_values, X_t<T>& x_1, Y_t<T>& y_1) {
  X_t<T_t<T,n>> x; Y_t<T_t<T,n>> y;
  
  for (int i=0;i<n;++i) {
    x(i).value()=x_values(i);
    x(i).tangent(i)=x_1(i);
  }

  F(x,y);

  for (int j=0;j<m;++j) {
    y_1(j) = 0;
    for (int i=0;i<n;i++) {
      y_1(j) += y(j).tangent(i);
    }
  }
}


template<typename T>
void Formula_F_x(Y_X_t<T>& y_x, X_t<T>& x_1, Y_t<T>& y_1) {
  for (int j = 0; j < m; j++) {
    y_1(j) = y_x(j).dot(x_1);
  }
}


template <typename T>
bool Validate_vt() {
  X_t<T> x_values = X_t<T>::Random();

  T tol = std::sqrt(std::numeric_limits<T>::epsilon());

  Y_X_t<T> y_x;
  Y_t<T> AD_y_1;
  Y_t<T> Formula_y_1;
  X_t<T> x_1 = X_t<T>::Random();

  // Populate y_x
  vt_F_x(x_values, y_x);

  // Run the AD version
  AD_F_x(x_values, x_1, AD_y_1);

  // Run the Formula version
  Formula_F_x(y_x, x_1, Formula_y_1);

  // Print the results
  std::cout << "AD Y^{(1)}\n" << AD_y_1 << "\n\n";
  std::cout << "Formula Y^{(1)\n" << Formula_y_1 << "\n\n";

  // Validate
  for (int j = 0; j < m; j++) {
    T diff = std::abs(Formula_y_1(j) - AD_y_1(j));
    if (diff > tol) {
      std::cout << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }
  }  
  return true;
}