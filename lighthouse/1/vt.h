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


template<typename T, int K>
void AD_F_x(const X_t<T>& x_values, const Eigen::Matrix<T, n, K>& x_1, 
  Y_t<T>& y_values, Eigen::Matrix<T, m, K>& y_1) {
  X_t<T_t<T, K>> x;
  Y_t<T_t<T, K>> y;
  
  for (int i = 0; i < n; ++i) {
    x(i).value() = x_values(i);
    for (int k = 0; k < K; ++k) {
      x(i).tangent(k) = x_1(i, k);
    }
  }

  F(x, y);

  for (int j = 0; j < m; ++j) {
    // Extract primal values
    y_values(j) = y(j).value();
    for (int k = 0; k < K; ++k) {
      // Extract Y_{(1)}
      y_1(j, k) = y(j).tangent(k);
    }
  }
}


template<typename T, int K>
void Formula_F_x(X_t<T>& x_values, Y_X_t<T>& y_x, const Eigen::Matrix<T, n, K>& x_1, 
  Y_t<T>& y_values, Eigen::Matrix<T, m, K>& y_1) {
  
  // y = f(x);
  F(x_values, y_values);

  // Y_{(1)} = F' * X_{(1)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < K; v1++) {
      y_1(j, v1) = 0;
      for (int i = 0; i < n; i++) {
        y_1(j, v1) += y_x(j)(i) * x_1(i, v1);
      }
    }
  }
}


template <typename T, int K>
bool Validate_vt() {
  X_t<T> x_values = X_t<T>::Random();

  T tol = std::sqrt(std::numeric_limits<T>::epsilon());

  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Y_X_t<T> y_x;
  Eigen::Matrix<T, m, K> AD_y_1;
  Eigen::Matrix<T, m, K> Formula_y_1;
  Eigen::Matrix<T, n, K> x_1 = Eigen::Matrix<T, n, K>::Random().cwiseAbs();

  // Show the seed
  std::cout << "Seed for X^({1}): \n" << x_1 << "\n\n";

  // Populate y_x
  vt_F_x(x_values, y_x);

  // Run the AD version
  AD_F_x(x_values, x_1, AD_y_values, AD_y_1);

  // Run the Formula version
  Formula_F_x(x_values, y_x, x_1, Formula_y_values, Formula_y_1);

  // Print the results
  std::cout << "AD y\n" << AD_y_values << "\n\n";
  std::cout << "AD Y^{(1)}\n" << AD_y_1 << "\n\n";

  std::cout << "Formula y\n" << Formula_y_values << "\n\n";
  std::cout << "Formula Y^{(1)}\n" << Formula_y_1 << "\n\n";

  // Validate
  T diff;
  T maxDiff = 0;
  for (int j = 0; j < m; j++) {
    diff = std::abs(Formula_y_values(j) - AD_y_values(j));
    if (diff > tol) {
      std::cout << "Validation for y Failed\n";
      std::cout << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }

    maxDiff = std::max(maxDiff, diff);

    diff = std::abs(Formula_y_1(j) - AD_y_1(j));
    if (diff > tol) {
      std::cout << "Validation for Y^{(1)} Failed\n";
      std::cout << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }

    maxDiff = std::max(maxDiff, diff);
  }
  std::cout << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}