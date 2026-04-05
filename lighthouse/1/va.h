// Jacobian by vector adjoint AD

#include "F.h"
#include "ad.h"
#include <cmath>
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


template<typename T, int K>
void AD_F_x(const X_t<T>& x_values, Eigen::Matrix<T, K, m>& y_1, 
  Y_t<T>& y_values, Eigen::Matrix<T, K, n>& x_1) {
  X_t<A_t<T, K>> x; 
  Y_t<A_t<T, K>> y;

  for (int i=0;i<n;++i) {
    x(i).value()=x_values(i);
    x(i).register_input(); 
  }

  F(x,y);

  A_t<T,m>::tape::init_adjoints();

  for (int u = 0; u < K; u++) {
    for (int j = 0; j < m; j++) {
      y(j).adjoint(u) = y_1(u, j);
    }
  }

  A_t<T,m>::tape::interpret();

  // Extract primal value
  for (int j = 0; j < m; j++) {
    y_values(j) = y(j).value();
  }

  // Extract X_{(1)}
  for (int u = 0; u < K; u++) {
    for (int i = 0; i < n; i++) {
      x_1(u, i) = x(i).adjoint(u);
    }
  }

  A_t<T,m>::tape::reset();
}


template<typename T, int K>
void Formula_F_x(const X_t<T>& x_values, Y_X_t<T>& y_x, Eigen::Matrix<T, K, m>& y_1, 
  Y_t<T>& y_values, Eigen::Matrix<T, K, n>& x_1) {

  // y = f(x)
  F(x_values, y_values);

  // X_{(1)} = F' * Y_{(1)}
  for (int u = 0; u < K; u++) {
    for (int i = 0; i < n; i++) {
      x_1(u, i) = 0;
      for (int j = 0; j < m; j++) {
        x_1(u, i) += y_x(j)(i) * y_1(u, j);
      }
    }
  }
}


template<typename T, int K>
bool Validate_va() {
  T tol = std::sqrt(std::numeric_limits<T>::epsilon());

  Y_X_t<T> y_x;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, K, n> AD_x_1;
  Eigen::Matrix<T, K, n> Formula_x_1;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, K, m> y_1 = Eigen::Matrix<T, K, m>::Random().cwiseAbs();

  // Show the seed
  std::cout << "Seed for x: \n" << x_values << "\n\n";
  std::cout << "Seed for Y_({1}): \n" << y_1 << "\n\n";

  // Populate y_x
  va_F_x(x_values, y_x);

  // Run the AD version
  AD_F_x(x_values, y_1, AD_y_values, AD_x_1);

  // Run the Formula version
  Formula_F_x(x_values,y_x, y_1, Formula_y_values, Formula_x_1);

  // Print the results
  std::cout << "AD y\n" << AD_y_values << "\n\n";
  std::cout << "AD X_{(1)}\n" << AD_x_1 << "\n\n";

  std::cout << "Formula y\n" << Formula_y_values << "\n\n";
  std::cout << "Formula X_{(1)}\n" << Formula_x_1 << "\n\n";

  // Validate
  T diff;
  T maxDiff = 0;
  for (int j = 0; j < m; j++) {
    diff = std::abs(Formula_y_values(j) - AD_y_values(j));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      std::cout << "Validation for y Failed\n";
      std::cout << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }
  }

  for (int i = 0; i < n; i++) {
    diff = std::abs(Formula_x_1(i) - AD_x_1(i));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      std::cout << "Validation for X_{(1)} Failed\n";
      std::cout << "Validation failed at index " << i << " Diff: " << diff << "\n";
      return false;
    }

    maxDiff = std::max(maxDiff, diff);
  }
  std::cout << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}