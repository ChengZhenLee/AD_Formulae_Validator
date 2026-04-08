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


template<typename T, int V>
void AD_F_x(const X_t<T>& x_values, const Eigen::Matrix<T, n, V>& x_1, 
  Y_t<T>& y_values, Eigen::Matrix<T, m, V>& y_1) {
  X_t<T_t<T, V>> x;
  Y_t<T_t<T, V>> y;
  
  for (int i = 0; i < n; ++i) {
    x(i).value() = x_values(i);
    for (int v = 0; v < V; ++v) {
      x(i).tangent(v) = x_1(i, v);
    }
  }

  F(x, y);

  for (int j = 0; j < m; ++j) {
    // Extract primal values
    y_values(j) = y(j).value();
    for (int v = 0; v < V; ++v) {
      // Extract Y_{(1)}
      y_1(j, v) = y(j).tangent(v);
    }
  }
}


template<typename T, int V>
void Formula_F_x(X_t<T>& x_values, Y_X_t<T>& y_x, const Eigen::Matrix<T, n, V>& x_1, 
  Y_t<T>& y_values, Eigen::Matrix<T, m, V>& y_1) {
  
  // y = f(x);
  F(x_values, y_values);

  // Y_{(1)} = F' * X_{(1)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V; v1++) {
      y_1(j, v1) = 0;
      for (int i = 0; i < n; i++) {
        y_1(j, v1) += y_x(j)(i) * x_1(i, v1);
      }
    }
  }
}


template <typename T, int V>
bool Validate_vt(std::ofstream& out) {

  T tol = std::sqrt(std::numeric_limits<T>::epsilon());

  Y_X_t<T> y_x;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, m, V> AD_y_1;
  Eigen::Matrix<T, m, V> Formula_y_1;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, n, V> x_1 = Eigen::Matrix<T, n, V>::Random().cwiseAbs();

  out << "\n=== Testing First Derivative (Tangent Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for X^({1}): \n" << x_1 << "\n\n";

  // Populate y_x
  vt_F_x(x_values, y_x);

  // Run the AD version
  AD_F_x(x_values, x_1, AD_y_values, AD_y_1);

  // Run the Formula version
  Formula_F_x(x_values, y_x, x_1, Formula_y_values, Formula_y_1);

  // Print the results
  out << "AD y\n" << AD_y_values << "\n\n";
  out << "AD Y^{(1)}\n" << AD_y_1 << "\n\n";

  out << "Formula y\n" << Formula_y_values << "\n\n";
  out << "Formula Y^{(1)}\n" << Formula_y_1 << "\n\n";

  // Validate
  T diff;
  T maxDiff = 0;
  for (int j = 0; j < m; j++) {
    diff = std::abs(Formula_y_values(j) - AD_y_values(j));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      out << "Validation for y Failed\n";
      out << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }

    diff = std::abs(Formula_y_1(j) - AD_y_1(j));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      out << "Validation for Y^{(1)} Failed\n";
      out << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}