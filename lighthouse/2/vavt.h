// Hessian by vector adjoint of vector tangent AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void vavt_F_xx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx) {
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
    for (int m2=0;m2<m;++m2) {
      // populate Jacobian
      y_x(m2)(n1) = y(m2).tangent(n1).value();
      for (int i2=0;i2<n;++i2)
        y_xx(m2)(i2)(n1)=x(i2).value().adjoint(m2);
    }
    // reset tape for next n1
    A_t<T,m>::tape::reset();
  }
}	


template<typename T, int V1, int U2>
void AD_F_xx(const X_t<T>& x_values, 
  const Eigen::Matrix<T, n, V1>& x_1, const Eigen::Matrix<T, U2, m>& y_2, 
  const Eigen::Vector<Eigen::Matrix<T, m, V1>, U2>& y_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, U2, n>& x_2,
  Eigen::Vector<Eigen::Matrix<T, n, V1>, U2>& x_1_2) {

  X_t<T_t<A_t<T,U2>,V1>> x; Y_t<T_t<A_t<T,U2>,V1>> y;

  for (int v1 = 0; v1 < V1; v1++) {
    for (int i=0;i<n;++i) {
      x(i).value().value()=x_values(i);

      // set X^{(1)}
      x(i).tangent(v1).value() = x_1(i, v1);

      x(i).value().register_input();
    }
  }
    
  F(x,y);

  A_t<T,U2>::tape::init_adjoints();
  
  for (int v1 = 0; v1 < V1; v1++) {
    for (int u2 = 0; u2 < U2; u2++) {
      for (int j = 0; j < m; j++) {
        // set Y^{(1)}_{(2)}
        y(j).tangent(v1).adjoint(u2) = y_1_2(u2)(j, v1);

        // Set Y_{(2)}
        y(j).value().adjoint(u2) = y_2(u2, j);
      }
    }
  }

  A_t<T,U2>::tape::interpret();

  // Extract y (only once)
  for (int j = 0; j < m; j++) {
    y_values(j) = y(j).value().value();
  }
  
  // Extract Y^{(1)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      y_1(j, v1) = y(j).tangent(v1).value();
    }
  }

  for (int u2 = 0; u2 < U2; u2++) {
    for (int i = 0; i < n; i++) {
      // Extract X_{(2)}
      x_2(u2, i) += x(i).value().adjoint(u2);

      // Extract X^{(1)}_{(2)}
      for (int v1 = 0; v1 < V1; v1++) {
        x_1_2(u2)(i, v1) = x(i).tangent(v1).adjoint(u2);
      }
    }
  }

  // reset tape
  A_t<T,U2>::tape::reset();
}


template<typename T, int V1, int U2>
void Formula_F_xx(const Y_X_t<T>& y_x, const Y_XX_t<T>& y_xx,
  const X_t<T>& x_values, 
  const Eigen::Matrix<T, n, V1>& x_1, const Eigen::Matrix<T, U2, m>& y_2, 
  const Eigen::Vector<Eigen::Matrix<T, m, V1>, U2>& y_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, U2, n>& x_2,
  Eigen::Vector<Eigen::Matrix<T, n, V1>, U2>& x_1_2) {
  
  // y = f(x)
  F(x_values, y_values);

  // Y^{(1)} = F' * X^{(1)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      T sum = 0;
      for (int i = 0; i < n; i++) {
        sum += y_x(j)(i) * x_1(i, v1);
      }
      y_1(j, v1) = sum;
    }
  }

  // X_{(2)} = F'' * X^{(1)} * Y^{(1)}_{(2)} + F' * Y_{(2)}
  for (int u2 = 0; u2 < U2; u2++) {
    for (int i2 = 0; i2 < n; i2++) {
      T term1 = 0;
      T term2 = 0;
      for (int j = 0; j < m; j++) {
        term1 += y_x(j)(i2) * y_2(u2, j);
        for (int v1 = 0; v1 < V1; v1++) {
          for (int i1 = 0; i1 < n; i1++) {
            term2 += y_xx(j)(i1)(i2) * x_1(i1, v1) * y_1_2(u2)(j, v1);
          }
        }
      }
      x_2(u2, i2) += term1 + term2;
    }
  }

  // X^{(1)}_{(2)} = F' * Y^{(1)}_{(2)}
  for (int u2 = 0; u2 < U2; u2++) {
    for (int i = 0; i < n; i++) {
      for (int v1 = 0; v1 < V1; v1++) {
        T sum = 0;
        for (int j = 0; j < m; j++) {
          sum += y_x(j)(i) * y_1_2(u2)(j, v1);
        }
        x_1_2(u2)(i, v1) += sum;
      }
    }
  }
}

template<typename T, int V1, int U2>
bool Validate_vavt(std::ofstream& out) {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 4.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, m, V1> AD_y_1;
  Eigen::Matrix<T, m, V1> Formula_y_1;
  Eigen::Matrix<T, U2, n> AD_x_2;
  Eigen::Matrix<T, U2, n> Formula_x_2;
  Eigen::Vector<Eigen::Matrix<T, n, V1>, U2> AD_x_1_2;
  Eigen::Vector<Eigen::Matrix<T, n, V1>, U2> Formula_x_1_2;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, n, V1> x_1 = Eigen::Matrix<T, n, V1>::Random().cwiseAbs();
  Eigen::Matrix<T, U2, m> y_2 = Eigen::Matrix<T, U2, m>::Random().cwiseAbs();
  Eigen::Vector<Eigen::Matrix<T, m, V1>, U2> y_1_2;
  for (int u2 = 0; u2 < U2; u2++) {
    y_1_2(u2) = Eigen::Matrix<T, m, V1>::Random().cwiseAbs();
  }

  out << "\n=== Testing Second Derivative (Adjoint over Tangent Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for Y_({1}): \n" << x_1 << "\n\n";
  out << "Seed for X^({2}): \n" << y_2 << "\n\n";
  out << "Seed for Y_{(1)}^{(2)}: (nested matrices not printed)\n\n";

  // Populate y_x and y_xx
  vavt_F_xx(x_values, y_x, y_xx);

  // Run the AD version
  AD_F_xx(x_values, x_1, y_2, y_1_2, 
          AD_y_values, AD_y_1, AD_x_2, AD_x_1_2);

  // Run the Formula version
  Formula_F_xx(y_x, y_xx,
    x_values, x_1, y_2, y_1_2, 
    Formula_y_values, Formula_y_1, Formula_x_2, Formula_x_1_2);

  T diff;
  T maxDiff = 0;

  // Compare y
  for (int j = 0; j < m; j++) {
    diff = std::abs(AD_y_values(j) - Formula_y_values(j));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      out << "Validation for y Failed\n";
      out << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }

    // Compare Y^{(1)}
    for (int v1 = 0; v1 < V1; v1++) {
      diff = std::abs(AD_y_1(j, v1) - Formula_y_1(j, v1));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y^{(1)} Failed\n";
        out << "Validation failed at index " << j << "," << v1 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }

  // Compare X_{(2)}
  for (int u2 = 0; u2 < U2; u2++) {
    for (int i = 0; i < n; i++) {
      diff = std::abs(AD_x_2(u2, i) - Formula_x_2(u2, i));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for X_{(2)} Failed\n";
        out << "Validation failed at index " << u2 << "," << i << " Diff: " << diff << "\n";
        return false;
      }

      // Compare X^{(1)}_{(2)}
      for (int v1 = 0; v1 < V1; v1++) {
        diff = std::abs(AD_x_1_2(u2)(i, v1) - Formula_x_1_2(u2)(i, v1));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for X^{(1)}_{(2)} Failed\n";
          out << "Validation failed at index " << u2 << "," << i << "," << v1 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}