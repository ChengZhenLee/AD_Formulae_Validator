// Hessian by vector tangent of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void vtva_F_xx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx) {
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

  // extract Jacobian from X_{(1)}
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      y_x(j)(i) = x(i).adjoint(j).value();
    }
  }

  A_t<T_t<T,n>,m>::tape::reset();
}	


template<typename T, int U1, int V2>
void AD_F_xx(const X_t<T>& x_values, 
  const Eigen::Matrix<T, U1, m>& y_1, const Eigen::Matrix<T, n, V2>& x_2, 
  const Eigen::Vector<Eigen::Matrix<T, m, V2>, U1>& y_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, U1, n>& x_1, Eigen::Matrix<T, m, V2>& y_2,
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1>& x_1_2) {

  X_t<A_t<T_t<T,V2>,U1>> x; Y_t<A_t<T_t<T,V2>,U1>> y;

  for (int i = 0; i < n; i++) {
    x(i).value().value() = x_values(i);

    // Set X^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      x(i).value().tangent(v2) = x_2(i, v2);
    }
    x(i).register_input(); 
  }

  F(x,y);

  A_t<T_t<T,V2>,U1>::tape::init_adjoints();

  // set Y_{(1)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int j = 0; j < m; j++) {
      y(j).adjoint(u1).value() = y_1(u1, j);
      // set Y_{(1, 2)}
      for (int v2 = 0; v2 < V2; v2++) {
        y(j).adjoint(u1).tangent(v2) = y_1_2(u1)(j, v2);
      }
    }
  }

  A_t<T_t<T,V2>,U1>::tape::interpret();

  // Extract y
  for (int j = 0; j < m; j++) {
    y_values(j) = y(j).value().value();
  }

  // Extract Y^{(2)}
  for (int j = 0; j < m; j++) {
    for (int v2 = 0; v2 < V2; v2++) {
      y_2(j, v2) = y(j).value().tangent(v2);
    }
  }

  // Extract X_{(1)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i = 0; i < n; i++) {
      x_1(u1, i) = x(i).adjoint(u1).value();

      // Extract X_{(1)}^{(2)}
      for (int v2 = 0; v2 < V2; v2++) {
        x_1_2(u1)(i, v2) = x(i).adjoint(u1).tangent(v2);
      }
    }
  }

  A_t<T_t<T,V2>,U1>::tape::reset();
}


template<typename T, int U1, int V2>
void Formula_F_xx(const Y_X_t<T>& y_x, const Y_XX_t<T>& y_xx,
  const X_t<T>& x_values, 
  const Eigen::Matrix<T, U1, m>& y_1, const Eigen::Matrix<T, n, V2>& x_2, 
  const Eigen::Vector<Eigen::Matrix<T, m, V2>, U1>& y_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, U1, n>& x_1, Eigen::Matrix<T, m, V2>& y_2,
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1>& x_1_2) {

  // y = f(x)
  F(x_values, y_values);

  // Y^{(2)} = F' * X^{(2)}
  for (int j = 0; j < m; j++) {
    for (int v2 = 0; v2 < V2; v2++) {
      T sum = 0;
      for (int i = 0; i < n; i++) {
        sum += y_x(j)(i) * x_2(i, v2); 
      }
      y_2(j, v2) = sum;
    }
  }

  // X_{(1)} = F' * Y_{(1)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i = 0; i < n; i++) {
      T sum = 0;
      for (int j = 0; j < m; j++) {
        sum += y_x(j)(i) * y_1(u1, j);
      }
      x_1(u1, i) = sum;
    }
  }

  // X_{(1)}^{(2)} = F'' * Y_{(1)} * X^{(2)} + F' * Y_{(1)}^{(2)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i1 = 0; i1 < n; i1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        T term1 = 0;
        T term2 = 0;
        for (int j = 0; j < m; j++) {
          term2 += y_x(j)(i1) * y_1_2(u1)(j, v2);
          for (int i2 = 0; i2 < n; i2++) {
            term1 += y_xx(j)(i1)(i2) * y_1(u1, j) * x_2(i2, v2);
          }
        }
        x_1_2(u1)(i1, v2) = term1 + term2;
      }
    }
  }
}


template<typename T, int U1, int V2>
bool Validate_vtva(std::ofstream& out) {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 4.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, U1, n> AD_x_1;
  Eigen::Matrix<T, U1, n> Formula_x_1;
  Eigen::Matrix<T, m, V2> AD_y_2;
  Eigen::Matrix<T, m, V2> Formula_y_2;
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1> AD_x_1_2;
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1> Formula_x_1_2;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, U1, m> y_1 = Eigen::Matrix<T, U1, m>::Random().cwiseAbs();
  Eigen::Matrix<T, n, V2> x_2 = Eigen::Matrix<T, n, V2>::Random().cwiseAbs();
  Eigen::Vector<Eigen::Matrix<T, m, V2>, U1> y_1_2;
  for (int u1 = 0; u1 < U1; u1++) {
    y_1_2(u1) = Eigen::Matrix<T, m, V2>::Random().cwiseAbs();
  }

  out << "\n=== Testing Second Derivative (Tangent over Adjoint Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for Y_({1}): \n" << y_1 << "\n\n";
  out << "Seed for X^({2}): \n" << x_2 << "\n\n";
  out << "Seed for Y_{(1)}^{(2)}: (nested matrices not printed)\n\n";

  // Populate y_x and y_xx
  vtva_F_xx(x_values, y_x, y_xx);

  // Run the AD version
  AD_F_xx(x_values, y_1, x_2, y_1_2, 
          AD_y_values, AD_x_1, AD_y_2, AD_x_1_2);

  // Run the Formula version
  Formula_F_xx(y_x, y_xx,
    x_values, y_1, x_2, y_1_2, 
    Formula_y_values, Formula_x_1, Formula_y_2, Formula_x_1_2);

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
    // Compare Y^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      diff = std::abs(AD_y_2(j, v2) - Formula_y_2(j, v2));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y^{(2)} Failed\n";
        out << "Validation failed at index " << j << "," << v2 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }
    
  // Compare X_{(1)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i = 0; i < m; i++) {
      diff = std::abs(AD_x_1(u1, i) - Formula_x_1(u1, i));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for X_{(1)} Failed\n";
        out << "Validation failed at index " << u1 << "," << i << " Diff: " << diff << "\n";
        return false;
      }

      // Compare X_{(1)}^{(2)}
      for (int v2 = 0; v2 < V2; v2++) {
        diff = std::abs(AD_x_1_2(u1)(i, v2) - Formula_x_1_2(u1)(i, v2));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for X_{(1)}^{(2)} Failed\n";
          out << "Validation failed at index " << u1 << "," << i << "," << v2 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}