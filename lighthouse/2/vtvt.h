// Hessian by vector tangent of vector tangent AD

#include "F.h"
#include "ad_types.h"

template<typename T>
void vtvt_F_xx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx) {
  // activate x and y
  X_t<T_t<T_t<T,n>,n>> x; Y_t<T_t<T_t<T,n>,n>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value().value()=x_values(i);
    // set X^{(1)}
    x(i).tangent(i).value()=1;
    // set X^{(2)}
    x(i).value().tangent(i)=1;
  }
  // run overloaded F
  F(x,y);
  // extract Hessian from Y^{(1,2)}
  for (int j=0;j<m;++j) 
    for (int n1=0;n1<n;++n1) 
      for (int n2=0;n2<n;++n2) 
        y_xx(j)(n1)(n2)=y(j).tangent(n1).tangent(n2);

  // extract Jacobian from Y^{(1)}
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      y_x(j)(i) = y(j).tangent(i).value();
    }
  }
}


template<typename T, int V1, int V2>
void AD_F_xx(const X_t<T>& x_values, 
  const Eigen::Matrix<T, n, V1>& x_1, const Eigen::Matrix<T, n, V2>& x_2, 
  const Eigen::Vector<Eigen::Matrix<T, V1, V2>, n>& x_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, m, V2>& y_2,
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m>& y_1_2) {

  X_t<T_t<T_t<T,V2>,V1>> x; Y_t<T_t<T_t<T,V2>,V1>> y;
  for (int i=0;i<n;++i) {
    // Seed x
    x(i).value().value()=x_values(i);

    // Seed X^{(1)}
    for (int v1 = 0; v1 < V1; v1++) {
      x(i).tangent(v1).value() = x_1(i, v1);
    }

    // Seed X^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      x(i).value().tangent(v2) = x_2(i, v2);
    }

    // Seed X^{(1, 2)}
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        x(i).tangent(v1).tangent(v2) = x_1_2(i)(v1, v2);
      }
    }
  }

  F(x,y);

  // Extract y
  for (int j = 0; j < m; j++) {
    y_values(j) = y(j).value().value();
  }

  // Extract Y^({1})
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      y_1(j, v1) = y(j).tangent(v1).value();
    }
  }

  // Extract Y^({2})
  for (int j = 0; j < m; j++) {
    for (int v2 = 0; v2 < V2; v2++) {
      y_2(j, v2) = y(j).value().tangent(v2);
    }
  }

  // Extract Y^({1, 2})
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        y_1_2(j)(v1, v2) = y(j).tangent(v1).tangent(v2);
      }
    }
  }
}


template<typename T, int V1, int V2>
void Formula_F_xx(const Y_X_t<T>& y_x, const Y_XX_t<T>& y_xx,
  const X_t<T>& x_values, 
  const Eigen::Matrix<T, n, V1>& x_1, const Eigen::Matrix<T, n, V2>& x_2, 
  const Eigen::Vector<Eigen::Matrix<T, V1, V2>, n>& x_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, m, V2>& y_2,
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m>& y_1_2) {
  
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

  // Y^{(1, 2)} = F'' * X^{(1)} * X^{(2)} + F' * X^{(1, 2)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        T term1 = 0;
        T term2 = 0;
        for (int i1 = 0; i1 < n; i1++) {
          term2 += y_x(j)(i1) * x_1_2(i1)(v1, v2);
          for (int i2 = 0; i2 < n; i2++) {
            term1 += y_xx(j)(i1)(i2) * x_1(i1, v1) * x_2(i2, v2);
          }
        }
        y_1_2(j)(v1, v2) = term1 + term2;
      }
    }
  }
}


template<typename T, int V1, int V2>
bool Validate_vtvt() {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 4.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, m, V1> AD_y_1;
  Eigen::Matrix<T, m, V1> Formula_y_1;
  Eigen::Matrix<T, m, V2> AD_y_2;
  Eigen::Matrix<T, m, V2> Formula_y_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m> AD_y_1_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m> Formula_y_1_2;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, n, V1> x_1 = Eigen::Matrix<T, n, V1>::Random().cwiseAbs();
  Eigen::Matrix<T, n, V2> x_2 = Eigen::Matrix<T, n, V2>::Random().cwiseAbs();
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, n> x_1_2;
  for (int i = 0; i < n; i++) {
    x_1_2(i) = Eigen::Matrix<T, V1, V2>::Random().cwiseAbs();
  }

  // Show the seeds
  std::cout << "Seed for x: \n" << x_values << "\n\n";
  std::cout << "Seed for X^({1}): \n" << x_1 << "\n\n";
  std::cout << "Seed for X^({2}): \n" << x_2 << "\n\n";
  std::cout << "Seed for X^{(1, 2)}: (nested matrices not printed)\n\n";

  // Populate y_x and y_xx
  vtvt_F_xx(x_values, y_x, y_xx);

  // Run the AD version
  AD_F_xx(x_values, x_1, x_2, x_1_2, 
          AD_y_values, AD_y_1, AD_y_2, AD_y_1_2);

  // Run the Formula version
  Formula_F_xx(y_x, y_xx, 
    x_values, x_1, x_2, x_1_2, 
              Formula_y_values, Formula_y_1, Formula_y_2, Formula_y_1_2);

  T diff;
  T maxDiff = 0;

  // Compare y
  for (int j = 0; j < m; j++) {
    diff = std::abs(AD_y_values(j) - Formula_y_values(j));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      std::cout << "Validation for y Failed\n";
      std::cout << "Validation failed at index " << j << " Diff: " << diff << "\n";
      return false;
    }

    // Compare y_1
    for (int v1 = 0; v1 < V1; v1++) {
      diff = std::abs(AD_y_1(j, v1) - Formula_y_1(j, v1));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        std::cout << "Validation for Y_{(1)} Failed\n";
        std::cout << "Validation failed at index " << j << "," << v1 << " Diff: " << diff << "\n";
        return false;
      }

      // Compare y_1_2
      for (int v2 = 0; v2 < V2; v2++) {
        diff = std::abs(AD_y_1_2(j)(v1, v2) - Formula_y_1_2(j)(v1, v2));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          std::cout << "Validation for Y_{(1, 2)} Failed\n";
          std::cout << "Validation failed at index " << j << "," << v1 << "," << v2 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }
    // Compare y_2
    for (int v2 = 0; v2 < V2; v2++) {
      diff = std::abs(AD_y_2(j, v2) - Formula_y_2(j, v2));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        std::cout << "Validation for Y_{(2)} Failed\n";
        std::cout << "Validation failed at index " << j << "," << v2 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }
  std::cout << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}