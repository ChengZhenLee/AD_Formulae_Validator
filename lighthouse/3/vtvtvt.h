// Hessian by vector tangent of vector tangent of vector tangent AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void vtvtvt_F_xxx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx, Y_XXX_t<T>& y_xxx) {
  // activate x and y
  X_t<T_t<T_t<T_t<T,n>,n>,n>> x; Y_t<T_t<T_t<T_t<T,n>,n>,n>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value().value()=x_values(i);
    // set X^{(1)}
    x(i).tangent(i).value().value()=1;
    // set X^{(2)}
    x(i).value().tangent(i).value()=1;
    // set X^{(3)}
    x(i).value().value().tangent(i)=1;
  }
  // run overloaded F
  F(x,y);
  for (int j=0;j<m;++j) 
    for (int n1=0;n1<n;++n1) {
      // Extract Jacobian
      y_x(j)(n1) = y(j).tangent(n1).value().value();
      for (int n2=0;n2<n;++n2) {
        // Extract Hessian
        y_xx(j)(n1)(n2) = y(j).tangent(n1).tangent(n2).value();
        for (int n3=0;n3<n;++n3) 
          // extract F^{[3]} from Y^{(1,2,3)}
          y_xxx(j)(n1)(n2)(n3)=y(j).tangent(n1).tangent(n2).tangent(n3);
    }
  }
}	


template<typename T, int V1, int V2, int V3> 
void AD_F_xxx(
  const X_t<T>& x_values,
  const Eigen::Matrix<T, n, V1>& x_1, const Eigen::Matrix<T, n, V2>& x_2, const Eigen::Matrix<T, n, V3>& x_3,
  const Eigen::Vector<Eigen::Matrix<T, V1, V2>, n>& x_1_2, const Eigen::Vector<Eigen::Matrix<T, V1, V3>, n>& x_1_3,
  const Eigen::Vector<Eigen::Matrix<T, V2, V3>, n>& x_2_3, const Eigen::Matrix<Eigen::Matrix<T, V2, V3>, n, V1>& x_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, m, V2>& y_2, Eigen::Matrix<T, m, V3>& y_3,
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m>& y_1_2, Eigen::Vector<Eigen::Matrix<T, V1, V3>, m>& y_1_3, 
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m>& y_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1>& y_1_2_3
) { 

  X_t<T_t<T_t<T_t<T,V3>,V2>,V1>> x; Y_t<T_t<T_t<T_t<T,V3>,V2>,V1>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value().value() = x_values(i);
    // set X^{(1)}
    for (int v1 = 0; v1 < V1; v1++) {
      x(i).tangent(v1).value().value() = x_1(i, v1);
    }
    // set X^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      x(i).value().tangent(v2).value() = x_2(i, v2);
    }
    // set X^{(3)}
    for (int v3 = 0; v3 < V3; v3++) {
      x(i).value().value().tangent(v3) = x_3(i, v3);
    }
    // set X^{(1, 2)}
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        x(i).tangent(v1).tangent(v2).value() = x_1_2(i)(v1, v2);
        // set X^{(1, 2, 3)}
        for (int v3 = 0; v3 < V3; v3++) {
          x(i).tangent(v1).tangent(v2).tangent(v3) = x_1_2_3(i, v1)(v2, v3);
        }
      }
    }
    // set X^{(1, 3)}
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v3 = 0; v3 < V3; v3++) {
        x(i).tangent(v1).value().tangent(v3) = x_1_3(i)(v1, v3);
      }
    }
    // set X^{(2, 3)}
    for (int v2 = 0; v2 < V2; v2++) {
      for (int v3 = 0; v3 < V3; v3++) {
        x(i).value().tangent(v2).tangent(v3) = x_2_3(i)(v2, v3);
      }
    }
  }

  F(x,y);

  for (int j = 0; j < m; j++) {
    // Extract y
    y_values(j) = y(j).value().value().value();
    // Extract Y^{(1)}
    for (int v1 = 0; v1 < V1; v1++) {
      y_1(j, v1) = y(j).tangent(v1).value().value();
    }
    // Extract Y^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      y_2(j, v2) = y(j).value().tangent(v2).value();
    }
    // Extract Y^{(3)}
    for (int v3 = 0; v3 < V3; v3++) {
      y_3(j, v3) = y(j).value().value().tangent(v3);
    }
    // Extract Y^{(1, 2)}
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        y_1_2(j)(v1, v2) = y(j).tangent(v1).tangent(v2).value();
        // Extract Y^{(1, 2, 3)}
        for (int v3 = 0; v3 < V3; v3++) {
          y_1_2_3(j, v1)(v2, v3) = y(j).tangent(v1).tangent(v2).tangent(v3);
        }
      }
    }
    // Extract Y^{(1, 3)}
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v3 = 0; v3 < V3; v3++) { 
        y_1_3(j)(v1, v3) = y(j).tangent(v1).value().tangent(v3);
      }
    }
    // Extract Y^{(2, 3)}
    for (int v2 = 0; v2 < V2; v2++) {
      for (int v3 = 0; v3 < V3; v3++) {
        y_2_3(j)(v2, v3) = y(j).value().tangent(v2).tangent(v3);
      }
    }
  }
}


template<typename T, int V1, int V2, int V3>
void Formula_F_xxx(const Y_X_t<T>& y_x, Y_XX_t<T>& y_xx, Y_XXX_t<T>& y_xxx,
  const X_t<T>& x_values,
  const Eigen::Matrix<T, n, V1>& x_1, const Eigen::Matrix<T, n, V2>& x_2, const Eigen::Matrix<T, n, V3>& x_3,
  const Eigen::Vector<Eigen::Matrix<T, V1, V2>, n>& x_1_2, const Eigen::Vector<Eigen::Matrix<T, V1, V3>, n>& x_1_3,
  const Eigen::Vector<Eigen::Matrix<T, V2, V3>, n>& x_2_3, const Eigen::Matrix<Eigen::Matrix<T, V2, V3>, n, V1>& x_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, m, V2>& y_2, Eigen::Matrix<T, m, V3>& y_3,
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m>& y_1_2, Eigen::Vector<Eigen::Matrix<T, V1, V3>, m>& y_1_3, 
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m>& y_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1>& y_1_2_3
) {
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
  // Y^{(3)} = F' * X^{(3)}
  for (int j = 0; j < m; j++) {
    for (int v3 = 0; v3 < V3; v3++) {
      T sum = 0;
      for (int i = 0; i < n; i++) {
        sum += y_x(j)(i) * x_3(i, v3);
      }
      y_3(j, v3) = sum;
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
  // Y^{(1, 3)} = F'' * X^{(1)} * X^{(3)} + F' * X^{(1, 3)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v3 = 0; v3 < V3; v3++) {
        T term1 = 0;
        T term2 = 0;
        for (int i1 = 0; i1 < n; i1++) {
          term2 += y_x(j)(i1) * x_1_3(i1)(v1, v3);
          for (int i3 = 0; i3 < n; i3++) {
            term1 += y_xx(j)(i1)(i3) * x_1(i1, v1) * x_3(i3, v3);
          }
        }
        y_1_3(j)(v1, v3) = term1 + term2;
      }
    }
  }  
  // Y^{(2, 3)} = F'' * X^{(2)} * X^{(3)} + F' * X^{(2, 3)}
  for (int j = 0; j < m; j++) {
    for (int v2 = 0; v2 < V2; v2++) {
      for (int v3 = 0; v3 < V3; v3++) {
        T term1 = 0;
        T term2 = 0;
        for (int i2 = 0; i2 < n; i2++) {
          term2 += y_x(j)(i2) * x_2_3(i2)(v2, v3);
          for (int i3 = 0; i3 < n; i3++) {
            term1 += y_xx(j)(i2)(i3) * x_2(i2, v2) * x_3(i3, v3);
          }
        }
        y_2_3(j)(v2, v3) = term1 + term2;
      }
    }
  }

  // Y^{(1, 2, 3)} = F''' * X^{(1)} * X^{(2)} * X^{(3)}
  // + F'' * X^{(1)} * X^{(2, 3)}
  // + F'' * X^{(2)} * X^{(1, 3)} 
  // + F'' * X^{(3)} * X^{(1, 2)}
  // + F' * X^{(1, 2, 3)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        for (int v3 = 0; v3 < V3; v3++) {
          T term1 = 0;
          T term2 = 0;
          T term3 = 0;
          T term4 = 0;
          T term5 = 0;
          for (int i1 = 0; i1 < n; i1++) {
            term5 += y_x(j)(i1) * x_1_2_3(i1, v1)(v2, v3);
            for (int i2 = 0; i2 < n; i2++) {
              term2 += y_xx(j)(i1)(i2) * x_1(i1, v1) * x_2_3(i2)(v2, v3);
              term3 += y_xx(j)(i1)(i2) * x_2(i2, v2) * x_1_3(i1)(v1, v3);
              for (int i3 = 0; i3 < n; i3++) {
                term1 += y_xxx(j)(i1)(i2)(i3) * x_1(i1, v1) * x_2(i2, v2) * x_3(i3, v3);
              }
            }
            for (int i3 = 0; i3 < n; i3++) {
              term4 += y_xx(j)(i1)(i3) * x_3(i3, v3) * x_1_2(i1)(v1, v2);
            }
          }
          y_1_2_3(j, v1)(v2, v3) = term1 + term2 + term3 + term4 + term5;
        }
      }
    }
  }
}


template<typename T, int V1, int V2, int V3>
bool Validate_vtvtvt(std::ofstream& out) {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 8.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;
  Y_XXX_t<T> y_xxx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, m, V1> AD_y_1;
  Eigen::Matrix<T, m, V1> Formula_y_1;
  Eigen::Matrix<T, m, V2> AD_y_2;
  Eigen::Matrix<T, m, V2> Formula_y_2;
  Eigen::Matrix<T, m, V3> AD_y_3;
  Eigen::Matrix<T, m, V3> Formula_y_3;
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m> AD_y_1_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m> Formula_y_1_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V3>, m> AD_y_1_3;
  Eigen::Vector<Eigen::Matrix<T, V1, V3>, m> Formula_y_1_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m> AD_y_2_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m> Formula_y_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1> AD_y_1_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1> Formula_y_1_2_3;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, n, V1> x_1 = Eigen::Matrix<T, n ,V1>::Random().cwiseAbs();
  Eigen::Matrix<T, n, V2> x_2 = Eigen::Matrix<T, n ,V2>::Random().cwiseAbs();
  Eigen::Matrix<T, n, V3> x_3 = Eigen::Matrix<T, n ,V3>::Random().cwiseAbs();
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, n> x_1_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V3>, n> x_1_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, n> x_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, n, V1> x_1_2_3;
  for (int i = 0; i < n; i++) {
    x_1_2(i) = Eigen::Matrix<T, V1, V2>::Random().cwiseAbs();
    x_1_3(i) = Eigen::Matrix<T, V1, V3>::Random().cwiseAbs();
    x_2_3(i) = Eigen::Matrix<T, V2, V3>::Random().cwiseAbs();
    for (int v1 = 0; v1 < V1; v1++) {
      x_1_2_3(i, v1) = Eigen::Matrix<T, V2, V3>::Random().cwiseAbs();
    }
  }

  out << "\n=== Testing Third Derivative (Tangent over Tangent over Tangent Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for X^({1}): \n" << x_1 << "\n\n";
  out << "Seed for X^({2}): \n" << x_2 << "\n\n";
  out << "Seed for X^{(3)}: \n" << x_3 << "\n\n";
  out << "Seed for X^{(1, 2)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 3)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(2, 3)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 2, 3)}: (nested matrices not printed)\n\n";

  // Populate y_x, y_xx and y_xxx
  vtvtvt_F_xxx(x_values, y_x, y_xx, y_xxx);

  // Run the AD version
  AD_F_xxx(
    x_values, x_1, x_2, x_3, 
    x_1_2, x_1_3, x_2_3,
    x_1_2_3,
    AD_y_values, 
    AD_y_1, AD_y_2, AD_y_3,
    AD_y_1_2, AD_y_1_3, AD_y_2_3,
    AD_y_1_2_3
  );

  // Run the Formula version
  Formula_F_xxx(y_x, y_xx, y_xxx, 
    x_values, x_1, x_2, x_3, 
    x_1_2, x_1_3, x_2_3,
    x_1_2_3,
    Formula_y_values, 
    Formula_y_1, Formula_y_2, Formula_y_3,
    Formula_y_1_2, Formula_y_1_3, Formula_y_2_3,
    Formula_y_1_2_3);

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

    // Compare y_1
    for (int v1 = 0; v1 < V1; v1++) {
      diff = std::abs(AD_y_1(j, v1) - Formula_y_1(j, v1));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y_{(1)} Failed\n";
        out << "Validation failed at index " << j << "," << v1 << " Diff: " << diff << "\n";
        return false;
      }

      // Compare y_1_2
      for (int v2 = 0; v2 < V2; v2++) {
        diff = std::abs(AD_y_1_2(j)(v1, v2) - Formula_y_1_2(j)(v1, v2));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y_{(1, 2)} Failed\n";
          out << "Validation failed at index " << j << "," << v1 << "," << v2 << " Diff: " << diff << "\n";
          return false;
        }

        // Compare y_1_2_3
        for (int v3 = 0; v3 < V3; v3++) {
          diff = std::abs(AD_y_1_2_3(j, v1)(v2, v3) - Formula_y_1_2_3(j, v1)(v2, v3));
          maxDiff = std::max(maxDiff, diff);
          if (diff > tol) {
            out << "Validation for Y_{(1, 2, 3)} Failed\n";
            out << "Validation failed at index " << j << "," << v1 << "," << v2 << "," << v3 <<" Diff: " << diff << "\n";
            return false;
          }
        }
      }

      // Compare y_1_3
      for (int v3 = 0; v3 < V3; v3++) {
        diff = std::abs(AD_y_1_3(j)(v1, v3) - Formula_y_1_3(j)(v1, v3));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y_{(1, 3)} Failed\n";
          out << "Validation failed at index " << j << "," << v1 << "," << v3 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }

    // Compare y_2
    for (int v2 = 0; v2 < V2; v2++) {
      diff = std::abs(AD_y_2(j, v2) - Formula_y_2(j, v2));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y_{(2)} Failed\n";
        out << "Validation failed at index " << j << "," << v2 << " Diff: " << diff << "\n";
        return false;
      }
      
      // Compare y_2_3
      for (int v3 = 0; v3 < V3; v3++) {
        diff = std::abs(AD_y_2_3(j)(v2, v3) - Formula_y_2_3(j)(v2, v3));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y_{(2, 3)} Failed\n";
          out << "Validation failed at index " << j << "," << v2 << "," << v3 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }

    // Compare y_3
    for (int v3 = 0; v3 < V3; v3++) {
      diff = std::abs(AD_y_3(j, v3) - Formula_y_3(j, v3));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y_{(3)} Failed\n";
        out << "Validation failed at index " << j << "," << v3 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}