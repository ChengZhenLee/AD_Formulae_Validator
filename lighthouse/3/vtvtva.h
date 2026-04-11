// Third derivative by vector tangent of vector tangent of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void vtvtva_F_xxx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx, Y_XXX_t<T>& y_xxx) {
  // activate x and y
  X_t<A_t<T_t<T_t<T,n>,n>,m>> x; Y_t<A_t<T_t<T_t<T,n>,n>,m>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value().value().value()=x_values(i);
    // set X^{(2)}
    x(i).value().tangent(i).value()=1;
    // set X^{(3)}
    x(i).value().value().tangent(i)=1;
    // register x with A_t<T_t<T,n>,m>::tape
    x(i).register_input(); 
  }
  // run overloaded F
  F(x,y);
  // allocate vector of adjoints of A_t<T_t<T_t<T,n>,n>,m>::tape
  A_t<T_t<T_t<T,n>,n>,m>::tape::init_adjoints();
  // set Y_{(1)}
  for (int m1=0;m1<m;++m1) 
    y(m1).adjoint(m1).value().value()=1;
  // interpret A_t<T_t<T_t<T,n>,n>,m>::tape
  A_t<T_t<T_t<T,n>,n>,m>::tape::interpret();
  // extract third derivative from X_{(1)}^{(2,3)}
  for (int m1=0;m1<m;++m1) 
    for (int i1=0;i1<n;++i1) 
      for (int n2=0;n2<n;++n2) 
        for (int n3=0;n3<n;++n3) 
          y_xxx(m1)(i1)(n2)(n3)=x(i1).adjoint(m1).tangent(n2).tangent(n3);
  for (int m1 = 0; m1 < m; m1++) {
    for (int n2 = 0; n2 < n; n2++) {
      // Extract Jacobian
      y_x(m1)(n2) = y(m1).value().tangent(n2).value();
      for (int n3 = 0; n3 < n; n3++) {
        // Extract Hessian
        y_xx(m1)(n2)(n3) = y(m1).value().tangent(n2).tangent(n3);
      }
    }
  }
  A_t<T_t<T_t<T,n>,n>,m>::tape::reset();

}	


template<typename T, int U1, int V2, int V3>
void AD_F_xxx(
  const X_t<T>& x_values,
  const Eigen::Matrix<T, U1, m>& y_1, const Eigen::Matrix<T, n, V2>& x_2, const Eigen::Matrix<T, n, V3>& x_3,
  const Eigen::Vector<Eigen::Matrix<T, m, V2>, U1>& y_1_2, const Eigen::Vector<Eigen::Matrix<T, m, V3>, U1>& y_1_3,
  const Eigen::Vector<Eigen::Matrix<T, V2, V3>, n>& x_2_3, const Eigen::Matrix<Eigen::Matrix<T, V2, V3>, U1, m>& y_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, U1, n>& x_1, Eigen::Matrix<T, m, V2>& y_2, Eigen::Matrix<T, m, V3>& y_3,
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1>& x_1_2, Eigen::Vector<Eigen::Matrix<T, n, V3>, U1>& x_1_3, 
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m>& y_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, U1, n>& x_1_2_3
) {

  X_t<A_t<T_t<T_t<T, V3>, V2>, U1>> x; Y_t<A_t<T_t<T_t<T, V3>, V2>, U1>> y;

  for (int i = 0; i < n; i++) {
    // Set X^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      x(i).value().tangent(v2).value() = x_2(i, v2);
      // Set X^{(2, 3)}
      for (int v3 = 0; v3 < V3; v3++) {
        x(i).value().tangent(v2).tangent(v3) = x_2_3(i)(v2, v3);
      }
    }
    // Set X^{(3)}
    for (int v3 = 0; v3 < V3; v3++) {
      x(i).value().value().tangent(v3) = x_3(i, v3);
    }

    x(i).register_input();
  }

  F(x, y);

  A_t<T_t<T_t<T, V3>, V2>, U1>::tape::init_adjoints();

  // Set Y_{(1)}
  for (int j = 0; j < m; j++) {
    for (int u1 = 0; u1 < U1; u1++) {
      y(j).adjoint(u1).value().value() = y_1(u1, j);
      // Set Y_{(1)}^{(2)}
      for (int v2 = 0; v2 < V2; v2++) {
        y(j).adjoint(u1).tangent(v2).value() = y_1_2(u1)(j, v2);
        // Set Y_{(1)}^{(2, 3)}
        for (int v3 = 0; v3 < V3; v3++) {
          y(j).adjoint(u1).tangent(v2).tangent(v3) = y_1_2_3(u1, j)(v2, v3);
        }
      }
      // Set Y_{(1)}^{(3)}
      for (int v3 = 0; v3 < V3; v3++) {
        y(j).adjoint(u1).value().tangent(v3) = y_1_3(u1)(j, v3);
      }
    }
  }

  A_t<T_t<T_t<T, V3>, V2>, U1>::tape::interpret();

  // Extract X_{(1)}
  for (int i = 0; i < n; i++) {
    for (int u1 = 0; u1 < U1; u1++) {
      x_1(u1, i) = x(i).adjoint(u1).value().value();
      // Extract X_{(1)}^{(2)}
      for (int v2 = 0; v2 < V2; v2++) {
        x_1_2(u1)(i, v2) = x(i).adjoint(u1).tangent(v2).value();
        // Extract X_{(1)}^{(2, 3)}
        for (int v3 = 0; v3 < V3; v3++) {
          x_1_2_3(u1, i)(v2, v3) = x(i).adjoint(u1).tangent(v2).tangent(v3);
        }
      }
      // Extract X_{(1)}^{(3)}
      for (int v3 = 0; v3 < V3; v3++) {
        x_1_3(u1)(i, v3) = x(i).adjoint(u1).value().tangent(v3);
      }
    }
  }

  // Extract y
  for (int j = 0; j < m; j++) {
    y_values(j) = y(j).value().value().value();
    // Extract Y^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      y_2(j, v2) = y(j).value().tangent(v2).value();
      // Extract Y^{(2, 3)}
      for (int v3 = 0; v3 < V3; v3++) {
        y_2_3(j)(v2, v3) = y(j).value().tangent(v2).tangent(v3);
      }
    }
    // Extract Y^{(3)}
    for (int v3 = 0; v3 < V3; v3++) {
      y_3(j, v3) = y(j).value().value().tangent(v3);
    }
  }

  A_t<T_t<T_t<T, V3>, V2>, U1>::tape::reset();
}


template<typename T, int U1, int V2, int V3>
void Formula_F_xxx(
  const Y_X_t<T>& y_x, const Y_XX_t<T>& y_xx, const Y_XXX_t<T>&  y_xxx,
  const X_t<T>& x_values,
  const Eigen::Matrix<T, U1, m>& y_1, const Eigen::Matrix<T, n, V2>& x_2, const Eigen::Matrix<T, n, V3>& x_3,
  const Eigen::Vector<Eigen::Matrix<T, m, V2>, U1>& y_1_2, const Eigen::Vector<Eigen::Matrix<T, m, V3>, U1>& y_1_3,
  const Eigen::Vector<Eigen::Matrix<T, V2, V3>, n>& x_2_3, const Eigen::Matrix<Eigen::Matrix<T, V2, V3>, U1, m>& y_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, U1, n>& x_1, Eigen::Matrix<T, m, V2>& y_2, Eigen::Matrix<T, m, V3>& y_3,
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1>& x_1_2, Eigen::Vector<Eigen::Matrix<T, n, V3>, U1>& x_1_3, 
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m>& y_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, U1, n>& x_1_2_3
) {
  // y = f(x)
  F(x_values, y_values);

  for (int j = 0; j < m; j++) {
    // Y^{(3)} = F' * X^{(3)}
    for (int v3 = 0; v3 < V3; v3++) {
      T sum = 0;
      for (int i3 = 0; i3 < n; i3++) {
        sum += y_x(j)(i3) * x_3(i3, v3);
      }
      y_3(j, v3) = sum;
    }
    // Y^{(2)} = F' * X^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      T sum1 = 0;
      for (int i2 = 0; i2 < n; i2++) {
        sum1 += y_x(j)(i2) * x_2(i2, v2);
        // Y^{(2, 3)} = F'' * X^{(2)} * X^{(3)} + F' * X^{(2, 3)}
        for (int v3 = 0; v3 < V3; v3++) {
          T sum2 = 0;
          for (int i3 = 0; i3 < n; i3++) {
            sum2 += y_xx(j)(i2)(i3) * x_2(i2, v2) * x_3(i3, v3);
          }
          y_2_3(j)(v2, v3) = sum2;
        }
      }
      y_2(j, v2) = sum1;
    }
  }

  // X_{(1)} += F' * Y_{(1)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i = 0; i < n; i++) {
      T sum = 0;
      for (int j = 0; j < m; j++) {
        sum += y_x(j)(i) * y_1(u1, j);
      }
      x_1(u1, i) = sum;
    }
  }

  // X_{(1)}^{(3)} += F'' * Y_{(1)} * X^{(3)} + F' * Y_{(1)}^{(3)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i1 = 0; i1 < n; i1++) {
      for (int v3 = 0; v3 < V3; v3++) {
        T term1 = 0;
        T term2 = 0;
        for (int j = 0; j < m; j++) {
          term2 += y_x(j)(i1) * y_1_3(u1)(j, v3);
          for (int i3 = 0; i3 < n; i3++) {
            term1 += y_xx(j)(i1)(i3) * y_1(u1, j) * x_3(i3, v3);
          }
        }
        x_1_3(u1)(i1, v3) = term1 + term2;
      }
    }
  }

  // X_{(1)}^{(2)} += F'' * Y_{(1)} * X^{(2)} + F' * Y_{(1)}^{(2)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i1 = 0; i1 < n; i1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        T term1 = 0;
        T term2 = 0;
        for (int j = 0; j < m; j++) {
          term2 += y_x(j)(i1) * y_1_3(u1)(j, v2);
          for (int i2 = 0; i2 < n; i2++) {
            term1 += y_xx(j)(i1)(i2) * y_1(u1, j) * x_3(i2, v2);
          }
        }
        x_1_3(u1)(i1, v2) = term1 + term2;
      }
    }
  }

  // X_{(1)}^{(2, 3)} += F''' * Y_{(1)} * X^{(2)} * X^{(3)}
  // + F'' * Y_{(1)}^{(3)} * X^{(2)}
  // + F'' * Y_{(1)}^{(2)} * X^{(3)}
  // + F'' * Y_{(1)} * X^{(2, 3)}
  // + F' * Y_{(1)}^{(2, 3)}
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i1 = 0; i1 < n; i1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        for (int v3 = 0; v3 < V3; v3++) {
          T term1 = 0;
          T term2 = 0;
          T term3 = 0;
          T term4 = 0;
          T term5 = 0;
          for (int j = 0; j < m; j++) {
            term5 += y_x(j)(i1) * y_1_2_3(u1, j)(v2, v3);
            for (int i2 = 0; i2 < n; i2++) {
              for (int i3 = 0; i3 < n; i3++) {
                term1 += y_xxx(j)(i1)(i2)(i3) * y_1(u1, j) * x_2(i2, v2) * x_3(i3, v3);
              }
              term2 += y_xx(j)(i1)(i2) * y_1_3(u1)(j, v3) * x_2(i2, v2);
              term3 += y_xx(j)(i1)(i2) * y_1(u1, j) * x_2_3(i2)(v2, v3);
            }
            for (int i3 = 0; i3 < n; i3++) {
              term4 += y_xx(j)(i1)(i3) * y_1_2(u1)(j, v2) * x_3(i3, v3);
            }
          }
        }
      }
    }
  }
}


template<typename T, int U1, int V2, int V3>
bool Validate_vtvtva(std::ofstream& out) {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 8.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;
  Y_XXX_t<T> y_xxx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, U1, n> AD_x_1;
  Eigen::Matrix<T, U1, n> Formula_x_1;
  Eigen::Matrix<T, m, V2> AD_y_2;
  Eigen::Matrix<T, m, V2> Formula_y_2;
  Eigen::Matrix<T, m, V3> AD_y_3;
  Eigen::Matrix<T, m, V3> Formula_y_3;
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1> AD_x_1_2;
  Eigen::Vector<Eigen::Matrix<T, n, V2>, U1> Formula_x_1_2;
  Eigen::Vector<Eigen::Matrix<T, n, V3>, U1> AD_x_1_3;
  Eigen::Vector<Eigen::Matrix<T, n, V3>, U1> Formula_x_1_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m> AD_y_2_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m> Formula_y_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, U1, n> AD_x_1_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, U1, n> Formula_x_1_2_3;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, U1, m> y_1 = Eigen::Matrix<T, U1, m>::Random().cwiseAbs();
  Eigen::Matrix<T, n, V2> x_2 = Eigen::Matrix<T, n ,V2>::Random().cwiseAbs();
  Eigen::Matrix<T, n, V3> x_3 = Eigen::Matrix<T, n ,V3>::Random().cwiseAbs();
  Eigen::Vector<Eigen::Matrix<T, m, V2>, U1> y_1_2;
  Eigen::Vector<Eigen::Matrix<T, m, V3>, U1> y_1_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, n> x_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, U1, m> y_1_2_3;
  for (int u1 = 0; u1 < U1; u1++) {
    y_1_2(u1) = Eigen::Matrix<T, m, V2>::Random().cwiseAbs();
    y_1_3(u1) = Eigen::Matrix<T, m, V3>::Random().cwiseAbs();
    for (int j = 0; j < m; j++) {
      y_1_2_3(u1, j) = Eigen::Matrix<T, V2, V3>::Random().cwiseAbs();
    }
  }
  for (int i = 0; i < n; i++) {
    x_2_3(i) = Eigen::Matrix<T, V2, V3>::Random().cwiseAbs();
  }

  out << "\n=== Testing Third Derivative (Tangent over Tangent over Adjoint Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for Y_{(1)}: \n" << y_1 << "\n\n";
  out << "Seed for X^{(2)}: \n" << x_2 << "\n\n";
  out << "Seed for X^{(3)}: \n" << x_3 << "\n\n";
  out << "Seed for Y_{(1)}^{(2)}: (nested matrices not printed)\n\n";
  out << "Seed for Y_{(1)}^{(3)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(2, 3)}: (nested matrices not printed)\n\n";
  out << "Seed for Y_{(1)}^{(2, 3)}: (nested matrices not printed)\n\n";

  // Populate y_x, y_xx and y_xxx
  vtvtva_F_xxx(x_values, y_x, y_xx, y_xxx);

  // Run the AD version
  AD_F_xxx(
    x_values, y_1, x_2, x_3, 
    y_1_2, y_1_3, x_2_3,
    y_1_2_3,
    AD_y_values, 
    AD_x_1, AD_y_2, AD_y_3,
    AD_x_1_2, AD_x_1_3, AD_y_2_3,
    AD_x_1_2_3
  );

  // Run the Formula version
  Formula_F_xxx(y_x, y_xx, y_xxx, 
    x_values, y_1, x_2, x_3, 
    y_1_2, y_1_3, x_2_3,
    y_1_2_3,
    Formula_y_values, 
    Formula_x_1, Formula_y_2, Formula_y_3,
    Formula_x_1_2, Formula_x_1_3, Formula_y_2_3,
    Formula_x_1_2_3
  );

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

    // Compare y_2
    for (int v2 = 0; v2 < V2; v2++) {
      diff = std::abs(AD_y_2(j, v2) - Formula_y_2(j, v2));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y^{(2)} Failed\n";
        out << "Validation failed at index " << j << "," << v2 << " Diff: " << diff << "\n";
        return false;
      }

      // Compare y_2_3
      for (int v3 = 0; v3 < V3; v3++) {
        diff = std::abs(AD_y_2_3(j)(v2, v3) - Formula_y_2_3(j)(v2, v3));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y^{(2, 3)} Failed\n";
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
        out << "Validation for Y^{(3)} Failed\n";
        out << "Validation failed at index " << j << "," << v3 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }

  // Compare x_1
  for (int u1 = 0; u1 < U1; u1++) {
    for (int i1 = 0; i1 < n; i1++) {
      diff = std::abs(AD_x_1(u1, i1) - Formula_x_1(u1, i1));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for X_{(1)} Failed\n";
        out << "Validation failed at index " << u1 << "," << i1 << " Diff: " << diff << "\n";
        return false;
      }

      // Compare x_1_2
      for (int v2 = 0; v2 < V2; v2++) {
        diff = std::abs(AD_x_1_2(u1)(i1, v2) - Formula_x_1_2(u1)(i1, v2));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for X_{(1)}^{(2)} Failed\n";
          out << "Validation failed at index " << u1 << "," << i1 << "," << v2 << " Diff: " << diff << "\n";
          return false;
        }

        // Compare x_1_2_3
        for (int v3 = 0; v3 < V3; v3++) {
          diff = std::abs(AD_x_1_2_3(u1, i1)(v2, v3) - Formula_x_1_2_3(u1, i1)(v2, v3));
          maxDiff = std::max(maxDiff, diff);
          if (diff > tol) {
            out << "Validation for X_{(1)}^{(2, 3)} Failed\n";
            out << "Validation failed at index " << u1 << "," << i1 << "," << v2 << "," << v3 << " Diff: " << diff << "\n";
            return false;
          }
        }
      }

      // Compare x_1_3
      for (int v3 = 0; v3 < V3; v3++) {
        diff = std::abs(AD_x_1_2(u1)(i1, v3) - Formula_x_1_2(u1)(i1, v3));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for X_{(1)}^{(3)} Failed\n";
          out << "Validation failed at index " << u1 << "," << i1 << "," << v3 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}