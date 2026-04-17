// Third derivative by vector adjoint of vector tangent of vector adjoint AD

#include "F.h"
#include "ad_types.h"
#include <iostream>


template<typename T>
void vastsa_F_xxx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx, Y_XXX_t<T>& y_xxx) {
  for (int m1=0;m1<m;++m1) {
    for (int n2=0;n2<n;++n2) {
      // activate x and y
      X_t<A_t<T_t<A_t<T,n>,1>,1>> x; Y_t<A_t<T_t<A_t<T,n>,1>,1>> y;
      for (int i=0;i<n;++i) {
        // set x
        x(i).value().value().value()=x_values(i);
        // register x with A_t<T_t<A_t<T,n>,1>,1>::tape
        x(i).register_input(); 
        // register x with A_t<T,n>::tape
        x(i).value().value().register_input(); 
      }
      // set X^{(2)}
      x(n2).value().tangent().value()=1;
      // run overloaded F
      F(x,y);
      // allocate vector of adjoints of A_t<T_t<A_t<T,n>,n>,m>::tape
      A_t<T_t<A_t<T,n>,1>,1>::tape::init_adjoints();
      // set Y_{(1)}
      y(m1).adjoint().value().value()=1;
      // interpret A_t<T_t<A_t<T,n>,n>,m>::tape
      A_t<T_t<A_t<T,n>,1>,1>::tape::interpret();
      // set X_{(1,3)}^{(2)}
      A_t<T,n>::tape::init_adjoints();
      for (int m3=0;m3<n;++m3) 
        x(m3).adjoint().tangent().adjoint(m3)=1;
      A_t<T,n>::tape::interpret();
      // set X_{(1,3)}^{(2)}
      // extract third derivative from X_{(3)}
      for (int i1=0;i1<n;++i1) {
        for (int m3=0;m3<n;++m3) {
          y_xxx(m1)(i1)(n2)(m3)=x(i1).value().value().adjoint(m3);
        }
      }
      // extract Jacobian
      for (int i1 = 0; i1 < n; i1++) {
        y_x(m1)(i1) = x(i1).adjoint().value().value();
      }
      // extract Hessian
      for (int i1 = 0; i1 < n; i1++) {
        y_xx(m1)(i1)(n2) = x(i1).adjoint().tangent().value();
      }
      A_t<T,n>::tape::reset();
      A_t<T_t<A_t<T,n>,1>,1>::tape::reset();
    }	
  }
}


template<typename T, int U>
void TA_F_xxx(
  X_t<A_t<T_t<A_t<T, U>, 1>, 1>>& x, Y_t<A_t<T_t<A_t<T, U>, 1>, 1>>& y,
  const X_t<T>& x_values,
  const Eigen::Matrix<T, 1, m>& y_1, const Eigen::Matrix<T, n, 1>& x_2,
  const Eigen::Vector<Eigen::Matrix<T, m, 1>, 1>& y_1_2,
  Y_t<T>& y_values,
  Eigen::Matrix<T, 1, n>& x_1, Eigen::Matrix<T, m, 1>& y_2,
  Eigen::Vector<Eigen::Matrix<T, n, 1>, 1>& x_1_2
) {
  for (int i = 0; i < n; i++) {
    // Seed X^{(2)}
    x(i).value().tangent().value() = x_2(i, 0);

    x(i).register_input();
  }

  F(x, y);

  A_t<T_t<A_t<T, U>, 1>, 1>::tape::init_adjoints();

  // Seed Y_{(1)}
  for (int j = 0; j < m; j++) {
    y(j).adjoint().value().value() = y_1(0, j);

    // Seed Y_{(1)}^{(2)}
    y(j).adjoint().tangent().value() = y_1_2(0)(j, 0);
  }

  A_t<T_t<A_t<T, U>, 1>, 1>::tape::interpret();

  // Extract y
  for (int j = 0; j < m; j++) {
    y_values(j) = y(j).value().value().value();
  }

  // Extract Y^{(2)}
  for (int j = 0; j < m; j++) {
    y_2(j, 0) = y(j).value().tangent().value();
  }

  // Extract X_{(1)}
  for (int i = 0; i < n; i++) {
    x_1(0, i) = x(i).adjoint().value().value();

    // Extract X_{(1)}^{(2)}
    x_1_2(0)(i, 0) = x(i).adjoint().tangent().value();
  }

  A_t<T_t<A_t<T, U>, 1>, 1>::tape::reset();
}


template<typename T, int U>
void AD_F_xxx(
  const X_t<T>& x_values,
  const Eigen::Matrix<T, 1, m>& y_1, const Eigen::Matrix<T, n, 1>& x_2, const Eigen::Matrix<T, U, m>& y_3,
  const Eigen::Vector<Eigen::Matrix<T, m, 1>, 1>& y_1_2, const Eigen::Vector<Eigen::Matrix<T, m, 1>, U>& y_2_3,
  const Eigen::Vector<Eigen::Matrix<T, 1, n>, U>& x_1_3, const Eigen::Matrix<Eigen::Matrix<T, n, 1>, U, 1>& x_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, 1, n>& x_1, Eigen::Matrix<T, m, 1>& y_2, Eigen::Matrix<T, U, n>& x_3,
  Eigen::Vector<Eigen::Matrix<T, n, 1>, 1>& x_1_2, Eigen::Vector<Eigen::Matrix<T, 1, m>, U>& y_1_3, 
  Eigen::Vector<Eigen::Matrix<T, n, 1>, U>& x_2_3,
  Eigen::Matrix<Eigen::Matrix<T, m, 1>, U, 1>& y_1_2_3
) {

  X_t<A_t<T_t<A_t<T, U>, 1>, 1>> x; Y_t<A_t<T_t<A_t<T, U>, 1>, 1>> y;

  for (int i = 0; i < n; i++) {
    // Set x
    x(i).value().value().value() = x_values(i);
    x(i).value().value().register_input();
  }

  TA_F_xxx(x, y, x_values, y_1, x_2, y_1_2, y_values, x_1, y_2, x_1_2);

  A_t<T, U>::tape::init_adjoints();

  // Seed Y_{(3)}
  for (int j = 0; j < m; j++) {
    for (int u = 0; u < U; u++) {
      y(j).value().value().adjoint(u) = y_3(u, j);

      // Seed Y^{(2)}_{(3)}
      y(j).value().tangent().adjoint(u) = y_2_3(u)(j, 0);
    }
  }

  // Seed X_{(1, 3)}
  for (int i = 0; i < n; i++) {
    for (int u = 0; u < U; u++) {
      x(i).adjoint().value().adjoint(u) = x_1_3(u)(0, i);

      // Seed X_{(1, 3)}^{(2)}
      x(i).adjoint().tangent().adjoint(u) = x_1_2_3(u, 0)(i, 0);
    }
  }

  A_t<T, U>::tape::interpret();

  // Extract X_{(3)}
  for (int i = 0; i < n; i++) { 
    for (int u = 0; u < U; u++) {
      x_3(u, i) = x(i).value().value().adjoint(u);

      // Extract X^{(2)}_{(3)}
      x_2_3(u)(i, 0) = x(i).value().tangent().adjoint(u);
    }
  }

  // Extract Y_{(1, 3)}
  for (int j = 0; j < m; j++) {
    for (int u = 0; u < U; u++) {
      y_1_3(u)(0, j) = y(j).adjoint().value().adjoint(u);

      // Extract Y_{(1, 3)}^{(2)}
      y_1_2_3(u, 0)(j, 0) = y(j).adjoint().tangent().adjoint(u);
    }
  }
}


template<typename T, int U>
void Formula_F_xxx(
  Y_X_t<T>& y_x, Y_XX_t<T>& y_xx, Y_XXX_t<T>& y_xxx,
  const X_t<T>& x_values,
  const Eigen::Matrix<T, 1, m>& y_1, const Eigen::Matrix<T, n, 1>& x_2, const Eigen::Matrix<T, U, m>& y_3,
  const Eigen::Vector<Eigen::Matrix<T, m, 1>, 1>& y_1_2, const Eigen::Vector<Eigen::Matrix<T, m, 1>, U>& y_2_3,
  const Eigen::Vector<Eigen::Matrix<T, 1, n>, U>& x_1_3, const Eigen::Matrix<Eigen::Matrix<T, n, 1>, U, 1>& x_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, 1, n>& x_1, Eigen::Matrix<T, m, 1>& y_2, Eigen::Matrix<T, U, n>& x_3,
  Eigen::Vector<Eigen::Matrix<T, n, 1>, 1>& x_1_2, Eigen::Vector<Eigen::Matrix<T, 1, m>, U>& y_1_3, 
  Eigen::Vector<Eigen::Matrix<T, n, 1>, U>& x_2_3,
  Eigen::Matrix<Eigen::Matrix<T, m, 1>, U, 1>& y_1_2_3
) {
  // y = f(x)
  F(x_values, y_values);

  // Y^{(2)} = F' * X^{(2)}
  for (int j = 0; j < m; j++) {
    T sum = 0;
    for (int i = 0; i < n; i++) {
      sum += y_x(j)(i) * x_2(i, 0);
    }
    y_2(j, 0) = sum;
  }

  // X_{(1)} += F' * Y_{(1)}
  for (int i = 0; i < n; i++) {
    T sum = 0;
    for (int j = 0; j < m; j++) {
      sum += y_x(j)(i) * y_1(0, j);
    }
    x_1(0, i) = sum;
  }

  // X_{(1)}^{(2)} += F'' * Y_{(1)} * X^{(2)} + F' * Y_{(1)}^{(2)}
  for (int i1 = 0; i1 < n; i1++) {
    T term1 = 0;
    T term2 = 0;
    for (int j = 0; j < m; j++) {
      term2 += y_x(j)(i1) * y_1_2(0)(j, 0);
      for (int i2 = 0; i2 < n; i2++) {
        term1 += y_xx(j)(i1)(i2) * y_1(0, j) * x_2(i2, 0);
      }
    }
    x_1_2(0)(i1, 0) = term1 + term2;
  }

  // X^{(2)}_{(3)} += F'' * Y_{(1)} * X_{(1, 3)}^{(2)} + F' * Y^{(2)}_{(3)}
  for (int u = 0; u < U; u++) {
    for (int i2 = 0; i2 < n; i2++) {
      T term1 = 0;
      T term2 = 0;
      for (int j = 0; j < m; j++) {
        term2 += y_x(j)(i2) * y_2_3(u)(j, 0);
        for (int i1 = 0; i1 < n; i1++) {
          term1 += y_xx(j)(i1)(i2) * y_1(0, j) * x_1_2_3(u, 0)(i1, 0);
        }
      }
      x_2_3(u)(i2, 0) = term1 + term2;
    }
  }

  // Y_{(1, 3)} = F'' * X^{(2)} * X_{(1, 3)}^{(2)} + F' * X_{(1, 3)}
  for (int u = 0; u < U; u++) {
    for (int j = 0; j < m; j++) {
      T term1 = 0;
      T term2 = 0;
      for (int i1 = 0; i1 < n; i1++) {
        term2 += y_x(j)(i1) * x_1_3(u)(0, i1);
        for (int i2 = 0; i2 < n; i2++) {
          term1 += y_xx(j)(i1)(i2) * x_2(i2, 0) * x_1_2_3(u, 0)(i1, 0);
        }
      }
      y_1_3(u)(0, j) = term1 + term2;
    }
  }

  // Y_{(1, 3)}^{(2)} = F' * X_{(1, 3)}^{(2)}
  for (int u = 0; u < U; u++) {
    for (int j = 0; j < m; j++) {
      T sum = 0;
      for (int i = 0; i < n; i++) {
        sum += y_x(j)(i) * x_1_2_3(u, 0)(i, 0);
      }
      y_1_2_3(u, 0)(j, 0) = sum;
    }
  }

  // X_{(3)} += F' * Y_{(3)}
  // + F''' * Y_{(1)} * X^{(2)} * X_{(1, 3)}^{(2)}
  // + F'' * X^{(2)} * Y^{(2)}_{(3)}
  // + F'' * Y_{(1)} * X_{(1, 3)}
  // + F'' * Y_{(1)}^{(2)} * X_{(1, 3)}^{(2)}
  for (int u = 0; u < U; u++) {
    for (int i3 = 0; i3 < n; i3++) {
      T term1 = 0;
      T term2 = 0;
      T term3 = 0;
      T term4 = 0;
      T term5 = 0;
      for (int j = 0; j < m; j++) {
        term1 += y_x(j)(i3) * y_3(u, j);
        for (int i1 = 0; i1 < n; i1++) {
          term4 += y_xx(j)(i1)(i3) * y_1(0, j) * x_1_3(u)(0, i1);
          term5 += y_xx(j)(i1)(i3) * y_1_2(0)(j, 0) * x_1_2_3(u, 0)(j, 0);
          for (int i2 = 0; i2 < n; i2++) {
            term2 += y_xxx(j)(i1)(i2)(i3) * y_1(0, j) * x_2(i2, 0) * x_1_2_3(u, 0)(i1, 0);
          }
        }
        for (int i2 = 0; i2 < n; i2++) {
          term3 += y_xx(j)(i2)(i3) * x_2(i2, 0) * y_2_3(u)(j, 0);
        }
      }
      x_3(u, i3) = term1 + term2 + term3 + term4 + term5;
    }
  }
}


template<typename T, int U>
bool Validate_vastsa(std::ofstream& out) {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 8.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;
  Y_XXX_t<T> y_xxx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, 1, n> AD_x_1;
  Eigen::Matrix<T, 1, n> Formula_x_1;
  Eigen::Matrix<T, m, 1> AD_y_2;
  Eigen::Matrix<T, m, 1> Formula_y_2;
  Eigen::Matrix<T, U, n> AD_x_3;
  Eigen::Matrix<T, U, n> Formula_x_3;
  Eigen::Vector<Eigen::Matrix<T, n, 1>, 1> AD_x_1_2;
  Eigen::Vector<Eigen::Matrix<T, n, 1>, 1> Formula_x_1_2;
  Eigen::Vector<Eigen::Matrix<T, 1, m>, U> AD_y_1_3;
  Eigen::Vector<Eigen::Matrix<T, 1, m>, U> Formula_y_1_3;
  Eigen::Vector<Eigen::Matrix<T, n, 1>, U> AD_x_2_3;
  Eigen::Vector<Eigen::Matrix<T, n, 1>, U> Formula_x_2_3;
  Eigen::Matrix<Eigen::Matrix<T, m, 1>, U, 1> AD_y_1_2_3;
  Eigen::Matrix<Eigen::Matrix<T, m, 1>, U, 1> Formula_y_1_2_3;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, 1, m> y_1 = Eigen::Matrix<T, 1, m>::Random().cwiseAbs();
  Eigen::Matrix<T, n, 1> x_2 = Eigen::Matrix<T, n, 1>::Random().cwiseAbs();
  Eigen::Matrix<T, U, m> y_3 = Eigen::Matrix<T, U ,m>::Random().cwiseAbs();
  Eigen::Vector<Eigen::Matrix<T, m, 1>, 1> y_1_2;
  Eigen::Vector<Eigen::Matrix<T, m, 1>, U> y_2_3;
  Eigen::Vector<Eigen::Matrix<T, 1, n>, U> x_1_3;
  Eigen::Matrix<Eigen::Matrix<T, n, 1>, U, 1> x_1_2_3;
  y_1_2(0) = Eigen::Matrix<T, m, 1>::Random().cwiseAbs();
  for (int u = 0; u < U; u++) {
    y_2_3(u) = Eigen::Matrix<T, m, 1>::Random().cwiseAbs();
    x_1_3(u) = Eigen::Matrix<T, 1, n>::Random().cwiseAbs();
    x_1_2_3(u, 0) = Eigen::Matrix<T, n, 1>::Random().cwiseAbs();
  }

  out << "\n=== Testing Third Derivative (Adjoint over Tangent over Adjoint Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for Y_{(1)}: \n" << y_1 << "\n\n";
  out << "Seed for X^{(2)}: \n" << x_2 << "\n\n";
  out << "Seed for Y_{(3)}: \n" << y_3 << "\n\n";
  out << "Seed for Y_{(1)}^{(2)}: (nested matrices not printed)\n\n";
  out << "Seed for Y^{(2)}_{(3)}: (nested matrices not printed)\n\n";
  out << "Seed for X_{(1, 3)}: (nested matrices not printed)\n\n";
  out << "Seed for X_{(1, 3)}^{(2)}: (nested matrices not printed)\n\n";

  // Populate y_x, y_xx and y_xxx
  vastsa_F_xxx(x_values, y_x, y_xx, y_xxx);

  // Run the AD version
  AD_F_xxx(
    x_values, y_1, x_2, y_3, 
    y_1_2, y_2_3, x_1_3,
    x_1_2_3,
    AD_y_values, 
    AD_x_1, AD_y_2, AD_x_3,
    AD_x_1_2, AD_y_1_3, AD_x_2_3,
    AD_y_1_2_3
  );

  // Run the Formula version
  Formula_F_xxx(y_x, y_xx, y_xxx, 
    x_values, y_1, x_2, y_3, 
    y_1_2, y_2_3, x_1_3,
    x_1_2_3,
    Formula_y_values, 
    Formula_x_1, Formula_y_2, Formula_x_3,
    Formula_x_1_2, Formula_y_1_3, Formula_x_2_3,
    Formula_y_1_2_3
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
    diff = std::abs(AD_y_2(j, 0) - Formula_y_2(j, 0));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      out << "Validation for Y^{(2)} Failed\n";
      out << "Validation failed at index " << j << "," << 0 << " Diff: " << diff << "\n";
      return false;
    }

    for (int u = 0; u < U; u++) {
      // Compare y_1_3
      diff = std::abs(AD_y_1_3(u)(0, j) - Formula_y_1_3(u)(0, j));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y_{(1, 3)} Failed\n";
        out << "Validation failed at index " << u << "," << 0 << "," << j << " Diff: " << diff << "\n";
        return false;
      }

      // Compare y_1_2_3
      diff = std::abs(AD_y_1_2_3(u, 0)(j, 0) - Formula_y_1_2_3(u, 0)(j, 0));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y_{(1, 3)}^{(2)} Failed\n";
        out << "Validation failed at index " << u << "," << 0 << "," << j << "," << 0 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }

  for (int i = 0; i < n; i++) {
    // Compare x_1
    diff = std::abs(AD_x_1(0, i) - Formula_x_1(0, i));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      out << "Validation for X_{(1)} Failed\n";
      out << "Validation failed at index " << 0 << "," << i << " Diff: " << diff << "\n";
      return false;
    }

    // Compare x_1_2
    diff = std::abs(AD_x_1_2(0)(i, 0) - Formula_x_1_2(0)(i, 0));
    maxDiff = std::max(maxDiff, diff);
    if (diff > tol) {
      out << "Validation for X_{(1)}^{(2)} Failed\n";
      out << "Validation failed at index " << 0 << "," << i << "," << 0 << " Diff: " << diff << "\n";
      return false;
    }

    // Compare x_3
    for (int u = 0; u < U; u++ ) {
      diff = std::abs(AD_x_3(u, i) - Formula_x_3(u, i));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for X_{(3)} Failed\n";
        out << "Validation failed at index " << u << "," << i << " Diff: " << diff << "\n";
        return false;
      }

      // Compare x_2_3
      diff = std::abs(AD_x_2_3(u)(i, 0) - Formula_x_2_3(u)(i, 0));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for X^{(2)}_{(3)} Failed\n";
        out << "Validation failed at index " << u << "," << i << "," << 0 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}