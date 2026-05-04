// Hessian by vector adjoint of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void vava_F_xx(const X_t<T>& x_values, Y_X_t<T>& y_x, Y_XX_t<T>& y_xx) {
  for (int m1=0;m1<m;++m1) { // shared trailing index
    // activate x and y
    X_t<A_t<A_t<T,n>,m>> x; Y_t<A_t<A_t<T,n>,m>> y;
    for (int i=0;i<n;++i) {
      // set x
      x(i).value().value()=x_values(i);
      // register x with A_t<T,n>::tape
      x(i).value().register_input(); 
      // register x with A_t<A_t<T,n>,m>::tape
      x(i).register_input(); 
    }
    // run overloaded F
    F(x,y);
    // allocate vector of adjoints of A_t<A_t<T,n>,m>::tape
    A_t<A_t<T,n>,m>::tape::init_adjoints();
    // set Y_{(1)}
    y(m1).adjoint(m1).value()=1;
    // interpret A_t<A_t<T,n>,m>::tape
    A_t<A_t<T,n>,m>::tape::interpret();
    // allocate vector of adjoints of A_t<T,n>::tape
    A_t<T,n>::tape::init_adjoints();
    // set X_{(1,2)}
    for (int m2=0;m2<n;++m2) 
      x(m2).adjoint(m1).adjoint(m2)=1;
    // interpret A_t<T,n>::tape
    A_t<T,n>::tape::interpret();
    // extract (n x n)-slice of Hessian from X_{(2)}
    for (int m2=0;m2<n;++m2) 
      for (int i2=0;i2<n;++i2) 
        y_xx(m1)(m2)(i2)=x(i2).value().adjoint(m2);
      
    // extract Jacobian
    for (int i = 0; i < n; i++) {
      y_x(m1)(i) = x(i).adjoint(m1).value();
      std::cout << x(i).adjoint(m1).value() << "\n";
    }
    // reset both tapes for next m1
    A_t<A_t<T,n>,m>::tape::reset();
    A_t<T,n>::tape::reset();
  }
}	


template<typename T, int U1, int U2>
void A_F_xx(
 X_t<A_t<A_t<T, U2>, U1>>& x, Y_t<A_t<A_t<T, U2>, U1>>& y,
 Eigen::Matrix<T, U1, m>& y_1, 
 Y_t<T>& y_values, Eigen::Matrix<T, U1, n>& x_1
) {
  for (int i = 0; i < n; i++) {
    x(i).register_input();
  }

  F(x, y);

  A_t<A_t<T, U2>, U1>::tape::init_adjoints();

  for (int j = 0; j < m; j++) {
    for (int u1 = 0; u1 < U1; u1++) {
      // Seed Y_{(1)}
      y(j).adjoint(u1).value() = y_1(u1, j);
    }
  }

  A_t<A_t<T, U2>, U1>::tape::interpret();

  for (int j = 0; j < m; j++) {
    // Extract y
    y_values(j) = y(j).value().value();
  }

  for (int i = 0; i < n; i++) {
    for (int u1 = 0; u1 < U1; u1++) {
      // Extract X_{(1)}
      x_1(u1, i) = x(i).adjoint(u1).value();
    }
  }

}


template<typename T, int U1, int U2>
void AD_F_xx(const X_t<T>& x_values, 
  Eigen::Matrix<T, U1, m>& y_1, Eigen::Matrix<T, U2, m>& y_2, 
  Eigen::Vector<Eigen::Matrix<T, U1, n>, U2>& x_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, U1, n>& x_1, Eigen::Matrix<T, U2, n>& x_2,
  Eigen::Vector<Eigen::Matrix<T, U1, m>, U2>& y_1_2) {

  X_t<A_t<A_t<T, U2>, U1>> x; 
  Y_t<A_t<A_t<T, U2>, U1>> y;

  for (int i = 0; i < n; i++) {
    // Set x
    x(i).value().value() = x_values(i);

    x(i).value().register_input();
  }

  for (int j = 0; j < m; j++) {
    for (int u1 = 0; u1 < U1; u1++) {
      y(j).adjoint(u1).register_input();
    }
  }

  // First order adjoint
  A_F_xx<T, U1, U2>(x, y, y_1, y_values, x_1);

  A_t<T, U2>::tape::init_adjoints();

  for (int u2 = 0; u2 < U2; u2++) {
    for (int j = 0; j < m; j++) {
      // Seed Y_{(2)}
      y(j).value().adjoint(u2) = y_2(u2, j);
    }

    for (int u1 = 0; u1 < U1; u1++) {
      for (int i = 0; i < n; i++) {
        // Seed X_{(1, 2)}
        x(i).adjoint(u1).adjoint(u2) = x_1_2(u2)(u1, i);
      }
    }
  }

  A_t<T, U2>::tape::interpret();

  for (int u2 = 0; u2 < U2; u2++) {
    for (int i = 0; i < n; i++) {
      // Extract X_{(2)}
      x_2(u2, i) = x(i).value().adjoint(u2);
    }

    for (int j = 0; j < m; j++) {
      for (int u1 = 0; u1 < U1; u1++) {
        // Extract Y_{(1, 2)}
        y_1_2(u2)(u1, j) = y(j).adjoint(u1).adjoint(u2);
      }
    }
  }

  A_t<A_t<T, U2>, U1>::tape::reset();
  A_t<T, U2>::tape::reset();
}


template<typename T, int U1, int U2>
void Formula_F_xx(const Y_X_t<T>& y_x, const Y_XX_t<T>& y_xx,
  const X_t<T>& x_values, 
  const Eigen::Matrix<T, U1, m>& y_1, const Eigen::Matrix<T, U2, m>& y_2, 
  const Eigen::Vector<Eigen::Matrix<T, U1, n>, U2>& x_1_2, 
  Y_t<T>& y_values, 
  Eigen::Matrix<T, U1, n>& x_1, Eigen::Matrix<T, U2, n>& x_2,
  Eigen::Vector<Eigen::Matrix<T, U1, m>, U2>& y_1_2) {

  // y = f(x)
  F(x_values, y_values);

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

  // Y_{(1, 2)} = F' * X_{(1, 2)}
  for (int u2 = 0; u2 < U2; u2++) {
    for (int u1 = 0; u1 < U1; u1++) {
      for (int j = 0; j < m; j++) {
        T sum = 0;
        for (int i = 0; i < n; i++) {
          sum += y_x(j)(i) * x_1_2(u2)(u1, i);
        }
        y_1_2(u2)(u1, j) = sum;
      }
    }
  }

  // X_{(2)} += F' * Y_{(2)} + F'' * Y_{(1)} * X_{(1, 2)}
  for (int u2 = 0; u2 < U2; u2++) {
    for (int i2 = 0; i2 < n; i2++) {
      T term1 = 0;
      T term2 = 0;
      for (int j = 0; j < m; j++) {
        term1 += y_x(j)(i2) * y_2(u2, j);
        for (int i1 = 0; i1 < n; i1++) {
          for (int u1 = 0; u1 < U1; u1++) {
            term2 += y_xx(j)(i1)(i2) * y_1(u1, j) * x_1_2(u2)(u1, i1);
          }
        }
      }
      x_2(u2, i2) = term1 + term2;
    }
  }
}


template<typename T, int U1, int U2>
bool Validate_vava(std::ofstream& out) {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 4.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, U1, n> AD_x_1;
  Eigen::Matrix<T, U1, n> Formula_x_1;
  Eigen::Matrix<T, U2, n> AD_x_2;
  Eigen::Matrix<T, U2, n> Formula_x_2;
  Eigen::Vector<Eigen::Matrix<T, U1, m>, U2> AD_y_1_2;
  Eigen::Vector<Eigen::Matrix<T, U1, m>, U2> Formula_y_1_2;

  // Inputs
  X_t<T> x_values = X_t<T>::Random();
  Eigen::Matrix<T, U1, m> y_1 = Eigen::Matrix<T, U1, m>::Random().cwiseAbs();
  Eigen::Matrix<T, U2, m> y_2 = Eigen::Matrix<T, U2, m>::Random().cwiseAbs();
  Eigen::Vector<Eigen::Matrix<T, U1, n>, U2> x_1_2;
  for (int u2 = 0; u2 < U2; u2++) {
    x_1_2(u2) = Eigen::Matrix<T, U1, n>::Random().cwiseAbs();
  }

  out << "\n=== Testing Second Derivative (Adjoint over Adjoint Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for Y_({1}): \n" << y_1 << "\n\n";
  out << "Seed for Y_({2}): \n" << y_2 << "\n\n";
  out << "Seed for X_{(1, 2)}: (nested matrices not printed)\n\n";

  // Populate y_x and y_xx
  vava_F_xx(x_values, y_x, y_xx);

  // Run the AD version
  AD_F_xx(x_values, y_1, y_2, x_1_2, 
          AD_y_values, AD_x_1, AD_x_2, AD_y_1_2);

  // Run the Formula version
  Formula_F_xx(y_x, y_xx, 
    x_values, y_1, y_2, x_1_2, 
              Formula_y_values, Formula_x_1, Formula_x_2, Formula_y_1_2);


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
  }

  // Compare x_1
  for (int i = 0; i < n; i++) {
    for (int u1 = 0; u1 < U1; u1++) {
      diff = std::abs(AD_x_1(u1, i) - Formula_x_1(u1, i));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for X_{(1)} Failed\n";
        out << "Validation failed at index " << u1 << "," << i << " Diff: " << diff << "\n";
        return false;
      }
    }
  }
  
  // Compare x_2
  for (int i = 0; i < n; i++) {
    for (int u2 = 0; u2 < U2; u2++) {
      diff = std::abs(AD_x_2(u2, i) - Formula_x_2(u2, i));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for X_{(2)} Failed\n";
        out << "Validation failed at index " << u2 << "," << i << " Diff: " << diff << "\n";
        return false;
      }
    }
  }

  // Compare y_1_2
  for (int u2 = 0; u2 < U2; u2++) {
    for (int u1 = 0; u1 < U1; u1++) {
      for (int j = 0; j < m; j++) {
        std::cout << "AD_y_1_2 " << u1 << ", " << j << ": " << AD_y_1_2(u2)(u1, j) << "\n";
        std::cout << "Formula_y_1_2 " << u1 << ", " << j << ": " << Formula_y_1_2(u2)(u1, j) << "\n";
        std::cout << "\n";

        // diff = std::abs(AD_y_1_2(u2)(u1, j) - Formula_y_1_2(u2)(u1, j));
        // maxDiff = std::max(maxDiff, diff);
        // if (diff > tol) {
        //   out << "Validation for Y_{(1, 2)} Failed\n";
        //   out << "Validation failed at index " << u2 << "," << u1 << "," << j << " Diff: " << diff << "\n";
        //   return false;
        // }
      }
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}