// Third derivative by vector adjoint of vector tangent of vector adjoint AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

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
      for (int i1=0;i1<n;++i1) 
        for (int m3=0;m3<n;++m3) 
          y_xxx(m1)(i1)(n2)(m3)=x(i1).value().value().adjoint(m3);
      // extract Jacobian
      for (int i1 = 0; i1 < n; i1++) {
        y_x(m1)(i1) = x(i1).adjoint(m1).value().value();
      }
      // extract Hessian
      for (int i1 = 0; i1 < n; i1++) {
        y_xx(m1)(i1)(n2) = x(i1).adjoint(m1).tangent(n2).value();
      }
      A_t<T,n>::tape::reset();
      A_t<T_t<A_t<T,n>,1>,1>::tape::reset();
    }	
  }
}


template<typename T, int U>
void AD_F_xxx(
  const X_t<T>& x_values,
  const Eigen::Matrix<T, U, m>& y_1, const Eigen::Matrix<T, n, 1>& x_2, const Eigen::Matrix<T, 1, m>& y_3,
  const Eigen::Vector<Eigen::Matrix<T, m, 1>, U>& y_1_2, const Eigen::Vector<Eigen::Matrix<T, m, 1>, 1>& y_2_3,
  const Eigen::Vector<Eigen::Matrix<T, U, n>, 1>& x_1_3, const Eigen::Matrix<Eigen::Matrix<T, n, 1>, 1, U>& x_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, U, n>& x_1, Eigen::Matrix<T, m, 1>& y_2, Eigen::Matrix<T, U, n>& x_3,
  Eigen::Vector<Eigen::Matrix<T, n, 1>, U>& x_1_2, Eigen::Vector<Eigen::Matrix<T, U, m>, 1>& y_1_3, 
  Eigen::Vector<Eigen::Matrix<T, n, 1>, 1>& x_2_3,
  Eigen::Matrix<Eigen::Matrix<T, m, 1>, 1, U>& y_1_2_3
) {

}


template<typename T, int U>
void Formula_F_xxx(
  Y_X_t<T>& y_x, Y_XX_t<T>& y_xx, Y_XXX_t<T>& y_xxx,
  const X_t<T>& x_values,
  const Eigen::Matrix<T, U, m>& y_1, const Eigen::Matrix<T, n, 1>& x_2, const Eigen::Matrix<T, 1, m>& y_3,
  const Eigen::Vector<Eigen::Matrix<T, m, 1>, U>& y_1_2, const Eigen::Vector<Eigen::Matrix<T, m, 1>, 1>& y_2_3,
  const Eigen::Vector<Eigen::Matrix<T, U, n>, 1>& x_1_3, const Eigen::Matrix<Eigen::Matrix<T, n, 1>, 1, U>& x_1_2_3,
  Y_t<T>& y_values,
  Eigen::Matrix<T, U, n>& x_1, Eigen::Matrix<T, m, 1>& y_2, Eigen::Matrix<T, U, n>& x_3,
  Eigen::Vector<Eigen::Matrix<T, n, 1>, U>& x_1_2, Eigen::Vector<Eigen::Matrix<T, U, m>, 1>& y_1_3, 
  Eigen::Vector<Eigen::Matrix<T, n, 1>, 1>& x_2_3,
  Eigen::Matrix<Eigen::Matrix<T, m, 1>, 1, U>& y_1_2_3
) {

}


template<typename T, int U>
bool Validate_vastsa(std::ofstream& out) {
  return true;
}