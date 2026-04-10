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
}	


template<typename T, int U1, int V2, int V3>
void AD_F_xx(
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

}


template<typename T, int U1, int V2, int V3>
void Formula_F_xx(
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

}


template<typename T, int U1, int V2, int V3>
bool Validate_vtvtva(std::ofstream out) {
  return true;
}