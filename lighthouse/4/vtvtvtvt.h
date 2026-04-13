// Hessian by vector tangent of vector tangent of vector tangent AD

#include "F.h"
#include "ad.h"
#include "ad_types.h"

template<typename T>
void vtvtvtvt_F_xxxx(const X_t<T>& x_values, 
  Y_X_t<T>& y_x, 
  Y_XX_t<T>& y_xx, 
  Y_XXX_t<T>& y_xxx, 
  Y_XXXX_t<T>& y_xxxx) {
  // activate x and y
  X_t<T_t<T_t<T_t<T_t<T,n>,n>,n>,n>> x; 
  Y_t<T_t<T_t<T_t<T_t<T,n>,n>,n>,n>> y;
  for (int i=0;i<n;++i) {
    // set x
    x(i).value().value().value().value()=x_values(i);
    // set X^{(1)}
    x(i).tangent(i).value().value().value()=1;
    // set X^{(2)}
    x(i).value().tangent(i).value().value()=1;
    // set X^{(3)}
    x(i).value().value().tangent(i).value()=1;
    // set X^{(4)}
    x(i).value().value().value().tangent(i)=1;
  }
  // run overloaded F
  F(x,y);
  // extract F^{[4]} from Y^{(1,2,3,4)}
  for (int j=0;j<m;++j) 
    for (int n1=0;n1<n;++n1) {
      y_x(j)(n1) = y(j).tangent(n1).value().value().value();
      for (int n2=0;n2<n;++n2) {
        y_xx(j)(n1)(n2) = y(j).tangent(n1).tangent(n2).value().value();
        for (int n3=0;n3<n;++n3) {
          y_xxx(j)(n1)(n2)(n3) = y(j).tangent(n1).tangent(n2).tangent(n3).value();
          for (int n4=0;n4<n;++n4) 
            y_xxxx(j)(n1)(n2)(n3)(n4)=y(j).tangent(n1).tangent(n2).tangent(n3).tangent(n4);
        }
      }
    }
}	


template<typename T, int V1, int V2, int V3, int V4>
void AD_F_xxxx(
  const X_t<T>& x_values,
  Eigen::Matrix<T, n, V1>& x_1, Eigen::Matrix<T, n, V2>& x_2, Eigen::Matrix<T, n, V3>& x_3, Eigen::Matrix<T, n, V4>& x_4, 
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, n>& x_1_2, Eigen::Vector<Eigen::Matrix<T, V1, V3>, n>& x_1_3, Eigen::Vector<Eigen::Matrix<T, V1, V4>, n>& x_1_4, 
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, n>& x_2_3, Eigen::Vector<Eigen::Matrix<T, V2, V4>, n>& x_2_4, 
  Eigen::Vector<Eigen::Matrix<T, V3, V4>, n>& x_3_4, 
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, n, V1>& x_1_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V4>, n, V1>& x_1_2_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V3>, n, V1>& x_1_3_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, n, V2>& x_2_3_4,
  Eigen::Vector<Eigen::Matrix<Eigen::Matrix<T, V3, V4>, V1, V2>, n>& x_1_2_3_4,
  Y_t<T>& y_values,
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, m, V2>& y_2, Eigen::Matrix<T, m, V3>& y_3, Eigen::Matrix<T, m, V4>& y_4,
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m>& y_1_2, Eigen::Vector<Eigen::Matrix<T, V1, V3>, m>& y_1_3, Eigen::Vector<Eigen::Matrix<T, V1, V4>, m>& y_1_4,
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m>& y_2_3, Eigen::Vector<Eigen::Matrix<T, V2, V4>, m>& y_2_4,
  Eigen::Vector<Eigen::Matrix<T, V3, V4>, m>& y_3_4,
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1>& y_1_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V4>, m, V1>& y_1_2_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V1>& y_1_3_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V2>& y_2_3_4,
  Eigen::Vector<Eigen::Matrix<Eigen::Matrix<T, V3, V4>, V1, V2>, m>& y_1_2_3_4
) {
  X_t<T_t<T_t<T_t<T_t<T,V4>,V3>,V2>,V1>> x; 
  Y_t<T_t<T_t<T_t<T_t<T,V4>,V3>,V2>,V1>> y;

  for (int i = 0; i < n; i++) {
    // Seed x
    x(i).value().value().value().value() = x(i);

    // Seed X^{(1)}
    for (int v1 = 0; v1 < V1; v1++) {
      x(i).tangent(v1).value().value().value() = x_1(i, v1);

      // Seed X^{(1, 2)}
      for (int v2 = 0; v2 < V2; v2++) {
        x(i).tangent(v1).tangent(v2).value().value() = x_1_2(i)(v1, v2);  

        // Seed X^{(1, 2, 3)}
        for (int v3 = 0; v3 < V3; v3++) {
          x(i).tangent(v1).tangent(v2).tangent(v3).value() = x_1_2_3(i, v1)(v2, v3);

          // Seed X^{(1, 2, 3, 4)}
          for (int v4 = 0; v4 < V4; v4++) {
            x(i).tangent(v1).tangent(v2).tangent(v3).tangent(v4) = x_1_2_3_4(i)(v1, v2)(v3, v4);
          }
        }

        // Seed X^{(1, 2, 4)}
        for (int v4 = 0; v4 < V4; v4++) {
          x(i).tangent(v1).tangent(v2).value().tangent(v4) = x_1_2_4(i, v1)(v2, v4);
        }
      }

      // Seed X^{(1, 3)}
      for (int v3 = 0; v3 < V3; v3++) {
        x(i).tangent(v1).value().tangent(v3).value() = x_1_3(i)(v1, v3);

        // Seed X^{(1, 3, 4)}
        for (int v4 = 0; v4 < V4; v4++) {
          x(i).tangent(v1).value().tangent(v3).tangent(v4) = x_1_3_4(i, v1)(v3, v4);
        }
      }

      // Seed X^{(1, 4)}
      for (int v4 = 0; v4 < V4; v4++) {
        x(i).tangent(v1).value().value().tangent(v4) = x_1_4(i)(v1, v4);  
      }

    }
    // Seed X^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      x(i).value().tangent(v2).value().value() = x_2(i, v2);

      // Seed X^{(2, 3)}
      for (int v3 = 0; v3 < V3; v3++) {
        x(i).value().tangent(v2).tangent(v3).value() = x_2_3(i)(v2, v3);

        // Seed X^{(2, 3, 4)}
        for (int v4 = 0; v4 < V4; v4++) {
          x(i).value().tangent(v2).tangent(v3).tangent(v4) = x_2_3_4(i, v2)(v3, v4);
        }
      }

      // Seed X^{(2, 4)}
      for (int v4 = 0; v4 < V4; v4++) {
        x(i).value().tangent(v2).value().tangent(v4) = x_2_4(i)(v2, v4);
      }
    }
    // Seed X^{(3)}
    for (int v3 = 0; v3 < V3; v3++) {
      x(i).value().value().tangent(v3).value() = x_3(i, v3);

      // Seed X^{(3, 4)}
      for (int v4 = 0; v4 < V4; v4++) {
        x(i).value().value().tangent(v3).tangent(v4) = x_3_4(i)(v3, v4);
      }

    }

    // Seed X^{(4)}
    for (int v4 = 0; v4 < V4; v4++) {
      x(i).value().value().value().tangent(v4) = x_4(i, v4);
    }
  }

  F(x,y);

  // Extract y
  for (int j = 0; j < m; j++) {
    y_values(j) = y(j).value().value().value().value();

    // Extract Y^{(1)}
    for (int v1 = 0; v1 < V1; v1++) {
      y_1(j, v1) = y(j).tangent(v1).value().value().value();

      // Extract Y^{(1, 2)}
      for (int v2 = 0; v2 < V2; v2++) {
        y_1_2(j)(v1, v2) = y(j).tangent(v1).tangent(v2).value().value();

        // Extract Y^{(1, 2, 3)}
        for (int v3 = 0; v3 < V3; v3++) {
          y_1_2_3(j, v1)(v2, v3) = y(j).tangent(v1).tangent(v2).tangent(v3).value();

          // Extract Y^{(1, 2, 3, 4)}
          for (int v4 = 0; v4 < V4; v4++) {
            y_1_2_3_4(j)(v1, v2)(v3, v4) = y(j).tangent(v1).tangent(v2).tangent(v3).tangent(v4);
          }
        }

        // Extract Y^{(1, 2, 4)}
        for (int v4 = 0; v4 < V4; v4++) {
          y_1_2_4(j, v1)(v2, v4) = y(j).tangent(v1).tangent(v2).value().tangent(v4);
        }
      }

      // Extract Y^{(1, 3)}
      for (int v3 = 0; v3 < V3; v3++) {
        y_1_3(j)(v1, v3) = y(j).tangent(v1).value().tangent(v3).value();

        // Extract Y^{(1, 3, 4)}
        for (int v4 = 0; v4 < V4; v4++) {
          y_1_3_4(j, v1)(v3, v4) = y(j).tangent(v1).value().tangent(v3).tangent(v4);
        }
      }

      // Extract Y^{(1, 4)}
      for (int v4 = 0; v4 < V4; v4++) {
        y_1_4(j)(v1, v4) = y(j).tangent(v1).value().value().tangent(v4);
      }
    }

    // Extract Y^{(2)}
    for (int v2 = 0; v2 < V2; v2++) {
      y_2(j, v2) = y(j).value().tangent(v2).value().value();

      // Extract Y^{(2, 3)}
      for (int v3 = 0; v3 < V3; v3++) {
        y_2_3(j)(v2, v3) = y(j).value().tangent(v2).tangent(v3).value();

        // Extract Y^{(2, 3, 4)}
        for (int v4 = 0; v4 < V4; v4++) {
          y_2_3_4(j, v2)(v3, v4) = y(j).value().tangent(v2).tangent(v3).tangent(v4);
        }
      }
      // Extract Y^{(2, 4)}
      for (int v4 = 0; v4 < V4; v4++) {
        y_2_4(j)(v2, v4) = y(j).value().tangent(v2).value().tangent(v4);
      }
    }

    // Extract Y^{(3)}
        for (int v3 = 0; v3 < V3; v3++) {
      y_3(j, v3) = y(j).value().value().tangent(v3).value();

      // Extract Y^{(3, 4)}
      for (int v4 = 0; v4 < V4; v4++) {
        y_3_4(j)(v3, v4) = y(j).value().value().tangent(v3).tangent(v4);
      }
    }

    // Extract Y^{(4)}
        for (int v4 = 0; v4 < V4; v4++) {
      y_4(j, v4) = y(j).value().value().value().tangent(v4);
    }
  }
}


template<typename T, int V1, int V2, int V3, int V4>
void Formula_F_xxxx(
  Y_X_t<T>& y_x, Y_XX_t<T>& y_xx, Y_XXX_t<T>& y_xxx, Y_XXXX_t<T>& y_xxxx,
  const X_t<T>& x_values,
  Eigen::Matrix<T, n, V1>& x_1, Eigen::Matrix<T, n, V2>& x_2, Eigen::Matrix<T, n, V3>& x_3, Eigen::Matrix<T, n, V4>& x_4, 
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, n>& x_1_2, Eigen::Vector<Eigen::Matrix<T, V1, V3>, n>& x_1_3, Eigen::Vector<Eigen::Matrix<T, V1, V4>, n>& x_1_4, 
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, n>& x_2_3, Eigen::Vector<Eigen::Matrix<T, V2, V4>, n>& x_2_4, 
  Eigen::Vector<Eigen::Matrix<T, V3, V4>, n>& x_3_4, 
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, n, V1>& x_1_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V4>, n, V1>& x_1_2_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V3>, n, V1>& x_1_3_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, n, V2>& x_2_3_4,
  Eigen::Vector<Eigen::Matrix<Eigen::Matrix<T, V3, V4>, V1, V2>, n>& x_1_2_3_4,
  Y_t<T>& y_values,
  Eigen::Matrix<T, m, V1>& y_1, Eigen::Matrix<T, m, V2>& y_2, Eigen::Matrix<T, m, V3>& y_3, Eigen::Matrix<T, m, V4>& y_4,
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m>& y_1_2, Eigen::Vector<Eigen::Matrix<T, V1, V3>, m>& y_1_3, Eigen::Vector<Eigen::Matrix<T, V1, V4>, m>& y_1_4,
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m>& y_2_3, Eigen::Vector<Eigen::Matrix<T, V2, V4>, m>& y_2_4,
  Eigen::Vector<Eigen::Matrix<T, V3, V4>, m>& y_3_4,
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1>& y_1_2_3,
  Eigen::Matrix<Eigen::Matrix<T, V2, V4>, m, V1>& y_1_2_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V1>& y_1_3_4,
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V2>& y_2_3_4,
  Eigen::Vector<Eigen::Matrix<Eigen::Matrix<T, V3, V4>, V1, V2>, m>& y_1_2_3_4
) {
  // y = f(x)
  F(x_values, y_values);

  for (int j = 0; j < m; j++) {
    // Y^{(1)} = F' * X^{(1)}
    for (int v1 = 0; v1 < V1; v1++ ) {
      T sum = 0;
      for (int i1 = 0; i1 < n; i1++) {
          sum += y_x(j, i1) * x_1(i1, v1);
      }
      y_1(j, v1) = sum;
    }

    // Y^{(2)} = F' * X^{(2)}
    for (int v2 = 0; v2 < V2; v2++ ) {
      T sum = 0;
      for (int i2 = 0; i2 < n; i2++) {
          sum += y_x(j, i2) * x_2(i2, v2);
      }
      y_2(j, v2) = sum;
    }

    // Y^{(3)} = F' * X^{(3)}
    for (int v3 = 0; v3 < V3; v3++ ) {
      T sum = 0;
      for (int i3 = 0; i3 < n; i3++) {
          sum += y_x(j, i3) * x_3(i3, v3);
      }
      y_3(j, v3) = sum;
    }

    // Y^{(4)} = F' * X^{(4)}
    for (int v4 = 0; v4 < V4; v4++ ) {
      T sum = 0;
      for (int i4 = 0; i4 < n; i4++) {
          sum += y_x(j, i4) * x_4(i4, v4);
      }
      y_4(j, v4) = sum;
    }
  }

  // Y^{(1, 2)} = F'' * X^{(1)} * X^{(2)} + F' * X^{(1, 2)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        
      }
    }
  }

  // Y^{(1, 3)} = F'' * X^{(1)} * X^{(3)} + F' * X^{(1, 3)}

  // Y^{(1, 4)} = F'' * X^{(1)} * X^{(4)} + F' * X^{(1, 4)}

  // Y^{(2, 3)} = F'' * X^{(2)} * X^{(3)} + F' * X^{(2, 3)}

  // Y^{(2, 4)} = F'' * X^{(2)} * X^{(4)} + F' * X^{(2, 4)}

  // Y^{(3, 4)} = F'' * X^{(3)} * X^{(4)} + F' * X^{(3, 4)}

  // Y^{(1, 2, 3)} = F''' * X^{(1)} * X^{(2)} * X^{(3)}
  // + F'' * X^{(2)} * X^{(1, 3)}
  // + F'' * X^{(1)} * X^{(2, 3)}
  // + F'' * X^{(3)} * X^{(1, 2)}
  // + F' * X^{(1, 2, 3)}

  // Y^{(1, 2, 4)} = F''' * X^{(1)} * X^{(2)} * X^{(4)}
  // + F'' * X^{(2)} * X^{(1, 4)}
  // + F'' * X^{(1)} * X^{(2, 4)}
  // + F'' * X^{(4)} * X^{(1, 2)}
  // + F' * X^{(1, 2, 4)}

  // Y^{(2, 3, 4)} = F''' * X^{(2)} * X^{(3)} * X^{(4)}
  // + F'' * X^{(2)} * X^{(3, 4)}
  // + F'' * X^{(3)} * X^{(2, 4)}
  // + F'' * X^{(4)} * X^{(2, 3)}
  // + F' * X^{(2, 3, 4)}

  // Y^{(1, 2, 3, 4)}
  // = F'''' * X^{(1)} * X^{(2)} * X^{(3)} * X^{(4)}
  // + F''' * X^{(1)} * X^{(2)} * X^{(3, 4)}
  // + F''' * X^{(1)} * X^{(3)} * X^{(2, 4)}
  // + F''' * X^{(1)} * X^{(4)} * X^{(2, 3)}
  // + F''' * X^{(2)} * X^{(3)} * X^{(1, 4)}
  // + F''' * X^{(2)} * X^{(4)} * X^{(1, 3)}
  // + F''' * X^{(3)} * X^{(4)} * X^{(1, 4)}
  // + F'' * X^{(1)} * X^{(2, 3, 4)}
  // + F'' * X^{(1, 2)} * X^{(3, 4)}
  // + F'' * X^{(1, 4)} * X^{(2, 3)}
  // + F'' * X^{(1, 2, 3)} * X^{(4)}
  // + F'' * X^{(2)} * X^{(1, 3, 4)}
  // + F'' * X^{(2, 4)} * X^{(1, 3)}
  // + F'' * X^{(3)} * X^{(1, 2, 4)}
  // + F' * X^{(1, 2, 3, 4)}
}