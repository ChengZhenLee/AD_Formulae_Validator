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
    x(i).value().value().value().value() = x_values(i);

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
          sum += y_x(j)(i1) * x_1(i1, v1);
      }
      y_1(j, v1) = sum;
    }

    // Y^{(2)} = F' * X^{(2)}
    for (int v2 = 0; v2 < V2; v2++ ) {
      T sum = 0;
      for (int i2 = 0; i2 < n; i2++) {
          sum += y_x(j)(i2) * x_2(i2, v2);
      }
      y_2(j, v2) = sum;
    }

    // Y^{(3)} = F' * X^{(3)}
    for (int v3 = 0; v3 < V3; v3++ ) {
      T sum = 0;
      for (int i3 = 0; i3 < n; i3++) {
          sum += y_x(j)(i3) * x_3(i3, v3);
      }
      y_3(j, v3) = sum;
    }

    // Y^{(4)} = F' * X^{(4)}
    for (int v4 = 0; v4 < V4; v4++ ) {
      T sum = 0;
      for (int i4 = 0; i4 < n; i4++) {
          sum += y_x(j)(i4) * x_4(i4, v4);
      }
      y_4(j, v4) = sum;
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

  // Y^{(1, 4)} = F'' * X^{(1)} * X^{(4)} + F' * X^{(1, 4)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v4 = 0; v4 < V4; v4++) {
        T term1 = 0;
        T term2 = 0;
        for (int i1 = 0; i1 < n; i1++) {
          term2 += y_x(j)(i1) * x_1_4(i1)(v1, v4);
          for (int i4 = 0; i4 < n; i4++) {
            term1 += y_xx(j)(i1)(i4) * x_1(i1, v1) * x_4(i4, v4);
          }
        }
        y_1_4(j)(v1, v4) = term1 + term2;
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

  // Y^{(2, 4)} = F'' * X^{(2)} * X^{(4)} + F' * X^{(2, 4)}
  for (int j = 0; j < m; j++) {
    for (int v2 = 0; v2 < V2; v2++) {
      for (int v4 = 0; v4 < V4; v4++) {
        T term1 = 0;
        T term2 = 0;
        for (int i2 = 0; i2 < n; i2++) {
          term2 += y_x(j)(i2) * x_2_4(i2)(v2, v4);
          for (int i4 = 0; i4 < n; i4++) {
            term1 += y_xx(j)(i2)(i4) * x_2(i2, v2) * x_4(i4, v4);
          }
        }
        y_2_4(j)(v2, v4) = term1 + term2;
      }
    }
  }

  // Y^{(3, 4)} = F'' * X^{(3)} * X^{(4)} + F' * X^{(3, 4)}
  for (int j = 0; j < m; j++) {
    for (int v3 = 0; v3 < V3; v3++) {
      for (int v4 = 0; v4 < V4; v4++) {
        T term1 = 0;
        T term2 = 0;
        for (int i3 = 0; i3 < n; i3++) {
          term2 += y_x(j)(i3) * x_3_4(i3)(v3, v4);
          for (int i4 = 0; i4 < n; i4++) {
            term1 += y_xx(j)(i3)(i4) * x_3(i3, v3) * x_4(i4, v4);
          }
        }
        y_3_4(j)(v3, v4) = term1 + term2;
      }
    }
  }

  // Y^{(1, 2, 3)} = F''' * X^{(1)} * X^{(2)} * X^{(3)}
  // + F'' * X^{(2)} * X^{(1, 3)}
  // + F'' * X^{(1)} * X^{(2, 3)}
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
              term2 += y_xx(j)(i1)(i2) * x_2(i2, v2) * x_1_3(i1)(v1, v3);
              term3 += y_xx(j)(i1)(i2) * x_1(i1, v1) * x_2_3(i2)(v2, v3);
              for (int i3 = 0; i3 < n; i3++) {
                term1 += y_xxx(j)(i1)(i2)(i3) * x_1(i1, v1) * x_2(i2, v2) * x_3(i3, v3);
              }
            }
            for (int i3 = 0; i3 < n; i3++) {
              term4 += y_xx(j)(i1)(i3) * x_1_2(i1)(v1, v2) * x_3(i3, v3);
            }
          }
          y_1_2_3(j, v1)(v2, v3) = term1 + term2 + term3 + term4 + term5;
        }
      }
    }
  }

  // Y^{(1, 2, 4)} = F''' * X^{(1)} * X^{(2)} * X^{(4)}
  // + F'' * X^{(2)} * X^{(1, 4)}
  // + F'' * X^{(1)} * X^{(2, 4)}
  // + F'' * X^{(4)} * X^{(1, 2)}
  // + F' * X^{(1, 2, 4)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        for (int v4 = 0; v4 < V4; v4++) {
          T term1 = 0;
          T term2 = 0;
          T term3 = 0;
          T term4 = 0;
          T term5 = 0;
          for (int i1 = 0; i1 < n; i1++) {
            term5 += y_x(j)(i1) * x_1_2_4(i1, v1)(v2, v4);
            for (int i2 = 0; i2 < n; i2++) {
              term2 += y_xx(j)(i1)(i2) * x_2(i2, v2) * x_1_4(i1)(v1, v4);
              term3 += y_xx(j)(i1)(i2) * x_1(i1, v1) * x_2_4(i2)(v2, v4);
              for (int i4 = 0; i4 < n; i4++) {
                term1 += y_xxx(j)(i1)(i2)(i4) * x_1(i1, v1) * x_2(i2, v2) * x_4(i4, v4);
              }
            }
            for (int i4 = 0; i4 < n; i4++) {
              term4 += y_xx(j)(i1)(i4) * x_1_2(i1)(v1, v2) * x_4(i4, v4);
            }
          }
          y_1_2_4(j, v1)(v2, v4) = term1 + term2 + term3 + term4 + term5;
        }
      }
    }
  }

  // Y^{(1, 3, 4)} = F''' * X^{(1)} * X^{(3)} * X^{(4)}
  // + F'' * X^{(1)} * X^{(3, 4)}
  // + F'' * X^{(3)} * X^{(1, 4)}
  // + F'' * X^{(4)} * X^{(1, 3)}
  // + F' * X^{(1, 3, 4)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v3 = 0; v3 < V3; v3++) {
        for (int v4 = 0; v4 < V4; v4++) {
          T term1 = 0;
          T term2 = 0;
          T term3 = 0;
          T term4 = 0;
          T term5 = 0;
          for (int i1 = 0; i1 < n; i1++) {
            term5 += y_x(j)(i1) * x_1_3_4(i1, v1)(v3, v4);
            for (int i3 = 0; i3 < n; i3++) {
              term2 += y_xx(j)(i1)(i3) * x_1(i1, v1) * x_3_4(i3)(v3, v4);
              term3 += y_xx(j)(i1)(i3) * x_3(i3, v3) * x_1_4(i1)(v1, v4);
              for (int i4 = 0; i4 < n; i4++) {
                term1 += y_xxx(j)(i1)(i3)(i4) * x_1(i1, v1) * x_3(i3, v3) * x_4(i4, v4);
              }
            }
            for (int i4 = 0; i4 < n; i4++) {
              term4 += y_xx(j)(i1)(i4) * x_1_3(i1)(v1, v3) * x_4(i4, v4);
            }
          }
          y_1_3_4(j, v1)(v3, v4) = term1 + term2 + term3 + term4 + term5;
        }
      }
    }
  }

  // Y^{(2, 3, 4)} = F''' * X^{(2)} * X^{(3)} * X^{(4)}
  // + F'' * X^{(2)} * X^{(3, 4)}
  // + F'' * X^{(3)} * X^{(2, 4)}
  // + F'' * X^{(4)} * X^{(2, 3)}
  // + F' * X^{(2, 3, 4)}
  for (int j = 0; j < m; j++) {
    for (int v2 = 0; v2 < V2; v2++) {
      for (int v3 = 0; v3 < V3; v3++) {
        for (int v4 = 0; v4 < V4; v4++) {
          T term1 = 0;
          T term2 = 0;
          T term3 = 0;
          T term4 = 0;
          T term5 = 0;
          for (int i2 = 0; i2 < n; i2++) {
            term5 += y_x(j)(i2) * x_2_3_4(i2, v2)(v3, v4);
            for (int i3 = 0; i3 < n; i3++) {
              term2 += y_xx(j)(i2)(i3) * x_2(i2, v2) * x_3_4(i3)(v3, v4);
              term3 += y_xx(j)(i2)(i3) * x_3(i3, v3) * x_2_4(i2)(v2, v4);
              for (int i4 = 0; i4 < n; i4++) {
                term1 += y_xxx(j)(i2)(i3)(i4) * x_2(i2, v2) * x_3(i3, v3) * x_4(i4, v4);
              }
            }
            for (int i4 = 0; i4 < n; i4++) {
              term4 += y_xx(j)(i2)(i4) * x_2_3(i2)(v2, v3) * x_4(i4, v4);
            }
          }
          y_1_2_3(j, v2)(v3, v4) = term1 + term2 + term3 + term4 + term5;
        }
      }
    }
  }

  // Y^{(1, 2, 3, 4)}
  // = F'''' * X^{(1)} * X^{(2)} * X^{(3)} * X^{(4)}
  // + F''' * X^{(1)} * X^{(2)} * X^{(3, 4)}
  // + F''' * X^{(1)} * X^{(3)} * X^{(2, 4)}
  // + F''' * X^{(1)} * X^{(4)} * X^{(2, 3)}
  // + F''' * X^{(2)} * X^{(3)} * X^{(1, 4)}
  // + F''' * X^{(2)} * X^{(4)} * X^{(1, 3)}
  // + F''' * X^{(3)} * X^{(4)} * X^{(1, 2)}
  // + F'' * X^{(1)} * X^{(2, 3, 4)}
  // + F'' * X^{(1, 2)} * X^{(3, 4)}
  // + F'' * X^{(1, 4)} * X^{(2, 3)}
  // + F'' * X^{(1, 2, 3)} * X^{(4)}
  // + F'' * X^{(2)} * X^{(1, 3, 4)}
  // + F'' * X^{(2, 4)} * X^{(1, 3)}
  // + F'' * X^{(3)} * X^{(1, 2, 4)}
  // + F' * X^{(1, 2, 3, 4)}
  for (int j = 0; j < m; j++) {
    for (int v1 = 0; v1 < V1; v1++) {
      for (int v2 = 0; v2 < V2; v2++) {
        for (int v3 = 0; v3 < V3; v3++) {
          for (int v4 = 0; v4 < V4; v4++) {
            T term1, term2, term3, term4, term5, term6, term7, term8, term9,
            term10, term11, term12, term13, term14, term15 = 0;
            for (int i1 = 0; i1 < n; i1++) {
              term15 += y_x(j)(i1) * x_1_2_3_4(i1)(v1, v2)(v3, v4);
              for (int i2 = 0; i2 < n; i2++) {
                term8 += y_xx(j)(i1)(i2) * x_1(i1 ,v1) * x_2_3_4(i2, v2)(v3, v4);
                term10 += y_xx(j)(i1)(i2) * x_1_4(i1)(v1, v4) * x_2_3(i2)(v2, v3);
                term12 += y_xx(j)(i1)(i2) * x_2(i2, v2) * x_1_3_4(i1, v1)(v3, v4);
                term13 += y_xx(j)(i1)(i2) * x_2_4(i2)(v2, v4) * x_1_3(i1)(v1, v3);
                for (int i3 = 0; i3 < n; i3++) {
                  term2 += y_xxx(j)(i1)(i2)(i3) * x_1(i1, v1) * x_2(i1, v2) * x_3_4(i3)(v3, v4);
                  term3 += y_xxx(j)(i1)(i2)(i3) * x_1(i1, v1) * x_3(i3, v3) * x_2_4(i2)(v2, v4);
                  term5 += y_xxx(j)(i1)(i2)(i3) * x_2(i2, v2) * x_3(i3, v3) * x_1_4(i1)(v1, v4);
                  for (int i4 = 0; i4 < n; i4++) {
                    term1 += y_xxxx(j)(i1)(i2)(i3)(i4) * x_1(i1, v1) * x_2(i2, v2) * x_3(i3, v3) * x_4(i4, v4);
                  }
                }
                for (int i4 = 0; i4 < n; i4++) {
                  term4 += y_xxx(j)(i1)(i2)(i4) * x_1(i1, v1)  * x_4(i4, v4) * x_2_3(i1)(v1, v2);
                  term6 += y_xxx(j)(i1)(i2)(i4) * x_2(i2, v2)  * x_4(i4, v4) * x_1_3(i1)(v1, v3);
                }
              }
              for (int i3 = 0; i3 < n; i3++) {
                term9 += y_xx(j)(i1)(i3) * x_1_2(i1)(v1, v2) * x_3_4(i3)(v3, v4);
                term14 += y_xx(j)(i1)(i3) * x_3(i3, v3) * x_1_2_4(i1, v1)(v2, v4);
                for (int i4 = 0; i4 < n; i4++) {
                  term7 += y_xxx(j)(i1)(i3)(i4) * x_3(i3, v3) * x_4(i4, v4) * x_1_2(i1)(v1, v2);
                }
              }
              for (int i4 = 0; i4 < n; i4++) {
                term11 += y_xx(j)(i1)(i4) * x_1_2_3(i1, v1)(v2, v3) * x_4(i4, v4);
              }
            }
            y_1_2_3_4(j)(v1, v2)(v3, v4) = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + term11 + term12 + term13 + term14 + term15;
          }
        }
      }
    }
  }
}


template<typename T, int V1, int V2, int V3, int V4>
bool Validate_vtvtvtvt(std::ofstream& out) {
  T tol = std::pow(std::numeric_limits<T>::epsilon(), 1.0 / 8.0);

  Y_X_t<T> y_x;
  Y_XX_t<T> y_xx;
  Y_XXX_t<T> y_xxx;
  Y_XXXX_t<T> y_xxxx;

  // Outputs
  Y_t<T> AD_y_values;
  Y_t<T> Formula_y_values;
  Eigen::Matrix<T, m, V1> AD_y_1;
  Eigen::Matrix<T, m, V1> Formula_y_1;
  Eigen::Matrix<T, m, V2> AD_y_2; 
  Eigen::Matrix<T, m, V2> Formula_y_2; 
  Eigen::Matrix<T, m, V3> AD_y_3; 
  Eigen::Matrix<T, m, V3> Formula_y_3; 
  Eigen::Matrix<T, m, V4> AD_y_4;
  Eigen::Matrix<T, m, V4> Formula_y_4;
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m> AD_y_1_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, m> Formula_y_1_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V3>, m> AD_y_1_3; 
  Eigen::Vector<Eigen::Matrix<T, V1, V3>, m> Formula_y_1_3; 
  Eigen::Vector<Eigen::Matrix<T, V1, V4>, m> AD_y_1_4;
  Eigen::Vector<Eigen::Matrix<T, V1, V4>, m> Formula_y_1_4;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m> AD_y_2_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, m> Formula_y_2_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V4>, m> AD_y_2_4;
  Eigen::Vector<Eigen::Matrix<T, V2, V4>, m> Formula_y_2_4;
  Eigen::Vector<Eigen::Matrix<T, V3, V4>, m> AD_y_3_4;
  Eigen::Vector<Eigen::Matrix<T, V3, V4>, m> Formula_y_3_4;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1> AD_y_1_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, m, V1> Formula_y_1_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V4>, m, V1> AD_y_1_2_4;
  Eigen::Matrix<Eigen::Matrix<T, V2, V4>, m, V1> Formula_y_1_2_4;
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V1> AD_y_1_3_4;
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V1> Formula_y_1_3_4;
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V2> AD_y_2_3_4;
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, m, V2> Formula_y_2_3_4;
  Eigen::Vector<Eigen::Matrix<Eigen::Matrix<T, V3, V4>, V1, V2>, m> AD_y_1_2_3_4;
  Eigen::Vector<Eigen::Matrix<Eigen::Matrix<T, V3, V4>, V1, V2>, m> Formula_y_1_2_3_4;

  // Inputs
  X_t<T> x_values;
  Eigen::Matrix<T, n, V1> x_1;
  Eigen::Matrix<T, n, V2> x_2;
  Eigen::Matrix<T, n, V3> x_3;
  Eigen::Matrix<T, n, V4> x_4;
  Eigen::Vector<Eigen::Matrix<T, V1, V2>, n> x_1_2;
  Eigen::Vector<Eigen::Matrix<T, V1, V3>, n> x_1_3;
  Eigen::Vector<Eigen::Matrix<T, V1, V4>, n> x_1_4; 
  Eigen::Vector<Eigen::Matrix<T, V2, V3>, n> x_2_3;
  Eigen::Vector<Eigen::Matrix<T, V2, V4>, n> x_2_4; 
  Eigen::Vector<Eigen::Matrix<T, V3, V4>, n> x_3_4;
  Eigen::Matrix<Eigen::Matrix<T, V2, V3>, n, V1> x_1_2_3;
  Eigen::Matrix<Eigen::Matrix<T, V2, V4>, n, V1> x_1_2_4;
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, n, V1> x_1_3_4;
  Eigen::Matrix<Eigen::Matrix<T, V3, V4>, n, V2> x_2_3_4;
  Eigen::Vector<Eigen::Matrix<Eigen::Matrix<T, V3, V4>, V1, V2>, n> x_1_2_3_4;
  for (int i = 0; i < n; i++) {
    x_1_2(i) = Eigen::Matrix<T, V1, V2>::Random().cwiseAbs();
    x_1_3(i) = Eigen::Matrix<T, V1, V3>::Random().cwiseAbs();
    x_1_4(i) = Eigen::Matrix<T, V1, V4>::Random().cwiseAbs();
    x_2_3(i) = Eigen::Matrix<T, V2, V3>::Random().cwiseAbs();
    x_2_4(i) = Eigen::Matrix<T, V2, V4>::Random().cwiseAbs();
    x_3_4(i) = Eigen::Matrix<T, V3, V4>::Random().cwiseAbs();
    for (int v1 = 0; v1 < V1; v1++) {
      x_1_2_3(i, v1) = Eigen::Matrix<T, V2, V3>::Random().cwiseAbs();
      x_1_2_4(i, v1) = Eigen::Matrix<T, V2, V4>::Random().cwiseAbs();
      x_1_3_4(i, v1) = Eigen::Matrix<T, V3, V4>::Random().cwiseAbs();
      for (int v2 = 0; v2 < V2; v2++) {
        x_1_2_3_4(i)(v1, v2) = Eigen::Matrix<T, V3, V4>::Random().cwiseAbs();
      }
    }
    for (int v2 = 0; v2 < V2; v2++) {
      x_2_3_4(i, v2) = Eigen::Matrix<T, V3, V4>::Random().cwiseAbs();
    }
  }

  out << "\n=== Testing Fourth Derivative (Tangent over Tangent over Tangent over Tangent Mode) ===\n";

  // Show the seeds
  out << "Seed for x: \n" << x_values << "\n\n";
  out << "Seed for X^{(1)}: \n" << x_1 << "\n\n";
  out << "Seed for X^{(2)}: \n" << x_2 << "\n\n";
  out << "Seed for X^{(3)}: \n" << x_3 << "\n\n";
  out << "Seed for X^{(4)}: \n" << x_4 << "\n\n";
  out << "Seed for X^{(1, 2)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 3)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 4)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(2, 3)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(2, 4)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(3, 4)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 2, 3)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 2, 4)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 3, 4)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(2, 3, 4)}: (nested matrices not printed)\n\n";
  out << "Seed for X^{(1, 2, 3, 4)}: (nested matrices not printed)\n\n";

  // Populate y_x, y_xx and y_xxx
  vastsa_F_xxx(x_values, y_x, y_xx, y_xxx);

  // Run the AD version
  AD_F_xxxx(
    x_values,
    x_1, x_2, x_3, x_4, 
    x_1_2, x_1_3, x_1_4, x_2_3, x_2_4, x_3_4, 
    x_1_2_3, x_1_2_4, x_1_3_4, x_2_3_4,
    x_1_2_3_4,
    AD_y_values,
    AD_y_1, AD_y_2, AD_y_3, AD_y_4,
    AD_y_1_2, AD_y_1_3, AD_y_1_4, AD_y_2_3, AD_y_2_4, AD_y_3_4,
    AD_y_1_2_3, AD_y_1_2_4, AD_y_1_3_4, AD_y_2_3_4,
    AD_y_1_2_3_4
  );

  // Run the Formula version
  Formula_F_xxxx(
    y_x, y_xx, y_xxx, y_xxxx,
    x_values,
    x_1, x_2, x_3, x_4, 
    x_1_2, x_1_3, x_1_4, x_2_3, x_2_4, x_3_4, 
    x_1_2_3, x_1_2_4, x_1_3_4, x_2_3_4,
    x_1_2_3_4,
    Formula_y_values,
    Formula_y_1, Formula_y_2, Formula_y_3, Formula_y_4,
    Formula_y_1_2, Formula_y_1_3, Formula_y_1_4, Formula_y_2_3, Formula_y_2_4, Formula_y_3_4,
    Formula_y_1_2_3, Formula_y_1_2_4, Formula_y_1_3_4, Formula_y_2_3_4,
    Formula_y_1_2_3_4
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

    for (int v1 = 0; v1 < V1; v1++) {
      // Compare y_1
      diff = std::abs(AD_y_1(j, v1) - Formula_y_1(j, v1));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y^{(1)} Failed\n";
        out << "Validation failed at index " << j << "," << v1 << " Diff: " << diff << "\n";
        return false;
      }

      for (int v2 = 0; v2 < V2; v2++) {
        // Compare y_1_2
        diff = std::abs(AD_y_1_2(j)(v1, v2) - Formula_y_1_2(j)(v1, v2));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y^{(1, 2)} Failed\n";
          out << "Validation failed at index " << j << "," << v1 << "," << v2 << " Diff: " << diff << "\n";
          return false;
        }

        for (int v3 = 0; v3 < V3; v3++) {
          // Compare y_1_2_3
          diff = std::abs(AD_y_1_2_3(j, v1)(v2, v3) - Formula_y_1_2_3(j, v1)(v2, v3));
          maxDiff = std::max(maxDiff, diff);
          if (diff > tol) {
            out << "Validation for Y^{(1, 2, 3)} Failed\n";
            out << "Validation failed at index " << j << "," << v1 << "," << v2 << "," << v3 << " Diff: " << diff << "\n";
            return false;
          }

          for (int v4 = 0; v4 < V4; v4++) {
            // Compare y_1_2_3_4
            diff = std::abs(AD_y_1_2_3_4(j)(v1, v2)(v3, v4) - Formula_y_1_2_3_4(j)(v1, v2)(v3, v4));
            maxDiff = std::max(maxDiff, diff);
            if (diff > tol) {
              out << "Validation for Y^{(1, 2, 3, 4)} Failed\n";
              out << "Validation failed at index " << j << "," << v1 << "," << v2 << "," << v3 << "," << v4 << " Diff: " << diff << "\n";
              return false;
            }
          }
        }
      }

      for (int v3 = 0; v3 < V3; v3++) {
        // Compare y_1_3
        diff = std::abs(AD_y_1_3(j)(v1, v3) - Formula_y_1_3(j)(v1, v3));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y^{(1, 3)} Failed\n";
          out << "Validation failed at index " << j << "," << v1 << "," << v3 << " Diff: " << diff << "\n";
          return false;
        }

        for (int v4 = 0; v4 < V4; v4++) {
          // Compare y_1_3_4
          diff = std::abs(AD_y_1_3_4(j, v1)(v3, v4) - Formula_y_1_3_4(j, v1)(v3, v4));
          maxDiff = std::max(maxDiff, diff);
          if (diff > tol) {
            out << "Validation for Y^{(1, 3, 4)} Failed\n";
            out << "Validation failed at index " << j << "," << v1 << "," << v3 << "," << v4 << " Diff: " << diff << "\n";
            return false;
          }
        }
      }

      for (int v4 = 0; v4 < V4; v4++) {
        // Compare y_1_4
        diff = std::abs(AD_y_1_4(j)(v1, v4) - Formula_y_1_4(j)(v1, v4));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y^{(1, 4)} Failed\n";
          out << "Validation failed at index " << j << "," << v1 << "," << v4 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }

    for (int v2 = 0; v2 < V2; v2++) {
      // Compare y_2
      diff = std::abs(AD_y_2(j, v2) - Formula_y_2(j, v2));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y^{(2)} Failed\n";
        out << "Validation failed at index " << j << "," << v2 << " Diff: " << diff << "\n";
        return false;
      }

      for (int v3 = 0; v3 < V3; v3++) {
        // Compare y_2_3
        diff = std::abs(AD_y_2_3(j)(v2, v3) - Formula_y_2_3(j)(v2, v3));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y^{(2, 3)} Failed\n";
          out << "Validation failed at index " << j << "," << v2 << "," << v3 << " Diff: " << diff << "\n";
          return false;
        }

        for (int v4 = 0; v4 < V4; v4++) {
          // Compare y_2_3_4
          diff = std::abs(AD_y_2_3_4(j, v2)(v3, v4) - Formula_y_2_3_4(j, v2)(v3, v4));
          maxDiff = std::max(maxDiff, diff);
          if (diff > tol) {
            out << "Validation for Y^{(2, 3, 4)} Failed\n";
            out << "Validation failed at index " << j << "," << v2 << "," << v3 << "," << v4 << " Diff: " << diff << "\n";
            return false;
          }
        }
      }

      for (int v4 = 0; v4 < V4; v4++) {
        // Compare y_2_4
        diff = std::abs(AD_y_2_4(j)(v2, v4) - Formula_y_2_4(j)(v2, v4));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y^{(2, 4)} Failed\n";
          out << "Validation failed at index " << j << "," << v2 << "," << v4 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }

    for (int v3 = 0; v3 < V3; v3++) {
      // Compare y_3
      diff = std::abs(AD_y_3(j, v3) - Formula_y_3(j, v3));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y^{(3)} Failed\n";
        out << "Validation failed at index " << j << "," << v3 << " Diff: " << diff << "\n";
        return false;
      }

      for (int v4 = 0; v4 < V4; v4++) {
        // Compare y_3_4
        diff = std::abs(AD_y_3_4(j)(v3, v4) - Formula_y_3_4(j)(v3, v4));
        maxDiff = std::max(maxDiff, diff);
        if (diff > tol) {
          out << "Validation for Y^{(3, 4)} Failed\n";
          out << "Validation failed at index " << j << "," << v3 << "," << v4 << " Diff: " << diff << "\n";
          return false;
        }
      }
    }

    for (int v4 = 0; v4 < V4; v4++) {
      // Compare y_4
      diff = std::abs(AD_y_4(j, v4) - Formula_y_4(j, v4));
      maxDiff = std::max(maxDiff, diff);
      if (diff > tol) {
        out << "Validation for Y^{(4)} Failed\n";
        out << "Validation failed at index " << j << "," << v4 << " Diff: " << diff << "\n";
        return false;
      }
    }
  }
  out << "Maximum difference:\n" << maxDiff << "\n";

  return true;
}