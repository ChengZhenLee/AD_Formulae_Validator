#ifndef LIGHTHOUSE_TYPES_H
#define LIGHTHOUSE_TYPES_H

#include "Eigen/Dense"
#include <iostream>

constexpr int n=3, m=2;

// input type
template<typename T>
using X_t=Eigen::Vector<T,n>;

// output type
template<typename T>
using Y_t=Eigen::Vector<T,m>;

// first derivative (Jacobian) type
template<typename T>
using Y_X_t=Y_t<X_t<T>>;

template <typename T>
void print(const Y_X_t<T>& y_x) {
  std::cout << "y_x=\n";
  for (int j=0;j<m;++j) {
    std::cout << y_x(j).transpose() << '\n'; 
  }
}

// second derivative (Hessian) type
template<typename T>
using Y_XX_t=Y_t<X_t<X_t<T>>>;

template <typename T>
void print(const Y_XX_t<T>& y_xx) {
  for (int j=0;j<m;++j) {
    std::cout << "y(" << j << ")_xx=\n";
    for (int i1=0;i1<n;++i1) {
      std::cout << y_xx(j)(i1).transpose() << '\n';
    }
  }
}

// third derivative type
template<typename T>
using Y_XXX_t=Y_t<X_t<X_t<X_t<T>>>>;

template <typename T>
void print(const Y_XXX_t<T>& y_xxx) {
  for (int j=0;j<m;++j) {
    for (int i1=0;i1<n;++i1) {
      std::cout << "y(" << j << ")_x(" << i1 << ")_xx=\n";
      for (int i2=0;i2<n;++i2) {
        std::cout << y_xxx(j)(i1)(i2).transpose() << '\n';
      }
    }
  }
}

// fourth derivative type
template<typename T>
using Y_XXXX_t=Y_t<X_t<X_t<X_t<X_t<T>>>>>;

template <typename T>
void print(const Y_XXXX_t<T>& y_xxxx) {
  for (int j=0;j<m;++j) {
    for (int i1=0;i1<n;++i1) {
      for (int i2=0;i2<n;++i2) {
        std::cout << "y(" << j << ")_x(" << i1 << ")_x(" << i2 << ")_xx=\n";
        for (int i3=0;i3<n;++i3) {
          std::cout << y_xxxx(j)(i1)(i2)(i3).transpose() << '\n';
        }
      }
    }
  }
}

#endif
