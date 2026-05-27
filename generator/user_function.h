#ifndef USER_FUNCTION_H
#define USER_FUNCTION_H

#include <vector>
#include "structures.h"

template<typename T>
void f(X_t<T>&x, Y_t<T>&y) {
    T v=x(1)*x(2);
    v=tan(v);
    T w=x(0)-v;
    y(0)=v/w;
    y(1)=y(0)*x(0);
};


#endif