#ifndef USER_FUNCTION_H
#define USER_FUNCTION_H

#include <vector>
#include "structures.h"

// Replace the body of f() below with your own formula to validate it.
//
// Keep the signature exactly as-is: it's instantiated for both plain
// doubles and nested AD types, so it must stay generic over T.
//   - x.size() must equal the "x" value in generator/configs.txt.
//   - y.size() must equal the "y" value in generator/configs.txt.
//   - Only use operations with a generic AD-aware overload: the arithmetic
//     operators (+ - * / and unary -), plus sin, cos, tan, exp, sqrt,
//     pow(value, int n), log10, fabs. Anything from <cmath>/std:: (e.g.
//     std::sin) will not compile against the AD types and should be avoided.
//
// After editing, run the validator with the mode you want to check, e.g.:
//   ADValidator.exe --sequence=at
template<typename T>
void f(X_t<T>&x, Y_t<T>&y) {
    T v=x[1]*x[2];
    v=tan(v);
    T w=x[0]-v;
    y[0]=v/w;
    y[1]=y[0]*x[0];
};


#endif