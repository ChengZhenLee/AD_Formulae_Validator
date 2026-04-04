#include "ad.h"


template<typename T, int size>
using T_t=ad::tangent_t<T, size>;

template<typename T, int size>
using A_t=ad::adjoint_t<T, size>;
