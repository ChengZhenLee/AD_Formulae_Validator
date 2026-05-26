#ifndef VALIDATOR_H
#define VALIDATOR_H

#include "structures.h"
#include <algorithm>
#include <limits>


template<typename T>
bool validatorAll(std::vector<Param<T>> parametersA, std::vector<Param<T>> parametersB);

template<typename T>
bool validateParameter(Param<T> &parameterA, Param<T> &parameterB);


#endif