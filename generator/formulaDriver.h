#ifndef FORMULADRIVER_H
#define FORMULADRIVER_H

#include <format>
#include "structures.h"
#include "utils.h"
#include "configManager.h"
#include "user_function.h"


template<typename T>
void runFormulaDriver(std::deque<Param<T>> parameters);

template<typename T>
void formulaDriver(
    std::string sequence, 
    std::deque<Equation<T>>& equations, 
    const std::deque<Param<T>>& derivates, 
    const std::deque<Param<T>>& inputs, 
    std::deque<Param<T>>& outputs
);

template<typename T>
void tangentMode(
    size_t order, 
    std::deque<Equation<T>>& equations, 
    const std::deque<Param<T>>& inputs, 
    std::deque<Param<T>>& outputs,
    std::deque<Param<T>>& derivatives
);

template<typename T>
void adjointMode(
    size_t order, 
    std::deque<Equation<T>>& equations, 
    const std::deque<Param<T>>& inputs,
    std::deque<Param<T>>& outputs,
    std::deque<Param<T>>& derivatives
);

#endif