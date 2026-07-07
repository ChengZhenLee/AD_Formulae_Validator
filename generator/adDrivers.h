#ifndef ADDRIVER_H
#define ADDRIVER_H

#include "structures.h"
#include "utils.h"
#include "user_function.h"
#include "readWrite.h"
#include "configManager.h"

template <typename T>
std::vector<Param<T>> runADDrivers(X_t<T>& x); 

template <typename T>
void AD_F_1(X_t<T_t<A_t<double,2>,3>>& x, Y_t<T_t<A_t<double,2>,3>>& y, std::vector<Param<T>> parameters); 

template <typename T>
void AD_F_2(X_t<T_t<A_t<double,2>,3>>& x, Y_t<T_t<A_t<double,2>,3>>& y, std::vector<Param<T>> parameters); 


#endif
