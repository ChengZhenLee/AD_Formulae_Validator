#ifndef ADDRIVER_H
#define ADDRIVER_H

#include "structures.h"

template <typename T>
std::vector<Param<T>> runADDrivers(std::string mode, X_t<T>& x); 

template <typename T>
void AD_F_1(X_t<A_t<T_t<double,3>,2>>& x, Y_t<A_t<T_t<double,3>,2>>& y, std::vector<Param<T>> parameters); 

template <typename T>
void AD_F_2(X_t<A_t<T_t<double,3>,2>>& x, Y_t<A_t<T_t<double,3>,2>>& y, std::vector<Param<T>> parameters); 


#endif
