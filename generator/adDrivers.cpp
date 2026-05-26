#include "adDrivers.h"
#include "structures.h"


template <typename T>
void AD_F_1(X_t<A_t<T_t<double,3>,2>>& x, Y_t<A_t<T_t<double,3>,2>>& y, std::vector<Param<T>> parameters) {
	seedADForOrder(x, parameters, 1, "ta");
	f(x, y);
	extractADForOrder(y, parameters, 1, "ta");
	extractPrimal(y, parameters);
}

template <typename T>
void AD_F_2(X_t<A_t<T_t<double,3>,2>>& x, Y_t<A_t<T_t<double,3>,2>>& y, std::vector<Param<T>> parameters) {
	seedPrimal(x, parameters);
	for (size_t i = 0; i < 3; i++) {
		x[i].register_input();
	}
	AD_F_1(x, y, parameters);
	seedADForOrder(x, parameters, 2, "ta");
	A_t<T_t<double,3>,2>::tape::interpret();
	extractADForOrder(y, parameters, 2, "ta");
	A_t<T_t<double,3>,2>::tape::reset();
}

template<typename T> 
std::vector<Param<T>> runADDrivers(std::string mode, X_t<T>& x_input) {
	std::vector<Param<T>> parameters; 
	if (mode == "random") {
		parameters = generateRandom(2, x_input); 
	} else if (mode == "identity") { 
		parameters = generateIdentity(2, x_input); 
	}
	X_t<A_t<T_t<double,3>,2>> x;
	Y_t<A_t<T_t<double,3>,2>> y;
	AD_F_2(x, y, parameters);
	return parameters;
}
