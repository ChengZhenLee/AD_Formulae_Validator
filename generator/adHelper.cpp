#include "adHelper.h"


template <typename T>
void AD_F_2(X_t<T_t<T_t<double,3>,3>>& x, Y_t<T_t<T_t<double,3>,3>>& y, std::vector<Param<T>> parameters) {
	seedADForOrder(x, y, parameters, 2, "tt");
	AD_F_1(x, y, parameters);
	extractADForOrder(x, y, parameters, 2, "tt");
}

template <typename T>
void AD_F_1(X_t<T_t<T_t<double,3>,3>>& x, Y_t<T_t<T_t<double,3>,3>>& y, std::vector<Param<T>> parameters) {
	seedADForOrder(x, y, parameters, 1, "tt");
	f(x, y);
	extractADForOrder(x, y, parameters, 1, "tt");
	extractPrimal(y, parameters);
}

template<typename T> 
std::vector<Param<T>> runHelperDrivers(std::vector<Param<T>> &parameters) {
	ConfigManager& cm = ConfigManager::getInstance();
	X_t<T_t<T_t<double,3>,3>> x;
	Y_t<T_t<T_t<double,3>,3>> y;
	x.resize(cm.getXShape());
	y.resize(cm.getYShape());
	AD_F_2<double>(x, y, parameters);
	return parameters;
}

int main(int argc, char** argv) {
	ConfigManager::getInstance().load("generator/configs.txt");

	std::vector<Param<double>> parameters = readParameters<double>("generator/derivatives_seeds.bin");

	auto results = runHelperDrivers<double>(parameters);
	writeParameters<double>(results, "generator/derivatives.bin");

	return 0;
}
