#include "adDrivers.h"


template <typename T>
void AD_F_2(X_t<T_t<A_t<double,2>,3>>& x, Y_t<T_t<A_t<double,2>,3>>& y, std::vector<Param<T>> parameters) {
	AD_F_1(x, y, parameters);
	A_t<double,2>::tape::init_adjoints();
	seedADForOrder(x, y, parameters, 2, "at");
	A_t<double,2>::tape::interpret();
	extractADForOrder(x, y, parameters, 2, "at");
}

template <typename T>
void AD_F_1(X_t<T_t<A_t<double,2>,3>>& x, Y_t<T_t<A_t<double,2>,3>>& y, std::vector<Param<T>> parameters) {
	seedADForOrder(x, y, parameters, 1, "at");
	f(x, y);
	extractADForOrder(x, y, parameters, 1, "at");
	extractPrimal(y, parameters);
}

template<typename T> 
std::vector<Param<T>> runADDrivers(std::vector<Param<T>> &parameters) {
	ConfigManager& cm = ConfigManager::getInstance();
	X_t<T_t<A_t<double,2>,3>> x;
	Y_t<T_t<A_t<double,2>,3>> y;
	x.resize(cm.getXShape());
	y.resize(cm.getYShape());
	seedPrimal(x, parameters);

	for (size_t i = 0; i < 3; i++) {
		x[i].value().register_input();
	}
	AD_F_2<double>(x, y, parameters);
	// Clean up all recorded tapes before returning
	A_t<double,2>::tape::reset();

	return parameters;
}

int main(int argc, char** argv) {
	ConfigManager::getInstance().load("generator/configs.txt");

	std::vector<Param<double>> parameters = readParameters<double>("generator/parameters.bin");

	auto results = runADDrivers<double>(parameters);
	writeParameters<double>(results, "generator/results.bin");

	return 0;
}
