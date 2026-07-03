#include "adDrivers.h"


template <typename T>
void AD_F_2(X_t<A_t<A_t<double,2>,2>>& x, Y_t<A_t<A_t<double,2>,2>>& y, std::vector<Param<T>> parameters) {
	AD_F_1(x, y, parameters);
	A_t<double,2>::tape::init_adjoints();
	seedADForOrder(x, y, parameters, 2, "aa");
	A_t<double,2>::tape::interpret();
	extractADForOrder(x, y, parameters, 2, "aa");
}

template <typename T>
void AD_F_1(X_t<A_t<A_t<double,2>,2>>& x, Y_t<A_t<A_t<double,2>,2>>& y, std::vector<Param<T>> parameters) {
	f(x, y);
	A_t<A_t<double,2>,2>::tape::init_adjoints();
	seedADForOrder(x, y, parameters, 1, "aa");
	A_t<A_t<double,2>,2>::tape::interpret();
	extractADForOrder(x, y, parameters, 1, "aa");
	extractPrimal(y, parameters);
}

template<typename T> 
std::vector<Param<T>> runADDrivers(std::vector<Param<T>> &parameters) {
	ConfigManager& cm = ConfigManager::getInstance();
	X_t<A_t<A_t<double,2>,2>> x;
	Y_t<A_t<A_t<double,2>,2>> y;
	x.resize(cm.getXShape());
	y.resize(cm.getYShape());
	seedPrimal(x, parameters);

	for (size_t i = 0; i < 3; i++) {
		x[i].register_input();
		x[i].value().register_input();
	}
	AD_F_2<double>(x, y, parameters);
	// Clean up all recorded tapes before returning
	A_t<double,2>::tape::reset();
	A_t<A_t<double,2>,2>::tape::reset();

	return parameters;
}

int main(int argc, char** argv) {
	ConfigManager::getInstance().load("generator/configs.txt");

	std::vector<Param<double>> parameters = readParameters<double>("generator/parameters.bin");

	auto results = runADDrivers<double>(parameters);
	writeParameters<double>(results, "generator/results.bin");

	return 0;
}
