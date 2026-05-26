#include "validator.h"


// Validates the tensors of all parameters
template<typename T>
bool validatorAll(std::vector<Param<T>> parametersA, std::vector<Param<T>> parametersB) {
    for (Param<T> paramA : parametersA) {
        auto it = std::find_if(parametersB.begin(), parametersB.end(), 
        [&](const Param<T>& paramB) {
            return paramB.name == paramA.name;
        });

        if (it != parametersB.end()) {
            T tolerance = std::sqrt(std::numeric_limits<T>::epsilon, std::pow(2, paramA.highestOrder));

            if (!validateParameter(paramA, it)) return false;
        }
    }

    return true;
}


// Helper function to validate the tensor of a single parameter
template<typename T>
bool validateParameter(Param<T> &parameterA, Param<T> &parameterB) {
    if (parameterA.name != parameterB.name)
        return false;

    return Tensor<T>::compareTensors(parameterA.tensor, parameterB.tensor);
}