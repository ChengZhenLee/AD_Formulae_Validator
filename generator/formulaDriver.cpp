#include <format>
#include "structures.h"
#include "utils.h"


template<typename T>
void formulaDriver(std::string sequence, std::deque<Equation<T>>& equations, const std::deque<Param<T>>& derivates, const std::deque<Param<T>>& inputs, std::deque<Param<T>>& outputs) {
    if (sequence.length() == 0) return;

    for (size_t order = 0; order < sequence.length(); order++) {
        if (sequence[order] == 't') {
            formulaDriver(sequence.substr(1, sequence.length() - 1), equations, inputs, outputs);
            tangentMode(order, equations, inputs);
        } else if (sequence[order] == 'a') {
            formulaDriver(sequence.substr(1, sequence.length() - 1), equations, inputs, outputs);
            adjointMode(order, equations, inputs);
        }
    }
}


// Note: inputs also include the Derivative functions F_1, F_1_2 etc.
template<typename T>
void tangentMode(size_t order, std::deque<Equation<T>>& equations, const std::deque<Param<T>>& inputs, std::deque<Param<T>>& outputs) {
    std::deque<Equation<T>> newEquations = {};

    for (const Equation<T>& e : equations) {
        // 1. Identify the target for the current order (e.g., Y_2)
        std::string targetName = std::format("{}_{}", e.leftSide.name, order);
        Param<T>& targetParam = findParamByName(targetName, outputs);

        Equation<T> newEquation(targetParam);
        Tensor<T> equationSum(targetParam.tensor.shape);

        for (const Monomial<T>& m : e.rightSide) {
            // 2. Perform Product Rule: d(A*B*C) = dA*B*C + A*dB*C + A*B*dC
            
            for (size_t i = 0; i < m.parameters.size(); ++i) {
                std::deque<Tensor<T>> productChain = {};
                std::deque<Param<T>> paramChain = {};

                // Handle the derivative of the i-th parameter
                // The inputs should also include for example F_1_2 for the second order derivative Hessian
                std::string derivedName = std::format("{}_{}", m.parameters[i].name, order);
                Param<T> derivedParam = findParamByName(derivedName, inputs);

                // Build the Monomial: [F_higher, params..., derivedParam]
                for (size_t j = 0; j < m.parameters.size(); ++j) {
                    if (i == j) {
                        productChain.push_back(derivedParam.tensor);
                        paramChain.push_back(derivedParam);
                    } else {
                        productChain.push_back(m.parameters[j].tensor);
                        paramChain.push_back(m.parameters[j]);
                    }
                }

                // If deriving the Derivative function F_{order}, we need to add the derived primal e.g. X_2 for F_1_2
                if (derivedName.front() == 'F') {
                    Param<T> derivedPrimal = findParamByName(std::format("X_{}", order));
                    productChain.push_back(derivedPrimal.tensor);
                    paramChain.push_back(derivedPrimal);
                }

                // Compute product and store in new Monomial
                equationSum += Tensor<T>::recProduct(productChain);
                newEquation.rightSide.push_back(Monomial<T>(paramChain));
            }
        }
        
        // Update the actual data in the output registry
        targetParam.tensor = equationSum;
        newEquations.push_back(newEquation);
    }

    // Append all newly derived equations to the system
    for (auto& ne : newEquations) {
        equations.push_back(ne);
    }
}


template<typename T>
void adjointMode(size_t order, std::deque<Equation<T>>& equations, const std::deque<Param<T>>& inputs) {

}
