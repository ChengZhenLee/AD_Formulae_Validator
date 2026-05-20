#include "formulaDriver.h"
#include <unordered_set>


template<typename T>
void runFormulaDriver(std::deque<Param<T>> parameters) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string sequence = cm.getSequence();

    // split parameters into input and outputs
    std::deque<Param<T>> inputs = {};
    std::deque<Param<T>> outputs = {};

    for (Param<T>& p : parameters) {
        if (p.role == ParamRole::Input) {
            inputs.push_back(p);
        } else if (p.role == ParamRole::Output) {
            outputs.push_back(p);
        }
    }

    // calculate the required Derivative functions (F_1, F_1_2, ...)
    std::deque<Param<T>> derivatives = getDerivatives<T>(sequence);

    std::deque<Equation<T>> equations = {};

    formulaDriver(sequence, equations, derivatives, inputs, outputs);
}


template<typename T>
void formulaDriver(
    std::string sequence, 
    std::deque<Equation<T>>& equations, 
    const std::deque<Param<T>>& derivatives, 
    const std::deque<Param<T>>& inputs, 
    std::deque<Param<T>>& outputs) 
{
    if (sequence.length() == 0) {
        Param<T>& X = findParamByName("X", inputs);
        Param<T>& Y = findParamByName("Y", outputs);
        // Dummy F variable
        Param<T> F = Param<T>({}, {}, "F", ParamRole::Input); 

        // Calculate primal Y via the primal function
        f(X.data, Y.data);

        // Initialize equations with Y = F * X (using a dummy F value)
        Equation<T> newEquation = Equation<T>(Y);
        Monomial<T> newMonomial = Monomial<T>({Y, F});
        newEquation.rightSide.push_back(newMonomial);

        equations.push_back(newEquation);
    };

    for (size_t order = 0; order < sequence.length(); order++) {
        if (sequence[order] == 't') {
            formulaDriver(sequence.substr(1, sequence.length() - 1), equations, derivatives, inputs, outputs);
            tangentMode(order, equations, inputs, outputs, derivatives);
        } else if (sequence[order] == 'a') {
            formulaDriver(sequence.substr(1, sequence.length() - 1), equations, derivatives, inputs, outputs);
            adjointMode(order, equations, inputs, outputs, derivatives);
        }
    }
}


// Note: inputs also include the Derivative functions F_1, F_1_2 etc.
template<typename T>
void tangentMode(
    size_t order, 
    std::deque<Equation<T>>& equations, 
    const std::deque<Param<T>>& inputs, 
    std::deque<Param<T>>& outputs,
    std::deque<Param<T>>& derivatives)
{
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
                    Param<T> derivedPrimal = findParamByName(std::format("X_{}", order), derivatives);
                    productChain.push_back(derivedPrimal.tensor);
                    paramChain.push_back(derivedPrimal);
                }

                // Compute product and store in new Monomial
                equationSum += Tensor<T>::product(productChain);
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
void adjointMode(
    size_t order, 
    std::deque<Equation<T>>& equations, 
    const std::deque<Param<T>>& inputs,
    std::deque<Param<T>>& outputs,
    std::deque<Param<T>>& derivatives
) 
{
    std::deque<Equation<T>> newEquations = {};
    std::deque<std::string> newOutputNames = {};
    std::unordered_set<std::string> seen = {};

    // Find the names of the new outputs
    for (const Equation<T>& e : equations) {
        for (const Monomial<T>& m : e.rightSide) {
            for (const Param<T>& p : m.parameters) {
                if (p.name[0] != 'F' && seen.find(p.name) == seen.end()) {
                    seen.insert(p.name);
                    newOutputNames.push_back(std::format("{}_{}", p.name, order));
                }
            }
        }
    }

    // Iterate over all the new outputs
    for (std::string& targetName : newOutputNames) {
        Param<T>& targetParam = findParamByName(targetName, outputs);
        Equation<T> newEquation = Equation(targetParam);
        Tensor<T> equationSum(targetParam.tensor.shape);

        // If the name is X_{order}
        if (targetName == std::format("X_{}", order)) {
            // Look for all equations with monomials with "F", 
            // and multiply the derivative of it with the derivative of the output
            for (Equation<T>&e : equations) {
                Param<T> leftSide = e.leftSide;

                for (Monomial<T>& m : e.rightSide) {
                    std::deque<Tensor<T>> productChain = {};
                    std::deque<Param<T>> paramChain = {};
                    std::string derivedName;
                    Param<T> derivedParam;

                    for (Param<T>& p : m.parameters) {
                        
                        if (p.name[0] == 'F') {
                            derivedName = std::format("{}_{}", p.name, order);
                            derivedParam = findParamByName(derivedName, derivatives);

                            productChain.push_back(derivedParam.tensor);
                            paramChain.push_back(derivedParam);
                        }

                        else {
                            productChain.push_back(p.tensor);
                            paramChain.push_back(p);
                        }
                    }

                    // Multiply by the derivative of the previous output
                    derivedName = std::format("{}_{}", leftSide.name, order);
                    derivedParam = findParamByName(derivedName, inputs);
                    productChain.push_back(derivedParam.tensor);
                    paramChain.push_back(derivedParam);

                    // Sum it into the new equation
                    equationSum += Tensor<T>::product(productChain);
                    newEquation.rightSide.push_back(Monomial<T>(paramChain));
                }

            }
        }

        // Name is not "X_{order}"
        else {
            for (Equation<T>&e : equations) {
                Param<T> leftSide = e.leftSide;

                // Only consider monomials where we have the original name
                // e.g. for targetName == Y_1_2, we are only looking for Y_1 in the monomial
                for (Monomial<T>& m : e.rightSide) {
                    std::deque<Tensor<T>> productChain = {};
                    std::deque<Param<T>> paramChain = {};
                    std::string derivedName;
                    Param<T> derivedParam;
                    bool found = false;
                    size_t pos = targetName.find_last_of('_');
                    std::string originalName = (pos == std::string::npos)
                        ? targetName
                        : targetName.substr(0, pos);

                    for (Param<T>& p : m.parameters) {
                        if (p.name == originalName) {
                            found = true;
                        }
                    }

                    // Skip this monomial
                    if (!found) continue;

                    // If the original name is found, swap it with the derivative of the orignal output
                    for (Param<T>& p : m.parameters) {
                        if (p.name != originalName) {
                            productChain.push_back(p.tensor);
                            paramChain.push_back(p);
                        } else {
                            derivedName = std::format("{}_{}", leftSide.name, order);
                            derivedParam = findParamByName(derivedName, inputs);
                            productChain.push_back(derivedParam.tensor);
                            paramChain.push_back(derivedParam);
                        }
                    }

                    // Sum it into the new equation
                    equationSum += Tensor<T>::product(productChain);
                    newEquation.rightSide.push_back(Monomial<T>(paramChain));
                }
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