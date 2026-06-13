#include <cstddef>
#include <deque>
#include <unordered_set>
#include "structures.h"
#include "configManager.h"

template<typename T>
void runFormulaDriver(std::deque<Param<T>> parameters) {
    // === TEMPORARY DEBUG DUMP ===
    std::cout << "\n========================================\n";
    std::cout << "   DUMPING ALL PARAMETERS FROM FILE     \n";
    std::cout << "========================================\n";
    for (size_t i = 0; i < parameters.size(); ++i) {
        std::cout << "Index [" << i << "] "
                  << "Name: " << parameters[i].name << " | "
                  << "Role: " << (parameters[i].role == ParamRole::Input ? "Input" : "Output") << " | "
                  << "Shape: [";
        for (size_t s : parameters[i].tensor.shape) std::cout << s << " ";
        std::cout << "]\n";
    }
    std::cout << "========================================\n\n";
    // === END DEBUG DUMP ===

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

    // calculate the required Derivative functions (F_1, F_2, ...)
    // Extract primal X from inputs
    const Param<T>& Xparam = findParamByName<T>("X", inputs);
    std::vector<T> x0;
    x0.reserve(Xparam.tensor.data.size());
    for (const auto &v : Xparam.tensor.data) x0.push_back(v);

    std::vector<Param<T>> derivedVec = getDerivatives<T>(sequence, x0);
    std::deque<Param<T>> derivatives(derivedVec.begin(), derivedVec.end());
    // === TEMPORARY DEBUG DUMP ===
    std::cout << "\n========================================\n";
    std::cout << "   DUMPING ALL PARAMETERS FROM FILE     \n";
    std::cout << "========================================\n";
    for (size_t i = 0; i < derivatives.size(); ++i) {
        std::cout << "Index [" << i << "] "
                  << "Name: " << derivatives[i].name << " | "
                  << "Role: " << (derivatives[i].role == ParamRole::Input ? "Input" : "Output") << " | "
                  << "Shape: [";
        for (size_t s : derivatives[i].tensor.shape) std::cout << s << " ";
        std::cout << "]\n";
    }
    std::cout << "========================================\n\n";
    // === END DEBUG DUMP ===

    std::deque<Equation<T>> equations = {};

    formulaDriver(sequence, equations, derivatives, inputs, outputs);
}


template<typename T>
void formulaDriver(
    std::string sequence, 
    std::deque<Equation<T>>& equations, 
    std::deque<Param<T>>& derivatives, 
    std::deque<Param<T>>& inputs, 
    std::deque<Param<T>>& outputs) 
{

    // Initialize base primal equations and execute initial evaluation
    Param<T>& X = findParamByName("X", inputs);
    Param<T>& Y = findParamByName("Y", outputs);
    Param<T> F = Param<T>({}, {}, "F", ParamRole::Input); 
    inputs.push_back(F);

    // Calculate primal Y via the primal function
    f(X.tensor.data, Y.tensor.data);

    // Initialize equations with Y = F * X (using a dummy F value)
    Equation<T> newEquation = Equation<T>(Y);
    Monomial<T> newMonomial = Monomial<T>({F, X});
    newEquation.rightSide.push_back(newMonomial);
    equations.push_back(newEquation);

    // Linearly progress forward over differentiation orders
    for (size_t order = 1; order <= sequence.length(); order++) {
        if (sequence[order - 1] == 't') {
            tangentMode(order, equations, inputs, outputs, derivatives);
        } else if (sequence[order - 1] == 'a') {
            adjointMode(order, equations, inputs, outputs, derivatives);
        }
    }
}

// Note: inputs also include the Derivative functions F_1, F_1_2 etc.
template<typename T>
void tangentMode(
    size_t order, 
    std::deque<Equation<T>>& equations, 
    std::deque<Param<T>>& inputs, 
    std::deque<Param<T>>& outputs,
    std::deque<Param<T>>& derivatives)
{
    std::deque<Equation<T>> newEquations = {};

    for (const Equation<T>& e : equations) {
        // 1. Identify the target for the current order (e.g., Y_1)
        std::string targetName = makeParamName(e.leftSide.name, order);
        Param<T>& targetParam = findParamByName(targetName, outputs);

        Equation<T> newEquation(targetParam);
        Tensor<T> equationSum(targetParam.tensor.shape);

        for (const Monomial<T>& m : e.rightSide) {
            // 2. Perform Product Rule: d(A*B*C) = dA*B*C + A*dB*C + A*B*dC
            for (size_t i = 0; i < m.parameters.size(); ++i) {
                std::deque<Tensor<T>> productChain = {};
                std::deque<Param<T>> paramChain = {};

                // Handle the derivative of the i-th parameter
                std::string derivedName = makeParamName(m.parameters[i].name, order);
                Param<T> derivedParam = findParamByName(derivedName, derivatives);

                bool insertPrimalAfterDerived = false;
                Param<T> derivedPrimal;
                if (derivedName.front() == 'F') {
                    insertPrimalAfterDerived = true;
                    try {
                        derivedPrimal = findParamByName(makeParamName("X", order), derivatives);
                    } catch (...) {
                        derivedPrimal = findParamByName(makeParamName("X", order), inputs);
                    }
                }

                // Build the Monomial: [F_higher, params..., derivedParam]
                for (size_t j = 0; j < m.parameters.size(); ++j) {
                    if (i == j) {
                        productChain.push_back(derivedParam.tensor);
                        paramChain.push_back(derivedParam);
                        if (insertPrimalAfterDerived) {
                            productChain.push_back(derivedPrimal.tensor);
                            paramChain.push_back(derivedPrimal);
                        }
                    } else {
                        productChain.push_back(m.parameters[j].tensor);
                        paramChain.push_back(m.parameters[j]);
                    }
                }

                // Validate product chain order and resulting shape against targetParam
                std::deque<size_t> resultShape = Tensor<T>::productShape(productChain);
                if (resultShape.empty()) {
                    throw std::runtime_error(std::format("Invalid tensor chain order when building monomial in tangentMode: {} tensors; adjacent shapes were not contractable", productChain.size()));
                }
                if (resultShape != targetParam.tensor.shape) {
                    throw std::runtime_error("Tensor chain result shape does not match targetParam shape for tangentMode");
                }
                equationSum = equationSum + Tensor<T>::product(productChain);
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
    std::deque<Param<T>>& inputs,
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
                    newOutputNames.push_back(makeParamName(p.name, order));
                }
            }
        }
    }

    // Iterate over all the new outputs
    for (std::string& targetName : newOutputNames) {
        Param<T>& targetParam = findParamByName(targetName, outputs);
        Equation<T> newEquation = Equation<T>(targetParam);
        Tensor<T> equationSum(targetParam.tensor.shape);

        // If the name matches the expected X dependency identifier
        if (targetName == makeParamName("X", order)) {
            for (Equation<T>& e : equations) {
                Param<T> leftSide = e.leftSide;

                for (Monomial<T>& m : e.rightSide) {
                    std::deque<Tensor<T>> productChain = {};
                    std::deque<Param<T>> paramChain = {};
                    std::string derivedName;
                    Param<T> derivedParam;

                    for (Param<T>& p : m.parameters) {
                        if (p.name[0] == 'F' && p.name.length() > 1) {
                            derivedName = makeParamName(p.name, order);
                            derivedParam = findParamByName(derivedName, derivatives);

                            productChain.push_back(derivedParam.tensor);
                            paramChain.push_back(derivedParam);
                        } else {
                            productChain.push_back(p.tensor);
                            paramChain.push_back(p);
                        }
                    }

                    // Multiply by the derivative of the previous output
                    derivedName = makeParamName(leftSide.name, order);
                    derivedParam = findParamByName(derivedName, inputs);
                    productChain.push_back(derivedParam.tensor);
                    paramChain.push_back(derivedParam);

                    // Validate product chain shape alignment
                    std::deque<size_t> resultShape = Tensor<T>::productShape(productChain);
                    if (resultShape.empty()) {
                        throw std::runtime_error(std::format("Invalid tensor chain order when building X_{} in adjointMode: {} tensors; adjacent shapes were not contractable", order, productChain.size()));
                    }
                    if (resultShape != targetParam.tensor.shape) {
                        throw std::runtime_error(std::format("Tensor chain result shape does not match targetParam shape for X_{}", order));
                    }
                    equationSum = equationSum + Tensor<T>::product(productChain);
                    newEquation.rightSide.push_back(Monomial<T>(paramChain));
                }
            }
        }
        else {
            for (Equation<T>& e : equations) {
                Param<T> leftSide = e.leftSide;

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

                    if (!found) continue;

                    // Swap the original element with the derivative tracking parameter
                    for (Param<T>& p : m.parameters) {
                        if (p.name != originalName) {
                            productChain.push_back(p.tensor);
                            paramChain.push_back(p);
                        } else {
                            derivedName = makeParamName(leftSide.name, order);
                            derivedParam = findParamByName(derivedName, inputs);
                            productChain.push_back(derivedParam.tensor);
                            paramChain.push_back(derivedParam);
                        }
                    }

                    // Validate product chain shape alignments
                    std::deque<size_t> resultShape = Tensor<T>::productShape(productChain);
                    if (resultShape.empty()) {
                        throw std::runtime_error(std::format("Invalid tensor chain order when building {} in adjointMode: {} tensors; adjacent shapes were not contractable", targetName, productChain.size()));
                    }
                    if (resultShape != targetParam.tensor.shape) {
                        throw std::runtime_error(std::format("Tensor chain result shape does not match targetParam shape for {}", targetName));
                    }
                    equationSum = equationSum + Tensor<T>::product(productChain);
                    newEquation.rightSide.push_back(Monomial<T>(paramChain));
                }
            }
        }

        // Update target data output
        targetParam.tensor = equationSum;
        newEquations.push_back(newEquation);
    }

    // Finalize all new calculations back into global tracking deque
    for (auto& ne : newEquations) {
        equations.push_back(ne);
    }
}