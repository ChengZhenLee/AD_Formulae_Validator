#include <cstddef>
#include <deque>
#include <fstream>
#include <ostream>
#include <unordered_set>
#include "structures.h"
#include "configManager.h"
#include "utils.h"


// Renders a single monomial (a chain of contracted parameters) as "A * B * C".
template<typename T>
std::string monomialToString(const Monomial<T>& m) {
    std::string s;
    for (size_t i = 0; i < m.parameters.size(); ++i) {
        if (i > 0) s += " * ";
        s += m.parameters[i].name;
    }
    return s;
}

// Writes every equation accumulated during the symbolic derivation as
// "LEFTSIDE = monomial + monomial + ...", one per line, so a peer can inspect
// exactly which terms the formula driver derived for their sequence.
template<typename T>
void saveEquations(const std::deque<Equation<T>>& equations, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Unable to open file: " << filename << "\n";
        return;
    }

    for (const Equation<T>& e : equations) {
        outFile << e.leftSide.name << " = ";
        for (size_t i = 0; i < e.rightSide.size(); ++i) {
            if (i > 0) outFile << " + ";
            outFile << monomialToString(e.rightSide[i]);
        }
        outFile << "\n";
    }
}

template<typename T>
std::vector<Param<T>> runFormulaDriver(std::vector<Param<T>> parameters, const std::string& equationsPath = "") {
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

    std::deque<Equation<T>> equations = {};

    formulaDriver(sequence, equations, derivatives, inputs, outputs);

    if (!equationsPath.empty()) {
        saveEquations(equations, equationsPath);
    }

    return parameters;
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
    Param<T> F = Param<T>({}, {}, "F", ParamRole::Input, {"j", "i"}); 

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

template<typename T>
void tangentMode(
    size_t order, 
    std::deque<Equation<T>>& equations, 
    std::deque<Param<T>>& inputs, 
    std::deque<Param<T>>& outputs,
    std::deque<Param<T>>& derivatives)
{
    std::deque<Equation<T>> newEquations = {};

    for (Equation<T>& e : equations) {
        std::string targetName = makeParamName(e.leftSide.name, order);
        Param<T>& targetParam = findParamByName(targetName, outputs);

        Equation<T> newEquation(targetParam);
        Tensor<T> equationSum(targetParam.tensor.shape);

        for (const Monomial<T>& m : e.rightSide) {
            for (size_t i = 0; i < m.parameters.size(); ++i) {
                if (m.parameters[i].name == "X") continue;

                std::deque<Param<T>> collectedParams = {};

                // Process the active derivative element for the current slot
                Param<T> derivedParam;
                std::string derivedName = makeParamName(m.parameters[i].name, order);
                if (derivedName.front() == 'F') {
                    derivedParam = findParamByName(derivedName, derivatives);
                } else {
                    derivedParam = findParamByName(derivedName, inputs);
                }

                bool insertPrimalAfterDerived = false;
                Param<T> derivedPrimal;
                if (derivedName.front() == 'F' && derivedName.length() > 1) {
                    insertPrimalAfterDerived = true;
                    derivedPrimal = findParamByName(makeParamName("X", order), inputs);
                }

                // Build the execution queue
                for (size_t j = 0; j < m.parameters.size(); ++j) {
                    if (i == j) {
                        collectedParams.push_back(derivedParam);
                        if (insertPrimalAfterDerived) {
                            collectedParams.push_back(derivedPrimal);
                        }
                    } else {
                        if (m.parameters[j].name == "X") continue;
                        collectedParams.push_back(m.parameters[j]);
                    }
                }

                if (collectedParams.empty()) continue;

                std::deque<Param<T>> paramChainCopy = collectedParams;

                Tensor<T> finalMonomialTensor = evaluateMonomialSequence(collectedParams, targetParam.indexNames);

                equationSum = equationSum + finalMonomialTensor;

                newEquation.rightSide.push_back(Monomial<T>(paramChainCopy));
            }
        }
        
        targetParam.tensor = equationSum;
        newEquations.push_back(newEquation);
    }

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

    for (std::string& targetName : newOutputNames) {
        Param<T>& targetParam = findParamByName(targetName, outputs);
        Equation<T> newEquation = Equation<T>(targetParam);
        Tensor<T> equationSum(targetParam.tensor.shape);

        // CASE 1: X tracking path dependencies
        if (targetName == makeParamName("X", order)) {
            for (Equation<T>& e : equations) {
                Param<T> leftSide = e.leftSide;

                for (Monomial<T>& m : e.rightSide) {
                    std::deque<Param<T>> collectedParams = {};
                    for (const Param<T>& p : m.parameters) {
                        // Standard filter for primal X or its tangent/higher-order mutations
                        if (p.name == "X") {
                            continue;
                        }
                        
                        if (p.name[0] == 'F') {
                            std::string lookupName = makeParamName(p.name, p.highestOrder + 1);
                            try {
                                Param<T> found = findParamByName(lookupName, derivatives);
                                collectedParams.push_back(found);
                            } catch (const std::exception&) {
                            }
                        } else {
                            try {
                                Param<T> found = findParamByName(p.name, inputs);
                                collectedParams.push_back(found);
                            } catch (const std::exception&) {
                            }
                        }
                    }

                    // Append the incoming adjoint seed node tracker
                    std::string seedName = makeParamName(leftSide.name, order);
                    try {
                        Param<T> foundSeed = findParamByName(seedName, inputs);
                        collectedParams.push_back(foundSeed);
                    } catch (const std::exception&) {
                    }

                    if (collectedParams.empty()) {
                        continue;
                    }

                    std::deque<Param<T>> paramChainCopy = collectedParams;

                    Tensor<T> finalMonomialTensor = evaluateMonomialSequence(collectedParams, targetParam.indexNames);

                    equationSum = equationSum + finalMonomialTensor;

                    // Reassemble param chains for symbolic registration updates
                    std::deque<Param<T>> paramChain;
                    for (const auto& cp : collectedParams) paramChain.push_back(cp);
                    newEquation.rightSide.push_back(Monomial<T>(paramChainCopy));
                }
            }
        }
        // CASE 2: Companion matrix parameter updates
        else {
            for (Equation<T>& e : equations) {
                Param<T> leftSide = e.leftSide;

                for (Monomial<T>& m : e.rightSide) {
                    bool found = false;
                    size_t pos = targetName.find_last_of('_');
                    std::string originalName = (pos == std::string::npos) ? targetName : targetName.substr(0, pos);

                    for (const Param<T>& p : m.parameters) {
                        if (p.name == originalName) found = true;
                    }
                    if (!found) continue;

                    std::deque<Param<T>> collectedParams = {};

                    for (const Param<T>& p : m.parameters) {
                        if (p.name != originalName) {
                            collectedParams.push_back(p);
                        } else {
                            std::string derivedName = makeParamName(leftSide.name, order);
                            if (derivedName.front() == 'F') {
                                collectedParams.push_back(findParamByName(derivedName, derivatives));
                            } else {
                                collectedParams.push_back(findParamByName(derivedName, inputs));
                            }
                        }
                    }

                    if (collectedParams.empty()) continue;

                    std::deque<Param<T>> paramChainCopy = collectedParams;

                    Tensor<T> finalMonomialTensor = evaluateMonomialSequence(collectedParams, targetParam.indexNames);

                    equationSum = equationSum + finalMonomialTensor;

                    std::deque<Param<T>> paramChain;
                    for (const auto& cp : collectedParams) paramChain.push_back(cp);
                    newEquation.rightSide.push_back(Monomial<T>(paramChainCopy));
                }
            }
        }

        targetParam.tensor = equationSum;
        newEquations.push_back(newEquation);
    }

    for (auto& ne : newEquations) {
        equations.push_back(ne);
    }
}


// --- Helpers ---
// 1. Contract tensors by finding contractable tensors in the remaining array
//    via equal active orders ergo equal indices
// 2. Realign the final product to the expected shape by permutating the indices

template <typename T>
static Tensor<T> evaluateMonomialSequence(
    std::deque<Param<T>>& collectedParams, 
    const std::deque<std::string>& targetIndexNames // Changed from std::vector to std::deque
) {
    if (collectedParams.empty()) {
        throw std::invalid_argument("Monomial evaluation failed: Parameter pool is empty.");
    }

    Param<T> runningParam = collectedParams.front();
    collectedParams.pop_front();

    // DEBUG PARAMETER --------------
    int step = 1;

    while (!collectedParams.empty()) {
        bool foundMatch = false;
        std::deque<std::string> nextIndexNames;
        Tensor<T> stepResult;

        auto targetIt = collectedParams.begin();
        for (; targetIt != collectedParams.end(); ++targetIt) {
            for (const std::string& runIdx : runningParam.indexNames) {
                if (std::find(targetIt->indexNames.begin(), targetIt->indexNames.end(), runIdx) != targetIt->indexNames.end()) {
                    foundMatch = true;
                    break;
                }
            }
            if (foundMatch) break; 
        }

        // ------------------- DEBUG
        Param<T> rightParam;

        if (foundMatch) {
            // -------------------
            rightParam = *targetIt;
            //---------------------

            Param<T> matchingParam = *targetIt;
            collectedParams.erase(targetIt);
            stepResult = contractByMetadata(runningParam, matchingParam, nextIndexNames);
        } else {
            //-----------------------------------
            rightParam = collectedParams.front();
            //-----------------------------------

            Param<T> standaloneParam = collectedParams.front();
            collectedParams.pop_front();
            stepResult = contractByMetadata(runningParam, standaloneParam, nextIndexNames);
        }

        // Execute the contraction
        stepResult = contractByMetadata(runningParam, rightParam, nextIndexNames);

        step++;


        runningParam.tensor = stepResult;
        runningParam.indexNames = nextIndexNames;
    }

    return alignToTarget(runningParam.tensor, runningParam.indexNames, targetIndexNames);
}

template<typename T>
static Tensor<T> contractByMetadata(
    const Param<T>& paramA, 
    const Param<T>& paramB,
    std::deque<std::string>& outIndexNames
) {
    std::vector<size_t> axesA;
    std::vector<size_t> axesB;
    outIndexNames.clear();
    std::vector<std::string> seen = {};

    // Find physical axis index offsets where string index names intersect
    // Also skip over indices that have been seen before (e.g. there are multiple 'i')
    for (size_t posA = 0; posA < paramA.indexNames.size(); ++posA) {
        for (size_t posB = 0; posB < paramB.indexNames.size(); ++posB) {
            if (paramA.indexNames[posA] == paramB.indexNames[posB] && 
                std::find(seen.begin(), seen.end(), paramA.indexNames[posA]) == seen.end()) {
                seen.push_back(paramA.indexNames[posA]);
                axesA.push_back(posA);
                axesB.push_back(posB);
            }
        }
    }

    // Determine the lingering uncontracted index names layout
    for (size_t i = 0; i < paramA.indexNames.size(); ++i) {
        if (std::find(axesA.begin(), axesA.end(), i) == axesA.end()) {
            outIndexNames.push_back(paramA.indexNames[i]);
        }
    }
    for (size_t j = 0; j < paramB.indexNames.size(); ++j) {
        if (std::find(axesB.begin(), axesB.end(), j) == axesB.end()) {
            outIndexNames.push_back(paramB.indexNames[j]);
        }
    }

    // Multi-dimensional contract execution engine
    return Tensor<T>::productGeneralContraction(paramA.tensor, paramB.tensor, axesA, axesB);
}

template<typename T>
static Tensor<T> alignToTarget(
    const Tensor<T>& sourceTensor,
    const std::deque<std::string>& sourceIndexNames,
    const std::deque<std::string>& targetIndexNames
) {
    if (std::equal(sourceIndexNames.begin(), sourceIndexNames.end(), targetIndexNames.begin(), targetIndexNames.end())) {
        return sourceTensor;
    }

    std::vector<size_t> permutationMap;
    for (const std::string& targetIdx : targetIndexNames) {
        auto it = std::find(sourceIndexNames.begin(), sourceIndexNames.end(), targetIdx);
        if (it == sourceIndexNames.end()) {
            throw std::runtime_error("Critical axis alignment error: Missing structural index name dimension during target mapping.");
        }
        permutationMap.push_back(std::distance(sourceIndexNames.begin(), it));
    }
    return Tensor<T>::permuteAxes(sourceTensor, permutationMap);
}