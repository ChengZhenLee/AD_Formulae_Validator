#include <cstddef>
#include <deque>
#include <unordered_set>
#include "structures.h"
#include "configManager.h"
#include "utils.h"

// ---------- TODO USED FOR TESTING
template<typename T>
void printIndexNames(const std::deque<std::string>& indices) {
    std::cout << "[";
    for (size_t i = 0; i < indices.size(); ++i) {
        std::cout << indices[i] << (i + 1 < indices.size() ? ", " : "");
    }
    std::cout << "]";
}

template<typename T>
std::vector<Param<T>> runFormulaDriver(std::vector<Param<T>> parameters) {
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
        std::cout << "]\n"
                  << "Index Names: [";
        for (std::string s : parameters[i].indexNames) std::cout << s << " ";
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
    std::cout << "   DUMPING ALL DERIVATIVES FROM FILE     \n";
    std::cout << "========================================\n";
    for (size_t i = 0; i < derivatives.size(); ++i) {
        std::cout << "Index [" << i << "] "
                  << "Name: " << derivatives[i].name << " | "
                  << "Role: " << (derivatives[i].role == ParamRole::Input ? "Input" : "Output") << " | "
                  << "Active Orders: [";
        for (std::string s : parameters[i].indexNames) std::cout << s << " ";
        std::cout << "]\n";
    }
    std::cout << "========================================\n\n";
    // === END DEBUG DUMP ===

    std::deque<Equation<T>> equations = {};

    formulaDriver(sequence, equations, derivatives, inputs, outputs);

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

    debugDumpEquations<T>(equations);
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
                        if (m.parameters[i].name.front() == 'F') continue;
                        collectedParams.push_back(m.parameters[j]);
                    }
                }

                if (collectedParams.empty()) continue;

                Tensor<T> finalMonomialTensor = evaluateMonomialSequence(collectedParams, targetParam.indexNames);
                
                equationSum = equationSum + finalMonomialTensor;
                
                std::deque<Param<T>> paramChain;
                for (const auto& cp : collectedParams) paramChain.push_back(cp);
                newEquation.rightSide.push_back(Monomial<T>(paramChain));
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
                    std::cout << "\n--------------------------------------------------\n";
                    std::cout << "[ADJOINT MODE - CASE 1] Processing Monomial for: " << targetName << "\n";
                    std::cout << "Originating Equation Left Side: " << leftSide.name << "\n";
                    std::cout << "--------------------------------------------------\n";

                    for (const Param<T>& p : m.parameters) {
                        std::cout << "  • Examining Monomial Param: " << p.name << " | Indices: [ ";
                        for (auto& s : p.indexNames) std::cout << s << " "; std::cout << "]\n";

                        // Standard filter for primal X or its tangent/higher-order mutations
                        if (p.name == "X" || p.name.rfind("X_", 0) == 0) {
                            std::cout << "    --> SKIP: Matches differentiation target group.\n";
                            continue;
                        }
                        
                        if (p.name[0] == 'F') {
                            std::string lookupName = makeParamName(p.name, p.highestOrder + 1);
                            std::cout << "    --> ACTION: Treating as Derivative. Generated Name: " << lookupName << "\n";
                            try {
                                Param<T> found = findParamByName(lookupName, derivatives);
                                collectedParams.push_back(found);
                                std::cout << "        SUCCESS: Collected from 'derivatives'.\n";
                            } catch (const std::exception& e) {
                                std::cout << "        ⚠️ FAILURE: Could not find " << lookupName << " in 'derivatives' pool!\n";
                            }
                        } else {
                            std::cout << "    --> ACTION: Treating as Input/Intermediate baseline: " << p.name << "\n";
                            try {
                                Param<T> found = findParamByName(p.name, inputs);
                                collectedParams.push_back(found);
                                std::cout << "        SUCCESS: Collected from 'inputs'.\n";
                            } catch (const std::exception& e) {
                                std::cout << "        ⚠️ FAILURE: Could not find " << p.name << " in 'inputs' pool!\n";
                            }
                        }
                    }

                    // Append the incoming adjoint seed node tracker
                    std::string seedName = makeParamName(leftSide.name, order);
                    std::cout << "  • Appending Adjoint Seed Tracker. Generated Name: " << seedName << "\n";
                    try {
                        Param<T> foundSeed = findParamByName(seedName, inputs);
                        collectedParams.push_back(foundSeed);
                        std::cout << "    SUCCESS: Collected seed from 'inputs'.\n";
                    } catch (const std::exception& e) {
                        std::cout << "    ⚠️ FAILURE: Could not find seed " << seedName << " in 'inputs' pool!\n";
                    }

                    std::cout << "Total parameters pushed to pipeline queue: " << collectedParams.size() << "\n";
                    std::cout << "--------------------------------------------------\n";

                    if (collectedParams.empty()) {
                        std::cout << "🛑 WARNING: collectedParams is completely EMPTY. Skipping execution step.\n";
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
    // === ENHANCED INITIAL POOL DIAGNOSTIC ===
    std::cout << "\n======================================================\n";
    std::cout << " INITIAL PARAMETER POOL IN COLLECTED_PARAMS (Size: " << collectedParams.size() << ")\n";
    std::cout << "======================================================\n";
    
    for (size_t idx = 0; idx < collectedParams.size(); ++idx) {
        const auto& param = collectedParams[idx];
        std::cout << "  [" << idx << "] Name: " << (param.name.empty() ? "(intermediate)" : param.name)
                  << " | Index Names: [ ";
        for (const std::string& s : param.indexNames) {
            std::cout << s << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "======================================================\n\n";

    if (collectedParams.empty()) {
        throw std::invalid_argument("Monomial evaluation failed: Parameter pool is empty.");
    }

    std::cout << "\n======================================================\n";
    std::cout << "        STARTING MONOMIAL EVALUATION PIPELINE          \n";
    std::cout << "======================================================\n";


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

        // --- PIPELINE TELEMETRY PRE-CONTRACTION ---
        std::cout << "\n--- [STEP " << step << "] ---\n";
        std::cout << "Left Operand Name:  " << (runningParam.name.empty() ? "(intermediate)" : runningParam.name) << "\n";
        std::cout << "  Tracked Strings: [ "; for (auto& s : runningParam.indexNames) std::cout << s << " "; std::cout << "]\n";
        std::cout << "  Tensor Shape:    [ "; for (size_t d : runningParam.tensor.shape) std::cout << d << " "; std::cout << "] (Rank " << runningParam.tensor.shape.size() << ")\n\n";

        std::cout << "Right Operand Name: " << rightParam.name << "\n";
        std::cout << "  Tracked Strings: [ "; for (auto& s : rightParam.indexNames) std::cout << s << " "; std::cout << "]\n";
        std::cout << "  Tensor Shape:    [ "; for (size_t d : rightParam.tensor.shape) std::cout << d << " "; std::cout << "] (Rank " << rightParam.tensor.shape.size() << ")\n";

        // Execute the contraction
        stepResult = contractByMetadata(runningParam, rightParam, nextIndexNames);

        // --- PIPELINE TELEMETRY POST-CONTRACTION ---
        std::cout << "\nResulting Output:\n";
        std::cout << "  Tracked Strings: [ "; for (auto& s : nextIndexNames) std::cout << s << " "; std::cout << "]\n";
        std::cout << "  Tensor Shape:    [ "; for (size_t d : stepResult.shape) std::cout << d << " "; std::cout << "] (Rank " << stepResult.shape.size() << ")\n";
        std::cout << "------------------------------------------------------\n";

        // CRITICAL CHECK: Verify string tracking size matches physical rank
        if (nextIndexNames.size() != stepResult.shape.size()) {
            std::cerr << "\n⚠️ [DESYNCHRONIZATION DETECTED AT STEP " << step << "]:\n"
                      << "   The metadata string tracker has " << nextIndexNames.size() << " indices,\n"
                      << "   but the physical underlying Tensor has rank " << stepResult.shape.size() << "!\n"
                      << "   Check your multi-axis contraction logic inside contractByMetadata.\n\n";
        }
        step++;


        runningParam.tensor = stepResult;
        runningParam.indexNames = nextIndexNames;
    }

    std::cout << "\nFinal alignment check against target names: [ ";
    for (auto& s : targetIndexNames) std::cout << s << " "; std::cout << "]\n";
    std::cout << "======================================================\n\n";

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
    // ------------ DEBUG
    // Add this debug block:
    if (targetIndexNames.size() != sourceTensor.shape.size()) {
        printIndexNames<T>(sourceIndexNames);
        std::cout << "\n";
        printIndexNames<T>(targetIndexNames);
        std::cerr << "\n[RANK MISMATCH CRASH]:\n"
                  << " -> Target Expected Rank (Size of targetIndexNames): " << targetIndexNames.size() << "\n"
                  << " -> Actual Source Tensor Rank (sourceTensor.rank()): " << sourceTensor.shape.size() << "\n"
                  << " -> Source Index Names Size: " << sourceIndexNames.size() << "\n";
    }


    if (std::equal(sourceIndexNames.begin(), sourceIndexNames.end(), targetIndexNames.begin(), targetIndexNames.end())) {
        return sourceTensor;
    }

    std::vector<size_t> permutationMap;
    for (const std::string& targetIdx : targetIndexNames) {
        auto it = std::find(sourceIndexNames.begin(), sourceIndexNames.end(), targetIdx);
        if (it == sourceIndexNames.end()) {
            std::cerr << "\n========== ALIGNMENT CRASH DEBUG TAPE ==========\n";
            std::cerr << "[CRITICAL ERROR]: Missing expected structural index name dimension.\n";
            std::cerr << "[MISSING TARGET INDEX]: " << targetIdx << "\n\n";
            
            std::cerr << "[SOURCE STATE]:\n";
            std::cerr << "  -> Index Names:  [ ";
            for (const std::string& val : sourceIndexNames) std::cerr << val << " ";
            std::cerr << "]\n";
            
            std::cerr << "[EXPECTED TARGET STATE]:\n";
            std::cerr << "  -> Index Names:  [ ";
            for (const std::string& val : targetIndexNames) std::cerr << val << " ";
            std::cerr << "]\n";
            std::cerr << "================================================\n\n";
            
            throw std::runtime_error("Critical axis alignment error: Missing structural index name dimension during target mapping.");
        }
        permutationMap.push_back(std::distance(sourceIndexNames.begin(), it));
    }
    return Tensor<T>::permuteAxes(sourceTensor, permutationMap);
}

// TODO:
// Used for testing

template<typename T>
void debugDumpEquations(const std::deque<Equation<T>>& equations) {
    std::cout << "\n==================================================\n";
    std::cout << "          SYMBOLIC EQUATION METADATA DUMP         \n";
    std::cout << "==================================================\n";
    
    for (size_t idx = 0; idx < equations.size(); ++idx) {
        const auto& eq = equations[idx];
        
        std::cout << "Eq [" << idx << "]: " << eq.leftSide.name;
        printIndexNames<T>(eq.leftSide.indexNames);
        std::cout << " = ";

        for (size_t m_idx = 0; m_idx < eq.rightSide.size(); ++m_idx) {
            const auto& monomial = eq.rightSide[m_idx];
            
            // Print terms within a single monomial chain
            for (size_t p_idx = 0; p_idx < monomial.parameters.size(); ++p_idx) {
                const auto& param = monomial.parameters[p_idx];
                std::cout << param.name;
                printIndexNames<T>(param.indexNames);
                
                if (p_idx + 1 < monomial.parameters.size()) {
                    std::cout << " * ";
                }
            }
            
            if (m_idx + 1 < eq.rightSide.size()) {
                std::cout << " + ";
            }
        }
        std::cout << "\n--------------------------------------------------\n";
    }
    std::cout << "==================================================\n\n";
}