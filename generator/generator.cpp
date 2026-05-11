#include <format>
#include "generator.h"
#include "utils.h"


std::string generateDrivers() {
}


std::string generateTangent(size_t curOrder, std::string sequence, std::string ADNested) {
    std::string result = "";
    std::string functionName = getCurrentLayerFunctionName(curOrder);
    std::string subFunctionName = getCurrentLayerFunctionName(curOrder - 1);

    // Function signature
    result += std::format("void {0}({1}& x, {1}& y, std::vector<Params> parameters) {{\n", 
        functionName, ADNested);

    // Seed the primal values if it is the outer most order
    if (curOrder == sequence.length()) {
        // Seed the primal values
        result += std::format("\tseedPrimal(x, parameters);\n");
    }
    
    // Seed values
    result += std::format("\tseedADForOrder(x, parameters, {}, {});\n",
        curOrder, sequence);

    // Run the lower layer function
    // Don't pass in parameters to the primal function
    if (curOrder == 1) {
        result += std::format("\t{}(x, y);\n", 
            subFunctionName);
    }
    result += std::format("\t{}(x, y, parameters);\n",
        subFunctionName);

    // Extract
    result += std::format("\textractADForOrder({}& y, parameters, {}, {});\n",
        ADNested, curOrder, sequence);
    
    // Extract primal values if it is the last layer
    if (curOrder == 1) {
        //extract primal
        result += std::format("\textractPrimal(y, parameters);\n");
    }

    // Close the function
    result += std::format("}}\n");

    return result;
}


std::string generateAdjoint(size_t curOrder, std::string sequence, std::string ADNested) {
    std::string result = "";
    std::string functionName = getCurrentLayerFunctionName(curOrder);
    std::string subFunctionName = getCurrentLayerFunctionName(curOrder - 1);
    std::string curADType = getCurrentLayerADType(ADNested, curOrder, sequence);

    // Function signature
    result += std::format("void {0}({1}& x, {1}& y, std::vector<Params> parameters) {{\n", 
        functionName, ADNested);

    // Seed the primal values if it is the outer most order
    if (curOrder == sequence.length()) {
        // Seed the primal values
        result += std::format("\tseedPrimal(x, parameters);\n");
    }

    // Register values
    result += generateRegisterInputString(curOrder);

    // Run the lower layer function
    // Don't pass in parameters to the primal function
    if (curOrder == 1) {
        result += std::format("\t{}(x, y);\n", 
            subFunctionName);
    }
    result += std::format("\t{}(x, y, parameters);\n",
        subFunctionName);

    // Seed values
    result += std::format("\tseedADForOrder(x, parameters, {}, {});\n",
        curOrder, sequence);

    // Interpret
    result += std::format("\t{}::tape::interpret();\n", curADType);
    
    // Extract
    result += std::format("\textractADForOrder({}& y, parameters, {}, {});\n",
        ADNested, curOrder, sequence);
    
    // Extract primal values if it is the last layer
    if (curOrder == 1) {
        //extract primal
        result += std::format("\textractPrimal(y, parameters);\n");
    }

    // TODO: Reset all tapes only at the top level
    if (curOrder == sequence.length()) {
        result += generateResetTapeString(sequence);
    }

    result += "}}\n";

    return result;
}