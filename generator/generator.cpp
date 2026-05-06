#include <format>
#include "generator.h"
#include "utils.h"
#include "configManager.h"


std::string generateDrivers() {
}


std::string generateTangent(size_t curOrder, std::string sequence, std::string ADNested) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t xShape = cm.getXShape();
    size_t yShape = cm.getYShape();

    std::string result = "";

    std::string functionName = getCurrentLayerFunctionName(curOrder);
    std::string subFunctionName = getCurrentLayerFunctionName(curOrder - 1);

    // Function signature
    result += std::format("void {0}({1}& x, {1}& y, std::vector<Params> parameters) {{\n", 
        functionName, ADNested);

    // Seed the primal values if it is the outer most order
    if (curOrder == sequence.length()) {
        // Seed the primal values

        result += std::format("\tfor (int i = 0; i < {}; i++) {{\n", xShape);
        result += std::format("\t\tseedPrimal(x(i), parameters, {{i}});\n");
        result += std::format("\t}}");
    }
    
    // Seed values
    result += std::format("\tseedADForOrder({}& x, parameters, {}, {{i}});\n",
        ADNested, curOrder, sequence);

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
    }

    return result;
}


std::string generateAdjoint() {
    std::string result = "";

    // If curOrder == totalOrder, seed x (top layer)

    // Register current layer inputs

    // Run the lower layer function

    // Seed

    // Interpret

    // Extract

    // If curOrder == 1, extract y (bottom layer)

    //

    return result;
}