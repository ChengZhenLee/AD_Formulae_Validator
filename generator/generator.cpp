#include <format>
#include <fstream>
#include "generator.h"
#include "utils.h"
#include "configManager.h"


// Generate the header file
void generateADHeader() {
    std::ofstream outFile("adDrivers.h");
    ConfigManager cm = ConfigManager::getInstance();
    std::string sequence = cm.getSequence();
    size_t order = sequence.length();
    std::string XADNested = generateXNestedADType(order, sequence);
    std::string YADNested = generateYNestedADType(order, sequence);
    std::string result = "";

    // Include the defines
    result += "#ifndef ADDRIVER_H\n";
    result += "#define ADDRIVER_H\n";
    result += "\n";

    // Include the required headers
    result += "#include \"structures.h\"\n";
    result += "#include \"ad.h\"\n";
    result += "#include \"types.h\"\n";
    result += "#include \"utils.h\"\n";
    result += "\n";

    // Include the functions
    result += "std::vector<Param<T>> runDrivers(std::string mode, X_t<T>& x); \n";
    for (size_t i = 0; i < order; i++) {
        result += std::format(
            "void {0}({1}& x, {2}& y, std::vector<Param<T>> parameters); \n", 
            getCurrentLayerFunctionName(i + 1), XADNested, YADNested);
    }

    result += "\n#endif\n";

    outFile << result;

    outFile.close();
}


// Generate and writes all the necessary drivers into an output file
void generateADDrivers() {
    std::ofstream outFile("adDrivers.cpp");

    ConfigManager cm = ConfigManager::getInstance();
    std::string sequence = cm.getSequence();

    size_t order = sequence.length();
    std::string XADNested = generateXNestedADType(order, sequence);
    std::string YADNested = generateYNestedADType(order, sequence);

    // Write the includes
    outFile << "#include \"adDriver.h\"\n";

    for (int i = 0; i < order; i++) {
        if (sequence[i] == 't') {
            outFile << generateTangent(i + 1, sequence, XADNested, YADNested); 
        }
        else if (sequence[i] == 'a') {
            outFile << generateAdjoint(i + 1, sequence, XADNested, YADNested);
        }
        outFile << "\n";
    }

    outFile << generateInterface(sequence, XADNested, YADNested);

    outFile.close();
}


std::string generateInterface(std::string sequence, std::string XADNested, std::string YADNested) {
    // void runDrivers(std::string mode)
    size_t order = sequence.length();
    std::string functionName = getCurrentLayerFunctionName(order);
    std::string result = "";

    // Code to seed the parameters
    result += "template<typename T> \n";
    result += "std::vector<Param<T>> runDrivers(std::string mode, X_t<T>& x) {{\n";
    result += "\tstd::vector<Param<T>> parameters; \n";
    result += "\tif (mode == \"random\") {{\n";
    result += std::format("\t\tparameters = generateRandom({}, x); \n", sequence);
    result += "\t}} else if (mode == \"identity\") {{ \n";
    result += std::format("\t\tparameters = generateIdentity({}, x); \n", sequence);
    result += "\t}}\n";

    // Initialize the AD types 
    result += std::format("\t{} x;\n", XADNested);
    result += std::format("\t{} y;\n", YADNested);

    // Run the drivers
    result += std::format("\t{}(x, y, parameters);\n", functionName);

    // Return the parameters
    result += "\treturn parameters;\n";
    result += "}}\n";

    return result;
}


std::string generateTangent(size_t curOrder, std::string sequence, std::string XADNested, std::string YADNested) {
    std::string result = "";
    std::string functionName = getCurrentLayerFunctionName(curOrder);
    std::string subFunctionName = getCurrentLayerFunctionName(curOrder - 1);

    // Function signature
    result += "template <typename T>";
    result += std::format("void {0}({1}& x, {2}& y, std::vector<Param<T>> parameters) {{\n", 
        functionName, XADNested, YADNested);

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
        YADNested, curOrder, sequence);
    
    // Extract primal values if it is the last layer
    if (curOrder == 1) {
        //extract primal
        result += std::format("\textractPrimal(y, parameters);\n");
    }

    // Close the function
    result += std::format("}}\n");

    return result;
}


std::string generateAdjoint(size_t curOrder, std::string sequence, std::string XADNested, std::string YADNested) {
    std::string result = "";
    std::string functionName = getCurrentLayerFunctionName(curOrder);
    std::string subFunctionName = getCurrentLayerFunctionName(curOrder - 1);
    std::string curADType = getCurrentLayerADType(curOrder, sequence);

    // Function signature
    result += "template <typename T>";
    result += std::format("void {0}({1}& x, {2}& y, std::vector<Param<T>> parameters) {{\n", 
        functionName, XADNested, YADNested);

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
        YADNested, curOrder, sequence);
    
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


// Generate formula-based (manual) calculation driver for validation
std::string generateFormulaDriver(std::ofstream& outFile, std::string sequence) {
}


// TODO: Generate main driver function that:
// 1. generates the AD functions and calls them
// 2. generates the Formula functions and calls them
// 3. compares the outputs and prints out the results
std::string generateMainDriver(std::ofstream& outFile, std::string sequence) {
    // TODO: Implement comparison logic
}


// TODO: Generate validation function that compares AD results with formula results