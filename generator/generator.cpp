#include <fstream>
#include "generator.h"
#include "utils.h"
#include "configManager.h"
#include <format>


// Generate the header file
void generateADHeader(std::string filename) {
    std::ofstream outFile(filename);
    ConfigManager cm = ConfigManager::getInstance();
    std::string sequence = cm.getSequence();
    size_t order = sequence.length();
    std::string XADNested = generateXNestedADType(sequence);
    std::string YADNested = generateYNestedADType(sequence);
    std::string result = "";

    // Include the defines
    result += "#ifndef ADDRIVER_H\n";
    result += "#define ADDRIVER_H\n";
    result += "\n";

    // Include the required headers
    result += "#include \"structures.h\"\n";
    result += "\n";

    // Include the functions
    result += "template <typename T>\n";
    result += "std::vector<Param<T>> runADDrivers(std::string mode, X_t<T>& x); \n\n";
    for (size_t i = 0; i < order; i++) {
            result += "template <typename T>\n";
            result += std::format(
                "void {0}({1}& x, {2}& y, std::vector<Param<T>> parameters); \n\n", 
                getCurrentLayerFunctionName(i + 1), XADNested, YADNested);
    }

    result += "\n#endif\n";

    outFile << result;

    outFile.close();
}


// Generate and writes all the necessary drivers into an output file
void generateADDrivers(std::string filename) {
    std::ofstream outFile(filename);

    ConfigManager cm = ConfigManager::getInstance();
    std::string sequence = cm.getSequence();

    size_t order = sequence.length();
    std::string XADNested = generateXNestedADType(sequence);
    std::string YADNested = generateYNestedADType(sequence);

    // Write the includes
    outFile << "#include \"adDrivers.h\"\n";
    outFile << "#include \"structures.h\"\n";
    outFile << "#include \"user_function.h\"";
    outFile << "#include \"readWrite.h\"\n";
    outFile << "\n\n";

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
    outFile << "\n";

    outFile << generateMain(sequence);
    outFile.close();
}


// TODO: Generate main function that:
// 1. reads a file of parameters 
// 2. calls the AD functions
// 3. writes the parameters to the output file
std::string generateMain(std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string type = cm.getType();

    std::string result = "";

    result += "int main(int argc, char** argv) {\n";

    // Read the parameters
    result += std::format(
        "\tstd::vector<Param<{}>> parameters = readParameters<{}>(\"generator/parameters.txt\");\n\n",
        type, type
    );

    // Run and write the results
    result += std::format(
        "\tauto results = runADDrivers<{}>(parameters);\n"
        "\twriteParameters<{}>(results, \"generator/results.txt\");\n\n", 
        type, type
    );

    result += "\treturn 0;\n";
    result += "}\n";

    return result;
}


// Generates the function to be called by the main function
std::string generateInterface(std::string sequence, std::string XADNested, std::string YADNested) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string type = cm.getType();

    size_t order = sequence.length();
    std::string functionName = getCurrentLayerFunctionName(order);
    std::string result = "";

    // Code to seed the parameters
    result += "template<typename T> \n";
    result += "std::vector<Param<T>> runADDrivers(std::vector<Param<T>> &parameters) {\n";

    // Initialize the AD types 
    result += std::format("\t{} x;\n", XADNested);
    result += format("\t{} y;\n", YADNested);

    // Run the drivers
    result += std::format("\t{}<{}>(x, y, parameters);\n", functionName, type);

    // Return the parameters
    result += "\treturn parameters;\n";
    result += "}\n";

    return result;
}


std::string generateTangent(size_t curOrder, std::string sequence, std::string XADNested, std::string YADNested) {
    std::string result = "";
    std::string functionName = getCurrentLayerFunctionName(curOrder);
    std::string subFunctionName = getCurrentLayerFunctionName(curOrder - 1);

    // Function signature
    result += "template <typename T>\n";
    result += std::format("void {0}({1}& x, {2}& y, std::vector<Param<T>> parameters) {{\n", 
        functionName, XADNested, YADNested);

    // Seed the primal values if it is the outer most order
    if (curOrder == sequence.length()) {
        // Seed the primal values
        result += std::format("\tseedPrimal(x, parameters);\n");
    }
    
    // Seed values
    result += std::format("\tseedADForOrder(x, parameters, {}, \"{}\");\n",
        curOrder, sequence);

    // Run the lower layer function
    // Don't pass in parameters to the primal function
    if (curOrder == 1) {
        result += std::format("\t{}(x, y);\n", 
            subFunctionName);
    } else {
        result += std::format("\t{}(x, y, parameters);\n",
            subFunctionName);
    }

    // Extract
    result += std::format("\textractADForOrder(y, parameters, {}, \"{}\");\n",
        curOrder, sequence);
    
    // Extract primal values if it is the last layer
    if (curOrder == 1) {
        //extract primal
        result += std::format("\textractPrimal(y, parameters);\n");
    }

    // Close the function
    result += "}\n";

    return result;
}


std::string generateAdjoint(size_t curOrder, std::string sequence, std::string XADNested, std::string YADNested) {
    std::string result = "";
    std::string functionName = getCurrentLayerFunctionName(curOrder);
    std::string subFunctionName = getCurrentLayerFunctionName(curOrder - 1);
    std::string curADType = getCurrentLayerADType(curOrder, sequence);

    // Function signature
    result += "template <typename T>\n";
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
    } else {
        result += std::format("\t{}(x, y, parameters);\n",
            subFunctionName);
    }

    // Seed values
    result += std::format("\tseedADForOrder(x, parameters, {}, \"{}\");\n",
        curOrder, sequence);

    // Interpret
    result += std::format("\t{}::tape::interpret();\n", curADType);
    
    // Extract
    result += std::format("\textractADForOrder(y, parameters, {}, \"{}\");\n",
        curOrder, sequence);
    
    // Extract primal values if it is the last layer
    if (curOrder == 1) {
        //extract primal
        result += std::format("\textractPrimal(y, parameters);\n");
    }

    // TODO: Reset all tapes only at the top level
    if (curOrder == sequence.length()) {
        result += generateResetTapeString(sequence);
    }

    result += "}\n";

    return result;
}