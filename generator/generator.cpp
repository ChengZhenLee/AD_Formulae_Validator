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
    result += "#include \"utils.h\"\n";
    result += "#include \"user_function.h\"\n";
    result += "#include \"readWrite.h\"\n";
    result += "#include \"configManager.h\"\n";
    result += "\n";

    // Include the functions
    result += "template <typename T>\n";
    result += "std::vector<Param<T>> runADDrivers(X_t<T>& x); \n\n";
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
    outFile << "\n\n";

    for (int i = 1; i <= order; i++) {
        char mode = sequence[i - 1];
        if (mode == 't') {
            outFile << generateTangent(order - i + 1, sequence, XADNested, YADNested); 
        }
        else if (mode == 'a') {
            outFile << generateAdjoint(order - i + 1, sequence, XADNested, YADNested);
        }
        outFile << "\n";
    }

    outFile << generateInterface(sequence, XADNested, YADNested);
    outFile << "\n";

    outFile << generateMain(sequence);
    outFile.close();
}


// Generate main function that:
// 1. reads a file of parameters 
// 2. calls the AD functions
// 3. writes the parameters to the output file
std::string generateMain(std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string type = cm.getType();

    std::string result = "";

    result += "int main(int argc, char** argv) {\n";

    // Load the configs
    result += "\tConfigManager::getInstance().load(\"generator/configs.txt\");\n\n";

    // Read the parameters
    result += std::format(
        "\tstd::vector<Param<{}>> parameters = readParameters<{}>(\"generator/parameters.bin\");\n\n",
        type, type
    );

    // Run and write the results
    result += std::format(
        "\tauto results = runADDrivers<{}>(parameters);\n"
        "\twriteParameters<{}>(results, \"generator/results.bin\");\n\n", 
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

    // Get the config manager
    result += "\tConfigManager& cm = ConfigManager::getInstance();\n";

    // Initialize the AD types 
    result += std::format("\t{} x;\n", XADNested);
    result += std::format("\t{} y;\n", YADNested);

    // Resize x and y deques
    result += "\tx.resize(cm.getXShape());\n";
    result += "\ty.resize(cm.getYShape());\n";

    // Seed primal values first
    result += "\tseedPrimal(x, parameters);\n\n";

    // Register inputs for all adjoint tapes
    result += generateRegisterInputString(sequence);

    // Run the drivers
    result += std::format("\t{}<{}>(x, y, parameters);\n", functionName, type);

    // Reset all adjoint tapes
    result += "\t// Clean up all recorded tapes before returning\n";
    result += generateResetTapeString(sequence); 
    result += "\n";

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

    // Seed values
    result += std::format("\tseedADForOrder(x, y, parameters, {}, \"{}\");\n",
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
    result += std::format("\textractADForOrder(x, y, parameters, {}, \"{}\");\n",
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
    size_t outerAdjointOrder = sequence.length() - sequence.find('a');

    // Function signature
    result += "template <typename T>\n";
    result += std::format("void {0}({1}& x, {2}& y, std::vector<Param<T>> parameters) {{\n", 
        functionName, XADNested, YADNested);

    // Run the lower layer function (executes the underlying math or lower AD layers)
    if (curOrder == 1) {
        result += std::format("\t{}(x, y);\n", subFunctionName);
    } else {
        result += std::format("\t{}(x, y, parameters);\n", subFunctionName);
    }

    // Initialize the tape
    result += std::format("\t{}::tape::init_adjoints();\n", curADType);

    // Seed values
    result += std::format("\tseedADForOrder(x, y, parameters, {}, \"{}\");\n",
        curOrder, sequence);

    // Interpret the tape
    result += std::format("\t{}::tape::interpret();\n", curADType);
    
    // Extract derivatives
    result += std::format("\textractADForOrder(x, y, parameters, {}, \"{}\");\n",
        curOrder, sequence);
    
    // Extract primal values if it is the base layer (curOrder == 1)
    if (curOrder == 1) {
        result += std::format("\textractPrimal(y, parameters);\n");
    }

    result += "}\n";

    return result;
}


// Generate helper drivers to calculate the derivatives
void generateHelperHeader(std::string filename) {
    std::ofstream outFile(filename);
    ConfigManager cm = ConfigManager::getInstance();
    size_t order = cm.getSequence().length();
    std::string sequence = "";
    for (size_t i = 0; i < order; i++) {
        sequence += 't';
    }
    std::string XADNested = generateXNestedADType(sequence);
    std::string YADNested = generateYNestedADType(sequence);
    std::string result = "";

    // Include the defines
    result += "#ifndef ADHELPER_H\n";
    result += "#define ADHELPER_H\n";
    result += "\n";

    // Include the required headers
    result += "#include \"structures.h\"\n";
    result += "#include \"utils.h\"\n";
    result += "#include \"user_function.h\"\n";
    result += "#include \"readWrite.h\"\n";
    result += "#include \"configManager.h\"\n";
    result += "\n";

    // Include the functions
    result += "template <typename T>\n";
    result += "std::vector<Param<T>> runADHelper(X_t<T>& x); \n\n";
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


// Generate the Helper drivers
void generateHelperDrivers(std::string filename) {
    std::ofstream outFile(filename);

    ConfigManager cm = ConfigManager::getInstance();
    size_t order = cm.getSequence().length();
    std::string sequence = "";
    for (size_t i = 0; i < order; i++) {
        sequence += 't';
    }
    std::string XADNested = generateXNestedADType(sequence);
    std::string YADNested = generateYNestedADType(sequence);

    // Write the includes
    outFile << "#include \"adHelper.h\"\n";
    outFile << "\n\n";

    // The helpers only use tangent mode
    for (int i = 1; i <= order; i++) {
        outFile << generateTangent(order - i + 1, sequence, XADNested, YADNested); 
        outFile << "\n";
    }

    outFile << generateHelperInterface(sequence, XADNested, YADNested);
    outFile << "\n";

    outFile << generateHelperMain(sequence);
    outFile.close();
}


// Generates the helper function's interface
std::string generateHelperInterface(std::string sequence, std::string XADNested, std::string YADNested) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string type = cm.getType();

    size_t order = sequence.length();
    std::string functionName = getCurrentLayerFunctionName(order);
    std::string result = "";

    // Code to seed the parameters
    result += "template<typename T> \n";
    result += "std::vector<Param<T>> runHelperDrivers(std::vector<Param<T>> &parameters) {\n";

    // Get the config manager
    result += "\tConfigManager& cm = ConfigManager::getInstance();\n";

    // Initialize the AD types 
    result += std::format("\t{} x;\n", XADNested);
    result += format("\t{} y;\n", YADNested);

    // Resize x and y deques
    result += std::format("\tx.resize(cm.getXShape());\n");
    result += std::format("\ty.resize(cm.getYShape());\n");

    // Run the drivers
    result += std::format("\t{}<{}>(x, y, parameters);\n", functionName, type);

    // Return the parameters
    result += "\treturn parameters;\n";
    result += "}\n";

    return result;
}


// Generate helper's main function
std::string generateHelperMain(std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string type = cm.getType();

    std::string result = "";

    result += "int main(int argc, char** argv) {\n";

    // Load the configs
    result += "\tConfigManager::getInstance().load(\"generator/configs.txt\");\n\n";

    // Read the parameters
    result += std::format(
        "\tstd::vector<Param<{}>> parameters = readParameters<{}>(\"generator/derivatives_seeds.bin\");\n\n",
        type, type
    );

    // Run and write the results
    result += std::format(
        "\tauto results = runHelperDrivers<{}>(parameters);\n"
        "\twriteParameters<{}>(results, \"generator/derivatives.bin\");\n\n", 
        type, type
    );

    result += "\treturn 0;\n";
    result += "}\n";

    return result;
}


// ----------------
// Helper Functions
// ----------------

// Generate the complex nested AD types
std::string generateNestedADType(std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    std::string V = std::to_string(cm.getTangentShape());
    std::string U = std::to_string(cm.getAdjointShape());

    std::string result = cm.getType();
    for (int i = 0; i < sequence.length(); i++) {
        if (sequence[i] == 't') {
            result = std::format("T_t<{},{}>", result, V);
        } else if (sequence[i] == 'a') {
            result = std::format("A_t<{},{}>", result, U);
        }
    }
    return result;
}

// Generate the x variable for complex nested AD types
std::string generateXNestedADType(std::string sequence) {
    std::string complexType = generateNestedADType(sequence);
    return std::format("X_t<{}>", complexType);
}

// Generate the y variable for complex nested AD types
std::string generateYNestedADType(std::string sequence) {
    std::string complexType = generateNestedADType(sequence);
    return std::format("Y_t<{}>", complexType);
}

// Determine the current layer's AD type
std::string getCurrentLayerADType(size_t curOrder, std::string sequence) {
    std::string subsequence = sequence.substr(0, sequence.length() - curOrder + 1);
    return generateNestedADType(subsequence);
}

// Determine the function name of the current layer
std::string getCurrentLayerFunctionName(size_t curOrder) {
    ConfigManager cm = ConfigManager::getInstance();
    if (curOrder == 0) return cm.getPrimalFunctionName();
    return std::format("AD_F_{}", curOrder);
}

// Generate the string required to register an input for all adjoint layers
std::string generateRegisterInputString(std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t xShape = cm.getXShape();
    std::string result = "\tfor (size_t i = 0; i < " + std::to_string(xShape) + "; i++) {\n";
    
    for (int i = sequence.length() - 1; i >= 0; i--) {
        if (sequence[i] == 't') continue;

        std::string xVariable = "x[i]";
        size_t curOrder = sequence.length() - i;
        for (int i = 0; i < curOrder - 1; i++) {
            xVariable += ".value()";
        }
        result += std::format("\t\t{}.register_input();\n", xVariable);
    }
    result += "\t}\n";
    return result;
}

// Generate the string required to reset all tapes in a (highest level) adjoint mode driver
std::string generateResetTapeString(std::string sequence) {
    std::string result;
    std::string curSubstring;
    for (int i = 0; i < sequence.length(); i++) {
        if (sequence[i] == 'a') {
            curSubstring = sequence.substr(0, i + 1);
            result += std::format("\t{}::tape::reset();\n", generateNestedADType(curSubstring));
        }
    }
    return result;
}


// Generate the string required to initialize the adjoints of all tapes in a (highest level) adjoint mode driver
std::string generateInitAdjointsString(std::string sequence) {
    std::string result;
    std::string curSubstring;
    for (int i = 0; i < sequence.length(); i++) {
        if (sequence[i] == 'a') {
            curSubstring = sequence.substr(0, i + 1);
            result = std::format("\t{}::tape::init_adjoints();\n", generateNestedADType(curSubstring)) + result;
        }
    }
    return result;
}
