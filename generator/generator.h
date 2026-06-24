#ifndef GENERATOR_H
#define GENERATOR_H

#include <string>

// Main AD driver generation functions
void generateADHeader(std::string filename);
void generateADDrivers(std::string filename);

// Helper driver generation functions (for computing derivatives)
void generateHelperHeader(std::string filename);
void generateHelperDrivers(std::string filename);

// Code generation utilities
std::string generateInterface(
    std::string sequence, 
    std::string XADNested, 
    std::string YADNested
);

std::string generateHelperInterface(
    std::string sequence, 
    std::string XADNested, 
    std::string YADNested
);

std::string generateTangent(
    size_t curOrder, 
    std::string sequence,
    std::string XADNested,
    std::string YADNested
);

std::string generateAdjoint(
    size_t curOrder, 
    std::string sequence,
    std::string XADNested,
    std::string YADNested
);

std::string generateMain(std::string sequence);
std::string generateHelperMain(std::string sequence);

// AD type generation helpers
std::string generateNestedADType(std::string sequence);
std::string generateXNestedADType(std::string sequence);
std::string generateYNestedADType(std::string sequence);
std::string getCurrentLayerADType(size_t curOrder, std::string sequence);
std::string getCurrentLayerFunctionName(size_t curOrder);

// Code string generation helpers
std::string generateRegisterInputString(std::string sequence);
std::string generateResetTapeString(std::string sequence);
std::string generateInitAdjointsString(std::string sequence);

#endif