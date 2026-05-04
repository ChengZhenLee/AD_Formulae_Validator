#include "generator.h"
#include "configManager.h"

// Generate the parameters needed for a given order and sequence of AD
std::vector<Param> generateParameters(int order, std::string sequence) {
    std::vector<Param> result = {};
    char curMode = sequence.front();

    std::vector<Param> recResult = generateParameters(--order, sequence.substr(1));

    if (curMode == 't') {
        // Append the new parameters
    } else if (curMode == 'a') {
        // Append the new parameters
    }

    return result;
}