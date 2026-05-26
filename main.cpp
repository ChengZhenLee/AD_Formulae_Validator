#include <iostream>
#include "generator/generator.h"
#include "generator/configManager.h"

int main(int argc, char** argv) {
    std::cout << "Generator-focused AD Validator runner\n";

    // Load generator config and produce driver/header outputs
    ConfigManager::getInstance().load("generator/configs.txt");
    try {
        generateADHeader();
        generateADDrivers();
        std::cout << "Generated AD drivers and header in generator/\n";
    } catch (const std::exception &e) {
        std::cerr << "Generator failed: " << e.what() << "\n";
        return 2;
    }

    return 0;
}
