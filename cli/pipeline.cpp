#include "pipeline.h"
#include "console.h"
#include "globals.h"
#include <filesystem>
#include <format>
#include <cstdlib>
#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <climits>
#endif
#include "generator/codegen/generator.h"
#include "generator/utils.h"
#include "generator/readWrite.h"
#include "generator/symbolic/formulaDriver.hpp"
#include "generator/symbolic/validator.h"

namespace fs = std::filesystem;

#ifdef _WIN32
const char* NULL_DEVICE = "NUL";
#else
const char* NULL_DEVICE = "/dev/null";
#endif

#ifdef _WIN32
const char* EXE_EXT = ".exe";
#else
const char* EXE_EXT = "";
#endif

std::string getExecutableDir() {
#if defined(_WIN32)
    char buffer[MAX_PATH];
    DWORD len = GetModuleFileNameA(NULL, buffer, MAX_PATH);
    if (len == 0 || len == MAX_PATH) {
        return fs::current_path().string();
    }
    return fs::path(buffer).parent_path().string();
#elif defined(__APPLE__)
    char buffer[PATH_MAX];
    uint32_t size = sizeof(buffer);
    if (_NSGetExecutablePath(buffer, &size) != 0) {
        return fs::current_path().string();
    }
    std::error_code ec;
    fs::path resolved = fs::canonical(buffer, ec);
    if (ec) return fs::path(buffer).parent_path().string();
    return resolved.parent_path().string();
#elif defined(__linux__)
    std::error_code ec;
    fs::path resolved = fs::canonical("/proc/self/exe", ec);
    if (ec) return fs::current_path().string();
    return resolved.parent_path().string();
#else
    return fs::current_path().string();
#endif
}

std::string quotePath(const std::string& path) {
    return "\"" + path + "\"";
}

bool handleADDriverCompilation(const std::string& current_sequence, const PipelineCache& cache) {
    std::string exePath = g_genDir + "/adDrivers" + EXE_EXT;
    bool needs_rebuild = (current_sequence != cache.last_sequence) || !fs::exists(exePath);

    if (!needs_rebuild) {
        announceSkip(1, "sequence unchanged");
        return true;
    }

    announceStage(1, "Compiling AD driver", "Generating & compiling AD driver", "sequence changed");
    try {
        generateADHeader(g_genDir + "/adDrivers.h");
        generateADDrivers(g_genDir + "/adDrivers.cpp");
    } catch (const std::exception &e) {
        printError(std::format("Generation failed: {}\n", e.what()));
        return false;
    }

    std::string logPath = g_genDir + "/adDrivers_build.log";
    std::string compileCmd = "g++ -std=c++20 -g -O0 " + quotePath(g_genDir + "/adDrivers.cpp") +
        " -I " + quotePath(g_genDir) +
        " -I " + quotePath(g_resourceDir + "/include/ad") +
        " -I " + quotePath(g_resourceDir + "/include") +
        " -o " + quotePath(exePath) +
        " > " + quotePath(logPath) + " 2>&1";
    if (std::system(compileCmd.c_str()) != 0) {
        printError("\nYour function in generator/user_function.h failed to compile:\n\n");
        printCompileErrors(logPath);
        return false;
    }
    return true;
}

bool handleHelperDriverCompilation(size_t current_length, const PipelineCache& cache) {
    std::string exePath = g_genDir + "/adHelper" + EXE_EXT;
    bool needs_rebuild = (current_length != cache.last_length) || !fs::exists(exePath);

    if (!needs_rebuild) {
        announceSkip(2, "sequence length unchanged");
        return true;
    }

    announceStage(2, "Compiling helper driver", "Generating & compiling helper driver", "sequence length changed");
    try {
        generateHelperHeader(g_genDir + "/adHelper.h");
        generateHelperDrivers(g_genDir + "/adHelper.cpp");
    } catch (const std::exception &e) {
        printError(std::format("Helper driver generation failed: {}\n", e.what()));
        return false;
    }

    std::string logPath = g_genDir + "/adHelper_build.log";
    std::string helperCompileCmd = "g++ -std=c++20 -g -O0 " + quotePath(g_genDir + "/adHelper.cpp") +
        " -I " + quotePath(g_genDir) +
        " -I " + quotePath(g_resourceDir + "/include/ad") +
        " -I " + quotePath(g_resourceDir + "/include") +
        " -o " + quotePath(exePath) +
        " > " + quotePath(logPath) + " 2>&1";
    if (std::system(helperCompileCmd.c_str()) != 0) {
        printError("\nYour function in generator/user_function.h failed to compile:\n\n");
        printCompileErrors(logPath);
        return false;
    }
    return true;
}

ParameterContainer generateSeeds(const std::string& target_type, ConfigManager& cm, X_t<double>& x_input) {
    announceStage(3, "Generating seeds", "Generating input & derivative seeds");
    ParameterContainer parameters;

    if (target_type == "double") {
        parameters = generateRandomSeeds<double>(cm.getSequence(), x_input);
        writeParameters<double>(std::get<std::vector<Param<double>>>(parameters), g_genDir + "/parameters.bin");

        auto derivSeeds = generateDerivativeSeeds<double>(cm.getSequence(), x_input);
        writeParameters<double>(derivSeeds, g_genDir + "/derivatives_seeds.bin");
    }
    else if (target_type == "float") {
        X_t<float> x_input_float(x_input.begin(), x_input.end());
        parameters = generateRandomSeeds<float>(cm.getSequence(), x_input_float);
        writeParameters<float>(std::get<std::vector<Param<float>>>(parameters), g_genDir + "/parameters.bin");

        auto derivSeeds = generateRandomSeeds<float>(cm.getSequence(), x_input_float);
        writeParameters<float>(derivSeeds, g_genDir + "/derivatives_seeds.bin");
    }
    else if (target_type == "int") {
        X_t<int> x_input_int(x_input.begin(), x_input.end());
        parameters = generateRandomSeeds<int>(cm.getSequence(), x_input_int);
        writeParameters<int>(std::get<std::vector<Param<int>>>(parameters), g_genDir + "/parameters.bin");

        auto derivSeeds = generateRandomSeeds<int>(cm.getSequence(), x_input_int);
        writeParameters<int>(derivSeeds, g_genDir + "/derivatives_seeds.bin");
    }
    return parameters;
}

bool runExecutables() {
    announceStage(4, "Running drivers", "Running compiled drivers");

    if (std::system(quotePath(g_genDir + "/adHelper" + EXE_EXT).c_str()) != 0) {
        printError("Execution of helper driver failed!\n");
        return false;
    }
    if (std::system(quotePath(g_genDir + "/adDrivers" + EXE_EXT).c_str()) != 0) {
        printError("Execution of AD driver failed!\n");
        return false;
    }
    return true;
}

ParameterContainer runFormula(const std::string& target_type, const ParameterContainer& parameters) {
    announceStage(5, "Running formula driver", "Running formula driver");
    std::string equationsPath = g_genDir + "/equations.txt";
    if (target_type == "double") {
        return runFormulaDriver<double>(std::get<std::vector<Param<double>>>(parameters), equationsPath);
    } else if (target_type == "float") {
        return runFormulaDriver<float>(std::get<std::vector<Param<float>>>(parameters), equationsPath);
    } else {
        return runFormulaDriver<int>(std::get<std::vector<Param<int>>>(parameters), equationsPath);
    }
}

bool validate(const std::string& target_type, const ParameterContainer& formulaResults) {
    announceStage(6, "Validating", "Validating results");
    bool allValid = false;
    std::string resultsPath = g_genDir + "/results.bin";
    std::string outPath = g_genDir + "/validation_results.txt";
    if (target_type == "double") {
        allValid = validateParameters<double>(std::get<std::vector<Param<double>>>(formulaResults), resultsPath, outPath);
    } else if (target_type == "float") {
        allValid = validateParameters<float>(std::get<std::vector<Param<float>>>(formulaResults), resultsPath, outPath);
    } else if (target_type == "int") {
        allValid = validateParameters<int>(std::get<std::vector<Param<int>>>(formulaResults), resultsPath, outPath);
    }
    return allValid;
}
