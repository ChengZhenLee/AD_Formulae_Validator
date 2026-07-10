#include <iostream>
#include <filesystem>
#include <format>
#include <cstdlib>
#include "cli/console.h"
#include "cli/pipeline.h"
#include "cli/globals.h"
#include "generator/configManager.h"

namespace fs = std::filesystem;

// Definitions for the globals declared extern in cli/globals.h.
bool g_verbose = false;
std::string g_resourceDir;
std::string g_genDir;

// ============================================================================
// MAIN PIPELINE EXECUTION
//
// This file only parses arguments and orchestrates the six pipeline stages
// in order (see cli/pipeline.h) plus the couple of one-time setup steps
// that precede them (locating the resource folder, checking for g++). The
// stages themselves, and all terminal-output formatting, live in cli/.
// ============================================================================

int main(int argc, char** argv) {
    initConsole();

    std::string sequenceOverride = "";

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg.rfind("--sequence=", 0) == 0) {
            sequenceOverride = arg.substr(std::string("--sequence=").length());
        } else if (arg == "--verbose" || arg == "-v") {
            g_verbose = true;
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else {
            printError(std::format("Unrecognized argument: {}\n\n", arg));
            printUsage(argv[0]);
            return 1;
        }
    }

    g_resourceDir = getExecutableDir();
    g_genDir = g_resourceDir + "/generator";

    // The generated AD/helper driver programs are compiled with relative
    // paths like "generator/parameters.bin" baked in (see generator.cpp),
    // resolved against whatever directory they're launched from. Since they
    // run as child processes of this one and inherit its working directory,
    // pin that directory to the resource folder now so the drivers find
    // their files regardless of where the user invoked this executable from.
    std::error_code chdirErr;
    fs::current_path(g_resourceDir, chdirErr);
    if (chdirErr) {
        printError(std::format("Error: could not switch to resource directory \"{}\": {}\n",
            g_resourceDir, chdirErr.message()));
        return 1;
    }

    std::string gppCheckCmd = std::string("g++ --version > ") + NULL_DEVICE + " 2>&1";
    if (std::system(gppCheckCmd.c_str()) != 0) {
        printError("Error: g++ was not found on PATH.\n"
                    "This tool compiles generated code at runtime and requires a\n"
                    "C++20-capable g++ (e.g. MinGW-w64 / w64devkit) to be installed and on PATH.\n");
        return 1;
    }

    std::cout << colorize(Color::BOLD, "Starting AD Validator...") << "\n" << std::flush;

    std::string configsPath = g_genDir + "/configs.txt";
    ConfigManager::getInstance().load(configsPath);

    if (!sequenceOverride.empty()) {
        std::string configuredSequence = ConfigManager::getInstance().getRawSequence();
        ConfigManager::getInstance().set("sequence", sequenceOverride);

        // A --sequence override becomes the new default for future runs,
        // rather than reverting to whatever configs.txt said before.
        if (sequenceOverride != configuredSequence) {
            ConfigManager::getInstance().persistSequence(configsPath, sequenceOverride);
            std::cout << colorize(Color::DIM, std::format(
                "Updated default sequence in generator/configs.txt: {} -> {}\n",
                configuredSequence, sequenceOverride));
        }
    }
    ConfigManager cm = ConfigManager::getInstance();

    std::string current_sequence = cm.getSequence();
    size_t current_length = current_sequence.length();
    std::string target_type = cm.getType();
    std::string display_sequence = cm.getRawSequence();

    // Load caching information
    PipelineCache cache;
    const std::string cache_file = g_genDir + "/.pipeline_cache";
    cache.load(cache_file);

    // Stage 1 & 2: Conditional compilation
    if (!handleADDriverCompilation(current_sequence, cache)) return 1;
    if (!handleHelperDriverCompilation(current_length, cache)) return 2;

    // Update and save cache right away if compilations succeeded
    cache.last_sequence = current_sequence;
    cache.last_length = current_length;
    cache.save(cache_file);

    // Stage 3: Input generation
    X_t<double> x_input(cm.getXShape());
    fillWithRandomValues(x_input);

    ParameterContainer parameters;
    try {
        parameters = generateSeeds(target_type, cm, x_input);
    } catch (const std::exception &e) {
        printError(std::format("Seed generation failed: {}\n", e.what()));
        return 3;
    }

    // Stage 4: Execute drivers
    if (!runExecutables()) return 4;

    // Stage 5 & 6: Formula & validation
    bool allValid = false;
    try {
        ParameterContainer formulaResults = runFormula(target_type, parameters);
        allValid = validate(target_type, formulaResults);
    } catch (const std::exception &e) {
        printError(std::format("Pipeline execution/validation failure: {}\n", e.what()));
        return 5;
    }

    std::string equationsPath = g_genDir + "/equations.txt";
    std::string validationPath = g_genDir + "/validation_results.txt";

    std::cout << "\n" << std::string(40, '-') << "\n";
    std::cout << colorize(Color::BOLD, "Sequence: ") << display_sequence << "\n";
    std::cout << colorize(Color::BOLD, "Result:   ")
               << (allValid ? colorize(Color::GREEN, "VALID") : colorize(Color::RED, "INVALID"))
               << "\n";
    std::cout << std::string(40, '-') << "\n";
    std::cout << "Equations written to " << equationsPath << "\n";
    std::cout << "Validation results for individual parameters written to " << validationPath << "\n";
    if (!allValid) {
        std::cout << colorize(Color::YELLOW, std::format("See {} for per-parameter details.\n", validationPath));
    }

    return allValid ? 0 : 6;
}
