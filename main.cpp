#include <iostream>
#include <random>
#include <variant>
#include <fstream>
#include <filesystem>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <climits>
#endif
#include "generator/generator.h"
#include "generator/configManager.h"
#include "generator/structures.h"
#include "generator/readWrite.h"
#include "generator/formulaDriver.hpp"
#include "generator/validator.h"

namespace fs = std::filesystem;

using ParameterContainer = std::variant<
    std::vector<Param<double>>,
    std::vector<Param<float>>,
    std::vector<Param<int>>
>;

// Quiet by default: only short stage messages and the final verdict print
// unless --verbose is set.
static bool g_verbose = false;

// All generator inputs/outputs and the g++ invocation live under
// <directory containing this executable>/generator, not the process's
// current working directory. This lets the tool run correctly whether it's
// launched from the build tree during development or from a standalone
// redistributable folder (see cmake dist/ packaging) handed to a peer.
static std::string g_resourceDir;
static std::string g_genDir;

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

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

// Null device to discard command output when probing for a tool on PATH.
#ifdef _WIN32
static const char* NULL_DEVICE = "NUL";
#else
static const char* NULL_DEVICE = "/dev/null";
#endif

// Extension for compiler-produced executables: none on Linux/macOS.
#ifdef _WIN32
static const char* EXE_EXT = ".exe";
#else
static const char* EXE_EXT = "";
#endif

// Wraps a path in double quotes so it survives std::system()'s shell
// invocation even when it contains spaces (e.g. a "AD Validator" folder).
std::string quotePath(const std::string& path) {
    return "\"" + path + "\"";
}

template<typename T>
void fillWithRandomValues(X_t<T>& x_input) {
    std::random_device rd;
    std::mt19937 gen(rd());

    using DistributionType = std::conditional_t<
        std::is_floating_point_v<T>,
        std::uniform_real_distribution<T>,
        std::uniform_int_distribution<T>
    >;
    DistributionType dist(static_cast<T>(0), static_cast<T>(100));

    for (auto& element : x_input) {
        element = dist(gen);
    }
}

// Struct to track our build states
struct PipelineCache {
    std::string last_sequence = "";
    size_t last_length = 0;

    void load(const std::string& filepath) {
        std::ifstream f(filepath);
        if (f.is_open()) {
            std::getline(f, last_sequence);
            f >> last_length;
        }
    }

    void save(const std::string& filepath) const {
        std::ofstream f(filepath);
        if (f.is_open()) {
            f << last_sequence << "\n" << last_length;
        }
    }
};

void printUsage(const char* argv0) {
    std::cout <<
        "Usage: " << argv0 << " [--sequence=<mode>] [--verbose] [--help]\n\n"
        "Validates a user-supplied primal function (generator/user_function.h) by\n"
        "comparing a numerically-computed AD driver against a symbolically-derived\n"
        "formula, for an arbitrary combination of tangent/adjoint differentiation\n"
        "orders. Prints a single VALID / INVALID verdict.\n\n"
        "Options:\n"
        "  --sequence=<mode>  Overrides the 'sequence' value from generator/configs.txt.\n"
        "                     Each character is 't' (tangent) or 'a' (adjoint), read\n"
        "                     left to right as outermost-to-innermost: the FIRST\n"
        "                     character is the last-applied, outermost differentiation\n"
        "                     (wraps everything before it); the LAST character is the\n"
        "                     first-applied, innermost one.\n"
        "                       e.g. --sequence=ta means \"tangent over adjoint\":\n"
        "                       adjoint runs first (innermost), tangent wraps it.\n"
        "  --verbose, -v      Show pipeline stage-by-stage progress.\n"
        "  --help, -h         Show this message.\n\n"
        "To test your own function: edit the body of f() in generator/user_function.h\n"
        "(keep the signature `template<typename T> void f(X_t<T>& x, Y_t<T>& y)`\n"
        "unchanged), set input/output sizes and precision in generator/configs.txt\n"
        "(x, y, T), then run this tool with the sequence you want to check.\n"
        "Only these operations are AD-aware inside f(): + - * / (unary and binary),\n"
        "sin, cos, tan, exp, sqrt, pow(x, int n), log10, fabs. Avoid std::sin etc.\n\n"
        "This tool compiles generated code at runtime, so a C++20-capable g++\n"
        "(MinGW-w64 on Windows, GCC/Clang on Linux/macOS) must be installed and on PATH.\n\n"
        "generator/ and include/ next to this executable are its resource folder -\n"
        "keep them alongside it if you copy this tool elsewhere.\n";
}

// Reads a build log and prints just the compiler error lines (or, if none are
// found, the tail of the log) so a peer's own mistake in user_function.h is
// easy to spot without wading through the full g++ output.
void printCompileErrors(const std::string& logPath) {
    std::ifstream log(logPath);
    if (!log.is_open()) {
        std::cerr << "  (could not open build log at " << logPath << ")\n";
        return;
    }

    std::vector<std::string> allLines;
    std::vector<std::string> errorLines;
    std::string line;
    while (std::getline(log, line)) {
        allLines.push_back(line);
        if (line.find(" error:") != std::string::npos) {
            errorLines.push_back(line);
        }
    }

    const std::vector<std::string>& toShow = errorLines.empty() ? allLines : errorLines;
    size_t start = toShow.size() > 20 ? toShow.size() - 20 : 0;
    for (size_t i = start; i < toShow.size(); i++) {
        std::cerr << "  " << toShow[i] << "\n";
    }
    if (start > 0) {
        std::cerr << "  ... (" << start << " earlier line(s) omitted)\n";
    }
    std::cerr << "\nFull compiler output: " << logPath << "\n";
}

// ============================================================================
// PIPELINE STAGES
// ============================================================================

bool handleADDriverCompilation(const std::string& current_sequence, const PipelineCache& cache) {
    std::string exePath = g_genDir + "/adDrivers" + EXE_EXT;
    bool needs_rebuild = (current_sequence != cache.last_sequence) || !fs::exists(exePath);

    if (!needs_rebuild) {
        if (g_verbose) std::cout << "--> Sequence unchanged. Skipping AD driver generation & compilation.\n";
        return true;
    }

    std::cout << "Compiling AD drivers...\n" << std::flush;
    if (g_verbose) std::cout << "\n[1/6] Generating & Compiling AD drivers (Sequence Changed)...\n";
    try {
        generateADHeader(g_genDir + "/adDrivers.h");
        generateADDrivers(g_genDir + "/adDrivers.cpp");
    } catch (const std::exception &e) {
        std::cerr << "Generation failed: " << e.what() << "\n";
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
        std::cerr << "\nYour function in generator/user_function.h failed to compile:\n\n";
        printCompileErrors(logPath);
        return false;
    }
    return true;
}

bool handleHelperDriverCompilation(size_t current_length, const PipelineCache& cache) {
    std::string exePath = g_genDir + "/adHelper" + EXE_EXT;
    bool needs_rebuild = (current_length != cache.last_length) || !fs::exists(exePath);

    if (!needs_rebuild) {
        if (g_verbose) std::cout << "--> Sequence length unchanged. Skipping Helper driver generation & compilation.\n";
        return true;
    }

    std::cout << "Compiling helper drivers...\n" << std::flush;
    if (g_verbose) std::cout << "\n[1.5/6] Generating & Compiling helper drivers (Length Changed)...\n";
    try {
        generateHelperHeader(g_genDir + "/adHelper.h");
        generateHelperDrivers(g_genDir + "/adHelper.cpp");
    } catch (const std::exception &e) {
        std::cerr << "Helper driver generation failed: " << e.what() << "\n";
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
        std::cerr << "\nYour function in generator/user_function.h failed to compile:\n\n";
        printCompileErrors(logPath);
        return false;
    }
    return true;
}

ParameterContainer generateSeeds(const std::string& target_type, ConfigManager& cm, X_t<double>& x_input) {
    std::cout << "Generating seeds...\n" << std::flush;
    if (g_verbose) std::cout << "\n[2/6] Generating Input & Derivative Seeds...\n";
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
    std::cout << "Running drivers...\n" << std::flush;
    if (g_verbose) std::cout << "\n[4/6] Running compiled drivers...\n";

    if (std::system(quotePath(g_genDir + "/adHelper" + EXE_EXT).c_str()) != 0) {
        std::cerr << "Execution of helper driver failed!\n";
        return false;
    }
    if (std::system(quotePath(g_genDir + "/adDrivers" + EXE_EXT).c_str()) != 0) {
        std::cerr << "Execution of AD driver failed!\n";
        return false;
    }
    return true;
}

ParameterContainer runFormula(const std::string& target_type, const ParameterContainer& parameters) {
    std::cout << "Running formula driver...\n" << std::flush;
    if (g_verbose) std::cout << "\n[5/6] Running formula driver...\n";
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
    std::cout << "Validating...\n" << std::flush;
    if (g_verbose) std::cout << "\n[6/6] Validating results...\n";
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

// ============================================================================
// MAIN PIPELINE EXECUTION
// ============================================================================

int main(int argc, char** argv) {
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
            std::cerr << "Unrecognized argument: " << arg << "\n\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    g_resourceDir = getExecutableDir();
    g_genDir = g_resourceDir + "/generator";

    std::string gppCheckCmd = std::string("g++ --version > ") + NULL_DEVICE + " 2>&1";
    if (std::system(gppCheckCmd.c_str()) != 0) {
        std::cerr << "Error: g++ was not found on PATH.\n"
                     "This tool compiles generated code at runtime and requires a\n"
                     "C++20-capable g++ (e.g. MinGW-w64 / w64devkit) to be installed and on PATH.\n";
        return 1;
    }

    std::cout << "Starting AD Validator...\n" << std::flush;
    if (g_verbose) std::cout << "=== AD Validator Pipeline Started ===\n";

    ConfigManager::getInstance().load(g_genDir + "/configs.txt");
    if (!sequenceOverride.empty()) {
        ConfigManager::getInstance().set("sequence", sequenceOverride);
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

    // Phase 1: Conditional Compilation
    if (!handleADDriverCompilation(current_sequence, cache)) return 1;
    if (!handleHelperDriverCompilation(current_length, cache)) return 2;

    // Update and save cache right away if compilations succeeded
    cache.last_sequence = current_sequence;
    cache.last_length = current_length;
    cache.save(cache_file);

    // Phase 2: Input Generation
    X_t<double> x_input(cm.getXShape());
    fillWithRandomValues(x_input);

    ParameterContainer parameters;
    try {
        parameters = generateSeeds(target_type, cm, x_input);
    } catch (const std::exception &e) {
        std::cerr << "Seed generation failed: " << e.what() << "\n";
        return 3;
    }

    // Phase 3 & 4: Execute drivers
    if (!runExecutables()) return 4;

    // Phase 5 & 6: Formula & Validation
    bool allValid = false;
    try {
        ParameterContainer formulaResults = runFormula(target_type, parameters);
        allValid = validate(target_type, formulaResults);
    } catch (const std::exception &e) {
        std::cerr << "Pipeline execution/validation failure: " << e.what() << "\n";
        return 5;
    }

    std::cout << "\nSequence: " << display_sequence << "\n";
    std::cout << "Result: " << (allValid ? "VALID" : "INVALID") << "\n";
    std::cout << "Equations written to " << g_genDir << "/equations.txt\n";
    std::cout << "Validation results for individual parameters written to " << g_genDir << "/validation_results.txt\n";
    if (!allValid) {
        std::cout << "See " << g_genDir << "/validation_results.txt for per-parameter details.\n";
    }

    return allValid ? 0 : 6;
}
