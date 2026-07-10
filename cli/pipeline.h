#ifndef CLI_PIPELINE_H
#define CLI_PIPELINE_H

#include <string>
#include <variant>
#include <vector>
#include <random>
#include <fstream>
#include <type_traits>
#include "generator/structures.h"
#include "generator/configManager.h"

// The three possible shapes the pipeline runs in, selected at runtime by
// the "T" value in generator/configs.txt.
using ParameterContainer = std::variant<
    std::vector<Param<double>>,
    std::vector<Param<float>>,
    std::vector<Param<int>>
>;

// Null device to discard command output when probing for a tool on PATH.
extern const char* NULL_DEVICE;

// Extension for compiler-produced executables: none on Linux/macOS.
extern const char* EXE_EXT;

// Resolves the directory containing the running executable (not the
// process's current working directory - see globals.h).
std::string getExecutableDir();

// Wraps a path in double quotes so it survives std::system()'s shell
// invocation even when it contains spaces (e.g. a "AD Validator" folder).
std::string quotePath(const std::string& path);

// Fills `x_input` with independent uniform random values in [0, 100), used
// to seed the primal input the pipeline validates against.
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

// Tracks the sequence/length the AD and helper drivers were last generated
// and compiled for, persisted to generator/.pipeline_cache, so an unchanged
// sequence between runs can skip regenerating and recompiling them.
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

// ============================================================================
// PIPELINE STAGES
//
// Each function below is one numbered stage of the pipeline (see
// console.h's announceStage/TOTAL_STAGES) and reads/writes g_genDir
// (globals.h) for its inputs/outputs. Defined in pipeline.cpp; called in
// order from main() in main.cpp.
// ============================================================================

// Stage 1: (re)generates and compiles generator/adDrivers.cpp for
// `current_sequence` if it differs from `cache`'s last sequence, or if the
// compiled driver is missing.
bool handleADDriverCompilation(const std::string& current_sequence, const PipelineCache& cache);

// Stage 2: (re)generates and compiles generator/adHelper.cpp if
// `current_length` differs from `cache`'s last length, or if the compiled
// helper is missing.
bool handleHelperDriverCompilation(size_t current_length, const PipelineCache& cache);

// Stage 3: builds and writes the random primal/derivative seeds both
// drivers will run on.
ParameterContainer generateSeeds(const std::string& target_type, ConfigManager& cm, X_t<double>& x_input);

// Stage 4: runs the compiled helper and AD driver executables.
bool runExecutables();

// Stage 5: runs the symbolic formula driver on the same seeds.
ParameterContainer runFormula(const std::string& target_type, const ParameterContainer& parameters);

// Stage 6: compares the AD driver's and formula driver's results.
bool validate(const std::string& target_type, const ParameterContainer& formulaResults);

#endif
