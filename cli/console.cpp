#include "console.h"
#include "globals.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <format>
#include <cstdlib>
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#else
#include <unistd.h>
#endif

namespace Color {
    bool enabled = false;
}

std::string colorize(const char* code, const std::string& text) {
    if (!Color::enabled) return text;
    return std::string(code) + text + Color::RESET;
}

void initConsole() {
#ifdef _WIN32
    bool isTTY = _isatty(_fileno(stdout)) != 0;
#else
    bool isTTY = isatty(fileno(stdout)) != 0;
#endif
    Color::enabled = isTTY && std::getenv("NO_COLOR") == nullptr;

#ifdef _WIN32
    if (Color::enabled) {
        HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
        DWORD mode = 0;
        if (hOut == INVALID_HANDLE_VALUE || !GetConsoleMode(hOut, &mode) ||
            !SetConsoleMode(hOut, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING)) {
            Color::enabled = false;
        }
    }
#endif
}

void printError(const std::string& message) {
    std::cerr << colorize(Color::RED, message);
}

void announceStage(int stageNum, const std::string& quietLabel, const std::string& verboseLabel, const std::string& reason) {
    if (g_verbose) {
        std::cout << "\n" << colorize(Color::CYAN, std::format("[{}/{}]", stageNum, TOTAL_STAGES)) << " " << verboseLabel;
        if (!reason.empty()) std::cout << " " << colorize(Color::DIM, std::format("({})", reason));
        std::cout << "...\n";
    } else {
        std::cout << quietLabel << "...\n" << std::flush;
    }
}

void announceSkip(int stageNum, const std::string& reason) {
    if (!g_verbose) return;
    std::cout << "\n" << colorize(Color::CYAN, std::format("[{}/{}]", stageNum, TOTAL_STAGES))
               << " " << colorize(Color::DIM, "Skipped - " + reason) << "\n";
}

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
        "                     Also becomes the new default written back to\n"
        "                     generator/configs.txt, so later runs without\n"
        "                     --sequence reuse it instead of the old value.\n"
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
        "keep them alongside it if you copy this tool elsewhere.\n\n"
        "Output is colored automatically when running in a terminal, and left as\n"
        "plain text when redirected to a file/pipe or when NO_COLOR is set.\n";
}

void printCompileErrors(const std::string& logPath) {
    std::ifstream log(logPath);
    if (!log.is_open()) {
        printError(std::format("  (could not open build log at {})\n", logPath));
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
        std::cerr << "  " << colorize(Color::RED, toShow[i]) << "\n";
    }
    if (start > 0) {
        std::cerr << colorize(Color::DIM, std::format("  ... ({} earlier line(s) omitted)\n", start));
    }
    std::cerr << "\nFull compiler output: " << logPath << "\n";
}
