#ifndef CLI_CONSOLE_H
#define CLI_CONSOLE_H

#include <string>

// ============================================================================
// TERMINAL OUTPUT
//
// Color is opt-out, not opt-in: enabled automatically when stdout is a real
// terminal, disabled automatically when it isn't (e.g. `ADValidator.exe >
// log.txt`, or piped into another program) so redirected output stays
// plain, grep-friendly text. Also disabled if the NO_COLOR environment
// variable is set (https://no-color.org), regardless of terminal. On
// Windows, color additionally requires turning on VT100 escape-sequence
// processing for the current console (off by default on older consoles);
// if that fails, color is disabled rather than leaving raw escape codes in
// the output.
// ============================================================================
namespace Color {
    // Set by initConsole(); read by colorize().
    extern bool enabled;

    constexpr const char* RESET  = "\x1b[0m";
    constexpr const char* BOLD   = "\x1b[1m";
    constexpr const char* DIM    = "\x1b[2m";
    constexpr const char* RED    = "\x1b[31m";
    constexpr const char* GREEN  = "\x1b[32m";
    constexpr const char* YELLOW = "\x1b[33m";
    constexpr const char* CYAN   = "\x1b[36m";
}

// Total number of pipeline stages shown in verbose "[n/6]" headers. Keep in
// sync with the announceStage()/announceSkip() call sites in pipeline.cpp.
constexpr int TOTAL_STAGES = 6;

// Wraps `text` in the given color code, or returns it unchanged if color
// output isn't enabled/available.
std::string colorize(const char* code, const std::string& text);

// Detects whether stdout supports color and, on Windows, enables VT100
// escape-sequence processing so ANSI codes actually render instead of
// printing as literal garbage. Must run before any colorize() calls -
// main() calls this first thing.
void initConsole();

// Prints an error message to stderr in red (if color is enabled). Callers
// still supply their own trailing "\n" so multi-line messages read exactly
// as before.
void printError(const std::string& message);

// Announces the start of pipeline stage `stageNum`. Quiet mode gets one
// short, present-tense line (e.g. "Compiling AD driver..."); verbose mode
// additionally gets a numbered "[n/6]" header and, if given, a short reason
// this stage is running rather than being skipped (e.g. "sequence changed").
void announceStage(int stageNum, const std::string& quietLabel, const std::string& verboseLabel, const std::string& reason = "");

// Announces (verbose mode only) that pipeline stage `stageNum` was skipped
// because its cached output is still valid for this run.
void announceSkip(int stageNum, const std::string& reason);

// Prints --help usage text.
void printUsage(const char* argv0);

// Reads a build log and prints just the compiler error lines (or, if none
// are found, the tail of the log) so a peer's own mistake in
// user_function.h is easy to spot without wading through the full g++
// output.
void printCompileErrors(const std::string& logPath);

#endif
