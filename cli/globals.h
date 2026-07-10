#ifndef CLI_GLOBALS_H
#define CLI_GLOBALS_H

#include <string>

// Process-wide pipeline state, owned and set once by main() in main.cpp.
// console.cpp and pipeline.cpp both read these rather than threading them
// through every function call.
extern bool g_verbose;

// <directory containing this executable>/generator and its parent. All
// generator inputs/outputs and the g++ invocations live here, not the
// process's current working directory - see getExecutableDir() in
// pipeline.cpp for why.
extern std::string g_resourceDir;
extern std::string g_genDir;

#endif
