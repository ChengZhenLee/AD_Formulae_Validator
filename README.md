# AD Formulae Validator

This README describes the runtime flow that is driven by `main.cpp` and the `generator/` folder.

## What this project does

`main.cpp` controls an automatic code-generation pipeline:

- loads generator settings from `generator/configs.txt`
- uses `generator/generator.cpp` to create `generator/adDrivers.h` and `generator/adDrivers.cpp`
- writes `generator/parameters.txt` from manual or random input
- compiles the generated AD driver as `generator/adDrivers.exe`
- runs the compiled AD driver and the formula validation driver

## Key files in `generator/`

- `generator/generator.cpp` / `generator/generator.h`
  - generate AD driver source code and header declarations
  - build the generated `main()` that reads parameters and invokes AD drivers
- `generator/configManager.h` / `generator/configManager.cpp`
  - load runtime settings such as data type, input shape, and AD sequence
- `generator/readWrite.h`
  - read/write parameter vectors to `generator/parameters.txt`
- `generator/formulaDriver.hpp` / `generator/formulaDriver.h`
  - run the formula driver for validation using the same parameters
- `generator/structures.h`, `generator/utils.h`, `generator/user_function.h`
  - shared support code for AD type creation, seeding, and extraction

## Runtime behavior

When `ADValidator` runs, `main.cpp` performs these phases:

1. Load configuration from `generator/configs.txt`
2. Generate the AD driver header and implementation files
3. Populate input parameters (manual or random) and save them to `generator/parameters.txt`
4. Compile the generated driver source into `generator/adDrivers.exe`
5. Execute the generated AD driver and validate results with `runFormulaDriver`

## Input and output files

- `generator/configs.txt` — configuration for type, input size, and AD sequence
- `generator/parameters.txt` — generated seed/input values used by the AD driver
- `generator/adDrivers.h` / `generator/adDrivers.cpp` — generated AD driver code
- `generator/results.txt` — output produced by the generated AD driver

## Notes

- This workflow is centered on the generator code and the runtime orchestration in `main.cpp`.
- The AD driver is emitted dynamically at runtime, compiled, and executed as part of the validator pipeline.
