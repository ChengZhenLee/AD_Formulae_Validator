# AD Formulae Validator

A small C++ project that builds an `ADValidator` executable using CMake and runs automatic validation of adjoint/tangent operations.

## Prerequisites

- `cmake` version 3.10 or newer
- A C++20-compatible compiler (e.g. `g++`, `clang++`)
- `make` or another CMake generator
- Standard build tools (`build-essential` on Debian/Ubuntu)

## Build Instructions

1. Open a terminal and go to the project root:

```bash
cd "AD_Formulae_Validator"
```

2. Create a build directory and enter it:

```bash
mkdir -p build
cd build
```

3. Generate the build files with CMake:

```bash
cmake ..
```

4. Build the executable:

```bash
make
```

## Run the Validator

From the `build` directory:

```bash
./ADValidator
```

The program prints a summary to the terminal and writes detailed validation output to `validation_results.txt` in the current directory.

## Notes

- The project uses C++20.
- If you prefer a different generator, replace `make` with `cmake --build .`.
- The generated executable is named `ADValidator`.
