# AD Formulae Validator

A command-line tool for validating a hand-written algorithmic differentiation
(AD) formula against a numerically-generated AD driver, for any combination
of tangent (forward) and adjoint (reverse) differentiation orders.

You write a primal function `f(x, y)` once. The tool generates and compiles
an AD driver for the differentiation mode you ask for, runs it, independently
re-derives the expected derivatives symbolically, compares the two, and
prints a single `VALID` / `INVALID` verdict.

The default primal function is the Lighthouse function, and the default AD sequence is "ta", which is "tangent-over-adjoint" mode.

## How it works

The validator computes the same higher-order derivatives two ways, from two
independent code paths, and checks that they agree within tolerance:

1. **Numeric AD driver** (`generator/codegen/generator.cpp` → generated
   `generator/adDrivers.cpp`). For a requested sequence such as `ta`, the
   primal type `T` is wrapped in nested AD types from `include/ad/ad.h`
   (`ad::tangent_t<T,V>` for each `t`, `ad::adjoint_t<T,U>` for each `a`,
   composed outermost-to-innermost). `f()` is then run *once*, unmodified,
   through this nested type. Each layer seeds its own tangent directions or
   records/interprets its own adjoint tape, so the nested type mechanically
   produces every derivative of the requested order — this is standard
   operator-overloading AD, and it is the thing being validated.
2. **Symbolic formula driver** (`generator/symbolic/formulaDriver.hpp`). This is a
   second, independent implementation that never touches the AD library. It
   starts from the base equation `y = f(x)` and applies the multivariate
   chain/product rule order-by-order: a `t` step differentiates every term
   of the current equation set with respect to the next tangent direction
   (sum rule over monomials, product rule within a monomial); an `a` step
   propagates adjoint sensitivities backward through the same equation set.
   Each application is expressed as a **tensor contraction over named
   indices** (`i` = input, `j` = output, `v_k`/`u_k` = the tangent/adjoint
   seed index introduced at order `k`), i.e. plain Einstein summation, so
   the result only depends on the shapes and the chain rule — not on any AD
   machinery. `generator/symbolic/formulaDriver.hpp` documents this in detail.

Both drivers are run on the *same* random seed values (`generator/parameters.bin`,
`generator/derivatives_seeds.bin`), and `generator/symbolic/validator.h` compares
their output tensors element-by-element. Because the two derivations share
no code, agreement is strong evidence that the operator-overloading AD
library computes the derivative the chain rule says it should for the given
formula and sequence — which is the property a hand-written AD formula
needs for a paper to rely on it.

## Requirements

- A C++20-capable `g++` on `PATH` (MinGW-w64 / w64devkit on Windows, GCC or
  Clang's `g++` on Linux/macOS). The tool compiles generated code at runtime,
  so this is a hard requirement, not just a build-time one — it checks for
  `g++` at startup and fails fast with a clear message if it's missing.
- CMake (3.10+) and a C++20 compiler, only if you're building the tool
  itself from source.

## Quick start

```
ADValidator.exe --sequence=at
```

```
Starting AD Validator...
Compiling AD drivers...
Generating seeds...
Running drivers...
Running formula driver...
Validating...

Sequence: at
Result: VALID
Equations written to <path>/generator/equations.txt
```

Exit code is `0` when the result is `VALID`, non-zero otherwise (see
[Exit codes](#exit-codes)).

## Validating your own function

1. Open `generator/user_function.h` and replace the body of `f()` with your
   formula. Keep the signature exactly as-is:
   ```cpp
   template<typename T>
   void f(X_t<T>& x, Y_t<T>& y);
   ```
   It gets instantiated for both plain numeric types and nested AD types, so
   it must stay generic over `T`.
2. In `generator/configs.txt`, set:
   - `x` — number of input vector elements (must match what `f()` expects)
   - `y` — number of output vector elements (must match what `f()` expects)
   - `V` — tangent (forward) seed size
   - `U` — adjoint (reverse) seed size
   - `T` — numeric type: `double`, `float`, or `int`
   - `sequence` — the default AD sequence to validate when `--sequence` isn't
     passed on the command line (see below)
3. Only operations with a generic AD-aware overload are safe to use inside
   `f()`: the arithmetic operators (`+ - * /`, unary and binary), and
   `sin`, `cos`, `tan`, `exp`, `sqrt`, `pow(x, int n)`, `log10`, `fabs`.
   Anything called via `std::` (e.g. `std::sin`) will not compile against the
   AD types — avoid it.
4. Run the tool with the sequence you want to check (see below). If `f()`
   doesn't compile, the tool prints the relevant compiler error lines
   directly instead of the raw build log.

## Choosing a differentiation sequence

`--sequence=<mode>` overrides the `sequence` value in `generator/configs.txt`.
Each character is `t` (tangent) or `a` (adjoint).

The sequence is read left-to-right as **outermost-to-innermost**: the first
character is the last-applied (outermost) differentiation, wrapping
everything before it; the last character is the first-applied (innermost)
one. This matches standard "X over Y" AD terminology — e.g. `--sequence=ta`
means "tangent-over-adjoint": adjoint runs first (innermost), tangent wraps
it as the second, outer derivative.

Sequences can be any length and any mix of `t`/`a`, e.g. `a`, `t`, `at`,
`ta`, `aa`, `tt`, `tat`, `aat`, etc.

Passing `--sequence` also **persists**: it rewrites the `sequence` line in
`generator/configs.txt`, so it becomes the new default for subsequent runs
that omit `--sequence`, instead of reverting to whatever was configured
before. Every other line in `configs.txt` is left untouched.

## Command-line options

| Option | Description |
|---|---|
| `--sequence=<mode>` | Overrides `sequence` from `generator/configs.txt` for this run, and persists as the new default in that file. |
| `--verbose`, `-v` | Show pipeline stage-by-stage progress instead of the short one-line status messages. |
| `--help`, `-h` | Show usage and exit. |

## Output files

All written under the `generator/` folder next to the executable:

- `generator/validation_results.txt` — per-parameter pass/fail breakdown
  (written on every run, most useful when the result is `INVALID`).
- `generator/equations.txt` — every symbolic equation the formula driver
  derived for the current sequence, one per line as
  `LEFTSIDE = monomial + monomial + ...`, so you can inspect exactly what
  was derived rather than only the final verdict.
- `generator/*_build.log` — full compiler output for the generated AD/helper
  drivers, kept even on success; only surfaced to the console on failure.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Result is VALID |
| 1 | AD driver failed to generate/compile, `g++` not found, or a bad CLI argument |
| 2 | Helper driver failed to generate/compile |
| 3 | Seed generation failed |
| 4 | Running a compiled driver failed |
| 5 | Formula driver or validation threw an exception |
| 6 | Ran successfully but the result is INVALID |

## Distributing the tool

Don't hand someone the whole `build/` directory — it's full of CMake/Ninja
build machinery they don't need. Building `ADValidator` also assembles a
`build/dist/` folder alongside it:

```
dist/
  ADValidator.exe        (or ADValidator on Linux/macOS)
  generator/              — full generator/ source tree (see Project layout)
  include/                — the AD library, JSON, and Eigen headers it depends on
```

`generator/` is copied wholesale, so a peer reviewing `dist/` sees the same
source as the repository, including `generator/codegen/` and
`generator/symbolic/` — only `ADValidator.exe` needs any of it, and only the
flat files (`structures.h`, `utils.h`/`.hpp`, `readWrite.h`,
`configManager.h`, `user_function.h`, `configs.txt`) at runtime, when it
compiles the AD/helper drivers it generates.

Zip up `dist/` and send it as-is. The executable locates `generator/` and
`include/` relative to its own location, not the current working directory,
so `dist/` can be copied anywhere and run from any folder — it isn't tied to
this repository. A peer only needs to edit `generator/user_function.h` and
`generator/configs.txt` inside their copy of `dist/`, then run the
executable.

## Project layout

- `main.cpp` — CLI entry point: argument parsing plus calling the six
  pipeline stages (compile AD driver, compile helper driver, seed, run,
  formula, validate) in order. The stages themselves and all terminal
  output live in `cli/`, not here.
- `cli/pipeline.h` / `pipeline.cpp` — the pipeline stage implementations,
  plus the executable-location/build-caching plumbing they share.
- `cli/console.h` / `console.cpp` — all terminal output: colored
  stage/verdict printing (auto-disabled when not writing to a terminal, or
  when `NO_COLOR` is set), `--help` text, and compiler-error reporting.
- `cli/globals.h` — the small set of process-wide state (verbose flag,
  resource/generator directory paths) shared between `main.cpp` and `cli/`.
- `generator/codegen/generator.cpp` / `generator.h` — generates the AD
  driver and helper driver source code for a given sequence.
- `generator/symbolic/formulaDriver.hpp` — independently re-derives the
  expected derivatives symbolically, for cross-validation against the
  numeric AD driver.
- `generator/symbolic/validator.h` — compares AD driver output against
  formula driver output within tolerance.
- `generator/configManager.h` — loads `generator/configs.txt` and handles
  the `--sequence` CLI override.
- `generator/structures.h` — shared types (`Param`, `Tensor`, `Equation`,
  `Monomial`, AD type aliases).
- `generator/readWrite.h` — binary (de)serialization of parameters between
  the AD driver, helper driver, and formula driver.
- `generator/utils.h` / `utils.hpp` — parameter-tree generation (one `Param`
  per input/output per differentiation order) and the recursive
  seed/extract routines that map tensor data into/out of the nested AD
  types at each order.
- `generator/user_function.h` — where you define the formula under test.
- `include/ad/ad.h` — the tangent/adjoint AD type library.

`generator/codegen/` and `generator/symbolic/` group source by role;
`structures.h`, `utils.h`/`utils.hpp`, `readWrite.h`, `configManager.h`,
and `user_function.h` stay flat directly under `generator/` because the
*generated* AD/helper driver source `#include`s them unqualified and is
compiled with `generator/` itself as its only project-local `-I` path (see
`handleADDriverCompilation` in `cli/pipeline.cpp`) — moving one of those
into a subfolder would break that compile.

`generator/adDrivers.{h,cpp}` and `generator/adHelper.{h,cpp}` are **not**
source files to read or edit and are not checked into git (see
`.gitignore`) — they're `generator/codegen/generator.cpp`'s output,
written fresh into `generator/` on whichever run needs them.
