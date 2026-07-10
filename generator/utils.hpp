#ifndef UTILS_HPP
#define UTILS_HPP

#include "structures.h"
#include "configManager.h"
#include "readWrite.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <format>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <type_traits>
#include "user_function.h"

// ============================================================================
// This file has two independent jobs:
//   1. Build the *parameter tree*: one Param per input/output per
//      differentiation order, with the tensor shapes and Einstein index
//      names the symbolic formula driver needs (generateParameters and the
//      seed generators below it).
//   2. Walk the *nested AD type* the generated driver builds for a given
//      sequence (X_t<...>/A_t<...> composed per generator.cpp) and copy
//      data in and out of it, one differentiation order/layer at a time
//      (seedAD/extractAD and friends). A nested AD object is unwrapped one
//      layer per recursive call — x.tangent(i)/x.adjoint(i) to descend into
//      an order that's active for the Param being seeded, x.value() to skip
//      through an order that isn't — until the innermost primal scalar is
//      reached, at which point the accumulated per-order coordinates are
//      assembled into one flat Tensor index.
// ============================================================================

// Nested AD types must have a .value() member
template<typename T, typename = void>
struct has_value_member : std::false_type {};

template<typename T>
struct has_value_member<T, std::void_t<decltype(std::declval<T>().value())>> : std::true_type {};

template<typename T>
inline constexpr bool isADNestedType = has_value_member<T>::value;

// Recursively builds one Param per input/output per differentiation order
// implied by `sequence` (order-1-is-innermost numbering, see
// ConfigManager::getSequence). Base case (empty sequence) is the bare
// primal X ("i" index) and Y ("j" index). Each recursive step takes every
// Param built for the shorter (order-1) sequence and derives one more from
// it for the current trailing order: a tangent order appends a new "v_k"
// axis of width V to the *right* (seed lives with the input, so this stays
// an Input-role Param); an adjoint order prepends a new "u_k" axis of
// width U to the *left* and flips Input<->Output role (an adjoint's seed
// is supplied on the *output* side and its result read back on the
// *input* side, standard reverse-mode convention). Names accumulate
// "_<order>" suffixes as this recurses, which is what Param's constructor
// parses back out into activeOrders/highestOrder.
template <typename T>
std::vector<Param<T>> generateParameters(std::string sequence) {
    ConfigManager cm = ConfigManager::getInstance();
    size_t xShape = cm.getXShape();
    size_t yShape = cm.getYShape();
    size_t V = cm.getTangentShape();
    size_t U = cm.getAdjointShape();
    size_t order = sequence.size();

    // Base case
    if (order == 0 || sequence.empty()) {
        std::vector<Param<T>> base;
        base.emplace_back(
            std::map<size_t, size_t>{}, 
            std::deque<size_t>{xShape}, 
            "X", 
            ParamRole::Input, 
            std::deque<std::string>{"i"});

        base.emplace_back(
            std::map<size_t, size_t>{}, 
            std::deque<size_t>{yShape}, 
            "Y", 
            ParamRole::Output, 
            std::deque<std::string>{"j"});

        return base;
    }

    char curMode = sequence[order - 1];
    std::vector<Param<T>> result = generateParameters<T>(sequence.substr(0, order - 1));

    size_t currentSize = result.size();
    for (size_t i = 0; i < currentSize; i++) {
        Param p = result[i];
        std::map<size_t, size_t> newOShape = p.orderedShape;
        std::deque<size_t> newShapes = p.tensor.shape;
        ParamRole newRole = p.role;
        std::string newName = std::format("{}_{}", p.name, std::to_string(order));
        std::deque<std::string> newIndexNames = p.indexNames;

        if (curMode == 't') {
            newShapes.push_back(V);
            newOShape[order] = V;
            newIndexNames.push_back(std::format("v_{}", order));
        } else if (curMode == 'a') {
            newShapes.push_front(U);
            newOShape[order] = U;
            newIndexNames.push_front(std::format("u_{}", order));
            newRole = (p.role == ParamRole::Input) ? ParamRole::Output : ParamRole::Input;
        }

        result.emplace_back(newOShape, newShapes, newName, newRole, newIndexNames);
    }

    return result;
}

// Builds the full parameter tree for `sequence` and fills it with the
// random data both drivers will run on: the given primal `x` values go
// into the "X" Param, and every other Param whose name starts with "X"
// (i.e. every tangent/adjoint seed derived from X) gets independent random
// values in the same range as `x`. Written to generator/parameters.bin and
// read back by both the generated AD driver and the formula driver, so
// both operate on identical seed data.
template<typename T>
std::vector<Param<T>> generateRandomSeeds(std::string sequence, const X_t<T>& x) {
    std::vector<Param<T>> parameters = generateParameters<T>(sequence);
    std::random_device rd;
    std::mt19937 gen(rd());

    using DistributionType = std::conditional_t<
        std::is_floating_point_v<T>,
        std::uniform_real_distribution<T>,
        std::uniform_int_distribution<T>
    >;

    DistributionType dist(
        static_cast<T>(*std::min_element(x.begin(), x.end())),
        static_cast<T>(*std::max_element(x.begin(), x.end()))
    );

    for (Param<T>& p : parameters) {
        if (p.name == "X") {
            p.tensor.data.resize(x.size());
            for (int i = 0; i < x.size(); i++) {
                p.tensor.data[i] = x[i];
            }
        } else if (p.name.rfind("X", 0) == 0 && p.name.length() > 1) {
            for (int i = 0; i < p.tensor.data.size(); i++) {
                p.tensor.data[i] = dist(gen);
            }
        }
    }

    return parameters;
}


// Builds the seeds for the *helper* driver (generator/adHelper.cpp), which
// always runs in pure tangent mode regardless of the requested sequence.
// Directional (tangent) derivatives of f are only equal to a true partial
// derivative when the seed direction is a standard basis vector, so X_1,
// X_2, ... here are seeded with identity matrices instead of random data:
// running f() forward through nested-tangent types with an identity-seeded
// direction extracts the literal Jacobian (and higher tangent-only
// derivative tensors) of f at the random primal point used elsewhere. The
// resulting Y_k outputs are renamed to F_k so runFormulaDriver/getDerivatives
// can find them as the ground-truth values for the symbolic "F" (dF/dX)
// placeholders that formulaDriver.hpp's monomials reference — i.e. these
// are the elementary partials the chain rule composes; validating the
// requested sequence checks that composition, not these values themselves
// (they come from ordinary, well-established forward-mode AD).
template<typename T>
std::vector<Param<T>> generateDerivativeSeeds(std::string sequence, const X_t<T>& x) {
    size_t order = sequence.length();
    std::string derivativeSequence = "";
    for (size_t i = 0; i < order; i++) {
        derivativeSequence += 't';
    }
    std::vector<Param<T>> parameters = generateParameters<T>(derivativeSequence);

    for (Param<T>& p : parameters) {
        // Seed X_1, X_2, X_3... with identity matrices
        if (p.name.rfind("X", 0) == 0 && p.name.length() == 3) {
            // Seed with identity
            size_t rows = p.tensor.shape[0];
            size_t cols = p.tensor.shape.size() > 1 ? p.tensor.shape[1] : rows;

            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    size_t idx = p.tensor.getIndex({i, j});
                    p.tensor.data[idx] = (i == j) ? static_cast<T>(1) :static_cast<T>(0);
                }
            }
        }

        // It's easier to change all the "Y" to "F"
        // But be careful not to change the original y=f(x)
        if (p.name.rfind("Y", 0) == 0 && p.name.length() > 1) {
            p.name[0] = 'F';

            // Populate with the correct indices
            std::deque<std::string> newIndexNames = {"j"};
            for (int i = 0; i < p.activeOrders.size(); i++) {
                newIndexNames.push_back("i");
            }
            p.indexNames = newIndexNames;
        }
    }

    return parameters;
}


// Reads back the elementary partial-derivative tensors (F, F_1, F_2, ...)
// computed by the helper driver into generator/derivatives.bin (see
// generateDerivativeSeeds), filtering to just the "F"-named entries.
// `sequence`/`x0`/`h` are unused — kept for interface symmetry with the
// other seed generators; the values were already computed by the time this
// runs, this just loads them.
template<typename T>
std::vector<Param<T>> getDerivatives(std::string sequence, const std::vector<T>& x0, T h) {
    // Read back the results and return the computed derivative parameter values.
    std::vector<Param<T>> results = readParameters<T>("generator/derivatives.bin");
    std::vector<Param<T>> derivatives;
    for (auto& p : results) {
        if (p.name.starts_with("F")) {
            derivatives.push_back(p);
        }
    }

    return derivatives;
}

// Copies data from Param `p`'s tensor into the nested AD object `x`, one
// wrapper layer per recursive call. `x` starts as the outermost/full nested
// type; `curOrder` is a recursion depth counter (1 at the outermost call),
// and `orderLabel` below converts it to the actual differentiation order
// number that outer layer corresponds to (order-1-is-innermost numbering,
// so orderLabel counts *down* from sequence.length() as curOrder — and
// recursion depth — increases). At each layer: if `p` doesn't carry data
// for this order (`!p.isActive`), just unwrap through `.value()` and
// recurse with the same coordinates; if it does, walk every seed index at
// this order via `.tangent(i)`/`.adjoint(i)`, recording each `i` onto
// `rightCoords` (tangent axes were appended on the right by
// generateParameters) or `leftCoords` (adjoint axes were prepended).
// Recursion bottoms out once every layer has been unwound to the bare
// primal type, at which point `leftCoords`+`primalIndex`+`rightCoords` is
// the full multi-index into `p.tensor`, and that scalar is copied into `x`.
template<typename ADNested, typename T>
void seedAD(ADNested& x, const Param<T>& p, size_t curOrder, const std::string& sequence,
    std::deque<size_t>& leftCoords, std::deque<size_t>& rightCoords, size_t& primalIndex) {
    if constexpr(!isADNestedType<ADNested>) {
        if (curOrder > sequence.length()) {
            std::deque<size_t> finalCoords = leftCoords;
            std::reverse(finalCoords.begin(), finalCoords.end());
            finalCoords.push_back(primalIndex);
            finalCoords.insert(finalCoords.end(), rightCoords.begin(), rightCoords.end());
            const size_t flatIndex = p.tensor.getIndex(finalCoords);
            if (flatIndex >= p.tensor.data.size()) {
                throw std::out_of_range("Seed coordinate exceeds tensor storage size.");
            }
            x = p.tensor.data[flatIndex];
            return;
        }
    }

    size_t orderLabel = sequence.length() - curOrder + 1;
    char curType = sequence[orderLabel - 1];
    if (p.isActive(orderLabel)) {
        size_t curShape = p.getShapeAt(orderLabel);
        if (curType == 't') {
            // COMPILER GUARD: Only compile this loop if x actually has a .tangent() method
            if constexpr (requires { x.tangent(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    rightCoords.push_back(i);
                    seedAD(x.tangent(i), p, curOrder + 1, sequence, 
                    leftCoords, rightCoords, primalIndex);
                    rightCoords.pop_back();
                }
            }
        } else if (curType == 'a') {
            // COMPILER GUARD: Only compile this loop if x actually has an .adjoint() method
            if constexpr (requires { x.adjoint(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    leftCoords.push_back(i);
                    seedAD(x.adjoint(i), p, curOrder + 1, sequence, 
                    leftCoords, rightCoords, primalIndex);
                    leftCoords.pop_back();
                }
            }
        }
    } else {
        if constexpr(requires { x.value(); })
        seedAD(x.value(), p, curOrder + 1, sequence, 
        leftCoords, rightCoords, primalIndex);
    }
}

// The read-back counterpart of seedAD: same recursive walk down the nested
// AD object `y`, but copies each reached primal scalar out into `p.tensor`
// instead of in. Coordinates must land in the same per-axis order the
// tensor's shape/indexNames were built in by generateParameters, which is
// why leftCoords is accumulated with push_front/pop_front here versus
// push_back/pop_back in seedAD — the two walk the same recursion shape but
// feed a differently-shaped (mirror-image) coordinate list into the shared
// getIndex() lookup.
template<typename ADNested, typename T>
void extractAD(ADNested& y, Param<T>& p, size_t curOrder, const std::string& sequence,
    std::deque<size_t>& leftCoords, std::deque<size_t>& rightCoords, size_t& primalIndex) {
    if constexpr(!isADNestedType<ADNested>) {
        if (curOrder > sequence.length()) {
            std::deque<size_t> finalCoords = leftCoords;
            std::reverse(finalCoords.begin(), finalCoords.end());
            finalCoords.push_back(primalIndex);
            finalCoords.insert(finalCoords.end(), rightCoords.begin(), rightCoords.end());
            const size_t flatIndex = p.tensor.getIndex(finalCoords);
            if (flatIndex >= p.tensor.data.size()) {
                throw std::out_of_range("Extract coordinate exceeds tensor storage size.");
            }
            p.tensor.data[flatIndex] = y;
            return;
        }
    }

    size_t orderLabel = sequence.length() - curOrder + 1;
    char curType = sequence[orderLabel - 1];
    if (p.isActive(orderLabel)) {
        size_t curShape = p.getShapeAt(orderLabel);
        if (curType == 't') {
            // COMPILER GUARD: Only compile this loop if y actually has a .tangent() method
            if constexpr(requires { y.tangent(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    rightCoords.push_back(i);
                    extractAD(y.tangent(i), p, curOrder + 1, sequence,
                    leftCoords, rightCoords, primalIndex);
                    rightCoords.pop_back();
                }
            }
            
        } else if (curType == 'a') {
            // COMPILER GUARD: Only compile this loop if y actually has a .adjoint() method
            if constexpr(requires { y.adjoint(0); }) {
                for (size_t i = 0; i < curShape; i++) {
                    leftCoords.push_front(i);
                    extractAD(y.adjoint(i), p, curOrder + 1, sequence, 
                    leftCoords, rightCoords, primalIndex);
                    leftCoords.pop_front();
                }
            }
        }
    } else {
        if constexpr(requires { y.value(); })
            extractAD(y.value(), p, curOrder + 1, sequence, 
            leftCoords, rightCoords, primalIndex);
    }
}

// Seeds every Param "owned" by `targetOrder`'s generated driver layer (see
// generateTangent/generateAdjoint in generator.cpp, which call this once
// per layer with their own order number). A Param can only be seeded when
// its data is actually needed at that point in the call stack: an adjoint
// order's seed is the incoming sensitivity, only meaningful once that
// order's tape exists (after f() has run for that layer) — and if a Param
// is active at multiple adjoint orders, it must wait for the *outermost*
// one, since inner adjoint tapes no longer exist by the time the outer
// tape is being seeded. `ownedHere` encodes exactly that: for an adjoint
// target order, this Param belongs here iff it's active at this order and
// has no later (outer) adjoint order still to come; for a tangent target
// order, it belongs here iff this is literally its highest order and it
// has no adjoint order at all (tangent seeds are needed up front, before
// f() runs, so they're owned by the single non-adjoint layer that created
// them).
template<typename ADNestedX, typename ADNestedY, typename T>
void seedADForOrder(ADNestedX& x, ADNestedY& y, std::vector<Param<T>>& parameters, size_t targetOrder, const std::string& sequence) {
    for (Param<T>& p : parameters) {
        bool hasLaterAdjoint = std::any_of(p.activeOrders.begin(), p.activeOrders.end(),
            [&](int o) { return sequence[o - 1] == 'a' && (size_t)o > targetOrder; });
        bool hasAnyAdjoint = std::any_of(p.activeOrders.begin(), p.activeOrders.end(),
            [&](int o) { return sequence[o - 1] == 'a'; });
        bool ownedHere = sequence[targetOrder - 1] == 'a'
            ? p.isActive(targetOrder) && !hasLaterAdjoint
            : p.highestOrder == targetOrder && !hasAnyAdjoint;
        if (ownedHere && p.role == ParamRole::Input) {
            if (p.name[0] == 'X') {
                for (size_t i = 0; i < x.size(); i++) {
                    std::deque<size_t> leftCoords = {};
                    std::deque<size_t> rightCoords = {};
                    seedAD(x[i], p, 1, sequence,
                    leftCoords, rightCoords, i);
                }
            }
            else if (p.name[0] == 'Y' || p.name[0] == 'F') {
                for (size_t j = 0; j < y.size(); j++) {
                    std::deque<size_t> leftCoords = {};
                    std::deque<size_t> rightCoords = {};
                    seedAD(y[j], p, 1, sequence, 
                    leftCoords, rightCoords, j);
                }
            }
        }
    }
}

// Read-back counterpart of seedADForOrder: same "ownedHere" rule, just
// applied to Output-role Params (a layer's results) instead of Input-role
// ones (a layer's seeds).
template<typename ADNestedX, typename ADNestedY, typename T>
void extractADForOrder(ADNestedX& x, ADNestedY& y, std::vector<Param<T>>& parameters, size_t targetOrder, const std::string& sequence) {
    for (Param<T>& p : parameters) {
        bool hasLaterAdjoint = std::any_of(p.activeOrders.begin(), p.activeOrders.end(),
            [&](int o) { return sequence[o - 1] == 'a' && (size_t)o > targetOrder; });
        bool hasAnyAdjoint = std::any_of(p.activeOrders.begin(), p.activeOrders.end(),
            [&](int o) { return sequence[o - 1] == 'a'; });
        bool ownedHere = sequence[targetOrder - 1] == 'a'
            ? p.isActive(targetOrder) && !hasLaterAdjoint
            : p.highestOrder == targetOrder && !hasAnyAdjoint;
        if (ownedHere && p.role == ParamRole::Output) {
            if (p.name[0] == 'X') {
                for (size_t i = 0; i < x.size(); i++) {
                    std::deque<size_t> leftCoords = {};
                    std::deque<size_t> rightCoords = {};
                    extractAD(x[i], p, 1, sequence, 
                    leftCoords, rightCoords, i);
                }
            }
            else if (p.name[0] == 'Y' || p.name[0] == 'F') {
                for (size_t j = 0; j < y.size(); j++) {
                    std::deque<size_t> leftCoords = {};
                    std::deque<size_t> rightCoords = {};
                    extractAD(y[j], p, 1, sequence, 
                    leftCoords, rightCoords, j);
                }
            }
        }
    }
}

// Seed the primal value for the nested type
template<typename ADNested, typename T>
void seedPrimal(ADNested& x, const std::vector<Param<T>>& parameters) {
    const size_t xSize = x.size();

    for (const auto& p : parameters) {
        if (p.name.size() == 1 && p.name[0] == 'X') {
            for (size_t i = 0; i < xSize; ++i) {
                seedElement(x[i], p, i);
            }
        }
    }
}

// Extract the primal value for the nested type
template<typename ADNested, typename T>
void extractPrimal(ADNested& y, std::vector<Param<T>>& parameters) {
    const size_t ySize = y.size();

    for (Param<T>& p : parameters) {
        if (p.name.length() == 1 && p.name[0] == 'Y') {
            for (size_t j = 0; j < ySize; j++) {
                extractElement(y[j], p, j);
            }
        }
    }
}

// Helper functions for seeding and extracting
template<typename ADNested, typename T>
void seedElement(ADNested& element, const Param<T> &p, size_t index) {
    if constexpr(isADNestedType<ADNested>) {
        seedElement(element.value(), p, index);
    } else {
        element = p.tensor.data[p.tensor.getIndex({index})];
    }
}


template<typename ADNested, typename T>
void extractElement(ADNested& element, Param<T> &p, size_t index) {
    if constexpr(isADNestedType<ADNested>) {
        extractElement(element.value(), p, index);
    } else {
        p.tensor.data[p.tensor.getIndex({index})] = element;
    }
}


// Find a specific parameter by name from a list of parameters
// Handle std::deque<Param<T>> and const std::deque<Param<T>>
template <typename T>
Param<T>& findParamByName(const std::string& targetName, std::deque<Param<T>>& parameters) {
    for (auto& p : parameters) {
        if (p.name == targetName) {
            return p;
        }
    }
    throw std::runtime_error(std::format("Parameter not found: {}", targetName));
}

template <typename T>
Param<T>& findParamByName(const std::string& targetName, const std::deque<Param<T>>& parameters) {
    for (auto& p : parameters) {
        if (p.name == targetName) {
            return p;
        }
    }
    throw std::runtime_error(std::format("Parameter not found: {}", targetName));
}

#endif
