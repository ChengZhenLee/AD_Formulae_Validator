#include <iostream>
#include <vector>
#include <deque>
#include "structures.h"

using namespace std;

int main() {
    using T = double;
    // Test productShape simple merge of (2,3) * (3,4) -> (2,4)
    deque<size_t> a = {2,3};
    deque<size_t> b = {3,4};
    auto r = Tensor<T>::productShape(a,b);
    if (r.size() != 2 || r[0] != 2 || r[1] != 4) {
        cerr << "productShape test failed" << endl;
        return 1;
    }

    // Test validateTensorChain: chain [(2,3),(3,4),(4,5)] should be valid
    Tensor<T> A({2,3});
    Tensor<T> B({3,4});
    Tensor<T> C({4,5});
    deque<Tensor<T>> chain = {A,B,C};
    if (!Tensor<T>::validateTensorChain(chain)) {
        cerr << "validateTensorChain test failed" << endl;
        return 2;
    }

    // Test product numeric: multiply A*B -> shape (2,4)
    for (size_t i=0;i<A.data.size();++i) A.data[i] = 1.0;
    for (size_t i=0;i<B.data.size();++i) B.data[i] = 2.0;
    Tensor<T> P = Tensor<T>::product(A,B);
    if (P.shape.size()!=2 || P.shape[0]!=2 || P.shape[1]!=4) {
        cerr << "numeric product shape mismatch" << endl;
        return 3;
    }
    // numeric value check: each entry should be sum_{k=0..2} 1*2 = 6
    for (size_t i = 0; i < P.data.size(); ++i) {
        if (std::abs(P.data[i] - 6.0) > 1e-12) {
            cerr << "numeric product value mismatch at idx " << i << ": " << P.data[i] << endl;
            return 4;
        }
    }

    // Test invalid chain detection
    Tensor<T> D({2,3});
    Tensor<T> E({4,5});
    deque<Tensor<T>> badChain = {D,E};
    if (Tensor<T>::validateTensorChain(badChain)) {
        cerr << "validateTensorChain false-negative" << endl;
        return 5;
    }

    // Test scalar/empty-shape handling: productShape({}, {3,4}) -> {3,4}
    deque<size_t> emptyShape = {};
    deque<size_t> v = {3,4};
    auto r2 = Tensor<T>::productShape(emptyShape, v);
    if (r2 != v) {
        cerr << "productShape empty-shape failed" << endl;
        return 6;
    }

    // Test triple chain numeric: (2,3)*(3,4)*(4,2) -> (2,2)
    Tensor<T> C2({4,2});
    for (size_t i=0;i<C2.data.size();++i) C2.data[i] = 3.0; // fill with 3
    deque<Tensor<T>> chain3 = {A,B,C2};
    Tensor<T> R3 = Tensor<T>::product(chain3);
    if (R3.shape.size()!=2 || R3.shape[0]!=2 || R3.shape[1]!=2) {
        cerr << "triple product shape mismatch" << endl;
        return 7;
    }

    cout << "All generator tests passed" << endl;
    return 0;
}
