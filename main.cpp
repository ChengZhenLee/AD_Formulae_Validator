#include "1/va.h"
#include "1/vt.h"
#include <iostream>


int main() {
  std::cout << "=== Testing First Derivative (Adjoint Mode) ===\n";
  bool va_pass = Validate_va<double>();
  
  std::cout << "\n=== Testing First Derivative (Tangent Mode) ===\n";
  bool vt_pass = Validate_vt<double>();
  
  std::cout << "\n=== Summary ===\n";
  std::cout << "Adjoint validation: " << (va_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent validation: " << (vt_pass ? "PASSED" : "FAILED") << "\n";

  return 0;
}
