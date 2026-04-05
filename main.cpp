#include "1/va.h"
#include "1/vt.h"
#include "2/vtvt.h"
#include <iostream>
#include <ctime>


int main() {
  // Seed random
  std::srand(static_cast<unsigned int>(std::time(nullptr)));

  std::cout << "=== Testing First Derivative (Adjoint Mode) ===\n";
  bool va_pass = Validate_va<double, 2>();
  
  std::cout << "\n=== Testing First Derivative (Tangent Mode) ===\n";
  bool vt_pass = Validate_vt<double, 3>();

  std::cout << "\n=== Testing Second Derivative (Tangent over Tangent Mode) ===\n";
  bool vtvt_pass = Validate_vtvt<double, 3, 3>();
  
  std::cout << "\n=== Summary ===\n";
  std::cout << "Adjoint validation: " << (va_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent validation: " << (vt_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent over Tangent validation: " << (vtvt_pass ? "PASSED" : "FAILED") << "\n";

  return 0;
}
