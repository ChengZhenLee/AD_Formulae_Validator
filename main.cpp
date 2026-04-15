#include "1/va.h"
#include "1/vt.h"
#include "2/vtvt.h"
#include "2/vtva.h"
#include "2/vavt.h"
#include "2/vava.h"
#include "3/vtvtvt.h"
#include "3/vtvtva.h"
#include "3/vastsa.h"
#include "4/vtvtvtvt.h"
#include <iostream>
#include <ctime>


int main() {
  // Seed random
  std::srand(static_cast<unsigned int>(std::time(nullptr)));

  // Prepare log file
  std::ofstream logFile("validation_results.txt");

  if (!logFile.is_open()) {
    std::cout << "Could not open log file";
    return 1;
  }

  bool va_pass = Validate_va<double, 2>(logFile);
  
  bool vt_pass = Validate_vt<double, 3>(logFile);

  bool vtvt_pass = Validate_vtvt<double, 3, 3>(logFile);

  bool vtva_pass = Validate_vtva<double, 2, 3>(logFile);

  bool vavt_pass = Validate_vavt<double, 3, 2>(logFile);

  bool vava_pass = Validate_vava<double, 2, 3>(logFile);

  bool vtvtvt_pass = Validate_vtvtvt<double, 3, 3, 3>(logFile);

  bool vtvtva_pass = Validate_vtvtva<double, 2, 3, 3>(logFile);

  bool vastsa_pass = Validate_vastsa<double, 2>(logFile);

  bool vtvtvtvt_pass = Validate_vtvtvtvt<double, 3, 3, 3, 3>(logFile);

  logFile.close();
  
  std::cout << "\n=== Summary ===\n";
  std::cout << "Adjoint validation: " << (va_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent validation: " << (vt_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent over Tangent validation: " << (vtvt_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent over Adjoint validation: " << (vtva_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Adjoint over Tangent validation: " << (vavt_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Adjoint over Adjoint validation: " << (vava_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent over Tangent over Tangent validation: " << (vtvtvt_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent over Tangent over Adjoint validation: " << (vtvtva_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Adjoint over Tangent over Adjoint validation: " << (vastsa_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "Tangent over Tangent over Tangent over Tangent validation: " << (vtvtvtvt_pass ? "PASSED" : "FAILED") << "\n";
  std::cout << "You can find the log file at validations_results.txt\n";

  return 0;
}
