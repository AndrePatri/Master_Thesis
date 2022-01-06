//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 17:55:01
//

// Include Files
#include "ModelController.h"

// Function Definitions
//
// B_CHECK_INV_CALC
//     B_CHECK_INV = B_CHECK_INV_CALC(H1)
//
// Arguments    : double h1
//                double B_check_inv[4]
// Return Type  : void
//
void ModelController::B_check_inv_calc(double h1, double B_check_inv[4])
{
  //     This function was generated by the Symbolic Math Toolbox version 8.7.
  //     29-Jul-2021 16:16:22
  B_check_inv[0] = 0.5;
  B_check_inv[1] = h1 / 2.0;
  B_check_inv[2] = 0.0;
  B_check_inv[3] = 0.5;
}

//
// Arguments    : void
// Return Type  : void
//
ModelController::ModelController()
{
}

//
// Arguments    : void
// Return Type  : void
//
ModelController::~ModelController()
{
  // (no terminate code required)
}

//
// File trailer for ModelController.cpp
//
// [EOF]
//
