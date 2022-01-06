//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 15:29:23
//

// Include Files
#include "ModelController.h"
#include <cstring>

// Function Definitions
//
// A_STATE_LIM_CALC
//     A_STATE_LIM = A_STATE_LIM_CALC(DELTA)
//
// Arguments    : double delta
//                double A_state_lim[28]
// Return Type  : void
//
void ModelController::A_state_lim_calc(double delta, double A_state_lim[28])
{
  double t3;
  //     This function was generated by the Symbolic Math Toolbox version 8.7.
  //     29-Jul-2021 16:01:06
  t3 = delta * delta / 2.0;
  std::memset(&A_state_lim[0], 0, 12U * sizeof(double));
  A_state_lim[12] = t3;
  A_state_lim[13] = 0.0;
  A_state_lim[14] = -t3;
  A_state_lim[15] = 0.0;
  A_state_lim[16] = 0.0;
  A_state_lim[17] = t3;
  A_state_lim[18] = 0.0;
  A_state_lim[19] = -t3;
  std::memset(&A_state_lim[20], 0, 8U * sizeof(double));
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
