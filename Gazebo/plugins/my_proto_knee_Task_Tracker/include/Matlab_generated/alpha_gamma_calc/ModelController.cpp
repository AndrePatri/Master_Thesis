//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 00:05:51
//

// Include Files
#include "ModelController.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include <algorithm>
#include <cmath>

// Function Declarations
static double rt_atan2d_snf(double u0, double u1);

// Function Definitions
//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int b_u0;
    int b_u1;
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }
    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }
    y = std::atan2(static_cast<double>(b_u0), static_cast<double>(b_u1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }
  return y;
}

//
// Arguments    : void
// Return Type  : void
//
ModelController::ModelController()
{
  static const double dv[18]{0.1,
                             0.145,
                             0.236,
                             0.256,
                             0.355,
                             -1.0,
                             0.0,
                             0.78539816339744828,
                             2.5132741228718345,
                             0.53376787572416418,
                             0.048333333333333332,
                             0.58800260354756761,
                             -0.58800260354756761,
                             -2.5535900500422253,
                             2.5535900500422253,
                             1.2566370614359172,
                             0.21749999999999997,
                             0.145};
  std::copy(&dv[0], &dv[18], &this->robot_dim[0]);
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
// DESCRIPTION:
//  computes, given a linkage beta, the corresponding alpha and gamma,
//  employing proper kinematic relationships and controls
//  PARAMETERS:
//  beta-->linkage angle (see the system's drawing for clarifications)
//
// Arguments    : double beta
//                double *alpha
//                double *b_gamma
// Return Type  : void
//
void ModelController::alpha_gamma_calc(double beta, double *alpha,
                                       double *b_gamma)
{
  double A;
  double A_tmp;
  double B;
  double B_tmp;
  double cos_gamma;
  double sin_gamma;
  A_tmp = std::cos(beta);
  A = ((((this->robot_dim[0] * this->robot_dim[0] -
          this->robot_dim[1] * this->robot_dim[1]) -
         this->robot_dim[2] * this->robot_dim[2]) -
        this->robot_dim[3] * this->robot_dim[3]) +
       2.0 * this->robot_dim[1] * this->robot_dim[3] * A_tmp) /
      (2.0 * this->robot_dim[2] *
       (this->robot_dim[3] * A_tmp - this->robot_dim[1]));
  B_tmp = std::sin(beta);
  B = this->robot_dim[3] * B_tmp /
      (this->robot_dim[3] * std::cos(beta) - this->robot_dim[1]);
  sin_gamma = B * B + 1.0;
  cos_gamma = (-A * B + std::sqrt(sin_gamma - A * A)) / sin_gamma;
  sin_gamma = A + B * cos_gamma;
  *b_gamma = rt_atan2d_snf(sin_gamma, cos_gamma);
  A = this->robot_dim[3] / this->robot_dim[0];
  B = this->robot_dim[2] / this->robot_dim[0];
  *alpha = rt_atan2d_snf(A * B_tmp - B * cos_gamma,
                         (this->robot_dim[1] / this->robot_dim[0] - A * A_tmp) -
                             B * sin_gamma);
}

//
// File trailer for ModelController.cpp
//
// [EOF]
//
