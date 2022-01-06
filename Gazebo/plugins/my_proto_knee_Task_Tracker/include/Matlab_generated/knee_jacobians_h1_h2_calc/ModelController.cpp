//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 01:07:28
//

// Include Files
#include "ModelController.h"
#include "rt_nonfinite.h"
#include <algorithm>
#include <cmath>

// Function Declarations
static double rt_powd_snf(double u0, double u1);

// Function Definitions
//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = std::abs(u0);
    d1 = std::abs(u1);
    if (std::isinf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = std::pow(u0, u1);
    }
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
// Arguments    : double beta
//                double knee_jacobians[2]
// Return Type  : void
//
void ModelController::knee_jacobians_h1_h2_calc(double beta,
                                                double knee_jacobians[2])
{
  double a_tmp;
  double a_tmp_tmp;
  double b_a_tmp;
  double b_a_tmp_tmp;
  double b_knee_jacobians_tmp;
  double c_a_tmp_tmp;
  double c_knee_jacobians_tmp;
  double d_a_tmp_tmp;
  double d_knee_jacobians_tmp;
  double e_a_tmp_tmp;
  double e_knee_jacobians_tmp;
  double f_a_tmp_tmp;
  double f_knee_jacobians_tmp;
  double g_a_tmp_tmp;
  double g_knee_jacobians_tmp;
  double h1_tmp_tmp;
  double h1_tmp_tmp_tmp;
  double h_a_tmp_tmp;
  double h_knee_jacobians_tmp;
  double i_a_tmp_tmp;
  double i_knee_jacobians_tmp;
  double j_a_tmp_tmp;
  double j_knee_jacobians_tmp;
  double k_a_tmp_tmp;
  double k_knee_jacobians_tmp;
  double knee_jacobians_tmp;
  double knee_jacobians_tmp_tmp_tmp;
  double l_a_tmp_tmp;
  h1_tmp_tmp = std::cos(beta);
  h1_tmp_tmp_tmp = std::sin(beta);
  a_tmp_tmp = this->robot_dim[3] * h1_tmp_tmp;
  a_tmp = this->robot_dim[1] - a_tmp_tmp;
  b_a_tmp_tmp = this->robot_dim[3] * this->robot_dim[3];
  c_a_tmp_tmp = 2.0 * h1_tmp_tmp;
  d_a_tmp_tmp = this->robot_dim[0] * this->robot_dim[0];
  e_a_tmp_tmp = this->robot_dim[2] * this->robot_dim[2];
  f_a_tmp_tmp = c_a_tmp_tmp * this->robot_dim[1];
  g_a_tmp_tmp = this->robot_dim[1] * this->robot_dim[1];
  h_a_tmp_tmp = f_a_tmp_tmp * this->robot_dim[3];
  b_a_tmp = (((-d_a_tmp_tmp + g_a_tmp_tmp) - h_a_tmp_tmp) + e_a_tmp_tmp) +
            b_a_tmp_tmp;
  i_a_tmp_tmp = a_tmp * a_tmp;
  j_a_tmp_tmp = h1_tmp_tmp_tmp * h1_tmp_tmp_tmp;
  k_a_tmp_tmp = b_a_tmp_tmp * j_a_tmp_tmp;
  l_a_tmp_tmp = k_a_tmp_tmp / i_a_tmp_tmp;
  knee_jacobians_tmp = 2.0 * d_a_tmp_tmp;
  b_knee_jacobians_tmp = rt_powd_snf(a_tmp, 3.0);
  c_knee_jacobians_tmp = 4.0 * e_a_tmp_tmp;
  h_a_tmp_tmp = (g_a_tmp_tmp - h_a_tmp_tmp) + b_a_tmp_tmp;
  d_knee_jacobians_tmp = 2.0 * b_a_tmp_tmp;
  e_knee_jacobians_tmp = this->robot_dim[1] * b_a_tmp_tmp;
  f_knee_jacobians_tmp = this->robot_dim[2] * b_knee_jacobians_tmp;
  g_knee_jacobians_tmp = c_knee_jacobians_tmp * b_knee_jacobians_tmp;
  knee_jacobians_tmp_tmp_tmp = this->robot_dim[3] * h1_tmp_tmp_tmp;
  h_knee_jacobians_tmp = this->robot_dim[2] * i_a_tmp_tmp;
  i_knee_jacobians_tmp = 2.0 * this->robot_dim[2] * i_a_tmp_tmp;
  j_knee_jacobians_tmp = std::sqrt(
      (l_a_tmp_tmp - b_a_tmp * b_a_tmp / (c_knee_jacobians_tmp * i_a_tmp_tmp)) +
      1.0);
  k_knee_jacobians_tmp = j_knee_jacobians_tmp + knee_jacobians_tmp_tmp_tmp *
                                                    b_a_tmp /
                                                    i_knee_jacobians_tmp;
  knee_jacobians[0] =
      -(2.0 * this->robot_dim[0] *
        ((h_knee_jacobians_tmp *
              (((a_tmp_tmp * b_a_tmp / i_knee_jacobians_tmp +
                 e_knee_jacobians_tmp * j_a_tmp_tmp / h_knee_jacobians_tmp) -
                k_a_tmp_tmp * b_a_tmp / f_knee_jacobians_tmp) +
               knee_jacobians_tmp_tmp_tmp *
                   ((((((((((rt_powd_snf(this->robot_dim[0], 4.0) +
                             c_a_tmp_tmp * d_a_tmp_tmp * this->robot_dim[1] *
                                 this->robot_dim[3]) -
                            knee_jacobians_tmp * e_a_tmp_tmp) -
                           knee_jacobians_tmp * b_a_tmp_tmp) -
                          rt_powd_snf(this->robot_dim[1], 4.0)) +
                         c_a_tmp_tmp * rt_powd_snf(this->robot_dim[1], 3.0) *
                             this->robot_dim[3]) +
                        f_a_tmp_tmp * e_a_tmp_tmp * this->robot_dim[3]) -
                       f_a_tmp_tmp * rt_powd_snf(this->robot_dim[3], 3.0)) +
                      rt_powd_snf(this->robot_dim[2], 4.0)) -
                     2.0 * e_a_tmp_tmp * b_a_tmp_tmp) +
                    rt_powd_snf(this->robot_dim[3], 4.0)) /
                   (g_knee_jacobians_tmp * j_knee_jacobians_tmp)) /
              (this->robot_dim[0] * h_a_tmp_tmp) -
          a_tmp_tmp / this->robot_dim[0]) +
         2.0 * this->robot_dim[2] * b_a_tmp_tmp * h1_tmp_tmp_tmp *
             (this->robot_dim[3] - this->robot_dim[1] * h1_tmp_tmp) *
             k_knee_jacobians_tmp /
             (this->robot_dim[0] * b_knee_jacobians_tmp *
              ((l_a_tmp_tmp + 1.0) * (l_a_tmp_tmp + 1.0)))) *
        h_a_tmp_tmp) /
      (a_tmp * (((((((d_knee_jacobians_tmp * (h1_tmp_tmp * h1_tmp_tmp) +
                      d_knee_jacobians_tmp * j_a_tmp_tmp) +
                     d_a_tmp_tmp) +
                    g_a_tmp_tmp) -
                   e_a_tmp_tmp) -
                  b_a_tmp_tmp) -
                 2.0 * this->robot_dim[1] * this->robot_dim[3] * h1_tmp_tmp) +
                2.0 * this->robot_dim[2] * this->robot_dim[3] * h1_tmp_tmp_tmp *
                    j_knee_jacobians_tmp));
  knee_jacobians[1] =
      -(i_a_tmp_tmp *
            (((this->robot_dim[3] * std::cos(beta) *
                   ((((-(this->robot_dim[0] * this->robot_dim[0]) +
                       this->robot_dim[1] * this->robot_dim[1]) -
                      2.0 * std::cos(beta) * this->robot_dim[1] *
                          this->robot_dim[3]) +
                     this->robot_dim[2] * this->robot_dim[2]) +
                    this->robot_dim[3] * this->robot_dim[3]) /
                   (2.0 * this->robot_dim[2] * (a_tmp * a_tmp)) +
               e_knee_jacobians_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                   (this->robot_dim[2] * (a_tmp * a_tmp))) -
              b_a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) * b_a_tmp /
                  f_knee_jacobians_tmp) +
             this->robot_dim[3] * std::sin(beta) *
                 ((((((((((rt_powd_snf(this->robot_dim[0], 4.0) +
                           2.0 * std::cos(beta) *
                               (this->robot_dim[0] * this->robot_dim[0]) *
                               this->robot_dim[1] * this->robot_dim[3]) -
                          2.0 * (this->robot_dim[0] * this->robot_dim[0]) *
                              (this->robot_dim[2] * this->robot_dim[2])) -
                         2.0 * (this->robot_dim[0] * this->robot_dim[0]) *
                             (this->robot_dim[3] * this->robot_dim[3])) -
                        rt_powd_snf(this->robot_dim[1], 4.0)) +
                       2.0 * std::cos(beta) *
                           rt_powd_snf(this->robot_dim[1], 3.0) *
                           this->robot_dim[3]) +
                      2.0 * std::cos(beta) * this->robot_dim[1] *
                          (this->robot_dim[2] * this->robot_dim[2]) *
                          this->robot_dim[3]) -
                     2.0 * std::cos(beta) * this->robot_dim[1] *
                         rt_powd_snf(this->robot_dim[3], 3.0)) +
                    rt_powd_snf(this->robot_dim[2], 4.0)) -
                   2.0 * (this->robot_dim[2] * this->robot_dim[2]) *
                       (this->robot_dim[3] * this->robot_dim[3])) +
                  rt_powd_snf(this->robot_dim[3], 4.0)) /
                 (g_knee_jacobians_tmp *
                  std::sqrt((b_a_tmp_tmp * (h1_tmp_tmp_tmp * h1_tmp_tmp_tmp) /
                                 (a_tmp * a_tmp) -
                             b_a_tmp * b_a_tmp /
                                 (c_knee_jacobians_tmp * (a_tmp * a_tmp))) +
                            1.0))) /
            h_a_tmp_tmp +
        d_knee_jacobians_tmp * h1_tmp_tmp_tmp *
            (this->robot_dim[3] - this->robot_dim[1] * std::cos(beta)) *
            k_knee_jacobians_tmp /
            (b_knee_jacobians_tmp *
             ((l_a_tmp_tmp + 1.0) * (l_a_tmp_tmp + 1.0)))) /
      (b_a_tmp / (2.0 * this->robot_dim[2] * a_tmp) -
       knee_jacobians_tmp_tmp_tmp * a_tmp * k_knee_jacobians_tmp / h_a_tmp_tmp);
}

//
// File trailer for ModelController.cpp
//
// [EOF]
//
