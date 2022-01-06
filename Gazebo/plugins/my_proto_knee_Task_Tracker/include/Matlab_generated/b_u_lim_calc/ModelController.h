//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ModelController.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 06-Aug-2021 10:29:52
//

#ifndef MODELCONTROLLER_H
#define MODELCONTROLLER_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
class ModelController {
public:
  ModelController();
  ~ModelController();
  void b_u_lim_calc(double alpha, double b_gamma, double h1, double h2,
                    const double u_lim[4], double phi_r, double beta,
                    double phi_r_dot, double beta_dot, double b_u_lim[8]);
  double robot_knee_el_par[7];
  double robot_dim[7];
  double robot_inertial_par[15];
  double g;
};

#endif
//
// File trailer for ModelController.h
//
// [EOF]
//
