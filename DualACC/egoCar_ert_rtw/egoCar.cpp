//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: egoCar.cpp
//
// Code generated for Simulink model 'egoCar'.
//
// Model version                  : 1.6
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Wed Jun 18 12:33:19 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Apple->ARM64
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "egoCar.h"
#include "rtwtypes.h"
#include <cstring>
#include <cmath>
#include "egoCar_private.h"
#include "cmath"

// Named constants for MATLAB Function: '<S34>/optimizer'
const real_T egoCar_RMDscale{ 0.02 };

const real_T egoCar_RMVscale{ 0.2 };

const int32_T egoCar_degrees{ 4 };

const int32_T egoCar_ny{ 2 };

const int32_T egoCar_p{ 30 };

const real_T egoCar_voff{ 0.4 };

//
// This function updates continuous states using the ODE3 fixed-step
// solver algorithm
//
void egoCar::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  // Solver Matrices
  static const real_T rt_ODE3_A[3]{
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3]{
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t { rtsiGetT(si) };

  time_T tnew { rtsiGetSolverStopTime(si) };

  time_T h { rtsiGetStepSize(si) };

  real_T *x { rtsiGetContStates(si) };

  ODE3_IntgData *id { static_cast<ODE3_IntgData *>(rtsiGetSolverData(si)) };

  real_T *y { id->y };

  real_T *f0 { id->f[0] };

  real_T *f1 { id->f[1] };

  real_T *f2 { id->f[2] };

  real_T hB[3];
  int_T i;
  int_T nXc { 3 };

  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  // Save the state values at time t in y, we'll use x as ynew.
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  // Assumes that rtsiSetT and ModelOutputs are up-to-date
  // f0 = f(t,y)
  rtsiSetdX(si, f0);
  egoCar_derivatives();

  // f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*));
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  this->step();
  egoCar_derivatives();

  // f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*));
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  this->step();
  egoCar_derivatives();

  // tnew = t + hA(3);
  // ynew = y + f*hB(:,3);
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

//
// Output and update for atomic system:
//    '<S1>/DataTypeConversion_L0'
//    '<S1>/DataTypeConversion_amax'
//    '<S1>/DataTypeConversion_amin'
//    '<S1>/DataTypeConversion_atrack'
//    '<S1>/DataTypeConversion_vset'
//
void egoCar::egoCar_DataTypeConversion_L0(real_T rtu_u, real_T *rty_y)
{
  *rty_y = rtu_u;
}

// Function for MATLAB Function: '<S34>/optimizer'
real_T egoCar::egoCar_norm(const real_T x[4])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[3]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

// Function for MATLAB Function: '<S34>/optimizer'
real_T egoCar::egoCar_maximum(const real_T x[4])
{
  real_T ex;
  int32_T idx;
  int32_T k;
  if (!std::isnan(x[0])) {
    idx = 1;
  } else {
    boolean_T exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 5)) {
      if (!std::isnan(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    for (k = idx + 1; k < 5; k++) {
      real_T x_0;
      x_0 = x[k - 1];
      if (ex < x_0) {
        ex = x_0;
      }
    }
  }

  return ex;
}

// Function for MATLAB Function: '<S34>/optimizer'
real_T egoCar::egoCar_xnrm2(int32_T n, const real_T x[16], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (int32_T k{ix0}; k < kend; k++) {
        real_T absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = std::sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = std::sqrt(b * b + 1.0) * a;
  } else if (std::isnan(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

// Function for MATLAB Function: '<S34>/optimizer'
void egoCar::egoCar_xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T
  ia0, const real_T x[16], int32_T ix0, real_T y[4])
{
  if ((b_m != 0) && (n != 0)) {
    int32_T b;
    if (n - 1 >= 0) {
      std::memset(&y[0], 0, static_cast<uint32_T>(n) * sizeof(real_T));
    }

    b = ((n - 1) << 2) + ia0;
    for (int32_T b_iy{ia0}; b_iy <= b; b_iy += 4) {
      real_T c;
      int32_T d;
      int32_T ia;
      c = 0.0;
      d = (b_iy + b_m) - 1;
      for (ia = b_iy; ia <= d; ia++) {
        c += x[((ix0 + ia) - b_iy) - 1] * b_A[ia - 1];
      }

      ia = (b_iy - ia0) >> 2;
      y[ia] += c;
    }
  }
}

// Function for MATLAB Function: '<S34>/optimizer'
void egoCar::egoCar_xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0,
  const real_T y[4], real_T b_A[16], int32_T ia0)
{
  if (!(alpha1 == 0.0)) {
    int32_T jA;
    jA = ia0;
    for (int32_T j{0}; j < n; j++) {
      real_T temp;
      temp = y[j];
      if (temp != 0.0) {
        int32_T b;
        temp *= alpha1;
        b = (b_m + jA) - 1;
        for (int32_T ijA{jA}; ijA <= b; ijA++) {
          b_A[ijA - 1] += b_A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }

      jA += 4;
    }
  }
}

// Function for MATLAB Function: '<S34>/optimizer'
void egoCar::egoCar_KWIKfactor(const real_T b_Ac[384], const int32_T iC[96],
  int32_T nA, const real_T b_Linv[16], real_T D[16], real_T b_H[16], int32_T n,
  real_T RLinv[16], real_T *Status)
{
  real_T Q[16];
  real_T R[16];
  real_T TL[16];
  real_T b_A[16];
  real_T tau[4];
  real_T work[4];
  real_T atmp;
  real_T b_A_0;
  real_T beta1;
  int32_T b_coltop;
  int32_T b_lastv;
  int32_T e_tmp;
  int32_T exitg1;
  int32_T ii;
  int32_T k_i;
  int32_T knt;
  boolean_T exitg2;
  *Status = 1.0;
  std::memset(&RLinv[0], 0, sizeof(real_T) << 4U);
  for (ii = 0; ii < nA; ii++) {
    b_lastv = iC[ii];
    for (k_i = 0; k_i < 4; k_i++) {
      RLinv[k_i + (ii << 2)] = ((b_Ac[b_lastv - 1] * b_Linv[k_i] + b_Linv[k_i +
        4] * b_Ac[b_lastv + 95]) + b_Linv[k_i + 8] * b_Ac[b_lastv + 191]) +
        b_Linv[k_i + 12] * b_Ac[b_lastv + 287];
    }
  }

  std::memcpy(&b_A[0], &RLinv[0], sizeof(real_T) << 4U);
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  tau[2] = 0.0;
  work[2] = 0.0;
  tau[3] = 0.0;
  work[3] = 0.0;
  for (k_i = 0; k_i < 4; k_i++) {
    ii = (k_i << 2) + k_i;
    if (k_i + 1 < 4) {
      atmp = b_A[ii];
      b_lastv = ii + 2;
      tau[k_i] = 0.0;
      beta1 = egoCar_xnrm2(3 - k_i, b_A, ii + 2);
      if (beta1 != 0.0) {
        b_A_0 = b_A[ii];
        beta1 = rt_hypotd_snf(b_A_0, beta1);
        if (b_A_0 >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          e_tmp = (ii - k_i) + 4;
          do {
            knt++;
            for (b_coltop = b_lastv; b_coltop <= e_tmp; b_coltop++) {
              b_A[b_coltop - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while ((std::abs(beta1) < 1.0020841800044864E-292) && (knt < 20));

          beta1 = rt_hypotd_snf(atmp, egoCar_xnrm2(3 - k_i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[k_i] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          for (b_coltop = b_lastv; b_coltop <= e_tmp; b_coltop++) {
            b_A[b_coltop - 1] *= atmp;
          }

          for (b_lastv = 0; b_lastv < knt; b_lastv++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[k_i] = (beta1 - b_A_0) / beta1;
          atmp = 1.0 / (b_A_0 - beta1);
          b_coltop = (ii - k_i) + 4;
          for (knt = b_lastv; knt <= b_coltop; knt++) {
            b_A[knt - 1] *= atmp;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 4 - k_i;
        knt = (ii - k_i) + 3;
        while ((b_lastv > 0) && (b_A[knt] == 0.0)) {
          b_lastv--;
          knt--;
        }

        knt = 3 - k_i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          b_coltop = (((knt - 1) << 2) + ii) + 4;
          e_tmp = b_coltop;
          do {
            exitg1 = 0;
            if (e_tmp + 1 <= b_coltop + b_lastv) {
              if (b_A[e_tmp] != 0.0) {
                exitg1 = 1;
              } else {
                e_tmp++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        knt = 0;
      }

      if (b_lastv > 0) {
        egoCar_xgemv(b_lastv, knt, b_A, ii + 5, b_A, ii + 1, work);
        egoCar_xgerc(b_lastv, knt, -tau[k_i], ii + 1, work, b_A, ii + 5);
      }

      b_A[ii] = atmp;
    } else {
      tau[3] = 0.0;
    }
  }

  for (k_i = 0; k_i < 4; k_i++) {
    for (ii = 0; ii <= k_i; ii++) {
      b_lastv = k_i << 2;
      R[ii + b_lastv] = b_A[b_lastv + ii];
    }

    for (ii = k_i + 2; ii < 5; ii++) {
      R[(ii + (k_i << 2)) - 1] = 0.0;
    }

    work[k_i] = 0.0;
  }

  for (k_i = 3; k_i >= 0; k_i--) {
    ii = ((k_i << 2) + k_i) + 5;
    if (k_i + 1 < 4) {
      b_A[ii - 5] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 4 - k_i;
        knt = ii - k_i;
        while ((b_lastv > 0) && (b_A[knt - 2] == 0.0)) {
          b_lastv--;
          knt--;
        }

        knt = 3 - k_i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          b_coltop = ((knt - 1) << 2) + ii;
          e_tmp = b_coltop;
          do {
            exitg1 = 0;
            if (e_tmp <= (b_coltop + b_lastv) - 1) {
              if (b_A[e_tmp - 1] != 0.0) {
                exitg1 = 1;
              } else {
                e_tmp++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        knt = 0;
      }

      if (b_lastv > 0) {
        egoCar_xgemv(b_lastv, knt, b_A, ii, b_A, ii - 4, work);
        egoCar_xgerc(b_lastv, knt, -tau[k_i], ii - 4, work, b_A, ii);
      }

      knt = ii - k_i;
      for (b_lastv = ii - 3; b_lastv < knt; b_lastv++) {
        b_A[b_lastv - 1] *= -tau[k_i];
      }
    }

    b_A[ii - 5] = 1.0 - tau[k_i];
    for (b_lastv = 0; b_lastv < k_i; b_lastv++) {
      b_A[(ii - b_lastv) - 6] = 0.0;
    }
  }

  for (k_i = 0; k_i < 4; k_i++) {
    ii = k_i << 2;
    Q[ii] = b_A[ii];
    Q[ii + 1] = b_A[ii + 1];
    Q[ii + 2] = b_A[ii + 2];
    Q[ii + 3] = b_A[ii + 3];
  }

  k_i = 0;
  do {
    exitg1 = 0;
    if (k_i <= nA - 1) {
      if (std::abs(R[(k_i << 2) + k_i]) < 1.0E-12) {
        *Status = -2.0;
        exitg1 = 1;
      } else {
        k_i++;
      }
    } else {
      for (k_i = 0; k_i < n; k_i++) {
        for (ii = 0; ii < n; ii++) {
          b_lastv = k_i << 2;
          knt = ii << 2;
          TL[k_i + knt] = ((b_Linv[b_lastv + 1] * Q[knt + 1] + b_Linv[b_lastv] *
                            Q[knt]) + b_Linv[b_lastv + 2] * Q[knt + 2]) +
            b_Linv[b_lastv + 3] * Q[knt + 3];
        }
      }

      std::memset(&RLinv[0], 0, sizeof(real_T) << 4U);
      for (k_i = nA; k_i >= 1; k_i--) {
        knt = (k_i - 1) << 2;
        b_coltop = (k_i + knt) - 1;
        RLinv[b_coltop] = 1.0;
        for (ii = k_i; ii <= nA; ii++) {
          e_tmp = (((ii - 1) << 2) + k_i) - 1;
          RLinv[e_tmp] /= R[b_coltop];
        }

        if (k_i > 1) {
          for (ii = 0; ii <= k_i - 2; ii++) {
            for (b_lastv = k_i; b_lastv <= nA; b_lastv++) {
              b_coltop = (b_lastv - 1) << 2;
              e_tmp = b_coltop + ii;
              RLinv[e_tmp] -= RLinv[(b_coltop + k_i) - 1] * R[knt + ii];
            }
          }
        }
      }

      for (b_lastv = 0; b_lastv < n; b_lastv++) {
        for (knt = b_lastv + 1; knt <= n; knt++) {
          k_i = ((knt - 1) << 2) + b_lastv;
          b_H[k_i] = 0.0;
          for (b_coltop = nA + 1; b_coltop <= n; b_coltop++) {
            ii = (b_coltop - 1) << 2;
            b_H[k_i] -= TL[(ii + knt) - 1] * TL[ii + b_lastv];
          }

          b_H[(knt + (b_lastv << 2)) - 1] = b_H[k_i];
        }
      }

      for (b_lastv = 0; b_lastv < nA; b_lastv++) {
        for (knt = 0; knt < n; knt++) {
          k_i = (b_lastv << 2) + knt;
          D[k_i] = 0.0;
          for (b_coltop = b_lastv + 1; b_coltop <= nA; b_coltop++) {
            ii = (b_coltop - 1) << 2;
            D[k_i] += TL[ii + knt] * RLinv[ii + b_lastv];
          }
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S34>/optimizer'
void egoCar::egoCar_DropConstraint(int32_T kDrop, boolean_T iA[96], int32_T *nA,
  int32_T iC[96])
{
  if (kDrop > 0) {
    iA[iC[kDrop - 1] - 1] = false;
    if (kDrop < *nA) {
      for (int32_T i{kDrop}; i < *nA; i++) {
        iC[i - 1] = iC[i];
      }
    }

    iC[*nA - 1] = 0;
    (*nA)--;
  }
}

// Function for MATLAB Function: '<S34>/optimizer'
void egoCar::egoCar_qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16],
  const real_T f[4], const real_T b_Ac[384], const real_T b[96], boolean_T iA[96],
  int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[96], int32_T
  *status)
{
  real_T cTol[96];
  real_T D[16];
  real_T RLinv[16];
  real_T U[16];
  real_T b_H[16];
  real_T Opt[8];
  real_T Rhs[8];
  real_T r[4];
  real_T z[4];
  real_T Xnorm0;
  real_T cMin;
  real_T cVal;
  real_T rMin;
  real_T t;
  real_T zTa;
  real_T zTa_tmp;
  int32_T iC[96];
  int32_T U_tmp;
  int32_T U_tmp_0;
  int32_T b_exponent;
  int32_T exitg1;
  int32_T exitg3;
  int32_T exponent;
  int32_T i;
  int32_T iSave;
  int32_T nA;
  int32_T r_tmp;
  int32_T tmp;
  boolean_T ColdReset;
  boolean_T DualFeasible;
  boolean_T cTolComputed;
  boolean_T exitg2;
  boolean_T exitg4;
  boolean_T guard1;
  boolean_T guard2;
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  *status = 1;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 96; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (i = 0; i < 96; i++) {
    if (iA[i]) {
      nA++;
      iC[nA - 1] = i + 1;
    }
  }

  guard1 = false;
  if (nA > 0) {
    std::memset(&Opt[0], 0, sizeof(real_T) << 3U);
    Rhs[0] = f[0];
    Rhs[4] = 0.0;
    Rhs[1] = f[1];
    Rhs[5] = 0.0;
    Rhs[2] = f[2];
    Rhs[6] = 0.0;
    Rhs[3] = f[3];
    Rhs[7] = 0.0;
    DualFeasible = false;
    tmp = static_cast<int32_T>(std::round(0.3 * static_cast<real_T>(nA)));
    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (*status <= maxiter)) {
        egoCar_KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, egoCar_degrees, RLinv,
                          &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            std::memset(&iA[0], 0, 96U * sizeof(boolean_T));
            std::memset(&iC[0], 0, 96U * sizeof(int32_T));
            ColdReset = true;
          }
        } else {
          for (i = 0; i < nA; i++) {
            Rhs[i + 4] = b[iC[i] - 1];
            for (r_tmp = i + 1; r_tmp <= nA; r_tmp++) {
              U_tmp_0 = ((i << 2) + r_tmp) - 1;
              U[U_tmp_0] = 0.0;
              for (iSave = 0; iSave < nA; iSave++) {
                U_tmp = iSave << 2;
                U[U_tmp_0] += RLinv[(U_tmp + r_tmp) - 1] * RLinv[U_tmp + i];
              }

              U[i + ((r_tmp - 1) << 2)] = U[U_tmp_0];
            }
          }

          for (i = 0; i < 4; i++) {
            Opt[i] = ((b_H[i + 4] * Rhs[1] + b_H[i] * Rhs[0]) + b_H[i + 8] *
                      Rhs[2]) + b_H[i + 12] * Rhs[3];
            for (r_tmp = 0; r_tmp < nA; r_tmp++) {
              Opt[i] += D[(r_tmp << 2) + i] * Rhs[r_tmp + 4];
            }
          }

          Xnorm0 = -1.0E-12;
          i = -1;
          for (r_tmp = 0; r_tmp < nA; r_tmp++) {
            iSave = r_tmp << 2;
            Opt[r_tmp + 4] = ((D[iSave + 1] * Rhs[1] + D[iSave] * Rhs[0]) +
                              D[iSave + 2] * Rhs[2]) + D[iSave + 3] * Rhs[3];
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[r_tmp + 4] += U[(iSave << 2) + r_tmp] * Rhs[iSave + 4];
            }

            cMin = Opt[r_tmp + 4];
            lambda[iC[r_tmp] - 1] = cMin;
            if ((cMin < Xnorm0) && (r_tmp + 1 <= nA)) {
              i = r_tmp;
              Xnorm0 = cMin;
            }
          }

          if (i + 1 <= 0) {
            DualFeasible = true;
            x[0] = Opt[0];
            x[1] = Opt[1];
            x[2] = Opt[2];
            x[3] = Opt[3];
          } else {
            (*status)++;
            if (tmp <= 5) {
              r_tmp = 5;
            } else {
              r_tmp = tmp;
            }

            if (*status > r_tmp) {
              nA = 0;
              std::memset(&iA[0], 0, 96U * sizeof(boolean_T));
              std::memset(&iC[0], 0, 96U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              egoCar_DropConstraint(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          std::memset(&lambda[0], 0, 96U * sizeof(real_T));
          Xnorm0 = f[1];
          cMin = f[0];
          cVal = f[2];
          t = f[3];
          for (tmp = 0; tmp < 4; tmp++) {
            x[tmp] = ((-b_Hinv[tmp + 4] * Xnorm0 + -b_Hinv[tmp] * cMin) +
                      -b_Hinv[tmp + 8] * cVal) + -b_Hinv[tmp + 12] * t;
          }
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    Xnorm0 = f[1];
    cMin = f[0];
    cVal = f[2];
    t = f[3];
    for (tmp = 0; tmp < 4; tmp++) {
      x[tmp] = ((-b_Hinv[tmp + 4] * Xnorm0 + -b_Hinv[tmp] * cMin) + -b_Hinv[tmp
                + 8] * cVal) + -b_Hinv[tmp + 12] * t;
    }

    guard1 = true;
  }

  if (guard1) {
    Xnorm0 = egoCar_norm(x);
    exitg2 = false;
    while ((!exitg2) && (*status <= maxiter)) {
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 96; i++) {
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 96] * x[1]);
          z[2] = std::abs(b_Ac[i + 192] * x[2]);
          z[3] = std::abs(b_Ac[i + 288] * x[3]);
          cTol[i] = std::fmax(cTol[i], egoCar_maximum(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 96] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 192] * x[2])
                   + b_Ac[i + 288] * x[3]) - b[i]) / cTol[i];
          if (cVal < cMin) {
            cMin = cVal;
            tmp = i;
          }
        }
      }

      cTolComputed = true;
      if (tmp + 1 <= 0) {
        exitg2 = true;
      } else if (*status == maxiter) {
        *status = 0;
        exitg2 = true;
      } else {
        do {
          exitg1 = 0;
          if ((tmp + 1 > 0) && (*status <= maxiter)) {
            guard2 = false;
            if (nA == 0) {
              for (r_tmp = 0; r_tmp < 4; r_tmp++) {
                z[r_tmp] = ((b_Hinv[r_tmp + 4] * b_Ac[tmp + 96] + b_Hinv[r_tmp] *
                             b_Ac[tmp]) + b_Hinv[r_tmp + 8] * b_Ac[tmp + 192]) +
                  b_Hinv[r_tmp + 12] * b_Ac[tmp + 288];
              }

              guard2 = true;
            } else {
              egoCar_KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, egoCar_degrees,
                                RLinv, &cMin);
              if (cMin <= 0.0) {
                *status = -2;
                exitg1 = 1;
              } else {
                for (r_tmp = 0; r_tmp < 16; r_tmp++) {
                  U[r_tmp] = -b_H[r_tmp];
                }

                for (r_tmp = 0; r_tmp < 4; r_tmp++) {
                  z[r_tmp] = ((U[r_tmp + 4] * b_Ac[tmp + 96] + U[r_tmp] *
                               b_Ac[tmp]) + U[r_tmp + 8] * b_Ac[tmp + 192]) +
                    U[r_tmp + 12] * b_Ac[tmp + 288];
                }

                for (i = 0; i < nA; i++) {
                  r_tmp = i << 2;
                  r[i] = ((D[r_tmp + 1] * b_Ac[tmp + 96] + D[r_tmp] * b_Ac[tmp])
                          + D[r_tmp + 2] * b_Ac[tmp + 192]) + D[r_tmp + 3] *
                    b_Ac[tmp + 288];
                }

                guard2 = true;
              }
            }

            if (guard2) {
              i = 0;
              cMin = 0.0;
              DualFeasible = true;
              ColdReset = true;
              if (nA > 0) {
                r_tmp = 0;
                exitg4 = false;
                while ((!exitg4) && (r_tmp <= nA - 1)) {
                  if (r[r_tmp] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    r_tmp++;
                  }
                }
              }

              if ((nA != 0) && (!ColdReset)) {
                for (r_tmp = 0; r_tmp < nA; r_tmp++) {
                  cVal = r[r_tmp];
                  if (cVal > 1.0E-12) {
                    cVal = lambda[iC[r_tmp] - 1] / cVal;
                    if ((i == 0) || (cVal < rMin)) {
                      rMin = cVal;
                      i = r_tmp + 1;
                    }
                  }
                }

                if (i > 0) {
                  cMin = rMin;
                  DualFeasible = false;
                }
              }

              cVal = b_Ac[tmp + 96];
              t = b_Ac[tmp + 192];
              zTa_tmp = b_Ac[tmp + 288];
              zTa = ((cVal * z[1] + z[0] * b_Ac[tmp]) + t * z[2]) + zTa_tmp * z
                [3];
              if (zTa <= 0.0) {
                cVal = 0.0;
                ColdReset = true;
              } else {
                cVal = (b[tmp] - (((cVal * x[1] + b_Ac[tmp] * x[0]) + t * x[2])
                                  + zTa_tmp * x[3])) / zTa;
                ColdReset = false;
              }

              if (DualFeasible && ColdReset) {
                *status = -1;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  t = cMin;
                } else if (DualFeasible) {
                  t = cVal;
                } else if (cMin < cVal) {
                  t = cMin;
                } else {
                  t = cVal;
                }

                for (r_tmp = 0; r_tmp < nA; r_tmp++) {
                  iSave = iC[r_tmp];
                  lambda[iSave - 1] -= t * r[r_tmp];
                  if ((iSave <= 96) && (lambda[iSave - 1] < 0.0)) {
                    lambda[iSave - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                std::frexp(1.0, &exponent);
                if (std::abs(t - cMin) < 2.2204460492503131E-16) {
                  egoCar_DropConstraint(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  x[3] += t * z[3];
                  std::frexp(1.0, &b_exponent);
                  if (std::abs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == egoCar_degrees) {
                      *status = -1;
                      exitg1 = 1;
                    } else {
                      nA++;
                      iC[nA - 1] = tmp + 1;
                      i = nA - 1;
                      exitg4 = false;
                      while ((!exitg4) && (i + 1 > 1)) {
                        r_tmp = iC[i - 1];
                        if (iC[i] > r_tmp) {
                          exitg4 = true;
                        } else {
                          iSave = iC[i];
                          iC[i] = r_tmp;
                          iC[i - 1] = iSave;
                          i--;
                        }
                      }

                      iA[tmp] = true;
                      tmp = -1;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            cMin = egoCar_norm(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 96; tmp++) {
                cTol[tmp] = std::fmax(std::abs(b[tmp]), 1.0);
              }

              cTolComputed = false;
            }

            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

// Model step function
void egoCar::step()
{
  // local block i/o variables
  real_T rtb_y;
  real_T rtb_y_n;
  real_T Bc[96];
  real_T a__1[96];
  real_T vseq[62];
  real_T rseq[60];
  real_T f[4];
  real_T rtb_xest[4];
  real_T xk[4];
  real_T rtb_TmpSignalConversionAtSFun_f[2];
  real_T y_innov[2];
  real_T ymax_incr[2];
  real_T ymin_incr[2];
  real_T rtb_y_f3;
  real_T rtb_y_p;
  real_T xk_0;
  real_T xk_1;
  real_T y;
  real_T y_e;
  real_T y_h;
  real_T y_innov_0;
  real_T y_innov_1;
  real_T y_o;
  int32_T d_i;
  int32_T i;
  int8_T rtb_TmpSignalConversionAtSFun_c[2];
  uint8_T b_Mrows;
  boolean_T ymax_incr_flag[2];
  boolean_T ymin_incr_flag[2];
  boolean_T b_Del_Save_Flag0;
  boolean_T tmp;
  boolean_T umax_incr_flag;
  boolean_T umin_incr_flag;
  static const real_T f_a[8]{ -2.8913158861098121E-20, 0.0063245553203367423,
    -2.1779195409073748E-5, -0.012649106890738636, 0.0282842628623524,
    -9.7399522859880982E-6, 0.0, 1.0 };

  static const real_T e_a[8]{ -0.28621232262067842, 0.59329074470859378,
    3.6639065311833567, 0.0089491794894280257, 0.0032542098114138676,
    -0.056579892681977516, -0.014230625940085229, 0.094687375766569978 };

  static const real_T b_Mlim[96]{ 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.4, 0.4, 0.4, 0.6,
    0.6, 0.6 };

  static const real_T b_Mv[5952]{ -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -6.7762635780344027E-20, -6.0986372202309624E-20, -5.4210108624275222E-20,
    -6.0986372202309624E-20, -5.4210108624275222E-20, -6.7762635780344027E-20,
    -7.453889935837843E-20, -6.7762635780344027E-20, -6.7762635780344027E-20,
    -8.1315162936412833E-20, -8.1315162936412833E-20, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.099999999999999936,
    2.7105054312137611E-20, 0.099999999999999936, 3.3881317890172014E-20,
    0.099999999999999922, 4.0657581468206416E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.0999999999999999, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999853, 5.4210108624275222E-20,
    0.099999999999999839, 5.4210108624275222E-20, 0.099999999999999825,
    5.4210108624275222E-20, 0.099999999999999811, 5.4210108624275222E-20,
    0.0999999999999998, 5.4210108624275222E-20, 0.099999999999999784,
    6.0986372202309624E-20, 0.09999999999999977, 6.7762635780344027E-20,
    0.099999999999999756, 6.0986372202309624E-20, 0.099999999999999756,
    5.4210108624275222E-20, 0.099999999999999742, 6.0986372202309624E-20,
    0.099999999999999728, 5.4210108624275222E-20, 0.099999999999999714,
    6.7762635780344027E-20, 0.099999999999999714, 7.453889935837843E-20,
    0.0999999999999997, 6.7762635780344027E-20, 0.099999999999999686,
    6.7762635780344027E-20, 0.099999999999999659, 8.1315162936412833E-20,
    0.099999999999999659, 8.1315162936412833E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -3.3881317890172014E-20, -4.0657581468206416E-20,
    -4.7433845046240819E-20, -4.7433845046240819E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -6.0986372202309624E-20, -6.7762635780344027E-20, -6.0986372202309624E-20,
    -5.4210108624275222E-20, -6.0986372202309624E-20, -5.4210108624275222E-20,
    -6.7762635780344027E-20, -7.453889935837843E-20, -6.7762635780344027E-20,
    -6.7762635780344027E-20, -8.1315162936412833E-20, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999853,
    5.4210108624275222E-20, 0.099999999999999839, 5.4210108624275222E-20,
    0.099999999999999825, 5.4210108624275222E-20, 0.099999999999999811,
    5.4210108624275222E-20, 0.0999999999999998, 5.4210108624275222E-20,
    0.099999999999999784, 6.0986372202309624E-20, 0.09999999999999977,
    6.7762635780344027E-20, 0.099999999999999756, 6.0986372202309624E-20,
    0.099999999999999756, 5.4210108624275222E-20, 0.099999999999999742,
    6.0986372202309624E-20, 0.099999999999999728, 5.4210108624275222E-20,
    0.099999999999999714, 6.7762635780344027E-20, 0.099999999999999714,
    7.453889935837843E-20, 0.0999999999999997, 6.7762635780344027E-20,
    0.099999999999999686, 6.7762635780344027E-20, 0.099999999999999659,
    8.1315162936412833E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -6.7762635780344027E-20, -6.0986372202309624E-20, -5.4210108624275222E-20,
    -6.0986372202309624E-20, -5.4210108624275222E-20, -6.7762635780344027E-20,
    -7.453889935837843E-20, -6.7762635780344027E-20, -6.7762635780344027E-20,
    0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.099999999999999922,
    4.0657581468206416E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0999999999999999,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999853, 5.4210108624275222E-20, 0.099999999999999839,
    5.4210108624275222E-20, 0.099999999999999825, 5.4210108624275222E-20,
    0.099999999999999811, 5.4210108624275222E-20, 0.0999999999999998,
    5.4210108624275222E-20, 0.099999999999999784, 6.0986372202309624E-20,
    0.09999999999999977, 6.7762635780344027E-20, 0.099999999999999756,
    6.0986372202309624E-20, 0.099999999999999756, 5.4210108624275222E-20,
    0.099999999999999742, 6.0986372202309624E-20, 0.099999999999999728,
    5.4210108624275222E-20, 0.099999999999999714, 6.7762635780344027E-20,
    0.099999999999999714, 7.453889935837843E-20, 0.0999999999999997,
    6.7762635780344027E-20, 0.099999999999999686, 6.7762635780344027E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -3.3881317890172014E-20,
    -4.0657581468206416E-20, -4.7433845046240819E-20, -4.7433845046240819E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -6.0986372202309624E-20, -6.7762635780344027E-20,
    -6.0986372202309624E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -5.4210108624275222E-20, -6.7762635780344027E-20, -7.453889935837843E-20,
    -6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.099999999999999936,
    2.7105054312137611E-20, 0.099999999999999936, 3.3881317890172014E-20,
    0.099999999999999922, 4.0657581468206416E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.0999999999999999, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999853, 5.4210108624275222E-20,
    0.099999999999999839, 5.4210108624275222E-20, 0.099999999999999825,
    5.4210108624275222E-20, 0.099999999999999811, 5.4210108624275222E-20,
    0.0999999999999998, 5.4210108624275222E-20, 0.099999999999999784,
    6.0986372202309624E-20, 0.09999999999999977, 6.7762635780344027E-20,
    0.099999999999999756, 6.0986372202309624E-20, 0.099999999999999756,
    5.4210108624275222E-20, 0.099999999999999742, 6.0986372202309624E-20,
    0.099999999999999728, 5.4210108624275222E-20, 0.099999999999999714,
    6.7762635780344027E-20, 0.099999999999999714, 7.453889935837843E-20,
    0.0999999999999997, 6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -2.0328790734103208E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -3.3881317890172014E-20,
    -4.0657581468206416E-20, -4.7433845046240819E-20, -4.7433845046240819E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -6.0986372202309624E-20, -6.7762635780344027E-20,
    -6.0986372202309624E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -5.4210108624275222E-20, -6.7762635780344027E-20, -7.453889935837843E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.099999999999999936,
    2.7105054312137611E-20, 0.099999999999999936, 3.3881317890172014E-20,
    0.099999999999999922, 4.0657581468206416E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.0999999999999999, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999853, 5.4210108624275222E-20,
    0.099999999999999839, 5.4210108624275222E-20, 0.099999999999999825,
    5.4210108624275222E-20, 0.099999999999999811, 5.4210108624275222E-20,
    0.0999999999999998, 5.4210108624275222E-20, 0.099999999999999784,
    6.0986372202309624E-20, 0.09999999999999977, 6.7762635780344027E-20,
    0.099999999999999756, 6.0986372202309624E-20, 0.099999999999999756,
    5.4210108624275222E-20, 0.099999999999999742, 6.0986372202309624E-20,
    0.099999999999999728, 5.4210108624275222E-20, 0.099999999999999714,
    6.7762635780344027E-20, 0.099999999999999714, 7.453889935837843E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -6.7762635780344027E-20, -6.0986372202309624E-20, -5.4210108624275222E-20,
    -6.0986372202309624E-20, -5.4210108624275222E-20, -6.7762635780344027E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.099999999999999936,
    2.7105054312137611E-20, 0.099999999999999936, 3.3881317890172014E-20,
    0.099999999999999922, 4.0657581468206416E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.0999999999999999, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999853, 5.4210108624275222E-20,
    0.099999999999999839, 5.4210108624275222E-20, 0.099999999999999825,
    5.4210108624275222E-20, 0.099999999999999811, 5.4210108624275222E-20,
    0.0999999999999998, 5.4210108624275222E-20, 0.099999999999999784,
    6.0986372202309624E-20, 0.09999999999999977, 6.7762635780344027E-20,
    0.099999999999999756, 6.0986372202309624E-20, 0.099999999999999756,
    5.4210108624275222E-20, 0.099999999999999742, 6.0986372202309624E-20,
    0.099999999999999728, 5.4210108624275222E-20, 0.099999999999999714,
    6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -3.3881317890172014E-20, -4.0657581468206416E-20,
    -4.7433845046240819E-20, -4.7433845046240819E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -6.0986372202309624E-20, -6.7762635780344027E-20, -6.0986372202309624E-20,
    -5.4210108624275222E-20, -6.0986372202309624E-20, -5.4210108624275222E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999853,
    5.4210108624275222E-20, 0.099999999999999839, 5.4210108624275222E-20,
    0.099999999999999825, 5.4210108624275222E-20, 0.099999999999999811,
    5.4210108624275222E-20, 0.0999999999999998, 5.4210108624275222E-20,
    0.099999999999999784, 6.0986372202309624E-20, 0.09999999999999977,
    6.7762635780344027E-20, 0.099999999999999756, 6.0986372202309624E-20,
    0.099999999999999756, 5.4210108624275222E-20, 0.099999999999999742,
    6.0986372202309624E-20, 0.099999999999999728, 5.4210108624275222E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -6.7762635780344027E-20, -6.0986372202309624E-20, -5.4210108624275222E-20,
    -6.0986372202309624E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.099999999999999922,
    4.0657581468206416E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0999999999999999,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999853, 5.4210108624275222E-20, 0.099999999999999839,
    5.4210108624275222E-20, 0.099999999999999825, 5.4210108624275222E-20,
    0.099999999999999811, 5.4210108624275222E-20, 0.0999999999999998,
    5.4210108624275222E-20, 0.099999999999999784, 6.0986372202309624E-20,
    0.09999999999999977, 6.7762635780344027E-20, 0.099999999999999756,
    6.0986372202309624E-20, 0.099999999999999756, 5.4210108624275222E-20,
    0.099999999999999742, 6.0986372202309624E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -6.7762635780344027E-20, -6.0986372202309624E-20, -5.4210108624275222E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999853,
    5.4210108624275222E-20, 0.099999999999999839, 5.4210108624275222E-20,
    0.099999999999999825, 5.4210108624275222E-20, 0.099999999999999811,
    5.4210108624275222E-20, 0.0999999999999998, 5.4210108624275222E-20,
    0.099999999999999784, 6.0986372202309624E-20, 0.09999999999999977,
    6.7762635780344027E-20, 0.099999999999999756, 6.0986372202309624E-20,
    0.099999999999999756, 5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -6.0986372202309624E-20,
    -6.7762635780344027E-20, -6.0986372202309624E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999853,
    5.4210108624275222E-20, 0.099999999999999839, 5.4210108624275222E-20,
    0.099999999999999825, 5.4210108624275222E-20, 0.099999999999999811,
    5.4210108624275222E-20, 0.0999999999999998, 5.4210108624275222E-20,
    0.099999999999999784, 6.0986372202309624E-20, 0.09999999999999977,
    6.7762635780344027E-20, 0.099999999999999756, 6.0986372202309624E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -3.3881317890172014E-20, -4.0657581468206416E-20,
    -4.7433845046240819E-20, -4.7433845046240819E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -6.0986372202309624E-20, -6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999853,
    5.4210108624275222E-20, 0.099999999999999839, 5.4210108624275222E-20,
    0.099999999999999825, 5.4210108624275222E-20, 0.099999999999999811,
    5.4210108624275222E-20, 0.0999999999999998, 5.4210108624275222E-20,
    0.099999999999999784, 6.0986372202309624E-20, 0.09999999999999977,
    6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -2.0328790734103208E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -3.3881317890172014E-20,
    -4.0657581468206416E-20, -4.7433845046240819E-20, -4.7433845046240819E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -6.0986372202309624E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999853,
    5.4210108624275222E-20, 0.099999999999999839, 5.4210108624275222E-20,
    0.099999999999999825, 5.4210108624275222E-20, 0.099999999999999811,
    5.4210108624275222E-20, 0.0999999999999998, 5.4210108624275222E-20,
    0.099999999999999784, 6.0986372202309624E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -3.3881317890172014E-20, -4.0657581468206416E-20,
    -4.7433845046240819E-20, -4.7433845046240819E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.099999999999999936,
    2.7105054312137611E-20, 0.099999999999999936, 3.3881317890172014E-20,
    0.099999999999999922, 4.0657581468206416E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.0999999999999999, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999853, 5.4210108624275222E-20,
    0.099999999999999839, 5.4210108624275222E-20, 0.099999999999999825,
    5.4210108624275222E-20, 0.099999999999999811, 5.4210108624275222E-20,
    0.0999999999999998, 5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -3.3881317890172014E-20, -4.0657581468206416E-20,
    -4.7433845046240819E-20, -4.7433845046240819E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.099999999999999922,
    4.0657581468206416E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0999999999999999,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999853, 5.4210108624275222E-20, 0.099999999999999839,
    5.4210108624275222E-20, 0.099999999999999825, 5.4210108624275222E-20,
    0.099999999999999811, 5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -3.3881317890172014E-20, -4.0657581468206416E-20,
    -4.7433845046240819E-20, -4.7433845046240819E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.099999999999999922,
    4.0657581468206416E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0999999999999999,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999853, 5.4210108624275222E-20, 0.099999999999999839,
    5.4210108624275222E-20, 0.099999999999999825, 5.4210108624275222E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -3.3881317890172014E-20,
    -4.0657581468206416E-20, -4.7433845046240819E-20, -4.7433845046240819E-20,
    -4.0657581468206416E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -5.4210108624275222E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.099999999999999922,
    4.0657581468206416E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0999999999999999,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999867, 5.4210108624275222E-20,
    0.099999999999999853, 5.4210108624275222E-20, 0.099999999999999839,
    5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    -5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.099999999999999936,
    2.7105054312137611E-20, 0.099999999999999936, 3.3881317890172014E-20,
    0.099999999999999922, 4.0657581468206416E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.0999999999999999, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999867, 5.4210108624275222E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999853, 5.4210108624275222E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, -5.4210108624275222E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.099999999999999922,
    4.0657581468206416E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0999999999999999,
    4.0657581468206416E-20, 0.099999999999999881, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999867,
    5.4210108624275222E-20, 0.099999999999999867, 5.4210108624275222E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, -5.4210108624275222E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.099999999999999867, 5.4210108624275222E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.099999999999999881,
    4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -3.3881317890172014E-20, -4.0657581468206416E-20,
    -4.7433845046240819E-20, -4.7433845046240819E-20, -4.0657581468206416E-20,
    -4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.0999999999999999, 4.0657581468206416E-20,
    0.099999999999999881, 4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, -4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.099999999999999936,
    2.7105054312137611E-20, 0.099999999999999936, 3.3881317890172014E-20,
    0.099999999999999922, 4.0657581468206416E-20, 0.099999999999999908,
    4.7433845046240819E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.0999999999999999, 4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    -4.7433845046240819E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.099999999999999922,
    4.0657581468206416E-20, 0.099999999999999908, 4.7433845046240819E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, -4.7433845046240819E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20,
    0.099999999999999908, 4.7433845046240819E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -2.0328790734103208E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -3.3881317890172014E-20, -4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.099999999999999936,
    3.3881317890172014E-20, 0.099999999999999922, 4.0657581468206416E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -2.0328790734103208E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, -2.7105054312137611E-20, -3.3881317890172014E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.09999999999999995,
    2.7105054312137611E-20, 0.099999999999999936, 2.7105054312137611E-20,
    0.099999999999999936, 3.3881317890172014E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    -2.7105054312137611E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.099999999999999978, 2.0328790734103208E-20, 0.099999999999999964,
    2.7105054312137611E-20, 0.09999999999999995, 2.7105054312137611E-20,
    0.099999999999999936, 2.7105054312137611E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, -2.7105054312137611E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978,
    2.0328790734103208E-20, 0.099999999999999964, 2.7105054312137611E-20,
    0.09999999999999995, 2.7105054312137611E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, -2.7105054312137611E-20, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20,
    0.099999999999999964, 2.7105054312137611E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -2.0328790734103208E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.099999999999999978, 2.0328790734103208E-20, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mx[384]{ -0.0051781079403026538, -0.0042394762134830653,
    -0.0034709895529211774, -0.0028418058905889595, -0.0023266738769033403,
    -0.0019049194554039356, -0.0015596161402757595, -0.0012769056970405443,
    -0.0010454419629475661, -0.00085593548562338053, -0.00070078070473059454,
    -0.0005737507141265946, -0.00046974735425589281, -0.00038459660510631273,
    -0.00031488106812992269, -0.000257802814040007, -0.00021107109208459395,
    -0.00017281039417540775, -0.00014148518416293054, -0.00011583827137908915,
    -9.4840355161449356E-5, -7.7648715403512756E-5, -6.3573391237851953E-5,
    -5.2049490483883769E-5, -4.2614518541191467E-5, -3.4889816857281371E-5,
    -2.8565366030310904E-5, -2.3387343641940707E-5, -1.9147937472455725E-5,
    -1.5677005266709844E-5, -0.00057322369001704143, 0.0051781079403026538,
    -0.0010425395534268329, 0.0042394762134830653, -0.0014267828837077742,
    0.0034709895529211774, -0.0017413747148738808, 0.0028418058905889595,
    -0.001998940721716688, 0.0023266738769033403, -0.0022098179324663877,
    0.0019049194554039356, -0.0023824695900304732, 0.0015596161402757595,
    -0.0025238248116480782, 0.0012769056970405443, -0.0026395566786945645,
    0.0010454419629475661, -0.0027343099173566545, 0.00085593548562338053,
    -0.0028118873078030451, 0.00070078070473059454, -0.0028754023031050428,
    0.0005737507141265946, -0.0029274039830403916, 0.00046974735425589281,
    -0.0029699793576151794, 0.00038459660510631273, -0.0030048371261033719,
    0.00031488106812992269, -0.0030333762531483273, 0.000257802814040007,
    -0.003056742114126031, 0.00021107109208459395, -0.0030758724630806219,
    0.00017281039417540775, -0.0030915350680868576, 0.00014148518416293054,
    -0.0031043585244787762, 0.00011583827137908915, -0.0031148574825875935,
    9.4840355161449356E-5, -0.0031234533024665592, 7.7648715403512756E-5,
    -0.0031304909645493871, 6.3573391237851953E-5, -0.003136252914926369,
    5.2049490483883769E-5, -0.0031409704008977125, 4.2614518541191467E-5,
    -0.0031448327517396651, 3.4889816857281371E-5, -0.0031479949771531478,
    2.8565366030310904E-5, -0.00315058398834733, 2.3387343641940707E-5,
    -0.00315270369143207, 1.9147937472455725E-5, -0.00315443915753494,
    1.5677005266709844E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.013222330410818706,
    0.013691646135095874, 0.014075889351464657, 0.014390481089367378,
    0.014648047019852584, 0.014858924168085968, 0.015031575774466021,
    0.015172930954177686, 0.01528866278691449, 0.01538341599748619,
    0.015460993364934115, 0.01552450834140656, 0.015576510005925576,
    0.015619085367878538, 0.015653943126032854, 0.015682482244617147,
    0.015705848098667846, 0.015724978441951087, 0.015740641042314017,
    0.015753464494904314, 0.015763963449900626, 0.015772559267231291,
    0.015779596927227746, 0.015785358875896553, 0.015790076360469363,
    0.015793938710166291, 0.015797100934642306, 0.015799689945068959,
    0.015801809647525297, 0.015803545113113673, 0.0012727474058932024,
    -0.013222330410818706, 0.0026192279050969297, -0.013691646135095874,
    0.0040082446582548484, -0.014075889351464657, 0.0054320871506457992,
    -0.014390481089367378, 0.0068844425467455067, -0.014648047019852584,
    0.008360142333971126, -0.014858924168085968, 0.009854954892123409,
    -0.015031575774466021, 0.011365415663609886, -0.015172930954177686,
    0.012888688108583794, -0.01528866278691449, 0.014422449864640252,
    -0.01538341599748619, 0.015964799542258597, -0.015460993364934115,
    0.017514180415364684, -0.01552450834140656, 0.019069317944447484,
    -0.015576510005925576, 0.020629168624813315, -0.015619085367878538,
    0.022192878107078467, -0.015653943126032854, 0.023759746909128628,
    -0.015682482244617147, 0.025329202343445587, -0.015705848098667846,
    0.026900775533146275, -0.015724978441951087, 0.028474082594307119,
    -0.015740641042314017, 0.030048809229354274, -0.015753464494904314,
    0.031624698113198423, -0.015763963449900626, 0.0332015385658754,
    -0.015772559267231291, 0.034779158097219473, -0.015779596927227746,
    0.036357415484227378, -0.015785358875896553, 0.037936195103283277,
    -0.015790076360469363, 0.039515402289777117, -0.015793938710166291,
    0.041094959538881412, -0.015797100934642306, 0.04267480339501039,
    -0.015799689945068959, 0.044254881905124527, -0.015801809647525297,
    0.045835152533672645, -0.015803545113113673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0181340739972458E-5, 1.0542719041302999E-5, 1.0838590570097395E-5,
    1.1080829689681565E-5, 1.1279158306483658E-5, 1.144153604427495E-5,
    1.1574479691819914E-5, 1.1683324744491336E-5, 1.1772439536433822E-5,
    1.1845400557151281E-5, 1.1905135988588624E-5, 1.1954043223354757E-5,
    1.1994085080505796E-5, 1.2026868580365706E-5, 1.2053709439894543E-5,
    1.2075684877029846E-5, 1.2093676843224852E-5, 1.2108407419257041E-5,
    1.2120467794865149E-5, 1.2130341995269178E-5, 1.2138426306802011E-5,
    1.2145045181271405E-5, 1.2150464257350259E-5, 1.2154901021589286E-5,
    1.2158533536915936E-5, 1.2161507588924891E-5, 1.2163942536765875E-5,
    1.2165936103445428E-5, 1.2167568297794291E-5, 1.2168904625502704E-5,
    0.028285259662161152, -1.0181340739972458E-5, 0.028286296467046231,
    -1.0542719041302999E-5, 0.028287366025317579, -1.0838590570097395E-5,
    0.028288462399793535, -1.1080829689681565E-5, 0.028289580729520879,
    -1.1279158306483658E-5, 0.028290717034687728, -1.144153604427495E-5,
    0.0282918680568997, -1.1574479691819914E-5, 0.028293031128409113,
    -1.1683324744491336E-5, 0.028294204065048886, -1.1772439536433822E-5,
    0.028295385078574275, -1.1845400557151281E-5, 0.0282965727048943,
    -1.1905135988588624E-5, 0.028297765745312665, -1.1954043223354757E-5,
    0.028298963218419834, -1.1994085080505796E-5, 0.028300164320705649,
    -1.2026868580365706E-5, 0.028301368394311629, -1.2053709439894543E-5,
    0.028302574900628808, -1.2075684877029846E-5, 0.028303783398681455,
    -1.2093676843224852E-5, 0.028304993527429184, -1.2108407419257041E-5,
    0.028306204991277124, -1.2120467794865149E-5, 0.028307417548212668,
    -1.2130341995269178E-5, 0.028308631000092649, -1.2138426306802011E-5,
    0.028309845184691158, -1.2145045181271405E-5, 0.028311059969188864,
    -1.2150464257350259E-5, 0.028312275244842491, -1.2154901021589286E-5,
    0.028313490922620574, -1.2158533536915936E-5, 0.028314706929630312,
    -1.2161507588924891E-5, 0.028315923206192136, -1.2163942536765875E-5,
    0.02831713970344454, -1.2165936103445428E-5, 0.028318356381383112,
    -1.2167568297794291E-5, 0.028319573207255002, -1.2168904625502704E-5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0 };

  static const real_T b_Mu1[96]{ -0.00093653765389914228, -0.0035160023017820649,
    -0.0074405818047014721, -0.012466448205861281, -0.018393972058572367,
    -0.025059710595610407, -0.032329848197080679, -0.040094825899733177,
    -0.048264944411079785, -0.056766764161831143, -0.065540157918117251,
    -0.074535897664471232, -0.083713678910717348, -0.0930405031312616,
    -0.10248935341839394, -0.1120381101989191, -0.12166866349801714,
    -0.1313661861223655, -0.14111853859280921, -0.15091578194443769,
    -0.16074977884102493, -0.17061386699515452, -0.18050259178723282,
    -0.1904114873524522, -0.20033689734995552, -0.21027582822103932,
    -0.22022582904713195, -0.2301848931858255, -0.24015137773727019,
    -0.25012393760883478, -3.1731173050451781E-5, 0.00093653765389914228,
    -0.00024199884910901012, 0.0035160023017820649, -0.000779709097649322,
    0.0074405818047014721, -0.0017667758970694289, 0.012466448205861281,
    -0.0033030139707138926, 0.018393972058572367, -0.0054701447021948749,
    0.025059710595610407, -0.0083350759014597375, 0.032329848197080679,
    -0.011952587050133481, 0.040094825899733177, -0.016367527794460168,
    0.048264944411079785, -0.021616617919084473, 0.056766764161831143,
    -0.027729921040941397, 0.065540157918117251, -0.034732051167764379,
    0.074535897664471232, -0.04264316054464129, 0.083713678910717348,
    -0.051479748434369123, 0.0930405031312616, -0.0612553232908029,
    0.10248935341839394, -0.071980944900540272, 0.1120381101989191,
    -0.0836656682509912, 0.12166866349801714, -0.096316906938816949,
    0.1313661861223655, -0.10994073070359503, 0.14111853859280921,
    -0.12454210902778072, 0.15091578194443769, -0.14012511057948704,
    0.16074977884102493, -0.15669306650242215, 0.17061386699515452,
    -0.17424870410638291, 0.18050259178723282, -0.19279425632377312,
    0.1904114873524522, -0.21233155132502135, 0.20033689734995552,
    -0.23286208588947932, 0.21027582822103932, -0.25438708547643285,
    0.22022582904713195, -0.27690755340708595, 0.2301848931858255,
    -0.30042431113136348, 0.24015137773727019, -0.324938031195581,
    0.25012393760883478, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0 };

  static const uint8_T b_Mrows_0[96]{ 2U, 4U, 6U, 8U, 10U, 12U, 14U, 16U, 18U,
    20U, 22U, 24U, 26U, 28U, 30U, 32U, 34U, 36U, 38U, 40U, 42U, 44U, 46U, 48U,
    50U, 52U, 54U, 56U, 58U, 60U, 61U, 62U, 63U, 64U, 65U, 66U, 67U, 68U, 69U,
    70U, 71U, 72U, 73U, 74U, 75U, 76U, 77U, 78U, 79U, 80U, 81U, 82U, 83U, 84U,
    85U, 86U, 87U, 88U, 89U, 90U, 91U, 92U, 93U, 94U, 95U, 96U, 97U, 98U, 99U,
    100U, 101U, 102U, 103U, 104U, 105U, 106U, 107U, 108U, 109U, 110U, 111U, 112U,
    113U, 114U, 115U, 116U, 117U, 118U, 119U, 120U, 121U, 122U, 123U, 151U, 152U,
    153U };

  static const real_T b_Linv[16]{ 7.9998080541108729, -2.9530049756913286,
    -2.1896728723561849, 0.0, 0.0, 8.7034491930031166, -2.0693880092142272, 0.0,
    0.0, 0.0, 9.0725349962825437, 0.0, 0.0, 0.0, 0.0, 0.001 };

  static const real_T b_Hinv[16]{ 77.511834577007519, -21.170045986459318,
    -19.865883764862005, 0.0, -21.170045986459318, 80.032394587866222,
    -18.77459513448354, 0.0, -19.865883764862005, -18.77459513448354,
    82.310891258771491, 0.0, 0.0, 0.0, 0.0, 1.0E-6 };

  static const real_T b_Ac[384]{ -0.00093653765389914228, -0.0035160023017820649,
    -0.0074405818047014721, -0.012466448205861281, -0.018393972058572367,
    -0.025059710595610407, -0.032329848197080679, -0.040094825899733177,
    -0.048264944411079785, -0.056766764161831143, -0.065540157918117251,
    -0.074535897664471232, -0.083713678910717348, -0.0930405031312616,
    -0.10248935341839394, -0.1120381101989191, -0.12166866349801714,
    -0.1313661861223655, -0.14111853859280921, -0.15091578194443769,
    -0.16074977884102493, -0.17061386699515452, -0.18050259178723282,
    -0.1904114873524522, -0.20033689734995552, -0.21027582822103932,
    -0.22022582904713195, -0.2301848931858255, -0.24015137773727019,
    -0.25012393760883478, -3.1731173050451781E-5, 0.00093653765389914228,
    -0.00024199884910901012, 0.0035160023017820649, -0.000779709097649322,
    0.0074405818047014721, -0.0017667758970694289, 0.012466448205861281,
    -0.0033030139707138926, 0.018393972058572367, -0.0054701447021948749,
    0.025059710595610407, -0.0083350759014597375, 0.032329848197080679,
    -0.011952587050133481, 0.040094825899733177, -0.016367527794460168,
    0.048264944411079785, -0.021616617919084473, 0.056766764161831143,
    -0.027729921040941397, 0.065540157918117251, -0.034732051167764379,
    0.074535897664471232, -0.04264316054464129, 0.083713678910717348,
    -0.051479748434369123, 0.0930405031312616, -0.0612553232908029,
    0.10248935341839394, -0.071980944900540272, 0.1120381101989191,
    -0.0836656682509912, 0.12166866349801714, -0.096316906938816949,
    0.1313661861223655, -0.10994073070359503, 0.14111853859280921,
    -0.12454210902778072, 0.15091578194443769, -0.14012511057948704,
    0.16074977884102493, -0.15669306650242215, 0.17061386699515452,
    -0.17424870410638291, 0.18050259178723282, -0.19279425632377312,
    0.1904114873524522, -0.21233155132502135, 0.20033689734995552,
    -0.23286208588947932, 0.21027582822103932, -0.25438708547643285,
    0.22022582904713195, -0.27690755340708595, 0.2301848931858255,
    -0.30042431113136348, 0.24015137773727019, -0.324938031195581,
    0.25012393760883478, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.0,
    -0.00093653765389914228, -0.0035160023017820649, -0.0074405818047014721,
    -0.012466448205861281, -0.018393972058572367, -0.025059710595610407,
    -0.032329848197080679, -0.040094825899733177, -0.048264944411079785,
    -0.056766764161831143, -0.065540157918117251, -0.074535897664471232,
    -0.083713678910717348, -0.0930405031312616, -0.10248935341839394,
    -0.1120381101989191, -0.12166866349801714, -0.1313661861223655,
    -0.14111853859280921, -0.15091578194443769, -0.16074977884102493,
    -0.17061386699515452, -0.18050259178723282, -0.1904114873524522,
    -0.20033689734995552, -0.21027582822103932, -0.22022582904713195,
    -0.2301848931858255, -0.24015137773727019, 0.0, 0.0, -3.1731173050451781E-5,
    0.00093653765389914228, -0.00024199884910901012, 0.0035160023017820649,
    -0.000779709097649322, 0.0074405818047014721, -0.0017667758970694289,
    0.012466448205861281, -0.0033030139707138926, 0.018393972058572367,
    -0.0054701447021948749, 0.025059710595610407, -0.0083350759014597375,
    0.032329848197080679, -0.011952587050133481, 0.040094825899733177,
    -0.016367527794460168, 0.048264944411079785, -0.021616617919084473,
    0.056766764161831143, -0.027729921040941397, 0.065540157918117251,
    -0.034732051167764379, 0.074535897664471232, -0.04264316054464129,
    0.083713678910717348, -0.051479748434369123, 0.0930405031312616,
    -0.0612553232908029, 0.10248935341839394, -0.071980944900540272,
    0.1120381101989191, -0.0836656682509912, 0.12166866349801714,
    -0.096316906938816949, 0.1313661861223655, -0.10994073070359503,
    0.14111853859280921, -0.12454210902778072, 0.15091578194443769,
    -0.14012511057948704, 0.16074977884102493, -0.15669306650242215,
    0.17061386699515452, -0.17424870410638291, 0.18050259178723282,
    -0.19279425632377312, 0.1904114873524522, -0.21233155132502135,
    0.20033689734995552, -0.23286208588947932, 0.21027582822103932,
    -0.25438708547643285, 0.22022582904713195, -0.27690755340708595,
    0.2301848931858255, -0.30042431113136348, 0.24015137773727019, -0.0, -1.0,
    -1.0, 0.0, 1.0, 1.0, -0.0, -0.0, -0.00093653765389914228,
    -0.0035160023017820649, -0.0074405818047014721, -0.012466448205861281,
    -0.018393972058572367, -0.025059710595610407, -0.032329848197080679,
    -0.040094825899733177, -0.048264944411079785, -0.056766764161831143,
    -0.065540157918117251, -0.074535897664471232, -0.083713678910717348,
    -0.0930405031312616, -0.10248935341839394, -0.1120381101989191,
    -0.12166866349801714, -0.1313661861223655, -0.14111853859280921,
    -0.15091578194443769, -0.16074977884102493, -0.17061386699515452,
    -0.18050259178723282, -0.1904114873524522, -0.20033689734995552,
    -0.21027582822103932, -0.22022582904713195, -0.2301848931858255, 0.0, 0.0,
    0.0, 0.0, -3.1731173050451781E-5, 0.00093653765389914228,
    -0.00024199884910901012, 0.0035160023017820649, -0.000779709097649322,
    0.0074405818047014721, -0.0017667758970694289, 0.012466448205861281,
    -0.0033030139707138926, 0.018393972058572367, -0.0054701447021948749,
    0.025059710595610407, -0.0083350759014597375, 0.032329848197080679,
    -0.011952587050133481, 0.040094825899733177, -0.016367527794460168,
    0.048264944411079785, -0.021616617919084473, 0.056766764161831143,
    -0.027729921040941397, 0.065540157918117251, -0.034732051167764379,
    0.074535897664471232, -0.04264316054464129, 0.083713678910717348,
    -0.051479748434369123, 0.0930405031312616, -0.0612553232908029,
    0.10248935341839394, -0.071980944900540272, 0.1120381101989191,
    -0.0836656682509912, 0.12166866349801714, -0.096316906938816949,
    0.1313661861223655, -0.10994073070359503, 0.14111853859280921,
    -0.12454210902778072, 0.15091578194443769, -0.14012511057948704,
    0.16074977884102493, -0.15669306650242215, 0.17061386699515452,
    -0.17424870410638291, 0.18050259178723282, -0.19279425632377312,
    0.1904114873524522, -0.21233155132502135, 0.20033689734995552,
    -0.23286208588947932, 0.21027582822103932, -0.25438708547643285,
    0.22022582904713195, -0.27690755340708595, 0.2301848931858255, -0.0, -0.0,
    -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Kr[180]{ 0.0, -9.3653765389914245E-6, 0.0,
    -3.5160023017820657E-5, 0.0, -7.4405818047014734E-5, 0.0,
    -0.00012466448205861283, 0.0, -0.0001839397205857237, 0.0,
    -0.00025059710595610411, 0.0, -0.00032329848197080687, 0.0,
    -0.00040094825899733187, 0.0, -0.00048264944411079792, 0.0,
    -0.00056766764161831154, 0.0, -0.00065540157918117265, 0.0,
    -0.00074535897664471248, 0.0, -0.00083713678910717364, 0.0,
    -0.00093040503131261616, 0.0, -0.0010248935341839395, 0.0,
    -0.0011203811019891911, 0.0, -0.0012166866349801716, 0.0,
    -0.0013136618612236554, 0.0, -0.0014111853859280923, 0.0,
    -0.0015091578194443773, 0.0, -0.0016074977884102497, 0.0,
    -0.0017061386699515455, 0.0, -0.0018050259178723285, 0.0,
    -0.0019041148735245224, 0.0, -0.0020033689734995554, 0.0,
    -0.0021027582822103937, 0.0, -0.00220225829047132, 0.0,
    -0.0023018489318582555, 0.0, -0.0024015137773727023, 0.0,
    -0.0025012393760883481, -0.0, -0.0, 0.0, -9.3653765389914245E-6, 0.0,
    -3.5160023017820657E-5, 0.0, -7.4405818047014734E-5, 0.0,
    -0.00012466448205861283, 0.0, -0.0001839397205857237, 0.0,
    -0.00025059710595610411, 0.0, -0.00032329848197080687, 0.0,
    -0.00040094825899733187, 0.0, -0.00048264944411079792, 0.0,
    -0.00056766764161831154, 0.0, -0.00065540157918117265, 0.0,
    -0.00074535897664471248, 0.0, -0.00083713678910717364, 0.0,
    -0.00093040503131261616, 0.0, -0.0010248935341839395, 0.0,
    -0.0011203811019891911, 0.0, -0.0012166866349801716, 0.0,
    -0.0013136618612236554, 0.0, -0.0014111853859280923, 0.0,
    -0.0015091578194443773, 0.0, -0.0016074977884102497, 0.0,
    -0.0017061386699515455, 0.0, -0.0018050259178723285, 0.0,
    -0.0019041148735245224, 0.0, -0.0020033689734995554, 0.0,
    -0.0021027582822103937, 0.0, -0.00220225829047132, 0.0,
    -0.0023018489318582555, 0.0, -0.0024015137773727023, -0.0, -0.0, -0.0, -0.0,
    0.0, -9.3653765389914245E-6, 0.0, -3.5160023017820657E-5, 0.0,
    -7.4405818047014734E-5, 0.0, -0.00012466448205861283, 0.0,
    -0.0001839397205857237, 0.0, -0.00025059710595610411, 0.0,
    -0.00032329848197080687, 0.0, -0.00040094825899733187, 0.0,
    -0.00048264944411079792, 0.0, -0.00056766764161831154, 0.0,
    -0.00065540157918117265, 0.0, -0.00074535897664471248, 0.0,
    -0.00083713678910717364, 0.0, -0.00093040503131261616, 0.0,
    -0.0010248935341839395, 0.0, -0.0011203811019891911, 0.0,
    -0.0012166866349801716, 0.0, -0.0013136618612236554, 0.0,
    -0.0014111853859280923, 0.0, -0.0015091578194443773, 0.0,
    -0.0016074977884102497, 0.0, -0.0017061386699515455, 0.0,
    -0.0018050259178723285, 0.0, -0.0019041148735245224, 0.0,
    -0.0020033689734995554, 0.0, -0.0021027582822103937, 0.0,
    -0.00220225829047132, 0.0, -0.0023018489318582555 };

  static const real_T b_Kv[186]{ 2.1048813364274074E-21, 0.0,
    2.0367999231504753E-21, 0.0, 1.963101434528115E-21, 0.0,
    1.9185804023058977E-21, 0.0, 1.868966093870934E-21, 0.0,
    1.796328274089597E-21, 0.0, 1.7349911122253774E-21, 0.0,
    1.7030514264606644E-21, 0.0, 1.6488191057351622E-21, 0.0,
    1.6067314995395007E-21, 0.0, 1.5421026063446421E-21, 0.0,
    1.4541645260462805E-21, 0.0, 1.3774160449568385E-21, 0.0,
    1.3124717843253368E-21, 0.0, 1.2423326250080741E-21, 0.0,
    1.1669575798210663E-21, 0.0, 1.0863130912963038E-21, 0.0,
    1.0003716849027306E-21, 0.0, 9.0911086639684616E-22, 0.0,
    8.1251221904972838E-22, 0.0, 7.45810987081722E-22, 0.0,
    6.750972627362992E-22, 0.0, 6.003613905675536E-22, 0.0,
    5.0397030408364459E-22, 0.0, 4.0286631237260206E-22, 0.0,
    3.1466905176009051E-22, 0.0, 2.40050378749355E-22, 0.0,
    1.7968319574083117E-22, 0.0, 1.1661610016180726E-22, 0.0,
    5.08471718523989E-23, 0.0, 0.0, 0.0, 1.9645136911104807E-21, 0.0,
    1.9014926490178118E-21, 0.0, 1.8334112357408797E-21, 0.0,
    1.7936108616867852E-21, 0.0, 1.749089829464568E-21, 0.0,
    1.6825264637454713E-21, 0.0, 1.6268377012482672E-21, 0.0,
    1.5993986539523138E-21, 0.0, 1.5505099109034678E-21, 0.0,
    1.5132266474620986E-21, 0.0, 1.4541899839823039E-21, 0.0,
    1.3726120335033124E-21, 0.0, 1.3016230104890837E-21, 0.0,
    1.2418235866837748E-21, 0.0, 1.1768793260522731E-21, 0.0,
    1.1067401667350104E-21, 0.0, 1.0313651215480027E-21, 0.0,
    9.5072063302324014E-22, 0.0, 8.64779226629667E-22, 0.0,
    7.7351840812378241E-22, 0.0, 7.1081787534493059E-22, 0.0,
    6.4411664337692421E-22, 0.0, 5.734029190315014E-22, 0.0,
    4.8171798957862287E-22, 0.0, 3.8532690309471381E-22, 0.0,
    3.0117196866780426E-22, 0.0, 2.299237653394257E-22, 0.0,
    1.7225414961282313E-22, 0.0, 1.1188696660429932E-22, 0.0,
    4.8819871025275388E-23, 0.0, 0.0, 0.0, 1.8276572918198257E-21, 0.0,
    1.7692342070093791E-21, 0.0, 1.7062131649167103E-21, 0.0,
    1.6706783323232952E-21, 0.0, 1.6308779582692007E-21, 0.0,
    1.5700836357052248E-21, 0.0, 1.5197935603278867E-21, 0.0,
    1.4966513785141996E-21, 0.0, 1.4529390408764875E-21, 0.0,
    1.4203235881694001E-21, 0.0, 1.3667670343862724E-21, 0.0,
    1.2914570805647193E-21, 0.0, 1.2261524204274862E-21, 0.0,
    1.171436687755016E-21, 0.0, 1.1116372639497072E-21, 0.0,
    1.0466930033182054E-21, 0.0, 9.7655384400094273E-22, 0.0,
    9.01178798813935E-22, 0.0, 8.2053431028917238E-22, 0.0,
    7.3459290389559928E-22, 0.0, 6.7587866607323165E-22, 0.0,
    6.1317813329437983E-22, 0.0, 5.4647690132637346E-22, 0.0,
    4.5948988663919219E-22, 0.0, 3.6780495718631362E-22, 0.0,
    2.8768716104416306E-22, 0.0, 2.1980551695901197E-22, 0.0,
    1.6483060397239186E-22, 0.0, 1.071609882457893E-22, 0.0,
    4.6793805237265468E-23, 0.0, 0.0, 0.0 };

  static const real_T b_Kx[12]{ 7.689433750374233E-6, -0.00052983264547518268,
    -4.0797700035767348E-7, 0.033752729948155841, 6.2634718616695827E-6,
    -0.00049099757092383356, -3.7807356318853408E-7, 0.03125149057206749,
    5.0972730132946579E-6, -0.00045360941458819537, -3.4928426905767272E-7,
    0.028849976794694789 };

  static const real_T b_Ku1[3]{ 0.0056257498156150479, 0.0053016816587513772,
    0.0049805833709413213 };

  static const real_T d_a[16]{ 0.81873075307798182, 1.5605438447738586E-5,
    -0.020266511909212012, 0.0, -0.090634596591614353, 0.99996475781523508,
    0.045768413341926564, 0.0, -6.97896423745501E-5, -2.7136872273192444E-8,
    1.0000352421847649, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T c_a[4]{ -1.1168314590257722, -0.63245448068765309,
    -0.0016088636636616162, 0.0 };

  static const real_T b_a[8]{ -4.0924793231964868E-18, -0.0027223994261342177,
    3.53553285779405, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T a[8]{ -0.28835920047413488, 0.59326526995080475,
    3.6969901567350294, 0.0089491794894280968, 0.0077934105482969848,
    -0.056577847513402488, -0.016886650855456559, 0.09468737576656984 };

  if ((&egoCar_M)->isMajorTimeStep()) {
    // set solver stop time
    rtsiSetSolverStopTime(&(&egoCar_M)->solverInfo,(((&egoCar_M)
      ->Timing.clockTick0+1)*(&egoCar_M)->Timing.stepSize0));
  }                                    // end MajorTimeStep

  // Update absolute time of base rate at minor time step
  if ((&egoCar_M)->isMinorTimeStep()) {
    (&egoCar_M)->Timing.t[0] = rtsiGetT(&(&egoCar_M)->solverInfo);
  }

  // Sum: '<S2>/Sum1' incorporates:
  //   Constant: '<Root>/v0 host'
  //   TransferFcn: '<S2>/Transfer Fcn'

  egoCar_B.Sum1 = 2.0 * egoCar_X.TransferFcn_CSTATE + 20.0;

  // MATLAB Function: '<S1>/DataTypeConversion_dmin' incorporates:
  //   Constant: '<Root>/Time gap'
  //   Constant: '<S1>/Default spacing constant'
  //   Product: '<S1>/Product2'
  //   Sum: '<S1>/Sum1'

  egoCar_DataTypeConversion_L0(egoCar_B.Sum1 * 1.4 + 10.0, &y_e);
  tmp = ((&egoCar_M)->isMajorTimeStep());

  // Sum: '<Root>/Sum' incorporates:
  //   Constant: '<Root>/x0 host'
  //   Inport: '<Root>/d_lead'
  //   Integrator: '<S2>/Integrator1'
  //   Sum: '<S2>/Sum'

  egoCar_Y.d_rel = egoCar_U.d_lead - (egoCar_X.Integrator1_CSTATE + 10.0);

  // MATLAB Function: '<S1>/DataTypeConversion_reldist'
  egoCar_DataTypeConversion_L0(egoCar_Y.d_rel, &y_h);

  // MATLAB Function: '<S1>/DataTypeConversion_vego'
  egoCar_DataTypeConversion_L0(egoCar_B.Sum1, &y_o);
  if (tmp) {
    // MATLAB Function: '<S1>/DataTypeConversion_L0' incorporates:
    //   Constant: '<S1>/Default spacing constant'

    egoCar_DataTypeConversion_L0(10.0, &rtb_y_n);

    // MATLAB Function: '<S1>/DataTypeConversion_vset' incorporates:
    //   Constant: '<Root>/Set velocity'

    egoCar_DataTypeConversion_L0(30.0, &rtb_y);
  }

  // Sum: '<Root>/Sum1' incorporates:
  //   Inport: '<Root>/v_lead'

  egoCar_Y.v_rel = egoCar_U.v_lead - egoCar_B.Sum1;

  // MATLAB Function: '<S1>/DataTypeConversion_vlead' incorporates:
  //   Sum: '<S1>/Sum6'

  egoCar_DataTypeConversion_L0(egoCar_B.Sum1 + egoCar_Y.v_rel, &y);
  if (tmp) {
    // MATLAB Function: '<S1>/DataTypeConversion_amin' incorporates:
    //   Constant: '<S1>/Minimum longitudinal acceleration constant'

    egoCar_DataTypeConversion_L0(-3.0, &rtb_y_f3);

    // MATLAB Function: '<S1>/DataTypeConversion_amax' incorporates:
    //   Constant: '<S1>/Maximum longitudinal acceleration constant'

    egoCar_DataTypeConversion_L0(2.0, &rtb_y_p);

    // SignalConversion generated from: '<S35>/ SFunction ' incorporates:
    //   Constant: '<S1>/Minimum velocity constant'
    //   MATLAB Function: '<S34>/optimizer'

    rtb_TmpSignalConversionAtSFun_f[0] = y_e;
    rtb_TmpSignalConversionAtSFun_f[1] = 0.0;

    // SignalConversion generated from: '<S35>/ SFunction ' incorporates:
    //   Constant: '<S1>/Maximum velocity constant'
    //   Constant: '<S1>/Unconstrained'
    //   MATLAB Function: '<S34>/optimizer'

    rtb_TmpSignalConversionAtSFun_c[0] = 0;
    rtb_TmpSignalConversionAtSFun_c[1] = 50;

    // MATLAB Function: '<S34>/optimizer' incorporates:
    //   SignalConversion generated from: '<S35>/ SFunction '

    std::memset(&vseq[0], 0, 62U * sizeof(real_T));
    for (i = 0; i < 31; i++) {
      vseq[(i << 1) + 1] = 1.0;
    }

    for (i = 0; i < 30; i++) {
      d_i = i << 1;
      rseq[d_i] = rtb_y_n * 0.02 - 0.76;
      rseq[d_i + 1] = rtb_y * 0.02 - 0.4;
    }

    for (i = 0; i < 31; i++) {
      vseq[i << 1] = egoCar_RMDscale * y - egoCar_voff;
    }

    y_e = vseq[0];
    y = vseq[1];
    xk[0] = egoCar_DW.last_x_PreviousInput[0];
    xk[1] = egoCar_DW.last_x_PreviousInput[1];
    xk[2] = egoCar_DW.last_x_PreviousInput[2];
    xk[3] = egoCar_DW.last_x_PreviousInput[3];

    // SignalConversion generated from: '<S35>/ SFunction ' incorporates:
    //   MATLAB Function: '<S34>/optimizer'

    ymax_incr[0] = y_h * 0.02 - 0.76;
    ymax_incr[1] = y_o * 0.02 - 0.4;

    // MATLAB Function: '<S34>/optimizer' incorporates:
    //   Memory: '<S14>/Memory'
    //   UnitDelay: '<S14>/last_mv'

    y_h = egoCar_DW.last_x_PreviousInput[1];
    y_o = egoCar_DW.last_x_PreviousInput[0];
    xk_0 = egoCar_DW.last_x_PreviousInput[2];
    xk_1 = egoCar_DW.last_x_PreviousInput[3];
    for (i = 0; i < 2; i++) {
      y_innov[i] = ymax_incr[i] - ((((f_a[i + 2] * y_h + f_a[i] * y_o) + f_a[i +
        4] * xk_0) + f_a[i + 6] * xk_1) + (0.0 * y_e + 0.0 * y));
    }

    y_innov_0 = y_innov[1];
    y_innov_1 = y_innov[0];
    for (i = 0; i < 4; i++) {
      rtb_xest[i] = (e_a[i + 4] * y_innov_0 + e_a[i] * y_innov_1) + xk[i];
    }

    ymax_incr_flag[0] = false;
    ymax_incr[0] = 0.0;
    ymin_incr_flag[0] = false;
    ymin_incr[0] = 0.0;
    ymax_incr_flag[1] = false;
    ymax_incr[1] = 0.0;
    ymin_incr_flag[1] = false;
    ymin_incr[1] = 0.0;
    umax_incr_flag = false;
    y_h = 0.0;
    umin_incr_flag = false;
    y_o = 0.0;
    for (d_i = 0; d_i < 96; d_i++) {
      xk_0 = 0.0;
      for (i = 0; i < 62; i++) {
        xk_0 += b_Mv[96 * i + d_i] * vseq[i];
      }

      xk_1 = b_Mlim[d_i];
      xk_0 = -((((((b_Mx[d_i + 96] * rtb_xest[1] + b_Mx[d_i] * rtb_xest[0]) +
                   b_Mx[d_i + 192] * rtb_xest[2]) + b_Mx[d_i + 288] * rtb_xest[3])
                 + xk_1) + b_Mu1[d_i] * egoCar_DW.last_mv_DSTATE) + xk_0);
      Bc[d_i] = xk_0;
      b_Mrows = b_Mrows_0[d_i];
      if (b_Mrows <= 60) {
        i = (b_Mrows - (((b_Mrows - 1) / egoCar_ny) << 1)) - 1;
        b_Del_Save_Flag0 = ymax_incr_flag[i];
        if (!ymax_incr_flag[i]) {
          xk_1 = -(0.02 * static_cast<real_T>(rtb_TmpSignalConversionAtSFun_c[i])
                   - (-0.36 * static_cast<real_T>(i) + 0.76)) - (-xk_1);
          b_Del_Save_Flag0 = true;
        } else {
          xk_1 = ymax_incr[i];
        }

        ymax_incr[i] = xk_1;
        ymax_incr_flag[i] = b_Del_Save_Flag0;
        Bc[d_i] = xk_0 + xk_1;
      } else if (b_Mrows <= 120) {
        i = (b_Mrows - (((b_Mrows - 61) >> 1) << 1)) - 61;
        b_Del_Save_Flag0 = ymin_incr_flag[i];
        if (!ymin_incr_flag[i]) {
          xk_1 = (0.02 * rtb_TmpSignalConversionAtSFun_f[i] - (-0.36 *
                   static_cast<real_T>(i) + 0.76)) - (-xk_1);
          b_Del_Save_Flag0 = true;
        } else {
          xk_1 = ymin_incr[i];
        }

        ymin_incr[i] = xk_1;
        ymin_incr_flag[i] = b_Del_Save_Flag0;
        Bc[d_i] = xk_0 + xk_1;
      } else if (b_Mrows <= 150) {
        if (!umax_incr_flag) {
          y_h = -(egoCar_RMVscale * rtb_y_p) - (-xk_1);
          umax_incr_flag = true;
        }

        Bc[d_i] = xk_0 + y_h;
      } else {
        if (!umin_incr_flag) {
          y_o = egoCar_RMVscale * rtb_y_f3 - (-xk_1);
          umin_incr_flag = true;
        }

        Bc[d_i] = xk_0 + y_o;
      }
    }

    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    f[3] = 0.0;
    for (d_i = 0; d_i < 3; d_i++) {
      rtb_y_f3 = 0.0;
      for (i = 0; i < 60; i++) {
        rtb_y_f3 += b_Kr[60 * d_i + i] * rseq[i];
      }

      rtb_y_p = 0.0;
      for (i = 0; i < 62; i++) {
        rtb_y_p += b_Kv[62 * d_i + i] * vseq[i];
      }

      i = d_i << 2;
      f[d_i] = (((((b_Kx[i + 1] * rtb_xest[1] + b_Kx[i] * rtb_xest[0]) + b_Kx[i
                   + 2] * rtb_xest[2]) + b_Kx[i + 3] * rtb_xest[3]) + rtb_y_f3)
                + b_Ku1[d_i] * egoCar_DW.last_mv_DSTATE) + rtb_y_p;
    }

    std::memcpy(&egoCar_B.iAout[0], &egoCar_DW.Memory_PreviousInput[0], 96U *
                sizeof(boolean_T));
    egoCar_qpkwik(b_Linv, b_Hinv, f, b_Ac, Bc, egoCar_B.iAout, 400, 1.0E-6,
                  rtb_xest, a__1, &i);
    if ((i < 0) || (i == 0)) {
      rtb_xest[0] = 0.0;
    }

    egoCar_B.u = egoCar_DW.last_mv_DSTATE + rtb_xest[0];
    y_h = egoCar_DW.last_x_PreviousInput[1];
    y_o = egoCar_DW.last_x_PreviousInput[0];
    xk_0 = egoCar_DW.last_x_PreviousInput[2];
    xk_1 = egoCar_DW.last_x_PreviousInput[3];
    y_innov_0 = y_innov[1];
    y_innov_1 = y_innov[0];
    for (i = 0; i < 4; i++) {
      egoCar_B.xk1[i] = (((((d_a[i + 4] * y_h + d_a[i] * y_o) + d_a[i + 8] *
                            xk_0) + d_a[i + 12] * xk_1) + c_a[i] * egoCar_B.u) +
                         (b_a[i] * y_e + 0.0 * y)) + (a[i + 4] * y_innov_0 + a[i]
        * y_innov_1);
    }

    // Outport: '<Root>/a_ego' incorporates:
    //   Gain: '<S14>/umin_scale1'

    egoCar_Y.a_ego = 5.0 * egoCar_B.u;
  }

  // Outport: '<Root>/v_ego'
  egoCar_Y.v_ego = egoCar_B.Sum1;
  if (tmp) {
    // MATLAB Function: '<S1>/DataTypeConversion_atrack' incorporates:
    //   Constant: '<S1>/External control signal constant'

    egoCar_DataTypeConversion_L0(0.0, &y_e);
  }

  // Integrator: '<S2>/Integrator'
  egoCar_B.Integrator = egoCar_X.Integrator_CSTATE;
  if ((&egoCar_M)->isMajorTimeStep()) {
    if ((&egoCar_M)->isMajorTimeStep()) {
      // Update for Memory: '<S14>/last_x'
      egoCar_DW.last_x_PreviousInput[0] = egoCar_B.xk1[0];
      egoCar_DW.last_x_PreviousInput[1] = egoCar_B.xk1[1];
      egoCar_DW.last_x_PreviousInput[2] = egoCar_B.xk1[2];
      egoCar_DW.last_x_PreviousInput[3] = egoCar_B.xk1[3];

      // Update for UnitDelay: '<S14>/last_mv'
      egoCar_DW.last_mv_DSTATE = egoCar_B.u;

      // Update for Memory: '<S14>/Memory'
      std::memcpy(&egoCar_DW.Memory_PreviousInput[0], &egoCar_B.iAout[0], 96U *
                  sizeof(boolean_T));
    }
  }                                    // end MajorTimeStep

  if ((&egoCar_M)->isMajorTimeStep()) {
    rt_ertODEUpdateContinuousStates(&(&egoCar_M)->solverInfo);

    // Update absolute time for base rate
    // The "clockTick0" counts the number of times the code of this task has
    //  been executed. The absolute time is the multiplication of "clockTick0"
    //  and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
    //  overflow during the application lifespan selected.

    ++(&egoCar_M)->Timing.clockTick0;
    (&egoCar_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&egoCar_M)->solverInfo);

    {
      // Update absolute timer for sample time: [0.1s, 0.0s]
      // The "clockTick1" counts the number of times the code of this task has
      //  been executed. The resolution of this integer timer is 0.1, which is the step size
      //  of the task. Size of "clockTick1" ensures timer will not overflow during the
      //  application lifespan selected.

      (&egoCar_M)->Timing.clockTick1++;
    }
  }                                    // end MajorTimeStep
}

// Derivatives for root system: '<Root>'
void egoCar::egoCar_derivatives()
{
  egoCar::XDot_egoCar_T *_rtXdot;
  _rtXdot = ((XDot_egoCar_T *) (&egoCar_M)->derivs);

  // Derivatives for TransferFcn: '<S2>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE = -2.0 * egoCar_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += egoCar_B.Integrator;

  // Derivatives for Integrator: '<S2>/Integrator1'
  _rtXdot->Integrator1_CSTATE = egoCar_B.Sum1;

  // Derivatives for Integrator: '<S2>/Integrator' incorporates:
  //   Outport: '<Root>/a_ego'

  _rtXdot->Integrator_CSTATE = egoCar_Y.a_ego;
}

// Model initialize function
void egoCar::initialize()
{
  // Registration code
  {
    // Setup solver object
    rtsiSetSimTimeStepPtr(&(&egoCar_M)->solverInfo, &(&egoCar_M)
                          ->Timing.simTimeStep);
    rtsiSetTPtr(&(&egoCar_M)->solverInfo, (&egoCar_M)->getTPtrPtr());
    rtsiSetStepSizePtr(&(&egoCar_M)->solverInfo, &(&egoCar_M)->Timing.stepSize0);
    rtsiSetdXPtr(&(&egoCar_M)->solverInfo, &(&egoCar_M)->derivs);
    rtsiSetContStatesPtr(&(&egoCar_M)->solverInfo, (real_T **) &(&egoCar_M)
                         ->contStates);
    rtsiSetNumContStatesPtr(&(&egoCar_M)->solverInfo, &(&egoCar_M)
      ->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&egoCar_M)->solverInfo, &(&egoCar_M)
      ->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&egoCar_M)->solverInfo, &(&egoCar_M)
      ->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&egoCar_M)->solverInfo, &(&egoCar_M)
      ->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&(&egoCar_M)->solverInfo, (boolean_T**)
      &(&egoCar_M)->contStateDisabled);
    rtsiSetErrorStatusPtr(&(&egoCar_M)->solverInfo, (&egoCar_M)
                          ->getErrorStatusPtr());
    rtsiSetRTModelPtr(&(&egoCar_M)->solverInfo, (&egoCar_M));
  }

  rtsiSetSimTimeStep(&(&egoCar_M)->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&(&egoCar_M)->solverInfo, false);
  rtsiSetIsContModeFrozen(&(&egoCar_M)->solverInfo, false);
  (&egoCar_M)->intgData.y = (&egoCar_M)->odeY;
  (&egoCar_M)->intgData.f[0] = (&egoCar_M)->odeF[0];
  (&egoCar_M)->intgData.f[1] = (&egoCar_M)->odeF[1];
  (&egoCar_M)->intgData.f[2] = (&egoCar_M)->odeF[2];
  (&egoCar_M)->contStates = ((X_egoCar_T *) &egoCar_X);
  (&egoCar_M)->contStateDisabled = ((XDis_egoCar_T *) &egoCar_XDis);
  (&egoCar_M)->Timing.tStart = (0.0);
  rtsiSetSolverData(&(&egoCar_M)->solverInfo, static_cast<void *>(&(&egoCar_M)
    ->intgData));
  rtsiSetSolverName(&(&egoCar_M)->solverInfo,"ode3");
  (&egoCar_M)->setTPtr(&(&egoCar_M)->Timing.tArray[0]);
  (&egoCar_M)->Timing.stepSize0 = 0.1;

  // InitializeConditions for TransferFcn: '<S2>/Transfer Fcn'
  egoCar_X.TransferFcn_CSTATE = 0.0;

  // InitializeConditions for Integrator: '<S2>/Integrator1'
  egoCar_X.Integrator1_CSTATE = 0.0;

  // InitializeConditions for Integrator: '<S2>/Integrator'
  egoCar_X.Integrator_CSTATE = 0.0;
}

// Model terminate function
void egoCar::terminate()
{
  // (no terminate code required)
}

time_T** egoCar::RT_MODEL_egoCar_T::getTPtrPtr()
{
  return &(Timing.t);
}

time_T* egoCar::RT_MODEL_egoCar_T::getTPtr() const
{
  return (Timing.t);
}

void egoCar::RT_MODEL_egoCar_T::setTPtr(time_T* aTPtr)
{
  (Timing.t = aTPtr);
}

boolean_T egoCar::RT_MODEL_egoCar_T::isMinorTimeStep() const
{
  return ((Timing.simTimeStep) == MINOR_TIME_STEP);
}

boolean_T egoCar::RT_MODEL_egoCar_T::getStopRequested() const
{
  return (Timing.stopRequestedFlag);
}

void egoCar::RT_MODEL_egoCar_T::setStopRequested(boolean_T aStopRequested)
{
  (Timing.stopRequestedFlag = aStopRequested);
}

boolean_T egoCar::RT_MODEL_egoCar_T::isMajorTimeStep() const
{
  return ((Timing.simTimeStep) == MAJOR_TIME_STEP);
}

boolean_T* egoCar::RT_MODEL_egoCar_T::getStopRequestedPtr()
{
  return (&(Timing.stopRequestedFlag));
}

const char_T** egoCar::RT_MODEL_egoCar_T::getErrorStatusPtr()
{
  return &errorStatus;
}

time_T egoCar::RT_MODEL_egoCar_T::getTStart() const
{
  return (Timing.tStart);
}

const char_T* egoCar::RT_MODEL_egoCar_T::getErrorStatus() const
{
  return (errorStatus);
}

void egoCar::RT_MODEL_egoCar_T::setErrorStatus(const char_T* const aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

// Constructor
egoCar::egoCar() :
  egoCar_U(),
  egoCar_Y(),
  egoCar_B(),
  egoCar_DW(),
  egoCar_X(),
  egoCar_XDis(),
  egoCar_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
egoCar::~egoCar() = default;

// Real-Time Model get method
egoCar::RT_MODEL_egoCar_T * egoCar::getRTM()
{
  return (&egoCar_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
