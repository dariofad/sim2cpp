//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: egoCar.h
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
#ifndef egoCar_h_
#define egoCar_h_
#include <cmath>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_nonfinite.h"
#include "egoCar_types.h"

extern "C"
{

#include "rtGetNaN.h"

}

#include <cstring>
#ifndef ODE3_INTG
#define ODE3_INTG

// ODE3 Integration Data
struct ODE3_IntgData {
  real_T *y;                           // output
  real_T *f[3];                        // derivatives
};

#endif

// Class declaration for model egoCar
class egoCar final
{
  // public data and function members
 public:
  // Block signals (default storage)
  struct B_egoCar_T {
    real_T Sum1;                       // '<S2>/Sum1'
    real_T Integrator;                 // '<S2>/Integrator'
    real_T xk1[4];                     // '<S34>/optimizer'
    real_T u;                          // '<S34>/optimizer'
    boolean_T iAout[96];               // '<S34>/optimizer'
  };

  // Block states (default storage) for system '<Root>'
  struct DW_egoCar_T {
    real_T last_mv_DSTATE;             // '<S14>/last_mv'
    real_T last_x_PreviousInput[4];    // '<S14>/last_x'
    boolean_T Memory_PreviousInput[96];// '<S14>/Memory'
  };

  // Continuous states (default storage)
  struct X_egoCar_T {
    real_T TransferFcn_CSTATE;         // '<S2>/Transfer Fcn'
    real_T Integrator1_CSTATE;         // '<S2>/Integrator1'
    real_T Integrator_CSTATE;          // '<S2>/Integrator'
  };

  // State derivatives (default storage)
  struct XDot_egoCar_T {
    real_T TransferFcn_CSTATE;         // '<S2>/Transfer Fcn'
    real_T Integrator1_CSTATE;         // '<S2>/Integrator1'
    real_T Integrator_CSTATE;          // '<S2>/Integrator'
  };

  // State disabled
  struct XDis_egoCar_T {
    boolean_T TransferFcn_CSTATE;      // '<S2>/Transfer Fcn'
    boolean_T Integrator1_CSTATE;      // '<S2>/Integrator1'
    boolean_T Integrator_CSTATE;       // '<S2>/Integrator'
  };

  // Invariant block signals (default storage)
  struct ConstB_egoCar_T {
    real_T MathFunction[2];            // '<S14>/Math Function'
    real_T MathFunction1;              // '<S14>/Math Function1'
    real_T MathFunction2;              // '<S14>/Math Function2'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU_egoCar_T {
    real_T d_lead;                     // '<Root>/d_lead'
    real_T v_lead;                     // '<Root>/v_lead'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_egoCar_T {
    real_T d_rel;                      // '<Root>/d_rel'
    real_T v_rel;                      // '<Root>/v_rel'
    real_T v_ego;                      // '<Root>/v_ego'
    real_T a_ego;                      // '<Root>/a_ego'
  };

  // Real-time Model Data Structure
  using odeFSubArray = real_T[3];
  struct RT_MODEL_egoCar_T {
    const char_T *errorStatus;
    RTWSolverInfo solverInfo;
    X_egoCar_T *contStates;
    int_T *periodicContStateIndices;
    real_T *periodicContStateRanges;
    real_T *derivs;
    XDis_egoCar_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T CTOutputIncnstWithState;
    real_T odeY[3];
    real_T odeF[3][3];
    ODE3_IntgData intgData;

    //
    //  Sizes:
    //  The following substructure contains sizes information
    //  for many of the model attributes such as inputs, outputs,
    //  dwork, sample times, etc.

    struct {
      int_T numContStates;
      int_T numPeriodicContStates;
      int_T numSampTimes;
    } Sizes;

    //
    //  Timing:
    //  The following substructure contains information regarding
    //  the timing information for the model.

    struct {
      uint32_T clockTick0;
      time_T stepSize0;
      uint32_T clockTick1;
      time_T tStart;
      SimTimeStep simTimeStep;
      boolean_T stopRequestedFlag;
      time_T *t;
      time_T tArray[2];
    } Timing;

    time_T** getTPtrPtr();
    time_T* getTPtr() const;
    void setTPtr(time_T* aTPtr);
    boolean_T isMinorTimeStep() const;
    boolean_T getStopRequested() const;
    void setStopRequested(boolean_T aStopRequested);
    boolean_T isMajorTimeStep() const;
    boolean_T* getStopRequestedPtr();
    const char_T** getErrorStatusPtr();
    time_T getTStart() const;
    const char_T* getErrorStatus() const;
    void setErrorStatus(const char_T* const aErrorStatus);
  };

  // Copy Constructor
  egoCar(egoCar const&) = delete;

  // Assignment Operator
  egoCar& operator= (egoCar const&) & = delete;

  // Move Constructor
  egoCar(egoCar &&) = delete;

  // Move Assignment Operator
  egoCar& operator= (egoCar &&) = delete;

  // Real-Time Model get method
  egoCar::RT_MODEL_egoCar_T * getRTM();

  // Root inports set method
  void setExternalInputs(const ExtU_egoCar_T *pExtU_egoCar_T)
  {
    egoCar_U = *pExtU_egoCar_T;
  }

  // Root outports get method
  const ExtY_egoCar_T &getExternalOutputs() const
  {
    return egoCar_Y;
  }

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  egoCar();

  // Destructor
  ~egoCar();

  // private data and function members
 private:
  // External inputs
  ExtU_egoCar_T egoCar_U;

  // External outputs
  ExtY_egoCar_T egoCar_Y;

  // Block signals
  B_egoCar_T egoCar_B;

  // Block states
  DW_egoCar_T egoCar_DW;

  // Block continuous states
  X_egoCar_T egoCar_X;

  // Block Continuous state disabled vector
  XDis_egoCar_T egoCar_XDis;

  // private member function(s) for subsystem '<S1>/DataTypeConversion_L0'
  static void egoCar_DataTypeConversion_L0(real_T rtu_u, real_T *rty_y);

  // private member function(s) for subsystem '<Root>'
  real_T egoCar_norm(const real_T x[4]);
  real_T egoCar_maximum(const real_T x[4]);
  real_T egoCar_xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  void egoCar_xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0,
                    const real_T x[16], int32_T ix0, real_T y[4]);
  void egoCar_xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const
                    real_T y[4], real_T b_A[16], int32_T ia0);
  void egoCar_KWIKfactor(const real_T b_Ac[384], const int32_T iC[96], int32_T
    nA, const real_T b_Linv[16], real_T D[16], real_T b_H[16], int32_T n, real_T
    RLinv[16], real_T *Status);
  void egoCar_DropConstraint(int32_T kDrop, boolean_T iA[96], int32_T *nA,
    int32_T iC[96]);
  void egoCar_qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16], const
                     real_T f[4], const real_T b_Ac[384], const real_T b[96],
                     boolean_T iA[96], int32_T maxiter, real_T FeasTol, real_T
                     x[4], real_T lambda[96], int32_T *status);

  // Global mass matrix

  // Continuous states update member function
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  // Derivatives member function
  void egoCar_derivatives();

  // Real-Time Model
  RT_MODEL_egoCar_T egoCar_M;
};

extern const egoCar::ConstB_egoCar_T egoCar_ConstB;// constant block i/o

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S14>/Constant' : Unused code path elimination
//  Block '<S14>/Floor' : Unused code path elimination
//  Block '<S14>/Floor1' : Unused code path elimination
//  Block '<S15>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S16>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S17>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S18>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S19>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S20>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S21>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S22>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S23>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S24>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S25>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S26>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S27>/Vector Dimension Check' : Unused code path elimination
//  Block '<S28>/Vector Dimension Check' : Unused code path elimination
//  Block '<S29>/Vector Dimension Check' : Unused code path elimination
//  Block '<S30>/Vector Dimension Check' : Unused code path elimination
//  Block '<S31>/Vector Dimension Check' : Unused code path elimination
//  Block '<S32>/Vector Dimension Check' : Unused code path elimination
//  Block '<S14>/Min' : Unused code path elimination
//  Block '<S14>/constant' : Unused code path elimination
//  Block '<S33>/Vector Dimension Check' : Unused code path elimination
//  Block '<S14>/umin_scale2' : Unused code path elimination
//  Block '<S14>/umin_scale3' : Unused code path elimination
//  Block '<S14>/umin_scale5' : Unused code path elimination
//  Block '<S14>/ym_zero' : Unused code path elimination
//  Block '<S13>/m_zero' : Unused code path elimination
//  Block '<S13>/p_zero' : Unused code path elimination
//  Block '<S14>/Reshape' : Reshape block reduction
//  Block '<S14>/Reshape1' : Reshape block reduction
//  Block '<S14>/Reshape2' : Reshape block reduction
//  Block '<S14>/Reshape3' : Reshape block reduction
//  Block '<S14>/Reshape4' : Reshape block reduction
//  Block '<S14>/Reshape5' : Reshape block reduction


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'egoCar'
//  '<S1>'   : 'egoCar/Adaptive Cruise Control System'
//  '<S2>'   : 'egoCar/Ego Car'
//  '<S3>'   : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_L0'
//  '<S4>'   : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_amax'
//  '<S5>'   : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_amin'
//  '<S6>'   : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_atrack'
//  '<S7>'   : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_dmin'
//  '<S8>'   : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_optsgn'
//  '<S9>'   : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_reldist'
//  '<S10>'  : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_vego'
//  '<S11>'  : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_vlead'
//  '<S12>'  : 'egoCar/Adaptive Cruise Control System/DataTypeConversion_vset'
//  '<S13>'  : 'egoCar/Adaptive Cruise Control System/MPC'
//  '<S14>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC'
//  '<S15>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check'
//  '<S16>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check1'
//  '<S17>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check2'
//  '<S18>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check'
//  '<S19>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check1'
//  '<S20>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check2'
//  '<S21>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check3'
//  '<S22>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check4'
//  '<S23>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check5'
//  '<S24>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check6'
//  '<S25>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check7'
//  '<S26>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check8'
//  '<S27>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check'
//  '<S28>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check1'
//  '<S29>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check2'
//  '<S30>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check'
//  '<S31>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check1'
//  '<S32>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check6'
//  '<S33>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/moorx'
//  '<S34>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/optimizer'
//  '<S35>'  : 'egoCar/Adaptive Cruise Control System/MPC/MPC/optimizer/optimizer'

#endif                                 // egoCar_h_

//
// File trailer for generated code.
//
// [EOF]
//
