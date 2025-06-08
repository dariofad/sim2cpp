//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: mpcACCsystem.h
//
// Code generated for Simulink model 'mpcACCsystem'.
//
// Model version                  : 14.5
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Sun Jun  8 15:50:16 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Apple->ARM64
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef mpcACCsystem_h_
#define mpcACCsystem_h_
#include <cmath>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_nonfinite.h"
#include "mpcACCsystem_types.h"

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

// Class declaration for model mpcACCsystem
class mpcACCsystem final
{
  // public data and function members
 public:
  // Block signals (default storage)
  struct B_mpcACCsystem_T {
    real_T Sum1;                       // '<S2>/Sum1'
    real_T Sum1_m;                     // '<S3>/Sum1'
    real_T Integrator;                 // '<S2>/Integrator'
    real_T Integrator_c;               // '<S3>/Integrator'
    real_T xk1[4];                     // '<S35>/optimizer'
    real_T u;                          // '<S35>/optimizer'
    boolean_T iAout[96];               // '<S35>/optimizer'
  };

  // Block states (default storage) for system '<Root>'
  struct DW_mpcACCsystem_T {
    real_T last_mv_DSTATE;             // '<S15>/last_mv'
    real_T last_x_PreviousInput[4];    // '<S15>/last_x'
    real_T lastSin;                    // '<Root>/Sine Wave'
    real_T lastCos;                    // '<Root>/Sine Wave'
    int32_T systemEnable;              // '<Root>/Sine Wave'
    boolean_T Memory_PreviousInput[96];// '<S15>/Memory'
  };

  // Continuous states (default storage)
  struct X_mpcACCsystem_T {
    real_T TransferFcn_CSTATE;         // '<S2>/Transfer Fcn'
    real_T Integrator1_CSTATE;         // '<S3>/Integrator1'
    real_T Integrator1_CSTATE_p;       // '<S2>/Integrator1'
    real_T TransferFcn_CSTATE_i;       // '<S3>/Transfer Fcn'
    real_T Integrator_CSTATE;          // '<S2>/Integrator'
    real_T Integrator_CSTATE_a;        // '<S3>/Integrator'
  };

  // State derivatives (default storage)
  struct XDot_mpcACCsystem_T {
    real_T TransferFcn_CSTATE;         // '<S2>/Transfer Fcn'
    real_T Integrator1_CSTATE;         // '<S3>/Integrator1'
    real_T Integrator1_CSTATE_p;       // '<S2>/Integrator1'
    real_T TransferFcn_CSTATE_i;       // '<S3>/Transfer Fcn'
    real_T Integrator_CSTATE;          // '<S2>/Integrator'
    real_T Integrator_CSTATE_a;        // '<S3>/Integrator'
  };

  // State disabled
  struct XDis_mpcACCsystem_T {
    boolean_T TransferFcn_CSTATE;      // '<S2>/Transfer Fcn'
    boolean_T Integrator1_CSTATE;      // '<S3>/Integrator1'
    boolean_T Integrator1_CSTATE_p;    // '<S2>/Integrator1'
    boolean_T TransferFcn_CSTATE_i;    // '<S3>/Transfer Fcn'
    boolean_T Integrator_CSTATE;       // '<S2>/Integrator'
    boolean_T Integrator_CSTATE_a;     // '<S3>/Integrator'
  };

  // Invariant block signals (default storage)
  struct ConstB_mpcACCsystem_T {
    real_T MathFunction[2];            // '<S15>/Math Function'
    real_T MathFunction1;              // '<S15>/Math Function1'
    real_T MathFunction2;              // '<S15>/Math Function2'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_mpcACCsystem_T {
    real_T a_lead;                     // '<Root>/a_lead'
    real_T d_rel;                      // '<Root>/d_rel'
    real_T v_rel;                      // '<Root>/v_rel'
    real_T v_ego;                      // '<Root>/v_ego'
    real_T a_ego;                      // '<Root>/a_ego'
  };

  // Real-time Model Data Structure
  using odeFSubArray = real_T[6];
  struct RT_MODEL_mpcACCsystem_T {
    const char_T *errorStatus;
    RTWSolverInfo solverInfo;
    X_mpcACCsystem_T *contStates;
    int_T *periodicContStateIndices;
    real_T *periodicContStateRanges;
    real_T *derivs;
    XDis_mpcACCsystem_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T CTOutputIncnstWithState;
    real_T odeY[6];
    real_T odeF[3][6];
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
  mpcACCsystem(mpcACCsystem const&) = delete;

  // Assignment Operator
  mpcACCsystem& operator= (mpcACCsystem const&) & = delete;

  // Move Constructor
  mpcACCsystem(mpcACCsystem &&) = delete;

  // Move Assignment Operator
  mpcACCsystem& operator= (mpcACCsystem &&) = delete;

  // Real-Time Model get method
  mpcACCsystem::RT_MODEL_mpcACCsystem_T * getRTM();

  // Root outports get method
  const ExtY_mpcACCsystem_T &getExternalOutputs() const
  {
    return mpcACCsystem_Y;
  }

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  mpcACCsystem();

  // Destructor
  ~mpcACCsystem();

  // private data and function members
 private:
  // External outputs
  ExtY_mpcACCsystem_T mpcACCsystem_Y;

  // Block signals
  B_mpcACCsystem_T mpcACCsystem_B;

  // Block states
  DW_mpcACCsystem_T mpcACCsystem_DW;

  // Block continuous states
  X_mpcACCsystem_T mpcACCsystem_X;

  // Block Continuous state disabled vector
  XDis_mpcACCsystem_T mpcACCsystem_XDis;

  // private member function(s) for subsystem '<S1>/DataTypeConversion_L0'
  static void mpcACCsys_DataTypeConversion_L0(real_T rtu_u, real_T *rty_y);

  // private member function(s) for subsystem '<Root>'
  real_T mpcACCsystem_norm(const real_T x[4]);
  real_T mpcACCsystem_maximum(const real_T x[4]);
  real_T mpcACCsystem_xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  void mpcACCsystem_xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T
    ia0, const real_T x[16], int32_T ix0, real_T y[4]);
  void mpcACCsystem_xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0,
    const real_T y[4], real_T b_A[16], int32_T ia0);
  void mpcACCsystem_KWIKfactor(const real_T b_Ac[384], const int32_T iC[96],
    int32_T nA, const real_T b_Linv[16], real_T D[16], real_T b_H[16], int32_T n,
    real_T RLinv[16], real_T *Status);
  void mpcACCsystem_DropConstraint(int32_T kDrop, boolean_T iA[96], int32_T *nA,
    int32_T iC[96]);
  void mpcACCsystem_qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16],
    const real_T f[4], const real_T b_Ac[384], const real_T b[96], boolean_T iA
    [96], int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[96],
    int32_T *status);

  // Global mass matrix

  // Continuous states update member function
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  // Derivatives member function
  void mpcACCsystem_derivatives();

  // Real-Time Model
  RT_MODEL_mpcACCsystem_T mpcACCsystem_M;
};

extern const mpcACCsystem::ConstB_mpcACCsystem_T mpcACCsystem_ConstB;// constant block i/o 

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S15>/Constant' : Unused code path elimination
//  Block '<S15>/Floor' : Unused code path elimination
//  Block '<S15>/Floor1' : Unused code path elimination
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
//  Block '<S27>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S28>/Vector Dimension Check' : Unused code path elimination
//  Block '<S29>/Vector Dimension Check' : Unused code path elimination
//  Block '<S30>/Vector Dimension Check' : Unused code path elimination
//  Block '<S31>/Vector Dimension Check' : Unused code path elimination
//  Block '<S32>/Vector Dimension Check' : Unused code path elimination
//  Block '<S33>/Vector Dimension Check' : Unused code path elimination
//  Block '<S15>/Min' : Unused code path elimination
//  Block '<S15>/constant' : Unused code path elimination
//  Block '<S34>/Vector Dimension Check' : Unused code path elimination
//  Block '<S15>/umin_scale2' : Unused code path elimination
//  Block '<S15>/umin_scale3' : Unused code path elimination
//  Block '<S15>/umin_scale5' : Unused code path elimination
//  Block '<S15>/ym_zero' : Unused code path elimination
//  Block '<S14>/m_zero' : Unused code path elimination
//  Block '<S14>/p_zero' : Unused code path elimination
//  Block '<S15>/Reshape' : Reshape block reduction
//  Block '<S15>/Reshape1' : Reshape block reduction
//  Block '<S15>/Reshape2' : Reshape block reduction
//  Block '<S15>/Reshape3' : Reshape block reduction
//  Block '<S15>/Reshape4' : Reshape block reduction
//  Block '<S15>/Reshape5' : Reshape block reduction


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
//  '<Root>' : 'mpcACCsystem'
//  '<S1>'   : 'mpcACCsystem/Adaptive Cruise Control System'
//  '<S2>'   : 'mpcACCsystem/Ego Car'
//  '<S3>'   : 'mpcACCsystem/Lead Car'
//  '<S4>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_L0'
//  '<S5>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_amax'
//  '<S6>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_amin'
//  '<S7>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_atrack'
//  '<S8>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_dmin'
//  '<S9>'   : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_optsgn'
//  '<S10>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_reldist'
//  '<S11>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vego'
//  '<S12>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vlead'
//  '<S13>'  : 'mpcACCsystem/Adaptive Cruise Control System/DataTypeConversion_vset'
//  '<S14>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC'
//  '<S15>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC'
//  '<S16>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check'
//  '<S17>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check1'
//  '<S18>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Matrix Signal Check2'
//  '<S19>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check'
//  '<S20>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check1'
//  '<S21>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check2'
//  '<S22>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check3'
//  '<S23>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check4'
//  '<S24>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check5'
//  '<S25>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check6'
//  '<S26>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check7'
//  '<S27>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Preview Signal Check8'
//  '<S28>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check'
//  '<S29>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check1'
//  '<S30>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Scalar Signal Check2'
//  '<S31>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check'
//  '<S32>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check1'
//  '<S33>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/MPC Vector Signal Check6'
//  '<S34>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/moorx'
//  '<S35>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/optimizer'
//  '<S36>'  : 'mpcACCsystem/Adaptive Cruise Control System/MPC/MPC/optimizer/optimizer'

#endif                                 // mpcACCsystem_h_

//
// File trailer for generated code.
//
// [EOF]
//
