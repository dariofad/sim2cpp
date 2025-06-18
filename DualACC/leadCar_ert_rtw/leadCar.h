//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: leadCar.h
//
// Code generated for Simulink model 'leadCar'.
//
// Model version                  : 1.4
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Wed Jun 18 12:30:29 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Apple->ARM64
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef leadCar_h_
#define leadCar_h_
#include <cmath>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "leadCar_types.h"
#include <cstring>
#ifndef ODE3_INTG
#define ODE3_INTG

// ODE3 Integration Data
struct ODE3_IntgData {
  real_T *y;                           // output
  real_T *f[3];                        // derivatives
};

#endif

// Class declaration for model leadCar
class leadCar final
{
  // public data and function members
 public:
  // Block signals (default storage)
  struct B_leadCar_T {
    real_T Sum1;                       // '<S1>/Sum1'
    real_T Integrator;                 // '<S1>/Integrator'
  };

  // Block states (default storage) for system '<Root>'
  struct DW_leadCar_T {
    real_T lastSin;                    // '<Root>/Sine Wave'
    real_T lastCos;                    // '<Root>/Sine Wave'
    int32_T systemEnable;              // '<Root>/Sine Wave'
  };

  // Continuous states (default storage)
  struct X_leadCar_T {
    real_T Integrator1_CSTATE;         // '<S1>/Integrator1'
    real_T TransferFcn_CSTATE;         // '<S1>/Transfer Fcn'
    real_T Integrator_CSTATE;          // '<S1>/Integrator'
  };

  // State derivatives (default storage)
  struct XDot_leadCar_T {
    real_T Integrator1_CSTATE;         // '<S1>/Integrator1'
    real_T TransferFcn_CSTATE;         // '<S1>/Transfer Fcn'
    real_T Integrator_CSTATE;          // '<S1>/Integrator'
  };

  // State disabled
  struct XDis_leadCar_T {
    boolean_T Integrator1_CSTATE;      // '<S1>/Integrator1'
    boolean_T TransferFcn_CSTATE;      // '<S1>/Transfer Fcn'
    boolean_T Integrator_CSTATE;       // '<S1>/Integrator'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_leadCar_T {
    real_T a_lead;                     // '<Root>/a_lead'
    real_T d_lead;                     // '<Root>/d_lead'
    real_T v_lead;                     // '<Root>/v_lead'
  };

  // Real-time Model Data Structure
  using odeFSubArray = real_T[3];
  struct RT_MODEL_leadCar_T {
    const char_T *errorStatus;
    RTWSolverInfo solverInfo;
    X_leadCar_T *contStates;
    int_T *periodicContStateIndices;
    real_T *periodicContStateRanges;
    real_T *derivs;
    XDis_leadCar_T *contStateDisabled;
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
  leadCar(leadCar const&) = delete;

  // Assignment Operator
  leadCar& operator= (leadCar const&) & = delete;

  // Move Constructor
  leadCar(leadCar &&) = delete;

  // Move Assignment Operator
  leadCar& operator= (leadCar &&) = delete;

  // Real-Time Model get method
  leadCar::RT_MODEL_leadCar_T * getRTM();

  // Root outports get method
  const ExtY_leadCar_T &getExternalOutputs() const
  {
    return leadCar_Y;
  }

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  leadCar();

  // Destructor
  ~leadCar();

  // private data and function members
 private:
  // External outputs
  ExtY_leadCar_T leadCar_Y;

  // Block signals
  B_leadCar_T leadCar_B;

  // Block states
  DW_leadCar_T leadCar_DW;

  // Block continuous states
  X_leadCar_T leadCar_X;

  // Block Continuous state disabled vector
  XDis_leadCar_T leadCar_XDis;

  // Global mass matrix

  // Continuous states update member function
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  // Derivatives member function
  void leadCar_derivatives();

  // Real-Time Model
  RT_MODEL_leadCar_T leadCar_M;
};

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
//  '<Root>' : 'leadCar'
//  '<S1>'   : 'leadCar/Lead Car'

#endif                                 // leadCar_h_

//
// File trailer for generated code.
//
// [EOF]
//
