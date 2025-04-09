//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: automatic_transmission.h
//
// Code generated for Simulink model 'automatic_transmission'.
//
// Model version                  : 15.5
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Wed Apr  9 13:54:16 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Apple->ARM64
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
// Validation result: Not run
//
#ifndef automatic_transmission_h_
#define automatic_transmission_h_
#include <cmath>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include <cstring>
#ifndef ODE3_INTG
#define ODE3_INTG

// ODE3 Integration Data
struct ODE3_IntgData {
  real_T *y;                           // output
  real_T *f[3];                        // derivatives
};

#endif

extern "C"
{
  static real_T rtGetInf(void);
  static real32_T rtGetInfF(void);
  static real_T rtGetMinusInf(void);
  static real32_T rtGetMinusInfF(void);
}                                      // extern "C"

extern "C"
{
  static real_T rtGetNaN(void);
  static real32_T rtGetNaNF(void);
}                                      // extern "C"

extern "C"
{
  extern real_T rtInf;
  extern real_T rtMinusInf;
  extern real_T rtNaN;
  extern real32_T rtInfF;
  extern real32_T rtMinusInfF;
  extern real32_T rtNaNF;
  static boolean_T rtIsInf(real_T value);
  static boolean_T rtIsInfF(real32_T value);
  static boolean_T rtIsNaN(real_T value);
  static boolean_T rtIsNaNF(real32_T value);
}                                      // extern "C"

// Class declaration for model automatic_transmission
class automatic_transmission final
{
  // public data and function members
 public:
  // Block signals and states (default storage) for system '<Root>'
  struct DW {
    real_T VehicleSpeed;               // '<S5>/mph'
    real_T LookUpTable;                // '<S7>/Look-Up Table'
    real_T TransmissionRPM;            // '<S5>/FinalDriveRatio2'
    real_T engineimpellerinertia;      // '<S1>/engine + impeller inertia'
    real_T OutputTorque;               // '<S7>/Product'
    real_T VehicleInertia;             // '<S5>/Vehicle Inertia'
    real_T gear;                       // '<Root>/ShiftLogic'
    uint32_T temporalCounter_i1;       // '<Root>/ShiftLogic'
    uint8_T is_active_c1_automatic_transmis;// '<Root>/ShiftLogic'
    uint8_T is_active_gear_state;      // '<Root>/ShiftLogic'
    uint8_T is_gear_state;             // '<Root>/ShiftLogic'
    uint8_T is_active_selection_state; // '<Root>/ShiftLogic'
    uint8_T is_selection_state;        // '<Root>/ShiftLogic'
  };

  // Continuous states (default storage)
  struct X {
    real_T WheelSpeed_CSTATE;          // '<S5>/Wheel Speed'
    real_T Integrator_CSTATE;          // '<S1>/Integrator'
  };

  // State derivatives (default storage)
  struct XDot {
    real_T WheelSpeed_CSTATE;          // '<S5>/Wheel Speed'
    real_T Integrator_CSTATE;          // '<S1>/Integrator'
  };

  // State disabled
  struct XDis {
    boolean_T WheelSpeed_CSTATE;       // '<S5>/Wheel Speed'
    boolean_T Integrator_CSTATE;       // '<S1>/Integrator'
  };

  // Constant parameters (default storage)
  struct ConstP {
    // Expression: downth
    //  Referenced by: '<S3>/interp_down'

    real_T interp_down_RowIdx[6];

    // Pooled Parameter (Mixed Expressions)
    //  Referenced by:
    //    '<S3>/interp_down'
    //    '<S3>/interp_up'
    //    '<S7>/Look-Up Table'

    real_T pooled1[4];

    // Expression: downtab
    //  Referenced by: '<S3>/interp_down'

    real_T interp_down_Table[24];

    // Expression: upth
    //  Referenced by: '<S3>/interp_up'

    real_T interp_up_RowIdx[6];

    // Expression: uptab
    //  Referenced by: '<S3>/interp_up'

    real_T interp_up_Table[24];

    // Expression: thvec
    //  Referenced by: '<S1>/EngineTorque'

    real_T EngineTorque_RowIdx[10];

    // Expression: nevec
    //  Referenced by: '<S1>/EngineTorque'

    real_T EngineTorque_ColIdx[11];

    // Expression: emap
    //  Referenced by: '<S1>/EngineTorque'

    real_T EngineTorque_Table[110];

    // Expression: [2.393 1.450 1.000 0.677]
    //  Referenced by: '<S7>/Look-Up Table'

    real_T LookUpTable_YData[4];

    // Pooled Parameter (Expression: speedratio)
    //  Referenced by:
    //    '<S6>/FactorK'
    //    '<S6>/TorqueRatio'

    real_T pooled3[21];

    // Expression: Kfactor
    //  Referenced by: '<S6>/FactorK'

    real_T FactorK_YData[21];

    // Expression: Torkratio
    //  Referenced by: '<S6>/TorqueRatio'

    real_T TorqueRatio_YData[21];
  };

  // External inputs (root inport signals with default storage)
  struct ExtU {
    real_T throttle;                   // '<Root>/throttle'
    real_T brake;                      // '<Root>/brake'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY {
    real_T speed;                      // '<Root>/speed'
    real_T RPM;                        // '<Root>/RPM'
    real_T gear;                       // '<Root>/gear'
  };

  // Real-time Model Data Structure
  using odeFSubArray = real_T[2];
  struct RT_MODEL {
    const char_T *errorStatus;
    RTWSolverInfo solverInfo;
    X *contStates;
    int_T *periodicContStateIndices;
    real_T *periodicContStateRanges;
    real_T *derivs;
    XDis *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T CTOutputIncnstWithState;
    real_T odeY[2];
    real_T odeF[3][2];
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
      struct {
        uint8_T TID[3];
      } TaskCounters;

      time_T tStart;
      SimTimeStep simTimeStep;
      boolean_T stopRequestedFlag;
      time_T *t;
      time_T tArray[3];
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
  automatic_transmission(automatic_transmission const&) = delete;

  // Assignment Operator
  automatic_transmission& operator= (automatic_transmission const&) & = delete;

  // Move Constructor
  automatic_transmission(automatic_transmission &&) = delete;

  // Move Assignment Operator
  automatic_transmission& operator= (automatic_transmission &&) = delete;

  // Real-Time Model get method
  automatic_transmission::RT_MODEL * getRTM();

  // External inputs
  ExtU rtU;

  // External outputs
  ExtY rtY;

  // model initialize function
  void initialize();

  // model step function
  void step();

  // Constructor
  automatic_transmission();

  // Destructor
  ~automatic_transmission();

  // private data and function members
 private:
  // Block states
  DW rtDW;

  // Block continuous states
  X rtX;

  // Block Continuous state disabled vector
  XDis rtXDis;

  // private member function(s) for subsystem '<Root>'
  void gear_state(const int32_T *sfEvent);

  // Global mass matrix

  // Continuous states update member function
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  // Derivatives member function
  void automatic_transmission_derivatives();

  // Real-Time Model
  RT_MODEL rtM;
};

// Constant parameters (default storage)
extern const automatic_transmission::ConstP rtConstP;

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<Root>/Scope2' : Unused code path elimination
//  Block '<Root>/Scope1' : Unused code path elimination


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
//  '<Root>' : 'automatic_transmission'
//  '<S1>'   : 'automatic_transmission/Engine'
//  '<S2>'   : 'automatic_transmission/ShiftLogic'
//  '<S3>'   : 'automatic_transmission/ThresholdCalculation'
//  '<S4>'   : 'automatic_transmission/Transmission'
//  '<S5>'   : 'automatic_transmission/Vehicle'
//  '<S6>'   : 'automatic_transmission/Transmission/TorqueConverter'
//  '<S7>'   : 'automatic_transmission/Transmission/TransmissionRatio'

#endif                                 // automatic_transmission_h_

//
// File trailer for generated code.
//
// [EOF]
//
