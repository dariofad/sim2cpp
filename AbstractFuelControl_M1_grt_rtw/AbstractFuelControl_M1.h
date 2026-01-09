/*
 * AbstractFuelControl_M1.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "AbstractFuelControl_M1".
 *
 * Model version              : 21.0
 * Simulink Coder version : 25.1 (R2025a) 21-Nov-2024
 * C++ source code generated on : Fri Jan  9 07:42:46 2026
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: AMD->x86-64 (Linux 64)
 * Code generation objective: Debugging
 * Validation result: Not run
 */

#ifndef AbstractFuelControl_M1_h_
#define AbstractFuelControl_M1_h_
#include <cmath>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#include "rt_nonfinite.h"
#include "AbstractFuelControl_M1_types.h"
#include "rt_zcfcn.h"

extern "C"
{

#include "rtGetNaN.h"

}

#include <cfloat>
#include <cstring>
#include "zero_crossing_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

/* Block signals (default storage) */
struct B_AbstractFuelControl_M1_T {
  real_T Gain;                         /* '<S19>/Gain' */
  real_T fuelinjectortolerance095105;
                                /* '<S1>/fuel injector tolerance [0.95 1.05]' */
  real_T airbyfuel;                    /* '<S3>/Divide' */
  real_T Gain1;                        /* '<S18>/Gain1' */
  real_T delays;                       /* '<S3>/delay (s)' */
  real_T RTVm;                         /* '<S4>/RT//Vm' */
  real_T Sum;                          /* '<S6>/Sum' */
  real_T Pwon;                         /* '<S2>/Pwon' */
  real_T PulseGenerator_10ms;          /* '<S2>/PulseGenerator_10ms' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_AbstractFuelControl_M1_T {
  real_T fuelsystemtransportdelay_RWORK[61441];/* '<S3>/fuel system transport delay' */
  void *fuelsystemtransportdelay_PWORK[3];/* '<S3>/fuel system transport delay' */
  real32_T UnitDelay2_DSTATE;          /* '<S14>/Unit Delay2' */
  real32_T UnitDelay1_DSTATE;          /* '<S12>/UnitDelay1' */
  real32_T UnitDelay1_DSTATE_d;        /* '<S11>/UnitDelay1' */
  real32_T commanded_fuel;             /* '<S2>/commanded_fuel' */
  real32_T airbyfuel_ref;              /* '<S2>/mode_fb1' */
  real32_T engine_speed;               /* '<S7>/DataStoreMemory' */
  real32_T throttle_flow;              /* '<S7>/DataStoreMemory1' */
  real32_T airbyfuel_meas;             /* '<S7>/DataStoreMemory2' */
  real32_T throttle_angle;             /* '<S7>/DataStoreMemory3' */
  int32_T clockTickCounter;            /* '<S2>/PulseGenerator_10ms' */
  int_T fuelsystemtransportdelay_IWORK[4];/* '<S3>/fuel system transport delay' */
  boolean_T UnitDelay_DSTATE;          /* '<S16>/Unit Delay' */
  boolean_T UnitDelay1_DSTATE_a;       /* '<S15>/Unit Delay1' */
  boolean_T UnitDelay1_DSTATE_e;       /* '<S14>/Unit Delay1' */
  boolean_T controller_mode;           /* '<S2>/mode_fb' */
};

/* Continuous states (default storage) */
struct X_AbstractFuelControl_M1_T {
  real_T Integrator_CSTATE;            /* '<S19>/Integrator' */
  real_T Throttledelay_CSTATE;         /* '<S1>/Throttle delay' */
  real_T p00543bar_CSTATE;             /* '<S4>/p0 = 0.543 (bar)' */
  real_T Integrator_CSTATE_h;          /* '<S18>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S6>/Integrator' */
  real_T fuelsystemtransportdelay_CSTATE;/* '<S3>/fuel system transport delay' */
};

/* State derivatives (default storage) */
struct XDot_AbstractFuelControl_M1_T {
  real_T Integrator_CSTATE;            /* '<S19>/Integrator' */
  real_T Throttledelay_CSTATE;         /* '<S1>/Throttle delay' */
  real_T p00543bar_CSTATE;             /* '<S4>/p0 = 0.543 (bar)' */
  real_T Integrator_CSTATE_h;          /* '<S18>/Integrator' */
  real_T Integrator_CSTATE_c;          /* '<S6>/Integrator' */
  real_T fuelsystemtransportdelay_CSTATE;/* '<S3>/fuel system transport delay' */
};

/* State disabled  */
struct XDis_AbstractFuelControl_M1_T {
  boolean_T Integrator_CSTATE;         /* '<S19>/Integrator' */
  boolean_T Throttledelay_CSTATE;      /* '<S1>/Throttle delay' */
  boolean_T p00543bar_CSTATE;          /* '<S4>/p0 = 0.543 (bar)' */
  boolean_T Integrator_CSTATE_h;       /* '<S18>/Integrator' */
  boolean_T Integrator_CSTATE_c;       /* '<S6>/Integrator' */
  boolean_T fuelsystemtransportdelay_CSTATE;/* '<S3>/fuel system transport delay' */
};

/* Zero-crossing (trigger) state */
struct PrevZCX_AbstractFuelControl_M_T {
  ZCSigState fuel_controller_pwon_Trig_ZCE;/* '<S7>/fuel_controller_pwon' */
  ZCSigState fuel_controller_mode_10ms_Trig_;/* '<S7>/fuel_controller_mode_10ms' */
  ZCSigState fuel_controller_10ms_Trig_ZCE;/* '<S7>/fuel_controller_10ms' */
};

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
struct ODE3_IntgData {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
};

#endif

/* External inputs (root inport signals with default storage) */
struct ExtU_AbstractFuelControl_M1_T {
  real_T PedalAngle;                   /* '<Root>/Pedal Angle' */
  real_T EngineSpeed;                  /* '<Root>/Engine Speed' */
};

/* External outputs (root outports fed by signals with default storage) */
struct ExtY_AbstractFuelControl_M1_T {
  real_T AFref;                        /* '<Root>/AFref' */
  real_T AF;                           /* '<Root>/AF' */
  real_T controller_mode;              /* '<Root>/controller_mode' */
};

/* Parameters (default storage) */
struct P_AbstractFuelControl_M1_T_ {
  real_T AF_sensor_tol;                /* Variable: AF_sensor_tol
                                        * Referenced by: '<S1>/A//F sensor tolerance [0.99 1.01]'
                                        */
  real_T MAF_sensor_tol;               /* Variable: MAF_sensor_tol
                                        * Referenced by: '<S1>/MAF sensor tolerance [0.95 1.05]'
                                        */
  real_T fault_time;                   /* Variable: fault_time
                                        * Referenced by: '<S1>/A//F Sensor Fault Injection'
                                        */
  real_T fuel_inj_tol;                 /* Variable: fuel_inj_tol
                                        * Referenced by: '<S1>/fuel injector tolerance [0.95 1.05]'
                                        */
  real_T kappa_tol;                    /* Variable: kappa_tol
                                        * Referenced by: '<S6>/Kappa tolerance [0.9 1.1]'
                                        */
  real_T tau_ww_tol;                   /* Variable: tau_ww_tol
                                        * Referenced by: '<S6>/tau_ww tolerance [0.9 1.1]'
                                        */
  real32_T ki;                         /* Variable: ki
                                        * Referenced by: '<S12>/Gain1'
                                        */
  real32_T kp;                         /* Variable: kp
                                        * Referenced by: '<S12>/Gain'
                                        */
  real32_T pump_tol;                   /* Variable: pump_tol
                                        * Referenced by: '<S4>/Gain2'
                                        */
  real_T Pwon_Time;                    /* Expression: 0.001
                                        * Referenced by: '<S2>/Pwon'
                                        */
  real_T Pwon_Y0;                      /* Expression: 0
                                        * Referenced by: '<S2>/Pwon'
                                        */
  real_T Pwon_YFinal;                  /* Expression: 1
                                        * Referenced by: '<S2>/Pwon'
                                        */
  real_T PulseGenerator_10ms_Amp;      /* Expression: 1
                                        * Referenced by: '<S2>/PulseGenerator_10ms'
                                        */
  real_T PulseGenerator_10ms_Period;
                               /* Computed Parameter: PulseGenerator_10ms_Period
                                * Referenced by: '<S2>/PulseGenerator_10ms'
                                */
  real_T PulseGenerator_10ms_Duty;
                                 /* Computed Parameter: PulseGenerator_10ms_Duty
                                  * Referenced by: '<S2>/PulseGenerator_10ms'
                                  */
  real_T PulseGenerator_10ms_PhaseDelay;/* Expression: 0.01
                                         * Referenced by: '<S2>/PulseGenerator_10ms'
                                         */
  real_T EngineSpeed9001100_UpperSat;  /* Expression: 3100
                                        * Referenced by: '<S1>/Engine Speed [900 1100]'
                                        */
  real_T EngineSpeed9001100_LowerSat;  /* Expression: 900
                                        * Referenced by: '<S1>/Engine Speed [900 1100]'
                                        */
  real_T rpmtorads_Gain;               /* Expression: pi/30
                                        * Referenced by: '<S1>/(rpm) to (rad//s)'
                                        */
  real_T AFSensorFaultInjection_Y0;    /* Expression: 0
                                        * Referenced by: '<S1>/A//F Sensor Fault Injection'
                                        */
  real_T AFSensorFaultInjection_YFinal;/* Expression: 1
                                        * Referenced by: '<S1>/A//F Sensor Fault Injection'
                                        */
  real_T FaultySensorOutput_Value;     /* Expression: -1
                                        * Referenced by: '<S3>/FaultySensorOutput'
                                        */
  real_T Integrator_IC;                /* Expression: 14.7
                                        * Referenced by: '<S19>/Integrator'
                                        */
  real_T AF_sensor_Gain;               /* Expression: 1
                                        * Referenced by: '<S17>/A//F_sensor'
                                        */
  real_T Switch_Threshold;             /* Expression: 0.5
                                        * Referenced by: '<S3>/Switch'
                                        */
  real_T Throttledelay_A;              /* Computed Parameter: Throttledelay_A
                                        * Referenced by: '<S1>/Throttle delay'
                                        */
  real_T Throttledelay_C;              /* Computed Parameter: Throttledelay_C
                                        * Referenced by: '<S1>/Throttle delay'
                                        */
  real_T Baseopeningangle_Value;       /* Expression: 8.8
                                        * Referenced by: '<S1>/Base opening angle'
                                        */
  real_T theta090_UpperSat;            /* Expression: 90
                                        * Referenced by: '<S1>/theta [0 90]'
                                        */
  real_T theta090_LowerSat;            /* Expression: 0
                                        * Referenced by: '<S1>/theta [0 90]'
                                        */
  real_T p00543bar_IC;                 /* Expression: 0.982
                                        * Referenced by: '<S4>/p0 = 0.543 (bar)'
                                        */
  real_T AtmosphericPressurebar_Value; /* Expression: 1
                                        * Referenced by: '<S1>/Atmospheric Pressure (bar)'
                                        */
  real_T SonicFlow_Value;              /* Expression: 1.0
                                        * Referenced by: '<S5>/Sonic Flow '
                                        */
  real_T Switch_Threshold_g;           /* Expression: 0.5
                                        * Referenced by: '<S5>/Switch'
                                        */
  real_T Integrator_IC_l;              /* Expression: 14.7
                                        * Referenced by: '<S18>/Integrator'
                                        */
  real_T Gain_Gain;                    /* Expression: 50
                                        * Referenced by: '<S19>/Gain'
                                        */
  real_T radstorpm_Gain;               /* Expression: 60/(2*pi)
                                        * Referenced by: '<S6>/rad//s to rpm'
                                        */
  real_T Gain_Gain_l;                  /* Expression: 4*pi/4
                                        * Referenced by: '<S3>/Gain'
                                        */
  real_T uKappa_tableData[20];
  /* Expression: reshape([0.8,.7,.7,.8,.9,.7,.66,.65,.73,.85,.66,.66,.63,.66,.8,.6,.6,.6,.6,.7],5,4)
   * Referenced by: '<S6>/1-Kappa'
   */
  real_T uKappa_bp01Data[5];           /* Expression: [1000,1500,2000,2500,3000]
                                        * Referenced by: '<S6>/1-Kappa'
                                        */
  real_T uKappa_bp02Data[4];           /* Expression: [.1,.2,.3,.4]
                                        * Referenced by: '<S6>/1-Kappa'
                                        */
  real_T Integrator_IC_m;              /* Expression: .0112
                                        * Referenced by: '<S6>/Integrator'
                                        */
  real_T tau_ww_tableData[20];
  /* Expression: reshape([.4,.3,.35,.3,.2,.22,.22,.4,.35,.5,.20,.22,.5,.4,.35,.35,.3,.45,.5,.4],5,4)
   * Referenced by: '<S6>/tau_ww'
   */
  real_T tau_ww_bp01Data[5];           /* Expression: [1000,1500,2000,2500,3000]
                                        * Referenced by: '<S6>/tau_ww'
                                        */
  real_T tau_ww_bp02Data[4];           /* Expression: [0.1,0.2,0.3,0.4]
                                        * Referenced by: '<S6>/tau_ww'
                                        */
  real_T fuelsystemtransportdelay_MaxDel;/* Expression: 10
                                          * Referenced by: '<S3>/fuel system transport delay'
                                          */
  real_T fuelsystemtransportdelay_InitOu;/* Expression: 14.7
                                          * Referenced by: '<S3>/fuel system transport delay'
                                          */
  real_T Gain1_Gain;                   /* Expression: 10
                                        * Referenced by: '<S18>/Gain1'
                                        */
  real_T radstorpm_Gain_p;             /* Expression: 60/(2*pi)
                                        * Referenced by: '<S3>/rad//s to rpm'
                                        */
  real_T delays_tableData[20];
  /* Expression: reshape([0.8,0.6,0.4,0.3,0.2,0.4,0.3,0.2,0.2,0.2,0.3,0.25,0.2,0.2,0.2,0.25,0.2,0.2,0.2,0.2],5,4)
   * Referenced by: '<S3>/delay (s)'
   */
  real_T delays_bp01Data[5];           /* Expression: [800,1000,1500,2000,3000]
                                        * Referenced by: '<S3>/delay (s)'
                                        */
  real_T delays_bp02Data[4];           /* Expression: [0.05,0.15,0.2,0.25]
                                        * Referenced by: '<S3>/delay (s)'
                                        */
  real_T RTVm_Gain;                    /* Expression: 0.41328
                                        * Referenced by: '<S4>/RT//Vm'
                                        */
  real_T Gain_Gain_m;                  /* Expression: -1
                                        * Referenced by: '<S6>/Gain'
                                        */
  real_T Constant_Value;               /* Expression: 1
                                        * Referenced by: '<S6>/Constant'
                                        */
  real32_T Constant3_Value;            /* Computed Parameter: Constant3_Value
                                        * Referenced by: '<S8>/Constant3'
                                        */
  real32_T Constant2_Value;            /* Computed Parameter: Constant2_Value
                                        * Referenced by: '<S8>/Constant2'
                                        */
  real32_T fb_fuel_saturation_UpperSat;
                              /* Computed Parameter: fb_fuel_saturation_UpperSat
                               * Referenced by: '<S8>/fb_fuel_saturation'
                               */
  real32_T fb_fuel_saturation_LowerSat;
                              /* Computed Parameter: fb_fuel_saturation_LowerSat
                               * Referenced by: '<S8>/fb_fuel_saturation'
                               */
  real32_T Constant1_Value;            /* Computed Parameter: Constant1_Value
                                        * Referenced by: '<S11>/Constant1'
                                        */
  real32_T Constant2_Value_d;          /* Computed Parameter: Constant2_Value_d
                                        * Referenced by: '<S11>/Constant2'
                                        */
  real32_T Constant3_Value_c;          /* Computed Parameter: Constant3_Value_c
                                        * Referenced by: '<S11>/Constant3'
                                        */
  real32_T Constant4_Value;            /* Computed Parameter: Constant4_Value
                                        * Referenced by: '<S11>/Constant4'
                                        */
  real32_T Constant5_Value;            /* Computed Parameter: Constant5_Value
                                        * Referenced by: '<S11>/Constant5'
                                        */
  real32_T UnitDelay1_InitialCondition;
                              /* Computed Parameter: UnitDelay1_InitialCondition
                               * Referenced by: '<S11>/UnitDelay1'
                               */
  real32_T Gain_Gain_j;                /* Computed Parameter: Gain_Gain_j
                                        * Referenced by: '<S11>/Gain'
                                        */
  real32_T Constant1_Value_h;          /* Computed Parameter: Constant1_Value_h
                                        * Referenced by: '<S12>/Constant1'
                                        */
  real32_T UnitDelay1_InitialCondition_l;
                            /* Computed Parameter: UnitDelay1_InitialCondition_l
                             * Referenced by: '<S12>/UnitDelay1'
                             */
  real32_T fuel_saturation_UpperSat;
                                 /* Computed Parameter: fuel_saturation_UpperSat
                                  * Referenced by: '<S8>/fuel_saturation'
                                  */
  real32_T fuel_saturation_LowerSat;
                                 /* Computed Parameter: fuel_saturation_LowerSat
                                  * Referenced by: '<S8>/fuel_saturation'
                                  */
  real32_T airbyfuel_reference_power_Value;
                          /* Computed Parameter: airbyfuel_reference_power_Value
                           * Referenced by: '<S9>/airbyfuel_reference_power'
                           */
  real32_T airbyfuel_reference_Value;
                                /* Computed Parameter: airbyfuel_reference_Value
                                 * Referenced by: '<S9>/airbyfuel_reference'
                                 */
  real32_T UnitDelay2_InitialCondition;
                              /* Computed Parameter: UnitDelay2_InitialCondition
                               * Referenced by: '<S14>/Unit Delay2'
                               */
  real32_T sampling_sec_Value;         /* Computed Parameter: sampling_sec_Value
                                        * Referenced by: '<S14>/sampling_sec'
                                        */
  real32_T normal_mode_start_sec_Value;
                              /* Computed Parameter: normal_mode_start_sec_Value
                               * Referenced by: '<S14>/normal_mode_start_sec'
                               */
  real32_T Constant1_Value_f;          /* Computed Parameter: Constant1_Value_f
                                        * Referenced by: '<S15>/Constant1'
                                        */
  real32_T Constant_Value_d;           /* Computed Parameter: Constant_Value_d
                                        * Referenced by: '<S15>/Constant'
                                        */
  real32_T threshold_Value;            /* Computed Parameter: threshold_Value
                                        * Referenced by: '<S16>/threshold'
                                        */
  real32_T Constant1_Value_l;          /* Computed Parameter: Constant1_Value_l
                                        * Referenced by: '<S10>/Constant1'
                                        */
  real32_T Constant2_Value_k;          /* Computed Parameter: Constant2_Value_k
                                        * Referenced by: '<S10>/Constant2'
                                        */
  real32_T Constant3_Value_o;          /* Computed Parameter: Constant3_Value_o
                                        * Referenced by: '<S10>/Constant3'
                                        */
  real32_T DataStoreMemory_InitialValue;
                             /* Computed Parameter: DataStoreMemory_InitialValue
                              * Referenced by: '<S7>/DataStoreMemory'
                              */
  real32_T DataStoreMemory1_InitialValue;
                            /* Computed Parameter: DataStoreMemory1_InitialValue
                             * Referenced by: '<S7>/DataStoreMemory1'
                             */
  real32_T DataStoreMemory2_InitialValue;
                            /* Computed Parameter: DataStoreMemory2_InitialValue
                             * Referenced by: '<S7>/DataStoreMemory2'
                             */
  real32_T DataStoreMemory3_InitialValue;
                            /* Computed Parameter: DataStoreMemory3_InitialValue
                             * Referenced by: '<S7>/DataStoreMemory3'
                             */
  real32_T commanded_fuel_InitialValue;
                              /* Computed Parameter: commanded_fuel_InitialValue
                               * Referenced by: '<S2>/commanded_fuel'
                               */
  real32_T mode_fb1_InitialValue;   /* Computed Parameter: mode_fb1_InitialValue
                                     * Referenced by: '<S2>/mode_fb1'
                                     */
  uint32_T uKappa_maxIndex[2];         /* Computed Parameter: uKappa_maxIndex
                                        * Referenced by: '<S6>/1-Kappa'
                                        */
  uint32_T tau_ww_maxIndex[2];         /* Computed Parameter: tau_ww_maxIndex
                                        * Referenced by: '<S6>/tau_ww'
                                        */
  uint32_T delays_maxIndex[2];         /* Computed Parameter: delays_maxIndex
                                        * Referenced by: '<S3>/delay (s)'
                                        */
  boolean_T UnitDelay1_InitialCondition_c;
                            /* Computed Parameter: UnitDelay1_InitialCondition_c
                             * Referenced by: '<S14>/Unit Delay1'
                             */
  boolean_T UnitDelay1_InitialCondition_f;
                            /* Computed Parameter: UnitDelay1_InitialCondition_f
                             * Referenced by: '<S15>/Unit Delay1'
                             */
  boolean_T UnitDelay_InitialCondition;
                               /* Computed Parameter: UnitDelay_InitialCondition
                                * Referenced by: '<S16>/Unit Delay'
                                */
  boolean_T mode_fb_InitialValue;    /* Computed Parameter: mode_fb_InitialValue
                                      * Referenced by: '<S2>/mode_fb'
                                      */
};

/* Real-time Model Data Structure */
struct tag_RTM_AbstractFuelControl_M_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_AbstractFuelControl_M1_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_AbstractFuelControl_M1_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[6];
  real_T odeF[3][6];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    time_T tStart;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[3];
  } Timing;
};

/* Class declaration for model AbstractFuelControl_M1 */
class AbstractFuelControl_M1 final
{
  /* public data and function members */
 public:
  /* Copy Constructor */
  AbstractFuelControl_M1(AbstractFuelControl_M1 const&) = delete;

  /* Assignment Operator */
  AbstractFuelControl_M1& operator= (AbstractFuelControl_M1 const&) & = delete;

  /* Move Constructor */
  AbstractFuelControl_M1(AbstractFuelControl_M1 &&) = delete;

  /* Move Assignment Operator */
  AbstractFuelControl_M1& operator= (AbstractFuelControl_M1 &&) = delete;

  /* Real-Time Model get method */
  RT_MODEL_AbstractFuelControl__T * getRTM();

  /* Root inports set method */
  void setExternalInputs(const ExtU_AbstractFuelControl_M1_T
    *pExtU_AbstractFuelControl_M1_T)
  {
    AbstractFuelControl_M1_U = *pExtU_AbstractFuelControl_M1_T;
  }

  /* Root outports get method */
  const ExtY_AbstractFuelControl_M1_T &getExternalOutputs() const
  {
    return AbstractFuelControl_M1_Y;
  }

  void ModelPrevZCStateInit();

  /* model start function */
  void start();

  /* Initial conditions function */
  void initialize();

  /* model step function */
  void step();

  /* model terminate function */
  static void terminate();

  /* Constructor */
  AbstractFuelControl_M1();

  /* Destructor */
  ~AbstractFuelControl_M1();

  /* private data and function members */
 private:
  /* External inputs */
  ExtU_AbstractFuelControl_M1_T AbstractFuelControl_M1_U;

  /* External outputs */
  ExtY_AbstractFuelControl_M1_T AbstractFuelControl_M1_Y;

  /* Block signals */
  B_AbstractFuelControl_M1_T AbstractFuelControl_M1_B;

  /* Block states */
  DW_AbstractFuelControl_M1_T AbstractFuelControl_M1_DW;

  /* Tunable parameters */
  static P_AbstractFuelControl_M1_T AbstractFuelControl_M1_P;

  /* Block continuous states */
  X_AbstractFuelControl_M1_T AbstractFuelControl_M1_X;

  /* Block Continuous state disabled vector */
  XDis_AbstractFuelControl_M1_T AbstractFuelControl_M1_XDis;

  /* Triggered events */
  PrevZCX_AbstractFuelControl_M_T AbstractFuelControl_M1_PrevZCX;

  /* Global mass matrix */

  /* Continuous states update member function*/
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  /* Derivatives member function */
  void AbstractFuelControl_M1_derivatives();

  /* Real-Time Model */
  RT_MODEL_AbstractFuelControl__T AbstractFuelControl_M1_M;
};

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'AbstractFuelControl_M1'
 * '<S1>'   : 'AbstractFuelControl_M1/Model 1'
 * '<S2>'   : 'AbstractFuelControl_M1/Model 1/AF_Controller'
 * '<S3>'   : 'AbstractFuelControl_M1/Model 1/Cylinder and Exhaust'
 * '<S4>'   : 'AbstractFuelControl_M1/Model 1/Intake Manifold'
 * '<S5>'   : 'AbstractFuelControl_M1/Model 1/Throttle'
 * '<S6>'   : 'AbstractFuelControl_M1/Model 1/Wall wetting'
 * '<S7>'   : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller'
 * '<S8>'   : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_10ms'
 * '<S9>'   : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_mode_10ms'
 * '<S10>'  : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_pwon'
 * '<S11>'  : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_10ms/air_estimation'
 * '<S12>'  : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_10ms/feedback_PI_controller'
 * '<S13>'  : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_10ms/feedforward_controller'
 * '<S14>'  : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_mode_10ms/normal_mode_detection'
 * '<S15>'  : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_mode_10ms/power_mode_detection'
 * '<S16>'  : 'AbstractFuelControl_M1/Model 1/AF_Controller/fuel_controller/fuel_controller_mode_10ms/sensor_failure_detection'
 * '<S17>'  : 'AbstractFuelControl_M1/Model 1/Cylinder and Exhaust/A//F_sensor'
 * '<S18>'  : 'AbstractFuelControl_M1/Model 1/Cylinder and Exhaust/Filter'
 * '<S19>'  : 'AbstractFuelControl_M1/Model 1/Cylinder and Exhaust/A//F_sensor/Filter'
 */
#endif                                 /* AbstractFuelControl_M1_h_ */
