//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: leadCar.cpp
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
#include "leadCar.h"
#include <cmath>
#include "rtwtypes.h"

//
// This function updates continuous states using the ODE3 fixed-step
// solver algorithm
//
void leadCar::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
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
  leadCar_derivatives();

  // f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*));
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  this->step();
  leadCar_derivatives();

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
  leadCar_derivatives();

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

// Model step function
void leadCar::step()
{
  if ((&leadCar_M)->isMajorTimeStep()) {
    // set solver stop time
    rtsiSetSolverStopTime(&(&leadCar_M)->solverInfo,(((&leadCar_M)
      ->Timing.clockTick0+1)*(&leadCar_M)->Timing.stepSize0));
  }                                    // end MajorTimeStep

  // Update absolute time of base rate at minor time step
  if ((&leadCar_M)->isMinorTimeStep()) {
    (&leadCar_M)->Timing.t[0] = rtsiGetT(&(&leadCar_M)->solverInfo);
  }

  {
    real_T lastSin_tmp;
    if ((&leadCar_M)->isMajorTimeStep()) {
      // Sin: '<Root>/Sine Wave'
      if (leadCar_DW.systemEnable != 0) {
        lastSin_tmp = (((&leadCar_M)->Timing.clockTick1) * 0.1);
        leadCar_DW.lastSin = std::sin(0.2 * lastSin_tmp);
        leadCar_DW.lastCos = std::cos(0.2 * lastSin_tmp);
        leadCar_DW.systemEnable = 0;
      }

      // Outport: '<Root>/a_lead' incorporates:
      //   Sin: '<Root>/Sine Wave'

      leadCar_Y.a_lead = ((leadCar_DW.lastSin * 0.99980000666657776 +
                           leadCar_DW.lastCos * -0.019998666693333084) *
                          0.99980000666657776 + (leadCar_DW.lastCos *
        0.99980000666657776 - leadCar_DW.lastSin * -0.019998666693333084) *
                          0.019998666693333084) * 0.6;
    }

    // Outport: '<Root>/d_lead' incorporates:
    //   Constant: '<Root>/x0 lead'
    //   Integrator: '<S1>/Integrator1'
    //   Sum: '<S1>/Sum'

    leadCar_Y.d_lead = leadCar_X.Integrator1_CSTATE + 50.0;

    // Sum: '<S1>/Sum1' incorporates:
    //   Constant: '<Root>/v0 lead'
    //   TransferFcn: '<S1>/Transfer Fcn'

    leadCar_B.Sum1 = 2.0 * leadCar_X.TransferFcn_CSTATE + 25.0;

    // Outport: '<Root>/v_lead'
    leadCar_Y.v_lead = leadCar_B.Sum1;

    // Integrator: '<S1>/Integrator'
    leadCar_B.Integrator = leadCar_X.Integrator_CSTATE;
  }

  if ((&leadCar_M)->isMajorTimeStep()) {
    real_T HoldSine;
    if ((&leadCar_M)->isMajorTimeStep()) {
      // Update for Sin: '<Root>/Sine Wave'
      HoldSine = leadCar_DW.lastSin;
      leadCar_DW.lastSin = leadCar_DW.lastSin * 0.99980000666657776 +
        leadCar_DW.lastCos * 0.019998666693333084;
      leadCar_DW.lastCos = leadCar_DW.lastCos * 0.99980000666657776 - HoldSine *
        0.019998666693333084;
    }
  }                                    // end MajorTimeStep

  if ((&leadCar_M)->isMajorTimeStep()) {
    rt_ertODEUpdateContinuousStates(&(&leadCar_M)->solverInfo);

    // Update absolute time for base rate
    // The "clockTick0" counts the number of times the code of this task has
    //  been executed. The absolute time is the multiplication of "clockTick0"
    //  and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
    //  overflow during the application lifespan selected.

    ++(&leadCar_M)->Timing.clockTick0;
    (&leadCar_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&leadCar_M)->solverInfo);

    {
      // Update absolute timer for sample time: [0.1s, 0.0s]
      // The "clockTick1" counts the number of times the code of this task has
      //  been executed. The resolution of this integer timer is 0.1, which is the step size
      //  of the task. Size of "clockTick1" ensures timer will not overflow during the
      //  application lifespan selected.

      (&leadCar_M)->Timing.clockTick1++;
    }
  }                                    // end MajorTimeStep
}

// Derivatives for root system: '<Root>'
void leadCar::leadCar_derivatives()
{
  leadCar::XDot_leadCar_T *_rtXdot;
  _rtXdot = ((XDot_leadCar_T *) (&leadCar_M)->derivs);

  // Derivatives for Integrator: '<S1>/Integrator1'
  _rtXdot->Integrator1_CSTATE = leadCar_B.Sum1;

  // Derivatives for TransferFcn: '<S1>/Transfer Fcn'
  _rtXdot->TransferFcn_CSTATE = -2.0 * leadCar_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += leadCar_B.Integrator;

  // Derivatives for Integrator: '<S1>/Integrator' incorporates:
  //   Outport: '<Root>/a_lead'

  _rtXdot->Integrator_CSTATE = leadCar_Y.a_lead;
}

// Model initialize function
void leadCar::initialize()
{
  // Registration code
  {
    // Setup solver object
    rtsiSetSimTimeStepPtr(&(&leadCar_M)->solverInfo, &(&leadCar_M)
                          ->Timing.simTimeStep);
    rtsiSetTPtr(&(&leadCar_M)->solverInfo, (&leadCar_M)->getTPtrPtr());
    rtsiSetStepSizePtr(&(&leadCar_M)->solverInfo, &(&leadCar_M)
                       ->Timing.stepSize0);
    rtsiSetdXPtr(&(&leadCar_M)->solverInfo, &(&leadCar_M)->derivs);
    rtsiSetContStatesPtr(&(&leadCar_M)->solverInfo, (real_T **) &(&leadCar_M)
                         ->contStates);
    rtsiSetNumContStatesPtr(&(&leadCar_M)->solverInfo, &(&leadCar_M)
      ->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&leadCar_M)->solverInfo, &(&leadCar_M)
      ->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&leadCar_M)->solverInfo, &(&leadCar_M
      )->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&leadCar_M)->solverInfo, &(&leadCar_M)
      ->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&(&leadCar_M)->solverInfo, (boolean_T**)
      &(&leadCar_M)->contStateDisabled);
    rtsiSetErrorStatusPtr(&(&leadCar_M)->solverInfo, (&leadCar_M)
                          ->getErrorStatusPtr());
    rtsiSetRTModelPtr(&(&leadCar_M)->solverInfo, (&leadCar_M));
  }

  rtsiSetSimTimeStep(&(&leadCar_M)->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&(&leadCar_M)->solverInfo, false);
  rtsiSetIsContModeFrozen(&(&leadCar_M)->solverInfo, false);
  (&leadCar_M)->intgData.y = (&leadCar_M)->odeY;
  (&leadCar_M)->intgData.f[0] = (&leadCar_M)->odeF[0];
  (&leadCar_M)->intgData.f[1] = (&leadCar_M)->odeF[1];
  (&leadCar_M)->intgData.f[2] = (&leadCar_M)->odeF[2];
  (&leadCar_M)->contStates = ((X_leadCar_T *) &leadCar_X);
  (&leadCar_M)->contStateDisabled = ((XDis_leadCar_T *) &leadCar_XDis);
  (&leadCar_M)->Timing.tStart = (0.0);
  rtsiSetSolverData(&(&leadCar_M)->solverInfo, static_cast<void *>(&(&leadCar_M
    )->intgData));
  rtsiSetSolverName(&(&leadCar_M)->solverInfo,"ode3");
  (&leadCar_M)->setTPtr(&(&leadCar_M)->Timing.tArray[0]);
  (&leadCar_M)->Timing.stepSize0 = 0.1;

  // InitializeConditions for Integrator: '<S1>/Integrator1'
  leadCar_X.Integrator1_CSTATE = 0.0;

  // InitializeConditions for TransferFcn: '<S1>/Transfer Fcn'
  leadCar_X.TransferFcn_CSTATE = 0.0;

  // InitializeConditions for Integrator: '<S1>/Integrator'
  leadCar_X.Integrator_CSTATE = 0.0;

  // Enable for Sin: '<Root>/Sine Wave'
  leadCar_DW.systemEnable = 1;
}

// Model terminate function
void leadCar::terminate()
{
  // (no terminate code required)
}

time_T** leadCar::RT_MODEL_leadCar_T::getTPtrPtr()
{
  return &(Timing.t);
}

time_T* leadCar::RT_MODEL_leadCar_T::getTPtr() const
{
  return (Timing.t);
}

void leadCar::RT_MODEL_leadCar_T::setTPtr(time_T* aTPtr)
{
  (Timing.t = aTPtr);
}

boolean_T leadCar::RT_MODEL_leadCar_T::isMinorTimeStep() const
{
  return ((Timing.simTimeStep) == MINOR_TIME_STEP);
}

boolean_T leadCar::RT_MODEL_leadCar_T::getStopRequested() const
{
  return (Timing.stopRequestedFlag);
}

void leadCar::RT_MODEL_leadCar_T::setStopRequested(boolean_T aStopRequested)
{
  (Timing.stopRequestedFlag = aStopRequested);
}

boolean_T leadCar::RT_MODEL_leadCar_T::isMajorTimeStep() const
{
  return ((Timing.simTimeStep) == MAJOR_TIME_STEP);
}

boolean_T* leadCar::RT_MODEL_leadCar_T::getStopRequestedPtr()
{
  return (&(Timing.stopRequestedFlag));
}

const char_T** leadCar::RT_MODEL_leadCar_T::getErrorStatusPtr()
{
  return &errorStatus;
}

time_T leadCar::RT_MODEL_leadCar_T::getTStart() const
{
  return (Timing.tStart);
}

const char_T* leadCar::RT_MODEL_leadCar_T::getErrorStatus() const
{
  return (errorStatus);
}

void leadCar::RT_MODEL_leadCar_T::setErrorStatus(const char_T* const
  aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

// Constructor
leadCar::leadCar() :
  leadCar_Y(),
  leadCar_B(),
  leadCar_DW(),
  leadCar_X(),
  leadCar_XDis(),
  leadCar_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
leadCar::~leadCar() = default;

// Real-Time Model get method
leadCar::RT_MODEL_leadCar_T * leadCar::getRTM()
{
  return (&leadCar_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
