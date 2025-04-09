//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: automatic_transmission.cpp
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
#include "automatic_transmission.h"
#include "rtwtypes.h"
#include <cmath>
// #include "rt_look.h"
#include "cmath"
#include "limits"

// Named constants for Chart: '<Root>/ShiftLogic'
const int32_T CALL_EVENT{ -1 };

const uint8_T IN_downshifting{ 1U };

const uint8_T IN_first{ 1U };

const uint8_T IN_fourth{ 2U };

const uint8_T IN_second{ 3U };

const uint8_T IN_steady_state{ 2U };

const uint8_T IN_third{ 4U };

const uint8_T IN_upshifting{ 3U };

const int32_T event_DOWN{ 30 };

const int32_T event_UP{ 31 };

extern real_T rt_powd_snf(real_T u0, real_T u1);

// private model entry point functions
extern void automatic_transmission_derivatives();
static void rate_scheduler(automatic_transmission::RT_MODEL *const rtM);
extern "C"
{

#ifndef INTERP
#define INTERP(x,x1,x2,y1,y2)          ( (y1)+(((y2) - (y1))/((x2) - (x1)))*((x)-(x1)) )
#endif

#ifndef ZEROTECHNIQUE
#define ZEROTECHNIQUE

  enum class ZeroTechnique : int_T {
    NORMAL_INTERP,
    AVERAGE_VALUE,
    MIDDLE_VALUE
  };

#endif

  static int_T rt_GetLookupIndex(const real_T *x, int_T xlen, real_T u) ;
}                                      // extern "C"

extern "C"
{
  static real_T rt_Lookup(const real_T *x, int_T xlen, real_T u, const real_T *y);
}                                      // extern "C"

extern "C"
{
  static real_T rt_Lookup2D_Normal (const real_T *xVals, const int_T numX,
    const real_T *yVals, const int_T numY,
    const real_T *zVals,
    const real_T x, const real_T y);
}                                      // extern "C"

extern "C"
{
  real_T rtNaN { -std::numeric_limits<real_T>::quiet_NaN() };

  real_T rtInf { std::numeric_limits<real_T>::infinity() };

  real_T rtMinusInf { -std::numeric_limits<real_T>::infinity() };

  real32_T rtNaNF { -std::numeric_limits<real32_T>::quiet_NaN() };

  real32_T rtInfF { std::numeric_limits<real32_T>::infinity() };

  real32_T rtMinusInfF { -std::numeric_limits<real32_T>::infinity() };
}

extern "C"
{
  // Return rtInf needed by the generated code.
  static real_T rtGetInf(void)
  {
    return rtInf;
  }

  // Get rtInfF needed by the generated code.
  static real32_T rtGetInfF(void)
  {
    return rtInfF;
  }

  // Return rtMinusInf needed by the generated code.
  static real_T rtGetMinusInf(void)
  {
    return rtMinusInf;
  }

  // Return rtMinusInfF needed by the generated code.
  static real32_T rtGetMinusInfF(void)
  {
    return rtMinusInfF;
  }
}

extern "C"
{
  // Return rtNaN needed by the generated code.
  static real_T rtGetNaN(void)
  {
    return rtNaN;
  }

  // Return rtNaNF needed by the generated code.
  static real32_T rtGetNaNF(void)
  {
    return rtNaNF;
  }
}

extern "C"
{
  //
  // Routine to get the index of the input from a table using binary or
  // interpolation search.
  //
  // Inputs:
  // *x   : Pointer to table, x[0] ....x[xlen-1]
  // xlen : Number of values in xtable
  // u    : input value to look up
  //
  // Output:
  // idx  : the index into the table such that:
  // if u is negative
  // x[idx] <= u < x[idx+1]
  // else
  // x[idx] < u <= x[idx+1]
  //
  // Interpolation Search: If the table contains a large number of nearly
  // uniformly spaced entries, i.e., x[n] vs n is linear then the index
  // corresponding to the input can be found in one shot using the linear
  // interpolation formula. Therefore if you have a look-up table block with
  // many data points, using interpolation search might speed up the code.
  // Compile the generated code with the following flag:
  //
  // make_rtw OPTS=-DDOINTERPSEARCH
  //
  // to enable interpolation search.
  //
  static int_T rt_GetLookupIndex(const real_T *x, int_T xlen, real_T u)
  {
    int_T idx{ 0 };

    int_T bottom{ 0 };

    int_T top{ xlen-1 };

    int_T retValue{ 0 };

    boolean_T returnStatus{ 0U };

#ifdef DOINTERPSEARCH

    real_T offset{ 0 };

#endif

    //
    //  Deal with the extreme cases first:
    //    if u <= x[bottom] then return idx = bottom
    //    if u >= x[top]    then return idx = top-1

    if (u <= x[bottom]) {
      retValue = bottom;
      returnStatus = 1U;
    } else if (u >= x[top]) {
      retValue = top-1;
      returnStatus = 1U;
    } else {
      // else required to ensure safe programming, even *
      //  if it's expected that it will never be reached
    }

    if (returnStatus == 0U) {
      if (u < 0) {
        // For negative input find index such that: x[idx] <= u < x[idx+1]
        for (;;) {

#ifdef DOINTERPSEARCH

          offset = (u-x[bottom])/(x[top]-x[bottom]);
          idx = bottom + (int_T)((top-bottom)*(offset-DBL_EPSILON));

#else

          idx = (bottom + top)/2;

#endif

          if (u < x[idx]) {
            top = idx - 1;
          } else if (u >= x[idx+1]) {
            bottom = idx + 1;
          } else {
            // we have x[idx] <= u < x[idx+1], return idx
            retValue = idx;
            break;
          }
        }
      } else {
        // For non-negative input find index such that: x[idx] < u <= x[idx+1]
        for (;;) {

#ifdef DOINTERPSEARCH

          offset = (u-x[bottom])/(x[top]-x[bottom]);
          idx = bottom + (int_T)((top-bottom)*(offset-DBL_EPSILON));

#else

          idx = (bottom + top)/2;

#endif

          if (u <= x[idx]) {
            top = idx - 1;
          } else if (u > x[idx+1]) {
            bottom = idx + 1;
          } else {
            // we have x[idx] < u <= x[idx+1], return idx
            retValue = idx;
            break;
          }
        }
      }
    }

    return retValue;
  }
}                                      // extern "C"

extern "C"
{
  // 1D lookup routine for data type of real_T.
  static real_T rt_Lookup(const real_T *x, int_T xlen, real_T u, const real_T *y)
  {
    int_T idx{ rt_GetLookupIndex(x, xlen, u) };

    real_T num{ y[idx+1] - y[idx] };

    real_T den{ x[idx+1] - x[idx] };

    // Due to the way the binary search is implemented
    // in rt_look.c (rt_GetLookupIndex), den cannot be
    // 0.  Equivalently, m cannot be inf or nan.
    real_T m{ num/den };

    return (y[idx] + (m * (u - x[idx])));
  }
}                                      // extern "C"

extern "C"
{
  // 2D normal lookup routine for data type of real_T.
  static real_T rt_Lookup2D_Normal(const real_T *xVals, const int_T numX,
    const real_T *yVals, const int_T numY,
    const real_T *zVals,
    const real_T x, const real_T y)
  {
    int_T xIdx, yIdx;
    real_T ylo, yhi;
    real_T Zx0yhi, Zx0ylo, xlo, xhi;
    real_T corner1, corner2;
    xIdx = rt_GetLookupIndex(xVals,numX,x);
    xlo = xVals[xIdx];
    xhi = xVals[xIdx+1];
    yIdx = rt_GetLookupIndex(yVals,numY,y);
    ylo = yVals[yIdx];
    yhi = yVals[yIdx+1];
    corner1 = *(zVals + xIdx + (numX * yIdx));
    corner2 = *(zVals + (xIdx+1) + (numX * yIdx));
    Zx0ylo = INTERP(x, xlo, xhi, corner1, corner2);
    corner1 = *(zVals + xIdx + (numX * (yIdx+1)));
    corner2 = *(zVals + (xIdx+1) + (numX*(yIdx+1)));
    Zx0yhi = INTERP(x, xlo, xhi, corner1, corner2);
    return (INTERP(y,ylo,yhi,Zx0ylo,Zx0yhi));
  }
}                                      // extern "C"

extern "C"
{
  // Test if value is infinite
  static boolean_T rtIsInf(real_T value)
  {
    return std::isinf(value);
  }

  // Test if single-precision value is infinite
  static boolean_T rtIsInfF(real32_T value)
  {
    return std::isinf(value);
  }

  // Test if value is not a number
  static boolean_T rtIsNaN(real_T value)
  {
    return std::isnan(value);
  }

  // Test if single-precision value is not a number
  static boolean_T rtIsNaNF(real32_T value)
  {
    return std::isnan(value);
  }
}

//
//         This function updates active task flag for each subrate.
//         The function is called at model base rate, hence the
//         generated code self-manages all its subrates.
//
static void rate_scheduler(automatic_transmission::RT_MODEL *const rtM)
{
  // Compute which subrates run during the next base time step.  Subrates
  //  are an integer multiple of the base rate counter.  Therefore, the subtask
  //  counter is reset when it reaches its limit (zero means run).

  (rtM->Timing.TaskCounters.TID[2])++;
  if ((rtM->Timing.TaskCounters.TID[2]) > 3) {// Sample time: [0.04s, 0.0s]
    rtM->Timing.TaskCounters.TID[2] = 0;
  }
}

//
// This function updates continuous states using the ODE3 fixed-step
// solver algorithm
//
void automatic_transmission::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
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
  int_T nXc { 2 };

  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  // Save the state values at time t in y, we'll use x as ynew.
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  // Assumes that rtsiSetT and ModelOutputs are up-to-date
  // f0 = f(t,y)
  rtsiSetdX(si, f0);
  automatic_transmission_derivatives();

  // f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*));
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  this->step();
  automatic_transmission_derivatives();

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
  automatic_transmission_derivatives();

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

// Function for Chart: '<Root>/ShiftLogic'
void automatic_transmission::gear_state(const int32_T *sfEvent)
{
  switch (rtDW.is_gear_state) {
   case IN_first:
    rtDW.gear = 1.0;
    if (*sfEvent == event_UP) {
      rtDW.is_gear_state = IN_second;
      rtDW.gear = 2.0;
    }
    break;

   case IN_fourth:
    rtDW.gear = 4.0;
    if (*sfEvent == event_DOWN) {
      rtDW.is_gear_state = IN_third;
      rtDW.gear = 3.0;
    }
    break;

   case IN_second:
    rtDW.gear = 2.0;
    switch (*sfEvent) {
     case event_UP:
      rtDW.is_gear_state = IN_third;
      rtDW.gear = 3.0;
      break;

     case event_DOWN:
      rtDW.is_gear_state = IN_first;
      rtDW.gear = 1.0;
      break;
    }
    break;

   case IN_third:
    rtDW.gear = 3.0;
    switch (*sfEvent) {
     case event_UP:
      rtDW.is_gear_state = IN_fourth;
      rtDW.gear = 4.0;
      break;

     case event_DOWN:
      rtDW.is_gear_state = IN_second;
      rtDW.gear = 2.0;
      break;
    }
    break;
  }
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
    if (std::isinf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
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
      y = (rtNaN);
    } else {
      y = std::pow(u0, u1);
    }
  }

  return y;
}

// Model step function
void automatic_transmission::step()
{
  real_T interp_down;
  real_T interp_up;
  int32_T sfEvent;
  if ((&rtM)->isMajorTimeStep()) {
    // set solver stop time
    rtsiSetSolverStopTime(&(&rtM)->solverInfo,(((&rtM)->Timing.clockTick0+1)*
      (&rtM)->Timing.stepSize0));
  }                                    // end MajorTimeStep

  // Update absolute time of base rate at minor time step
  if ((&rtM)->isMinorTimeStep()) {
    (&rtM)->Timing.t[0] = rtsiGetT(&(&rtM)->solverInfo);
  }

  // Gain: '<S5>/mph' incorporates:
  //   Gain: '<S5>/LinearSpeed'
  //   Integrator: '<S5>/Wheel Speed'

  rtDW.VehicleSpeed = 6.2831853071795862 * rtX.WheelSpeed_CSTATE *
    0.011363636363636364;

  // Outport: '<Root>/speed'
  rtY.speed = rtDW.VehicleSpeed;

  // Integrator: '<S1>/Integrator'
  // Limited  Integrator
  if (rtX.Integrator_CSTATE >= 6000.0) {
    rtX.Integrator_CSTATE = 6000.0;
  } else if (rtX.Integrator_CSTATE <= 600.0) {
    rtX.Integrator_CSTATE = 600.0;
  }

  // Outport: '<Root>/RPM' incorporates:
  //   Integrator: '<S1>/Integrator'

  rtY.RPM = rtX.Integrator_CSTATE;
  if ((&rtM)->isMajorTimeStep() &&
      (&rtM)->Timing.TaskCounters.TID[2] == 0) {
    // Chart: '<Root>/ShiftLogic'
    sfEvent = CALL_EVENT;
    if (rtDW.temporalCounter_i1 < MAX_uint32_T) {
      rtDW.temporalCounter_i1++;
    }

    if (rtDW.is_active_c1_automatic_transmis == 0) {
      rtDW.is_active_c1_automatic_transmis = 1U;
      rtDW.is_active_gear_state = 1U;
      rtDW.is_gear_state = IN_first;
      rtDW.gear = 1.0;
      rtDW.is_active_selection_state = 1U;
      rtDW.is_selection_state = IN_steady_state;
    } else {
      if (rtDW.is_active_gear_state != 0) {
        gear_state(&sfEvent);
      }

      if (rtDW.is_active_selection_state != 0) {
        // Outputs for Function Call SubSystem: '<Root>/ThresholdCalculation'
        // Lookup2D: '<S3>/interp_down' incorporates:
        //   Inport: '<Root>/throttle'

        interp_down = rt_Lookup2D_Normal(&rtConstP.interp_down_RowIdx[0], 6,
          &rtConstP.pooled1[0], 4, &rtConstP.interp_down_Table[0], rtU.throttle,
          rtDW.gear);

        // Lookup2D: '<S3>/interp_up' incorporates:
        //   Inport: '<Root>/throttle'

        interp_up = rt_Lookup2D_Normal(&rtConstP.interp_up_RowIdx[0], 6,
          &rtConstP.pooled1[0], 4, &rtConstP.interp_up_Table[0], rtU.throttle,
          rtDW.gear);

        // End of Outputs for SubSystem: '<Root>/ThresholdCalculation'
        switch (rtDW.is_selection_state) {
         case IN_downshifting:
          if ((rtDW.temporalCounter_i1 >= 2U) && (rtDW.VehicleSpeed <=
               interp_down)) {
            sfEvent = event_DOWN;
            if (rtDW.is_active_gear_state != 0) {
              gear_state(&sfEvent);
            }

            rtDW.is_selection_state = IN_steady_state;
          } else if (rtDW.VehicleSpeed > interp_down) {
            rtDW.is_selection_state = IN_steady_state;
          }
          break;

         case IN_steady_state:
          if (rtDW.VehicleSpeed > interp_up) {
            rtDW.temporalCounter_i1 = 0U;
            rtDW.is_selection_state = IN_upshifting;
          } else if (rtDW.VehicleSpeed < interp_down) {
            rtDW.temporalCounter_i1 = 0U;
            rtDW.is_selection_state = IN_downshifting;
          }
          break;

         case IN_upshifting:
          if ((rtDW.temporalCounter_i1 >= 2U) && (rtDW.VehicleSpeed >= interp_up))
          {
            sfEvent = event_UP;
            if (rtDW.is_active_gear_state != 0) {
              gear_state(&sfEvent);
            }

            rtDW.is_selection_state = IN_steady_state;
          } else if (rtDW.VehicleSpeed < interp_up) {
            rtDW.is_selection_state = IN_steady_state;
          }
          break;
        }
      }
    }

    // End of Chart: '<Root>/ShiftLogic'

    // Outport: '<Root>/gear'
    rtY.gear = rtDW.gear;

    // Lookup: '<S7>/Look-Up Table'
    rtDW.LookUpTable = rt_Lookup(&rtConstP.pooled1[0], 4, rtDW.gear,
      &rtConstP.LookUpTable_YData[0]);
  }

  // Gain: '<S5>/FinalDriveRatio2' incorporates:
  //   Integrator: '<S5>/Wheel Speed'

  rtDW.TransmissionRPM = 3.23 * rtX.WheelSpeed_CSTATE;

  // Product: '<S6>/SpeedRatio' incorporates:
  //   Integrator: '<S1>/Integrator'
  //   Product: '<S7>/Product1'

  interp_down = rtDW.LookUpTable * rtDW.TransmissionRPM / rtX.Integrator_CSTATE;

  // Fcn: '<S6>/Impeller' incorporates:
  //   Integrator: '<S1>/Integrator'
  //   Lookup: '<S6>/FactorK'
  //   Product: '<S6>/Quotient'

  interp_up = rt_powd_snf(rtX.Integrator_CSTATE / rt_Lookup(&rtConstP.pooled3[0],
    21, interp_down, &rtConstP.FactorK_YData[0]), 2.0);

  // Gain: '<S1>/engine + impeller inertia' incorporates:
  //   Fcn: '<S6>/Impeller'
  //   Inport: '<Root>/throttle'
  //   Integrator: '<S1>/Integrator'
  //   Lookup2D: '<S1>/EngineTorque'
  //   Sum: '<S1>/Sum'

  rtDW.engineimpellerinertia = (rt_Lookup2D_Normal
    (&rtConstP.EngineTorque_RowIdx[0], 10, &rtConstP.EngineTorque_ColIdx[0], 11,
     &rtConstP.EngineTorque_Table[0], rtU.throttle, rtX.Integrator_CSTATE) -
    interp_up) * 45.472138452209627;

  // Product: '<S7>/Product' incorporates:
  //   Fcn: '<S6>/Impeller'
  //   Lookup: '<S6>/TorqueRatio'
  //   Product: '<S6>/Turbine'

  rtDW.OutputTorque = interp_up * rt_Lookup(&rtConstP.pooled3[0], 21,
    interp_down, &rtConstP.TorqueRatio_YData[0]) * rtDW.LookUpTable;

  // Signum: '<S5>/Sign'
  if (std::isnan(rtDW.VehicleSpeed)) {
    interp_down = (rtNaN);
  } else if (rtDW.VehicleSpeed < 0.0) {
    interp_down = -1.0;
  } else {
    interp_down = (rtDW.VehicleSpeed > 0.0);
  }

  // Gain: '<S5>/Vehicle Inertia' incorporates:
  //   Fcn: '<S5>/RoadLoad'
  //   Gain: '<S5>/Final Drive Ratio1'
  //   Inport: '<Root>/brake'
  //   Product: '<S5>/SignedLoad'
  //   Signum: '<S5>/Sign'
  //   Sum: '<S5>/Sum'
  //   Sum: '<S5>/Sum1'

  rtDW.VehicleInertia = (3.23 * rtDW.OutputTorque - ((0.02 * rt_powd_snf
    (rtDW.VehicleSpeed, 2.0) + 40.0) + rtU.brake) * interp_down) *
    0.082684618362373577;
  if ((&rtM)->isMajorTimeStep()) {
    rt_ertODEUpdateContinuousStates(&(&rtM)->solverInfo);

    // Update absolute time for base rate
    // The "clockTick0" counts the number of times the code of this task has
    //  been executed. The absolute time is the multiplication of "clockTick0"
    //  and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
    //  overflow during the application lifespan selected.

    ++(&rtM)->Timing.clockTick0;
    (&rtM)->Timing.t[0] = rtsiGetSolverStopTime(&(&rtM)->solverInfo);

    {
      // Update absolute timer for sample time: [0.01s, 0.0s]
      // The "clockTick1" counts the number of times the code of this task has
      //  been executed. The resolution of this integer timer is 0.01, which is the step size
      //  of the task. Size of "clockTick1" ensures timer will not overflow during the
      //  application lifespan selected.

      (&rtM)->Timing.clockTick1++;
    }

    rate_scheduler((&rtM));
  }                                    // end MajorTimeStep
}

// Derivatives for root system: '<Root>'
void automatic_transmission::automatic_transmission_derivatives()
{
  automatic_transmission::XDot *_rtXdot;
  boolean_T lsat;
  boolean_T usat;
  _rtXdot = ((XDot *) (&rtM)->derivs);

  // Derivatives for Integrator: '<S5>/Wheel Speed'
  _rtXdot->WheelSpeed_CSTATE = rtDW.VehicleInertia;

  // Derivatives for Integrator: '<S1>/Integrator'
  lsat = (rtX.Integrator_CSTATE <= 600.0);
  usat = (rtX.Integrator_CSTATE >= 6000.0);
  if (((!lsat) && (!usat)) || (lsat && (rtDW.engineimpellerinertia > 0.0)) ||
      (usat && (rtDW.engineimpellerinertia < 0.0))) {
    _rtXdot->Integrator_CSTATE = rtDW.engineimpellerinertia;
  } else {
    // in saturation
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  // End of Derivatives for Integrator: '<S1>/Integrator'
}

// Model initialize function
void automatic_transmission::initialize()
{
  // Registration code
  {
    // Setup solver object
    rtsiSetSimTimeStepPtr(&(&rtM)->solverInfo, &(&rtM)->Timing.simTimeStep);
    rtsiSetTPtr(&(&rtM)->solverInfo, (&rtM)->getTPtrPtr());
    rtsiSetStepSizePtr(&(&rtM)->solverInfo, &(&rtM)->Timing.stepSize0);
    rtsiSetdXPtr(&(&rtM)->solverInfo, &(&rtM)->derivs);
    rtsiSetContStatesPtr(&(&rtM)->solverInfo, (real_T **) &(&rtM)->contStates);
    rtsiSetNumContStatesPtr(&(&rtM)->solverInfo, &(&rtM)->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&rtM)->solverInfo, &(&rtM)
      ->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&rtM)->solverInfo, &(&rtM)
      ->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&rtM)->solverInfo, &(&rtM)
      ->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&(&rtM)->solverInfo, (boolean_T**) &(&rtM)
      ->contStateDisabled);
    rtsiSetErrorStatusPtr(&(&rtM)->solverInfo, (&rtM)->getErrorStatusPtr());
    rtsiSetRTModelPtr(&(&rtM)->solverInfo, (&rtM));
  }

  rtsiSetSimTimeStep(&(&rtM)->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&(&rtM)->solverInfo, false);
  rtsiSetIsContModeFrozen(&(&rtM)->solverInfo, false);
  (&rtM)->intgData.y = (&rtM)->odeY;
  (&rtM)->intgData.f[0] = (&rtM)->odeF[0];
  (&rtM)->intgData.f[1] = (&rtM)->odeF[1];
  (&rtM)->intgData.f[2] = (&rtM)->odeF[2];
  (&rtM)->contStates = ((X *) &rtX);
  (&rtM)->contStateDisabled = ((XDis *) &rtXDis);
  (&rtM)->Timing.tStart = (0.0);
  rtsiSetSolverData(&(&rtM)->solverInfo, static_cast<void *>(&(&rtM)->intgData));
  rtsiSetSolverName(&(&rtM)->solverInfo,"ode3");
  (&rtM)->setTPtr(&(&rtM)->Timing.tArray[0]);
  (&rtM)->Timing.stepSize0 = 0.01;

  // InitializeConditions for Integrator: '<S5>/Wheel Speed'
  rtX.WheelSpeed_CSTATE = 0.0;

  // InitializeConditions for Integrator: '<S1>/Integrator'
  rtX.Integrator_CSTATE = 1000.0;
}

time_T** automatic_transmission::RT_MODEL::getTPtrPtr()
{
  return &(Timing.t);
}

time_T* automatic_transmission::RT_MODEL::getTPtr() const
{
  return (Timing.t);
}

void automatic_transmission::RT_MODEL::setTPtr(time_T* aTPtr)
{
  (Timing.t = aTPtr);
}

boolean_T automatic_transmission::RT_MODEL::isMinorTimeStep() const
{
  return ((Timing.simTimeStep) == MINOR_TIME_STEP);
}

boolean_T automatic_transmission::RT_MODEL::getStopRequested() const
{
  return (Timing.stopRequestedFlag);
}

void automatic_transmission::RT_MODEL::setStopRequested(boolean_T aStopRequested)
{
  (Timing.stopRequestedFlag = aStopRequested);
}

boolean_T automatic_transmission::RT_MODEL::isMajorTimeStep() const
{
  return ((Timing.simTimeStep) == MAJOR_TIME_STEP);
}

boolean_T* automatic_transmission::RT_MODEL::getStopRequestedPtr()
{
  return (&(Timing.stopRequestedFlag));
}

const char_T** automatic_transmission::RT_MODEL::getErrorStatusPtr()
{
  return &errorStatus;
}

time_T automatic_transmission::RT_MODEL::getTStart() const
{
  return (Timing.tStart);
}

const char_T* automatic_transmission::RT_MODEL::getErrorStatus() const
{
  return (errorStatus);
}

void automatic_transmission::RT_MODEL::setErrorStatus(const char_T* const
  aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

// Constructor
automatic_transmission::automatic_transmission() :
  rtU(),
  rtY(),
  rtDW(),
  rtX(),
  rtXDis(),
  rtM()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
automatic_transmission::~automatic_transmission() = default;

// Real-Time Model get method
automatic_transmission::RT_MODEL * automatic_transmission::getRTM()
{
  return (&rtM);
}

//
// File trailer for generated code.
//
// [EOF]
//
