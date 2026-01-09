/*
 * AbstractFuelControl_M1.cpp
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

#include "AbstractFuelControl_M1.h"
#include <cmath>
#include "rtwtypes.h"
#include "AbstractFuelControl_M1_private.h"
#include "cmath"
#include "zero_crossing_types.h"

static void rate_scheduler(RT_MODEL_AbstractFuelControl__T *const
  AbstractFuelControl_M1_M);

/* For variable transport delay block, find the real delay time */
real_T rt_VTDelayfindtDInterpolate(
  real_T x,real_T* uBuf,int_T bufSz,int_T head,int_T tail,int_T* pLast,real_T t,
  real_T tStart,boolean_T discrete,boolean_T minorStepAndTAtLastMajorOutput,
  real_T initOutput,real_T* appliedDelay)
{
  int_T n, k;
  real_T f;
  int_T kp1;
  real_T tminustD = 0;
  real_T tL = 0;
  real_T tR = 0;
  real_T uL = 0;
  real_T uR = 0;
  real_T uD, fU;
  real_T* tBuf = uBuf + bufSz;
  real_T* xBuf = uBuf + 2* bufSz;
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the entry at head has not been added */
    if (*pLast == head) {
      *pLast = (*pLast == 0) ? bufSz-1 : *pLast-1;
    }

    head = (head == 0) ? bufSz-1 : head-1;
  }

  /*
   * The loop below finds k such that:
   *      x(t)-x(tminustD) =1 or
   *      x - xBuf[k+1] <= 1.0 < x - xBuf[k]
   *
   * Note that we have:
   *
   * tStart = tBuf[0] < tBuf[1] < ... < tBuf[tail] < ... tBuf[head] <= t
   *      0 = xBuf[0] < xBuf[1] < ... < xBuf[tail] < ... xBuf[head] <  x
   *
   * This is true if we assume the direction of transport is always positive
   * such as a flow goes through a pipe from one end to another. However, for
   * model such as convey belt, the transportation can change direction. For
   * this case, there will be more than one solution to x(t)-x(tminustD) = 1,
   * should found the minimum tminustD and tminustD > 0. The search will not
   * be as efficient as the following code.
   */

  /*
   * when x<=1, physically it means the flow didn't reach the output yet,
   * t-tD will be less then zero, so force output to be the initial output
   */
  if (x <= 1) {
    return initOutput;
  }

  /*
   * if the x is monoton increase, only one solution. use k=pLast for now
   */
  k= *pLast;
  n = 0;
  for (;;) {
    n++;
    if (n>bufSz)
      break;
    if (x - xBuf[k] > 1.0) {
      /* move k forward, unless k = head */
      if (k == head) {
        /* xxx this situation means tD= appliedDelay = 0
         *
         * linearly interpolate using (tBuf[head], xBuf[head])
         * and (t,x) to find (tD,xD) such that: x - xD = 1.0
         */
        int_T km1;
        f = (x - 1.0 - xBuf[k]) / (x - xBuf[k]);
        tminustD = (1.0-f)*tBuf[k] + f*t;
        km1 = k-1;
        if (km1 < 0)
          km1 = bufSz-1;
        tL = tBuf[km1];
        tR = tBuf[k];
        uL = uBuf[km1];
        uR = uBuf[k];
        break;
      }

      kp1 = k+1;
      if (kp1 == bufSz)
        kp1 = 0;
      if (x - xBuf[kp1] <= 1.0) {
        /*
         * linearly interpolate using (tBuf[k], xBuf[k])
         * and  (tBuf[k+1], xBuf[k+1]) to find (tminustD,xD)
         * such that: x - xD = 1.0
         */
        f = (x - 1.0 - xBuf[k]) / (xBuf[kp1] - xBuf[k]);
        tL = tBuf[k];
        tR = tBuf[kp1];
        uL = uBuf[k];
        uR = uBuf[kp1];
        tminustD = (1.0-f)*tL + f*tR;
        break;
      }

      k = kp1;
    } else {
      /* moved k backward, unless k = tail */
      if (k == tail) {
        /* This situation means tminustD <= Ttail*/
        f = (x - 1.0)/xBuf[k];
        if (discrete) {
          return(uBuf[tail]);
        }

        kp1 = k+1;
        if (kp1 == bufSz)
          kp1 = 0;

        /* * linearly interpolate using (tStart, 0)
         * and  (tBuf[tail], xBuf[tail]) to find (tminustD,xD)
         * such that: x - xD = 1.0
         */

        /* Here it is better to use Tstart because since x>1, tminustD
         * must > 0. Since x is monotone increase, its linearity is
         * better.
         */
        tminustD = (1-f)*tStart + f*tBuf[k];

        /* linearly interpolate using (t[tail], x[tail])
         * and  (tBuf[tail+1], xBuf[tail+1]) to find (tminustD,xD)
         * such that: x - xD = 1.0.
         * For time delay block, use t[tail] and t[tail+1], not good
         * for transport delay block since it may give tminstD < 0
         */

        /*  f = (tBuf[kp1]-tBuf[k])/(xBuf[kp1]-xBuf[k]);
         *  tminustD = tBuf[kp1]-f*(1+xBuf[kp1]-x);
         */
        tL = tBuf[k];
        tR = tBuf[kp1];
        uL = uBuf[k];
        uR = uBuf[kp1];
        break;
      }

      k = k - 1;
      if (k < 0)
        k = bufSz-1;
    }
  }

  *pLast = k;
  if (tR == tL) {
    fU = 1.0;
  } else {
    fU = (tminustD-tL)/(tR-tL);
  }

  /* for discrete signal, no interpolation, use either uL or uR
   * depend on wehre tminustD is.
   */
  if (discrete) {
    uD= (fU > (1.0-fU))? uR: uL;
  } else {
    uD = (1.0-fU)*uL + fU*uR;
  }

  /* we want return tD= t-(t-tD);*/
  *appliedDelay = t-tminustD;
  return uD;
}

real_T look2_binlxpw(real_T u0, real_T u1, const real_T bp0[], const real_T bp1[],
                     const real_T table[], const uint32_T maxIndex[], uint32_T
                     stride)
{
  real_T fractions[2];
  real_T frac;
  real_T yL_0d0;
  real_T yL_0d1;
  uint32_T bpIndices[2];
  uint32_T bpIdx;
  uint32_T iLeft;
  uint32_T iRght;

  /* Column-major Lookup 2-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex[0U]]) {
    /* Binary Search */
    bpIdx = maxIndex[0U] >> 1U;
    iLeft = 0U;
    iRght = maxIndex[0U];
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex[0U] - 1U;
    frac = (u0 - bp0[maxIndex[0U] - 1U]) / (bp0[maxIndex[0U]] - bp0[maxIndex[0U]
      - 1U]);
  }

  fractions[0U] = frac;
  bpIndices[0U] = iLeft;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u1 <= bp1[0U]) {
    iLeft = 0U;
    frac = (u1 - bp1[0U]) / (bp1[1U] - bp1[0U]);
  } else if (u1 < bp1[maxIndex[1U]]) {
    /* Binary Search */
    bpIdx = maxIndex[1U] >> 1U;
    iLeft = 0U;
    iRght = maxIndex[1U];
    while (iRght - iLeft > 1U) {
      if (u1 < bp1[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u1 - bp1[iLeft]) / (bp1[iLeft + 1U] - bp1[iLeft]);
  } else {
    iLeft = maxIndex[1U] - 1U;
    frac = (u1 - bp1[maxIndex[1U] - 1U]) / (bp1[maxIndex[1U]] - bp1[maxIndex[1U]
      - 1U]);
  }

  /* Column-major Interpolation 2-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  bpIdx = iLeft * stride + bpIndices[0U];
  yL_0d0 = table[bpIdx];
  yL_0d0 += (table[bpIdx + 1U] - yL_0d0) * fractions[0U];
  bpIdx += stride;
  yL_0d1 = table[bpIdx];
  return (((table[bpIdx + 1U] - yL_0d1) * fractions[0U] + yL_0d1) - yL_0d0) *
    frac + yL_0d0;
}

/*
 *         This function updates active task flag for each subrate.
 *         The function is called at model base rate, hence the
 *         generated code self-manages all its subrates.
 */
static void rate_scheduler(RT_MODEL_AbstractFuelControl__T *const
  AbstractFuelControl_M1_M)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (AbstractFuelControl_M1_M->Timing.TaskCounters.TID[2])++;
  if ((AbstractFuelControl_M1_M->Timing.TaskCounters.TID[2]) > 9) {/* Sample time: [0.01s, 0.0s] */
    AbstractFuelControl_M1_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
void AbstractFuelControl_M1::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
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
  int_T nXc { 6 };

  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  AbstractFuelControl_M1_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  this->step();
  AbstractFuelControl_M1_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  this->step();
  AbstractFuelControl_M1_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void AbstractFuelControl_M1::step()
{
  if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M))) {
    /* set solver stop time */
    if (!((&AbstractFuelControl_M1_M)->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&(&AbstractFuelControl_M1_M)->solverInfo,
                            (((&AbstractFuelControl_M1_M)->Timing.clockTickH0 +
        1) * (&AbstractFuelControl_M1_M)->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&(&AbstractFuelControl_M1_M)->solverInfo,
                            (((&AbstractFuelControl_M1_M)->Timing.clockTick0 + 1)
        * (&AbstractFuelControl_M1_M)->Timing.stepSize0 +
        (&AbstractFuelControl_M1_M)->Timing.clockTickH0 *
        (&AbstractFuelControl_M1_M)->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep((&AbstractFuelControl_M1_M))) {
    (&AbstractFuelControl_M1_M)->Timing.t[0] = rtsiGetT
      (&(&AbstractFuelControl_M1_M)->solverInfo);
  }

  {
    real_T appliedDelay;
    real_T rtb_AF_sensor;
    real_T rtb_Gain2;
    real_T rtb_Integrator;
    real_T rtb_Kappatolerance0911;
    real_T rtb_c_air_chrggcyl;
    real_T rtb_fuel_puddle_evap;
    real_T rtb_fuelsystemtransportdelay;
    real_T rtb_pratio;
    real_T rtb_radstorpm;
    real32_T rtb_DataStoreRead;
    real32_T rtb_DataStoreRead2;
    real32_T rtb_Sum2;
    real32_T rtb_Sum3;
    real32_T rtb_Sum_o;
    boolean_T rtb_DataStoreRead1_l;
    boolean_T rtb_LogicalOperator;
    boolean_T rtb_LogicalOperator_c;
    boolean_T rtb_RelationalOperator1;
    ZCEventType zcEvent;

    /* Outputs for Atomic SubSystem: '<Root>/Model 1' */
    /* Saturate: '<S1>/Engine Speed [900 1100]' incorporates:
     *  Inport: '<Root>/Engine Speed'
     */
    if (AbstractFuelControl_M1_U.EngineSpeed >
        AbstractFuelControl_M1_P.EngineSpeed9001100_UpperSat) {
      rtb_Integrator = AbstractFuelControl_M1_P.EngineSpeed9001100_UpperSat;
    } else if (AbstractFuelControl_M1_U.EngineSpeed <
               AbstractFuelControl_M1_P.EngineSpeed9001100_LowerSat) {
      rtb_Integrator = AbstractFuelControl_M1_P.EngineSpeed9001100_LowerSat;
    } else {
      rtb_Integrator = AbstractFuelControl_M1_U.EngineSpeed;
    }

    /* Gain: '<S1>/(rpm) to (rad//s)' incorporates:
     *  Saturate: '<S1>/Engine Speed [900 1100]'
     */
    rtb_radstorpm = AbstractFuelControl_M1_P.rpmtorads_Gain * rtb_Integrator;

    /* Gain: '<S17>/A//F_sensor' incorporates:
     *  Integrator: '<S19>/Integrator'
     */
    rtb_AF_sensor = AbstractFuelControl_M1_P.AF_sensor_Gain *
      AbstractFuelControl_M1_X.Integrator_CSTATE;

    /* Sum: '<S1>/Sum' incorporates:
     *  Constant: '<S1>/Base opening angle'
     *  TransferFcn: '<S1>/Throttle delay'
     */
    rtb_fuelsystemtransportdelay = AbstractFuelControl_M1_P.Throttledelay_C *
      AbstractFuelControl_M1_X.Throttledelay_CSTATE +
      AbstractFuelControl_M1_P.Baseopeningangle_Value;

    /* Saturate: '<S1>/theta [0 90]' */
    if (rtb_fuelsystemtransportdelay >
        AbstractFuelControl_M1_P.theta090_UpperSat) {
      rtb_fuelsystemtransportdelay = AbstractFuelControl_M1_P.theta090_UpperSat;
    } else if (rtb_fuelsystemtransportdelay <
               AbstractFuelControl_M1_P.theta090_LowerSat) {
      rtb_fuelsystemtransportdelay = AbstractFuelControl_M1_P.theta090_LowerSat;
    }

    /* End of Saturate: '<S1>/theta [0 90]' */

    /* MinMax: '<S5>/MinMax' incorporates:
     *  Constant: '<S1>/Atmospheric Pressure (bar)'
     *  Integrator: '<S4>/p0 = 0.543 (bar)'
     *  Product: '<S5>/Product1'
     *  Product: '<S5>/Product2'
     */
    rtb_pratio = std::fmin(AbstractFuelControl_M1_X.p00543bar_CSTATE /
      AbstractFuelControl_M1_P.AtmosphericPressurebar_Value, 1.0 /
      AbstractFuelControl_M1_X.p00543bar_CSTATE *
      AbstractFuelControl_M1_P.AtmosphericPressurebar_Value);

    /* Switch: '<S5>/Switch' incorporates:
     *  Constant: '<S5>/Sonic Flow '
     *  Fcn: '<S5>/g(pratio)'
     */
    if (rtb_pratio >= AbstractFuelControl_M1_P.Switch_Threshold_g) {
      /* Fcn: '<S5>/g(pratio)' */
      rtb_pratio -= rtb_pratio * rtb_pratio;
      if (rtb_pratio < 0.0) {
        rtb_Integrator = -std::sqrt(-rtb_pratio);
      } else {
        rtb_Integrator = std::sqrt(rtb_pratio);
      }

      rtb_pratio = 2.0 * rtb_Integrator;
    } else {
      rtb_pratio = AbstractFuelControl_M1_P.SonicFlow_Value;
    }

    /* End of Switch: '<S5>/Switch' */

    /* Sum: '<S5>/Sum' incorporates:
     *  Constant: '<S1>/Atmospheric Pressure (bar)'
     *  Integrator: '<S4>/p0 = 0.543 (bar)'
     */
    rtb_Integrator = AbstractFuelControl_M1_P.AtmosphericPressurebar_Value -
      AbstractFuelControl_M1_X.p00543bar_CSTATE;

    /* Signum: '<S5>/flow direction' */
    if (std::isnan(rtb_Integrator)) {
      rtb_Integrator = (rtNaN);
    } else if (rtb_Integrator < 0.0) {
      rtb_Integrator = -1.0;
    } else {
      rtb_Integrator = (rtb_Integrator > 0.0);
    }

    /* Product: '<S5>/Product' incorporates:
     *  Fcn: '<S5>/f(theta)'
     *  Signum: '<S5>/flow direction'
     */
    rtb_pratio = (((2.821 - 0.05231 * rtb_fuelsystemtransportdelay) + 0.10299 *
                   rtb_fuelsystemtransportdelay * rtb_fuelsystemtransportdelay)
                  - 0.00063 * rtb_fuelsystemtransportdelay *
                  rtb_fuelsystemtransportdelay * rtb_fuelsystemtransportdelay) *
      rtb_pratio * rtb_Integrator;

    /* Outputs for Atomic SubSystem: '<S1>/AF_Controller' */
    /* Step: '<S2>/Pwon' incorporates:
     *  Step: '<S1>/A//F Sensor Fault Injection'
     */
    rtb_Integrator = (&AbstractFuelControl_M1_M)->Timing.t[0];
    if (rtb_Integrator < AbstractFuelControl_M1_P.Pwon_Time) {
      /* Step: '<S2>/Pwon' */
      AbstractFuelControl_M1_B.Pwon = AbstractFuelControl_M1_P.Pwon_Y0;
    } else {
      /* Step: '<S2>/Pwon' */
      AbstractFuelControl_M1_B.Pwon = AbstractFuelControl_M1_P.Pwon_YFinal;
    }

    /* End of Step: '<S2>/Pwon' */
    if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M)) &&
        (&AbstractFuelControl_M1_M)->Timing.TaskCounters.TID[1] == 0) {
      /* DiscretePulseGenerator: '<S2>/PulseGenerator_10ms' */
      AbstractFuelControl_M1_B.PulseGenerator_10ms =
        (AbstractFuelControl_M1_DW.clockTickCounter <
         AbstractFuelControl_M1_P.PulseGenerator_10ms_Duty) &&
        (AbstractFuelControl_M1_DW.clockTickCounter >= 0) ?
        AbstractFuelControl_M1_P.PulseGenerator_10ms_Amp : 0.0;

      /* DiscretePulseGenerator: '<S2>/PulseGenerator_10ms' */
      if (AbstractFuelControl_M1_DW.clockTickCounter >=
          AbstractFuelControl_M1_P.PulseGenerator_10ms_Period - 1.0) {
        AbstractFuelControl_M1_DW.clockTickCounter = 0;
      } else {
        AbstractFuelControl_M1_DW.clockTickCounter++;
      }
    }

    /* Outputs for Atomic SubSystem: '<S2>/fuel_controller' */
    /* DataStoreWrite: '<S7>/DataStoreWrite' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion1'
     */
    AbstractFuelControl_M1_DW.engine_speed = static_cast<real32_T>(rtb_radstorpm);

    /* DataStoreWrite: '<S7>/DataStoreWrite3' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion4'
     */
    AbstractFuelControl_M1_DW.throttle_angle = static_cast<real32_T>
      (rtb_fuelsystemtransportdelay);

    /* DataStoreWrite: '<S7>/DataStoreWrite1' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion2'
     *  Gain: '<S1>/MAF sensor tolerance [0.95 1.05]'
     */
    AbstractFuelControl_M1_DW.throttle_flow = static_cast<real32_T>
      (AbstractFuelControl_M1_P.MAF_sensor_tol * rtb_pratio);

    /* End of Outputs for SubSystem: '<S2>/fuel_controller' */
    /* End of Outputs for SubSystem: '<S1>/AF_Controller' */

    /* Step: '<S1>/A//F Sensor Fault Injection' */
    if (rtb_Integrator < AbstractFuelControl_M1_P.fault_time) {
      rtb_Integrator = AbstractFuelControl_M1_P.AFSensorFaultInjection_Y0;
    } else {
      rtb_Integrator = AbstractFuelControl_M1_P.AFSensorFaultInjection_YFinal;
    }

    /* Switch: '<S3>/Switch' incorporates:
     *  Constant: '<S3>/FaultySensorOutput'
     *  Step: '<S1>/A//F Sensor Fault Injection'
     */
    if (rtb_Integrator >= AbstractFuelControl_M1_P.Switch_Threshold) {
      rtb_Integrator = AbstractFuelControl_M1_P.FaultySensorOutput_Value;
    } else {
      rtb_Integrator = rtb_AF_sensor;
    }

    /* Outputs for Atomic SubSystem: '<S1>/AF_Controller' */
    /* Outputs for Atomic SubSystem: '<S2>/fuel_controller' */
    /* DataStoreWrite: '<S7>/DataStoreWrite2' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion3'
     *  Gain: '<S1>/A//F sensor tolerance [0.99 1.01]'
     *  Switch: '<S3>/Switch'
     */
    AbstractFuelControl_M1_DW.airbyfuel_meas = static_cast<real32_T>
      (AbstractFuelControl_M1_P.AF_sensor_tol * rtb_Integrator);
    if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M)) &&
        (&AbstractFuelControl_M1_M)->Timing.TaskCounters.TID[1] == 0) {
      /* Outputs for Triggered SubSystem: '<S7>/fuel_controller_10ms' incorporates:
       *  TriggerPort: '<S8>/Function'
       */
      /* Outputs for Triggered SubSystem: '<S7>/fuel_controller_mode_10ms' incorporates:
       *  TriggerPort: '<S9>/Function'
       */
      /* Outputs for Triggered SubSystem: '<S7>/fuel_controller_pwon' incorporates:
       *  TriggerPort: '<S10>/Function'
       */
      if (rtsiIsModeUpdateTimeStep(&(&AbstractFuelControl_M1_M)->solverInfo)) {
        zcEvent = rt_ZCFcn(RISING_ZERO_CROSSING,
                           &AbstractFuelControl_M1_PrevZCX.fuel_controller_pwon_Trig_ZCE,
                           (AbstractFuelControl_M1_B.Pwon));
        if (zcEvent != NO_ZCEVENT) {
          /* DataStoreWrite: '<S10>/DataStoreWrite1' incorporates:
           *  Constant: '<S10>/Constant1'
           *  DataTypeConversion: '<S10>/Data Type Conversion'
           */
          AbstractFuelControl_M1_DW.controller_mode =
            (AbstractFuelControl_M1_P.Constant1_Value_l != 0.0F);

          /* DataStoreWrite: '<S10>/DataStoreWrite' incorporates:
           *  Constant: '<S10>/Constant2'
           */
          AbstractFuelControl_M1_DW.commanded_fuel =
            AbstractFuelControl_M1_P.Constant2_Value_k;

          /* DataStoreWrite: '<S10>/DataStoreWrite2' incorporates:
           *  Constant: '<S10>/Constant3'
           */
          AbstractFuelControl_M1_DW.airbyfuel_ref =
            AbstractFuelControl_M1_P.Constant3_Value_o;
        }

        zcEvent = rt_ZCFcn(RISING_ZERO_CROSSING,
                           &AbstractFuelControl_M1_PrevZCX.fuel_controller_mode_10ms_Trig_,
                           (AbstractFuelControl_M1_B.PulseGenerator_10ms));
        if (zcEvent != NO_ZCEVENT) {
          /* Outputs for Atomic SubSystem: '<S9>/sensor_failure_detection' */
          /* Logic: '<S16>/Logical Operator' incorporates:
           *  Constant: '<S16>/threshold'
           *  DataStoreRead: '<S9>/DataStoreRead2'
           *  RelationalOperator: '<S16>/Relational Operator'
           *  UnitDelay: '<S16>/Unit Delay'
           */
          rtb_LogicalOperator = ((AbstractFuelControl_M1_DW.airbyfuel_meas <=
            AbstractFuelControl_M1_P.threshold_Value) ||
            AbstractFuelControl_M1_DW.UnitDelay_DSTATE);

          /* Update for UnitDelay: '<S16>/Unit Delay' */
          AbstractFuelControl_M1_DW.UnitDelay_DSTATE = rtb_LogicalOperator;

          /* End of Outputs for SubSystem: '<S9>/sensor_failure_detection' */

          /* Outputs for Atomic SubSystem: '<S9>/normal_mode_detection' */
          /* Sum: '<S14>/Sum' incorporates:
           *  Constant: '<S14>/sampling_sec'
           *  UnitDelay: '<S14>/Unit Delay2'
           */
          rtb_Sum_o = AbstractFuelControl_M1_DW.UnitDelay2_DSTATE +
            AbstractFuelControl_M1_P.sampling_sec_Value;

          /* Logic: '<S14>/Logical Operator' incorporates:
           *  Constant: '<S14>/normal_mode_start_sec'
           *  RelationalOperator: '<S14>/Relational Operator'
           *  UnitDelay: '<S14>/Unit Delay1'
           */
          rtb_LogicalOperator_c = ((rtb_Sum_o >=
            AbstractFuelControl_M1_P.normal_mode_start_sec_Value) ||
            AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_e);

          /* Update for UnitDelay: '<S14>/Unit Delay2' */
          AbstractFuelControl_M1_DW.UnitDelay2_DSTATE = rtb_Sum_o;

          /* Update for UnitDelay: '<S14>/Unit Delay1' */
          AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_e = rtb_LogicalOperator_c;

          /* End of Outputs for SubSystem: '<S9>/normal_mode_detection' */

          /* Outputs for Atomic SubSystem: '<S9>/power_mode_detection' */
          /* Switch: '<S15>/Switch' incorporates:
           *  Constant: '<S15>/Constant'
           *  Constant: '<S15>/Constant1'
           *  Sum: '<S15>/Sum'
           *  UnitDelay: '<S15>/Unit Delay1'
           */
          if (AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_a) {
            rtb_Sum_o = AbstractFuelControl_M1_P.Constant_Value_d;
          } else {
            rtb_Sum_o = AbstractFuelControl_M1_P.Constant_Value_d +
              AbstractFuelControl_M1_P.Constant1_Value_f;
          }

          /* RelationalOperator: '<S15>/Relational Operator1' incorporates:
           *  DataStoreRead: '<S9>/DataStoreRead4'
           *  Switch: '<S15>/Switch'
           */
          rtb_RelationalOperator1 = (AbstractFuelControl_M1_DW.throttle_angle >=
            rtb_Sum_o);

          /* Update for UnitDelay: '<S15>/Unit Delay1' */
          AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_a =
            rtb_RelationalOperator1;

          /* End of Outputs for SubSystem: '<S9>/power_mode_detection' */

          /* DataStoreWrite: '<S9>/DataStoreWrite' incorporates:
           *  Logic: '<S9>/Logical Operator1'
           *  Logic: '<S9>/Logical Operator2'
           */
          AbstractFuelControl_M1_DW.controller_mode = (rtb_LogicalOperator ||
            (!rtb_LogicalOperator_c) || rtb_RelationalOperator1);

          /* Switch: '<S9>/Switch' incorporates:
           *  Logic: '<S9>/Logical Operator3'
           */
          if (rtb_LogicalOperator_c && rtb_RelationalOperator1) {
            /* DataStoreWrite: '<S9>/DataStoreWrite1' incorporates:
             *  Constant: '<S9>/airbyfuel_reference_power'
             */
            AbstractFuelControl_M1_DW.airbyfuel_ref =
              AbstractFuelControl_M1_P.airbyfuel_reference_power_Value;
          } else {
            /* DataStoreWrite: '<S9>/DataStoreWrite1' incorporates:
             *  Constant: '<S9>/airbyfuel_reference'
             */
            AbstractFuelControl_M1_DW.airbyfuel_ref =
              AbstractFuelControl_M1_P.airbyfuel_reference_Value;
          }

          /* End of Switch: '<S9>/Switch' */
        }

        zcEvent = rt_ZCFcn(RISING_ZERO_CROSSING,
                           &AbstractFuelControl_M1_PrevZCX.fuel_controller_10ms_Trig_ZCE,
                           (AbstractFuelControl_M1_B.PulseGenerator_10ms));
        if (zcEvent != NO_ZCEVENT) {
          /* Outputs for Atomic SubSystem: '<S8>/air_estimation' */
          /* Sum: '<S11>/Sum3' incorporates:
           *  Constant: '<S11>/Constant2'
           *  Constant: '<S11>/Constant3'
           *  Constant: '<S11>/Constant4'
           *  Constant: '<S11>/Constant5'
           *  DataStoreRead: '<S8>/DataStoreRead1'
           *  Product: '<S11>/Prod2'
           *  Product: '<S11>/Prod3'
           *  Product: '<S11>/Prod4'
           *  UnitDelay: '<S11>/UnitDelay1'
           */
          rtb_Sum3 = ((AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_d *
                       AbstractFuelControl_M1_DW.engine_speed *
                       AbstractFuelControl_M1_P.Constant3_Value_c +
                       AbstractFuelControl_M1_P.Constant2_Value_d) +
                      AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_d *
                      AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_d *
                      AbstractFuelControl_M1_DW.engine_speed *
                      AbstractFuelControl_M1_P.Constant4_Value) +
            AbstractFuelControl_M1_DW.engine_speed *
            AbstractFuelControl_M1_DW.engine_speed *
            AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_d *
            AbstractFuelControl_M1_P.Constant5_Value;

          /* Update for UnitDelay: '<S11>/UnitDelay1' incorporates:
           *  Constant: '<S11>/Constant1'
           *  DataStoreRead: '<S8>/DataStoreRead'
           *  Gain: '<S11>/Gain'
           *  Product: '<S11>/Prod1'
           *  Sum: '<S11>/Sum1'
           *  Sum: '<S11>/Sum2'
           */
          AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_d +=
            (AbstractFuelControl_M1_DW.throttle_flow - rtb_Sum3) *
            AbstractFuelControl_M1_P.Gain_Gain_j *
            AbstractFuelControl_M1_P.Constant1_Value;

          /* End of Outputs for SubSystem: '<S8>/air_estimation' */

          /* Outputs for Enabled SubSystem: '<S8>/feedback_PI_controller' incorporates:
           *  EnablePort: '<S12>/Enable'
           */
          /* Switch: '<S8>/Switch' incorporates:
           *  Constant: '<S8>/Constant3'
           *  DataStoreRead: '<S8>/DataStoreRead3'
           *  Logic: '<S8>/Logical Operator2'
           */
          if (!AbstractFuelControl_M1_DW.controller_mode) {
            /* Sum: '<S12>/Sum1' incorporates:
             *  DataStoreRead: '<S8>/DataStoreRead2'
             *  DataStoreRead: '<S8>/DataStoreRead4'
             */
            rtb_Sum_o = AbstractFuelControl_M1_DW.airbyfuel_meas -
              AbstractFuelControl_M1_DW.airbyfuel_ref;

            /* Sum: '<S12>/Sum2' incorporates:
             *  Constant: '<S12>/Constant1'
             *  Gain: '<S12>/Gain1'
             *  Product: '<S12>/Prod1'
             *  UnitDelay: '<S12>/UnitDelay1'
             */
            rtb_Sum2 = AbstractFuelControl_M1_P.ki * rtb_Sum_o *
              AbstractFuelControl_M1_P.Constant1_Value_h +
              AbstractFuelControl_M1_DW.UnitDelay1_DSTATE;

            /* Update for UnitDelay: '<S12>/UnitDelay1' */
            AbstractFuelControl_M1_DW.UnitDelay1_DSTATE = rtb_Sum2;

            /* Sum: '<S8>/Sum1' incorporates:
             *  Constant: '<S8>/Constant2'
             *  Gain: '<S12>/Gain'
             *  Sum: '<S12>/Sum3'
             */
            rtb_Sum_o = (AbstractFuelControl_M1_P.kp * rtb_Sum_o + rtb_Sum2) +
              AbstractFuelControl_M1_P.Constant2_Value;

            /* Saturate: '<S8>/fb_fuel_saturation' */
            if (rtb_Sum_o > AbstractFuelControl_M1_P.fb_fuel_saturation_UpperSat)
            {
              rtb_Sum_o = AbstractFuelControl_M1_P.fb_fuel_saturation_UpperSat;
            } else if (rtb_Sum_o <
                       AbstractFuelControl_M1_P.fb_fuel_saturation_LowerSat) {
              rtb_Sum_o = AbstractFuelControl_M1_P.fb_fuel_saturation_LowerSat;
            }

            /* End of Saturate: '<S8>/fb_fuel_saturation' */
          } else {
            rtb_Sum_o = AbstractFuelControl_M1_P.Constant3_Value;
          }

          /* End of Switch: '<S8>/Switch' */
          /* End of Outputs for SubSystem: '<S8>/feedback_PI_controller' */

          /* Outputs for Atomic SubSystem: '<S8>/feedforward_controller' */
          /* Product: '<S8>/Prod1' incorporates:
           *  DataStoreRead: '<S8>/DataStoreRead4'
           *  Product: '<S13>/Product'
           */
          rtb_Sum_o *= rtb_Sum3 / AbstractFuelControl_M1_DW.airbyfuel_ref;

          /* End of Outputs for SubSystem: '<S8>/feedforward_controller' */

          /* Saturate: '<S8>/fuel_saturation' */
          if (rtb_Sum_o > AbstractFuelControl_M1_P.fuel_saturation_UpperSat) {
            /* DataStoreWrite: '<S8>/DataStoreWrite' */
            AbstractFuelControl_M1_DW.commanded_fuel =
              AbstractFuelControl_M1_P.fuel_saturation_UpperSat;
          } else if (rtb_Sum_o <
                     AbstractFuelControl_M1_P.fuel_saturation_LowerSat) {
            /* DataStoreWrite: '<S8>/DataStoreWrite' */
            AbstractFuelControl_M1_DW.commanded_fuel =
              AbstractFuelControl_M1_P.fuel_saturation_LowerSat;
          } else {
            /* DataStoreWrite: '<S8>/DataStoreWrite' */
            AbstractFuelControl_M1_DW.commanded_fuel = rtb_Sum_o;
          }

          /* End of Saturate: '<S8>/fuel_saturation' */
        }
      }

      /* End of Outputs for SubSystem: '<S7>/fuel_controller_pwon' */
      /* End of Outputs for SubSystem: '<S7>/fuel_controller_mode_10ms' */
      /* End of Outputs for SubSystem: '<S7>/fuel_controller_10ms' */
    }

    /* End of Outputs for SubSystem: '<S2>/fuel_controller' */
    if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M)) &&
        (&AbstractFuelControl_M1_M)->Timing.TaskCounters.TID[2] == 0) {
      /* DataStoreRead: '<S2>/DataStoreRead' */
      rtb_DataStoreRead = AbstractFuelControl_M1_DW.commanded_fuel;

      /* DataStoreRead: '<S2>/DataStoreRead1' */
      rtb_DataStoreRead1_l = AbstractFuelControl_M1_DW.controller_mode;

      /* DataStoreRead: '<S2>/DataStoreRead2' */
      rtb_DataStoreRead2 = AbstractFuelControl_M1_DW.airbyfuel_ref;
    }

    /* End of Outputs for SubSystem: '<S1>/AF_Controller' */

    /* Integrator: '<S18>/Integrator' */
    rtb_Integrator = AbstractFuelControl_M1_X.Integrator_CSTATE_h;

    /* Gain: '<S19>/Gain' incorporates:
     *  Integrator: '<S18>/Integrator'
     *  Integrator: '<S19>/Integrator'
     *  Sum: '<S19>/Sum'
     */
    AbstractFuelControl_M1_B.Gain =
      (AbstractFuelControl_M1_X.Integrator_CSTATE_h -
       AbstractFuelControl_M1_X.Integrator_CSTATE) *
      AbstractFuelControl_M1_P.Gain_Gain;

    /* Gain: '<S4>/Gain2' incorporates:
     *  Fcn: '<S4>/Pumping'
     *  Integrator: '<S4>/p0 = 0.543 (bar)'
     */
    rtb_Gain2 = (((0.08979 * AbstractFuelControl_M1_X.p00543bar_CSTATE *
                   rtb_radstorpm - 0.366) - 0.0337 * rtb_radstorpm *
                  AbstractFuelControl_M1_X.p00543bar_CSTATE *
                  AbstractFuelControl_M1_X.p00543bar_CSTATE) + 0.0001 *
                 AbstractFuelControl_M1_X.p00543bar_CSTATE * rtb_radstorpm *
                 rtb_radstorpm) * AbstractFuelControl_M1_P.pump_tol;

    /* Gain: '<S6>/rad//s to rpm' */
    rtb_fuelsystemtransportdelay = AbstractFuelControl_M1_P.radstorpm_Gain *
      rtb_radstorpm;

    /* Gain: '<S3>/Gain' incorporates:
     *  Product: '<S3>/Product1'
     */
    rtb_c_air_chrggcyl = rtb_Gain2 / rtb_radstorpm *
      AbstractFuelControl_M1_P.Gain_Gain_l;

    /* Gain: '<S6>/Kappa tolerance [0.9 1.1]' incorporates:
     *  Gain: '<S3>/Gain'
     *  Lookup_n-D: '<S6>/1-Kappa'
     *  VariableTransportDelay: '<S3>/fuel system transport delay'
     */
    rtb_Kappatolerance0911 = AbstractFuelControl_M1_P.kappa_tol * look2_binlxpw
      (rtb_fuelsystemtransportdelay, rtb_c_air_chrggcyl,
       AbstractFuelControl_M1_P.uKappa_bp01Data,
       AbstractFuelControl_M1_P.uKappa_bp02Data,
       AbstractFuelControl_M1_P.uKappa_tableData,
       AbstractFuelControl_M1_P.uKappa_maxIndex, 5U);
    if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M)) &&
        (&AbstractFuelControl_M1_M)->Timing.TaskCounters.TID[2] == 0) {
      /* Gain: '<S1>/fuel injector tolerance [0.95 1.05]' incorporates:
       *  DataTypeConversion: '<S1>/Data Type Conversion'
       */
      AbstractFuelControl_M1_B.fuelinjectortolerance095105 =
        AbstractFuelControl_M1_P.fuel_inj_tol * rtb_DataStoreRead;
    }

    /* Product: '<S6>/Divide2' incorporates:
     *  Gain: '<S3>/Gain'
     *  Gain: '<S6>/tau_ww tolerance [0.9 1.1]'
     *  Integrator: '<S6>/Integrator'
     *  Lookup_n-D: '<S6>/tau_ww'
     *  VariableTransportDelay: '<S3>/fuel system transport delay'
     */
    rtb_fuel_puddle_evap = AbstractFuelControl_M1_X.Integrator_CSTATE_c /
      (AbstractFuelControl_M1_P.tau_ww_tol * look2_binlxpw
       (rtb_fuelsystemtransportdelay, rtb_c_air_chrggcyl,
        AbstractFuelControl_M1_P.tau_ww_bp01Data,
        AbstractFuelControl_M1_P.tau_ww_bp02Data,
        AbstractFuelControl_M1_P.tau_ww_tableData,
        AbstractFuelControl_M1_P.tau_ww_maxIndex, 5U));

    /* Product: '<S3>/Divide' incorporates:
     *  Product: '<S6>/Divide'
     *  Sum: '<S6>/Add'
     */
    AbstractFuelControl_M1_B.airbyfuel = rtb_Gain2 / (rtb_Kappatolerance0911 *
      AbstractFuelControl_M1_B.fuelinjectortolerance095105 +
      rtb_fuel_puddle_evap);

    /* VariableTransportDelay: '<S3>/fuel system transport delay' */
    rtb_fuelsystemtransportdelay = rt_VTDelayfindtDInterpolate
      (AbstractFuelControl_M1_X.fuelsystemtransportdelay_CSTATE,
       static_cast<real_T *>
       (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_PWORK[0]),
       AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[3],
       AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1],
       AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[0],
       &AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[2],
       (&AbstractFuelControl_M1_M)->Timing.t[0],
       AbstractFuelControl_M1_DW.fuelsystemtransportdelay_RWORK[0],false,
       rtmIsMinorTimeStep((&AbstractFuelControl_M1_M)) && ((static_cast<real_T *>
         (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_PWORK[0]))
        [AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1] +
        AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[3]] ==
        (&AbstractFuelControl_M1_M)->Timing.t[0]),
       AbstractFuelControl_M1_P.fuelsystemtransportdelay_InitOu,&appliedDelay);

    /* Gain: '<S18>/Gain1' incorporates:
     *  Sum: '<S18>/Sum'
     */
    AbstractFuelControl_M1_B.Gain1 = (rtb_fuelsystemtransportdelay -
      rtb_Integrator) * AbstractFuelControl_M1_P.Gain1_Gain;

    /* Lookup_n-D: '<S3>/delay (s)' incorporates:
     *  Gain: '<S3>/Gain'
     *  Gain: '<S3>/rad//s to rpm'
     */
    AbstractFuelControl_M1_B.delays = look2_binlxpw
      (AbstractFuelControl_M1_P.radstorpm_Gain_p * rtb_radstorpm,
       rtb_c_air_chrggcyl, AbstractFuelControl_M1_P.delays_bp01Data,
       AbstractFuelControl_M1_P.delays_bp02Data,
       AbstractFuelControl_M1_P.delays_tableData,
       AbstractFuelControl_M1_P.delays_maxIndex, 5U);

    /* Gain: '<S4>/RT//Vm' incorporates:
     *  Sum: '<S4>/Sum'
     */
    AbstractFuelControl_M1_B.RTVm = (rtb_pratio - rtb_Gain2) *
      AbstractFuelControl_M1_P.RTVm_Gain;

    /* Sum: '<S6>/Sum' incorporates:
     *  Constant: '<S6>/Constant'
     *  Gain: '<S6>/Gain'
     *  Product: '<S6>/Divide1'
     *  Sum: '<S6>/Add2'
     */
    AbstractFuelControl_M1_B.Sum = (AbstractFuelControl_M1_P.Gain_Gain_m *
      rtb_Kappatolerance0911 + AbstractFuelControl_M1_P.Constant_Value) *
      AbstractFuelControl_M1_B.fuelinjectortolerance095105 -
      rtb_fuel_puddle_evap;

    /* End of Outputs for SubSystem: '<Root>/Model 1' */
    if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M)) &&
        (&AbstractFuelControl_M1_M)->Timing.TaskCounters.TID[2] == 0) {
      /* Outport: '<Root>/AFref' incorporates:
       *  DataTypeConversion: '<Root>/Data Type Conversion'
       */
      AbstractFuelControl_M1_Y.AFref = rtb_DataStoreRead2;

      /* Outport: '<Root>/controller_mode' incorporates:
       *  DataTypeConversion: '<Root>/Data Type Conversion1'
       */
      AbstractFuelControl_M1_Y.controller_mode = rtb_DataStoreRead1_l;
    }

    /* Outport: '<Root>/AF' */
    AbstractFuelControl_M1_Y.AF = rtb_AF_sensor;
  }

  if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M))) {
    /* Matfile logging */
    rt_UpdateTXYLogVars((&AbstractFuelControl_M1_M)->rtwLogInfo,
                        ((&AbstractFuelControl_M1_M)->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M))) {
    boolean_T bufferFull;

    /* Update for Atomic SubSystem: '<Root>/Model 1' */
    /* Update for VariableTransportDelay: '<S3>/fuel system transport delay' */
    bufferFull = false;
    if (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1] <
        AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[3] - 1) {
      AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1]++;
    } else {
      AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1] = 0;
    }

    if (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[0] ==
        AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1]) {
      bufferFull = true;
      if (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[0] <
          AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[3] - 1) {
        AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[0]++;
      } else {
        AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[0] = 0;
      }
    }

    (static_cast<real_T *>
      (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_PWORK[0]))
      [AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1]] =
      AbstractFuelControl_M1_B.airbyfuel;
    (static_cast<real_T *>
      (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_PWORK[0]))
      [AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1] +
      AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[3]] =
      (&AbstractFuelControl_M1_M)->Timing.t[0];
    (static_cast<real_T *>
      (AbstractFuelControl_M1_DW.fuelsystemtransportdelay_PWORK[0]))
      [(AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[3] << 1) +
      AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1]] =
      AbstractFuelControl_M1_X.fuelsystemtransportdelay_CSTATE;
    if (bufferFull) {
      rtsiSetBlockStateForSolverChangedAtMajorStep(&(&AbstractFuelControl_M1_M
        )->solverInfo, true);
      rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
        (&(&AbstractFuelControl_M1_M)->solverInfo, true);
    }

    /* End of Update for VariableTransportDelay: '<S3>/fuel system transport delay' */
    /* End of Update for SubSystem: '<Root>/Model 1' */

    /* ContTimeOutputInconsistentWithStateAtMajorOutputFlag is set, need to run a minor output */
    if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M))) {
      if (rtsiGetContTimeOutputInconsistentWithStateAtMajorStep
          (&(&AbstractFuelControl_M1_M)->solverInfo)) {
        rtsiSetSimTimeStep(&(&AbstractFuelControl_M1_M)->solverInfo,
                           MINOR_TIME_STEP);
        rtsiSetContTimeOutputInconsistentWithStateAtMajorStep
          (&(&AbstractFuelControl_M1_M)->solverInfo, false);
        AbstractFuelControl_M1::step();
        rtsiSetSimTimeStep(&(&AbstractFuelControl_M1_M)->solverInfo,
                           MAJOR_TIME_STEP);
      }
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep((&AbstractFuelControl_M1_M))) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal((&AbstractFuelControl_M1_M))!=-1) &&
          !((rtmGetTFinal((&AbstractFuelControl_M1_M))-
             ((((&AbstractFuelControl_M1_M)->Timing.clockTick1+
                (&AbstractFuelControl_M1_M)->Timing.clockTickH1* 4294967296.0)) *
              0.001)) > ((((&AbstractFuelControl_M1_M)->Timing.clockTick1+
                           (&AbstractFuelControl_M1_M)->Timing.clockTickH1*
                           4294967296.0)) * 0.001) * (DBL_EPSILON))) {
        rtmSetErrorStatus((&AbstractFuelControl_M1_M), "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&(&AbstractFuelControl_M1_M)->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++(&AbstractFuelControl_M1_M)->Timing.clockTick0)) {
      ++(&AbstractFuelControl_M1_M)->Timing.clockTickH0;
    }

    (&AbstractFuelControl_M1_M)->Timing.t[0] = rtsiGetSolverStopTime
      (&(&AbstractFuelControl_M1_M)->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      (&AbstractFuelControl_M1_M)->Timing.clockTick1++;
      if (!(&AbstractFuelControl_M1_M)->Timing.clockTick1) {
        (&AbstractFuelControl_M1_M)->Timing.clockTickH1++;
      }
    }

    rate_scheduler((&AbstractFuelControl_M1_M));
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void AbstractFuelControl_M1::AbstractFuelControl_M1_derivatives()
{
  XDot_AbstractFuelControl_M1_T *_rtXdot;
  real_T instantDelay;
  _rtXdot = ((XDot_AbstractFuelControl_M1_T *) (&AbstractFuelControl_M1_M)
             ->derivs);

  /* Derivatives for Atomic SubSystem: '<Root>/Model 1' */
  /* Derivatives for Integrator: '<S19>/Integrator' */
  _rtXdot->Integrator_CSTATE = AbstractFuelControl_M1_B.Gain;

  /* Derivatives for TransferFcn: '<S1>/Throttle delay' incorporates:
   *  Inport: '<Root>/Pedal Angle'
   */
  _rtXdot->Throttledelay_CSTATE = AbstractFuelControl_M1_P.Throttledelay_A *
    AbstractFuelControl_M1_X.Throttledelay_CSTATE;
  _rtXdot->Throttledelay_CSTATE += AbstractFuelControl_M1_U.PedalAngle;

  /* Derivatives for Integrator: '<S4>/p0 = 0.543 (bar)' */
  _rtXdot->p00543bar_CSTATE = AbstractFuelControl_M1_B.RTVm;

  /* Derivatives for Integrator: '<S18>/Integrator' */
  _rtXdot->Integrator_CSTATE_h = AbstractFuelControl_M1_B.Gain1;

  /* Derivatives for Integrator: '<S6>/Integrator' */
  _rtXdot->Integrator_CSTATE_c = AbstractFuelControl_M1_B.Sum;

  /* Derivatives for VariableTransportDelay: '<S3>/fuel system transport delay' */
  instantDelay = AbstractFuelControl_M1_B.delays;
  if (AbstractFuelControl_M1_B.delays >
      AbstractFuelControl_M1_P.fuelsystemtransportdelay_MaxDel) {
    instantDelay = AbstractFuelControl_M1_P.fuelsystemtransportdelay_MaxDel;
  }

  if (instantDelay < 0.0) {
    _rtXdot->fuelsystemtransportdelay_CSTATE = 0.0;
  } else {
    _rtXdot->fuelsystemtransportdelay_CSTATE = 1.0 / instantDelay;
  }

  /* End of Derivatives for VariableTransportDelay: '<S3>/fuel system transport delay' */
  /* End of Derivatives for SubSystem: '<Root>/Model 1' */
}

/* Model initialize function */
void AbstractFuelControl_M1::initialize()
{
  /* Registration code */
  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
                          &(&AbstractFuelControl_M1_M)->Timing.simTimeStep);
    rtsiSetTPtr(&(&AbstractFuelControl_M1_M)->solverInfo, &rtmGetTPtr
                ((&AbstractFuelControl_M1_M)));
    rtsiSetStepSizePtr(&(&AbstractFuelControl_M1_M)->solverInfo,
                       &(&AbstractFuelControl_M1_M)->Timing.stepSize0);
    rtsiSetdXPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
                 &(&AbstractFuelControl_M1_M)->derivs);
    rtsiSetContStatesPtr(&(&AbstractFuelControl_M1_M)->solverInfo, (real_T **)
                         &(&AbstractFuelControl_M1_M)->contStates);
    rtsiSetNumContStatesPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
      &(&AbstractFuelControl_M1_M)->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&AbstractFuelControl_M1_M)->solverInfo, &(
      &AbstractFuelControl_M1_M)->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
      &(&AbstractFuelControl_M1_M)->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
      &(&AbstractFuelControl_M1_M)->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
      (boolean_T**) &(&AbstractFuelControl_M1_M)->contStateDisabled);
    rtsiSetErrorStatusPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
                          (&rtmGetErrorStatus((&AbstractFuelControl_M1_M))));
    rtsiSetRTModelPtr(&(&AbstractFuelControl_M1_M)->solverInfo,
                      (&AbstractFuelControl_M1_M));
  }

  rtsiSetSimTimeStep(&(&AbstractFuelControl_M1_M)->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&(&AbstractFuelControl_M1_M)->solverInfo,
    false);
  rtsiSetIsContModeFrozen(&(&AbstractFuelControl_M1_M)->solverInfo, false);
  (&AbstractFuelControl_M1_M)->intgData.y = (&AbstractFuelControl_M1_M)->odeY;
  (&AbstractFuelControl_M1_M)->intgData.f[0] = (&AbstractFuelControl_M1_M)->
    odeF[0];
  (&AbstractFuelControl_M1_M)->intgData.f[1] = (&AbstractFuelControl_M1_M)->
    odeF[1];
  (&AbstractFuelControl_M1_M)->intgData.f[2] = (&AbstractFuelControl_M1_M)->
    odeF[2];
  (&AbstractFuelControl_M1_M)->contStates = ((X_AbstractFuelControl_M1_T *)
    &AbstractFuelControl_M1_X);
  (&AbstractFuelControl_M1_M)->contStateDisabled =
    ((XDis_AbstractFuelControl_M1_T *) &AbstractFuelControl_M1_XDis);
  (&AbstractFuelControl_M1_M)->Timing.tStart = (0.0);
  rtsiSetSolverData(&(&AbstractFuelControl_M1_M)->solverInfo, static_cast<void *>
                    (&(&AbstractFuelControl_M1_M)->intgData));
  rtsiSetSolverName(&(&AbstractFuelControl_M1_M)->solverInfo,"ode3");
  rtmSetTPtr((&AbstractFuelControl_M1_M), &(&AbstractFuelControl_M1_M)
             ->Timing.tArray[0]);
  rtmSetTFinal((&AbstractFuelControl_M1_M), 40.0);
  (&AbstractFuelControl_M1_M)->Timing.stepSize0 = 0.001;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (nullptr);
    (&AbstractFuelControl_M1_M)->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo((&AbstractFuelControl_M1_M)->rtwLogInfo, (nullptr));
    rtliSetLogXSignalPtrs((&AbstractFuelControl_M1_M)->rtwLogInfo, (nullptr));
    rtliSetLogT((&AbstractFuelControl_M1_M)->rtwLogInfo, "tout");
    rtliSetLogX((&AbstractFuelControl_M1_M)->rtwLogInfo, "");
    rtliSetLogXFinal((&AbstractFuelControl_M1_M)->rtwLogInfo, "");
    rtliSetLogVarNameModifier((&AbstractFuelControl_M1_M)->rtwLogInfo, "rt_");
    rtliSetLogFormat((&AbstractFuelControl_M1_M)->rtwLogInfo, 0);
    rtliSetLogMaxRows((&AbstractFuelControl_M1_M)->rtwLogInfo, 1000);
    rtliSetLogDecimation((&AbstractFuelControl_M1_M)->rtwLogInfo, 1);

    /*
     * Set pointers to the data and signal info for each output
     */
    {
      static void * rt_LoggedOutputSignalPtrs[3];
      rt_LoggedOutputSignalPtrs[0] = &AbstractFuelControl_M1_Y.AFref;
      rt_LoggedOutputSignalPtrs[1] = &AbstractFuelControl_M1_Y.AF;
      rt_LoggedOutputSignalPtrs[2] = &AbstractFuelControl_M1_Y.controller_mode;
      rtliSetLogYSignalPtrs((&AbstractFuelControl_M1_M)->rtwLogInfo,
                            ((LogSignalPtrsType)rt_LoggedOutputSignalPtrs));
    }

    {
      static int_T rt_LoggedOutputWidths[] {
        1,
        1,
        1
      };

      static int_T rt_LoggedOutputNumDimensions[] {
        1,
        1,
        1
      };

      static int_T rt_LoggedOutputDimensions[] {
        1,
        1,
        1
      };

      static boolean_T rt_LoggedOutputIsVarDims[] {
        0,
        0,
        0
      };

      static void* rt_LoggedCurrentSignalDimensions[] {
        (nullptr),
        (nullptr),
        (nullptr)
      };

      static int_T rt_LoggedCurrentSignalDimensionsSize[] {
        4,
        4,
        4
      };

      static BuiltInDTypeId rt_LoggedOutputDataTypeIds[] {
        SS_DOUBLE,
        SS_DOUBLE,
        SS_DOUBLE
      };

      static int_T rt_LoggedOutputComplexSignals[] {
        0,
        0,
        0
      };

      static RTWPreprocessingFcnPtr rt_LoggingPreprocessingFcnPtrs[] {
        (nullptr),
        (nullptr),
        (nullptr)
      };

      static const char_T *rt_LoggedOutputLabels[]{
        "",
        "",
        "" };

      static const char_T *rt_LoggedOutputBlockNames[]{
        "AbstractFuelControl_M1/AFref",
        "AbstractFuelControl_M1/AF",
        "AbstractFuelControl_M1/controller_mode" };

      static RTWLogDataTypeConvert rt_RTWLogDataTypeConvert[] {
        { 0, SS_DOUBLE, SS_DOUBLE, 0, 0, 0, 1.0, 0, 0.0 },

        { 0, SS_DOUBLE, SS_DOUBLE, 0, 0, 0, 1.0, 0, 0.0 },

        { 0, SS_DOUBLE, SS_DOUBLE, 0, 0, 0, 1.0, 0, 0.0 }
      };

      static RTWLogSignalInfo rt_LoggedOutputSignalInfo[] {
        {
          3,
          rt_LoggedOutputWidths,
          rt_LoggedOutputNumDimensions,
          rt_LoggedOutputDimensions,
          rt_LoggedOutputIsVarDims,
          rt_LoggedCurrentSignalDimensions,
          rt_LoggedCurrentSignalDimensionsSize,
          rt_LoggedOutputDataTypeIds,
          rt_LoggedOutputComplexSignals,
          (nullptr),
          rt_LoggingPreprocessingFcnPtrs,

          { rt_LoggedOutputLabels },
          (nullptr),
          (nullptr),
          (nullptr),

          { rt_LoggedOutputBlockNames },

          { (nullptr) },
          (nullptr),
          rt_RTWLogDataTypeConvert
        }
      };

      rtliSetLogYSignalInfo((&AbstractFuelControl_M1_M)->rtwLogInfo,
                            rt_LoggedOutputSignalInfo);

      /* set currSigDims field */
      rt_LoggedCurrentSignalDimensions[0] = &rt_LoggedOutputWidths[0];
      rt_LoggedCurrentSignalDimensions[1] = &rt_LoggedOutputWidths[1];
      rt_LoggedCurrentSignalDimensions[2] = &rt_LoggedOutputWidths[2];
    }

    rtliSetLogY((&AbstractFuelControl_M1_M)->rtwLogInfo, "yout");
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime((&AbstractFuelControl_M1_M)->rtwLogInfo, 0.0,
    rtmGetTFinal((&AbstractFuelControl_M1_M)), (&AbstractFuelControl_M1_M)
    ->Timing.stepSize0, (&rtmGetErrorStatus((&AbstractFuelControl_M1_M))));

  {
    int_T j;

    /* Start for Atomic SubSystem: '<Root>/Model 1' */
    /* Start for Atomic SubSystem: '<S1>/AF_Controller' */
    /* Start for DiscretePulseGenerator: '<S2>/PulseGenerator_10ms' */
    AbstractFuelControl_M1_DW.clockTickCounter = -10;

    /* Start for Atomic SubSystem: '<S2>/fuel_controller' */
    /* Start for DataStoreMemory: '<S7>/DataStoreMemory' */
    AbstractFuelControl_M1_DW.engine_speed =
      AbstractFuelControl_M1_P.DataStoreMemory_InitialValue;

    /* Start for DataStoreMemory: '<S7>/DataStoreMemory1' */
    AbstractFuelControl_M1_DW.throttle_flow =
      AbstractFuelControl_M1_P.DataStoreMemory1_InitialValue;

    /* Start for DataStoreMemory: '<S7>/DataStoreMemory2' */
    AbstractFuelControl_M1_DW.airbyfuel_meas =
      AbstractFuelControl_M1_P.DataStoreMemory2_InitialValue;

    /* Start for DataStoreMemory: '<S7>/DataStoreMemory3' */
    AbstractFuelControl_M1_DW.throttle_angle =
      AbstractFuelControl_M1_P.DataStoreMemory3_InitialValue;

    /* End of Start for SubSystem: '<S2>/fuel_controller' */

    /* Start for DataStoreMemory: '<S2>/commanded_fuel' */
    AbstractFuelControl_M1_DW.commanded_fuel =
      AbstractFuelControl_M1_P.commanded_fuel_InitialValue;

    /* Start for DataStoreMemory: '<S2>/mode_fb' */
    AbstractFuelControl_M1_DW.controller_mode =
      AbstractFuelControl_M1_P.mode_fb_InitialValue;

    /* Start for DataStoreMemory: '<S2>/mode_fb1' */
    AbstractFuelControl_M1_DW.airbyfuel_ref =
      AbstractFuelControl_M1_P.mode_fb1_InitialValue;

    /* End of Start for SubSystem: '<S1>/AF_Controller' */

    /* Start for VariableTransportDelay: '<S3>/fuel system transport delay' */
    AbstractFuelControl_M1_DW.fuelsystemtransportdelay_RWORK[0] = 0.0;
    AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[0] = 0;
    AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[1] = 0;
    AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[2] = 0;
    AbstractFuelControl_M1_DW.fuelsystemtransportdelay_IWORK[3] = 20480;
    AbstractFuelControl_M1_DW.fuelsystemtransportdelay_PWORK[0] = (void *)
      &AbstractFuelControl_M1_DW.fuelsystemtransportdelay_RWORK[1];
    for (j = 0; j < 20480; j++) {
      AbstractFuelControl_M1_DW.fuelsystemtransportdelay_RWORK[j + 1] = 14.7;
      AbstractFuelControl_M1_DW.fuelsystemtransportdelay_RWORK[j + 20481] =
        (&AbstractFuelControl_M1_M)->Timing.t[0];
      AbstractFuelControl_M1_DW.fuelsystemtransportdelay_RWORK[40961] = 0.0;
    }

    /* End of Start for VariableTransportDelay: '<S3>/fuel system transport delay' */
    /* End of Start for SubSystem: '<Root>/Model 1' */
  }

  AbstractFuelControl_M1_PrevZCX.fuel_controller_10ms_Trig_ZCE =
    UNINITIALIZED_ZCSIG;
  AbstractFuelControl_M1_PrevZCX.fuel_controller_mode_10ms_Trig_ =
    UNINITIALIZED_ZCSIG;
  AbstractFuelControl_M1_PrevZCX.fuel_controller_pwon_Trig_ZCE =
    UNINITIALIZED_ZCSIG;

  /* SystemInitialize for Atomic SubSystem: '<Root>/Model 1' */
  /* InitializeConditions for Integrator: '<S19>/Integrator' */
  AbstractFuelControl_M1_X.Integrator_CSTATE =
    AbstractFuelControl_M1_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S1>/Throttle delay' */
  AbstractFuelControl_M1_X.Throttledelay_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S4>/p0 = 0.543 (bar)' */
  AbstractFuelControl_M1_X.p00543bar_CSTATE =
    AbstractFuelControl_M1_P.p00543bar_IC;

  /* InitializeConditions for Integrator: '<S18>/Integrator' */
  AbstractFuelControl_M1_X.Integrator_CSTATE_h =
    AbstractFuelControl_M1_P.Integrator_IC_l;

  /* InitializeConditions for Integrator: '<S6>/Integrator' */
  AbstractFuelControl_M1_X.Integrator_CSTATE_c =
    AbstractFuelControl_M1_P.Integrator_IC_m;

  /* InitializeConditions for VariableTransportDelay: '<S3>/fuel system transport delay' */
  AbstractFuelControl_M1_X.fuelsystemtransportdelay_CSTATE = 0.0;

  /* SystemInitialize for Atomic SubSystem: '<S1>/AF_Controller' */
  /* SystemInitialize for Atomic SubSystem: '<S2>/fuel_controller' */
  /* SystemInitialize for Triggered SubSystem: '<S7>/fuel_controller_mode_10ms' */
  /* SystemInitialize for Atomic SubSystem: '<S9>/sensor_failure_detection' */
  /* InitializeConditions for UnitDelay: '<S16>/Unit Delay' */
  AbstractFuelControl_M1_DW.UnitDelay_DSTATE =
    AbstractFuelControl_M1_P.UnitDelay_InitialCondition;

  /* End of SystemInitialize for SubSystem: '<S9>/sensor_failure_detection' */

  /* SystemInitialize for Atomic SubSystem: '<S9>/normal_mode_detection' */
  /* InitializeConditions for UnitDelay: '<S14>/Unit Delay2' */
  AbstractFuelControl_M1_DW.UnitDelay2_DSTATE =
    AbstractFuelControl_M1_P.UnitDelay2_InitialCondition;

  /* InitializeConditions for UnitDelay: '<S14>/Unit Delay1' */
  AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_e =
    AbstractFuelControl_M1_P.UnitDelay1_InitialCondition_c;

  /* End of SystemInitialize for SubSystem: '<S9>/normal_mode_detection' */

  /* SystemInitialize for Atomic SubSystem: '<S9>/power_mode_detection' */
  /* InitializeConditions for UnitDelay: '<S15>/Unit Delay1' */
  AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_a =
    AbstractFuelControl_M1_P.UnitDelay1_InitialCondition_f;

  /* End of SystemInitialize for SubSystem: '<S9>/power_mode_detection' */
  /* End of SystemInitialize for SubSystem: '<S7>/fuel_controller_mode_10ms' */

  /* SystemInitialize for Triggered SubSystem: '<S7>/fuel_controller_10ms' */
  /* SystemInitialize for Atomic SubSystem: '<S8>/air_estimation' */
  /* InitializeConditions for UnitDelay: '<S11>/UnitDelay1' */
  AbstractFuelControl_M1_DW.UnitDelay1_DSTATE_d =
    AbstractFuelControl_M1_P.UnitDelay1_InitialCondition;

  /* End of SystemInitialize for SubSystem: '<S8>/air_estimation' */

  /* SystemInitialize for Enabled SubSystem: '<S8>/feedback_PI_controller' */
  /* InitializeConditions for UnitDelay: '<S12>/UnitDelay1' */
  AbstractFuelControl_M1_DW.UnitDelay1_DSTATE =
    AbstractFuelControl_M1_P.UnitDelay1_InitialCondition_l;

  /* End of SystemInitialize for SubSystem: '<S8>/feedback_PI_controller' */
  /* End of SystemInitialize for SubSystem: '<S7>/fuel_controller_10ms' */
  /* End of SystemInitialize for SubSystem: '<S2>/fuel_controller' */
  /* End of SystemInitialize for SubSystem: '<S1>/AF_Controller' */
  /* End of SystemInitialize for SubSystem: '<Root>/Model 1' */
}

/* Model terminate function */
void AbstractFuelControl_M1::terminate()
{
  /* (no terminate code required) */
}

/* Constructor */
AbstractFuelControl_M1::AbstractFuelControl_M1() :
  AbstractFuelControl_M1_U(),
  AbstractFuelControl_M1_Y(),
  AbstractFuelControl_M1_B(),
  AbstractFuelControl_M1_DW(),
  AbstractFuelControl_M1_X(),
  AbstractFuelControl_M1_XDis(),
  AbstractFuelControl_M1_PrevZCX(),
  AbstractFuelControl_M1_M()
{
  /* Currently there is no constructor body generated.*/
}

/* Destructor */
/* Currently there is no destructor body generated.*/
AbstractFuelControl_M1::~AbstractFuelControl_M1() = default;

/* Real-Time Model get method */
RT_MODEL_AbstractFuelControl__T * AbstractFuelControl_M1::getRTM()
{
  return (&AbstractFuelControl_M1_M);
}
