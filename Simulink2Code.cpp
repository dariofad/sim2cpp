//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: Simulink2Code.cpp
//
// Code generated for Simulink model 'Simulink2Code'.
//
// Model version                  : 1.13
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Mon Apr  7 20:57:27 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Apple->ARM64
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "Simulink2Code.h"

// Model step function
void myPlus::step()
{
  // Outport: '<Root>/result' incorporates:
  //   Inport: '<Root>/x'
  //   Inport: '<Root>/y'
  //   Sum: '<Root>/Plus'

  rtY.result = rtU.x + rtU.y;
}

// Model initialize function
void myPlus::initialize()
{
  // (no initialization code required)
}

// Model terminate function
void myPlus::terminate()
{
  // (no terminate code required)
}

const char_T* myPlus::RT_MODEL::getErrorStatus() const
{
  return (errorStatus);
}

void myPlus::RT_MODEL::setErrorStatus(const char_T* const volatile aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

// Constructor
myPlus::myPlus() :
  rtU(),
  rtY(),
  rtM()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
myPlus::~myPlus() = default;

// Real-Time Model get method
myPlus::RT_MODEL * myPlus::getRTM()
{
  return (&rtM);
}

//
// File trailer for generated code.
//
// [EOF]
//
