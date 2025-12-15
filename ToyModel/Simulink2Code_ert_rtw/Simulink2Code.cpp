//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: Simulink2Code.cpp
//
// Code generated for Simulink model 'Simulink2Code'.
//
// Model version                  : 1.19
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Wed Apr  9 13:29:11 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Apple->ARM64
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "Simulink2Code.h"

// Model step function
void Simulink2Code::step()
{
  // Outport: '<Root>/result' incorporates:
  //   Inport: '<Root>/x'
  //   Inport: '<Root>/y'
  //   Sum: '<Root>/Plus'

  Simulink2Code_Y.result = Simulink2Code_U.x + Simulink2Code_U.y + model_offset;
  
}

// Model initialize function
void Simulink2Code::initialize()
{
  // (no initialization code required)
}

// Model terminate function
void Simulink2Code::terminate()
{
  // (no terminate code required)
}

const char_T* Simulink2Code::RT_MODEL_Simulink2Code_T::getErrorStatus() const
{
  return (errorStatus);
}

void Simulink2Code::RT_MODEL_Simulink2Code_T::setErrorStatus(const char_T* const
  volatile aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

// Constructor
Simulink2Code::Simulink2Code() :
  Simulink2Code_U(),
  Simulink2Code_Y(),
  Simulink2Code_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
Simulink2Code::~Simulink2Code() = default;

// Real-Time Model get method
Simulink2Code::RT_MODEL_Simulink2Code_T * Simulink2Code::getRTM()
{
  return (&Simulink2Code_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
