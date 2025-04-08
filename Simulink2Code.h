//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: Simulink2Code.h
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
#ifndef Simulink2Code_h_
#define Simulink2Code_h_
#include <math.h>
#include "rtwtypes.h"
#include "Simulink2Code_types.h"

// Class declaration for model Simulink2Code
class myPlus final
{
  // public data and function members
 public:
  // External inputs (root inport signals with default storage)
  struct ExtU {
    real_T x;                          // '<Root>/x'
    real_T y;                          // '<Root>/y'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY {
    real_T result;                     // '<Root>/result'
  };

  // Real-time Model Data Structure
  struct RT_MODEL {
    const char_T * volatile errorStatus;
    const char_T* getErrorStatus() const;
    void setErrorStatus(const char_T* const volatile aErrorStatus);
  };

  // Copy Constructor
  myPlus(myPlus const&) = delete;

  // Assignment Operator
  myPlus& operator= (myPlus const&) & = delete;

  // Move Constructor
  myPlus(myPlus &&) = delete;

  // Move Assignment Operator
  myPlus& operator= (myPlus &&) = delete;

  // Real-Time Model get method
  myPlus::RT_MODEL * getRTM();

  // External inputs
  ExtU rtU;

  // External outputs
  ExtY rtY;

  // model initialize function
  static void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  myPlus();

  // Destructor
  ~myPlus();

  // private data and function members
 private:
  // Real-Time Model
  RT_MODEL rtM;
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
//  '<Root>' : 'Simulink2Code'

#endif                                 // Simulink2Code_h_

//
// File trailer for generated code.
//
// [EOF]
//
