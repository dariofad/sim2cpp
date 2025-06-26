//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: Simulink2Code.h
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
#ifndef Simulink2Code_h_
#define Simulink2Code_h_
#include <cmath>
#include "rtwtypes.h"
#include "Simulink2Code_types.h"

// Class declaration for model Simulink2Code
class Simulink2Code final
{
  // public data and function members
 public:
  // External inputs (root inport signals with default storage)
  struct ExtU_Simulink2Code_T {
    int x;                          // '<Root>/x'
    int y;                          // '<Root>/y'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_Simulink2Code_T {
    int result;                     // '<Root>/result'
  };

  // Real-time Model Data Structure
  struct RT_MODEL_Simulink2Code_T {
    const char_T * volatile errorStatus;
    const char_T* getErrorStatus() const;
    void setErrorStatus(const char_T* const volatile aErrorStatus);
  };

  // Copy Constructor
  Simulink2Code(Simulink2Code const&) = delete;

  // Assignment Operator
  Simulink2Code& operator= (Simulink2Code const&) & = delete;

  // Move Constructor
  Simulink2Code(Simulink2Code &&) = delete;

  // Move Assignment Operator
  Simulink2Code& operator= (Simulink2Code &&) = delete;

  // Real-Time Model get method
  Simulink2Code::RT_MODEL_Simulink2Code_T * getRTM();

  // External inputs
  ExtU_Simulink2Code_T Simulink2Code_U;

  // External outputs
  ExtY_Simulink2Code_T Simulink2Code_Y;

  // model initialize function
  static void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  Simulink2Code();

  // Destructor
  ~Simulink2Code();

  // private data and function members
 private:
  // Real-Time Model
  RT_MODEL_Simulink2Code_T Simulink2Code_M;
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
