//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: main.cpp
//
// Code generated for Simulink model 'automatic_transmission'.
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
#include <stdio.h>         // This example main program uses printf/fflush
#include "automatic_transmission.h"    // Model header file

static automatic_transmission AT_Obj; // Instance of model class

//
// Associating rt_OneStep with a real-time clock or interrupt service routine
// is what makes the generated code "real-time".  The function rt_OneStep is
// always associated with the base rate of the model.  Subrates are managed
// by the base rate from inside the generated code.  Enabling/disabling
// interrupts and floating point context switches are target specific.  This
// example code indicates where these should take place relative to executing
// the generated code step function.  Overrun behavior should be tailored to
// your application needs.  This example simply sets an error status in the
// real-time model and returns from rt_OneStep.
//
void rt_OneStep(void);
void rt_OneStep(void)
{
  static boolean_T OverrunFlag{false};

  // Disable interrupts here

  // Check for overrun
  if (OverrunFlag)
  {
    AT_Obj.getRTM()->setErrorStatus("Overrun");
    return;
  }

  OverrunFlag = true;

  // Save FPU context here (if necessary)
  // Re-enable timer or interrupt here
  // Set model inputs here

  // Step the model
  AT_Obj.step();

  // Get model outputs here

  // Indicate task complete
  OverrunFlag = false;

  // Disable interrupts here
  // Restore FPU context here (if necessary)
  // Enable interrupts here
}

//
// The example main function illustrates what is required by your
// application code to initialize, execute, and terminate the generated code.
// Attaching rt_OneStep to a real-time clock is target specific. This example
// illustrates how you do this relative to initializing the model.
//
int_T main(int_T argc, const char *argv[])
{
  // Unused arguments
  (void)(argc);
  (void)(argv);

  // Initialize model
  AT_Obj.initialize();

  double throttle[10] = {2,2,2,2,2,20,20,20,20,20};
  double brake[10] = {1,1,1,1,1,10,10,10,10,10};

  for (int i = 0; i < 10; i++)
  {
    // Set input signals in steps
    AT_Obj.rtU.throttle = throttle[i];
    AT_Obj.rtU.brake = brake[i];

    // Perform a simulation step
    AT_Obj.step();

    // Print the output of the current step
    double Ts = 0.01; // sampling time
    printf("at time %f, input (throttle, brake): %f, %f; ", i * Ts, AT_Obj.rtU.throttle, AT_Obj.rtU.brake);
    printf("output (speed, RPM, gear): %f, %f, %f;\n", AT_Obj.rtY.speed, AT_Obj.rtY.RPM, AT_Obj.rtY.gear);
  }

  // Terminate model
  // AT_Obj.terminate();
  return 0;
}
