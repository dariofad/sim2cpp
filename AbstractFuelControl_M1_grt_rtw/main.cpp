#include "AbstractFuelControl_M1.h"
#include <stdio.h>
#include <chrono>
#include <thread>

static AbstractFuelControl_M1 model;

int main() {

    // slowdown the model to observe the perturbations
    using namespace std::chrono;
    auto interval = milliseconds(500);
    system_clock::time_point next_call_time = system_clock::now() + interval;

    model.initialize();

    ExtU_AbstractFuelControl_M1_T inputs{};
    double pedal_offset = 20.0;
    double rpm_offset = 900.0;

    const int num_steps = 20;
    double ts = 0.001;

    for (int i = 0; i < num_steps; ++i) {

        std::this_thread::sleep_until(next_call_time);	    

        inputs.PedalAngle = pedal_offset + (double)i * 2 / 10;
	inputs.EngineSpeed = rpm_offset + (double)i;
	    
        model.setExternalInputs(&inputs);
        model.step();

        const auto& out = model.getExternalOutputs();
        printf("time [s]: %1.3f, pedal: %f, rpm: %f, AFref: %f, AF: %f, mode: "
               "%f\n",
               i * ts, inputs.PedalAngle, inputs.EngineSpeed, out.AFref, out.AF,
               out.controller_mode);

        next_call_time += interval;
    }

    model.terminate();
    return 0;
}
