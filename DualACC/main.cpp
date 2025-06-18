#include <iostream>
#include "leadCar.h"
#include "egoCar.h"

int main() {
    leadCar lead;
    egoCar ego;

    const double Ts = 0.1;
    const int steps = 801;

    lead.initialize();
    ego.initialize();

    for (int i = 0; i < steps; ++i) {
        lead.step();
        const leadCar::ExtY_leadCar_T& leadOutput = lead.getExternalOutputs();

        egoCar::ExtU_egoCar_T egoInput;
        egoInput.d_lead = leadOutput.d_lead;
        egoInput.v_lead = leadOutput.v_lead;
        ego.setExternalInputs(&egoInput);

        ego.step();

        const egoCar::ExtY_egoCar_T& egoOutput = ego.getExternalOutputs();

        printf("t = %.2f, a_lead = %.3f, d_lead = %.3f, v_lead = %.3f; v_ego = %.3f, a_ego = %.3f; d_rel = %.3f, v_rel = %.3f\n",
               i * Ts,
               leadOutput.a_lead,
               leadOutput.d_lead,
               leadOutput.v_lead,
               egoOutput.a_ego,
               egoOutput.v_ego,
               egoOutput.d_rel,
               egoOutput.v_rel);
    }

    return 0;
}