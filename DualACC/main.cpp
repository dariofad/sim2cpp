#include <iostream>
#include "leadCar.h"
#include "egoCar.h"

#include <limits.h>

static leadCar lead;
static egoCar ego;

int main() {

    const double Ts = 0.1;
    const int steps = 801;

    lead.initialize();
    ego.initialize();

    for (int i = 0; i < INT_MAX; ++i) {
        lead.step();
        const leadCar::ExtY_leadCar_T& leadOutput = lead.getExternalOutputs();

        egoCar::ExtU_egoCar_T egoInput;
        egoInput.d_lead = leadOutput.d_lead;
        egoInput.v_lead = leadOutput.v_lead;
        ego.setExternalInputs(&egoInput);

        ego.step();

        const egoCar::ExtY_egoCar_T& egoOutput = ego.getExternalOutputs();

        printf("After step, t = %.2f, a_ego = %.3f, v_ego = %.3f, d_rel = %.3f\n",
               i * Ts,
               egoOutput.a_ego,
               egoOutput.v_ego,
               egoOutput.d_rel);
    }

    return 0;
}
