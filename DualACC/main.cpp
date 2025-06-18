#include <iostream>
#include "leadCar.h"
#include "egoCar.h"

int main() {
    leadCar lead;
    egoCar ego;

    lead.initialize();
    ego.initialize();

    for (int i = 0; i < 100; ++i) {
        lead.step();
        const leadCar::ExtY_leadCar_T& leadOutput = lead.getExternalOutputs();

        egoCar::ExtU_egoCar_T egoInput;
        egoInput.d_lead = leadOutput.d_lead;
        egoInput.v_lead = leadOutput.v_lead;
        ego.setExternalInputs(&egoInput);

        ego.step();

        const egoCar::ExtY_egoCar_T& egoOutput = ego.getExternalOutputs();

        std::cout << "Step " << i
          << ": a_lead = " << leadOutput.a_lead
          << ", d_lead = " << leadOutput.d_lead
          << ", v_lead = " << leadOutput.v_lead
          << " | a_ego = " << egoOutput.a_ego
          << ", v_ego = " << egoOutput.v_ego
          << " | d_rel = " << egoOutput.d_rel
          << ", v_rel = " << egoOutput.v_rel
          << std::endl;    }

    return 0;
}