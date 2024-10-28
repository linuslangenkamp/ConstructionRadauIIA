#include "integrator.h"
#include <iostream>

int main() {
    std::cout << "Hi" << std::endl;
    for (int i = 0; i <= 70; i++) {
        auto integr = Integrator::radauIIA((IntegratorSteps)i);
        std::cout << integr.steps << std::endl;
    }
    
    return 0;
}