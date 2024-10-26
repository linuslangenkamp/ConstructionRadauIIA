#include "integrator.h"
#include <iostream>

int main() {
    std::cout << "Hi" << std::endl;
    auto int78 = Integrator::radauIIA(IntegratorSteps::Steps78);
    return 0;
}