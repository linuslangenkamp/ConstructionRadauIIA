#include "integrator.h"
#include <iostream>

int main() {
    std::cout << "Hi" << std::endl;
    auto integr = Integrator::radauIIA(IntegratorSteps::Steps1);
    return 0;
}