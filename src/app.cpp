#include <iostream>
#include <string>
#include "lib/rocket1D.h"


//TODO: add function to export data

#pragma region declarations
void oneDAnalysis();
CombustionReactor configAnalysis();
void displayResult(RocketOneDimResolved resolved);
#pragma endregion

int main()
{
    std::cout << "Welcome to Rocket Engine Analyzer\n";
    oneDAnalysis();
    system("pause");

    return 0;

};


#pragma region rocket1D
void oneDAnalysis()
{
    CombustionReactor combParams = configAnalysis();
    RocketOneDimAnalyser rocketAnalyser(combParams);
    RocketOneDimResolved resolved = rocketAnalyser.resolve();
    displayResult(resolved);
};

CombustionReactor configAnalysis() {
    CombustionReactor combParams;
    std::cout << "Rocket Nozzle 1D analysis\n";
    std::cout << "Enter the propellant type: ";
    std::cin >> combParams.prop;
    std::cout << "Enter the thrust (kN): ";
    std::cin >> combParams.thrust;
    std::cout << "Enter the fuel to oxidizer ratio: ";
    std::cin >> combParams.fuelToOxidizerRatio;
    std::cout << "Enter the combustion pressure (MPa): ";
    std::cin >> combParams.combustionPressure;
    std::cout << "Enter the exit pressure (MPa): ";
    std::cin >> combParams.exitPressure;
    std::cout << "Enter the oxidizer temperature (K): ";
    std::cin >> combParams.oxidizerTemperature;
    std::cout << "Enter the fuel temperature (K): ";
    std::cin >> combParams.fuelTemperature;
    std::cout << "Enter the combustion temperature (K): ";
    std::cin >> combParams.combustionTemperature;
    std::cout << "Enter the resolved combustion ratio of specific heats: ";
    std::cin >> combParams.gamma;
    std::cout << "Enter the molecular weight (g/mol): ";
    std::cin >> combParams.molecularWeight;
    std::cout << "Enter the atmospheric pressure (MPa): ";
    std::cin >> combParams.atmopshericPressure;

    return combParams;
};

void displayResult(RocketOneDimResolved resolved) {
    std::cout << "\nRocket Nozzle 1D analysis results\n";
    std::cout << "Mass flow rate: " << resolved.massFlowRate << " kg/s\n";
    std::cout << "Throat area: " << resolved.throatArea << " m^2\n";
    std::cout << "Exit area to throat area ratio: " << resolved.exitAreaToThroatAreaRatio << "\n";
    std::cout << "Specific impulse: " << resolved.speficImpulse << " s\n";
    std::cout << "Exit Mach number: " << resolved.exitMachNumber << "\n";
    std::cout << "Thrust: " << resolved.thrust << " KN\n";
    std::cout << "C* (Characteristic velocity): " << resolved.cStar << " m/s\n";
    std::cout << "Thrust coefficient: " << resolved.thrustCoefficient << "\n";
};

#pragma endregion

