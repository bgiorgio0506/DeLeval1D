#include <iostream>
#include <string>
#include "lib/preliminaryCalculator.h"
#include "lib/engine.h"


//TODO: add function to export data

#pragma region declarations
void oneDAnalysis();
CombustionReactor configAnalysis();
void displayResult(engine eng);
#pragma endregion

int main()
{
    std::cout << "Welcome to Rocket Engine Analyzer\n";
	std::cout << "This program will help you on the first design analisys.\n The program can perform:\n 1. 1D Analysis\n 2. From the 1D Analysis can extrapolate the correct data to solve the geometry nozzle";
    std::cout << "\n\n In the future plans are also to add the ability to solve the combustion at equilibrium, preliminary design of the combustion chamber and cooling system.\n\n";
    oneDAnalysis();
    system("pause");

    return 0;

};


#pragma region rocket1D
void oneDAnalysis()
{
    CombustionReactor combParams = configAnalysis();
    int nLines = 0;
    std::cout << "Enter number of lines for generation: ";
    std::cin >> nLines;
    engine eng(combParams, nLines);
    eng.preliminaryAnalysis();
    displayResult(eng);
};

CombustionReactor configAnalysis() {
    CombustionReactor combParams;
    std::cout << "Rocket Nozzle Preliminary analysis\n";
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

void displayResult(engine eng) {
	bool wantGeometry;
	std::cout << "\n\n========================================\n";
    std::cout << "\nRocket Nozzle 1D analysis results\n";
    std::cout << "Mass flow rate: " << eng.resolved.massFlowRate << " kg/s\n";
    std::cout << "Throat area: " << eng.resolved.throatArea << " m^2\n";
    std::cout << "Exit area to throat area ratio: " << eng.resolved.exitAreaToThroatAreaRatio << "\n";
    std::cout << "Specific impulse: " << eng.resolved.speficImpulse << " s\n";
    std::cout << "Exit Mach number: " << eng.resolved.exitMachNumber << "\n";
    std::cout << "Thrust: " << eng.resolved.thrust << " KN\n";
    std::cout << "C* (Characteristic velocity): " << eng.resolved.cStar << " m/s\n";
    std::cout << "Thrust coefficient: " << eng.resolved.thrustCoefficient << "\n";
	std::cout << "========================================\n";

	std::cout << "Do you want to calculate the geometry of the nozzle? (1 for yes, 0 for no): ";
	std::cin >> wantGeometry;
	if (!wantGeometry) {
		//todo ask for export data
		std::cout << "Do you want to export the analysis data? (1 for yes, 0 for no): ";
		int exportData;
		std::cin >> exportData;
		if (exportData) {
			//todo implement export data
			std::cout << "Export data feature is not implemented yet.\n";
		}
		else {
			std::cout << "Exiting the program.\n";
            return;
		}
		
    }
    else {
		std::cout << "Generating Nozzle Geometry.\n";
        eng.generateGeometry();
        eng.exportResults("nozzleGeometry.csv");
		return;
    }
};

#pragma endregion

