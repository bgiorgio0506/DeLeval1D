/**
 * @file rocket1D.h
 * @author Giorgio Bella (bgiorgio0506@gmail.it)
 * @brief This file implements the 1 dimensional rocket analysis
 * @version 0.1
 * @date 2024-02-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <string>

using namespace std;

struct CombustionReactor
{
    string prop;
    double thrust;//kN
    double fuelToOxidizerRatio;
    double combustionPressure;//MPa
    double exitPressure;//MPa
    double oxidizerTemperature;//K
    double fuelTemperature;//K
    double combustionTemperature;//K
    double gamma;// resolved combustion ratio of specific heats
    double molecularWeight;//g/mol
    double atmopshericPressure;//MPa
};

struct preliminaryResults
{
    CombustionReactor intialConditions;
    double massFlowRate;
    double throatArea;
    double exitAreaToThroatAreaRatio;
    double speficImpulse;
    double exitMachNumber;
    double thrust;
    double cStar;
    double thrustCoefficient;
};



class PreliminaryAnalyser{
    public:
        PreliminaryAnalyser(CombustionReactor initializer){
            this->intialConditions = initializer;
            this->cStar = 0.0;
            this->massFlowRate = 0.0;
            this->exitMachNumber = 0.0;
            this->exitAreaToThroatAreaRatio = 0.0;
            this->speficImpulse = 0.0;
            this->throatArea = 0.0;
            this->thrustCoefficient = 0.0;

        };
        ~PreliminaryAnalyser() {};
        CombustionReactor intialConditions;
        double massFlowRate;
        double throatArea;
        double exitAreaToThroatAreaRatio;
        double speficImpulse;
        double exitMachNumber;
        double cStar;
        double thrustCoefficient;
        double getExitMachNumber();
        double getThroatArea();
        double getExitAreaToThroatAreaRatio();
        double getMassFlowRate();
        double getSpeficImpulse();
        double getCStar();
        double getThrustCoefficient();
        preliminaryResults resolve();

};