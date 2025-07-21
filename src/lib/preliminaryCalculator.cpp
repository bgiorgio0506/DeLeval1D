#include <cmath>
#include "preliminaryCalculator.h"


double PreliminaryAnalyser::getExitMachNumber(){
    const double gamma = this->intialConditions.gamma;
    const double exponent = (gamma - 1) /gamma;
    const double pressureRatio = pow(this->intialConditions.combustionPressure / this->intialConditions.exitPressure, exponent);
    const double specHeatRatio = 2 / (gamma - 1);


    return this->exitMachNumber = sqrt(specHeatRatio*(pressureRatio-1));
}

double PreliminaryAnalyser::getExitAreaToThroatAreaRatio(){ //todo check fo validity
    const double gamma = this->intialConditions.gamma;
    const double exponent = (gamma + 1) / (2 * (gamma - 1));
    const double mReciprocal = 1 / this->exitMachNumber;
    const double specHeatRatio = 2 / (gamma + 1);
    return this->exitAreaToThroatAreaRatio = mReciprocal * pow(specHeatRatio * (1 + (((gamma - 1) / 2) * pow(this->exitMachNumber, 2))), exponent);
}

double PreliminaryAnalyser::getCStar(){
    const double gamma = this->intialConditions.gamma;
    const double R = 8.314 / (this->intialConditions.molecularWeight/1000); //g/mol -> kg/mol
    const double T = this->intialConditions.combustionTemperature;
    const double exponent = -((gamma + 1) / (2 * (gamma - 1)));

     return this->cStar = (sqrt(gamma*R*T)/gamma)*pow(2/(gamma+1), exponent);
}

double PreliminaryAnalyser::getThrustCoefficient(){
    const double gamma = this->intialConditions.gamma;
    const double pressureRatio = this->intialConditions.exitPressure / this->intialConditions.combustionPressure;
    const double exitAreaToThroatAreaRatio = this->exitAreaToThroatAreaRatio;
    const double pressureCorrected = (this->intialConditions.exitPressure - this->intialConditions.atmopshericPressure) / (this->intialConditions.combustionPressure);
    const double exponent = (gamma + 1) / (gamma - 1);

    return this->thrustCoefficient = (gamma * sqrt(pow(2 / (gamma + 1), exponent) * ((2 / (gamma - 1)) * (1 - pow(pressureRatio, (gamma - 1) / gamma))))) + (pressureCorrected * exitAreaToThroatAreaRatio);
};

double PreliminaryAnalyser::getThroatArea(){
  //use the throat area ratio to calculate the throat area
    return this->throatArea = (this->intialConditions.thrust * 1000) / ((this->intialConditions.combustionPressure * 1e+6) * this->thrustCoefficient);
}

double PreliminaryAnalyser::getMassFlowRate(){
    return this->massFlowRate = ((this->intialConditions.combustionPressure*1e+6) * this->throatArea) / this->cStar;
}

double PreliminaryAnalyser::getSpeficImpulse(){
    return this->speficImpulse = (this->cStar * this->thrustCoefficient) / 9.81;
}

/**TODO: resolve()*/
preliminaryResults PreliminaryAnalyser::resolve(){
    preliminaryResults resolved;
    if(this->intialConditions.exitPressure == 0){
        this->intialConditions.exitPressure = 0.1;
    }
    resolved.intialConditions = this->intialConditions;
    resolved.exitMachNumber = this->getExitMachNumber();
    resolved.exitAreaToThroatAreaRatio = this->getExitAreaToThroatAreaRatio();
    resolved.cStar = this->getCStar();
    resolved.thrustCoefficient = this->getThrustCoefficient();
    resolved.throatArea = this->getThroatArea();
    if(this->throatArea != 0){
        resolved.massFlowRate = this->getMassFlowRate();
        resolved.speficImpulse = this->getSpeficImpulse();
        resolved.thrust = this->intialConditions.thrust;
    }else {
        resolved.massFlowRate = 0;
        resolved.speficImpulse = 0;
        resolved.thrust = 0;
    }

    return resolved;
    
}