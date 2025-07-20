#define _USE_MATH_DEFINES

#include "flowEq.h"
#include <math.h>


double FlowAnalyzer::getPmAngle(double mach, double gamma) {
	if (gamma <= 1) return 0;
	double angle = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt(((gamma - 1) * (pow(mach, 2) - 1)) / (gamma + 1))) - atan(sqrt(pow(mach, 2) - 1));
	return angle * (180 / M_PI);
}

double FlowAnalyzer::getMachFromPM(double M, double pm, double gamma){
    double nu = pm * (M_PI / 180);
    while (M <= 15) {
        double x = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt(((gamma - 1) * (pow(M, 2) - 1)) / (gamma + 1))) - atan(sqrt(pow(M, 2) - 1));
        double residual = nu - x;
        if (residual < 0.0001) {
            return M;
        }
        M += 0.001;
    }
    return 1.0;
}