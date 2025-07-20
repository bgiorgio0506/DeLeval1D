#ifndef ENGINE_H
#define ENGINE_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cmath>
#include <string>
#include <vector>
#include "geometryGenerator.h"
#include "csvExporter.h"
class engine
{
public:
	CombustionReactor intialConditions;
	double throatArea;
	double throatR;
	int nPoints;   // Number of points to use for geometry calculations
	int nSegments; // Number of segments to use for geometry calculations

	// algo calculated values
	GeometryGenerator geometryGenerator;
	preliminaryResults resolved;
	DataExporter exporter;
	double exitSoundSpeed;
	double exitTemperature;
	// bools
	bool isNozzleResolved = false;

	engine(CombustionReactor initialConditions, int nSegments) : intialConditions(initialConditions), throatArea(0), throatR(0), nSegments(nSegments)
	{
		this->throatR = sqrt(resolved.throatArea / M_PI); // Calculate the radius from the area
		this->throatArea = resolved.throatArea;
		geometryGenerator = GeometryGenerator(resolved);
	}
#pragma region dimensioning
	void preliminaryAnalysis();
	void generateGeometry();
	/*
	 * TODO: thermodynamic analysis
	 */
	void exportResults(const std::string &fileName);
#pragma region simulaion and performace
	/*
	 * TODO: methods for perfomance testing and simulations
	 */

#pragma region utils
	preliminaryResults getPreliminaryResults();
	geometryLine getNozzleContours();
	void reset();
};
#endif
