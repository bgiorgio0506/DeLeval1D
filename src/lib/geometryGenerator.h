
#pragma once
#include <cmath>
#include <vector>
#include "preliminaryCalculator.h"
#include "point.h"
#include "../thermo/flowEq.h"

struct GeometryParameters {
	int nLines; // Number of lines to generate
	double nozzleAreaRatio; // Area ratio for De Laval nozzle
	double throatRadious; //radius
};

typedef std::vector<cPoint> geometryLine;

class GeometryGenerator {
private:
	int nLines; // Number of lines to generate
	int nPoints; // Number of points per line
	CombustionReactor initialConditions; // Initial conditions for the combustion reactor
	preliminaryResults resolved1DimEngine;//resolved value from 1Dim analysis
	geometryLine geometryData;
	GeometryUtils geometryUtils;
	FlowAnalyzer flowAnalyzer;
	void initializeGeometryData(int nLines);
	tempPoint findNext(double theta, cPoint point, cPoint prev_point);
public:
	GeometryGenerator() = default;
	GeometryGenerator(preliminaryResults preliminaryRes);
	//TODO: data structure to hold the generated geometry data
	void generateGeometry(const GeometryParameters& params);
	void clearGeometry();
	//Get geomtry
	geometryLine getGeometry(PointType type = POINT_TYPE_UNKNOWN);

};

class GeometryUtils {
public:
	static std::vector<double> interpolateAngle(double theta, int n);
	static PlainPoint returnPlainPointIntersection(double x1, double y1, double tht1, double x2, double y2, double tht2);
};