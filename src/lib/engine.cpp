#define _USE_MATH_DEFINES

#include <stdexcept>
#include "engine.h"
#include "preliminaryCalculator.h"

void engine::preliminaryAnalysis()
{
	PreliminaryAnalyser analyser(this->intialConditions);
	this->resolved = analyser.resolve();
	this->geometryGenerator.updateConditions(this->resolved);
	this->isNozzleResolved = true;
}

void engine::generateGeometry()
{
	if (this->isNozzleResolved)
	{
		this->geometryGenerator.generateGeometry({nSegments, this->resolved.exitAreaToThroatAreaRatio, sqrt(this->resolved.throatArea / M_PI)});
	}
	else
		std::invalid_argument("Nozzle not yet resolved");
}

void engine::exportResults(const std::string &fileName)
{
	// Geometry Results
	const string path = "./data/" + fileName;
	std::vector<string> collumsNames = {"index", "x", "y", "mu", "flow_angle", "mach", "pm_angle", "type"};
	geometryLine geo = this->geometryGenerator.getGeometry();
	exporter.exportToCSV(geo, path, collumsNames);
}

#pragma region utils
preliminaryResults engine::getPreliminaryResults()
{
	return this->resolved;
}

geometryLine engine::getNozzleContours()
{
	return this->geometryGenerator.getGeometry(POINT_TYPE_WALL);
}