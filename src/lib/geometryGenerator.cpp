
#define _USE_MATH_DEFINES

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <vector>
#include "geometryGenerator.h"
#include "point.h"

#define DEG2RAD(angle_deg) ((angle_deg) * M_PI / 180.0)
#define RAD2DEG(angle_rad) ((angle_rad) * 180.0 / M_PI)

GeometryGenerator::GeometryGenerator(preliminaryResults preliminaryRes)
{
	// Initialize the geometry generator with the specified type and initial conditions
	// This could include setting up any necessary parameters or data structures
	this->initialConditions = preliminaryRes.intialConditions;
	this->resolved1DimEngine = preliminaryRes;
}

void GeometryGenerator::generateGeometry(const GeometryParameters &params)
{
	// Clear previous geometry data
	geometryData.clear();
	// Determines max wall angles and intermidiates
	const double max_wall_angle = this->flowAnalyzer.getPmAngle(this->resolved1DimEngine.exitMachNumber, this->initialConditions.gamma) / 2;
	std::vector<double> wallAngles = this->geometryUtils.interpolateAngle(max_wall_angle, nLines);
	// init geometry vector
	this->initializeGeometryData(params.nLines);
	// methods of charactestics
	// step 0
	const double throatRadius = params.throatRadious;
	cPoint throatPoint = cPoint(0.0, throatRadius, 0.19, 0);
	throatPoint.setFlowProps(1.026, 0.19, 0);
	throatPoint.setPointType(POINT_TYPE_WALL);
	double throatSlope = -75.23; // deg

	double x0_loc = throatRadius * tan(DEG2RAD(90 + throatSlope));
	double mach = this->flowAnalyzer.getMachFromPM(throatPoint.getProperties().mach, throatPoint.getProperties().pm_angle + throatPoint.mu, this->initialConditions.gamma);
	double mu = RAD2DEG(asin(1 / mach));
	this->geometryData[0].setPoint(x0_loc, 0, mu);
	this->geometryData[0].setFlowProps(mach, throatPoint.mu + throatPoint.getProperties().pm_angle, 0);
	this->geometryData[0].setPointType(POINT_TYPE_CENTERLINE);

	for (int i = 1; i <= nLines; i++)
	{
		cPoint point = this->geometryData[i];
		cPoint prev_point = this->geometryData[i - 1];
		tempPoint tmpPoint;
		if (point.getProperties().type != POINT_TYPE_WALL)
		{
			point.setFlowProps(
				this->flowAnalyzer.getMachFromPM(prev_point.getProperties().mach, wallAngles[i], this->initialConditions.gamma),
				wallAngles[i],
				wallAngles[i]);
			const double mu = RAD2DEG(asin(1 / point.getProperties().mach));
			point.setMu(mu);
			tmpPoint = this->findNext(wallAngles[i], point, prev_point);
		}
		else if (point.getProperties().type == POINT_TYPE_WALL)
		{
			point.setFlowProps(prev_point.getProperties().mach, prev_point.getProperties().pm_angle, prev_point.getProperties().flow_angle);
			point.setMu(prev_point.getPoint()[2]);
			tmpPoint = this->findNext(max_wall_angle, point, prev_point);
		}

		PlainPoint coords = this->geometryUtils.returnPlainPointIntersection(throatPoint.getPoint()[0], throatPoint.getPoint()[1], tmpPoint.rightDerivative, prev_point.getPoint()[0], throatPoint.getPoint()[1], tmpPoint.leftDerivative);
		point.setXY(coords.x, coords.y);
	}

	int j = 0;
	int k = 1;
	for (int i = 0; i <= nPoints; i++)
	{
		cPoint point = this->geometryData[i];
		cPoint prev_point = this->geometryData[i - 1];
		tempPoint tmpPoint;
		switch (point.getProperties().type)
		{
		case (POINT_TYPE_WALL):
		{
			point.setFlowProps(prev_point.getProperties().mach, prev_point.getProperties().pm_angle, prev_point.getProperties().flow_angle);
			point.setMu(prev_point.mu);
			tmpPoint = this->findNext(prev_point.getProperties().flow_angle, point, prev_point);
			PlainPoint coords = this->geometryUtils.returnPlainPointIntersection(this->geometryData[i - (nLines - j)].x, this->geometryData[i - (nLines - j)].y, tmpPoint.rightDerivative, prev_point.x, prev_point.y, tmpPoint.leftDerivative);
			point.setXY(coords.x, coords.y);
			k++;
			break;
		}

		case (POINT_TYPE_CENTERLINE):
		{
			double pm_angle = this->geometryData[i - (nLines - j)].getProperties().flow_angle + this->geometryData[i - (nLines - j)].getProperties().pm_angle;
			point.setFlowProps(
				this->flowAnalyzer.getMachFromPM(this->geometryData[i - (nLines - j)].getProperties().mach, pm_angle, this->initialConditions.gamma),
				pm_angle,
				0);
			point.setMu(RAD2DEG(1 / point.getProperties().mach));
			tmpPoint = {
				0.0, // not needed for next steps!!
				0.0,
				0.0,
				0.0,
				((this->geometryData[i - (nLines - j)].getProperties().flow_angle - this->geometryData[i - (nLines - j)].mu + point.getProperties().flow_angle - point.mu) / 2)};
			PlainPoint coords = this->geometryUtils.returnPlainPointIntersection(this->geometryData[i - (nLines - j)].x, this->geometryData[i - (nLines - j)].y, tmpPoint.rightDerivative, this->geometryData[i - (nLines - j + 1)].x, this->geometryData[i - (nLines - j + 1)].y, tmpPoint.leftDerivative);
			point.setXY(coords.x, coords.y);
			break;
		}
		// POINT UNKWON
		default:
		{
			double pm_angle = ((this->geometryData[i - (nLines - j)].getProperties().flow_angle + this->geometryData[i - (nLines - j)].mu) - (prev_point.getProperties().flow_angle + prev_point.mu)) / 2;
			point.setFlowProps(
				this->flowAnalyzer.getMachFromPM(this->geometryData[i - (nLines - j)].getProperties().mach, pm_angle, this->initialConditions.gamma),
				pm_angle,
				wallAngles[k]);
			point.setMu(RAD2DEG(1 / point.getProperties().mach));
			double tht_r = (this->geometryData[i - (nLines - j)].getProperties().flow_angle - this->geometryData[i - (nLines - j)].mu + point.getProperties().flow_angle - point.mu) / 2;
			double tht_l = (prev_point.getProperties().flow_angle + prev_point.mu + point.getProperties().flow_angle + point.mu) / 2;
			tmpPoint = {0.0, 0.0, 0.0, tht_l, tht_r};
			PlainPoint coords = this->geometryUtils.returnPlainPointIntersection(this->geometryData[i - (nLines - j)].x, this->geometryData[i - (nLines - j)].y, tmpPoint.rightDerivative, prev_point.x, prev_point.y, tmpPoint.leftDerivative);
			point.setXY(coords.x, coords.y);
			k = 1;
			j++;
			break;
		}
		}
	}
}

tempPoint GeometryGenerator::findNext(double theta, cPoint point, cPoint prev_point)
{
	tempPoint result;
	switch (point.getProperties().type)
	{
	case (POINT_TYPE_UNKNOWN || POINT_TYPE_CENTERLINE): {
		// find theta_ax = nu_ax and the Mach angle
		double pm_ax = (2 * (theta)+(prev_point.getProperties().flow_angle - prev_point.getProperties().pm_angle)) / 2;
		double M_ax = this->flowAnalyzer.getMachFromPM(1.001, pm_ax, this->initialConditions.gamma);
		double mu_ax = RAD2DEG(asin(1 / M_ax));
		// find left and right derivatives
		double tht_l = (prev_point.getProperties().flow_angle + prev_point.mu + point.getProperties().flow_angle + point.mu) / 2;
		double tht_r = (pm_ax - mu_ax + point.getProperties().flow_angle - point.mu) / 2;
		result = { pm_ax, M_ax, mu_ax, tht_l, tht_r };
		break;
	}
	case (POINT_TYPE_WALL): {
		float tht_r = theta;
		float tht_l = (prev_point.getProperties().flow_angle + prev_point.mu + point.getProperties().flow_angle + point.mu) / 2;
		result = { 0.0, 0.0, 0.0, tht_l, tht_r };
		break;
	}
	}
	return result;
}

void GeometryGenerator::initializeGeometryData(int nLines)
{
	// Initialize the geometry data structure with the specified number of lines and points
	nLines = nLines;
	// init characteristics method
	nPoints = nLines + nLines * (nLines + 1) / 2;
	if (nPoints <= 0)
	{
		throw std::invalid_argument("Number of points must be greater than zero.");
	}
	int j = 1 + nLines;
	int k = 0;
	for (int i = 0; i <= nPoints; i++)
	{
		cPoint point(0.0, 0.0, 0.0, cPointProperties{i, 0.0, 0.0, 0.0, POINT_TYPE_UNKNOWN});
		if (i == j + k)
		{
			point.setPointType(POINT_TYPE_WALL);
			k = k + 1;
			j = j + nLines - k;
		}
		this->geometryData.push_back(point);
	}

	for (int i = 0; i < nPoints; i++)
	{
		if (this->geometryData[i].getProperties().index == 1)
			this->geometryData[i].setPointType(POINT_TYPE_CENTERLINE);
		else if (i > 1 && this->geometryData[i - 1].getProperties().type == POINT_TYPE_WALL)
			this->geometryData[i].setPointType(POINT_TYPE_CENTERLINE);
	}
}

geometryLine GeometryGenerator::getGeometry(PointType type)
{
	geometryLine data;
	switch (type)
	{
	case POINT_TYPE_WALL:
		for (int i = 0; i < geometryData.size(); i++)
		{
			if (geometryData[i].getProperties().type == type)
				data.push_back(geometryData[i]);
		}
		break;
	case POINT_TYPE_CENTERLINE:
		for (int i = 0; i < geometryData.size(); i++)
		{
			if (geometryData[i].getProperties().type == type)
				data.push_back(geometryData[i]);
		}
		break;
	default:
		data = geometryData;
		break;
	}
	return data;
}

void GeometryGenerator::clearGeometry()
{
	this->geometryData.clear();
}

std::vector<double> GeometryUtils::interpolateAngle(double theta, int n)
{
	std::vector<double> result = {0};
	double temp = theta - (1 + fmod(theta, 1));
	double delta = temp / (n - 2);
	for (int i = 0; i < n - 2; i++)
	{
		result.push_back(result[i] + delta);
	}
	result.push_back(theta);
	return result;
}
/*
 * Note to self refactor this
 */
PlainPoint GeometryUtils::returnPlainPointIntersection(double x1, double y1, double tht1, double x2, double y2, double tht2)
{
	tht1 = DEG2RAD(tht1);
	tht2 = DEG2RAD(tht2);
	double x = (x1 * tan(tht1) - x2 * tan(tht2) + y2 - y1) / (tan(tht1) - tan(tht2));
	double y = (tan(tht1) * tan(tht2) * (x1 - x2) + tan(tht1) * y1 - tan(tht1) * y1) / (tan(tht1) - tan(tht2));
	return {x, y};
}
