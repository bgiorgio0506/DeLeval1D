
#ifndef POINT_H
#define POINT_H
#pragma once
#include <vector>

enum PointType
{
	POINT_TYPE_UNKNOWN = 0,
	POINT_TYPE_INLET,
	POINT_TYPE_OUTLET,
	POINT_TYPE_CENTERLINE,
	POINT_TYPE_WALL
};

struct cPointProperties
{
	int index;		   // Index of the point
	double flow_angle; // Flow angle in radians
	double pm_angle;
	double mach;
	PointType type;
};

struct PlainPoint
{
	double x;
	double y;
};

struct tempPoint
{
	double theta;
	double machAngle;
	double mu;
	double leftDerivative;
	double rightDerivative;
};

typedef std::vector<double> PointCoords; // Vector to hold point coordinates

class Point
{
public:
	double x;  // X coordinate
	double y;  // Y coordinate
	double mu; // angle in radians

	// Default constructor
	Point() : x(0.0), y(0.0), mu(0.0) {}

	// Parameterized constructor
	Point(double x, double y, double mu) : x(x), y(y), mu(mu) {}

	// Methods get & set
	PointCoords getPoint();

	void setPoint(double x, double y, double mu);
};

class cPoint : public Point
{
public:
	cPointProperties properties; // Properties of the point
	// Costruttore completo: inizializza tutto
	cPoint(double x, double y, double mu, cPointProperties props)
		: Point(x, y, mu), properties(props) {}

	cPoint(double x, double y, double mu, int index)
		: Point(x, y, mu)
	{
		properties.index = index;			  // Set index
		properties.type = POINT_TYPE_UNKNOWN; // Default type
	}

	cPoint() : Point(), properties()
	{
		// Default constructor initializes point and properties
		properties.index = -1;				  // Default index
		properties.flow_angle = 0.0;		  // Default flow angle
		properties.pm_angle = 0.0;			  // Default pm angle
		properties.mach = 0.0;				  // Default mach number
		properties.type = POINT_TYPE_UNKNOWN; // Default point type
	}

	// Metodo per ottenere le coordinate del punto
	PointCoords getPoint();

	// Metodo per impostare le coordinate del punto
	void setPoint(double x, double y, double mu);
	void setMu(double mu);
	void setXY(double x, double y);

	// Metodo per ottenere le propriet� del punto
	cPointProperties getProperties() const
	{
		return properties;
	}
	// Metodo per impostare le propriet� del punto
	void setPointType(PointType type);
	void setFlowProps(double mach, double pm_angle, double flow_angle);
};

#endif // !POINT_H
