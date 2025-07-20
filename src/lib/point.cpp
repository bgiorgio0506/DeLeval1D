#include "point.h"

PointCoords Point::getPoint() {
	return PointCoords{ x, y, mu };
}

void Point::setPoint(double x, double y, double mu) {
	this->x = x;
	this->y = y;
	this->mu = mu;
}

PointCoords cPoint::getPoint() {
	return Point::getPoint();
}

void cPoint::setPoint(double x, double y, double mu) {
	Point::setPoint(x, y, mu);
}

void cPoint::setFlowProps(double mach, double pm_angle, double flow_angle) {
	properties.mach = mach;
	properties.pm_angle = pm_angle;
	properties.flow_angle = flow_angle;
}

void cPoint::setPointType(PointType type) {
	properties.type = type;
}

void cPoint::setMu(double mu) {
	this->mu = mu;
}

void cPoint::setXY(double x, double y) {
	this->x = x;
	this->y = y;
}
