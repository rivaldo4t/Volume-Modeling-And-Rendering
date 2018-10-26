#pragma once
#include "Vector.hpp"
#include <iostream>
#include <algorithm>
#include <memory>
#include <vector>

class Triangle
{
public:
	lux::Vector p0, p1, p2;
	lux::Vector e1, e2, e3;
	lux::Vector n;

	Triangle(lux::Vector _p0, lux::Vector _p1, lux::Vector _p2)
	{
		p0 = _p0;
		p1 = _p1;
		p2 = _p2;
		e1 = p1 - p0;
		e2 = p2 - p0;
		e3 = p2 - p1;
		n = (e1 ^ e2);
	}

	void closestPointOnPlane(lux::Vector x, lux::Vector& p);

	void closestPointOnEdges(lux::Vector x, lux::Vector& p, double prevDist);

	void closestVertex(lux::Vector x, lux::Vector& p, double prevDist);

	lux::Vector closestPoint(lux::Vector x);
};

typedef std::vector<std::shared_ptr<Triangle>> Triangles;