#pragma once
#include "Vector.hpp"
#include <iostream>
#include <algorithm>
#include <memory>

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

	void closestPointOnPlane(lux::Vector x, lux::Vector& p)
	{
		double n_mag_sq = n.magnitudeSquared();
		lux::Vector n_unit = n.unitvector();

		lux::Vector r = x - p0 - (n_unit * (n_unit * (x - p0)));

		double v = ((e1 ^ r) * (n)) / n_mag_sq;
		double u = ((e2 ^ r) * (-n)) / n_mag_sq;
		double d = std::numeric_limits<double>::max();

		if ((0 <= u) && (u <= 1) && 
			(0 <= v) && (v <= 1) && 
			(0 <= (u + v)) && ((u + v) <= 1))
		{
			p = p0 + (u * e1) + (v * e2);
			d = (x - p).magnitude();
		}
		closestPointOnEdges(x, p, d);
	}

	void closestPointOnEdges(lux::Vector x, lux::Vector& p, double prevDist = std::numeric_limits<double>::max())
	{
		double e1_mag_sq = e1.magnitudeSquared();
		double q1 = (e1 * (x - p0)) / e1_mag_sq;
		
		double e2_mag_sq = e2.magnitudeSquared();
		double q2 = (e2 * (x - p0)) / e2_mag_sq;
		
		double e3_mag_sq = e3.magnitudeSquared();
		double q3 = (e3 * (x - p1)) / e3_mag_sq;
		
		double d = prevDist;

		if (((0 <= q1) && (q1 <= 1)) || 
			((0 <= q2) && (q2 <= 1)) || 
			((0 <= q3) && (q3 <= 1)))
		{
			if ((0 <= q1) && (q1 <= 1))
			{
				double dist = (p0 + (q1 * e1) - x).magnitude();
				if (dist < d)
				{
					d = dist;
					p = p0 + q1 * e1;
				}
			}
			if ((0 <= q2) && (q2 <= 1))
			{
				double dist = (p0 + (q2 * e2) - x).magnitude();
				if (dist < d)
				{
					d = dist;
					p = p0 + q2 * e2;
				}
			}

			if ((0 <= q3) && (q3 <= 1))
			{
				double dist = (p1 + (q3 * e3) - x).magnitude();
				if (dist < d)
				{
					d = dist;
					p = p1 + q3 * e3;
				}
			}
		}
		closestVertex(x, p, d);
	}

	void closestVertex(lux::Vector x, lux::Vector& p, double prevDist = std::numeric_limits<double>::max())
	{
		double d0 = (x - p0).magnitude();
		double d1 = (x - p1).magnitude();
		double d2 = (x - p2).magnitude();

		if (d0 < prevDist || d1 < prevDist || d2 < prevDist)
		{
			double d = std::min(d0, std::min(d1, d2));
			if (d == d0)
				p = p0;
			else if (d == d1)
				p = p1;
			else if (d == d2)
				p = p2;
		}
	}

	lux::Vector closestPoint(lux::Vector x)
	{
		lux::Vector p(0.0, 0.0, 0.0);
		closestPointOnPlane(x, p);
		return p;
	}
};

typedef std::vector<std::shared_ptr<Triangle>> Triangles;