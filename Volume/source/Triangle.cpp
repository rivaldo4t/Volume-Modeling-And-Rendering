#include "Triangle.hpp"

void Triangle::closestPointOnPlane(lux::Vector x, lux::Vector& p)
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

void Triangle::closestPointOnEdges(lux::Vector x, lux::Vector& p, double prevDist = std::numeric_limits<double>::max())
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

void Triangle::closestVertex(lux::Vector x, lux::Vector& p, double prevDist = std::numeric_limits<double>::max())
{
	double d0 = (x - p0).magnitude();
	double d1 = (x - p1).magnitude();
	double d2 = (x - p2).magnitude();

	if (d0 < prevDist || d1 < prevDist || d2 < prevDist)
	{
		if (d0 < d1)
		{
			if (d0 < d2)
				p = p0;
			else
				p = p2;
		}
		else if (d1 < d2)
			p = p1;
		else
			p = p2;
	}
}

lux::Vector Triangle::closestPoint(lux::Vector x)
{
	lux::Vector p(0.0, 0.0, 0.0);
	closestPointOnPlane(x, p);
	return p;
}