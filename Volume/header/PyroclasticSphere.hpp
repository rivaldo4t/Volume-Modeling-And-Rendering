#pragma once
#include "Vector.hpp"
#include "FractalSummedPerlinNoise.hpp"

class PyroClasticSphere
{
private:
	lux::Vector center;
	float radius;
public:
	double eval(lux::Vector x)
	{
		lux::Vector vec = x - center;
		float d = vec.magnitude();
		if (d < radius)
			return 1.0;
		
		float PNbound = radius * 0.1;
		if (d > radius + PNbound)
			return 0.0;
		
		lux::Vector surfacePoint = vec.unitvector();
		float scalingFact = 1.0;
		float gamma = 1.0;
		FSPN fspn; // params
		float r = abs(fspn.eval(radius * surfacePoint));
		r = pow(r, gamma);
		r *= scalingFact;
		if (d < radius + r)
			return 1.0;

		return 0.0;
	}
};