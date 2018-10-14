#pragma once
#include "Grid.hpp"
#include "FractalSummedPerlinNoise.hpp"

class StampedNoise : public Grid
{
private:
	lux::Vector p;
	float pScale;
	float fade;
	FSPN fspn; // params
public:
	StampedNoise(lux::Vector _p, float pScale) : Grid() {} // p+- pscale, pscale, pscale, fill grid params
	void computeNoise()
	{
		gridData.resize(Nx * Ny * Nz, 0);

#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector xijk = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					lux::Vector vec = xijk - p;
					if (vec.magnitude() <= pScale)
					{
						float fadeFact = vec.magnitude() / pScale;
						fadeFact = 1 - fadeFact;
						fadeFact = pow(fadeFact, fade);
						
						float noiseVal = fspn.eval(xijk);
						noiseVal *= fadeFact;

						int index = getIndex(i, j, k);
						if (noiseVal > gridData[index])
							gridData[index] = noiseVal;
					}
				}
			}
		}
	}
};
