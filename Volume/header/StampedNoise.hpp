#pragma once
#include "Grid.hpp"
#include "FractalSummedPerlinNoise.hpp"

class StampedNoise : public Grid
{
private:
	lux::Vector p;
	float pScale;
	float fade;
	FSPN fspn;
public:
	StampedNoise() : Grid () {}
	StampedNoise(lux::Vector _p, float pS, float f, FSPN& fractalNoise ,lux::Vector l, int nx, unsigned int ny, unsigned int nz, double d) : 
		Grid(l, nx, ny, nz, d), p(_p), pScale(pS), fade(f), fspn(fractalNoise) { fspn = FSPN(); }

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
					float dist = vec.magnitude();
					if (dist <= pScale)
					{
						float fadeFact = dist / pScale;
						fadeFact = 1 - fadeFact;
						fadeFact = pow(fadeFact, fade);
						
						float noiseVal = fspn.eval(xijk.unitvector()); //xijk or unit vector? pn gives zero if not unit vec
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
