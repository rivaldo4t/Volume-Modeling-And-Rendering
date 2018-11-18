#include "StampedNoise.hpp"
using namespace lux;

void StampedNoise::computeNoise(NoiseParams& param)
{
	gridData.resize(Nx * Ny * Nz, 0);

	fspn = FSPN(param.octaves, param.freq, param.fJump, 2);
	fade = param.wedgeSpecific;
	
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

					float noiseVal = fspn.eval(xijk);
					noiseVal *= fadeFact;

					int index = getIndex(i, j, k);
					if (noiseVal > gridData[index])
						gridData[index] = noiseVal;
				}
			}
		}
	}

	std::cout << "Noise Stamped\n";
}