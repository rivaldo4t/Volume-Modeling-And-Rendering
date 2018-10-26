#include "StampedNoise.hpp"

void StampedNoise::computeNoise()
{
	gridData.resize(Nx * Ny * Nz, 0);

	static int _oct, _freq, _fjump, _fade;
	std::vector<int> oct = { 1, 2, 3, 4, 5 };
	std::vector<double> freq = { 0.5, 1, 2.5, 3, 4.3, 5.6, 10 };
	std::vector<double> fjump = { 1, 1.5, 2 };
	std::vector<double> fad = { 0.5, 1, 1.5, 2, 3 };

	fspn = FSPN(oct[_oct], freq[_freq], fjump[_fjump], 2.0);
	fade = fad[_fade];
	//fspn = FSPN(1, 10, 5, 5, lux::Vector(0, 0.5, 0.3).unitvector());
	//fspn = FSPN(5, 2, 2, 2);
	//fade = 1;

	std::cout << "-------------------\n";
	std::cout << "octaves:\t" << fspn.octaves << std::endl;
	std::cout << "freq:\t\t" << fspn.freq << std::endl;
	std::cout << "fjump:\t\t" << fspn.fJump << std::endl;
	std::cout << "roughness:\t" << fspn.roughness << std::endl;
	std::cout << "fade:\t\t" << fade << std::endl;
	std::cout << "-------------------\n";

	_fade++;
	if (_fade == fad.size())
		_fjump++;
	if (_fjump == fjump.size())
		_freq++;
	if (_freq == freq.size())
		_oct++;

	_fade %= fad.size();
	_fjump %= fjump.size();
	_freq %= freq.size();
	_oct %= oct.size();

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