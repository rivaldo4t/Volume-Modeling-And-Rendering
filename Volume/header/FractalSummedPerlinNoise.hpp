#pragma once
#include "PerlinNoise.hpp"

class FSPN
{
private:
	lux::PerlinNoise PN;
	lux::Vector translate;
	int octaves;
	float freq;
	float fJump;
	float roughness;
public:
	FSPN()
	{ 
		PN = lux::PerlinNoise();
		octaves = 1;
		freq = 10;
		fJump = 2;
		roughness = 2;
	}

	FSPN(int o, float fq, float fj, float r, lux::Vector t = lux::Vector())
	{
		PN = lux::PerlinNoise();
		octaves = o;
		freq = fq;
		fJump = fj;
		roughness = r;
	}

	float eval(lux::Vector x)
	{
		float coeff = roughness == 1.0 ? (1.0 / float(octaves)) : (1 - roughness) / (1 - pow(roughness, octaves));
		float sumPN = 0;
		for (int i = 0; i < octaves; ++i)
		{
			lux::Vector evalAt = (x - translate) * pow(fJump, i) * freq;
			float pn = PN.eval(evalAt);
			sumPN += pow(roughness, i) * pn;
		}
		return coeff * sumPN;
	}
};