#pragma once
#include "PerlinNoise.hpp"

class FSPN
{
//private:
public:
	lux::PerlinNoise PN;
	lux::Vector translate;
	int octaves;
	float freq;
	float fJump;
	float roughness;
//public:
	FSPN()
	{ 
		PN = lux::PerlinNoise();
		octaves = 2;
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

	float eval(lux::Vector x);
};