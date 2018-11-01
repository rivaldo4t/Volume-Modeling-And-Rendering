#pragma once
#include "PerlinNoise.hpp"
#include <iostream>
#include <vector>

class NoiseParams
{
//private:
public:
	int octaves;
	float freq;
	float fJump;
	float wedgeSpecific; // fade, gamma, clump
	int o, fr, fj, w;
	bool randomVal;
	std::vector<int> octaves_values = { /*1, 2, 3,*/ 4, 5 };
	std::vector<double> freq_values = { /*0.5, 1,*/ 2.5, 3, 4.3, 5.6, 10 };
	std::vector<double> fJump_values = { /*1,*/ 1.5, 2 };
	std::vector<double> wedgeSpecific_values = { /*0.5, 1,*/ 1.5, 2, 3 };
//public:
	NoiseParams() { o = fr = fj = w = 0; }
	NoiseParams(bool randomInit)  : randomVal(randomInit)
	{
		int r = rand() % 10;
		w = r % wedgeSpecific_values.size();
		fj = r % fJump_values.size();
		fr = r % freq_values.size();
		o = r % octaves_values.size();
	}
	void updateParams();
	void setParams(int o, float fr, float fj, float ws);
};

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
	float eval(lux::Vector x) const;
};