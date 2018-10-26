#pragma once
#include "Grid.hpp"
#include "Fspn.hpp"

class StampedNoise : public Grid
{
private:
	lux::Vector p;
	float pScale;
	float fade;
	FSPN fspn;
public:
	StampedNoise() : Grid () {}
	StampedNoise(lux::Vector _p, float pS,lux::Vector l, int nx, unsigned int ny, unsigned int nz, double d) :
		Grid(l, nx, ny, nz, d), p(_p), pScale(pS) {}
	void computeNoise();
};
