#pragma once
#include "ScalarGrid.hpp"
#include "Fspn.hpp"

namespace lux
{
	class StampedNoise : public ScalarGrid
	{
	private:
		lux::Vector p;
		float pScale;
		float fade;
		FSPN fspn;
	public:
		StampedNoise() : ScalarGrid() {}
		StampedNoise(lux::Vector _p, float pS, lux::Vector l, unsigned int nx, unsigned int ny, unsigned int nz, double d) :
			ScalarGrid(l, nx, ny, nz, d), p(_p), pScale(pS) {}
		void computeNoise(NoiseParams& param);
	};
}
