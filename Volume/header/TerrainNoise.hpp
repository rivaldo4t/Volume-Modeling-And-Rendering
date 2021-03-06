#pragma once
#include "VectorField.hpp"
#include "Fspn.hpp"

namespace lux
{
	class TerrainNoise : public Field<double>
	{
	private:
		FSPN fspn;
		float up, down;
		float gammaUp, gammaDown;
	public:
		TerrainNoise() {}
		TerrainNoise(FSPN f = FSPN(), float u = 1.f, float d = 1.f, float gU = 1.f, float gD = 1.f) : 
		fspn(f), up(u), down(d), gammaUp(gU), gammaDown(gD) {}
		TerrainNoise(const NoiseParams& param, float u = 1.f, float d = 1.f, float gU = 1.f, float gD = 1.f) :
		up(u), down(d), gammaUp(gU), gammaDown(gD) { fspn = FSPN(param.octaves, param.freq, param.fJump, param.wedgeSpecific); }
		
		const FieldDataType eval(const Vector& p) const
		{
			float n = fspn.eval(p);
			if (n >= 0)
				return up * pow(n, gammaUp);
			else
				return -down * pow(-n, gammaDown);
		}
	};
}