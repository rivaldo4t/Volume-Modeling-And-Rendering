#pragma once
#include "VectorField.hpp"
#include "Fspn.hpp"

namespace lux
{
	class TerrainNoise : public VectorField
	{
	private:
		FSPN fspn;
		Vector up, down;
		float gammaUp, gammaDown;
	public:
		TerrainNoise() {}
		TerrainNoise(FSPN f = FSPN(), Vector u = Vector(), Vector d = Vector(), float gU = 1.f, float gD = 1.f) : 
		fspn(f), up(u), down(d), gammaUp(gU), gammaDown(gD) {}
		TerrainNoise(const NoiseParams& param, Vector u = Vector(), Vector d = Vector(), float gU = 1.f, float gD = 1.f) :
		up(u), down(d), gammaUp(gU), gammaDown(gD)
		{
			fspn = FSPN(param.octaves, param.freq, param.fJump, 2);
		}
		virtual std::unique_ptr<VectorField> clone() const override
		{
			return std::make_unique<TerrainNoise>(*this);
		}
		Vector eval(const Vector& p) const
		{
			float n = fspn.eval(p);
			if (n >= 0)
				return up * pow(n, gammaUp);
			else
				return -down * pow(-n, gammaDown);
		}
		Matrix grad(const Vector& p) const
		{
			return Matrix();
		}
	};
}