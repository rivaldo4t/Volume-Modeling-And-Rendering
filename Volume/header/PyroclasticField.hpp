#pragma once
#include "ScalarField.hpp"
#include "Fspn.hpp"

namespace lux
{
	class PyroclasticField : public ScalarField
	{
	private:
		SField f;
		FSPN fspn;
		float gamma;
		float scalingFact = 0.1;
	public:
		PyroclasticField(SField _f, FSPN _fspn = FSPN()) : f(_f), fspn(_fspn) {}
		PyroclasticField(SField _f, const NoiseParams& param) : f(_f)
		{ 
			fspn = FSPN(param.octaves, param.freq, param.fJump, 2);
			gamma = param.wedgeSpecific;
		}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<PyroclasticField>(*this);
		}
		double eval(const Vector& p) const
		{
			return f->eval(p) + scalingFact * pow(abs(fspn.eval(p.unitvector())), gamma);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};
}