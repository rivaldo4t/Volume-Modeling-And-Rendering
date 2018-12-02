#pragma once
#include "ScalarField.hpp"
#include "Fspn.hpp"
#include "ClosestPointTransform.hpp"

namespace lux
{
	class PyroclasticField : public Field<float>
	{
	private:
		SField f;
		FSPN fspn;
		float gamma;
		float scalingFact;
		CPT cpt;
	public:
		PyroclasticField(SField _f, FSPN _fspn = FSPN()) : f(_f), fspn(_fspn) {}
		PyroclasticField(SField _f, const NoiseParams& param, float sc = 0.5) : f(_f), scalingFact(sc)
		{ 
			fspn = FSPN(param.octaves, param.freq, param.fJump, 4);
			gamma = param.wedgeSpecific;
			cpt = CPT(f);
		}
		const FieldDataType eval(const Vector& p) const
		{
			return f->eval(p) + scalingFact * pow(abs(fspn.eval(p.unitvector())), gamma);
			
			// using cpt
			return f->eval(p) + scalingFact * pow(abs(fspn.eval(cpt.eval(p))), gamma);
		}
	};
}