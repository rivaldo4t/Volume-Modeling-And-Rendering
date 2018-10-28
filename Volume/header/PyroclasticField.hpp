#pragma once
#include "ScalarField.hpp"
#include "Fspn.hpp"

namespace lux
{
	class pyroclasticField : public ScalarField
	{
	private:
		SField f;
		FSPN fspn;
	public:
		pyroclasticField(SField _f, FSPN _fspn = FSPN()) : f(_f), fspn(_fspn) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<pyroclasticField>(*this);
		}
		double eval(const Vector& p) const
		{
			float scalingFact = 1;
			float gamma = 1;
			return f->eval(p) + scalingFact * pow(abs(fspn.eval(p.unitvector())), gamma);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};
}