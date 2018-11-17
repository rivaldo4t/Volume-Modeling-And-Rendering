#pragma once
#include "ScalarField.hpp"
#include "VectorField.hpp"

namespace lux
{
	class AdvectedField : public ScalarField
	{
		SField sf;
		VField vf;
		double deltaT = 0.2;
	public:
		AdvectedField(SField s, VField v) : sf(s), vf(v) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<AdvectedField>(*this);
		}
		double eval(const Vector& p) const
		{
			return sf->eval(p - vf->eval(p)*deltaT);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};
}