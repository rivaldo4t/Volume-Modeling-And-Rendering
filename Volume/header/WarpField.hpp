#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class WarpField : public ScalarField
	{
	private:
		SField noiseField;
		SField f;
		CPT cptF;
	public:
		WarpField() {}
		WarpField(SField s1, SField s2) : f(s1), noiseField(s2), cptF(CPT(f)) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<WarpField>(*this);
		}
		double eval(const Vector& p) const
		{
			return f->eval(p) + noiseField->eval(cptF.eval(p));
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};
}