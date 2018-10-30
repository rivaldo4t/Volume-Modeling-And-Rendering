#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class WarpField : public ScalarField
	{
	private:
		SField S;
		VField V;
	public:
		WarpField() {}
		WarpField(SField s, VField v) : S(s), V(v) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<WarpField>(*this);
		}
		double eval(const Vector& p) const
		{
			return S->eval(V->eval(p));
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};
}