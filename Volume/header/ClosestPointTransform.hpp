#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class CPT : public VectorField
	{
	private:
		SField f;
		IdentityVectorField I;
	public:
		CPT() {}
		CPT(SField _f, IdentityVectorField _I = IdentityVectorField()) : I(_I), f(_f) {}
		virtual std::unique_ptr<VectorField> clone() const override
		{
			return std::make_unique<CPT>(*this);
		}
		Vector eval(const Vector& p) const
		{
			return I.eval(p) - f->eval(p) * f->grad(p);
		}
		Matrix grad(const Vector& p) const
		{
			return Matrix();
		}
	};
}