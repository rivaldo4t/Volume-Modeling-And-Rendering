#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class NPT : public VectorField
	{
	private:
		VFIdentity I;
		SField f;
	public:
		NPT() {}
		NPT(VFIdentity _I, SField _f) : I(_I), f(_f) {}
		virtual std::unique_ptr<VectorField> clone() const override
		{
			return std::make_unique<NPT>(*this);
		}
		Vector eval(const Vector& p) const
		{
			Vector gradF = f->grad(p);
			return I.eval(p) - f->eval(p) * gradF / gradF.magnitudeSquared();
		}
		Matrix grad(const Vector& p) const
		{
			return Matrix();
		}
		Vector convergeEval(const Vector& p) const
		{
			float threshold = 0.1f;
			Vector x = eval(p), xNext;
			while (true)
			{
				xNext = eval(x);
				if ((xNext - x).magnitude < threshold)
					break;
				x = xNext;
			}
			return x;
		}
	};
}
