#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class NPT : public Field<Vector>
	{
	private:
		VFIdentity I;
		SField f;
	public:
		NPT() {}
		NPT(VFIdentity _I, SField _f) : I(_I), f(_f) {}
		const FieldDataType eval(const Vector& p) const
		{
			Vector gradF = f->grad(p);
			return I.eval(p) - f->eval(p) * gradF / gradF.magnitudeSquared();
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
