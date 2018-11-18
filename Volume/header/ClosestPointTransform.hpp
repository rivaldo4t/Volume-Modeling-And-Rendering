#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class CPT : public Field<Vector>
	{
	private:
		SField f;
		VFIdentity I;
	public:
		CPT() {}
		CPT(SField _f, VFIdentity _I = VFIdentity()) : I(_I), f(_f) {}
		const FieldDataType eval(const Vector& p) const { return I.eval(p) - f->eval(p) * f->grad(p); }
	};
}
