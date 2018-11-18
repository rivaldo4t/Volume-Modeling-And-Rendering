#pragma once
#include "ScalarField.hpp"
#include "VectorField.hpp"

namespace lux
{
	class AdvectedField : public Field<float>
	{
		SField sf;
		VField vf;
		double deltaT = 0.08;
	public:
		AdvectedField(SField s, VField v) : sf(s), vf(v) {}
		const FieldDataType eval(const Vector& p) const { return sf->eval(p - vf->eval(p)*deltaT); }
	};
}