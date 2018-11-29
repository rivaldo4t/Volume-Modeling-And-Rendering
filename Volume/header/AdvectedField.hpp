#pragma once
#include "ScalarField.hpp"
#include "VectorField.hpp"

namespace lux
{
	// create multiadvectedfields
	class AdvectedSField : public Field<float>
	{
		SField sf;
		VField adv;
		double deltaT;
	public:
		AdvectedSField(SField s, VField v, double d = 0.1) : sf(s), adv(v), deltaT(d) {}
		const FieldDataType eval(const Vector& p) const {
			/*if (p.Y() < 0.4)
				return sf->eval(p);
			else*/
				return sf->eval(p - adv->eval(p) * deltaT);
		}
	};

	class AdvectedVField : public Field<Vector>
	{
		VField vf;
		VField adv;
		double deltaT;
	public:
		AdvectedVField(VField v1, VField v2, double d = 0.1) : vf(v1), adv(v2), deltaT(d) {}
		const FieldDataType eval(const Vector& p) const { return vf->eval(p - adv->eval(p) * deltaT); }
	};
}