#pragma once
#include "ScalarField.hpp"
#include "VectorField.hpp"

namespace lux
{
	class AdvectedField : public Field<float>
	{
		SField sf;
		VField vf;
		double deltaT;
	public:
		AdvectedField(SField s, VField v, double d = 0.1) : sf(s), vf(v), deltaT(d) {}
		const FieldDataType eval(const Vector& p) const {
			/*if (p.Y() < 0.4)
				return sf->eval(p);
			else*/
				return sf->eval(p - vf->eval(p) * deltaT); 
		}
	};

	class AdvectedFieldCM : public Field<float>
	{
		SField sf;
		VField cm;
		double deltaT;
	public:
		AdvectedFieldCM(SField s, VField v, double d = 0.1) : sf(s), cm(v), deltaT(d) {}
		const FieldDataType eval(const Vector& p) const { return sf->eval(cm->eval(p)); }
	};
}