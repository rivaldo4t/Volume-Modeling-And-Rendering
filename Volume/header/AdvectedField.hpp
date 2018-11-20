#pragma once
#include "ScalarField.hpp"
#include "VectorField.hpp"

namespace lux
{
	class AdvectedField : public Field<float>
	{
		SField sf;
		VField vf;
		double deltaT = 0.1;
	public:
		AdvectedField(SField s, VField v) : sf(s), vf(v) {}
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
		double deltaT = 0.1;
	public:
		AdvectedFieldCM(SField s, VField v) : sf(s), cm(v) {}
		const FieldDataType eval(const Vector& p) const { return sf->eval(cm->eval(p)); }
	};
}