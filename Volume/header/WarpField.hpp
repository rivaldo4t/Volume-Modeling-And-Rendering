#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class WarpField : public Field<float>
	{
	private:
		SField f;
		CPT cptF;
		SField noiseField;
	public:
		WarpField() {}
		WarpField(SField s1, SField s2) : f(s1), noiseField(s2), cptF(CPT(f)) {}
		const FieldDataType eval(const Vector& p) const { return f->eval(p) + noiseField->eval(cptF.eval(p)); }
	};

	class WarpFieldVF : public Field<float>
	{
		SField f;
		VField v;
	public:
		WarpFieldVF() {}
		WarpFieldVF(SField s1, VField s2) : f(s1), v(s2) {}
		const FieldDataType eval(const Vector& p) const { return f->eval(v->eval(p)); }
	};
}