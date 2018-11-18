#pragma once
#include <iostream>
#include "VectorField.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class WarpField : public Field<double>
	{
	private:
		SField noiseField;
		SField f;
		CPT cptF;
	public:
		WarpField() {}
		WarpField(SField s1, SField s2) : f(s1), noiseField(s2), cptF(CPT(f)) {}
		const FieldDataType eval(const Vector& p) const { return f->eval(p) + noiseField->eval(cptF.eval(p)); }
	};
}