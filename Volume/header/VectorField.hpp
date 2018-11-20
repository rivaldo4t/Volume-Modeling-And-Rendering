#pragma once
#include "Field.hpp"
#include "Fspn.hpp"

namespace lux
{
	class VFIdentity : public Field<Vector>
	{
	public:
		VFIdentity() {}
		const FieldDataType eval(const Vector& p) const { return p; }
	};

	class VFRandom : public Field<Vector>
	{
	private:
		FSPN f;
		Vector delta = Vector(0.1, 0.1, 0.1);
	public:
		VFRandom() {}
		const FieldDataType eval(const Vector& p) const { return Vector(f.eval(p), f.eval(p - delta), f.eval(p + delta)); }
	};

	class VFCharMap : public Field<Vector>
	{
	private:
		VField vf;
		double deltaT = 0.1;
	public:
		VFCharMap(VField v) : vf(v) {}
		const FieldDataType eval(const Vector& p) const { return p - vf->eval(p) * deltaT; }
	};
}