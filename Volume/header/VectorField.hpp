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
		Vector delta;
	public:
		VFRandom(FSPN _f = FSPN(), Vector d = Vector(0.1, 0.1, 0.1)) : f(_f), delta(d) {}
		const FieldDataType eval(const Vector& p) const { return Vector(f.eval(p), f.eval(p - delta), f.eval(p + delta)); }
	};

	class VFCharMap : public Field<Vector>
	{
	private:
		VField vf;
		double deltaT;
	public:
		VFCharMap(VField v, double d = 0.1) : vf(v), deltaT(d) {}
		const FieldDataType eval(const Vector& p) const { return p - vf->eval(p) * deltaT; }
	};
}