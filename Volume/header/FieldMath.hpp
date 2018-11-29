#pragma once
#include "Field.hpp"

namespace lux
{
	class SFAdd : public Field<float>
	{
	private:
		SField s1, s2;
	public:
		SFAdd(SField a, SField b) : s1(a), s2(b) {}
		const FieldDataType eval(const Vector& p) const { return s1->eval(p) + s2->eval(p); }
	};

	class SFSubtract : public Field<float>
	{
	private:
		SField s1, s2;
	public:
		SFSubtract(SField a, SField b) : s1(a), s2(b) {}
		const FieldDataType eval(const Vector& p) const { return s1->eval(p) - s2->eval(p); }
	};

	class VFAdd : public Field<Vector>
	{
	private:
		VField v1, v2;
	public:
		VFAdd(VField a, VField b) : v1(a), v2(b) {}
		friend VField operator + (VField vf1, VField v2);
		const FieldDataType eval(const Vector& p) const { return v1->eval(p) + v2->eval(p); }
	};

	class VFSubtract : public Field<Vector>
	{
	private:
		VField v1, v2;
	public:
		VFSubtract(VField a, VField b) : v1(a), v2(b) {}
		friend VField operator - (VField vf1, VField vf2);
		const FieldDataType eval(const Vector& p) const { return v1->eval(p) - v2->eval(p); }
	};

	class VFMultiply : public Field<Vector>
	{
	private:
		VField v;
		SField s;
	public:
		VFMultiply(VField a, SField b) : v(a), s(b) {}
		friend VField operator * (VField vf, SField sf);
		const FieldDataType eval(const Vector& p) const { return v->eval(p) * s->eval(p); }
	};
}