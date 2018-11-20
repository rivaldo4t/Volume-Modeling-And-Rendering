#pragma once
#include "Field.hpp"

namespace lux 
{
	class SFSphere : public Field<float>
	{
		Vector center;
		double radius;
	public:
		SFSphere() {}
		SFSphere(Vector c, double r) : center(c), radius(r) {}
		const FieldDataType eval(const Vector& p) const { return radius - (p - center).magnitude(); }
		const FieldGradType grad(const Vector& p) const { return -(p - center) / (p-center).magnitude(); }
	};

	class SFPlane : public Field<float>
	{
		Vector anchor;
		Vector normal;
	public:
		SFPlane() {}
		SFPlane(Vector a, Vector n) : anchor(a), normal(n.unitvector()) {}
		const FieldDataType eval(const Vector& p) const { return (p - anchor)*normal; }
		const FieldGradType grad(const Vector& p) const { return eval(p) > 0 ? -normal : normal; }
	};

	class SFTorus : public Field<float>
	{
		double majorRadius, minorRadius;
		Vector center;
		Vector normal;
	public:
		SFTorus() {}
		SFTorus(double r1, double r2, Vector c, Vector n) : majorRadius(r1), minorRadius(r2), center(c), normal(n.unitvector()) {}
		const FieldDataType eval(const Vector& p) const
		{
			Vector perp = p - (p*normal)*normal;
			
			double fact1 = (perp - center).magnitude();
			fact1 *= 2 * majorRadius;
			fact1 *= fact1;
			
			double fact2 = (p - center).magnitude();
			fact2 = fact2*fact2 + majorRadius*majorRadius - minorRadius*minorRadius;
			fact2 *= fact2;
			
			return fact1 - fact2;
		}
	};

	class SFCone : public Field<float>
	{
		Vector center;
		Vector normal;
		double angle;
		double height;
	public:
		SFCone() {}
		SFCone(Vector c, Vector n, double a, double h) : center(c), normal(n.unitvector()), angle(a * M_PI / 180.0), height(h) {}
		const FieldDataType eval(const Vector& p) const
		{
			double temp = (p - center) * normal;
			if (p == center)
				return 0;
			else if (temp > height)
				return height - temp;
			else if (temp < 0)
				return temp;
			else
				return angle - acos(temp / (p - center).magnitude());

		}
	};

	class SFBox : public Field<float>
	{
		double boxLength;
	public:
		SFBox() {}
		SFBox(double l) : boxLength(l) {}
		const FieldDataType eval(const Vector& p) const
		{
			double x = p[0], y = p[1], z = p[2];
			double k = 6;
			return boxLength - pow(x, k) - pow(y, k) - pow(z, k);
		}
	};

	class SFIcosahedron : public Field<float>
	{
		Vector center;
	public:
		SFIcosahedron() {}
		SFIcosahedron(Vector c) : center(c) {}
		const FieldDataType eval(const Vector& p) const
		{
			double fact1 = 1.8 * M_PI;
			double fact2 = 1.61803399; // golden ration
			double dist = (p - center).magnitude();
			double x = p[0], y = p[1], z = p[2];
			if (dist > fact1)
				return -fact1;
			else
				return	cos(x + fact2 * y) + cos(x - fact2 * y) +
						cos(y + fact2 * z) + cos(y - fact2 * z) +
						cos(z + fact2 * x) + cos(z - fact2 * x) - 2.0;
		}
	};

	class SFEllipse : public Field<float>
	{
		double majorRadius, minorRadius;
		Vector center;
		Vector normal;
	public:
		SFEllipse() {}
		SFEllipse(double r1, double r2, Vector c, Vector n) : majorRadius(r1), minorRadius(r2), center(c), normal(n.unitvector()) {}
		const FieldDataType eval(const Vector& p) const
		{
			double fact1 = p*normal;
			Vector perp = p - fact1*normal;
			double fact2 = (perp - center).magnitude();
			fact1 *= fact1 / (majorRadius * majorRadius);
			fact2 *= fact2 / (minorRadius * minorRadius);

			return 1 - fact1 - fact2;
		}
	};

	class SFStienerPatch : public Field<float>
	{
	public:
		SFStienerPatch() {}
		const FieldDataType eval(const Vector& p) const
		{
			double x = p[0], y = p[1], z = p[2];
			return -(x*x*y*y + x*x*z*z + y*y*z*z - x*y*z);
		}
	};

	class SFCylinder : public Field<float>
	{
		Vector center;
		Vector normal;
		double height;
		double radius;
	public:
		SFCylinder() {}
		SFCylinder(Vector c, Vector n, double h, double r) : center(c), normal(n.unitvector()), height(h), radius(r) {}
		double getHeight()
		{ return height; }
		const FieldDataType eval(const Vector& p) const
		{
			Vector perp = p - (p*normal)*normal;
			return radius - (perp - center).magnitude();
		}
	};

	class SFIntersect : public Field<float>
	{
		SField _f, _g;
	public:
		SFIntersect() {}
		SFIntersect(const SField& f, const SField& g) : _f(f), _g(g){}
		const FieldDataType eval(const Vector& p) const { return std::min(_f->eval(p), _g->eval(p)); }
	};

	class SFUnion : public Field<float>
	{
		SField _f, _g;
	public:
		SFUnion() {}
		SFUnion(const SField& f, const SField& g) : _f(f), _g(g) {}
		const FieldDataType eval(const Vector& p) const { return std::max(_f->eval(p), _g->eval(p)); }
	};

	class SFCutout : public Field<float>
	{
		SField _f, _g;
	public:
		SFCutout() {}
		SFCutout(const SField& f, const SField& g) : _f(f), _g(g) {}
		const FieldDataType eval(const Vector& p) const { return std::min(_f->eval(p), -(_g->eval(p))); }
	};

	class SFShell : public Field<float>
	{
		SField _f;
		double thickness;
	public:
		SFShell() {}
		SFShell(const SField& f, double t) : _f(f), thickness(t) {}
		const FieldDataType eval(const Vector& p) const
		{
			double eval = _f->eval(p);
			return std::min(eval + thickness * 0.5, -(eval - thickness * 0.5));
		}
	};

	class SFTranslate : public Field<float>
	{
		SField _f;
		Vector _x;
	public:
		SFTranslate() {}
		SFTranslate(SField f, Vector x) : _f(f), _x(x) {}
		const FieldDataType eval(const Vector& p) const { return _f->eval(p - _x); }
	};

	class SFScale : public Field<float>
	{
		SField _f;
		double _s;
	public:
		SFScale() {}
		SFScale(SField f, double s) : _f(f), _s(s) {}
		const FieldDataType eval(const Vector& p) const { return _f->eval(p / _s); }
	};

	class SFRotate : public Field<float>
	{
		SField _f;
		Vector axis;
		double angle;
	public:
		SFRotate() {}
		SFRotate(SField f, Vector ax, double an) : _f(f), axis(ax.unitvector()), angle(-an * M_PI / 180.0) {}
		const FieldDataType eval(const Vector& p) const
		{
			double A = cos(angle);
			double BB = (axis * p) * (1 - A);
			double C = sin(angle);
			Vector rotated_p = (A * p) + (BB * axis) + (C * (p ^ axis));
			return _f->eval(rotated_p);
		}
	};

	class SFMask : public Field<float>
	{
		SField _f;
	public:
		SFMask() {}
		SFMask(SField f) : _f(f) {}
		const FieldDataType eval(const Vector& p) const { return _f->eval(p) > 0.0 ? 1.0 : 0.0; }
	};

	class SFClamp : public Field<float>
	{
		SField _f;
		double low, high;
	public:
		SFClamp() {}
		SFClamp(SField f, double l, double h) : _f(f), low(l), high(h) {}
		const FieldDataType eval(const Vector& p) const
		{
			double eval = _f->eval(p);
			return eval < low ? high : eval > high ? high : eval;
		}
	};
}