#pragma once
#include "Vector.hpp"
#include "Color.hpp"
#include <iostream>
#include <algorithm>
#include <memory>

namespace lux 
{
	class ScalarField 
	{
	public:
		virtual double eval(const Vector& p) const = 0;
		virtual Vector grad(const Vector& p) const = 0;
		virtual std::unique_ptr<ScalarField> clone() const = 0;
	};
	typedef std::shared_ptr<lux::ScalarField> SField;

	class SFEmpty : public ScalarField
	{
	public:
		SFEmpty() {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFEmpty>(*this);
		}
		double eval(const Vector& p) const
		{
			// caution
			return 0.0;
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFSphere : public ScalarField
	{
		Vector center;
		double radius;
	public:
		SFSphere() {}
		SFSphere(const Vector c, double r) : center(c), radius(r) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFSphere>(*this);
		}
		double eval(const Vector& p) const
		{
			return radius - (p - center).magnitude();
		}
		Vector grad(const Vector& p) const
		{
			Vector g = -(p - center);
			g.normalize();
			return g;
		}
	};

	class SFPlane : public ScalarField
	{
		Vector anchor;
		Vector normal;
	public:
		SFPlane() {}
		SFPlane(Vector a, Vector n) : anchor(a), normal(n.unitvector()) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFPlane>(*this);
		}
		double eval(const Vector& p) const
		{
			return (p - anchor)*normal;
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFTorus : public ScalarField
	{
		double majorRadius, minorRadius;
		Vector center;
		Vector normal;
	public:
		SFTorus() {}
		SFTorus(double r1, double r2, Vector c, Vector n) : majorRadius(r1), minorRadius(r2), center(c), normal(n.unitvector()) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFTorus>(*this);
		}
		double eval(const Vector& p) const
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
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFCone : public ScalarField
	{
		Vector center;
		Vector normal;
		double angle;
		double height;
	public:
		SFCone() {}
		SFCone(Vector c, Vector n, double a, double h) : center(c), normal(n.unitvector()), angle(a * M_PI / 180.0), height(h) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFCone>(*this);
		}
		double eval(const Vector& p) const
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
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFBox : public ScalarField
	{
		double boxLength;
	public:
		SFBox() {}
		SFBox(double l) : boxLength(l) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFBox>(*this);
		}
		double eval(const Vector& p) const
		{
			double x = p[0], y = p[1], z = p[2];
			double k = 6;
			return boxLength - pow(x, k) - pow(y, k) - pow(z, k);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFIcosahedron : public ScalarField
	{
		Vector center;
	public:
		SFIcosahedron() {}
		SFIcosahedron(Vector c) : center(c) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFIcosahedron>(*this);
		}
		double eval(const Vector& p) const
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
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFEllipse : public ScalarField
	{
		double majorRadius, minorRadius;
		Vector center;
		Vector normal;
	public:
		SFEllipse() {}
		SFEllipse(double r1, double r2, Vector c, Vector n) : majorRadius(r1), minorRadius(r2), center(c), normal(n.unitvector()) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFEllipse>(*this);
		}
		double eval(const Vector& p) const
		{
			double fact1 = p*normal;
			Vector perp = p - fact1*normal;
			double fact2 = (perp - center).magnitude();
			fact1 *= fact1 / (majorRadius * majorRadius);
			fact2 *= fact2 / (minorRadius * minorRadius);

			return 1 - fact1 - fact2;
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFStienerPatch : public ScalarField
	{
	public:
		SFStienerPatch() {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFStienerPatch>(*this);
		}
		double eval(const Vector& p) const
		{
			double x = p[0], y = p[1], z = p[2];
			return -(x*x*y*y + x*x*z*z + y*y*z*z - x*y*z);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFCylinder : public ScalarField
	{
		Vector center;
		Vector normal;
		double height;
		double radius;
	public:
		SFCylinder() {}
		SFCylinder(Vector c, Vector n, double h, double r) : center(c), normal(n.unitvector()), height(h), radius(r) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFCylinder>(*this);
		}
		double eval(const Vector& p) const
		{
			Vector perp = p - (p*normal)*normal;
			return radius - (perp - center).magnitude();
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
		double getHeight()
		{
			return height;
		}
	};

	class SFIntersect : public ScalarField
	{
		SField _f, _g;
	public:
		SFIntersect() {}
		SFIntersect(const SField& f, const SField& g) : _f(f), _g(g){}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFIntersect>(*this);
		}
		double eval(const Vector& p) const
		{
			return std::min(_f->eval(p), _g->eval(p));
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFUnion : public ScalarField
	{
		SField _f, _g;
	public:
		SFUnion() {}
		SFUnion(const SField& f, const SField& g) : _f(f), _g(g) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFUnion>(*this);
		}
		double eval(const Vector& p) const
		{
			return std::max(_f->eval(p), _g->eval(p));
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFCutout : public ScalarField
	{
		SField _f, _g;
	public:
		SFCutout() {}
		SFCutout(const SField& f, const SField& g) : _f(f), _g(g) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFCutout>(*this);
		}
		double eval(const Vector& p) const
		{
			return std::min(_f->eval(p), -(_g->eval(p)));
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFShell : public ScalarField
	{
		SField _f;
		double thickness;
	public:
		SFShell() {}
		SFShell(const SField& f, double t) : _f(f), thickness(t) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFShell>(*this);
		}
		double eval(const Vector& p) const
		{
			double eval = _f->eval(p);
			return std::min(eval + thickness * 0.5, -(eval - thickness * 0.5));
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFTranslate : public ScalarField
	{
		SField _f;
		Vector _x;
	public:
		SFTranslate() {}
		SFTranslate(SField f, Vector x) : _f(f), _x(x) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFTranslate>(*this);
		}
		double eval(const Vector& p) const
		{
			return _f->eval(p - _x);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFScale : public ScalarField
	{
		SField _f;
		double _s;
	public:
		SFScale() {}
		SFScale(SField f, double s) : _f(f), _s(s) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFScale>(*this);
		}
		double eval(const Vector& p) const
		{
			return _f->eval(p / _s);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFRotate : public ScalarField
	{
		SField _f;
		Vector axis;
		double angle;
	public:
		SFRotate() {}
		SFRotate(SField f, Vector ax, double an) : _f(f), axis(ax.unitvector()), angle(-an * M_PI / 180.0) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFRotate>(*this);
		}
		double eval(const Vector& p) const
		{
			double A = cos(angle);
			double B = (axis * p) * (1 - A);
			double C = sin(angle);
			Vector rotated_p = (A * p) + (B * axis) + (C * (p ^ axis));
			return _f->eval(rotated_p);
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFMask : public ScalarField
	{
		SField _f;
	public:
		SFMask() {}
		SFMask(SField f) : _f(f) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFMask>(*this);
		}
		double eval(const Vector& p) const
		{
			return _f->eval(p) > 0.0 ? 1.0 : 0.0;
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};

	class SFClamp : public ScalarField
	{
		SField _f;
		double low, high;
	public:
		SFClamp() {}
		SFClamp(SField f, double l, double h) : _f(f), low(l), high(h) {}
		virtual std::unique_ptr<ScalarField> clone() const override
		{
			return std::make_unique<SFClamp>(*this);
		}
		double eval(const Vector& p) const
		{
			double eval = _f->eval(p);
			return eval < low ? high : eval > high ? high : eval;
		}
		Vector grad(const Vector& p) const
		{
			return Vector();
		}
	};
}
