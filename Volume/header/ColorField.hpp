#pragma once
#include <iostream>
#include <memory>
#include "Vector"
#include "Color.hpp"
#include "ScalarField.hpp"

namespace lux
{
	class ColorField
	{
		SField _f, _g, _h;
	public:
		ColorField() {}
		ColorField(SField f, SField g, SField h) : _f(f), _g(g), _h(h) {}
		Color eval(const Vector& p) const
		{
			auto mask1 = std::make_shared<SFMask>(_f);
			auto mask2 = std::make_shared<SFMask>(_g);
			auto mask3 = std::make_shared<SFMask>(_h);

			Color red(1.0, 0.2, 0.2, 1.0);
			Color green(0.2, 1.0, 0.2, 1.0);
			Color blue(0.2, 0.2, 1.0, 1.0);

			Color c = red * mask1->eval(p);
			c += green * mask2->eval(p);
			c += blue * mask3->eval(p);

			return c;
		}
	};
	typedef std::shared_ptr<lux::ColorField> CField;

	class CFRotate : public ColorField
	{
		CField _f;
		Vector axis;
		double angle;
	public:
		CFRotate() {}
		CFRotate(CField f, Vector ax, double an) : _f(f), axis(ax.unitvector()), angle(-an * M_PI / 180.0) {}
		Color eval(const Vector& p) const
		{
			double A = cos(angle);
			double BB = (axis * p) * (1 - A);
			double C = sin(angle);
			Vector rotated_p = (A * p) + (BB * axis) + (C * (p ^ axis));
			return _f->eval(rotated_p);
		}
	};
}