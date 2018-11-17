#pragma once
#include <iostream>
#include <memory>
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Fspn.hpp"

namespace lux
{
	class VectorField
	{
	public:
		virtual Vector eval(const Vector& p) const = 0;
		virtual Matrix grad(const Vector& p) const = 0;
		virtual std::unique_ptr<VectorField> clone() const = 0;
	};
	typedef std::shared_ptr<lux::VectorField> VField;

	class VFIdentity : public VectorField
	{
	public:
		VFIdentity() {}
		virtual std::unique_ptr<VectorField> clone() const override
		{
			return std::make_unique<VFIdentity>(*this);
		}
		Vector eval(const Vector& p) const
		{
			return p;
		}
		Matrix grad(const Vector& p) const
		{
			return Matrix();
		}
	};

	class VFRandom : public VectorField
	{
	private:
		FSPN f;
		Vector delta = Vector(0.1, 0.1, 0.1);
	public:
		VFRandom() {}
		virtual std::unique_ptr<VectorField> clone() const override
		{
			return std::make_unique<VFRandom>(*this);
		}
		Vector eval(const Vector& p) const
		{
			return Vector(f.eval(p), f.eval(p - delta), f.eval(p + delta));
		}
		Matrix grad(const Vector& p) const
		{
			return Matrix();
		}
	};
}