#pragma once
#include <iostream>
#include <memory>
#include "Vector.hpp"
#include "Matrix.hpp"

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

	class IdentityVectorField : public VectorField
	{
	public:
		IdentityVectorField() {}
		virtual std::unique_ptr<VectorField> clone() const override
		{
			return std::make_unique<IdentityVectorField>(*this);
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
}