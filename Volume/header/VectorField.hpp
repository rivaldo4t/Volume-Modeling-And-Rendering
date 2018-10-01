#pragma once
#include <iostream>
#include <memory>
#include "Vector.hpp"
#include "Matrix.hpp"

namespace lux
{
	class VectorField
	{
		virtual Vector eval(const Vector& p) const = 0;
		virtual Matrix grad(const Vector& p) const = 0;
		virtual std::unique_ptr<VectorField> clone() const = 0;
	};
	typedef std::shared_ptr<lux::VectorField> VField;

	VField operator-(const VField v)
	{
		VField u;

	}
}