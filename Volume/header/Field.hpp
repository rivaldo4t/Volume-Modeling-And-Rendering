#pragma once
#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Color.hpp"
#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <algorithm>

namespace lux
{
	template <typename T>
	struct GradType
	{
	   typedef int GType;
	};

	template<>
	struct GradType<float>
	{
	   typedef Vector GType;
	};

	template<>
	struct GradType<Vector>
	{
	   typedef Matrix GType;
	};

	template <typename T>
	class Field 
	{
	  public:

		Field(){}

		virtual ~Field(){}

		typedef T FieldDataType;
		typedef typename GradType<T>::GType FieldGradType;

		virtual const FieldDataType eval(const Vector& P) const { FieldDataType base{}; return base; }
		virtual const FieldGradType grad(const Vector& P) const { FieldGradType base{}; return base; }
	};

	typedef std::shared_ptr<lux::Field<float>> SField;
	typedef std::shared_ptr<lux::Field<Vector>> VField;
	typedef std::shared_ptr<lux::Field<Matrix>> MField;
	//typedef std::shared_ptr<lux::Field<Color>> CField;
}