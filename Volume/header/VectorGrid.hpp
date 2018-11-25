#pragma once
#include "Grid.hpp"

namespace lux
{
	class VectorGrid : public Grid<Vector>
	{
	public:
		VectorGrid() : Grid<Vector>() {}
		VectorGrid(lux::Vector o, unsigned int x, unsigned int y, unsigned int z, double delta, Vector dVal = Vector(),
			std::vector<Vector>& data = std::vector<Vector>()) : Grid<Vector>(o, x, y, z, delta, dVal, data) {}

		void stamp(lux::VField v);
	};
}
