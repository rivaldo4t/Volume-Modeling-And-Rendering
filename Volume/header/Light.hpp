#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include "Vector.hpp"
#include "Color.hpp"
#include "ScalarField.hpp"
#include "ScalarGrid.hpp"

namespace lux
{
	class Light : public ScalarGrid
	{
	private:
		lux::Vector pos;
		lux::Color color;
		double delta_s;
	public:
		Light() {}
		Light(lux::Vector p, lux::Vector l, int nx, unsigned int ny, unsigned int nz, double d, double ds = 0.01, lux::Color c = lux::Color(0.8, 0.8, 0.8, 1.0))
			: ScalarGrid(l, nx, ny, nz, d), pos(p), delta_s(ds), color(c) {}
		lux::Vector getPosition() const { return pos; }
		lux::Color getColor() const { return color; }
		void setColor(lux::Color c) { color = c; }
		void computeDSM(lux::SField density);
		void writeDSM(std::string fileName);
		void readDSM(std::string fileName);
	};
}
