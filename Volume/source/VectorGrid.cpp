#include "VectorGrid.hpp"
using namespace lux;

void VectorGrid::stamp(lux::VField v)
{
	std::cout << "Stamping Field . . . . ";
	gridData.clear();
	gridData.resize(Nx * Ny * Nz, Vector());

#pragma omp parallel for
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				lux::Vector p = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
				gridData[getIndex(i, j, k)] = v->eval(p);
			}
		}
	}

	std::cout << "Done\n";
}
