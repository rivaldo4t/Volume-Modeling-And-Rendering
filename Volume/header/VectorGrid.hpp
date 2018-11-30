#pragma once
#include "Grid.hpp"
#include "ScalarGrid.hpp"

namespace lux
{
	class VectorGrid : public Grid<Vector>
	{
	public:
		VectorGrid() : Grid<Vector>() {}
		VectorGrid(lux::Vector o, unsigned int x, unsigned int y, unsigned int z, double delta, Vector dVal = Vector(),
			std::vector<Vector>& data = std::vector<Vector>()) : Grid<Vector>(o, x, y, z, delta, dVal, data) {}
		virtual ~VectorGrid() { std::cout << "------------VectorGrid destructor\n"; }

		void stamp(lux::VField v);
		void gsr()
		{
			std::cout << "Calculating Divergence . . . . ";
			ScalarGrid divergence(llc, Nx, Ny, Nz, delta_grid, 0.0);
			divergence.gridData.clear();
			divergence.gridData.resize(Nx * Ny * Nz, 0.0);
#pragma omp parallel for
			for (int i = 0; i < Nx; ++i)
			{
				for (int j = 0; j < Ny; ++j)
				{
					for (int k = 0; k < Nz; ++k)
					{
						lux::Vector p_i1_jk = llc + lux::Vector(float(i + 1) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
						lux::Vector p_1i_jk = llc + lux::Vector(float(i - 1) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
						lux::Vector pi_j1_k = llc + lux::Vector(float(i) * delta_grid, float(j + 1) * delta_grid, float(k) * delta_grid);
						lux::Vector pi_1j_k = llc + lux::Vector(float(i) * delta_grid, float(j - 1) * delta_grid, float(k) * delta_grid);
						lux::Vector pij_k1 = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k + 1) * delta_grid);
						lux::Vector pij_1k = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k - 1) * delta_grid);
						divergence.gridData[getIndex(i, j, k)] = (eval(p_i1_jk).X() - eval(p_1i_jk).X() +
							eval(pi_j1_k).Y() - eval(pi_1j_k).Y() +
							eval(pij_k1).Z() - eval(pij_1k).Z()) / (2 * delta_grid);
					}
				}
			}
			std::cout << "Done\n";

			std::cout << "Calculating Pressure . . . . ";
			ScalarGrid pressure(llc, Nx, Ny, Nz, delta_grid, 0.0);
			pressure.gridData.clear();
			pressure.gridData.resize(Nx * Ny * Nz, 0.0);
			int numiter = 10;
			while (numiter--)
			{
#pragma omp parallel for
				for (int i = 0; i < Nx; ++i)
				{
					for (int j = 0; j < Ny; ++j)
					{
						for (int k = 0; k < Nz; ++k)
						{
							pressure.gridData[getIndex(i, j, k)] =
								(delta_grid * delta_grid * divergence.gridData[getIndex(i, j, k)] +
									pressure.gridData[getIndex(i + 1, j, k)] + pressure.gridData[getIndex(i - 1, j, k)] +
									pressure.gridData[getIndex(i, j + 1, k)] + pressure.gridData[getIndex(i, j - 1, k)] +
									pressure.gridData[getIndex(i, j, k + 1)] + pressure.gridData[getIndex(i, j, k - 1)]) / 6.0;
						}
					}
				}
			}
			std::cout << "Done\n";
			//divergence.gridData.clear();

			std::cout << "Updating Velocity . . . . ";
#pragma omp parallel for
			for (int i = 0; i < Nx; ++i)
			{
				for (int j = 0; j < Ny; ++j)
				{
					for (int k = 0; k < Nz; ++k)
					{
						gridData[getIndex(i, j, k)] -=
							Vector(
								pressure.gridData[getIndex(i + 1, j, k)] - pressure.gridData[getIndex(i - 1, j, k)],
								pressure.gridData[getIndex(i, j + 1, k)] - pressure.gridData[getIndex(i, j - 1, k)],
								pressure.gridData[getIndex(i, j, k + 1)] - pressure.gridData[getIndex(i, j, k - 1)]
							) / (2 * delta_grid);
					}
				}
			}
			std::cout << "Done\n";
			//pressure.gridData.clear();
		}
	};
}
