#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include "Vector.hpp"
#include "Color.hpp"
#include "ScalarField.hpp"
#include "Grid.hpp"

class Light : public Grid
{
public:
	lux::Vector pos;
	lux::Color color;
	double delta_s;

	Light() {}
	Light(lux::Vector p, lux::Vector l, int nx, unsigned int ny, unsigned int nz, double d, double ds = 0.01, lux::Color c = lux::Color(0.8, 0.8, 0.8, 1.0))
		: Grid(l, nx, ny, nz, d), pos(p), delta_s(ds), color(c) {}

	void computeDSM(lux::SField density)
	{
		gridData.resize(Nx * Ny * Nz, 0);
		
#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector marchStartPos = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					lux::Vector nL = pos - marchStartPos;
					double sFar = nL.magnitude(), sNear = 0, s = 0;
					nL.normalize();
					lux::Vector X = marchStartPos + sNear * nL;

					if (density->eval(X) > 0)
					{
						while (s <= sFar)
						{
							X += delta_s * nL;
							double d = density->eval(X);
							if (d > 0)
								gridData[getIndex(i, j, k)] += d * delta_s;
							s += delta_s;
						}
					}
				}
			}
		}
	}

	void computeDSM(const Grid& g)
	{
		gridData.resize(Nx * Ny * Nz, 0);

#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector marchStartPos = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					lux::Vector nL = pos - marchStartPos;
					double sFar = nL.magnitude(), sNear = 0, s = 0;
					nL.normalize();
					lux::Vector X = marchStartPos + sNear * nL;

					if (g.eval(X) > 0)
					{
						while (s <= sFar)
						{
							X += delta_s * nL;
							double d = g.eval(X);
							if (d > 0)
								gridData[getIndex(i, j, k)] += d * delta_s;
							s += delta_s;
						}
					}
				}
			}
		}
	}

	void writeDSM(std::string fileName)
	{
		if (gridData.size() == 0)
			throw std::runtime_error("DSM empty");

		std::ofstream ofs(fileName, std::ios::out | std::ofstream::binary);
		if (!ofs)
			throw std::runtime_error("error opening file");

		ofs.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
		ofs.write(reinterpret_cast<const char*>(&delta_s), sizeof(delta_s));
		ofs.write(reinterpret_cast<const char*>(&color), sizeof(color));
		
		ofs.write(reinterpret_cast<const char*>(&llc), sizeof(llc));
		ofs.write(reinterpret_cast<const char*>(&Nx), sizeof(Nx));
		ofs.write(reinterpret_cast<const char*>(&Ny), sizeof(Ny));
		ofs.write(reinterpret_cast<const char*>(&Nz), sizeof(Nz));
		ofs.write(reinterpret_cast<const char*>(&delta_grid), sizeof(delta_grid));

		size_t size = gridData.size();
		ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
		ofs.write(reinterpret_cast<const char*>(&gridData[0]), size * sizeof(gridData[0]));
		ofs.close();
		std::cout << "file write: " << fileName << std::endl;
	}

	void readDSM(std::string fileName)
	{
		std::ifstream ifs(fileName, std::ios::in | std::ifstream::binary);
		if (!ifs)
			throw std::runtime_error("error opening file");

		ifs.read(reinterpret_cast<char*>(&pos), sizeof(pos));
		ifs.read(reinterpret_cast<char*>(&delta_s), sizeof(delta_s));
		ifs.read(reinterpret_cast<char*>(&color), sizeof(color));

		ifs.read(reinterpret_cast<char*>(&llc), sizeof(llc));
		ifs.read(reinterpret_cast<char*>(&Nx), sizeof(Nx));
		ifs.read(reinterpret_cast<char*>(&Ny), sizeof(Ny));
		ifs.read(reinterpret_cast<char*>(&Nz), sizeof(Nz));
		ifs.read(reinterpret_cast<char*>(&delta_grid), sizeof(delta_grid));

		urc = { llc.X() + (Nx - 1) * delta_grid,
			llc.Y() + (Ny - 1) * delta_grid,
			llc.Z() + (Nz - 1) * delta_grid };

		size_t size;
		ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
		gridData.clear();
		gridData.resize(size, 0);
		ifs.read(reinterpret_cast<char*>(&gridData[0]), size * sizeof(gridData[0]));
		ifs.close();
		std::cout << "file read: " << fileName << std::endl;
	}
};
