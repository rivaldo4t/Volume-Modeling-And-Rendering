#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include "Vector.hpp"
#include "Color.hpp"
#include "ScalarField.hpp"

class Light
{
public:
	lux::Vector pos;
	lux::Color color;
	std::vector<double> deepShadowMap;
	lux::Vector llc, urc;
	unsigned int Nx, Ny, Nz;
	double delta_grid;
	double delta_s;
	double defaultVal = 0.0;

	Light() {}
	Light(lux::Vector p, lux::Vector l, int nx, unsigned int ny, unsigned int nz, double d, double ds = 0.01, lux::Color c = lux::Color(0.8, 0.8, 0.8, 1.0))
		: pos(p), llc(l), Nx(nx), Ny(ny), Nz(nz), delta_grid(d), delta_s(ds), color(c) 
	{
		urc = { llc.X() + (Nx - 1) * delta_grid,
				llc.Y() + (Ny - 1) * delta_grid,
				llc.Z() + (Nz - 1) * delta_grid };
	}

	unsigned int getIndex(unsigned int i, unsigned int j, unsigned int k) const
	{
		i = std::min(i, Nx - 1);
		j = std::min(j, Ny - 1);
		k = std::min(k, Nz - 1);

		unsigned int index = i + (Nx * j) + (Nx * Ny * k);
		if (index > deepShadowMap.size())
			throw std::runtime_error("grid index out of range");
		return index;
	}

	bool withinGrid(lux::Vector p) const
	{
		return	p.X() >= llc.X() && p.X() <= urc.X() &&
			p.Y() >= llc.Y() && p.Y() <= urc.Y() &&
			p.Z() >= llc.Z() && p.Z() <= urc.Z();
	}

	void computeDSM(lux::SField density)
	{
		deepShadowMap.resize(Nx * Ny * Nz, 0);
		
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
								deepShadowMap[getIndex(i, j, k)] += d * delta_s;
							s += delta_s;
						}
					}
				}
			}
		}
	}

	void computeDSM2(const Grid& g)
	{
		deepShadowMap.resize(Nx * Ny * Nz, 0);

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
								deepShadowMap[getIndex(i, j, k)] += d * delta_s;
							s += delta_s;
						}
					}
				}
			}
		}
	}

	double eval(lux::Vector p) const
	{
		if (!withinGrid(p))
			return defaultVal;

		lux::Vector toPoint = p - llc;
		double x = toPoint.X();
		double y = toPoint.Y();
		double z = toPoint.Z();

		unsigned int i = floor(x / delta_grid);
		unsigned int j = floor(y / delta_grid);
		unsigned int k = floor(z / delta_grid);

		double wi = (x - i * delta_grid) / delta_grid;
		double wj = (y - j * delta_grid) / delta_grid;
		double wk = (z - k * delta_grid) / delta_grid;

		double c00 = deepShadowMap[getIndex(i, j, k)] * (1 - wi) + deepShadowMap[getIndex(i + 1, j, k)] * wi;
		double c01 = deepShadowMap[getIndex(i, j, k + 1)] * (1 - wi) + deepShadowMap[getIndex(i + 1, j, k + 1)] * wi;
		double c10 = deepShadowMap[getIndex(i, j + 1, k)] * (1 - wi) + deepShadowMap[getIndex(i + 1, j + 1, k)] * wi;
		double c11 = deepShadowMap[getIndex(i, j + 1, k + 1)] * (1 - wi) + deepShadowMap[getIndex(i + 1, j + 1, k + 1)] * wi;

		double c0 = c00 * (1 - wj) + c10 * wj;
		double c1 = c01 * (1 - wj) + c11 * wj;

		double c = c0 * (1 - wk) + c1 * wk;
		return c;
	}

	void writeDSM(std::string fileName)
	{
		if (deepShadowMap.size() == 0)
			throw std::runtime_error("DSM empty");

		std::ofstream ofs(fileName, std::ios::out | std::ofstream::binary);
		if (!ofs)
			throw std::runtime_error("error opening file");

		size_t size = deepShadowMap.size();
		ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
		ofs.write(reinterpret_cast<const char*>(&deepShadowMap[0]), size * sizeof(deepShadowMap[0]));
		ofs.close();
		std::cout << "file write: " << fileName << std::endl;
	}

	void readDSM(std::string fileName)
	{
		std::ifstream ifs(fileName, std::ios::in | std::ifstream::binary);
		if (!ifs)
			throw std::runtime_error("error opening file");

		size_t size;
		ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
		deepShadowMap.clear();
		deepShadowMap.resize(size, 0);
		ifs.read(reinterpret_cast<char*>(&deepShadowMap[0]), size * sizeof(deepShadowMap[0]));
		ifs.close();
		std::cout << "file read: " << fileName << std::endl;
	}
};
