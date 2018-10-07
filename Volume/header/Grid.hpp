#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <memory>
#include "ScalarField.hpp"
#include "Triangle.hpp"

class Grid
{
// private:
public:
	unsigned int Nx, Ny, Nz;
	double delta_grid = 0.01;
	double defaultVal = 0.0;
	lux::Vector llc, urc;
	std::vector<double> gridData;
// public:
	Grid() { Nx = 0; Ny = 0; Nz = 0; }
	Grid(lux::Vector o, unsigned int x, unsigned int y, unsigned int z, double delta, std::vector<double>& data = std::vector<double>())
		: llc(o), Nx(x), Ny(y), Nz(z), delta_grid(delta)
	{ 
		gridData = std::move(data); 
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
		if (index > gridData.size())
			throw std::runtime_error("grid index out of range");
		return index;
	}

	bool withinGrid(lux::Vector p) const
	{
		return	p.X() >= llc.X() && p.X() <= urc.X() &&
				p.Y() >= llc.Y() && p.Y() <= urc.Y() &&
				p.Z() >= llc.Z() && p.Z() <= urc.Z();
	}

	void stamp(lux::SField s)
	{
		gridData.resize(Nx * Ny * Nz, 0);

#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector p = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					gridData[getIndex(i, j, k)] = s->eval(p);
				}
			}
		}
	}

	void levelSet(Triangles& triangles)
	{
		gridData.clear();
		gridData.resize(Nx * Ny * Nz, 0);

#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector p = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					double dist = std::numeric_limits<double>::max();
					int closestPointIndex = -1;
					lux::Vector closestPoint;
					
					for (int i = 0; i < triangles.size(); ++i)
					{
						lux::Vector cp = triangles[i]->closestPoint(p);
						double mag = (cp - p).magnitude();

						if (mag <= dist)
						{
							dist = mag;
							closestPointIndex = i;
							closestPoint = cp;
						}
					}
					
					lux::Vector closestTriangleNormal = triangles[closestPointIndex]->n.unitvector();
					double dotProd = (p - closestPoint) * closestTriangleNormal;
					if (dotProd >= 0)
						dist = -dist;
					
					gridData[getIndex(i, j, k)] = dist;
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

		double c00 = gridData[getIndex(i, j, k)] * (1 - wi) + gridData[getIndex(i + 1, j, k)] * wi;
		double c01 = gridData[getIndex(i, j, k + 1)] * (1 - wi) + gridData[getIndex(i + 1, j, k + 1)] * wi;
		double c10 = gridData[getIndex(i, j + 1, k)] * (1 - wi) + gridData[getIndex(i + 1, j + 1, k)] * wi;
		double c11 = gridData[getIndex(i, j + 1, k + 1)] * (1 - wi) + gridData[getIndex(i + 1, j + 1, k + 1)] * wi;

		double c0 = c00 * (1 - wj) + c10 * wj;
		double c1 = c01 * (1 - wj) + c11 * wj;

		double c = c0 * (1 - wk) + c1 * wk;
		return c;
	}

	void writelevelSet(std::string fileName)
	{
		if (gridData.size() == 0)
			throw std::runtime_error("DSM empty");

		std::ofstream ofs(fileName, std::ios::out | std::ofstream::binary);
		if (!ofs)
			throw std::runtime_error("error opening file");

		size_t size = gridData.size();
		ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
		ofs.write(reinterpret_cast<const char*>(&gridData[0]), size * sizeof(gridData[0]));
		ofs.close();
		std::cout << "file write: " << fileName << std::endl;
	}

	void readlevelSet(std::string fileName)
	{
		std::ifstream ifs(fileName, std::ios::in | std::ifstream::binary);
		if (!ifs)
			throw std::runtime_error("error opening file");

		size_t size;
		ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
		gridData.clear();
		gridData.resize(size, 0);
		ifs.read(reinterpret_cast<char*>(&gridData[0]), size * sizeof(gridData[0]));
		ifs.close();
		std::cout << "file read: " << fileName << std::endl;
	}
};

class GridUnion : public Grid
{
	std::shared_ptr<Grid> f, g;
public:
	double eval(lux::Vector p)
	{
		return std::max(f->eval(p), g->eval(p));
	}
};
