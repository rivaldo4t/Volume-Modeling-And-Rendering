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
private:
	unsigned int Nx, Ny, Nz;
	lux::Vector llc, urc;
	double delta_grid = 0.01;
	std::vector<double> gridData;
	double defaultVal = 0.0;
public:
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
		unsigned int index = i + (Nx * j) + (Nx * Ny * k);
		if (index > gridData.size())
			throw std::runtime_error("grid index out of range");
		return index;
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

	bool withinGrid(lux::Vector p) const
	{
		return	p.X() >= llc.X() && p.X() <= urc.X() &&
				p.Y() >= llc.Y() && p.Y() <= urc.Y() &&
				p.Z() >= llc.Z() && p.Z() <= urc.Z();
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

		//// debugging
		//auto index = getIndex(i, j, k);
		//double d = gridData[index];
		//return d;

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
					//int closestPointIndex = -1;
					std::vector<std::pair<int, lux::Vector>> triangleIndexAndClosestPoint;
					//lux::Vector closestPoint;
					/*if (p.X() == -2.0 && p.Y() == 0.1 && p.Z() == 1.1)
					{
						double t = 5;
						t *= 5;
					}*/
					for (int i = 0; i < triangles.size(); ++i)
					{
						lux::Vector cp = triangles[i]->closestPoint(p);
						double mag = (cp - p).magnitude();
						//if (mag < dist)
						if (mag <= dist)
						{
							dist = mag;
							/*closestPointIndex = i;
							closestPoint = cp;*/
							// tempIndices.push_back(i);
							triangleIndexAndClosestPoint.push_back(std::make_pair(i, cp));
						}
					}
					//lux::Vector combinedNormal;
					bool atleastOnePositive = false;
					for (auto indexAndPoint : triangleIndexAndClosestPoint)
					{
						//combinedNormal = triangles[tri]->n;
						// auto n = triangles[indexAndPoint.first]->n;
						// auto dotpro = (p - indexAndPoint.second) * n;
						// if (dotpro >= 0)
						if ((p - indexAndPoint.second) * triangles[indexAndPoint.first]->n >= 0)
						{
							atleastOnePositive = true;
							break;
						}
					}
					if (atleastOnePositive)
						dist = -dist;
					//combinedNormal = triangles[closestPointIndex]->n;
					//double dotpro = (p - closestPoint) * combinedNormal;
					/*if (dotpro >= 0)
					{
						dist = -dist;
					}
					else
					{
						if (dist == -1.13342423423)
							dist = -1.13342423423;
					}*/
					gridData[getIndex(i, j, k)] = dist;
				}
			}
		}
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
