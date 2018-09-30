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
private:
//public:
	lux::Vector pos;
	lux::Color color;
	std::vector<double> deepShadowMap;
public:
	Light() {}
	Light(lux::Vector p, lux::Color c = lux::Color(1.0, 1.0, 1.0, 1.0)) : pos(p), color(c) {}

	void computeDSM(lux::SField density)
	{
		lux::Vector llc(-5, -5, -5);
		double delta_grid = 0.01;
		int Nx = 1000, Ny = 1000, Nz = 1000;
		double delta_s = 0.01;
		/*Nx = Ny = Nz = 100;
		delta_s = 0.1;*/
		deepShadowMap.resize(Nx * Ny * Nz, 0);
		
#pragma omp parallel for
		for (int i = 0; i < Nx; ++i)
		{
			//std::cout << i << std::endl;
			for (int j = 0; j < Ny; ++j)
			{
				for (int k = 0; k < Nz; ++k)
				{
					lux::Vector marchStartPos = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
					lux::Vector nL = pos - marchStartPos;
					double sFar = nL.magnitude(), sNear = 0, s = 0;
					nL.normalize();
					lux::Vector X = marchStartPos + sNear * nL;
					int index = i + (Nx * j) + (Nx * Ny * k);

					if (density->eval(X) > 0)
					{
						while (s <= sFar)
						{
							X += delta_s * nL;
							double d = density->eval(X);
							if (d > 0)
								deepShadowMap[index] += d * delta_s;
							s += delta_s;
						}
					}
				}
			}
		}

		/*for (auto i : deepShadowMap)
			if (i > 0)
				std::cout << i << std::endl;*/
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
