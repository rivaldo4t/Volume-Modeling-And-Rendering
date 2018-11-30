#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <memory>
#include "Field.hpp"

namespace lux
{
	template <typename T>
	class Grid : public Field<T>
	{
	public:
		unsigned int Nx, Ny, Nz;
		float delta_grid = 0.01f;
		T defaultVal;
		Vector llc, urc;
		std::vector<T> gridData;

		Grid() { Nx = 0; Ny = 0; Nz = 0; }
		Grid(Vector o, unsigned int x, unsigned int y, unsigned int z, float delta, T dVal,
			std::vector<T>& data = std::vector<T>())
			: llc(o), Nx(x), Ny(y), Nz(z), delta_grid(delta)
		{
			gridData = std::move(data);
			urc = { llc.X() + (Nx - 1) * delta_grid,
				llc.Y() + (Ny - 1) * delta_grid,
				llc.Z() + (Nz - 1) * delta_grid };
			defaultVal = dVal;
		}
		virtual ~Grid() { gridData.clear(); std::cout << "------------Grid destructor\n"; }

		unsigned int getIndex(unsigned int i, unsigned int j, unsigned int k) const;
		bool withinGrid(lux::Vector p) const;
		virtual const T eval(const Vector & p) const override;
		void writeGrid(std::string fileName);
		void readGrid(std::string fileName);
	};

	template<typename T>
	unsigned int Grid<T>::getIndex(unsigned int i, unsigned int j, unsigned int k) const
	{
		i = std::min(i, Nx - 1);
		j = std::min(j, Ny - 1);
		k = std::min(k, Nz - 1);

		unsigned int index = i + (Nx * j) + (Nx * Ny * k);
		if (index >= gridData.size())
			throw std::runtime_error("grid index out of range");
		return index;
	}

	template<typename T>
	bool Grid<T>::withinGrid(lux::Vector p) const
	{
		return	p.X() >= llc.X() && p.X() <= urc.X() &&
			p.Y() >= llc.Y() && p.Y() <= urc.Y() &&
			p.Z() >= llc.Z() && p.Z() <= urc.Z();
	}

	template<typename T>
	const T Grid<T>::eval(const Vector & p) const
	{
		if (!withinGrid(p))
			return defaultVal;

		lux::Vector toPoint = p - llc;
		float x = toPoint.X();
		float y = toPoint.Y();
		float z = toPoint.Z();

		unsigned int i = floor(x / delta_grid);
		unsigned int j = floor(y / delta_grid);
		unsigned int k = floor(z / delta_grid);

		float wi = (x - i * delta_grid) / delta_grid;
		float wj = (y - j * delta_grid) / delta_grid;
		float wk = (z - k * delta_grid) / delta_grid;

		T gijk = gridData[getIndex(i, j, k)];
		T gi_1jk = gridData[getIndex(i + 1, j, k)];
		T gijk_1 = gridData[getIndex(i, j, k + 1)];
		T gi_1jk_1 = gridData[getIndex(i + 1, j, k + 1)];
		T gij_1k = gridData[getIndex(i, j + 1, k)];
		T gi_1j_1k = gridData[getIndex(i + 1, j + 1, k)];
		T gij_1k_1 = gridData[getIndex(i, j + 1, k + 1)];
		T gi_1j_1k_1 = gridData[getIndex(i + 1, j + 1, k + 1)];

		T c00 = gijk		* (1 - wi) + gi_1jk		* wi;
		T c01 = gijk_1		* (1 - wi) + gi_1jk_1	* wi;
		T c10 = gij_1k		* (1 - wi) + gi_1j_1k	* wi;
		T c11 = gij_1k_1	* (1 - wi) + gi_1j_1k_1	* wi;

		T c0 = c00 * (1 - wj) + c10 * wj;
		T c1 = c01 * (1 - wj) + c11 * wj;

		T c = c0 * (1 - wk) + c1 * wk;
		return c;
	}

	template<typename T>
	void Grid<T>::writeGrid(std::string fileName)
	{
		std::cout << "Writing File: " << fileName << " . . . . ";
		if (gridData.size() == 0)
			throw std::runtime_error("Grid empty");

		std::ofstream ofs(fileName, std::ios::out | std::ofstream::binary);
		if (!ofs)
			throw std::runtime_error("error opening file");

		ofs.write(reinterpret_cast<const char*>(&llc), sizeof(llc));
		ofs.write(reinterpret_cast<const char*>(&Nx), sizeof(Nx));
		ofs.write(reinterpret_cast<const char*>(&Ny), sizeof(Ny));
		ofs.write(reinterpret_cast<const char*>(&Nz), sizeof(Nz));
		ofs.write(reinterpret_cast<const char*>(&delta_grid), sizeof(delta_grid));

		size_t size = gridData.size();
		ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
		ofs.write(reinterpret_cast<const char*>(&gridData[0]), size * sizeof(gridData[0]));
		ofs.close();
		std::cout << "Done\n";
	}

	template<typename T>
	void Grid<T>::readGrid(std::string fileName)
	{
		std::cout << "Reading File: " << fileName << " . . . . ";
		std::ifstream ifs(fileName, std::ios::in | std::ifstream::binary);
		if (!ifs)
			throw std::runtime_error("error opening file");

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
		std::cout << "Done\n";
	}
}