#include "ScalarGrid.hpp"
//remove lux
using namespace lux;

unsigned int ScalarGrid::getIndex(unsigned int i, unsigned int j, unsigned int k) const
{
	i = std::min(i, Nx - 1);
	j = std::min(j, Ny - 1);
	k = std::min(k, Nz - 1);

	unsigned int index = i + (Nx * j) + (Nx * Ny * k);
	if (index >= gridData.size())
		throw std::runtime_error("grid index out of range");
	//return std::numeric_limits<int>::max();
	return index;
}

bool ScalarGrid::withinGrid(lux::Vector p) const
{
	return	p.X() >= llc.X() && p.X() <= urc.X() &&
		p.Y() >= llc.Y() && p.Y() <= urc.Y() &&
		p.Z() >= llc.Z() && p.Z() <= urc.Z();
}

void ScalarGrid::writeGrid(std::string fileName)
{
	std::cout << "Writing File: " << fileName << " . . . . ";
	if (gridData.size() == 0)
		throw std::runtime_error("DSM empty");

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

void ScalarGrid::readGrid(std::string fileName)
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

//FieldDataType
const float lux::ScalarGrid::eval(const Vector & p) const
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

	float c00 = gridData[getIndex(i, j, k)] * (1 - wi) + gridData[getIndex(i + 1, j, k)] * wi;
	float c01 = gridData[getIndex(i, j, k + 1)] * (1 - wi) + gridData[getIndex(i + 1, j, k + 1)] * wi;
	float c10 = gridData[getIndex(i, j + 1, k)] * (1 - wi) + gridData[getIndex(i + 1, j + 1, k)] * wi;
	float c11 = gridData[getIndex(i, j + 1, k + 1)] * (1 - wi) + gridData[getIndex(i + 1, j + 1, k + 1)] * wi;

	float c0 = c00 * (1 - wj) + c10 * wj;
	float c1 = c01 * (1 - wj) + c11 * wj;

	float c = c0 * (1 - wk) + c1 * wk;
	return c;
}

void ScalarGrid::stamp(lux::SField s)
{
	gridData.resize(Nx * Ny * Nz, 0);

#pragma omp parallel for
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				lux::Vector p = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
				gridData[getIndex(i, j, k)] = s->eval(p);
			}
		}
	}

	std::cout << "Field stamped\n";
}

void ScalarGrid::stampWithDisplacement(lux::SField s, NoiseParams& param)
{
	gridData.resize(Nx * Ny * Nz, 0);

	FSPN fspn = FSPN(param.octaves, param.freq, param.fJump, 2);
	float gamma = param.wedgeSpecific;
	float scalingFact = 0.1f;

#pragma omp parallel for
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				lux::Vector p = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
				gridData[getIndex(i, j, k)] = s->eval(p) + scalingFact * pow(abs(fspn.eval(p.unitvector())), gamma);
			}
		}
	}

	std::cout << "Field stamped with displacement\n";
}

void ScalarGrid::pyroDisplace(NoiseParams& param)
{
	gridData.resize(Nx * Ny * Nz, 0);

	FSPN fspn = FSPN(param.octaves, param.freq, param.fJump, 2);
	float gamma = param.wedgeSpecific;
	float scalingFact = 0.2f;

#pragma omp parallel for
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				lux::Vector p = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
				gridData[getIndex(i, j, k)] += scalingFact * pow(abs(fspn.eval(p.unitvector())), gamma);
			}
		}
	}

	std::cout << "Field displaced\n";
}

void ScalarGrid::levelSet(Triangles& triangles)
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
				lux::Vector p = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
				float dist = std::numeric_limits<float>::max();
				int closestPointIndex = -1;
				lux::Vector closestPoint;

				for (int i = 0; i < triangles.size(); ++i)
				{
					lux::Vector cp = triangles[i]->closestPoint(p);
					float mag = (cp - p).magnitude();

					if (mag <= dist)
					{
						dist = mag;
						closestPointIndex = i;
						closestPoint = cp;
					}
				}

				lux::Vector closestTriangleNormal = triangles[closestPointIndex]->n.unitvector();
				float dotProd = (p - closestPoint) * closestTriangleNormal;
				if (dotProd >= 0)
					dist = -dist;

				gridData[getIndex(i, j, k)] = dist;
			}
		}
	}

	std::cout << "Levelset generated\n";
}