#include "Grid.hpp"

unsigned int Grid::getIndex(unsigned int i, unsigned int j, unsigned int k) const
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

bool Grid::withinGrid(lux::Vector p) const
{
	return	p.X() >= llc.X() && p.X() <= urc.X() &&
		p.Y() >= llc.Y() && p.Y() <= urc.Y() &&
		p.Z() >= llc.Z() && p.Z() <= urc.Z();
}

void Grid::stamp(lux::SField s)
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

	std::cout << "Field stamped\n";
}

void Grid::stampWithDisplacement(lux::SField s)
{
	gridData.resize(Nx * Ny * Nz, 0);

	static int _oct, _gam, _freq, _fjump, _rough;
	std::vector<double> oct = { 1, 2, 3, 5, 6 };
	std::vector<double> gam = { 0.33, 1, 1.6, 2 };
	std::vector<double> freq = { 1, 5, 10, 50 };
	std::vector<double> fjump = { 1.2, 2 };
	std::vector<double> rough = { 0.5, 1.0, 2.5, 4 };

	FSPN fspn = FSPN(oct[_oct], freq[_freq], fjump[_fjump], rough[_rough]);
	float gamma = gam[_gam];
	float scalingFact = 0.6;
	/*fspn = FSPN(1.6, 5, 2, 1.6);
	gamma = 1.6;*/

	std::cout << "--------------------\n";
	std::cout << "octaves:\t" << fspn.octaves << std::endl;
	std::cout << "freq:\t\t" << fspn.freq << std::endl;
	std::cout << "fjump:\t\t" << fspn.fJump << std::endl;
	std::cout << "roughness:\t" << fspn.roughness << std::endl;
	std::cout << "gamma:\t\t" << gamma << std::endl;
	std::cout << "--------------------\n";

	_rough++;
	if (_rough == rough.size())
		_fjump++;
	if (_fjump == fjump.size())
		_freq++;
	if (_freq == freq.size())
		_gam++;
	if (_gam == gam.size())
		_oct++;

	_rough %= rough.size();
	_fjump %= fjump.size();
	_freq %= freq.size();
	_gam %= gam.size();
	_oct %= oct.size();

#pragma omp parallel for
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				lux::Vector p = llc + lux::Vector(double(i) * delta_grid, double(j) * delta_grid, double(k) * delta_grid);
				gridData[getIndex(i, j, k)] = s->eval(p) +
					scalingFact * pow(abs(fspn.eval(p.unitvector())), gamma);
			}
		}
	}

	std::cout << "Field stamped\n";
}

void Grid::levelSet(Triangles& triangles)
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

	std::cout << "Levelset generated\n";
}

double Grid::eval(lux::Vector p) const
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

void Grid::writeGrid(std::string fileName)
{
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
	std::cout << "file write: " << fileName << std::endl;
}

void Grid::readGrid(std::string fileName)
{
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
	std::cout << "file read: " << fileName << std::endl;
}