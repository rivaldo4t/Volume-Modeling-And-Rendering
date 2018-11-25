#include "ScalarGrid.hpp"
//remove lux
using namespace lux;

void ScalarGrid::stamp(lux::SField s)
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
				gridData[getIndex(i, j, k)] = s->eval(p);
			}
		}
	}

	std::cout << "Field stamped\n";
}

void ScalarGrid::stampWithDisplacement(lux::SField s, NoiseParams& param)
{
	gridData.clear();
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
	gridData.clear();
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

void ScalarGrid::finiteDistanceGradient()
{
	static int N = 1;
	static float deltaX = 0.09, deltaY = 0.09, deltaZ = 0.09;
	static std::vector<float> coeff = { 0.5f };
	//static int N = 2;
	//static std::vector<float> coeff = { 0.57142857142857151f, -0.071428571428571438f };
	//static float deltaX = 0.09, deltaY = 0.09, deltaZ = 0.09;
	/*static int N = 4;
	static float deltaX = 0.09, deltaY = 0.09, deltaZ = 0.09;
	static std::vector<float> coeff = { 0.63039399624765435f, -0.15064623723160311f, 0.02101573900354391f, -0.00076349801959558095f };*/
	/*static int N = 7;
	static float deltaX = 0.08, deltaY = 0.08, deltaZ = 0.05;
	static std::vector<float> coeff = { 0.66418180822321082f, -0.20286437997462345f, 0.042119156095069443f,
		-0.0035797086604024858f, 0.00014591197349313843f, -2.8080248121598451e-06f, 2.0367364328490342e-08f };*/

	std::cout << "Computing Gradient of grid . . . . ";
	gridGradData.resize(Nx * Ny * Nz, Vector());

#pragma omp parallel for
	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < Nz; ++k)
			{
				lux::Vector p = llc + lux::Vector(float(i) * delta_grid, float(j) * delta_grid, float(k) * delta_grid);
				int index = getIndex(i, j, k);
				for (int n = -N; n <= N; ++n)
				{
					if (n == 0)
						continue;
					float a = coeff[abs(n) - 1];
					a = n < 0 ? -a : a;
					gridGradData[index] += a * Vector(	eval(p + n * deltaX * Vector(1, 0, 0)) / deltaX,
														eval(p + n * deltaY * Vector(0, 1, 0)) / deltaY,
														eval(p + n * deltaZ * Vector(0, 0, 1)) / deltaZ);
				}
			}
		}
	}

	std::cout << "Done\n";
}

const Vector ScalarGrid::grad(const Vector & p) const
{
	if (!withinGrid(p))
		//return Vector(); 
		return -(p - Vector()) / (p - Vector()).magnitude();

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

	Vector vijk			= gridGradData[getIndex(i, j, k)];
	Vector vi_1jk		= gridGradData[getIndex(i + 1, j, k)];
	Vector vijk_1		= gridGradData[getIndex(i, j, k + 1)];
	Vector vi_1jk_1		= gridGradData[getIndex(i + 1, j, k + 1)];
	Vector vij_1k		= gridGradData[getIndex(i, j + 1, k)];
	Vector vi_1j_1k		= gridGradData[getIndex(i + 1, j + 1, k)];
	Vector vij_1k_1		= gridGradData[getIndex(i, j + 1, k + 1)];
	Vector vi_1j_1k_1	= gridGradData[getIndex(i + 1, j + 1, k + 1)];

	Vector c00 = vijk		* (1 - wi) + vi_1jk		* wi;
	Vector c01 = vijk_1		* (1 - wi) + vi_1jk_1	* wi;
	Vector c10 = vij_1k		* (1 - wi) + vi_1j_1k	* wi;
	Vector c11 = vij_1k_1	* (1 - wi) + vi_1j_1k_1	* wi;

	Vector c0 = c00 * (1 - wj) + c10 * wj;
	Vector c1 = c01 * (1 - wj) + c11 * wj;

	Vector c = c0 * (1 - wk) + c1 * wk;
	return c;
}