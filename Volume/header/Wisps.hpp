#pragma once
#include <random>
#include "Grid.hpp"
#include "FractalSummedPerlinNoise.hpp"

class WispDot
{
private:
	lux::Vector p;
	lux::Vector offset;
	float pScale;
	float dScale;
	float clump;
	float density;
	FSPN fspn1, fspn2;
	lux::Vector generatedPosition;
public:
	WispDot()
	{ 
		p = lux::Vector(0.0, 0.0, 0.0);
		offset = lux::Vector(0.5, 0.0, 0.0);
		pScale = 0.5;
		dScale = 0.2;
		clump = 1;
		density = 1;
		fspn1 = FSPN(); fspn2 = FSPN();
	}

	lux::Vector getPos() { return generatedPosition; }
	float getDensity() { return density; }

	void generateDot()
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> distrib(-1.0, 1.0);

		lux::Vector r(distrib(gen), distrib(gen), distrib(gen));
		r.normalize();
		float q = pow(abs(fspn1.eval(r)), clump);
		lux::Vector r2 = r * q;
		lux::Vector p2 = p + r2 * pScale;
		lux::Vector D = lux::Vector(fspn2.eval(p2), fspn2.eval(p2 + offset), fspn2.eval(p2 - offset));

		generatedPosition = p2 + D * dScale;
	}
};

class Wisp : public Grid
{
private:
	int numberOfDots;
public:
	Wisp() {}
	Wisp(lux::Vector l, int nx, unsigned int ny, unsigned int nz, double d, int n) :
		Grid(l, nx, ny, nz, d), numberOfDots(n) {}

	void stampWispDot(const lux::Vector& p, const float& d)
	{
		if (!withinGrid(p))
			return;

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

		gridData[getIndex(i, j, k)] += d * (1 - wi) * (1 - wj) * (1 - wk);
		gridData[getIndex(i + 1, j, k)] += d * (wi) * (1 - wj) * (1 - wk);
		gridData[getIndex(i, j + 1, k)] += d * (1 - wi) * (wj) * (1 - wk);
		gridData[getIndex(i, j, k + 1)] += d * (1 - wi) * (1 - wj) * (wk);
		gridData[getIndex(i + 1, j + 1, k)] += d * (wi) * (wj) * (1 - wk);
		gridData[getIndex(i + 1, j, k + 1)] += d * (wi) * (1 - wj) * (wk);
		gridData[getIndex(i, j + 1, k + 1)] += d * (1 - wi) * (wj) * (wk);
		gridData[getIndex(i + 1, j + 1, k + 1)] += d * (wi) * (wj) * (wk);
	}

	void stampWisp()
	{
		gridData.resize(Nx * Ny * Nz, 0);

		for (int i = 0; i < numberOfDots; ++i)
		{
			WispDot dot;
			dot.generateDot();
			stampWispDot(dot.getPos(), dot.getDensity());
		}

		std::cout << "Wisp Dots Stamped\n";
	}
};
