#pragma once
#include <random>
#include "Grid.hpp"
#include "FractalSummedPerlinNoise.hpp"

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<double> distrib(0.0, 1.0);

class WispDot
{
private:
	lux::Vector guidePos;
	lux::Vector offset;
	float pScale;
	float dScale;
	float clump;
	float density;
	FSPN fspn1, fspn2;
	lux::Vector generatedPosition;
public:
	WispDot(lux::Vector pos = lux::Vector(), FSPN& f1 = FSPN(), FSPN& f2 = FSPN(), float clum = 2) :
		guidePos(pos), fspn1(f1), fspn2(f2), clump(clum)
	{ 
		offset = lux::Vector(0.1, 0.1, 0.1);
		pScale = 0.5;
		dScale = 2;
		density = 10;
	}

	lux::Vector getPos() { return generatedPosition; }
	float getDensity() { return density; }

	void generateDot()
	{
		float correlationFact = 0.7;

		lux::Vector r0(2 * distrib(gen) - 1, 2 * distrib(gen) - 1, 2 * distrib(gen) - 1);
		r0 = correlationFact * guidePos + (1 - correlationFact) * r0;
		lux::Vector r1 = r0.unitvector();
		float q = pow(abs(fspn1.eval(r0)), clump);
		lux::Vector r2 = r1 * q;
		lux::Vector p2 = guidePos + r2 * pScale;
		lux::Vector D = lux::Vector(fspn2.eval(r2), fspn2.eval(r2 + offset), fspn2.eval(r2 - offset));
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
		
		static int _oct, _freq, _fjump, _clump;
		std::vector<int> oct = { 1, 2, 3, 4, 5 };
		std::vector<double> freq = { 0.5, 1, 2.5, 3, 4.3, 5.6, 10 };
		std::vector<double> fjump = { 1, 1.5, 2 };
		std::vector<double> clu = { 0.5, 1, 1.5, 2, 3 };

		FSPN f1 = FSPN(oct[_oct], freq[_freq], fjump[_fjump], 2.0);
		FSPN f2 = FSPN(oct[oct.size() - 1 - _oct], freq[freq.size() - 1 - _freq], fjump[fjump.size() - 1 - _fjump], 2.0);
		float clum = clu[_clump];
		
		std::cout << "-------------------\n";
		std::cout << "f1 octaves:\t" << f1.octaves << std::endl;
		std::cout << "f1 freq:\t" << f1.freq << std::endl;
		std::cout << "f1 fjump:\t" << f1.fJump << std::endl;
		std::cout << "f1 roughness:\t" << f1.roughness << std::endl;
		std::cout << "f2 octaves:\t" << f2.octaves << std::endl;
		std::cout << "f2 freq:\t" << f2.freq << std::endl;
		std::cout << "f2 fjump:\t" << f2.fJump << std::endl;
		std::cout << "f2 roughness:\t" << f2.roughness << std::endl;
		std::cout << "clump:\t\t" << clum << std::endl;
		std::cout << "-------------------\n";

		_clump++;
		if (_clump == clu.size())
			_fjump++;
		if (_fjump == fjump.size())
			_freq++;
		if (_freq == freq.size())
			_oct++;

		_clump %= clu.size();
		_fjump %= fjump.size();
		_freq %= freq.size();
		_oct %= oct.size();

		WispDot dot(lux::Vector(), f1, f2, clum);
		for (int i = 0; i < numberOfDots; ++i)
		{
			dot.generateDot();
			stampWispDot(dot.getPos(), dot.getDensity());
			if (i % 500000 == 0)
				std::cout << ".";
		}

		std::cout << "\nWisp Dots Stamped\n";
	}
};
