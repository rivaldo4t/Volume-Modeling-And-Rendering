#pragma once
#include <random>
#include "Grid.hpp"
#include "Fspn.hpp"

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
	void generateDot();
};

class Wisp : public Grid
{
private:
	int numberOfDots;
public:
	Wisp() {}
	Wisp(lux::Vector l, int nx, unsigned int ny, unsigned int nz, double d, int n) :
		Grid(l, nx, ny, nz, d), numberOfDots(n) {}
	void stampWispDot(const lux::Vector& p, const float& d);
	void stampWisp();
};
