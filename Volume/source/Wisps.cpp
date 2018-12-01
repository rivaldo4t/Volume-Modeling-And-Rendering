#include "Wisps.hpp"
using namespace lux;

void WispDot::setPscale(float p) { pScale = p; }
void WispDot::setDscale(float d) { dScale = d; }
void WispDot::setDensity(float d) { density = d; }

void Wisp::stampWispDot(const lux::Vector& p, const float& d)
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

void Wisp::stampWisp(NoiseParams& param, lux::Vector guidePos, float psize, float dense)
{
	gridData.resize(Nx * Ny * Nz, 0);

	FSPN f1 = FSPN(param.octaves, param.freq, param.fJump, 2);
	FSPN f2 = FSPN(param.octaves_values[param.octaves_values.size() - 1 - param.o],
		param.freq_values[param.freq_values.size() - 1 - param.fr],
		param.fJump_values[param.fJump_values.size() - 1 - param.fj],
		2);
	float clum = param.wedgeSpecific;

	std::cout << "\n-------------------\n";
	std::cout << "f2 octaves:\t" << f2.octaves << std::endl;
	std::cout << "f2 freq:\t" << f2.freq << std::endl;
	std::cout << "f2 fjump:\t" << f2.fJump << std::endl;
	std::cout << "f2 roughness:\t" << f2.roughness << std::endl;
	std::cout << "-------------------\n\n";

	std::cout << "Stamping Wisp Dots . . . . ";

	WispDot dot(guidePos, f1, f2, clum);
	dot.setPscale(psize);
	dot.setDscale(psize * 0.5);
	dot.setDensity(dense);
	for (int i = 0; i < numberOfDots; ++i)
	{
		dot.generateDot();
		stampWispDot(dot.getPos(), dot.getDensity());
	}

	std::cout << "Done\n";
}

void WispDot::generateDot()
{
	float correlationFact = 0.7f;

	lux::Vector r0(2 * distrib(gen) - 1, 2 * distrib(gen) - 1, 2 * distrib(gen) - 1);
	r0 = correlationFact * guidePos + (1 - correlationFact) * r0;
	lux::Vector r1 = r0.unitvector();
	float q = pow(abs(fspn1.eval(r0)), clump);
	lux::Vector r2 = r1 * q;
	lux::Vector p2 = guidePos + r2 * pScale;
	lux::Vector D = lux::Vector(fspn2.eval(r2), fspn2.eval(r2 + offset), fspn2.eval(r2 - offset));
	generatedPosition = p2 + D * dScale;
}