#include "Fspn.hpp"

float FSPN::eval(lux::Vector x)
{
	float coeff = roughness == 1.0 ? (1.0 / float(octaves)) : (1 - roughness) / (1 - pow(roughness, octaves));
	float sumPN = 0;
	for (int i = 0; i < octaves; ++i)
	{
		lux::Vector evalAt = (x - translate) * pow(fJump, i) * freq;
		float pn = PN.eval(evalAt);
		sumPN += pow(roughness, i) * pn;
	}
	return coeff * sumPN;
}