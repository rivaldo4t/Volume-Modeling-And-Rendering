#include "Fspn.hpp"

void NoiseParams::updateParams()
{
	octaves = octaves_values[o];
	freq = freq_values[fr];
	fJump = fJump_values[fj];
	wedgeSpecific = wedgeSpecific_values[w];

	std::cout << "--------------------\n";
	std::cout << "octaves:\t" << octaves << std::endl;
	std::cout << "freq:\t\t" << freq << std::endl;
	std::cout << "fjump:\t\t" << fJump << std::endl;
	std::cout << "wedgeSpecific:\t\t" << wedgeSpecific << std::endl;
	std::cout << "--------------------\n";

	w++;
	if (w == wedgeSpecific_values.size())
		fj++;
	if (fj == fJump_values.size())
		fr++;
	if (fr == freq_values.size())
		o++;

	w %= wedgeSpecific_values.size();
	fj %= fJump_values.size();
	fr %= freq_values.size();
	o %= octaves_values.size();
}

float FSPN::eval(lux::Vector x) const
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