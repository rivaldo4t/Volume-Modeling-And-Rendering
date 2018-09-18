#include "RayMarch.hpp"

lux::Color marchRays(std::shared_ptr<Camera>& camera, double x, double y, lux::SFields& scalarFields)
{
	double u = (-1 + 2 * x) * camera->htanFOV();
	double v = (-1 + 2 * y) * camera->vtanFOV();

	lux::Vector q_ij = (u * camera->right()) + (v * camera->up());
	lux::Vector n_ij = (q_ij + camera->view()).unitvector();
	
	lux::Color c(1.0, 0.0, 0.0, 1.0);
	lux::Color L(0.0, 0.0, 0.0, 1.0);

	double sNear = 0.5, sFar = 5.0;
	double T = 1;
	double delta_s = 0.01;
	double delta_T;
	double kappa = 50;
	double s = sNear;

	lux::Vector X = camera->eye() + sNear * n_ij;
	auto field = scalarFields[0].get();

	while (s <= sFar)
	{
		X += delta_s * n_ij;
		double d = field->eval(X);
		if (d > 0)
		{
			delta_T = exp(-kappa * delta_s * d);
			L += c * ((1 - delta_T) * T);
			T *= delta_T;
		}
		s += delta_s;
	}

	return L;
}