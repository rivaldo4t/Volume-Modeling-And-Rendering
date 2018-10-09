#include "RayMarch.hpp"
#include<iostream>
#include <chrono>
#include <ctime>

#include <ImfRgbaFile.h>
namespace IMF = OPENEXR_IMF_NAMESPACE;

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, const lux::SField& sfield, const lux::CField& cfield)
{
	const int num_frames = 120 / 1;
	const double delta_rot = 360 / num_frames * M_PI / 180.0;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::vector<IMF::Rgba> exr;
	lux::Vector eye, view, up;

	std::cout << "\t\t    ...\n\t         Keep Calm\n\t\t    and\n\t    Let The Rays March\n\t\t    ...\n\n";
	for (int k = 0; k <= num_frames; k++)
	{
		start = std::chrono::system_clock::now();
		exr.clear();
		exr.resize(img_h * img_w);
		eye = lux::Vector(0, 0.2, 2) * cos(k * delta_rot) + lux::Vector(2, 0.2, 0) * sin(k * delta_rot);
		view = lux::Vector(0, 0.2, 0) - eye;
		up = lux::Vector(0, 1, 0);
		camera->setEyeViewUp(eye, view, up);

		std::cout << "|0%|==|==|==|==|==|==|==|==|==|==|==|100%|\n|0%|";
#pragma omp parallel for
		for (int j = 0; j < img_h; ++j)
		{
			if ((j) % (img_h / 10) == 0)
				std::cout << "==|";
			for (int i = 0; i < img_w; ++i)
			{
				double x = double(i) / (img_w - 1);
				double y = double(j) / (img_h - 1);
				double u = (-1 + 2 * x) * camera->htanFOV();
				double v = (-1 + 2 * y) * camera->vtanFOV();
				lux::Vector q_ij = (u * camera->right()) + (v * camera->up());
				lux::Vector n_ij = (q_ij + camera->view()).unitvector();

				lux::Color L = marchRays(camera->eye(), n_ij, sfield, cfield);
				exr[(img_h - 1 - j) * img_w + i] = IMF::Rgba(half(L[0]), half(L[1]), half(L[2]), half(L[3]));
			}
		}
		std::cout << "100%|\n";

		std::string fileName = "output/Render_1/frame_" + std::to_string(k + 1) + ".exr";
		IMF::RgbaOutputFile file(fileName.c_str(), img_w, img_h, IMF::WRITE_RGBA);
		file.setFrameBuffer(const_cast<IMF::Rgba*>(exr.data()), 1, img_w);
		file.writePixels(img_h);
		std::cout << "Frame Written : " << k + 1 << "\n";

		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n\t\t    ...\n\n";
	}
}

lux::Color marchRays(lux::Vector pos, lux::Vector dir, const lux::SField& density, const lux::CField& color)
{	
	lux::Color L(0.0, 0.0, 0.0, 1.0);
	lux::Color white(0.8, 0.8, 0.8, 1.0);
	lux::Color red(1.0, 0.2, 0.2, 1.0);

	double sNear = 0.5, sFar = 3.0;
	double T = 1;
	double delta_s = 0.01;
	double delta_T;
	double kappa = 10;
	double s = sNear;

	lux::Vector X = pos + sNear * dir;
	lux::Vector lightPos(0.0, 0.8, 0.0);

	while (s <= sFar)
	{
		X += delta_s * dir;
		double d = density->eval(X);
		lux::Color c = color->eval(X);
		c = c.isZero() ? white : c;
		//c *= marchRaysDSM(X, lightPos, density);

		// explicit color to objects
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

double marchRaysDSM(lux::Vector pos, lux::Vector lightPos, const lux::SField& density)
{
	lux::Vector nL = lightPos - pos;
	double sFar = nL.magnitude(), sNear = 0, s = 0;
	nL.normalize();
	lux::Vector X = pos + sNear * nL;
	double delta_s = 0.01;
	double dsm = 0.0;
	double kappa = 10;

	if (density->eval(X) > 0)
	{
		while (s <= sFar)
		{
			X += delta_s * nL;
			double d = density->eval(X);
			if (d > 0)
				dsm += d * delta_s;
			s += delta_s;
		}
	}
	return exp(-kappa * dsm);
}

void render2(const int img_w, const int img_h, std::shared_ptr<Camera> camera, const std::shared_ptr<Grid>& g, const std::vector<std::shared_ptr<Light>>& lights)
{
	const int num_frames = 120 / 1;
	const double delta_rot = 360 / num_frames * M_PI / 180.0;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::vector<IMF::Rgba> exr;
	lux::Vector eye, view, up;

	std::cout << "\t\t  ...\n\t       Keep Calm\n\t\t  and\n\t  Let The Rays March\n\t\t  ...\n\n";

	for (int k = 0; k <= num_frames; k++)
	{
		start = std::chrono::system_clock::now();
		exr.clear();
		exr.resize(img_h * img_w);
		
		eye = lux::Vector(0, 0.2, 2) * cos(k * delta_rot) + lux::Vector(2, 0.2, 0) * sin(k * delta_rot);
		view = lux::Vector(0, 0.2, 0) - eye;
		up = lux::Vector(0, 1, 0);
		eye = lux::Vector(0, 0, 2) * cos(k * delta_rot) + lux::Vector(2, 0, 0) * sin(k * delta_rot);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0);
		camera->setEyeViewUp(eye, view, up);

		std::cout << "|0%|==|==|==|==|==|==|==|==|==|==|100%|\n|0%|";

#pragma omp parallel for
		for (int j = 0; j < img_h; ++j)
		{
			if ((j) % (img_h / 10) == 0)
				std::cout << "==|";
			for (int i = 0; i < img_w; ++i)
			{
				double x = double(i) / (img_w - 1);
				double y = double(j) / (img_h - 1);
				double u = (-1 + 2 * x) * camera->htanFOV();
				double v = (-1 + 2 * y) * camera->vtanFOV();
				lux::Vector q_ij = (u * camera->right()) + (v * camera->up());
				lux::Vector n_ij = (q_ij + camera->view()).unitvector();

				lux::Color L = marchRays2(camera->eye(), n_ij, g, lights);
				exr[(img_h - 1 - j) * img_w + i] = IMF::Rgba(half(L[0]), half(L[1]), half(L[2]), half(L[3]));
			}
		}

		std::cout << "100%|\n";

		std::string fileName = "output/Render_1/frame_" + std::to_string(k + 1) + ".exr";
		IMF::RgbaOutputFile file(fileName.c_str(), img_w, img_h, IMF::WRITE_RGBA);
		file.setFrameBuffer(const_cast<IMF::Rgba*>(exr.data()), 1, img_w);
		file.writePixels(img_h);

		std::cout << "Frame Written : " << k + 1 << "\n";
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n\t\t    ...\n\n";
	}
}

lux::Color marchRays2(lux::Vector pos, lux::Vector dir, const std::shared_ptr<Grid>& g, const std::vector<std::shared_ptr<Light>>& lights)
{
	// should alpha be initialized to 1?
	lux::Color L(0.0, 0.0, 0.0, 0.0);
	lux::Color white(0.8, 0.8, 0.8, 1.0);

	double sNear = 0.2, sFar = 4.0;
	double T = 1;
	double delta_s = 0.01;
	double delta_T;
	double kappa = 10;
	double s = sNear;

	lux::Vector X = pos + sNear * dir;

	while (s <= sFar)
	{
		X += delta_s * dir;
		double d = g->eval(X);
		
		lux::Color c(0.0, 0.0, 0.0, 0.0);
		for (auto l : lights)
			c += white * l->getColor() * exp(-kappa * l->eval(X));
		c = c.isZero() ? white : c;

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

double marchRaysDSM2(lux::Vector pos, lux::Vector lightPos, const std::shared_ptr<Grid>& g)
{
	lux::Vector nL = lightPos - pos;
	double sFar = nL.magnitude(), sNear = 0, s = 0;
	nL.normalize();
	lux::Vector X = pos + sNear * nL;
	double delta_s = 0.01;
	double dsm = 0.0;
	double kappa = 10;

	if (g->eval(X) > 0)
	{
		while (s <= sFar)
		{
			X += delta_s * nL;
			double d = g->eval(X);
			if (d > 0)
				dsm += d * delta_s;
			s += delta_s;
		}
	}
	return exp(-kappa * dsm);
}