#include "RayMarch.hpp"
#include<iostream>
#include <chrono>
#include <ctime>

#include <ImfRgbaFile.h>
namespace IMF = OPENEXR_IMF_NAMESPACE;

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, lux::SField sfield, lux::CField cfield)
{
	const int num_frames = 120 / 1;
	const double delta_rot = 360 / num_frames;
	const double delta_rotr = delta_rot * M_PI / 180.0;

	std::vector<IMF::Rgba> exr;

	lux::Vector eye, view, up;
	std::chrono::time_point<std::chrono::system_clock> start, end;

	std::cout << "\t\t    ...\n\t         Keep Calm\n\t\t    and\n\t    Let The Rays March\n\t\t    ...\n\n";
	for (int k = 0; k <= num_frames; k++)
	{
		start = std::chrono::system_clock::now();

		exr.clear();
		exr.reserve(img_h * img_w);

		eye = lux::Vector(0, 0, 1) * cos(k * delta_rotr) + lux::Vector(1.4, 0, 0) * sin(k * delta_rotr);
		view = lux::Vector(0, 0.2, -0.4) - eye;
		up = lux::Vector(0, 1, 0);
		camera->setEyeViewUp(eye, view, up);

		std::cout << "|0%|==|==|==|==|==|==|==|==|==|==|==|100%|\n|0%|";
		for (int j = img_h; j >= 0; --j)
		{
			if ((img_h - j) % (img_h / 10) == 0)
				std::cout << "==|";
			for (int i = 0; i < img_w; ++i)
			{
				lux::Color L = marchRays(camera, double(i) / (img_w - 1), double(j) / (img_h - 1), sfield, cfield);
				exr.push_back(IMF::Rgba(half(L[0]), half(L[1]), half(L[2]), half(L[3])));
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

lux::Color marchRays(std::shared_ptr<Camera>& camera, double x, double y, lux::SField field, lux::CField color)
{
	double u = (-1 + 2 * x) * camera->htanFOV();
	double v = (-1 + 2 * y) * camera->vtanFOV();

	lux::Vector q_ij = (u * camera->right()) + (v * camera->up());
	lux::Vector n_ij = (q_ij + camera->view()).unitvector();
	
	lux::Color L(0.0, 0.0, 0.0, 1.0);
	lux::Color white(0.8, 0.8, 0.8, 1.0);
	lux::Color red(1.0, 0.2, 0.2, 1.0);

	double sNear = 0.5, sFar = 3.0;
	double T = 1;
	double delta_s = 0.01;
	double delta_T;
	double kappa = 1;
	double s = sNear;

	lux::Vector X = camera->eye() + sNear * n_ij;

	while (s <= sFar)
	{
		X += delta_s * n_ij;
		double d = field->eval(X);
		lux::Color c = color->eval(X);
		c = c.isZero() ? white : c;
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