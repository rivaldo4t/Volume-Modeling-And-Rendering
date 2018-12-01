#include "RayMarch.hpp"
#include<iostream>
#include <chrono>
#include <ctime>
#include <ImfRgbaFile.h>
namespace IMF = OPENEXR_IMF_NAMESPACE;

#define DSM_GRID

void roundTable(std::shared_ptr<Camera> camera, float stepDegrees)
{
	lux::Vector eye, view, up;
#if 1
	// circle around origin
	eye = lux::Vector(0, 0, 2) * cos(stepDegrees * M_PI / 180.f) + lux::Vector(2, 0, 0) * sin(stepDegrees * M_PI / 180.f);
	view = lux::Vector(0, 0, 0) - eye;
	up = lux::Vector(0, 1, 0);
#else
	// views at degrees steps
	float a = 0.5, b = 1.4;
	float theta = -atan(a / b);
	
	if (stepDegrees == 0.0)
	{
		eye = lux::Vector(0, a, b);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(1.0, 0.0, 0.0), theta);
	}
	
	else if (stepDegrees == 90.0)
	{
		eye = lux::Vector(b, a, 0);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(0.0, 0.0, -1.0), theta);
	}

	else if (stepDegrees == 180.0)
	{
		eye = lux::Vector(0, a, -b);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(-1.0, 0.0, 0.0), theta);
	}

	else if (stepDegrees == 270.0)
	{
		eye = lux::Vector(-b, a, 0);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(0.0, 0.0, 1.0), theta);
	}
#endif

	camera->setEyeViewUp(eye, view, up);
}

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, const lux::SField& sfield, const lux::CField& cfield)
{
	const int num_frames = 1;
	const float stepDegrees = 360 / num_frames;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::vector<IMF::Rgba> exr;
	lux::Vector eye, view, up;

	std::cout << "\t\t    ...\n\t         Keep Calm\n\t\t    and\n\t    Let The Rays March\n\t\t    ...\n\n";
	for (int k = 0; k < num_frames; k++)
	{
		start = std::chrono::system_clock::now();
		exr.clear();
		exr.resize(img_h * img_w);
		roundTable(camera, k * stepDegrees);

		std::cout << "|0%|==|==|==|==|==|==|==|==|==|==|==|100%|\n|0%|";
#pragma omp parallel for
		for (int j = 0; j < img_h; ++j)
		{
			if ((j) % (img_h / 10) == 0)
				std::cout << "==|";
			for (int i = 0; i < img_w; ++i)
			{
				float x = float(i) / (img_w - 1);
				float y = float(j) / (img_h - 1);
				float u = (-1 + 2 * x) * camera->htanFOV();
				float v = (-1 + 2 * y) * camera->vtanFOV();
				lux::Vector q_ij = (u * camera->right()) + (v * camera->up());
				lux::Vector n_ij = (q_ij + camera->view()).unitvector();

				lux::Color L = marchRays(camera->eye(), n_ij, sfield, cfield);
				exr[(img_h - 1 - j) * img_w + i] = IMF::Rgba(half(L[0]), half(L[1]), half(L[2]), half(L[3]));
			}
		}
		std::cout << "100%|\n";

		std::string fileName = "output/Render_4/frame_" + std::to_string(k + 1) + ".exr";
		IMF::RgbaOutputFile file(fileName.c_str(), img_w, img_h, IMF::WRITE_RGBA);
		file.setFrameBuffer(const_cast<IMF::Rgba*>(exr.data()), 1, img_w);
		file.writePixels(img_h);
		std::cout << "Frame Written : " << k + 1 << "\n";

		end = std::chrono::system_clock::now();
		std::chrono::duration<float> elapsed_seconds = end - start;
		std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n\t\t    ...\n\n";
	}

	std::cout << "Exit Render Loop\n";
}

lux::Color marchRays(lux::Vector pos, lux::Vector dir, const lux::SField& density, const lux::CField& color)
{	
	lux::Color L(0.0, 0.0, 0.0, 1.0);
	lux::Color white(0.8, 0.8, 0.8, 1.0);
	lux::Color red(1.0, 0.2, 0.2, 1.0);

	float sNear = 0.5, sFar = 3.0;
	float T = 1;
	float delta_s = 0.01;
	float delta_T;
	float kappa = 50;
	float s = sNear;

	lux::Vector X = pos + sNear * dir;
	lux::Vector lightPos(0.0, 0.8, 0.0);

	while (s <= sFar)
	{
		X += delta_s * dir;
		float d = density->eval(X);
		lux::Color c = white;// color->eval(X);
		//c = c.isZero() ? white : c;
		c *= marchRaysDSM(X, lightPos, density);

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

float marchRaysDSM(lux::Vector pos, lux::Vector lightPos, const lux::SField& density)
{
	lux::Vector nL = lightPos - pos;
	float sFar = nL.magnitude(), sNear = 0, s = 0;
	nL.normalize();
	lux::Vector X = pos + sNear * nL;
	float delta_s = 0.01;
	float dsm = 0.0;
	float kappa = 50;

	if (density->eval(X) > 0)
	{
		while (s <= sFar)
		{
			X += delta_s * nL;
			float d = density->eval(X);
			if (d > 0)
				dsm += d * delta_s;
			s += delta_s;
		}
	}
	return exp(-kappa * dsm);
}

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, lux::SField& g, std::vector<std::shared_ptr<lux::Light>>& lights)
{
	std::cout << "\t\t  ...\n\t       Keep Calm\n\t\t  and\n\t  Let The Rays March\n\t\t  ...\n\n";

	const int num_frames = 100;
	const float stepDegrees = 360.f / num_frames;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::vector<IMF::Rgba> exr;
	NoiseParams param;
	param.updateParams();

	lux::SField sp = std::make_shared<lux::SFSphere>(lux::Vector(), 0.005);
	lux::SField sp2 = std::make_shared<lux::SFSphere>(lux::Vector(0.0, 0.5, 0.0), 0.2);
	lux::SField to1 = std::make_shared<lux::SFTorus>(0.4, 0.04, lux::Vector(), lux::Vector(-1, 1, 1));
	to1 = std::make_shared<lux::SFTranslate>(to1, lux::Vector(-0.8, -0.3, 0.0));
	lux::SField to2 = std::make_shared<lux::SFTorus>(0.4, 0.04, lux::Vector(), lux::Vector(0, 1, 0));
	to2 = std::make_shared<lux::SFTranslate>(to2, lux::Vector(0.6, 0.7, 0.0));
	lux::SField pl = std::make_shared<lux::SFPlane>(lux::Vector(), lux::Vector(0.0, -0.2, 0.0));
	lux::SField bx = std::make_shared<lux::SFBox>(1);
	lux::SField g2 = std::make_shared<lux::SFUnion>(to1, to2);
	g = sp;

	float l = 500;
	float o = 1;
	std::shared_ptr<lux::Light> key = std::make_shared<lux::Light>(lux::Vector(0.0, 1.0, 1.0),
		lux::Vector(-o, -o, -o), 3 * l/5, 3 * l/5, 3 * l/5, 2 * o / (3 * l/5));
	key->setColor(lux::Color(0.9, 0.9, 0.9, 1.0));
	lights = { key };

	float dt = 0.0;
	for (int k = 0; k < num_frames; k++)
	{
		start = std::chrono::system_clock::now();
		exr.clear();
		exr.resize(img_h * img_w);
		//roundTable(camera, k * stepDegrees);
		
		auto w = std::make_shared<lux::Wisp>(lux::Vector(-1, -1, -1), l, l, l, 2.0 / l, 500000);
		w->stampWisp(param, lux::Vector(float(k) / 24.f - 0.5, 0.0, 0.0), 0.2 + float(k + 1) / 24.f, float(24-k)/6.f);
		auto g2 = std::make_shared<lux::SFUnion>(g, w);
		auto g3 = std::make_shared<lux::ScalarGrid>(lux::Vector(-1, -1, -1), l, l, l, 2.0 / l);
		g3->stamp(g2);
		g = g3;

		param.updateParams();
		dt += 0.0;

#ifdef DSM_GRID
		for (auto l : lights)
			l->computeDSM(g);
#endif

		std::cout << "\n|0%|==|==|==|==|==|==|==|==|==|==|100%|\n|0%|";

#pragma omp parallel for
		for (int j = 0; j < img_h; ++j)
		{
			if ((j) % (img_h / 10) == 0)
				std::cout << "==|";
			for (int i = 0; i < img_w; ++i)
			{
				float x = float(i) / (img_w - 1);
				float y = float(j) / (img_h - 1);
				float u = (-1 + 2 * x) * camera->htanFOV();
				float v = (-1 + 2 * y) * camera->vtanFOV();
				lux::Vector q_ij = (u * camera->right()) + (v * camera->up());
				lux::Vector n_ij = (q_ij + camera->view()).unitvector();

				lux::Color L = marchRays(camera->eye(), n_ij, g, lights);
				exr[(img_h - 1 - j) * img_w + i] = IMF::Rgba(half(L[0]), half(L[1]), half(L[2]), half(L[3]));
			}
		}

		std::cout << "100%|\n";

		std::string fileName = "output/Render_5/frame_" + std::to_string(k + 1) + ".exr";
		IMF::RgbaOutputFile file(fileName.c_str(), img_w, img_h, IMF::WRITE_RGBA);
		file.setFrameBuffer(const_cast<IMF::Rgba*>(exr.data()), 1, img_w);
		file.writePixels(img_h);
		std::cout << "Frame Written : " << k + 1 << "\n";

		end = std::chrono::system_clock::now();
		std::chrono::duration<float> elapsed_seconds = end - start;
		std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n\t\t    ...\n\n";
	}

	std::cout << "Exit Render Loop\n";
}

lux::Color marchRays(lux::Vector pos, lux::Vector dir, const lux::SField& g, const std::vector<std::shared_ptr<lux::Light>>& lights)
{
	lux::Color L(0.0, 0.0, 0.0, 0.0);
	lux::Color white(1.0, 1.0, 1.0, 1.0);

	float sNear = 0.2, sFar = 6.0;
	float T = 1;
	float delta_s = 0.01;
	float delta_T;
	float kappa = 1; //0.1 for wisp // 800 otherwise
	float s = sNear;

	lux::Vector X = pos + sNear * dir;

	while (s <= sFar)
	{
		X += delta_s * dir;
		float d = g->eval(X);

		if (d > 0)
		{
			lux::Color c;
			for (auto l : lights)
			{
				float dsm;
#ifdef DSM_GRID
				dsm = l->eval(X);
#else
				dsm = marchToLight(X, l->getPosition(), g);
#endif
				c += white * l->getColor() * exp(-kappa * dsm);
			}
			if (lights.size() == 0)
				c = white;

			delta_T = exp(-kappa * delta_s * d);
			L += c * ((1 - delta_T) * T);
			T *= delta_T;
		}

		s += delta_s;
	}

	return L;
}

float marchToLight(lux::Vector pos, lux::Vector lightPos, const lux::SField& g)
{
	lux::Vector nL = lightPos - pos;
	float sFar = nL.magnitude(), sNear = 0, s = 0;
	nL.normalize();
	lux::Vector X = pos + sNear * nL;
	float delta_s = 0.01;
	float dsm = 0.0;

	if (g->eval(X) > 0)
	{
		while (s <= sFar)
		{
			X += delta_s * nL;
			float d = g->eval(X);
			if (d > 0)
				dsm += d * delta_s;
			s += delta_s;
		}
	}
	return dsm;
}