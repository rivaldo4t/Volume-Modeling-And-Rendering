#include "RayMarch.hpp"
#include<iostream>
#include <chrono>
#include <ctime>

#include <ImfRgbaFile.h>
namespace IMF = OPENEXR_IMF_NAMESPACE;

#define DSM_GRID

void roundTable(lux::Vector& eye, lux::Vector& view, lux::Vector& up, double stepDegrees)
{
	//eye = lux::Vector(0, 0, 2) *cos(stepDegrees) + lux::Vector(2, 0, 0) * sin(stepDegrees);
	float a = 0.5, b = 1.4;
	//a = 0.2; b = 0.7;
	
	if (stepDegrees == 0.0)
	{
		eye = lux::Vector(0, a, b);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(1.0, 0.0, 0.0), atan(a / b));
	}
	
	else if (stepDegrees == 90.0)
	{
		eye = lux::Vector(b, a, 0);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(0.0, 0.0, -1.0), atan(a / b));
	}

	else if (stepDegrees == 180.0)
	{
		eye = lux::Vector(0, a, -b);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(-1.0, 0.0, 0.0), atan(a / b));
	}

	else if (stepDegrees == 270.0)
	{
		eye = lux::Vector(-b, a, 0);
		view = lux::Vector(0, 0, 0) - eye;
		up = lux::Vector(0, 1, 0).rodriguesRotation(lux::Vector(0.0, 0.0, 1.0), atan(a / b));
	}
}

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, const lux::SField& sfield, const lux::CField& cfield)
{
	const int num_frames = 1;
	const double delta_rot = 360 / num_frames * M_PI / 180.0;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::vector<IMF::Rgba> exr;
	lux::Vector eye, view, up;

	std::cout << "\t\t    ...\n\t         Keep Calm\n\t\t    and\n\t    Let The Rays March\n\t\t    ...\n\n";
	for (int k = 0; k < num_frames; k++)
	{
		start = std::chrono::system_clock::now();
		exr.clear();
		exr.resize(img_h * img_w);
		roundTable(eye, view, up, k * delta_rot);
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

		std::string fileName = "output/Render_3/frame_" + std::to_string(k + 1) + ".exr";
		IMF::RgbaOutputFile file(fileName.c_str(), img_w, img_h, IMF::WRITE_RGBA);
		file.setFrameBuffer(const_cast<IMF::Rgba*>(exr.data()), 1, img_w);
		file.writePixels(img_h);
		std::cout << "Frame Written : " << k + 1 << "\n";

		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n\t\t    ...\n\n";
	}

	std::cout << "Exit Render Loop\n";
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
	double kappa = 50;
	double s = sNear;

	lux::Vector X = pos + sNear * dir;
	lux::Vector lightPos(0.0, 0.8, 0.0);

	while (s <= sFar)
	{
		X += delta_s * dir;
		double d = density->eval(X);
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

double marchRaysDSM(lux::Vector pos, lux::Vector lightPos, const lux::SField& density)
{
	lux::Vector nL = lightPos - pos;
	double sFar = nL.magnitude(), sNear = 0, s = 0;
	nL.normalize();
	lux::Vector X = pos + sNear * nL;
	double delta_s = 0.01;
	double dsm = 0.0;
	double kappa = 50;

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

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, std::shared_ptr<Grid>& g, std::vector<std::shared_ptr<Light>>& lights)
{
	const int num_frames = 4;
	const double delta_rot = 360.f / num_frames;// *M_PI / 180.f;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::vector<IMF::Rgba> exr;
	lux::Vector eye, view, up;
	NoiseParams param;

#if 0
	// mountains and valleys
	param.octaves = 5;
	param.freq = 1.0f;
	param.fJump = 1.6f;
	param.wedgeSpecific = 1.5f;

	lux::SField tn = std::make_shared<lux::TerrainNoise>(param, 2.2f, 0.4f, 1.2f, 1.2f);
	lux::SField pl = std::make_shared<lux::SFPlane>(lux::Vector(0.0, 0.0, 0.0), lux::Vector(0.0, -1.0, 0.0));
	lux::SField wf = std::make_shared<lux::WarpField>(pl, tn);

	lux::SField box = std::make_shared<lux::SFBox>(0.6);
	wf = std::make_shared<lux::SFIntersect>(wf, box);
	
	//cauldron
	param.octaves = 4;
	param.freq = 1.8f;
	param.fJump = 1.1f;
	param.wedgeSpecific = 1.5f;
	lux::SField cauldron = std::make_shared<lux::SFSphere>(lux::Vector(), 0.5);
	cauldron = std::make_shared<lux::SFScale>(cauldron, 0.2);
	cauldron = std::make_shared<lux::PyroclasticField>(cauldron, param);
	cauldron = std::make_shared<lux::SFTranslate>(cauldron, lux::Vector(-0.2, 0.25, 0.7));
	wf = std::make_shared<lux::SFCutout>(wf, cauldron);

	//cave
	lux::SField cave = std::make_shared<lux::SFSphere>(lux::Vector(), 0.5);
	cave = std::make_shared<lux::SFScale>(cave, 0.08);
	cave = std::make_shared<lux::PyroclasticField>(cave, param);
	cave = std::make_shared<lux::SFTranslate>(cave, lux::Vector(0.48, 0.1, 0.82));
	wf = std::make_shared<lux::SFCutout>(wf, cave);

	//landbridge
	param.octaves = 3;
	param.freq = 1.0f;
	param.fJump = 1.8f;
	param.wedgeSpecific = 1.5f;
	lux::SField torus = std::make_shared<lux::SFTorus>(0.8, 0.17, lux::Vector(0.0, 0.0, 0.0), lux::Vector(0.0, 0.0, 1.0));
	torus = std::make_shared<lux::PyroclasticField>(torus, param);
	lux::SField plane = std::make_shared<lux::SFPlane>(lux::Vector(), lux::Vector(0.0, 1.0, 0.0));
	lux::SField landBridge = std::make_shared<lux::SFIntersect>(torus, plane);
	plane = std::make_shared<lux::SFPlane>(lux::Vector(), lux::Vector(0.3, 1.0, 0.0));
	landBridge = std::make_shared<lux::SFIntersect>(torus, plane);
	landBridge = std::make_shared<lux::SFScale>(landBridge, 0.25);
	landBridge = std::make_shared<lux::SFTranslate>(landBridge, lux::Vector(0.1, -0.1, 0.4));
	wf = std::make_shared<lux::SFUnion>(wf, landBridge);

	//g = std::make_shared<Grid>(lux::Vector(-1, -1, -1), 50, 50, 50, 0.04);
	g = std::make_shared<Grid>(lux::Vector(-1, -1, -1), 500, 500, 500, 0.004);
	g->stamp(wf);
	g->writeGrid("D:/temp/vol/exceptSmokeAndMonument.bin");

	std::shared_ptr<Light> key = std::make_shared<Light>(lux::Vector(0.0, 1.0, 0.0),
	//	lux::Vector(-1, -1, -1), 25, 25, 25, 0.08);
	lux::Vector(-1, -1, -1), 250, 250, 250, 0.008);
	//lux::Vector(-1, -1, -1), 500, 500, 500, 0.004);
	key->setColor(lux::Color(0.5, 0.5, 0.5, 1.0));
	lights = { key };
#endif

#ifdef DSM_GRID
	for (auto l : lights)
		l->computeDSM(g);
#endif
	//

	std::cout << "\t\t  ...\n\t       Keep Calm\n\t\t  and\n\t  Let The Rays March\n\t\t  ...\n\n";
	for (int k = 0; k < num_frames; k++)
	{
		start = std::chrono::system_clock::now();
		exr.clear();
		exr.resize(img_h * img_w);
		roundTable(eye, view, up, k * delta_rot);
		camera->setEyeViewUp(eye, view, up);
		//param.updateParams();
#if 1
		/*param.octaves = 4;
		param.freq = 0.8f; // 1.2
		param.fJump = 1.5f; // 1.6
		param.wedgeSpecific = 1.5f;*/

		/*g = std::make_shared<Grid>(lux::Vector(-1, -1, -1), 500, 500, 500, 0.004);
		lux::SField sp = std::make_shared<lux::SFSphere>(lux::Vector(), 0.5);
		lux::SField s = std::make_shared<lux::PyroclasticField>(sp, param);
		g->stamp(s);*/

		/*auto wisp1 = std::make_shared<Wisp>(lux::Vector(-1, -1, -1), 500, 500, 500, 0.004, 50000);
		wisp1->stampWisp(param);
		auto g1 = std::make_shared<Grid>(lux::Vector(-1, -1, -1), 500, 500, 500, 0.004);
		g1 = wisp1;
		g1 = std::make_shared<GridTranslate>(g1, lux::Vector(-0.5, 0.0, 0.0));
		
		auto wisp2 = std::make_shared<Wisp>(lux::Vector(-1, -1, -1), 500, 500, 500, 0.004, 50000);
		wisp2->stampWisp(param);
		auto g2 = std::make_shared<Grid>(lux::Vector(-1, -1, -1), 500, 500, 500, 0.004);
		g2 = wisp1;
		g2 = std::make_shared<GridTranslate>(g2, lux::Vector(0.5, 0.0, 0.0));

		std::shared_ptr<Grid> g3 = std::make_shared<GridUnion>(g1, g2);
		g = g3;*/

#if 0
		param.updateParams();
		auto noise = std::make_shared<StampedNoise>(lux::Vector(0, 0, 0), 0.1f, lux::Vector(-1, -1, -1), 100, 100, 100, 0.02);
		noise->computeNoise(param);
		g = noise;

		param.updateParams();
		auto noise2 = std::make_shared<StampedNoise>(lux::Vector(-0.1, 0.1, -0.1), 0.1f, lux::Vector(-1, -1, -1), 100, 100, 100, 0.02);
		noise2->computeNoise(param);
		std::shared_ptr<Grid> g2 = noise2;
		g = std::make_shared<GridUnion>(g, g2);

		param.updateParams();
		auto noise3 = std::make_shared<StampedNoise>(lux::Vector(0, 0.1, -0.1), 0.1f, lux::Vector(-1, -1, -1), 100, 100, 100, 0.02);
		noise3->computeNoise(param);
		g2 = noise3;
		g = std::make_shared<GridUnion>(g, g2);

		param.updateParams();
		auto noise4 = std::make_shared<StampedNoise>(lux::Vector(-0.1, 0.1, 0), 0.1f, lux::Vector(-1, -1, -1), 100, 100, 100, 0.02);
		noise4->computeNoise(param);
		g2 = noise4;
		g = std::make_shared<GridUnion>(g, g2);

		param.updateParams();
		auto noise5 = std::make_shared<StampedNoise>(lux::Vector(-0.2, 0.2, -0.2), 0.1f, lux::Vector(-1, -1, -1), 100, 100, 100, 0.02);
		noise5->computeNoise(param);
		g2 = noise5;
		g = std::make_shared<GridUnion>(g, g2);

		param.updateParams();
		auto noise6 = std::make_shared<StampedNoise>(lux::Vector(0.1, 0.2, -0.2), 0.1f, lux::Vector(-1, -1, -1), 100, 100, 100, 0.02);
		noise6->computeNoise(param);
		g2 = noise6;
		g = std::make_shared<GridUnion>(g, g2);

		param.updateParams();
		auto noise7 = std::make_shared<StampedNoise>(lux::Vector(-.2, 0.2, 0.1), 0.1f, lux::Vector(-1, -1, -1), 100, 100, 100, 0.02);
		noise7->computeNoise(param);
		g2 = noise7;
		g = std::make_shared<GridUnion>(g, g2);
#endif
#endif

//#ifdef DSM_GRID
//		for (auto l : lights)
//			l->computeDSM(g);
//#endif

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

				lux::Color L = marchRays(camera->eye(), n_ij, g, lights);
				exr[(img_h - 1 - j) * img_w + i] = IMF::Rgba(half(L[0]), half(L[1]), half(L[2]), half(L[3]));
			}
		}

		std::cout << "100%|\n";

		std::string fileName = "output/Render_3/frame_" + std::to_string(k + 1) + ".exr";
		IMF::RgbaOutputFile file(fileName.c_str(), img_w, img_h, IMF::WRITE_RGBA);
		file.setFrameBuffer(const_cast<IMF::Rgba*>(exr.data()), 1, img_w);
		file.writePixels(img_h);
		std::cout << "Frame Written : " << k + 1 << "\n";

		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n\t\t    ...\n\n";
	}

	std::cout << "Exit Render Loop\n";
}

lux::Color marchRays(lux::Vector pos, lux::Vector dir, const std::shared_ptr<Grid>& g, const std::vector<std::shared_ptr<Light>>& lights)
{
	lux::Color L(0.0, 0.0, 0.0, 0.0);
	lux::Color white(1.0, 1.0, 1.0, 1.0);

	double sNear = 0.2, sFar = 6.0;
	double T = 1;
	double delta_s = 0.01;
	double delta_T;
	double kappa = 800; //0.1 for wisp
	double s = sNear;

	lux::Vector X = pos + sNear * dir;

	while (s <= sFar)
	{
		X += delta_s * dir;
		double d = g->eval(X);
		
		lux::Color c;
		for (auto l : lights)
		{
			double dsm;
#ifdef DSM_GRID
			dsm = l->eval(X);
#else
			dsm = marchToLight(X, l->getPosition(), g);
#endif
			c += white * l->getColor() * exp(-kappa * dsm);
		}
		if (lights.size() == 0)
			c = white;

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

double marchToLight(lux::Vector pos, lux::Vector lightPos, const std::shared_ptr<Grid>& g)
{
	lux::Vector nL = lightPos - pos;
	double sFar = nL.magnitude(), sNear = 0, s = 0;
	nL.normalize();
	lux::Vector X = pos + sNear * nL;
	double delta_s = 0.01;
	double kappa = 50;
	double dsm = 0.0;

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
	return dsm;
}