#include <iostream>
#include <memory>
#include <vector>

#include "ScalarField.hpp"
#include "Camera.hpp"
#include "RayMarch.hpp"

#include <ImfRgbaFile.h>
namespace IMF = OPENEXR_IMF_NAMESPACE;

int main()
{
	// camera
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();
	
	// sphere
	std::shared_ptr<lux::ScalarField> sp = std::make_shared<lux::SFSphere>(lux::Vector(-0.2, 0.0, 0.0), 0.3);

	// plane
	std::shared_ptr<lux::ScalarField> p1 = std::make_shared<lux::SFPlane>(lux::Vector(0.0, 0.1, 0.0), lux::Vector(0.0, 1.0, 0.0));
	std::shared_ptr<lux::ScalarField> p2 = std::make_shared<lux::SFPlane>(lux::Vector(0.0, -0.1, 0.0), lux::Vector(0.0, 1.0, 0.0));

	// torus
	std::shared_ptr<lux::ScalarField> r = std::make_shared<lux::SFTorus>(0.2, 0.3, lux::Vector(0.0, 0.0, 0.0), lux::Vector(0.0, 1.0, 1.0));
	
	// cone
	std::shared_ptr<lux::ScalarField> co = std::make_shared<lux::SFCone>(lux::Vector(0.0, 0.0, 0.0), lux::Vector(1.0, 1.0, 0.0), 20.0, 0.5);

	// box
	std::shared_ptr<lux::ScalarField> b = std::make_shared<lux::SFBox>(0.001);

	// icosahedron
	std::shared_ptr<lux::ScalarField> i = std::make_shared<lux::SFIcosahedron>(lux::Vector(0.0, 0.0, 0.0));
	
	// ellipse
	std::shared_ptr<lux::ScalarField> e = std::make_shared<lux::SFEllipse>(0.2, 0.3, lux::Vector(0.1, 0.0, 0.0), lux::Vector(0.0, 1.0, 1.0));
	
	// stienerpatch
	std::shared_ptr<lux::ScalarField> st = std::make_shared<lux::SFStienerPatch>();

	// cylinder
	std::shared_ptr<lux::ScalarField> cy = std::make_shared<lux::SFCylinder>(lux::Vector(0.0, 0.0, 0.0), lux::Vector(0.0, 1.0, 0.0), 1.0, 0.2);

	std::shared_ptr<lux::ScalarField> f1 = std::make_shared<lux::SFIntersect>(sp, e);
	std::shared_ptr<lux::ScalarField> f2 = std::make_shared<lux::SFUnion>(f1, co);
	std::shared_ptr<lux::ScalarField> f3 = std::make_shared<lux::SFIntersect>(f2, e);
	std::shared_ptr<lux::ScalarField> f4 = std::make_shared<lux::SFIntersect>(cy, p1);
	std::shared_ptr<lux::ScalarField> f5 = std::make_shared<lux::SFIntersect>(f4, p2);
	

	lux::SFields scalarFields;
	scalarFields.push_back(f5->clone());

	std::vector<IMF::Rgba> exr;
	const int img_w = 1920/10;
	const int img_h = 1080/10;

	std::cout << "\t\t    ...\n\t         Keep Calm\n\t\t    and\n\t    Let The Rays March\n\t\t    ...\n\n|0%|==|==|==|==|==|==|==|==|==|==|==|100%|\n|0%|";
	for (int j = 0; j < img_h; ++j)
	{
		if (j % (img_h / 10) == 0)
			std::cout << "==|";
		for (int i = 0; i < img_w; ++i)
		{	
			lux::Color L = marchRays(camera, double(i) / (img_w - 1), double(j) / (img_h - 1), scalarFields);
			exr.push_back(IMF::Rgba(half(L[0]), half(L[1]), half(L[2]), half(L[3])));
		}
	}
	std::cout << "100%|\n";

	IMF::RgbaOutputFile file("out.exr", img_w, img_h, IMF::WRITE_RGBA);
	file.setFrameBuffer(const_cast<IMF::Rgba*>(exr.data()), 1, img_w);
	file.writePixels(img_h);

	return 0;
}