#include "Main.hpp"

int main()
{
	// camera
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();

	// render size
	const int img_w = 1920/1;
	const int img_h = 1080/1;
	
	// scalar and color field
	//auto fields = getHumanoid();

	//render loop
	//render(img_w, img_h, camera, fields.first, fields.second);
	
	Triangles triangles;
	auto obj = loadObjField("models/cleanbunny.obj", triangles);
	Grid g(lux::Vector(-0.8, -0.8, -0.8), 320, 320, 320, 0.005);
	g.levelSet(triangles);
	g.writeGrid("D:/temp/vol/level_cleanbunny_fine.dat");
	
	Light key(lux::Vector(0.0, 3.0, 3.0), lux::Vector(-0.8, -0.8, -0.8), 160, 160, 160, 0.01);
	key.color = lux::Color(0.8, 0.2, 0.2, 1.0);
	key.computeDSM(g);
	key.writeDSM("D:/temp/vol/key_level_cleanbunny_fine.dat");
	
	Light rim(lux::Vector(0.0, 0.0, -3.0), lux::Vector(-0.8, -0.8, -0.8), 160, 160, 160, 0.01);
	rim.color = lux::Color(0.2, 0.8, 0.2, 1.0);
	rim.computeDSM(g);
	rim.writeDSM("D:/temp/vol/rim_level_cleanbunny_fine.dat");
	
	Light fill(lux::Vector(0.0, -3.0, 0.0), lux::Vector(-0.8, -0.8, -0.8), 160, 160, 160, 0.01);
	fill.color = lux::Color(0.2, 0.2, 0.8, 1.0);
	fill.computeDSM(g);
	fill.writeDSM("D:/temp/vol/fill_level_cleanbunny_fine.dat");
	
	/*std::vector<std::shared_ptr<Light>> lights;
	lights.emplace_back(std::make_shared<Light>(key));
	lights.emplace_back(std::make_shared<Light>(rim));
	lights.emplace_back(std::make_shared<Light>(fill));
	
	render2(img_w, img_h, camera, g, lights);*/

	int t;
	std::cin >> t;
	return 0;
}