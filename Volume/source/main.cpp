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
	
	// lights
#if 0
	auto fields = getHumanoid();
	Grid g(lux::Vector(-2, -2, -2), 500, 500, 500, 0.007);
	g.stamp(fields.first);
	g.writelevelSet("D:/temp/vol/level_humanoid_2.dat");
	std::cout << "humanoid levelset genereted\n";
	
	Light key(lux::Vector(0.0, 3.0, 3.0), lux::Vector(-2, -2, -2), 400, 400, 400, 0.01);
	key.computeDSM2(g);
	key.writeDSM("D:/temp/vol/key_level_humanoid_2.dat");
	/*key.readDSM("D:/temp/vol/key_level_humanoid_2.dat");*/
	
	Light rim(lux::Vector(0.0, 0.0, -3.0), lux::Vector(-2, -2, -2), 400, 400, 400, 0.01);
	rim.computeDSM2(g);
	rim.writeDSM("D:/temp/vol/rim_level_humanoid_2.dat");
	//rim.readDSM("D:/temp/vol/rim_level_humanoid_2.dat");

	Light fill(lux::Vector(0.0, -3.0, 0.0), lux::Vector(-2, -2, -2), 400, 400, 400, 0.01);
	fill.computeDSM2(g);
	fill.writeDSM("D:/temp/vol/fill_level_humanoid_2.dat");
	//fill.readDSM("D:/temp/vol/fill_level_humanoid_2.dat");
	
	std::vector<std::shared_ptr<Light>> lights;
	lights.emplace_back(std::make_shared<Light>(key));
	lights.emplace_back(std::make_shared<Light>(rim));
	lights.emplace_back(std::make_shared<Light>(fill));
	
	render2(img_w, img_h, camera, g, lights);
#endif

#if 0
	{
		Triangles triangles;
		auto obj = loadObjField("models/cleanbunny.obj", triangles);
		Grid g(lux::Vector(-0.8, -0.8, -0.8), 500, 500, 500, 0.0032);
		g.levelSet(triangles);
		g.writelevelSet("D:/temp/vol/level_cleanbunny_2.dat");
		std::cout << "cleanbunny levelset genereted\n";
		//g.readlevelSet("D:/temp/vol/level1.dat");
	}

	{
		Triangles triangles;
		auto obj = loadObjField("models/cleanteapot.obj", triangles);
		Grid g(lux::Vector(-1.5, -1.5, -1.5), 500, 500, 500, 0.006);
		g.levelSet(triangles);
		g.writelevelSet("D:/temp/vol/level_cleanteapot_2.dat");
		std::cout << "cleanteapot levelset genereted\n";
		//g.readlevelSet("D:/temp/vol/level1.dat");
	}
#endif

	/*Grid g(lux::Vector(-0.8, -0.8, -0.8), 500, 500, 500, 0.0032);
	g.readlevelSet("D:/temp/vol/level_cleanbunny_2.dat");
	std::vector<std::shared_ptr<Light>> lights;
	render2(img_w, img_h, camera, g, lights);*/

	Grid g(lux::Vector(-2, -2, -2), 500, 500, 500, 0.007);
	g.readlevelSet("D:/temp/vol/level_humanoid_2.dat");

	int t;
	std::cin >> t;
	return 0;
}