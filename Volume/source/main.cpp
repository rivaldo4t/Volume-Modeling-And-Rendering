#include "Main.hpp"

int main()
{
	// camera
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();

	// render size
	const int img_w = 1920/10;
	const int img_h = 1080/10;
	
	// lights
#if 0
	{
		Light key(lux::Vector(0.0, 3.0, 2.0));
		Grid g(lux::Vector(-5, -5, -5), 800, 800, 800, 0.01);
		g.readlevelSet("D:/temp/vol/level_humanoid.dat");
		key.computeDSM2(g);
		key.writeDSM("D:/temp/vol/key_level_humanoid.dat");
		//key.readDSM("D:/temp/vol/key.dat");
	}

	{
		Light rim(lux::Vector(0.0, 0.0, -4.0));
		rim.computeDSM2(fields.first);
		rim.writeDSM("D:/temp/vol/rim_level.dat");
		// rim.readDSM("D:/temp/vol/rim.dat");
	}

	{
		Light fill(lux::Vector(0.0, -4.0, 0.0));
		fill.computeDSM2(fields.first);
		fill.writeDSM("D:/temp/vol/fill_level.dat");
		// fill.readDSM("D:/temp/vol/fill.dat");
	}
#endif
	// scalar and color field
	// auto fields = getHumanoid();

#if 0
	{
		Triangles triangles;
		auto obj = loadObjField("models/cleanteapot.obj", triangles);
		Grid g(lux::Vector(-2, -2, -2), 400, 400, 400, 0.01);
		g.levelSet(triangles);
		g.writelevelSet("D:/temp/vol/level_cleanteapot.dat");
		//g.readlevelSet("D:/temp/vol/level1.dat");
		std::cout << "cleanteapot levelset genereted\n";
	}
	{
		Triangles triangles;
		auto obj = loadObjField("models/cleanbunny.obj", triangles);
		Grid g(lux::Vector(-1, -1, -1), 400, 400, 400, 0.005);
		g.levelSet(triangles);
		g.writelevelSet("D:/temp/vol/level_cleanbunny.dat");
		//g.readlevelSet("D:/temp/vol/level1.dat");
		std::cout << "cleanbunny levelset genereted\n";
	}
	{
		auto fields = getHumanoid();
		Triangles triangles;
		Grid g(lux::Vector(-5, -5, -5), 800, 800, 800, 0.01);
		g.stamp(fields.first);
		g.writelevelSet("D:/temp/vol/level_humanoid.dat");
		//g.readlevelSet("D:/temp/vol/level1.dat");
		std::cout << "humanoid levelset genereted\n";
	}
#endif

	Grid g(lux::Vector(-5, -5, -5), 800, 800, 800, 0.01);
	g.readlevelSet("D:/temp/vol/level_humanoid.dat");

	//auto fields = getHumanoid();
	//Grid g(lux::Vector(-2, -2, -2), 40, 40, 40, 0.1);
	//g.stamp(fields.first);*/
	/*Triangles triangles;
	auto obj = loadObjField("models/cube1.obj", triangles);
	Grid g(lux::Vector(-1, -1, -1), 20, 20, 20, 0.1);
	g.levelSet(triangles);*/
	//g.readlevelSet("D:/temp/vol/level_humanoid.dat");
	//render(img_w, img_h, camera, fields.first, fields.second);
	render2(img_w, img_h, camera, g);
	int t;
	std::cin >> t;
	return 0;
}