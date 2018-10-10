#include "Main.hpp"

int main()
{
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();

	const int img_w = 1920/1;
	const int img_h = 1080/1;
	
	//auto fields = getHumanoid();
	//render(img_w, img_h, camera, fields.first, fields.second);

	/*std::shared_ptr<Grid> g = std::make_shared<Grid>(lux::Vector(-0.8, -0.8, -0.8), 320, 320, 320, 0.005);
	Triangles triangles;
	auto obj = loadObjField("models/cleanbunny.obj", triangles);
	g->levelSet(triangles);
	g->writeGrid("D:/temp/vol/level_cleanbunny.dat");

	std::shared_ptr<Light> key = std::make_shared<Light>(lux::Vector(-2.0, 1.0, 0.5), lux::Vector(-0.8, -0.8, -0.8), 320, 320, 320, 0.005);
	key->setColor(lux::Color(0.6, 0.2, 0.4, 0.8));
	key->computeDSM(g);
	key->writeDSM("D:/temp/vol/key_level_cleanbunny.dat");

	std::shared_ptr<Light> rim = std::make_shared<Light>(lux::Vector(0.5, 0.1, -3.0), lux::Vector(-0.8, -0.8, -0.8), 320, 320, 320, 0.005);
	rim->setColor(lux::Color(0.2, 0.4, 0.6, 0.2));
	rim->computeDSM(g);
	rim->writeDSM("D:/temp/vol/rim_level_cleanbunny.dat");

	std::shared_ptr<Light> fill = std::make_shared<Light>(lux::Vector(0.0, -2.0, 0.0), lux::Vector(-0.8, -0.8, -0.8), 320, 320, 320, 0.005);
	fill->setColor(lux::Color(0.2, 0.2, 0.2, 0.2));
	fill->computeDSM(g);
	fill->writeDSM("D:/temp/vol/fill_level_cleanbunny.dat");*/

	Grid g;
	g.readGrid("D:/temp/vol/level_cleanbunny.dat");
	std::shared_ptr<Grid> g1 = std::make_shared<Grid>(g);

	Light l;
	l.readDSM("D:/temp/vol/key_level_cleanbunny.dat");
	std::shared_ptr<Light> key = std::make_shared<Light>(l);

	Light l2;
	l2.readDSM("D:/temp/vol/rim_level_cleanbunny.dat");
	std::shared_ptr<Light> rim = std::make_shared<Light>(l2);

	Light l3;
	l3.readDSM("D:/temp/vol/fill_level_cleanbunny.dat");
	std::shared_ptr<Light> fill = std::make_shared<Light>(l3);

	std::vector<std::shared_ptr<Light>> lights = {key, rim, fill};
	render2(img_w, img_h, camera, g1, lights);

	int t;
	std::cin >> t;
	return 0;
}