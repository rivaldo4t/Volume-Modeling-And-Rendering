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

	//render(img_w, img_h, camera, fields.first, fields.second);
	
	/*Grid g2(lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	g2.stamp(getHumanoid().first);
	
	Triangles triangles;
	auto obj = loadObjField("models/cleanbunny.obj", triangles);
	Grid g(lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	g.levelSet(triangles);
	g.writeGrid("D:/temp/vol/level_cleanbunny.dat");*/
	
	/*Light key(lux::Vector(0.0, 3.0, 3.0), lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	key.color = lux::Color(0.8, 0.2, 0.2, 1.0);
	key.computeDSM(g);
	key.writeDSM("D:/temp/vol/key_level_cleanbunny.dat");
	
	Light rim(lux::Vector(0.0, 0.0, -3.0), lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	rim.color = lux::Color(0.2, 0.8, 0.2, 1.0);
	rim.computeDSM(g);
	rim.writeDSM("D:/temp/vol/rim_level_cleanbunny.dat");
	
	Light fill(lux::Vector(0.0, -3.0, 0.0), lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	fill.color = lux::Color(0.2, 0.2, 0.8, 1.0);
	fill.computeDSM(g);
	fill.writeDSM("D:/temp/vol/fill_level_cleanbunny.dat");*/

	/*Grid g;
	g.readGrid("D:/temp/vol/level_cleanbunny_fine.dat");
	std::shared_ptr<Grid> f = std::make_shared<Grid>(g);
	Light key;
	key.readDSM("D:/temp/vol/key_level_cleanbunny_fine.dat");
	Light rim;
	rim.readDSM("D:/temp/vol/rim_level_cleanbunny_fine.dat");
	Light fill;
	fill.readDSM("D:/temp/vol/fill_level_cleanbunny_fine.dat");*/
	
	std::vector<std::shared_ptr<Light>> lights;
	/*lights.emplace_back(std::make_shared<Light>(key));
	lights.emplace_back(std::make_shared<Light>(rim));
	lights.emplace_back(std::make_shared<Light>(fill));*/
	
	//std::shared_ptr<Grid> g1 = std::make_shared<Grid>(lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	/*std::shared_ptr<Grid> g1 = std::make_shared<Grid>(lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	g1->stamp(getHumanoid().first);

	std::shared_ptr<Grid> g2 = std::make_shared<Grid>(lux::Vector(-0.8, -0.8, -0.8), 16, 16, 16, 0.1);
	Triangles triangles;
	auto obj = loadObjField("models/cleanbunny.obj", triangles);
	g2->levelSet(triangles);

	std::shared_ptr<Grid> g4 = std::make_shared<GridScale>(g1, 0.7);
	std::shared_ptr<Grid> g5 = std::make_shared<GridTranslate>(g1, lux::Vector(1, 1, 0));
	std::shared_ptr<Grid> g3 = std::make_shared<GridUnion>(g5, g2);*/

	Grid g;
	g.readGrid("D:/temp/vol/level_cleanbunny_fine.dat");
	std::shared_ptr<Grid> g1 = std::make_shared<Grid>(g);
	std::shared_ptr<Grid> g2 = std::make_shared<GridScale>(g1, 0.2);
	std::shared_ptr<Grid> g3 = std::make_shared<GridTranslate>(g2, lux::Vector(-1.5, 0, 0));
	
	Grid f;
	f.readGrid("D:/temp/vol/level_cleanbunny_fine.dat");
	std::shared_ptr<Grid> g4 = std::make_shared<Grid>(f);
	std::shared_ptr<Grid> g5 = std::make_shared<GridScale>(g4, 0.2);
	std::shared_ptr<Grid> g6 = std::make_shared<GridTranslate>(g5, lux::Vector(1.5, 0, 0));

	std::shared_ptr<Grid> g7 = std::make_shared<GridUnion>(g3, g6);
	render2(img_w, img_h, camera, g7, lights);

	int t;
	std::cin >> t;
	return 0;
}