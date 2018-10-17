#include "Main.hpp"

int main()
{
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();

	const int img_w = 960 / 1;
	const int img_h = 540 / 1;
	
	std::shared_ptr<Grid> g;

	std::shared_ptr<Light> key = std::make_shared<Light>(lux::Vector(-0.0, 0.5, 2.0), lux::Vector(-1, -1, -1), 200, 200, 200, 0.01);
	key->setColor(lux::Color(0.6, 0.2, 0.4, 0.5));
	std::shared_ptr<Light> fill = std::make_shared<Light>(lux::Vector(0.0, -2.4, 0.2), lux::Vector(-1, -1, -1), 200, 200, 200, 0.01);
	fill->setColor(lux::Color(0.1, 0.2, 0.6, 0.2));

	std::vector<std::shared_ptr<Light>> lights;// = { key, fill };
	render(img_w, img_h, camera, g, lights);

	int t;
	std::cin >> t;
	return 0;
}