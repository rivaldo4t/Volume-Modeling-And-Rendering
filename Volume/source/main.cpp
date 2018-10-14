#include "Main.hpp"

int main()
{
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();

	const int img_w = 960 / 1;
	const int img_h = 540 / 1;
	
	FSPN f;
	std::shared_ptr<StampedNoise> g1 = std::make_shared<StampedNoise>(lux::Vector(0.0, 0.0, 0.0), 0.8, 1, f,
		lux::Vector(-1, -1, -1), 200, 200, 200, 0.01);
	//g1->computeNoise();

	std::shared_ptr<Grid> g = g1;
	std::vector<std::shared_ptr<Light>> lights;// = { key, rim };
	render(img_w, img_h, camera, g, lights);

	int t;
	std::cin >> t;
	return 0;
}