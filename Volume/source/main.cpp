#include "Main.hpp"

int main()
{
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();

	const int img_w = 960 / 1;
	const int img_h = 540 / 1;
	
	std::shared_ptr<Grid> g;
	std::vector<std::shared_ptr<Light>> lights;// = { key, rim };
	render(img_w, img_h, camera, g, lights);

	int t;
	std::cin >> t;
	return 0;
}