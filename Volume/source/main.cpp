#include <iostream>
#include <memory>
#include <vector>

#include "ScalarField.hpp"
#include "Camera.hpp"
#include "RayMarch.hpp"

std::pair<lux::SField, lux::CField> getHumanoid()
{
	// stienerpatch
	lux::SField st_temp = std::make_shared<lux::SFStienerPatch>();
	lux::SField st = std::make_shared<lux::SFScale>(st_temp, 0.5);

	// Humanoid
	lux::SField ring_temp = std::make_shared<lux::SFTorus>(0.6, 0.1, lux::Vector(0.0, 0.0, 0.0), lux::Vector(0.0, 1.0, 0.0));
	lux::SField ring_temp_2 = std::make_shared<lux::SFScale>(ring_temp, 0.3);
	lux::SField ring = std::make_shared<lux::SFTranslate>(ring_temp_2, lux::Vector(0.0, 0.3, 0.0));
	lux::SField head = std::make_shared<lux::SFSphere>(lux::Vector(0.0, 0.1, 0.0), 0.08);
	lux::SField horn_1 = std::make_shared<lux::SFCone>(lux::Vector(0.18, 0.27, 0.0), lux::Vector(-1.0, -1.0, 0.0), 25.0, 0.12);
	lux::SField horn_2 = std::make_shared<lux::SFCone>(lux::Vector(-0.18, 0.27, 0.0), lux::Vector(1.0, -1.0, 0.0), 25.0, 0.12);
	lux::SField horns = std::make_shared<lux::SFUnion>(horn_1, horn_2);
	lux::SField body_temp = std::make_shared<lux::SFEllipse>(0.3, 0.18, lux::Vector(0.0, 0.0, 0.0), lux::Vector(1.0, 0.0, 0.0));
	lux::SField body_temp_2 = std::make_shared<lux::SFTranslate>(body_temp, lux::Vector(0.0, 0.0, -0.15));
	lux::SField body = std::make_shared<lux::SFCutout>(body_temp_2, head);
	lux::SField body_head = std::make_shared<lux::SFUnion>(body, head);
	lux::SField body_head_horns = std::make_shared<lux::SFUnion>(body_head, horns);

	lux::SField torso_temp = std::make_shared<lux::SFBox>(1);
	lux::SField torso_temp_2 = std::make_shared<lux::SFScale>(torso_temp, 0.15);
	lux::SField torso_temp_3 = std::make_shared<lux::SFTranslate>(torso_temp_2, lux::Vector(0.0, -0.32, -0.15));
	lux::SField torso = std::make_shared<lux::SFShell>(torso_temp_3, 0.5);
	lux::SField body_head_horns_torso = std::make_shared<lux::SFUnion>(body_head_horns, torso_temp_3);

	double h = 0.6;
	lux::Vector cyAxis = lux::Vector(0, 1, 0).unitvector();
	lux::SField cy = std::make_shared<lux::SFCylinder>(lux::Vector(0.0, 0.0, 0.0), cyAxis, h, 0.1);
	lux::SField p1 = std::make_shared<lux::SFPlane>(cyAxis * h * 0.5, -cyAxis);
	lux::SField p2 = std::make_shared<lux::SFPlane>(-cyAxis * h * 0.5, cyAxis);
	lux::SField leg_temp = std::make_shared<lux::SFIntersect>(cy, p1);
	lux::SField leg_temp_2 = std::make_shared<lux::SFIntersect>(leg_temp, p2);
	lux::SField exp_1 = std::make_shared<lux::SFScale>(leg_temp_2, 0.3);

	cyAxis = lux::Vector(0, 1, 0).unitvector();
	cy = std::make_shared<lux::SFCylinder>(lux::Vector(0.0, 0.0, 0.0), cyAxis, h * 3, 0.1);
	p1 = std::make_shared<lux::SFPlane>(cyAxis * h * 3 * 0.5, -cyAxis);
	p2 = std::make_shared<lux::SFPlane>(-cyAxis * h * 3 * 0.5, cyAxis);
	leg_temp = std::make_shared<lux::SFIntersect>(cy, p1);
	leg_temp_2 = std::make_shared<lux::SFIntersect>(leg_temp, p2);
	lux::SField exp_2 = std::make_shared<lux::SFScale>(leg_temp_2, 0.3);
	lux::SField exp_3 = std::make_shared<lux::SFRotate>(exp_2, lux::Vector(0, 0, 1), 45);
	lux::SField exp_4 = std::make_shared<lux::SFTranslate>(exp_3, lux::Vector(-0.17, -0.265, 0));

	lux::SField exp_5 = std::make_shared<lux::SFUnion>(exp_1, exp_4);
	lux::SField exp_6 = std::make_shared<lux::SFRotate>(exp_5, lux::Vector(0, 0, 1), -120);

	lux::SField exp_leg_1 = std::make_shared<lux::SFRotate>(exp_6, lux::Vector(0, 1, 0), -45);
	lux::SField leg_1 = std::make_shared<lux::SFTranslate>(exp_leg_1, lux::Vector(0.3, 0.1, -0.3));

	lux::SField exp_leg_2 = std::make_shared<lux::SFRotate>(exp_6, lux::Vector(0, 1, 0), 45);
	lux::SField leg_2 = std::make_shared<lux::SFTranslate>(exp_leg_2, lux::Vector(0.3, 0.1, -0.1));

	lux::SField exp_leg_3 = std::make_shared<lux::SFRotate>(exp_6, lux::Vector(0, 1, 0), 180 - 45);
	lux::SField leg_3 = std::make_shared<lux::SFTranslate>(exp_leg_3, lux::Vector(-0.3, 0.1, -0.1));

	lux::SField exp_leg_4 = std::make_shared<lux::SFRotate>(exp_6, lux::Vector(0, 1, 0), 180 + 45);
	lux::SField leg_4 = std::make_shared<lux::SFTranslate>(exp_leg_4, lux::Vector(-0.3, 0.1, -0.3));

	lux::SField exp_7 = std::make_shared<lux::SFUnion>(leg_1, leg_2);
	lux::SField exp_8 = std::make_shared<lux::SFUnion>(leg_3, leg_4);
	lux::SField legs = std::make_shared<lux::SFUnion>(exp_7, exp_8);

	lux::SField robe_temp = std::make_shared<lux::SFCone>(lux::Vector(0.0, 0.0, 0.0), lux::Vector(0.0, -1.0, 0.0), 30.0, 0.5);
	lux::SField robe = std::make_shared<lux::SFTranslate>(robe_temp, lux::Vector(0.0, -0.12, -0.15));
	lux::SField robe_legs = std::make_shared<lux::SFUnion>(robe, legs);

	lux::SField model = std::make_shared<lux::SFUnion>(body_head_horns_torso, robe_legs);
	lux::SField humanoid_temp = std::make_shared<lux::SFUnion>(model, ring);
	lux::SField ic_temp = std::make_shared<lux::SFIcosahedron>(lux::Vector(0.0, 0.0, 0.0));
	lux::SField ic_temp_2 = std::make_shared<lux::SFScale>(ic_temp, 0.02);
	lux::SField ic = std::make_shared<lux::SFTranslate>(ic_temp_2, lux::Vector(0.0, 0.1, 0.0));
	lux::SField humanoid_temp_2 = std::make_shared<lux::SFUnion>(ic, humanoid_temp);

	lux::SField humanoid = std::make_shared<lux::SFTranslate>(humanoid_temp_2, lux::Vector(0.0, 0.2, -0.4));

	lux::SField head_c = std::make_shared<lux::SFTranslate>(head, lux::Vector(0.0, 0.2, -0.4));
	lux::SField body_c = std::make_shared<lux::SFTranslate>(body, lux::Vector(0.0, 0.2, -0.4));
	lux::SField horns_c = std::make_shared<lux::SFTranslate>(horns, lux::Vector(0.0, 0.2, -0.4));
	lux::CField colorField = std::make_shared<lux::ColorField>(head_c, body_c, horns_c);

	//humanoid = st;
	return std::make_pair(humanoid, colorField);
}

int main()
{
	// camera
	std::shared_ptr<Camera> camera = std::make_shared<Camera>();

	//scalar and color field
	auto fields = getHumanoid();

	const int img_w = 1920/10;
	const int img_h = 1080/10;

	render(img_w, img_h, camera, fields.first, fields.second);

	return 0;
}