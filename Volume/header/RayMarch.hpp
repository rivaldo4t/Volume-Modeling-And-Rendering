#pragma once
#include <vector>
#include <memory>
#include "Vector.hpp"
#include "Camera.hpp"
#include "ScalarField.hpp"
#include "ColorField.hpp"

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, const lux::SField& sfield, const lux::CField& cfield);

lux::Color marchRays(lux::Vector pos, lux::Vector dir, const lux::SField& field, const lux::CField& c);

double marchRaysDSM(lux::Vector pos, lux::Vector lightPos, const lux::SField& density);