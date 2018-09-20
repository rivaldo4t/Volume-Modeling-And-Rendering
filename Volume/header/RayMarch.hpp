#pragma once
#include <vector>
#include <memory>
#include "Vector.hpp"
#include "Camera.hpp"
#include "ScalarField.hpp"

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, lux::SField sfield, lux::CField cfield);

lux::Color marchRays(std::shared_ptr<Camera>& camera, double x, double y, lux::SField field, lux::CField c);