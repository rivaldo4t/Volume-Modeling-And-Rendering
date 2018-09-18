#pragma once
#include <vector>
#include <memory>
#include "Vector.hpp"
#include "Camera.hpp"
#include "ScalarField.hpp"

lux::Color marchRays(std::shared_ptr<Camera>& camera, double x, double y, lux::SFields& scalarFields);