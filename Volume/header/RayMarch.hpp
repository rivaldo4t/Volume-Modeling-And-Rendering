#pragma once
#include <vector>
#include <memory>
#include "Vector.hpp"
#include "Camera.hpp"
#include "ScalarField.hpp"
#include "ColorField.hpp"
#include "Grid.hpp"
#include "Light.hpp"
#include "StampedNoise.hpp"
#include "Wisps.hpp"
#include "PyroclasticField.hpp"

// for scalar fields
void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, const lux::SField& sfield, const lux::CField& cfield);
lux::Color marchRays(lux::Vector pos, lux::Vector dir, const lux::SField& field, const lux::CField& c);
double marchRaysDSM(lux::Vector pos, lux::Vector lightPos, const lux::SField& density);

// for grids
// removed const - fix it!
void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, std::shared_ptr<Grid>& g, std::vector<std::shared_ptr<Light>>& lights);
lux::Color marchRays(lux::Vector pos, lux::Vector dir, const std::shared_ptr<Grid>& g, const std::vector<std::shared_ptr<Light>>& lights);
double marchToLight(lux::Vector pos, lux::Vector lightPos, const std::shared_ptr<Grid>& g);