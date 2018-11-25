#pragma once
#include <vector>
#include <memory>
#include "Vector.hpp"
#include "Camera.hpp"
#include "ScalarField.hpp"
#include "ColorField.hpp"
#include "ScalarGrid.hpp"
#include "VectorGrid.hpp"
#include "Light.hpp"
#include "StampedNoise.hpp"
#include "Wisps.hpp"
#include "PyroclasticField.hpp"
#include "TerrainNoise.hpp"
#include "WarpField.hpp"
#include "AdvectedField.hpp"
#include "ObjFileLoaders.hpp"

void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, const lux::SField& sfield, const lux::CField& cfield);
lux::Color marchRays(lux::Vector pos, lux::Vector dir, const lux::SField& field, const lux::CField& c);
float marchRaysDSM(lux::Vector pos, lux::Vector lightPos, const lux::SField& density);

// with lights
// removed const - fix it!
void render(const int img_w, const int img_h, std::shared_ptr<Camera> camera, lux::SField& g, std::vector<std::shared_ptr<lux::Light>>& lights);
lux::Color marchRays(lux::Vector pos, lux::Vector dir, const lux::SField& g, const std::vector<std::shared_ptr<lux::Light>>& lights);
float marchToLight(lux::Vector pos, lux::Vector lightPos, const lux::SField& g);