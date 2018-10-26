#pragma once
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include "ScalarField.hpp"
#include "ColorField.hpp"
#include "Camera.hpp"
#include "RayMarch.hpp"
#include "Light.hpp"
#include "Triangle.hpp"
#include "Grid.hpp"
#include "StampedNoise.hpp"

#include "OBJ_Loader.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"


std::pair<lux::SField, lux::CField> getHumanoid();
// obj loader
lux::SField loadObjField(std::string fileName, Triangles& triangles);
// tiny obj
bool loadObj(std::string fileName, Triangles& triangles);