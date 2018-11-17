#pragma once
#include "ScalarField.hpp"
#include "Triangle.hpp"
#include <iostream>
#include <string>

#include "OBJ_Loader.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

// obj loader
lux::SField loadObjField(std::string fileName, Triangles& triangles);

// tiny obj
bool loadObj(std::string fileName, Triangles& triangles);
