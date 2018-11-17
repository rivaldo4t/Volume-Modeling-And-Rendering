#pragma once
#include "ScalarField.hpp"
#include "Triangle.hpp"
#include <iostream>
#include <string>

// OBJ_Loader
lux::SField loadObjField(std::string fileName, Triangles& triangles);

// tiny_obj_loader
bool loadObj(std::string fileName, Triangles& triangles);
