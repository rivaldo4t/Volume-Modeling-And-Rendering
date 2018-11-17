#include "ObjFileLoaders.hpp"
#include "OBJ_Loader.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

lux::SField loadObjField(std::string fileName, Triangles& triangles)
{
	objl::Loader Loader;
	bool loadout = Loader.LoadFile(fileName);
	if (loadout)
	{
		for (unsigned int i = 0; i < Loader.LoadedMeshes.size(); i++)
		{
			objl::Mesh curMesh = Loader.LoadedMeshes[i];
			for (unsigned int j = 0; j < curMesh.Indices.size(); j += 3)
			{
				lux::Vector p0 = { curMesh.Vertices[curMesh.Indices[j + 0]].Position.X,
					curMesh.Vertices[curMesh.Indices[j + 0]].Position.Y,
					curMesh.Vertices[curMesh.Indices[j + 0]].Position.Z };
				lux::Vector p1 = { curMesh.Vertices[curMesh.Indices[j + 1]].Position.X,
					curMesh.Vertices[curMesh.Indices[j + 1]].Position.Y,
					curMesh.Vertices[curMesh.Indices[j + 1]].Position.Z };
				lux::Vector p2 = { curMesh.Vertices[curMesh.Indices[j + 2]].Position.X,
					curMesh.Vertices[curMesh.Indices[j + 2]].Position.Y,
					curMesh.Vertices[curMesh.Indices[j + 2]].Position.Z };

				//triangles.push_back(std::make_shared<Triangle>(p0 / 20.0, p1 / 20.0, p2 / 20.0));
				triangles.push_back(std::make_shared<Triangle>(p0, p1, p2));
			}
		}
	}
	std::cout << "Number of triangles loaded: " << triangles.size() << std::endl;

	int i = 0;
	// negative normal
	//lux::SField obj = std::make_shared<lux::SFPlane>(triangles[i]->p0, -(triangles[i]->n).unitvector());
	/*for (i = 1; i < triangles.size(); ++i)
	{
	lux::SField p = std::make_shared<lux::SFPlane>(triangles[i]->p0, -(triangles[i]->n).unitvector());
	obj = std::make_shared<lux::SFIntersect>(obj, p);
	}
	obj = std::make_shared<lux::SFScale>(obj, 0.5);
	obj = std::make_shared<lux::SFRotate>(obj, lux::Vector(1.0, 1.0, 0.0), 30);*/

	lux::SField obj;
	return obj;
}

bool loadObj(std::string fileName, Triangles& triangles)
{
	std::cout << "Loading " << fileName << std::endl;

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, fileName.c_str(), "models/", true);
	if (!ret)
	{
		printf("Failed to load/parse .obj.\n");
		if (!err.empty())
			std::cerr << err << std::endl;
		return false;
	}

	// For each shape
	for (size_t i = 0; i < shapes.size(); i++)
	{
		size_t index_offset = 0;

		// For each face
		for (size_t f = 0; f < shapes[i].mesh.num_face_vertices.size(); f++)
		{
			size_t fnum = shapes[i].mesh.num_face_vertices[f];
			assert(fnum == 3); // triangulate = true

			lux::Vector p0, p1, p2;
			//lux::Vector n0, n1, n2;
			tinyobj::index_t idx;

			idx = shapes[i].mesh.indices[index_offset + 0];
			p0 = { attrib.vertices[3 * idx.vertex_index + 0],
				attrib.vertices[3 * idx.vertex_index + 1],
				attrib.vertices[3 * idx.vertex_index + 2] };
			/*n0 = { attrib.normals[3 * idx.normal_index + 0],
			attrib.normals[3 * idx.normal_index + 1],
			attrib.normals[3 * idx.normal_index + 2] };*/

			idx = shapes[i].mesh.indices[index_offset + 1];
			p1 = { attrib.vertices[3 * idx.vertex_index + 0],
				attrib.vertices[3 * idx.vertex_index + 1],
				attrib.vertices[3 * idx.vertex_index + 2] };
			/*n1 = { attrib.normals[3 * idx.normal_index + 0],
			attrib.normals[3 * idx.normal_index + 1],
			attrib.normals[3 * idx.normal_index + 2] };*/

			idx = shapes[i].mesh.indices[index_offset + 2];
			p2 = { attrib.vertices[3 * idx.vertex_index + 0],
				attrib.vertices[3 * idx.vertex_index + 1],
				attrib.vertices[3 * idx.vertex_index + 2] };
			/*n2 = { attrib.normals[3 * idx.normal_index + 0],
			attrib.normals[3 * idx.normal_index + 1],
			attrib.normals[3 * idx.normal_index + 2] };*/

			Triangle t(p0 / 2, p1 / 2, p2 / 2);
			//t.n = (n0 + n1 + n2) / 3;
			triangles.push_back(std::make_shared<Triangle>(t));

			index_offset += fnum;
		}
	}

	std::cout << "Load successful\n";

	return true;
}